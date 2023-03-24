#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <string>
#include <iostream>
#include <queue>
#include <algorithm>

#include "../include/PolarCode.h"
#include "../include/GaussianApproximation.h"
#include "../include/ScCrcAidedDecoder.h"
#include "../include/ScFlipProgDecoder.h"
#include "../include/Exceptions.h"
#include "../include/Domain.h"

#define DBL_MAX 1.7976931348623158e+308 
#define FROZEN_VALUE 0
#define ROOT_POSITION -1

ScFlipProgDecoder::ScFlipProgDecoder(PolarCode * codePtr,
	int level,
	double gammaLeft,
	double gammaRight,
	std::vector<double> omegaArr) : ScCrcAidedDecoder(codePtr) {

	std::size_t n = _codePtr->N();

	_gammaLeft = gammaLeft;
	_gammaRight = gammaRight;
	_levelMax = level;
	_omegaArr = omegaArr;

	_criticalSetTree = GetCriticalSetTree(_maskWithCrc, _levelMax);
	_subchannelsMeansGa = std::vector<double>(n, 0);
}

void ScFlipProgDecoder::SetSigma(double sigma) {
	BaseDecoder::SetSigma(sigma);

	GaussianApproximation ga(sigma);
	std::size_t n = _codePtr->N();
	for (std::size_t i = 0; i < n; i++)
	{
		_subchannelsMeansGa[i] = ga.GetMu(i + 1, n);
	}

	return;
}

// mask: 1 - black, 0 - white
std::vector<int> ScFlipProgDecoder::GetCriticalSet(std::vector<int> mask, int position) {
	std::vector<std::vector<int>> tree;
	std::size_t m = _codePtr->m();
	std::size_t n = _codePtr->N();
	std::vector<int> result;

	for (std::size_t i = 0; i < m + 1; i++)
	{
		std::vector<int> temp((int)1 << i, 0);
		tree.push_back(temp);
	}
	for (std::size_t i = position; i < n; i++)
	{
		tree[m][i] = mask[i];
	}

	for (int i = (int)m - 1; i >= 0; i--)
	{
		std::size_t levelSize = tree[i].size();
		for (int j = 0; j < levelSize; j++)
		{
			int left = j << 1;
			int right = left + 1;

			tree[i][j] = tree[i + 1][left] && tree[i + 1][right];
		}
	}

	for (std::size_t i = 0; i < m + 1; i++)
	{
		int levelSize = (int)tree[i].size();
		for (int j = 0; j < levelSize; j++)
		{
			if (!tree[i][j])
				continue;

			int criticalBit = j << (m - i);
			result.push_back(criticalBit);

			int left = j;
			int right = j;
			tree[i][j] = 0;
			for (std::size_t k = i + 1; k <= m; k++)
			{
				left = left << 1;
				right = (right << 1) + 1;

				for (std::size_t l = left; l <= right; l++)
					tree[k][l] = 0;
			}
		}
	}

	return result;
}

CriticalSetNode * ScFlipProgDecoder::GetCriticalSetTree(std::vector<int> mask, int levelMax) {
	std::queue<CriticalSetNode *> q;

	auto root = new CriticalSetNode();
	root->Bit = ROOT_POSITION;
	q.push(root);
	while (!q.empty()) {
		CriticalSetNode * cur = q.front();
		q.pop();

		auto criticalSet = GetCriticalSet(_maskWithCrc, cur->Bit + 1);
		for (std::size_t i = 0; i < criticalSet.size(); i++)
		{
			auto childNode = new CriticalSetNode();
			childNode->Bit = criticalSet[i];
			
			if (cur->Bit != ROOT_POSITION) {
				childNode->Path = cur->Path;
				childNode->Path.push_back(cur->Bit);
			}
				
			cur->Children.push_back(childNode);
			
			if (childNode->Path.size() < levelMax - 1) // indices of level from 1
				q.push(childNode);
		}
	}

	return root;
}

// with using means after gaussian approxiamtion procedure
std::vector<CriticalSetNode *> ScFlipProgDecoder::SortCriticalNodes(std::vector<CriticalSetNode *> criticalNodes, std::vector<double> llrs) {
	std::size_t length = criticalNodes.size();
	std::vector<CriticalSetNode *> result(length, 0);
	std::vector<double> sortingMetric(length, 0);
	
	for (std::size_t i = 0; i < length; i++)
	{
		int criticalInd = criticalNodes[i]->Bit;
		sortingMetric[i] = fabs(llrs[criticalInd]) / _subchannelsMeansGa[criticalInd];
	}

	for (std::size_t i = 0; i < length; i++)
	{
		auto minIt = std::min_element(sortingMetric.begin(), sortingMetric.end());
		int minInd = (int)std::distance(sortingMetric.begin(), minIt);
		result[i] = criticalNodes[minInd];

		criticalNodes.erase(criticalNodes.begin() + minInd);
		sortingMetric.erase(sortingMetric.begin() + minInd);
	}

	return result;
}

bool ScFlipProgDecoder::NoChild(CriticalSetNode * node, std::vector<double> inLlr) {
	int position = node->Bit;
	std::size_t n = _codePtr->N();
	std::size_t criticalSetSize = node->Children.size();
	
	if (position >= n)
		return true;

	if (criticalSetSize <= 0)
		return true;

	std::size_t level = node->Path.size() + 1;

	if (_omegaArr[level] >= 1.0)
		return false;

	std::vector<int> maskWithoutCriticalSet = _maskWithCrc;
	for (std::size_t i = 0; i < criticalSetSize; i++)
	{
		maskWithoutCriticalSet[node->Children[i]->Bit] = 0;
	}
	maskWithoutCriticalSet[node->Bit] = 0;

	int n1 = 0;
	int n2 = 0;
	for (std::size_t i = position; i < n; i++)
	{
		if (maskWithoutCriticalSet[i]) {
			n1++;

			double mu = _subchannelsMeansGa[i];
			double sigma = sqrt(2 * mu);
			if (inLlr[i] < mu - _gammaLeft * sigma)
				n2++;
		}
	}

	return ((double)n2 ) / n1 > _omegaArr[level];
}
bool ScFlipProgDecoder::NotSelect(CriticalSetNode * node, std::vector<double> inLlr) {
	int position = node->Bit;
	std::size_t level_dec = node->Path.size();

	if (_omegaArr[level_dec] >= 1.0)
		return false;

	//std::cout << level_dec << std::endl;

	double mu = _subchannelsMeansGa[position];
	double sigma = sqrt(2 * mu);
	
	return inLlr[position] > mu + _gammaRight * sigma;
}

std::vector<int> ScFlipProgDecoder::Decode(std::vector<double> inLlr) {
	std::size_t n = inLlr.size();
	std::size_t m = _codePtr->m();
	
	DecodeInternal(inLlr);

	if (!IsCrcPassed(_x)) {
		std::queue<CriticalSetNode *> q;
		q.push(_criticalSetTree);

		bool exitFlag = false;
		while (!q.empty() && !exitFlag) {
			auto currentNode = q.front();
			q.pop();

			auto criticalNodes = currentNode->Children;
			std::vector<CriticalSetNode *> suspectedNodes = SortCriticalNodes(criticalNodes, _beliefTree[m]);

			for (std::size_t i = 0; i < suspectedNodes.size(); i++)
			{
				int bitPosition = suspectedNodes[i]->Bit;
				if (NotSelect(suspectedNodes[i], inLlr))
					continue;

				std::size_t pathSize = suspectedNodes[i]->Path.size();
				for (std::size_t j = 0; j < pathSize; j++)
				{
					int prevBitPosition = suspectedNodes[i]->Path[j];
					int nextPrevBitPosition = (j != pathSize - 1) ? suspectedNodes[i]->Path[j + 1] : bitPosition;
					_x[prevBitPosition] = !_x[prevBitPosition];
					DecodeFromTo(prevBitPosition, nextPrevBitPosition);
				}
				
				_x[bitPosition] = !_x[bitPosition];
				DecodeFromTo(bitPosition, (int)n);
				
				if (IsCrcPassed(_x)) {
					exitFlag = true;
					break;
				}
				
				int minPosition = suspectedNodes[i]->Path.empty() ? bitPosition : suspectedNodes[i]->Path[0];
				_x[minPosition] = !_x[minPosition];
				DecodeFromTo(minPosition, (int)n);
				
				if (!NoChild(suspectedNodes[i], inLlr))
					q.push(suspectedNodes[i]);
			}
		}

	}

	return TakeResult();
}