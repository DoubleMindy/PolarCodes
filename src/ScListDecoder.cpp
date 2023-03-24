#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <string>
#include <iostream>
#include <queue>
#include <algorithm>

#include "../include/ScListDecoder.h"

#define DBL_MAX 1.7976931348623158e+308 
#define FROZEN_VALUE 0

ScListDecoder::ScListDecoder(PolarCode * codePtr, int L) : ScCrcAidedDecoder(codePtr) {
	_L = L;
	std::size_t m = _codePtr->m();
	std::size_t n = _codePtr->N();
	_treeHeight = m + 1;

	std::vector<double> b(n, -10000.0);
	std::vector<int> u(n, -1);

	std::vector<std::vector<double>> beliefTree;
	std::vector<std::vector<int>> uhatTree;

	for (std::size_t i = 0; i < _treeHeight; i++)
	{
		beliefTree.push_back(b);
		uhatTree.push_back(u);
	}

	for (std::size_t i = 0; i < _L; i++)
	{
		_beliefTrees.push_back(beliefTree);
		_uhatTrees.push_back(uhatTree);
	}

	for (std::size_t i = 0; i < _L; i++)
	{
		_candidates.push_back(u);
	}

	_metrics = std::vector<double>(_L, 0);
	
	// optimization allocation
	_areTakenOne = std::vector<bool>(_L, 0);
	_areTakenZero = std::vector<bool>(_L, 0);
}

double ScListDecoder::StepMetric(double belief, int decision) {
#ifdef DOMAIN_LLR
	
#ifdef MINSUM

	return (belief < 0 && decision == 0 || belief > 0 && decision == 1) ? -fabs(belief) : 0;
#else
	const double limit = 10000000.0;

	if (belief < -limit)
		belief = -limit;

	if (belief > limit)
		belief = limit;

	double p0_methric = -log(1 + exp(-belief));
	double p1_methric = -log(1 + exp(belief));
	return (decision) ? p1_methric : p0_methric;

#endif // MINSUM

#elif DOMAIN_P1
	return log((decision) ? belief : 1 - belief);
#endif // DOMAIN
}

void ScListDecoder::PassDownList(std::size_t iter) {
	std::size_t m = _codePtr->m();
	std::size_t n = _codePtr->N();

	std::size_t iterXor;
	std::size_t level;
	if (iter) {
		iterXor = iter ^ (iter - 1);
		level = m - FirstBitPos(iterXor);
	}
	else {
		level = 0;
	}

	//std::vector<int> binaryIter(m - level, 0);
	int size = (int)(m - level);
	std::size_t iterCopy = iter;
	for (int i = size - 1; i >= 0; i--)
	{
		_binaryIter[i] = iterCopy % 2;
		iterCopy = iterCopy >> 1;

	}

	std::size_t length = (std::size_t)1 << (m - level - 1);
	for (std::size_t i = level; i < m; i++)
	{
		std::size_t ones = ~0u;
		std::size_t offset = iter & (ones << (m - i));

		for (std::size_t j = 0; j < _L; j++)
		{
			if (!_binaryIter[i - level]) {
				FillLeftMessageInTree(_beliefTrees[j][i].begin() + offset,
					_beliefTrees[j][i].begin() + offset + length,
					_beliefTrees[j][i + 1].begin() + offset,
					length);
			}
			else {

				FillRightMessageInTree(_beliefTrees[j][i].begin() + offset,
					_beliefTrees[j][i].begin() + offset + length,
					_uhatTrees[j][i + 1].begin() + offset,
					_beliefTrees[j][i + 1].begin() + offset + length,
					length);
			}
		}
		
		length = length / 2;
	}

}

void ScListDecoder::PassUpList(std::size_t iter) {
	std::size_t m = _codePtr->m();
	std::size_t iterCopy = iter;

	std::size_t bit = iterCopy % 2;
	std::size_t length = 1;
	std::size_t level = m;
	std::size_t offset = iter;
	while (bit != 0)
	{
		offset -= length;
		for (std::size_t j = 0; j < _L; j++)
		{
			for (std::size_t i = 0; i < length; i++)
			{
				_uhatTrees[j][level - 1][offset + i] = _uhatTrees[j][level][offset + i] ^ _uhatTrees[j][level][offset + length + i];

			}
			for (std::size_t i = 0; i < length; i++)
			{
				_uhatTrees[j][level - 1][offset + length + i] = _uhatTrees[j][level][offset + length + i];

			}

		}
		
		iterCopy = iterCopy >> 1;
		bit = iterCopy % 2;
		length *= 2;
		level -= 1;

	}

}

void ScListDecoder::DecodeListInternal(std::vector<double> inLlr) {
	std::size_t n = inLlr.size();
	std::size_t m = _codePtr->m();

	// Fill each tree in the forrest with input llrs
	for (std::size_t j = 0; j < _L; j++) {
		for (std::size_t i = 0; i < n; i++)
			_beliefTrees[j][0][i] = inLlr[i];

		_metrics[j] = 0;

	}
		
	int logL = (int)FirstBitPos(_L) - 1;

	// first log(_L) bits
	std::size_t i_all = 0; // number of bit (j - number of condidate in list)
	std::size_t i_unfrozen = 0;
	while (i_unfrozen < logL)
	{
		PassDownList(i_all);

		if (_maskWithCrc[i_all]) {
			int value = 0;
			for (std::size_t j = 0; j < _L; j++) {
				_candidates[j][i_all] = value;
				if ((j + 1) % (1 << i_unfrozen) == 0) // all paths at the logL first steps 
					value = !value;

			}
			i_unfrozen++;

		}
		else {
			for (std::size_t j = 0; j < _L; j++) {
				_candidates[j][i_all] = FROZEN_VALUE;

			}
		}

		for (std::size_t j = 0; j < _L; j++) {
			_uhatTrees[j][m][i_all] = _candidates[j][i_all];
			_metrics[j] += StepMetric(_beliefTrees[j][m][i_all], _candidates[j][i_all]);

		}
			
		PassUpList(i_all);

		i_all++;

	}
	while (i_all < n)
	{
		PassDownList(i_all);

		if (_maskWithCrc[i_all]) {
			FillListMask(i_all);
			
			for (std::size_t j = 0, i = 0; i < _L; i++) {
				// add new path
				if (_areTakenOne[j] && _areTakenZero[j]) {
					_candidates[j][i_all] = 1;
					_beliefTrees.push_back(_beliefTrees[j]);
					_metrics.push_back(_metrics[j]);
					_uhatTrees.push_back(_uhatTrees[j]);
					_candidates.push_back(_candidates[j]);

					_candidates[j][i_all] = 0;
					_areTakenOne.push_back(false);
					_areTakenZero.push_back(false);

					j++;
					continue;
				}
				if (_areTakenZero[j]) {
					_candidates[j][i_all] = 0;
					j++;
					continue;
				}
				if (_areTakenOne[j]) {
					_candidates[j][i_all] = 1;
					j++;
					continue;
				}
				
				// delete path
				_beliefTrees.erase(_beliefTrees.begin() + j);
				_metrics.erase(_metrics.begin() + j);
				_uhatTrees.erase(_uhatTrees.begin() + j);
				_candidates.erase(_candidates.begin() + j);
				_areTakenZero.erase(_areTakenZero.begin() + j);
				_areTakenOne.erase(_areTakenOne.begin() + j);
			}
		}
		else {
			for (std::size_t j = 0; j < _L; j++) {
				_candidates[j][i_all] = FROZEN_VALUE;

			}
				
		}

		for (std::size_t j = 0; j < _L; j++) {
			_uhatTrees[j][m][i_all] = _candidates[j][i_all];
			_metrics[j] += StepMetric(_beliefTrees[j][m][i_all], _candidates[j][i_all]);

		}
			
		PassUpList(i_all);

		i_all++;

	}

	_operationsCount.Iterations += _L;
	_operationsCount.Normilizer++;

	return;
}

void ScListDecoder::FillListMask(std::size_t iter) {
	std::vector<int> indices(2 * _L, 0);
	std::vector<double> metricsNew(2 * _L, 0);

	std::size_t m = _codePtr->m();

	for (int j = 0; j < _L; j++)
	{
		_areTakenOne[j] = false;
		_areTakenZero[j] = false;

		indices[j] = j;
		indices[j + _L] = j + _L;

		double step0 = StepMetric(_beliefTrees[j][m][iter], 0);
		double step1 = StepMetric(_beliefTrees[j][m][iter], 1);

		metricsNew[j] = _metrics[j] + step0;
		metricsNew[j + _L] = _metrics[j] + step1;

	}

	for (std::size_t i = 0; i < _L; i++)
	{
		auto maxIt = std::max_element(metricsNew.begin(), metricsNew.end());
		int maxInd = (int)std::distance(metricsNew.begin(), maxIt);

		if (indices[maxInd] >= _L)
			_areTakenOne[indices[maxInd] - _L] = true;
		else
			_areTakenZero[indices[maxInd]] = true;

		indices.erase(indices.begin() + maxInd);
		
		metricsNew.erase(metricsNew.begin() + maxInd);

	}
}

std::vector<int> ScListDecoder::TakeListResult() {
	std::vector<int> result(_codePtr->k(), 0);
	std::vector<int> candidate(_codePtr->N(), 0);
	std::vector<int> codewordBits = _codePtr->UnfrozenBits();
	std::vector<double> metrics = _metrics;

	int maxInd = -1;
	std::size_t j = 0;
	for (; j < _L; j++)
	{
		auto maxIt = std::max_element(_metrics.begin(), _metrics.end());
		maxInd = (int)std::distance(_metrics.begin(), maxIt);

		if (IsCrcPassed(_candidates[maxInd]))
			break;
		
		_metrics[maxInd] = -100000000.0;

	}
	if (j < _L)
		candidate = _candidates[maxInd];

	for (std::size_t i = 0; i < codewordBits.size(); i++)
	{
		result[i] = candidate[codewordBits[i]];
	}
	return result;
}


std::vector<int> ScListDecoder::Decode(std::vector<double> inLlr) {

	DecodeListInternal(inLlr);

	return TakeListResult();
}