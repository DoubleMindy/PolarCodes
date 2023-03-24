#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <string>
#include <queue>
#include <iostream>
#include <algorithm>

#include "../include/PolarCode.h"
#include "../include/ScDecoder.h"
#include "../include/ScDecoderTreeMaker.h"
#include "../include/Exceptions.h"
#include "../include/Domain.h"
#include "../include/GaussianApproximation.h"

#define DBL_MAX 1.7976931348623158e+308 
#define FROZEN_VALUE 0

// TODO
//	input output into another cpp file
//  approximation sigma is fixed de-facto

ScDecoderTreeMaker::ScDecoderTreeMaker(PolarCode * codePtr, double approximationSigma) : ScDecoder(codePtr) {
	std::size_t n = _codePtr->N();

	approximationSigma = sqrt(pow(10, -1.0 / 10.0));
	GaussianApproximation ga(approximationSigma);
	_p = std::vector<double>(n, 0);
	for (std::size_t i = 0; i < n; i++)
	{
		_p[i] = ga.GetChannelErrorProbability(i + 1, n);
	}
}

std::size_t GetFirstDistinctBit(std::size_t first, std::size_t second) {
	std::size_t x = first ^ second;

	std::size_t m = 0;
	while (x > 0)
	{
		x = x >> 1;
		m++;
	}

	return m;
}

// Fill
void FillBinary(int number, std::vector<int> & binary, std::size_t lastDigitsCount) {
	std::size_t k = binary.size();
	for (std::size_t i = 0; i < lastDigitsCount; i++)
	{
		int digit = number % 2;
		number = number >> 1;
		binary[k - i - 1] = digit;
	}

	return;
}

struct PathNode {
	std::string binaryPath;
	std::string description;
};

std::string ScDecoderTreeMaker::GetPathInfo() {
	return _pathTrace;
}

bool CompareTree(PathNode first, PathNode second) {
	if (first.binaryPath.size() > second.binaryPath.size())
		return false;

	if (first.binaryPath.size() < second.binaryPath.size())
		return true;

	return first.binaryPath < second.binaryPath;
}

std::vector<int> ScDecoderTreeMaker::Decode(std::vector<double> inP1) {
	_pathTrace = "";
	std::vector<PathNode> unorderedTree;

	std::size_t n = inP1.size();
	std::size_t m = _codePtr->m();
	std::size_t k = _codePtr->k();
	for (std::size_t i = 0; i < n; i++)
	{
		_beliefTree[0][i] = inP1[i];
	}

	std::vector<double> beta(k, 0.0); // only for unfrozen bits
	std::vector<double> metrics(n, 0); // for all bits
	std::vector<int> A = _codePtr->UnfrozenBits(); // info set
	
	std::vector<int> binary(k, 0);

	// fill first frozen bits
	for (std::size_t i = 0; i < A[0]; i++)
	{
		PassDown(i);
		_x[i] = FROZEN_VALUE;
		_uhatTree[m][i] = _x[i];
		PassUp(i);

		double p0 = 1 - _beliefTree[m][i];
		double currentMetric = log(p0) - log(1 - _p[i]);
		// cumulative
		metrics[i] = currentMetric;
		if (i != 0)
			metrics[i] += metrics[i - 1];
	}

	int i = 0;
	int j = -1;

	std::size_t treePathCount = pow(2, k);

	std::size_t previous_number = treePathCount - 1;
	std::size_t current_j = 0;
	for (std::size_t number = 0; number < treePathCount; number++)
	{
		// 0 - if there is no distinctions, 1 - only last bit, therefore numerations started from 0
		current_j = GetFirstDistinctBit(previous_number, number);
		
		FillBinary(number, binary, current_j);

		std::cout << std::to_string(current_j) << std::endl;

		for (std::size_t l = 0; l < k; l++)
		{
			std::cout << std::to_string(binary[l]);
		}
		std::cout << std::endl;
		previous_number = number;

		current_j = k - current_j; // take into account numerations
		i = A[current_j];

		while (i < n)
		{
			PassDown(i); // get p1 metric in _beliefTree[m][i]

			double p = binary[current_j] ? _beliefTree[m][i] : 1 - _beliefTree[m][i];
			
			if (_maskWithCrc[i]) {

				double previous = (i == 0) ? 0 : metrics[i - 1];
				double currentMetric = previous + log(p / (1 - _p[i]));
				metrics[i] = currentMetric;

				_x[i] = binary[current_j];
				_uhatTree[m][i] = _x[i];
				PassUp(i);

				std::string binaryPath = "";
				for (std::size_t l = 0; l <= current_j; l++)
				{
					binaryPath += std::to_string(binary[l]);
				}
				std::string description = binaryPath + ", M=" + std::to_string(currentMetric) + ", " 
					+ std::to_string(p);

				PathNode newNode{binaryPath, description};
				unorderedTree.push_back(newNode);

				i++;
				current_j++;
			}
			else {
				_x[i] = FROZEN_VALUE;
				_uhatTree[m][i] = _x[i];
				PassUp(i);

				double p0 = 1 - _beliefTree[m][i];
				double currentMetric = log(p0) - log(1 - _p[i]);
				// cumulative
				metrics[i] = currentMetric;
				if (i != 0)
					metrics[i] += metrics[i - 1];

				i++;
			}

			// HERE trace
			// _path += std::to_string(i) + ": " + std::to_string(metrics[i]) + "\n";
			///////
		}
	}

	sort(unorderedTree.begin(), unorderedTree.end(), CompareTree);
	for (std::size_t i = 0; i < unorderedTree.size(); i++)
	{
		_pathTrace += unorderedTree[i].description + "\n";
	}
	
	return TakeResult();
}
