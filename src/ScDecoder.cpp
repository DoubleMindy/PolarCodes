#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <string>
#include <iostream>
#include <queue>
#include <algorithm>

#include "../include/PolarCode.h"
#include "../include/ScDecoder.h"
#include "../include/Exceptions.h"
#include "../include/Domain.h"

#define DBL_MAX 1.7976931348623158e+308 
#define FROZEN_VALUE 0

ScDecoder::ScDecoder(PolarCode * codePtr) : BaseDecoder(codePtr) {
	std::size_t m = _codePtr->m();
	std::size_t n = _codePtr->N();
	_treeHeight = m + 1;

	std::vector<double> b(n, -10000.0);
	std::vector<int> u(n, -1);

	for (std::size_t i = 0; i < _treeHeight; i++)
	{
		_beliefTree.push_back(b);
		_uhatTree.push_back(u);
	}

	_maskWithCrc = std::vector<int>(n, 0);
	auto maskInf = _codePtr->BitsMask();

	if (_codePtr->IsCrcUsed()) {
		auto maskCrc = _codePtr->CrcMask();
		for (std::size_t i = 0; i < n; i++)
			_maskWithCrc[i] = maskInf[i] || maskCrc[i];
	}
	else
		_maskWithCrc = maskInf;

	_x = std::vector<int>(n, -1);

	// optimization allocation
	_binaryIter = std::vector<int>(m, 0);

}

ScDecoder::~ScDecoder() {
}

void ScDecoder::FillLeftMessageInTree(std::vector<double>::iterator leftIt,
	std::vector<double>::iterator rightIt,
	std::vector<double>::iterator outIt,
	std::size_t n)
{
	for (std::size_t i = 0; i < n; i++, leftIt++, rightIt++, outIt++)
	{
		*outIt = f(*leftIt, *rightIt);
	}
}

void ScDecoder::FillRightMessageInTree(std::vector<double>::iterator leftIt,
	std::vector<double>::iterator rightIt,
	std::vector<int>::iterator uhatIt,
	std::vector<double>::iterator outIt,
	std::size_t n)
{
	for (std::size_t i = 0; i < n; i++, leftIt++, rightIt++, outIt++, uhatIt++)
	{
		*outIt = g(*leftIt, *rightIt, *uhatIt);
	}
}


std::size_t ScDecoder::FirstBitPos(std::size_t n) {
	std::size_t m = 0;
	while (n > 0)
	{
		n = n >> 1;
		m++;
	}

	return m;
}

void ScDecoder::PassDown(std::size_t iter) {
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

		if (!_binaryIter[i - level]) {
			FillLeftMessageInTree(_beliefTree[i].begin() + offset,
				_beliefTree[i].begin() + offset + length,
				_beliefTree[i + 1].begin() + offset,
				length);
		}
		else {

			FillRightMessageInTree(_beliefTree[i].begin() + offset,
				_beliefTree[i].begin() + offset + length,
				_uhatTree[i + 1].begin() + offset,
				_beliefTree[i + 1].begin() + offset + length,
				length);
		}

		length = length / 2;
	}
}

void ScDecoder::PassUp(std::size_t iter) {
	std::size_t m = _codePtr->m();
	std::size_t iterCopy = iter;

	std::size_t bit = iterCopy % 2;
	std::size_t length = 1;
	std::size_t level = m;
	std::size_t offset = iter;
	while (bit != 0)
	{
		offset -= length;
		for (std::size_t i = 0; i < length; i++)
		{
			_uhatTree[level - 1][offset + i] = _uhatTree[level][offset + i] ^ _uhatTree[level][offset + length + i];
		}
		for (std::size_t i = 0; i < length; i++)
		{
			_uhatTree[level - 1][offset + length + i] = _uhatTree[level][offset + length + i];
		}

		iterCopy = iterCopy >> 1;
		bit = iterCopy % 2;
		length *= 2;
		level -= 1;
	}
}

void ScDecoder::DecodeInternal(std::vector<double> inLlr) {
	std::size_t n = inLlr.size();
	std::size_t m = _codePtr->m();

	for (std::size_t i = 0; i < n; i++)
	{
		_beliefTree[0][i] = inLlr[i];
	}
	for (std::size_t i = 0; i < n; i++)
	{
		PassDown(i);
		if (_maskWithCrc[i]) {
			_x[i] = L(_beliefTree[m][i]);
		}
		else {
			_x[i] = FROZEN_VALUE;
		}
		_uhatTree[m][i] = _x[i];
		PassUp(i);
	}

	return;
}

std::vector<int> ScDecoder::TakeResult() {
	std::vector<int> result(_codePtr->k(), 0);
	std::vector<int> codewordBits = _codePtr->UnfrozenBits();
	for (std::size_t i = 0; i < codewordBits.size(); i++)
	{
		result[i] = _x[codewordBits[i]];
	}

	return result;
}

std::vector<int> ScDecoder::Decode(std::vector<double> inLlr) {
	
	DecodeInternal(inLlr);

	return TakeResult();
}