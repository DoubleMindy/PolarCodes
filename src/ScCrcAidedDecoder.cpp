#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <string>
#include <iostream>
#include <queue>
#include <algorithm>

#include "../include/PolarCode.h"
#include "../include/ScCrcAidedDecoder.h"
#include "../include/Exceptions.h"
#include "../include/Domain.h"

#define FROZEN_VALUE 0

ScCrcAidedDecoder::ScCrcAidedDecoder(PolarCode * codePtr) : ScDecoder(codePtr) {
	
}

void ScCrcAidedDecoder::DecodeFrom(int position) {
	std::size_t n = _codePtr->N();
	std::size_t m = _codePtr->m();

	_uhatTree[m][position] = _x[position];
	PassUp(position);
	for (std::size_t i = position + 1; i < n; i++)
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
}

void ScCrcAidedDecoder::DecodeFromTo(int startPosition, int endPosition) {
	std::size_t n = _codePtr->N();
	std::size_t m = _codePtr->m();

	_uhatTree[m][startPosition] = _x[startPosition];
	PassUp(startPosition);
	for (std::size_t i = startPosition + 1; i < endPosition; i++)
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

	if (endPosition != n) {
		PassDown(endPosition);
		if (_maskWithCrc[endPosition]) {
			_x[endPosition] = L(_beliefTree[m][endPosition]);
		}
		else {
			_x[endPosition] = FROZEN_VALUE;
		}
	}
}
