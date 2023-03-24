#include <math.h>

#include "../include/PolarCode.h"
#include "../include/Exceptions.h"
#include "../include/Encoder.h"
#include "../include/CRC.h"

Encoder::Encoder(PolarCode * codePtr) {
	_codePtr = codePtr;
	_crcPtr = new CRC(_codePtr->CrcPoly());
}

Encoder::~Encoder() {
	delete _crcPtr;
}

void PolarTransformRecursive(std::vector<int>::iterator begin, std::size_t n) {
	
	if (n == 1)
		return;

	std::size_t half_n = n / 2;
	for (auto it = begin; it < begin + half_n; it++)
		*it = (*it + *(it + half_n)) % 2;
	
	PolarTransformRecursive(begin, half_n);
	PolarTransformRecursive(begin + half_n, half_n);

	return;
}

std::vector<int> Encoder::PolarTransform(std::vector<int> x) {
	PolarTransformRecursive(x.begin(), x.size());
	return x;
}

std::vector<int> Encoder::Encode(std::vector<int> word) {
	std::size_t k = word.size();

	if (k != _codePtr->k())
		throw IncorrectVectorSizeException("Length of the word does not matched with code parameter k");

	std::size_t n = _codePtr->N();
	std::vector<int> codeword(n, 0);
	std::vector<int> codewordBits = _codePtr->UnfrozenBits();

	for (std::size_t i = 0; i < codewordBits.size(); i++)
	{ 
		codeword[codewordBits[i]] = word[i];
	}

	if (_codePtr->IsCrcUsed()) {
		std::vector<int> crc = _crcPtr->Calculate(word);
		std::vector<int> crcBits = _codePtr->CrcUnfrozenBits();

		for (std::size_t i = 0; i < crcBits.size(); i++)
		{
			codeword[crcBits[i]] = crc[i];
		}
	}

	PolarTransformRecursive(codeword.begin(), n);
	return codeword;
}