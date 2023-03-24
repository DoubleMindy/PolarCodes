
#include "../include/CRC.h"

CRC::CRC(std::vector<int> poly) {
	_poly = poly;
	_deg = poly.size();
	_poly.push_back(1); // add main degree
	_paddingSymbol = 0;
	_init = std::vector<int>(_deg, 0);
	// here exception if poly is empty
}

std::vector<int> CRC::Calculate(std::vector<int> bits) {
	std::size_t crcLength = _deg;
	std::size_t inputLength = bits.size();

	std::vector<int> crc = _init;

	std::size_t blocksCount = inputLength / crcLength;
	if (inputLength % crcLength) {
		blocksCount++;
	}
	
	std::vector<int> inputPadded(blocksCount * crcLength, _paddingSymbol);

	for (std::size_t i = 0; i < inputLength; i++)
	{
		inputPadded[i] = bits[i];
	}

	for (std::size_t i = 0; i < blocksCount; i++)
	{
		std::size_t offset = i * crcLength;


		for (std::size_t j = 0; j < crcLength; j++)
		{
			crc[j] ^= inputPadded[offset + j];

		}
		
		for (std::size_t j = 0; j < crcLength; j++)
		{
			int flag = crc[crcLength - 1] == 1;
			
			// shift
			for (std::size_t k = crcLength - 1; k > 0; k--)
			{
				crc[k] = crc[k - 1];
			}
			crc[0] = 0;

			// division
			if (flag)
				for (std::size_t k = 0; k < crcLength; k++)
				{
					crc[k] ^= _poly[k];
				}
		}
	}
	
	return crc;
}
