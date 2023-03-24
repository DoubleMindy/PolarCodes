#pragma once
#include <vector>
class PolarCode {

protected:
	
	std::size_t _m = 0;
	std::size_t _N = 0;
	std::size_t _k = 0;
	// with length of CRC hash
	std::size_t _k_extended = 0;

	std::vector<int> _bitsMask;
	std::vector<int> _crcMask;
	std::vector<int> _maskWithCrc;

	std::vector<int> _unfrozenPolarSequence;
	std::vector<int> _crcUnfrozenPolarSequence;
	std::vector<int> _unfrozenPolarSequenceWithCrc;

	std::vector<int> _unfrozenBits;
	std::vector<int> _unfrozenBitsWithCrc;
	std::vector<int> _crcUnfrozenBits;

	std::vector<int> _crcPoly;
	std::size_t _crcDeg;

public:
	PolarCode();
	PolarCode(int m, int k, std::vector<int> reliabilitySequence, std::vector<int> crcPoly);
	// Build code from k indeices of usedBits
	PolarCode(int m, std::vector<int> usedBits);

	void InitWithNewSequence(int m, std::vector<int> usedBits);

	std::size_t m();
	std::size_t N();
	std::size_t k();
	std::size_t kExt();
	
	std::vector<int> BitsMask();
	std::vector<int> CrcMask();
	std::vector<int> BitsMaskWithCrc();

	std::vector<int> UnfrozenPolarSequence();
	std::vector<int> UnfrozenPolarSequenceWithCrc();

	std::vector<int> UnfrozenBits();
	std::vector<int> UnfrozenBitsWithCrc();
	std::vector<int> CrcUnfrozenBits();

	std::vector<int> CrcPoly();
	std::size_t CrcDeg();

	bool IsCrcUsed();

	~PolarCode() {};
};