#pragma once

#include "BaseChannel.h"
#include <vector>
#include <random>

class BpskBscChannel : public BaseChannel {
protected:
	std::random_device _randomDevice;
	std::bernoulli_distribution _bernoulli_dist;
	double _coderate;
	double _fixedLlr;

public:
	BpskBscChannel();
	~BpskBscChannel() {};
	std::vector<double> Pass(std::vector<int> codeword) override;
	void SetSnr(double snr) override;
};