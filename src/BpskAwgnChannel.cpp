
#include <random>
#include "../include/BpskAwgnChannel.h"
#include "../include/CommonTransformations.h"

BpskAwgnChannel::BpskAwgnChannel() : BaseChannel() {
	_snr = 1.0; // trash value
	_sigma = snrToSigma(_snr);

	std::normal_distribution<double> normal_dist(0, _sigma);
	_normal_dist = normal_dist;
}

void BpskAwgnChannel::SetSnr(double snr) {
	BaseChannel::SetSnr(snr);
	_sigma = snrToSigma(snr);
	std::normal_distribution<double> normal_dist(0, _sigma);
	_normal_dist = normal_dist;
}

double InputToLlr(double input, double sigma) {
	return 2 * input / (sigma * sigma);
}

int ModulateBpsk(int input) {
	return 1 - 2 * input;
}

std::vector<double> BpskAwgnChannel::Pass(std::vector<int> codeword) {

	std::size_t n = codeword.size();
	std::vector<double> output(n, 0.0);
	for (std::size_t i = 0; i < n; i++) {
		output[i] = ModulateBpsk(codeword[i]);
		output[i] += _normal_dist(_randomDevice); // return LLRs
		output[i] = InputToLlr(output[i], _sigma);
	}

#ifdef DOMAIN_P1

	for (std::size_t i = 0; i < n; i++)
	{
		output[i] = LlrToP1(output[i]);
	}
	
#endif // DOMAIN_P1

	return output;

}