#include "../include/CommonTransformations.h"

#include <cmath>
#include <array>

double snrToSigma(double snr) {
	return sqrt(pow(10, -snr / 10) / 2);
}

double snrToEbN0(double snr, double coderate) {
	return snr - 10 * log10(coderate);
}

double ebnoToPErr(double sigma) {
	//return 1.0 / 2 * erfc(sqrt(ebno));
	return 1.0 / 2 * (1 + erf(-1.0 / sigma / sqrt(2)));
}

double LlrToP1(double llr) {
	if (llr > 300.0)
		return 0.0;

	if (llr < -300.0)
		return 1.0;

	return 1.0 / (1.0 + exp(llr));
}

std::array<double,2> LlrToLogP(double llr) {
	double denom = 0;
	if (llr > 300.0)
		denom = llr;
	else if (llr > -300.0)
		denom = std::log1p(exp(llr));
	return {llr-denom, -denom};
// function LogP = llr2logp(llr)
// if llr>308
//     denom=llr;
// elseif llr<-308
//     denom=0;
// else
//     denom=log1p(exp(llr));
// end
// logp0 = llr-denom;
// logp1 = -denom;
// LogP=[logp0 logp1];
// end
}
