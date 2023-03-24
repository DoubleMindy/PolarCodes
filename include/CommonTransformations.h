#pragma once


double snrToSigma(double snr);

double snrToEbN0(double snr, double coderate);

double ebnoToPErr(double sigma);

double LlrToP1(double llr);

#include <array>
extern std::array<double,2> LlrToLogP(double llr);