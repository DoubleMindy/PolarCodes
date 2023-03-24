#pragma once
#include <cstddef>

#define PI_HALF_SQRT 1.253314137315500251207882642405522626503493370304969158314

#define PHI_MAX exp(0.0218)
#define PHI_GAP_MAX exp(-0.4527*pow(10, 0.86) + 0.0218)
#define PHI_GAP_MIN PI_HALF_SQRT * sqrt(1.0 / 5.0) * (6.0 / 7.0 ) * exp(- 2.5)
#define PHI_MIN 0

class GaussianApproximation
{

public:
	double _sigma = 0;
	double phi(double x);
	double phiDerivative(double x);
	double NumericalSolve(double y);
	double phiInv(double y);

	double GetMu(std::size_t i, std::size_t n);
public:
	GaussianApproximation(double sigma);
	double GetChannelErrorProbability(std::size_t j, std::size_t n);
};
