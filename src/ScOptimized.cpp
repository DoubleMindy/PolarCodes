
#include "../include/PolarCode.h"
#include "../include/ScOptimized.h"

#define FROZEN_VALUE 0

ScOptimizedDecoder::ScOptimizedDecoder(PolarCode * codePtr) : BaseDecoder(codePtr) {
	_m = _codePtr->m();
	_n = _codePtr->N();
	_k = _codePtr->k();

	vector<vector<int>> _beta_temp;

	for (std::size_t i = 0; i <= _m; i++)
	{
		_alpha.push_back(vector<double>(1 << i, 0.0));
		_beta_temp.push_back(vector<int>(1 << i, 0.0));
	}

	_beta.push_back(_beta_temp);
	_beta.push_back(_beta_temp);

	_mask = _codePtr->BitsMaskWithCrc();
}


#ifdef DOMAIN_LLR

double ScOptimizedDecoder::f(double left, double right) {
	double sign = 1;

	if (left < 0) {
		left = -left;
		sign = -sign;

	}
	if (right < 0) {
		right = -right;
		sign = -sign;

	}


	return ((left < right) ? left : right) * sign;
}

double ScOptimizedDecoder::g(double left, double right, int left_hard) {


	return right + (left_hard == 1 ? -1 : 1) * left;
}

int ScOptimizedDecoder::HD(double llr) {

	return llr < 0;
}

#elif DOMAIN_P1

double ScOptimizedDecoder::f(double p1, double p2) {

	return p1 * (1 - p2) + p2 * (1 - p1);
}

double ScOptimizedDecoder::g(double p1, double p2, int b) {
	double p1_b = f(b, p1);


	if (p1_b == 0 && p2 == 1 || p2 == 0 && p1_b == 1)
		return 1 / 2;

	return p1_b * p2 / (p1_b * p2 + (1 - p1_b) * (1 - p2));
}

int ScOptimizedDecoder::HD(double p1) {

	return p1 >= 0.5;
}
#endif

void ScOptimizedDecoder::recursively_calc_alpha(std::size_t lambda, std::size_t phi) {

	if (lambda == _m)
		return;
	

	std::size_t lambda_big = 1 << lambda;
	std::size_t lambda_next = lambda + 1;
	std::size_t isPhiOdd = phi % 2;

	if (isPhiOdd) {
		for (std::size_t i = 0; i < lambda_big; i++) {
			_alpha[lambda][i] = g(_alpha[lambda_next][i], _alpha[lambda_next][i + lambda_big], _beta[!isPhiOdd][lambda][i]);
		}
	}
	else {
		recursively_calc_alpha(lambda_next, phi >> 1);

		for (std::size_t i = 0; i < lambda_big; i++) {
			_alpha[lambda][i] = f(_alpha[lambda_next][i], _alpha[lambda_next][i + lambda_big]);

		}
	}
		
	return;
}

void ScOptimizedDecoder::recursively_calc_beta(std::size_t lambda, std::size_t phi) {


	std::size_t isPhiOdd = phi % 2;
	std::size_t isPhiEven = !isPhiOdd;
	if (isPhiEven)
		return;


	std::size_t lambda_big = 1 << lambda;
	std::size_t lambda_next = lambda + 1;
	std::size_t isHalfPhiOdd = (phi >> 1) % 2;

	for (std::size_t i = 0; i < lambda_big; i++)
	{
		_beta[isHalfPhiOdd][lambda_next][i] = _beta[isPhiEven][lambda][i] ^ _beta[isPhiOdd][lambda][i];
		_beta[isHalfPhiOdd][lambda_next][lambda_big + i] = _beta[isPhiOdd][lambda][i];

	}

	recursively_calc_beta(lambda_next, phi >> 1);
}

std::vector<int> ScOptimizedDecoder::Decode(std::vector<double> llr) {
	std::vector<int> result(_codePtr->k(), 0);
	std::size_t i = 0;
	int iterationsCount = 0;
	for (std::size_t phi = 0; phi < _n; phi++)
	{
		_alpha[_m][phi] = llr[phi];

	}

	for (std::size_t phi = 0; phi < _n; phi++)
	{
		recursively_calc_alpha(0, phi);

		if (_mask[phi]) {
			int hard_desicion = HD(_alpha[0][0]);

			_beta[phi % 2][0][0] = hard_desicion;
			result[i++] = hard_desicion;

		}
		else {
			_beta[phi % 2][0][0] = FROZEN_VALUE;

		}

		recursively_calc_beta(0, phi);
		iterationsCount++;
	}

	_operationsCount.Iterations += (double)iterationsCount / _n;
	_operationsCount.Normilizer++;

	return result;
}