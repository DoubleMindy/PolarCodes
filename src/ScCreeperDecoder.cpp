#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>

#include "../include/ScCreeperDecoder.h"
#include "../include/GaussianApproximation.h"
#include "../include/CommonTransformations.h"

using std::vector;
using std::make_pair;
using std::max;
using std::cout;

#define LOW_INFINITY -1000000000.0
#define DUMMY_LENGTH -1

ScCreeperDecoder::ScCreeperDecoder(PolarCode * codePtr, double delta) : ScOptimizedDecoder(codePtr) {
	_delta = delta;

	_metric = vector<double>(_n, 0.0);
	_path = vector<int>(_n, 0);
}

void ScCreeperDecoder::recursively_calc_alpha_creeper(std::size_t lambda, std::size_t phi, bool isPathSwitched) {

	if (lambda == _m)
		return;

	std::size_t lambda_big = 1 << lambda;
	std::size_t lambda_next = lambda + 1;
	std::size_t isPhiOdd = phi % 2;

	if (!isPhiOdd || isPathSwitched)
		recursively_calc_alpha_creeper(lambda_next, phi >> 1, isPathSwitched);

	if (isPhiOdd) {
		for (std::size_t i = 0; i < lambda_big; i++) {
			_alpha[lambda][i] = g(_alpha[lambda_next][i], _alpha[lambda_next][i + lambda_big], _beta[!isPhiOdd][lambda][i]);
			//cout << _alpha[lambda][i] << " ";
		}
	}
	else {
		for (std::size_t i = 0; i < lambda_big; i++) {
			_alpha[lambda][i] = f(_alpha[lambda_next][i], _alpha[lambda_next][i + lambda_big]);
			//cout << _alpha[lambda][i] << " ";
		}
	}
	return;
}

double ScCreeperDecoder::calculate_step_metric_fano(double belief, int decision, double pe) {
#ifdef DOMAIN_LLR
	// double p1 = LlrToP1(belief);
	// double p0 = 1 - p1;
	// double abs = fabs(belief);
	//return (belief < 0 && decision == 0 || belief > 0 && decision == 1) ? -abs: 0.0;
	// return log(((decision) ? p1 : p0)) - log(1 - pe);
	auto logp = LlrToLogP(belief);
	return logp[static_cast<std::size_t>(decision)] - log(1 - pe);
#elif DOMAIN_P1
	double p0 = 1 - belief;
	double p1 = belief;
	return log(((decision) ? p1 : p0)) - log(1 - pe);
#endif
}

double ScCreeperDecoder::Q(double x, double delta) {
	//return x;
	return floor(x / delta) * delta;
}

void ScCreeperDecoder::LoadPath(int length, int lastBit) {
	std::size_t isPhiOdd;

	_path[length - 1] = lastBit;
	for (std::size_t phi = 0; phi < length; phi++)
	{
		isPhiOdd = phi % 2;

		_beta[isPhiOdd][0][0] = _path[phi];
		if (isPhiOdd)
			recursively_calc_beta(0, phi);
	}
}

std::vector<int> ScCreeperDecoder::Decode(std::vector<double> belief) {
	// local variables
	int next_bit_pos = 0;
	bool isPathSwitched = false;
	double m = 0.0;
	double m_max = -1000000000.0;
	double pm0 = 0.0;
	double pm1 = 0.0;
	double mx = 0.0;
	double my = 0.0;
	double T = 0.0;
	double max_temp = 0.0;
	int phi = 0;
	int next_bit = 0;
	std::pair<int, int> n_current = { 0, 0 };
	int iterationsCount = 0;
	vector<vector<int>> haveNotPassedCRC;

	vector<int> result(_k, 0);
	double approximationSigma = sqrt(pow(10, -1.0 / 10.0));
	GaussianApproximation ga(approximationSigma);
	std::vector<double> p = std::vector<double>(_n, 0);
	for (std::size_t i = 0; i < _n; i++)
	{
		p[i] = ga.GetChannelErrorProbability(i + 1, _n);
		_metric[i] = 0.0;
	}

	_NP.clear();
	_F.clear();
	_TP.clear();

	// alphas init
	for (std::size_t phi = 0; phi < _n; phi++)
	{
		_alpha[_m][phi] = belief[phi];
	}

	// init stack for the first step will be Rule One
	// dummy node
	_NP.push_front(make_pair(DUMMY_LENGTH, 1));
	_TP.push_front(LOW_INFINITY);
	_F.push_front(true);

	// root node
	std::pair<int, int> root = make_pair(DUMMY_LENGTH, 0);
	_NP.push_front(root);
	_TP.push_front(0.0);
	_F.push_front(true);

	// initial state of the variables
	T = LOW_INFINITY;
	m_max = LOW_INFINITY;
	n_current = root;

	std::string debug = "";
	int checksCount = 0;

	bool isFifthStep = false;// debug
	while (true) {
		iterationsCount++;
		isPathSwitched = false;
		next_bit_pos = n_current.first + 1;
		recursively_calc_alpha_creeper(0, next_bit_pos, isPathSwitched);
		pm0 = ((next_bit_pos != 0) ? _metric[next_bit_pos - 1] : 0)
			+ calculate_step_metric_fano(_alpha[0][0], 0, p[next_bit_pos]);

		if (_mask[next_bit_pos]) {
			T = Q(_TP[1], _delta);
			
			pm1 = ((next_bit_pos != 0) ? _metric[next_bit_pos - 1] : 0)
				+ calculate_step_metric_fano(_alpha[0][0], 1, p[next_bit_pos]);

			if (pm0 > pm1) {
				mx = pm0;
				my = pm1;
				next_bit = 0;
			}
			else {
				mx = pm1;
				my = pm0;
				next_bit = 1;
			}

			// Perform one of six rules
			if (T <= my) {
				// First Rule
				if (mx > m_max) {
					m_max = mx;
					_TP.push_front(my);
					_TP.push_front(LOW_INFINITY);

					_F.push_front(true);
					_F.push_front(true);
					debug += "1";
				}
				// Second Rule
				else {
					_F.push_front(false);
					_F.push_front(false);
					debug += "2";
				}

				_metric[next_bit_pos] = mx;
				n_current = { next_bit_pos, next_bit };
				_NP.push_front({ next_bit_pos, !next_bit }); // ny
				_NP.push_front({ next_bit_pos, next_bit }); // nx
			}
			else {
				// Third Rule
				if (mx >= T) {
					n_current = { next_bit_pos, next_bit };
					_TP[0] = max(my, _TP[0]);
					m_max = max(mx, m_max);
					debug += "3";
					_metric[next_bit_pos] = mx;
				}
				else {
					max_temp = max(mx, _TP[0]);
					n_current = _NP[1];
					isPathSwitched = true;

					if (_F[1]) {
						// Forth Rule
						if (max_temp >= Q(_TP[3], _delta)) {
							_TP[0] = max_temp;

							std::swap(_TP[0], _TP[1]);
							std::swap(_F[0], _F[1]);
							std::swap(_NP[0], _NP[1]);
							_TP[0] = LOW_INFINITY;
							debug += "4";
						}
						// Fifth Rule
						else {
							_TP[2] = max(max_temp, _TP[2]);
							_NP.pop_front();
							_NP.pop_front();
							_TP.pop_front();
							_TP.pop_front();
							_F.pop_front();
							_F.pop_front();
							debug += "5";
							isFifthStep = true;
						}
					}
					// Sixth Rule
					else {
						_TP[0] = max_temp;
						_NP.pop_front();
						_NP.pop_front();
						_F.pop_front();
						_F.pop_front();

						debug += "6";
					}
				}
			}
		}
		else {
			n_current = { next_bit_pos, 0 };
			_metric[next_bit_pos] = pm0;
		}
		
		if (isPathSwitched) {
			LoadPath(n_current.first + 1, n_current.second); // argument is length
			recursively_calc_alpha_creeper(0, n_current.first, isPathSwitched); // restore alpha for metric
			_metric[n_current.first] = ((n_current.first != 0) ? _metric[n_current.first - 1] : 0)
				+ calculate_step_metric_fano(_alpha[0][0], n_current.second, p[n_current.first]);
			isPathSwitched = false;
			continue;
		}

		phi = n_current.first;
		
		_path[phi] = _beta[phi % 2][0][0] = n_current.second;

		if (phi == _n - 1) {
			// here CRC check
			if (isFifthStep)
				std::cout << debug << "\n";
			break;
		}

		if (phi % 2 == 1)
			recursively_calc_beta(0, phi);
	}

	std::vector<int> codewordBits = _codePtr->UnfrozenBits();
	for (std::size_t i = 0; i < codewordBits.size(); i++)
	{
		result[i] = _path[codewordBits[i]];
	}

	_operationsCount.Iterations += (double)iterationsCount / _n;
	_operationsCount.Normilizer++;

	return result;
}
