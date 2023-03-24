#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

#include "../include/PolarCode.h"
#include "../include/ScDecoder.h"
#include "../include/ScFanoDecoder.h"
#include "../include/Exceptions.h"
#include "../include/Domain.h"
#include "../include/GaussianApproximation.h"

#define DBL_MAX 1.7976931348623158e+308 
#define FROZEN_VALUE 0

ScFanoDecoder::ScFanoDecoder(PolarCode * codePtr, double T, double delta) : ScDecoder(codePtr) {
	std::size_t n = _codePtr->N();

	_T = T;
	_delta = delta;

	double approximationSigma = sqrt(pow(10, -1.0 / 10.0));
	GaussianApproximation ga(approximationSigma);
	_p = std::vector<double>(n, 0);
	for (std::size_t i = 0; i < n; i++)
	{
		_p[i] = ga.GetChannelErrorProbability(i + 1, n);
	}
}

void UpdateT(double & T, double & delta, double & tau) {
	while (T + delta < tau)
		T += delta;
}

std::string ScFanoDecoder::GetPathInfo() {
	return _pathTrace;
}

void ScFanoDecoder::BackwardMove(std::vector<double> & beta, std::vector<bool> & gamma, std::vector<int> & mask, int & i, double & T, double & delta, bool & B, int & j) {

	while (true) {
		double mu = 0;

		if (j <= -1)
			mu = -1000;

		if (j >= 1)
			mu = beta[j - 1];

		if (mu >= T) {
			j--;
			// HERE trace
			_pathTrace += std::to_string(j) + ": " + std::to_string((j > -1) ? beta[j] : 0.0) + ", back, T=" + std::to_string(T) + "\n";
			///////
			if (!gamma[j + 1])
			{
				B = true;
				return;
			}
		}
		else {
			T -= delta;
			B = false;
			return;
		}
			
	}
}

std::vector<int> ScFanoDecoder::Decode(std::vector<double> inP1) {
	_pathTrace = "";

	std::size_t n = inP1.size();
	std::size_t m = _codePtr->m();
	std::size_t k = _codePtr->k();
	for (std::size_t i = 0; i < n; i++)
	{
		_beliefTree[0][i] = inP1[i];
	}

	int i = 0;
	int j = -1;
	bool B = false;
	std::vector<double> beta(k, 0.0); // only for frozen bits
	std::vector<double> metrics(n, 0); // for all bits
	std::vector<bool> gamma(k, 0);
	std::vector<int> A = _codePtr->UnfrozenBits(); // info set
	double T = _T;
	double delta = _delta;
	int iterationsCount = 0;
	//int checksCount = 0;

	std::size_t firstInfoBit = 0;
	for (std::size_t i = 0; i < n; i++)
	{
		if (_maskWithCrc[i]) {
			firstInfoBit = i;
			break;
		}
	}
	while (i < n)
	{
		iterationsCount++;
		PassDown(i); // get p1 metric in _beliefTree[m][i]
		double p0 = 1 - _beliefTree[m][i];
		double p1 = _beliefTree[m][i];
		
		if (_maskWithCrc[i]) {
			double previous = (i == 0) ? 0 : metrics[i - 1];
			double m0 = previous + log(p0 / (1 - _p[i]));
			double m1 = previous + log(p1 / (1 - _p[i]));

			double max = (m1 > m0) ? m1 : m0 ;
			int argmax = (m1 > m0) ? 1 : 0;

			double min = (m1 > m0) ? m0 : m1;
			int argmin = (m1 > m0) ? 0 : 1;

			/*if ((i == n - 1) && (max > T)) {
				if (!B) {
					_x[i] = argmax;
					if (IsCrcPassed(_x))
						break;
					else {
						max = -1000000.0;
					}
				}
				else {
					_x[i] = argmin;
					if (IsCrcPassed(_x))
						break;
					else {
						min = -1000000.0;
					}
				}
			}*/

			if (max > T) {
				if (!B) {
					_x[i] = argmax;
					_uhatTree[m][i] = _x[i];
					PassUp(i);

					metrics[i] = max;
					beta[j + 1] = max;
					gamma[j + 1] = false;

					double mu = 0;
					if (j != -1) 
						mu = beta[j];

					if (mu < T + delta) 
						UpdateT(T, delta, beta[j+1]);
					i++;
					j++;
					// HERE trace
					_pathTrace += std::to_string(j) + ": " + std::to_string(beta[j]) + ", right, T=" + std::to_string(T) 
						+ ", x[i]=" + std::to_string(argmax) + ", p=" + std::to_string(argmax == 1 ? p1 : p0) + "\n";
					///////
				}
				else {
					if (min > T) {
						_x[i] = argmin;
						_uhatTree[m][i] = _x[i];
						PassUp(i);

						metrics[i] = min;
						beta[j + 1] = min;
						gamma[j + 1] = true;
						B = false;

						i++;
						j++;

						// HERE trace
						_pathTrace += std::to_string(j) + ": " + std::to_string(beta[j]) + ", left, T=" + std::to_string(T) 
							+ ", x[i]=" + std::to_string(argmin) + ", p=" + std::to_string(argmin == 1 ? p1 : p0) + "\n";
						///////
					}
					else {
						if (j == -1) {
							T = T - delta;
							B = false;

							// HERE trace
							_pathTrace += std::to_string(j) + ": " + "0.0" + ", lower, T=" + std::to_string(T) + "\n";
							///////
						}
						else {
							BackwardMove(beta, gamma, _maskWithCrc, i, T, delta, B, j);
							i = A[j + 1];
						}
					}
				}
			}
			else {
				if (j == -1) {
					T = T - delta;
					// HERE trace
					_pathTrace += std::to_string(j) + ": " + "0.0" + ", lower, T=" + std::to_string(T) + "\n";
					///////
				}
				else {
					BackwardMove(beta, gamma, _maskWithCrc, i, T, delta, B, j);
					i = A[j + 1];
				}
					
			}
		}
		else {
			_x[i] = FROZEN_VALUE;
			_uhatTree[m][i] = _x[i];
			PassUp(i);

			double currentMetric = log(p0) - log(1 - _p[i]);
			// cumulative
			metrics[i] = currentMetric;
			if (i != 0)
				metrics[i] += metrics[i - 1];

			i++;
		}

		// HERE trace
		// _path += std::to_string(i) + ": " + std::to_string(metrics[i]) + "\n";
		///////
	}

	_operationsCount.Iterations += (double)iterationsCount / n;
	_operationsCount.Normilizer++;

	return TakeResult();
}