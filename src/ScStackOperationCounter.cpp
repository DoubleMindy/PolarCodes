#include <cmath>
#include "../include/PolarCode.h"
#include "../include/ScStackOperationCounter.h"

#define FROZEN_VALUE 0
#define MAX_METRIC 1000000.0

ScStackOperationCounter::ScStackOperationCounter(PolarCode * codePtr, int L, int D) : ScOptimizedDecoder(codePtr) {
	_L = L;
	_D = D;

	_kWithCrc = _codePtr->kExt();
	_crcPtr = new CRC(_codePtr->CrcPoly());

	_path_metrics = vector<double>(_D, 0.0);
	_path_lengths = vector<std::size_t>(_D, 0);
	_path_lengths_hits = vector<std::size_t>(_n + 1, 0);
	//_inactive_path_indices = vector<std::size_t>(_D, 0);
	_is_paths_active = vector<bool>(_D, false);

	vector<int> temp = vector<int>(_n, 0);
	for (std::size_t i = 0; i < _D; i++)
	{
		_paths.push_back(temp);
	}
}

double ScStackOperationCounter::f(double left, double right) {
	double sign = 1;

	_operationsCount.Comps += 3;
	_operationsCount.Muls += 1;

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

double ScStackOperationCounter::g(double left, double right, int left_hard) {

	_operationsCount.Comps += 1;
	_operationsCount.Muls += 1;
	_operationsCount.Sums += 1;

	return right + (left_hard == 1 ? -1 : 1) * left;
}

int ScStackOperationCounter::HD(double llr) {

	_operationsCount.Comps++;

	return llr < 0;
}

void ScStackOperationCounter::recursively_calc_alpha_stack(std::size_t lambda, std::size_t phi, bool isPathSwitched) {

	if (lambda == _m) {
		return;
	}
	_operationsCount.Comps++;

	std::size_t lambda_big = 1 << lambda;
	std::size_t lambda_next = lambda + 1;
	std::size_t isPhiOdd = phi % 2;

	_operationsCount.Sums++;
	_operationsCount.Comps+=2;

	if (!isPhiOdd || isPathSwitched) {
		recursively_calc_alpha_stack(lambda_next, phi >> 1, isPathSwitched);
	}

	if (isPhiOdd) {
		for (std::size_t i = 0; i < lambda_big; i++) {
			_alpha[lambda][i] = g(_alpha[lambda_next][i], _alpha[lambda_next][i + lambda_big], _beta[!isPhiOdd][lambda][i]);
			_operationsCount.Comps++;
			_operationsCount.Sums += 2;
		}
	}
	else {
		for (std::size_t i = 0; i < lambda_big; i++) {
			_alpha[lambda][i] = f(_alpha[lambda_next][i], _alpha[lambda_next][i + lambda_big]);
			_operationsCount.Comps++;
			_operationsCount.Sums += 2;
		}
	}

	return;
}

void ScStackOperationCounter::recursively_calc_beta(std::size_t lambda, std::size_t phi) {

	std::size_t isPhiOdd = phi % 2;
	std::size_t isPhiEven = !isPhiOdd;
	_operationsCount.Comps++;
	if (isPhiEven)
		return;

	_operationsCount.Sums ++;

	std::size_t lambda_big = 1 << lambda;
	std::size_t lambda_next = lambda + 1;
	std::size_t isHalfPhiOdd = (phi >> 1) % 2;

	for (std::size_t i = 0; i < lambda_big; i++)
	{
		_beta[isHalfPhiOdd][lambda_next][i] = _beta[isPhiEven][lambda][i] ^ _beta[isPhiOdd][lambda][i];
		_beta[isHalfPhiOdd][lambda_next][lambda_big + i] = _beta[isPhiOdd][lambda][i];

		_operationsCount.Comps++;
		_operationsCount.Sums += 2;
		_operationsCount.Xors ++;
	}

	recursively_calc_beta(lambda_next, phi >> 1);
}


double ScStackOperationCounter::calculate_step_metric(double newLlr, int decision) {
	// economize operations

	_operationsCount.Comps++;

	if (newLlr < 0) {
		_operationsCount.Comps++;

		if (decision == 0) {
			return fabs(newLlr);
		}
		else {
			return 0;
		}
	}
	else {
		_operationsCount.Comps++;

		if (decision == 1) {
			return fabs(newLlr);
		}
		else {
			return 0;
		}
	}
}

void ScStackOperationCounter::KillPath(std::size_t index, std::size_t & T) {
	_is_paths_active[index] = false;
	_inactive_path_indices.push(index);
	T--;

	_operationsCount.Sums++;
}
std::size_t ScStackOperationCounter::ClonePath(std::size_t index, std::size_t & T) {
	std::size_t new_index = _inactive_path_indices.top();
	_inactive_path_indices.pop();

	_is_paths_active[new_index] = true;
	_path_metrics[new_index] = _path_metrics[index];
	_path_lengths[new_index] = _path_lengths[index];
	T++;
	_operationsCount.Sums++;

	for (std::size_t i = 0; i < _path_lengths[index]; i++)
	{
		_paths[new_index][i] = _paths[index][i];
		_operationsCount.Sums++;
		_operationsCount.Comps++;
	}

	return new_index;
}

void ScStackOperationCounter::ExtendPath(std::size_t index, int bit, double pathMetric) {
	_paths[index][_path_lengths[index]] = bit;

	_path_lengths[index]++;
	_path_metrics[index] = pathMetric;
	_operationsCount.Sums++;
}

void ScStackOperationCounter::LoadPath(std::size_t index) {
	std::size_t isPhiOdd;

	for (std::size_t phi = 0; phi < _path_lengths[index]; phi++)
	{
		isPhiOdd = phi % 2;

		_beta[isPhiOdd][0][0] = _paths[index][phi];
		if (isPhiOdd)
			recursively_calc_beta(0, phi);

		_operationsCount.Sums++;
		_operationsCount.Comps += 2;
	}
}

std::vector<int> ScStackOperationCounter::Decode(std::vector<double> llr) {
	vector<int> result = vector<int>(_k, 0);
	// variables
	std::size_t T = 0;
	std::size_t min_index = 0;
	std::size_t max_index = 0;
	bool isPathSwitched = false;
	double pm0 = 0.0;
	double pm1 = 0.0;
	std::size_t phi = 0;
	bool isPhiEven = false;
	std::size_t passedLength = 0;
	bool isTruncateNeeded = false;
	int iterationsCount = 0;

	// init structures
	_inactive_path_indices = stack<std::size_t>();
	for (std::size_t i = 0; i < _D; i++)
	{
		_inactive_path_indices.push(i);
		_is_paths_active[i] = false;
		_path_metrics[i] = 0.0;
		_path_lengths[i] = 0;

		_operationsCount.Sums++;
		_operationsCount.Comps++;
	}
	for (std::size_t i = 0; i < _n + 1; i++)
	{
		_path_lengths_hits[i] = 0;
		_operationsCount.Sums++;
		_operationsCount.Comps++;
	}

	// set initial path
	min_index = _inactive_path_indices.top();
	_inactive_path_indices.pop();
	_is_paths_active[min_index] = true;
	_path_metrics[min_index] = 0.0;
	_path_lengths[min_index] = 0;
	T = 1;

	// load input llrs
	for (std::size_t phi = 0; phi < _n; phi++)
	{
		_alpha[_m][phi] = llr[phi];

		_operationsCount.Sums++;
		_operationsCount.Comps++;
	}
	
	while (true) {

		iterationsCount++;
		recursively_calc_alpha_stack(0, _path_lengths[min_index], isPathSwitched);

		pm0 = _path_metrics[min_index] + calculate_step_metric(_alpha[0][0], 0);
		pm1 = _path_metrics[min_index] + calculate_step_metric(_alpha[0][0], 1);

		_operationsCount.Sums += 2;

		// for checking poped up path on the condition of the width search L
		passedLength = _path_lengths[min_index];
		isTruncateNeeded = false;

		_operationsCount.Comps++;

		if (passedLength < _n) {
			_operationsCount.Comps++;

			if (_mask[_path_lengths[min_index]]) {

				_operationsCount.Comps++;

				if (T == _D) {
					if (_path_metrics[max_index] > (pm0 > pm1 ? pm0 : pm1)) {
						KillPath(max_index, T);
					}

					_operationsCount.Comps += 2;
				}

				_operationsCount.Comps++;

				if (pm0 < pm1) {
					_operationsCount.Comps++;

					if (T < _D) {
						max_index = ClonePath(min_index, T);
						ExtendPath(max_index, 1, pm1);
					}

					ExtendPath(min_index, 0, pm0);
				}
				else {
					_operationsCount.Comps++;

					if (T < _D) {
						max_index = ClonePath(min_index, T);
						ExtendPath(max_index, 0, pm0);
					}

					ExtendPath(min_index, 1, pm1);
				}

			}
			else {
				ExtendPath(min_index, FROZEN_VALUE, pm0);
			}

			// Update length info
			_path_lengths_hits[passedLength]++;
			_operationsCount.Sums++;
			_operationsCount.Comps++;
			if (_path_lengths_hits[passedLength] >= _L)
				isTruncateNeeded = true;
		}
		// Crc check
		else {
			if (IsCrcPassed(_paths[min_index]))
				break;

			_operationsCount.Sums++;
			_operationsCount.Comps++;

			KillPath(min_index, T);
			_path_lengths_hits[_n]++;
			if (_path_lengths_hits[_n] >= _L)
				break;
		}

		// loop for new indices
		std::size_t new_min_index = 0;
		std::size_t new_max_index = 0;
		double min_metric = MAX_METRIC;
		double max_metric = -MAX_METRIC;

		for (std::size_t i = 0; i < _D; i++)
		{
			_operationsCount.Comps += 2;
			_operationsCount.Sums++;

			if (!_is_paths_active[i])
				continue;

			_operationsCount.Comps += 1;

			if (isTruncateNeeded) {

				_operationsCount.Comps++;

				if (_path_lengths[i] <= passedLength) {
					KillPath(i, T);
					continue;
				}
			}

			_operationsCount.Comps += 1;

			if (_path_metrics[i] < min_metric) {
				min_metric = _path_metrics[i];
				new_min_index = i;
			}

			_operationsCount.Comps += 1;

			if (_path_metrics[i] > max_metric) {
				max_metric = _path_metrics[i];
				new_max_index = i;
			}
		}

		isPathSwitched = min_index != new_min_index;
		min_index = new_min_index;
		max_index = new_max_index;

		_operationsCount.Comps += 4;

		// stack is empty
		if (_inactive_path_indices.size() == _D)
		{
			min_index = 0;
			for (std::size_t i = 0; i < _n; i++)
			{
				_paths[min_index][i] = 0;
			}
			break;
		}

		if (isPathSwitched)
			LoadPath(min_index);

		phi = _path_lengths[min_index] - 1;
		_beta[phi % 2][0][0] = _paths[min_index][phi];
		if (phi % 2 == 1)
			recursively_calc_beta(0, phi);
	}

	std::vector<int> codewordBits = _codePtr->UnfrozenBits();
	for (std::size_t i = 0; i < codewordBits.size(); i++)
	{
		result[i] = _paths[min_index][codewordBits[i]];

		_operationsCount.Sums++;
		_operationsCount.Comps++;
	}

	_operationsCount.Iterations += (double) iterationsCount / _n;
	_operationsCount.Normilizer++;

	return result;
}
