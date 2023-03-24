#include <cmath>
#include "../include/PolarCode.h"
#include "../include/ScStackOptimizedDecoder.h"

#define FROZEN_VALUE 0
#define MAX_METRIC 1000000.0

ScStackOptimizedDecoder::ScStackOptimizedDecoder(PolarCode * codePtr, int L, int D) : ScOptimizedDecoder(codePtr) {
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

void ScStackOptimizedDecoder::recursively_calc_alpha_stack(std::size_t lambda, std::size_t phi, bool isPathSwitched) {

	if (lambda == _m)
		return;

	std::size_t lambda_big = 1 << lambda;
	std::size_t lambda_next = lambda + 1;
	std::size_t isPhiOdd = phi % 2;

	if (!isPhiOdd || isPathSwitched)
		recursively_calc_alpha_stack(lambda_next, phi >> 1, isPathSwitched);

	if (isPhiOdd) {
		for (std::size_t i = 0; i < lambda_big; i++) {
			_alpha[lambda][i] = g(_alpha[lambda_next][i], _alpha[lambda_next][i + lambda_big], _beta[!isPhiOdd][lambda][i]);
		}
	}
	else {
		for (std::size_t i = 0; i < lambda_big; i++) {
			_alpha[lambda][i] = f(_alpha[lambda_next][i], _alpha[lambda_next][i + lambda_big]);
		}
	}

	return;
}

double ScStackOptimizedDecoder::calculate_step_metric(double newLlr, int decision) {
	// economize operations
	return (newLlr < 0 && decision == 0 || newLlr > 0 && decision == 1) ? fabs(newLlr) : 0;
}

void ScStackOptimizedDecoder::KillPath(std::size_t index, std::size_t & T) {
	_is_paths_active[index] = false;
	_inactive_path_indices.push(index);
	T--;
}
std::size_t ScStackOptimizedDecoder::ClonePath(std::size_t index, std::size_t & T) {
	std::size_t new_index = _inactive_path_indices.top();
	_inactive_path_indices.pop();

	_is_paths_active[new_index] = true;
	_path_metrics[new_index] = _path_metrics[index];
	_path_lengths[new_index] = _path_lengths[index];
	T++;

	for (std::size_t i = 0; i < _path_lengths[index]; i++)
	{
		_paths[new_index][i] = _paths[index][i];
	}

	return new_index;
}

void ScStackOptimizedDecoder::ExtendPath(std::size_t index, int bit, double pathMetric) {
	_paths[index][_path_lengths[index]] = bit;
	
	_path_lengths[index]++;
	_path_metrics[index] = pathMetric;
}

void ScStackOptimizedDecoder::LoadPath(std::size_t index) {
	std::size_t isPhiOdd;

	for (std::size_t phi = 0; phi < _path_lengths[index]; phi++)
	{
		isPhiOdd = phi % 2;

		_beta[isPhiOdd][0][0] = _paths[index][phi];
		if (isPhiOdd)
			recursively_calc_beta(0, phi);
	}
}

std::vector<int> ScStackOptimizedDecoder::Decode(std::vector<double> llr) {
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

	// init structures
	_inactive_path_indices = stack<std::size_t>();
	for (std::size_t i = 0; i < _D; i++)
	{
		_inactive_path_indices.push(i);
		_is_paths_active[i] = false;
		_path_metrics[i] = 0.0;
		_path_lengths[i] = 0;
	}
	for (std::size_t i = 0; i < _n + 1; i++)
	{
		_path_lengths_hits[i] = 0;
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
	}
	int checks = 0;
	while (true) {
		
		recursively_calc_alpha_stack(0, _path_lengths[min_index], isPathSwitched);
				
		pm0 = _path_metrics[min_index] + calculate_step_metric(_alpha[0][0], 0);
		pm1 = _path_metrics[min_index] + calculate_step_metric(_alpha[0][0], 1);

		// for checking poped up path on the condition of the width search L
		passedLength = _path_lengths[min_index];
		isTruncateNeeded = false;

		if (passedLength < _n) {
			if (_mask[_path_lengths[min_index]]) {
				if (T == _D && _path_metrics[max_index] > (pm0 > pm1 ? pm0 : pm1)) {
					KillPath(max_index, T);
				}

				if (pm0 < pm1) {
					if (T < _D) {
						max_index = ClonePath(min_index, T);
						ExtendPath(max_index, 1, pm1);
					}

					ExtendPath(min_index, 0, pm0);
				}
				else {
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
			if (_path_lengths_hits[passedLength] >= _L)
				isTruncateNeeded = true;
		}
		// Crc check
		else {
			if (IsCrcPassed(_paths[min_index]))
				break;

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
			if (!_is_paths_active[i])
				continue;

			if (isTruncateNeeded && _path_lengths[i] <= passedLength) {
				KillPath(i, T);
				continue;
			}

			if (_path_metrics[i] < min_metric) {
				min_metric = _path_metrics[i];
				new_min_index = i;
			}
			if (_path_metrics[i] > max_metric) {
				max_metric = _path_metrics[i];
				new_max_index = i;
			}
		}		

		isPathSwitched = min_index != new_min_index;
		min_index = new_min_index;
		max_index = new_max_index;

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
	}

	return result;
}
