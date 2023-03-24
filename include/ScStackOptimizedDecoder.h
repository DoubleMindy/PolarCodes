#pragma once

#include <stack>

#include "CRC.h"
#include "BaseDecoder.h"
#include "ScOptimized.h"

using std::vector;
using std::stack;

class ScStackOptimizedDecoder : public ScOptimizedDecoder {

protected:

	std::size_t _L = 0;
	std::size_t _D = 0;

	std::size_t _kWithCrc;
	CRC * _crcPtr;

	vector<double> _path_metrics;
	vector<std::size_t> _path_lengths;
	vector<std::size_t> _path_lengths_hits;
	vector<vector<int>> _paths;
	stack<std::size_t> _inactive_path_indices;
	vector<bool> _is_paths_active;

	void recursively_calc_alpha_stack(std::size_t lambda, std::size_t phi, bool isPathSwitched);
	double calculate_step_metric(double newLlr, int decision);

	void KillPath(std::size_t index, std::size_t & T);
	std::size_t ClonePath(std::size_t index, std::size_t & T);
	void ExtendPath(std::size_t index, int bit, double pathMetric);
	void LoadPath(std::size_t index);

public:
	ScStackOptimizedDecoder(PolarCode * codePtr, int L, int D);

	std::vector<int> Decode(std::vector<double> llr) override;

	~ScStackOptimizedDecoder() {};
};