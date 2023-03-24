#pragma once

#include "BaseDecoder.h"
#include "ScOptimized.h"
#include <deque>
#include <utility>
#include <vector>

class ScCreeperDecoder : public ScOptimizedDecoder {

protected:
	double _delta;

	// nodes as pair of (bit number, bit value)
	std::deque<std::pair<int, int>> _NP;
	std::deque<double> _TP;
	std::deque<bool> _F;

	vector<int> _path;
	vector<double> _metric;

	void recursively_calc_alpha_creeper(std::size_t lambda, std::size_t phi, bool isPathSwitched);
	double calculate_step_metric_fano(double newLlr, int decision, double pe);
	double Q(double x, double T);
	void LoadPath(int length, int lastBit);
public:
	ScCreeperDecoder(PolarCode * codePtr, double delta);

	std::vector<int> Decode(std::vector<double> llr) override;

	~ScCreeperDecoder() = default;
};
