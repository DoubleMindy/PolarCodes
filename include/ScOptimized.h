#pragma once

#include "BaseDecoder.h"

using std::vector;

class ScOptimizedDecoder : public BaseDecoder {

protected:
	
	// _n = 2 ^ _m
	std::size_t _m = 0;
	std::size_t _n = 0;
	std::size_t _k = 0;

	vector<vector<double>> _alpha;
	vector<vector<vector<int>>> _beta;

	vector<int> _mask;

	double f(double left, double right);
	double g(double left, double right, int left_hard);
	int HD(double llr);

	void recursively_calc_alpha(std::size_t lambda, std::size_t phi);
	void recursively_calc_beta(std::size_t lambda, std::size_t phi);

public:
	ScOptimizedDecoder(PolarCode * codePtr);

	std::vector<int> Decode(std::vector<double> llr) override;

	~ScOptimizedDecoder() {};
};