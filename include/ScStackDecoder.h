#pragma once

#include "ScCrcAidedDecoder.h"

class ScStackDecoder : public ScCrcAidedDecoder {
protected:
	int _L;
	int _D;

	std::vector<std::vector<std::vector<double>>> _beliefTrees;
	std::vector<std::vector<std::vector<int>>> _uhatTrees;
	std::vector<std::vector<int>> _candidates;
	std::vector<std::pair<double, int>> _metrics; // metric and index of structures
	std::vector<int> _current_bits;
	std::vector<int> _paths_limits; // i-th element stands for a path of length i+1

	void PassDownAll(std::size_t iter);
	void PassUpAll(std::size_t iter);
	void PassDownStackLeader(std::size_t iter);
	void PassUpStackSelectevely(std::size_t iter, std::vector<int> elements);
	double StepMetric(double belief, int decision);
	void EliminatePaths(std::size_t iter);
public:
	ScStackDecoder(PolarCode * code, int L, int D);
	std::vector<int> Decode(std::vector<double> belief) override;
	~ScStackDecoder() {};
};