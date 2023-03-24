#pragma once



#include "ScDecoder.h"

class ScFanoDecoder : public ScDecoder {

private:
	double _T;
	double _delta;

	std::string _pathTrace;
	std::vector<double> _p; // channel error probabilities 

	void BackwardMove(std::vector<double> & beta, std::vector<bool> & gamma, std::vector<int> & mask, int & i, double & T, double & delta, bool & B, int & j);
protected:
	
public:
	ScFanoDecoder(PolarCode * code, double T, double delta);

	std::string GetPathInfo();
	std::vector<int> Decode(std::vector<double> llr) override;

	~ScFanoDecoder() {};
};