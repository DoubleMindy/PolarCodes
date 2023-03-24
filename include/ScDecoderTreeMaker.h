#pragma once

#include "ScDecoder.h"

class ScDecoderTreeMaker : public ScDecoder {

private:
	std::vector<double> _p;
	std::string _pathTrace = "";
public:
	ScDecoderTreeMaker(PolarCode * code, double approximationSigma);
	std::vector<int> Decode(std::vector<double> llr) override;
	std::string GetPathInfo();
	~ScDecoderTreeMaker() {};
};