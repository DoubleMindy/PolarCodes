#pragma once

#include "ScCrcAidedDecoder.h"

struct FanoState {
	std::vector<std::vector<double>> beliefTree;
	std::vector<std::vector<int>> uhatTree;
	std::vector<int> x;

	int i;
	int j;
	bool B;
	
	std::vector<double> betaJ;
	std::vector<double> betaI;
	std::vector<bool> gamma;

	// list mechanism
	double metric;
	bool isVisited;
};

class ScListFanoDecoder : public ScCrcAidedDecoder {

private:
	double _T;
	double _delta;
	int _L;

	std::vector<double> _p; // channel error probabilities 

	// list stuff
	std::vector<FanoState> _states;

	// current state
	int _i;
	int _j;
	bool _B;

	std::vector<double> _betaJ;
	std::vector<double> _betaI;
	std::vector<bool> _gamma;

	void UpdateT(double & T, double & delta, double & tau);
	void BackwardMove(std::vector<double> & beta, std::vector<bool> & gamma, int & i, double & T, double & delta, bool & B, int & j);

	FanoState SaveState();
	void LoadState(FanoState state);
public:
	ScListFanoDecoder(PolarCode * code, double T, double delta, double approximationSnr, int L);

	// All fano decoders works only with P1 domain
	std::vector<int> Decode(std::vector<double> p1) override;
	~ScListFanoDecoder() {};
};