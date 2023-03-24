#pragma once

#include <string>
#include <unordered_map>
#include <chrono>
#include <sstream>
#include <iostream>
#include <vector>
#include "DecoderType.h"
#include "SimulationType.h"
#include "OperationsCount.h"

struct SimulationParams {
    simulatorType simulator;
    std::unordered_map < std::string, std::string > simulatorParams;
    
	decoderType decoder;
	std::unordered_map < std::string, std::string > decoderParams;
	
	// Others code types are not forseen yet
	std::string code = "Polar";
	std::unordered_map < std::string, std::string > codeParams;

	std::string resultsFilename;
	std::string additionalFilename;
	std::vector<double> snrArray;

	std::string ToString() {
		std::stringstream ss;

		ss << std::string(60, '#') << "\n";
		ss << "Prameters of simulation run\n\n";

		ss << "ResultsFilename: " + resultsFilename + "\n\n";
		ss << "AdditionalFilename: " + additionalFilename + "\n\n";
		
		ss << "Simulation type: " + simulationTypeToString(simulator) + "\n";
		ss << "Simulation parameters:\n";
		for (auto simulationParam : simulatorParams)
		{
			ss << "\t" << simulationParam.first + ": " + simulationParam.second + "\n";
		}

		ss << "Decoder type: " + decoderTypeToString(decoder) + "\n";
		ss << "Decoder params:\n";
		for (auto decoderParam : decoderParams)
		{
			ss << "\t" << decoderParam.first + ": " + decoderParam.second + "\n";
		}

		ss << "Code type: " + code + "\n";
		ss << "Code params:\n";
		for (auto codeParam : codeParams)
		{
			ss << "\t" << codeParam.first + ": " + codeParam.second + "\n";
		}

		ss << "\n";
		ss << "SNR: ";
		for (std::size_t i = 0; i < snrArray.size(); i++)
		{
			ss <<  snrArray[i] << ", ";
		}
		ss << "\n";

		ss << std::string(60, '#') << "\n";

		return ss.str();
	}
};

struct SimulationIterationResults {
	double snr;
	double ebn0;
	double sigma;

	double fer;

	int rejectionsCount;
	int testsCount;
	OperationsCount operationsCount;
	std::chrono::milliseconds elapsedTime;
	static const bool includeOperations = true;
	
	static std::string GetHeader() {
		return (std::string("SNR, EbN0, sigma, FER, rejectionsCount, testsCount, time(ms), iter")
			+ ((includeOperations) ? ", sums, muls, comps, xors, total" : ""));
	}

	std::string ToString() {
		std::stringstream ss;
		
		double sums = (double)operationsCount.Sums / operationsCount.Normilizer / 1000;
		double muls = (double)operationsCount.Muls / operationsCount.Normilizer / 1000;
		double comps = (double)operationsCount.Comps / operationsCount.Normilizer / 1000;
		double xors = (double)operationsCount.Xors / operationsCount.Normilizer / 1000;
		double total = sums + muls + comps + xors;
		double iterations = operationsCount.Iterations / operationsCount.Normilizer;

		ss << snr << ", " << ebn0 << ", " << sigma << ", " << fer << ", " << rejectionsCount
			<< ", " << testsCount << ", " << elapsedTime.count() << ", " << iterations;

		if (includeOperations)
			ss << ", " << sums << ", " << muls << ", " << comps << ", " << xors << ", " << total;
		
		return ss.str();
	}
};