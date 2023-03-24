#include "../include/GaussianApproximation.h"
#include "../include/Exceptions.h"

#include <iostream>
#include <cstddef>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

using std::vector;

void SaveDumpGA(std::string filename, vector<int> leader) {
	std::ofstream file(filename, std::ios::trunc);

	if (!file.is_open())
		throw FileIsNotOpennedException("Could not open for writing dump file: " + filename);

	for (std::size_t i = 0; i < leader.size(); i++)
	{
		file << leader[i] << " ";
	}
}

std::vector<int> ReadSequenceFromFile1(std::string path) {
	std::vector<int> seq;
	std::string line;
	std::ifstream myFile(path);

	std::getline(myFile, line);

	int val;
	std::stringstream ss(line);

	while (ss >> val)
		seq.push_back(val);

	return seq;
}

std::vector<int> SortProb1(std::vector<double> p) {
	int n = (int)p.size();
	std::vector<int> indices(n);
	for (int i = 0; i < n; i++)
	{
		indices[i] = i;
	}

	for (int i = 0; i < n; i++)
	{
		double max = -1.0;
		int index = i;
		for (int j = i; j < n; j++)
		{
			if (p[j] > max) {
				max = p[j];
				index = j;
			}
		}

		if (index != i) {
			double p_temp = p[index];
			p[index] = p[i];
			p[i] = p_temp;

			int ind_temp = indices[index];
			indices[index] = indices[i];
			indices[i] = ind_temp;
		}
		
	}

	return indices;
}

// test correctness of GA
int main12(int argc, char* argv[]) {

	int n = 128;
	GaussianApproximation ga(1);
	
	//std::vector<int> seq = ReadSequenceFromFile1("C:\\Users\\ische\\source\\repos\\PolarCodes\\polar_sequences\\16.txt");
	std::vector<double> p(n, 0);
	for (std::size_t i = 0; i < n; i++)
	{
		p[i] = ga.GetChannelErrorProbability(i + 1, n);
	}
	std::vector<int> seq0 = SortProb1(p);
	SaveDumpGA("GaSeq.dump", seq0);
	return 0;
}
