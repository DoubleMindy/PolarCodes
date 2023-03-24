#include <vector>
#include <fstream>
#include <string>

#include "../include/ScDecoderTreeMaker.h"
#include "../include/PolarCode.h"
#include "../include/CommonTransformations.h"


void DumpInfo1(std::string filename, std::string info) {
	std::ofstream resultsFileStream;
	resultsFileStream.open(filename, std::fstream::out | std::fstream::app);
	resultsFileStream << info << std::endl;

	resultsFileStream.close();
}

std::string VecToStr1(std::vector<double> vec) {
	std::string result;

	for (std::size_t i = 0; i < vec.size(); i++)
	{
		result += std::to_string(vec[i]) + ", ";
	}

	return result;
}


int main0() {
	std::vector<double> output = { 0.446035, 0.997261, 0.0692609, 0.0672047, 0.994517, 0.71172, 0.948817, 0.000515609, 0.145983, 0.00175451, 0.771573, 0.975468, 0.987648, 0.289571, 0.00332316, 0.96251 };
	/*for (std::size_t i = 0; i < output.size(); i++)
	{
		output[i] = LlrToP1(output[i]);
	}*/
	auto code = new PolarCode(4, 8, { 0, 1, 2, 4, 8, 3, 5, 9, 6, 10, 12, 7, 11, 13, 14, 15 }, {});
	auto decoderPtr = new ScDecoderTreeMaker(code, 1.0);
	decoderPtr->Decode(output);
	auto fullTree = decoderPtr->GetPathInfo();
	auto filename = "fullTree.debug";
	DumpInfo1(filename, VecToStr1(output));
	DumpInfo1(filename, fullTree);
	return 0;
}