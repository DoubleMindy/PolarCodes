#include <random>
#include <algorithm>
#include <string>
#include <fstream>
#include <chrono>
#include "OperationsCount.h"
#include "CommonTransformations.h"
#include "SimulationParameters.h"
#include "MonteCarloSimulator.h"
#include "ScFanoDecoder.h"
#include "ScDecoder.h"
#include "ScListDecoder.h"
#include "ScStackDecoder.h"
#include "ScFanoDecoder.h"

MonteCarloSimulator::MonteCarloSimulator(int maxTestsCount,
	int maxRejectionsCount,
	PolarCode * codePtr,
	Encoder * encoderPtr,
	BaseChannel * channelPtr,
	BaseDecoder * decoderPtr,
	bool isSigmaDependOnR) : BaseSimulator(codePtr, encoderPtr, channelPtr, decoderPtr, isSigmaDependOnR)
{
	_maxTestsCount = maxTestsCount;
	_maxRejectionsCount = maxRejectionsCount;
}

void DumpInfo(std::string filename, std::string info) {
	std::ofstream resultsFileStream;
	resultsFileStream.open(filename, std::fstream::out | std::fstream::app);
	resultsFileStream << info << std::endl;

	resultsFileStream.close();
}

template<class T>
std::string VecToStr(std::vector<T> vec) {
	std::ostringstream streamObj;
	std::string result;

	for (std::size_t i = 0; i < vec.size(); i++)
	{
		streamObj << vec[i] << ", ";
		//result += std::to_string(vec[i]) + ", ";
	}

	return streamObj.str();
}

#ifdef PARALLEL_DECODER1
std::vector<int> ReadSequenceFromFileParallel(std::string path) {
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
#endif // PARALLEL_DECODER
SimulationIterationResults MonteCarloSimulator::Run(double snr)
{	
	SimulationIterationResults result;

	std::size_t n = _codePtr->N();
	std::size_t k = _codePtr->k();
	std::vector<int> word(k, 0);
	std::vector<int> codeword(n, 0);
	std::vector<double> channelOuput(n, 0);
	std::vector<int> decoded(n, 0);
	
	auto t1 = std::chrono::steady_clock::now();

	double sigma = GetSigma(snr, (double)k / n);
	int tests = 0;
	int wrong_dec = 0;
	//int wrong_bits = 0;

	std::random_device randomDevice;
	std::uniform_int_distribution<> uniform_discrete_dist(0, 1);

#ifdef PARALLEL_DECODER1
	// HERE parallel decoder
	std::vector<double> channelOuput1(n, 0);
	std::vector<int> parallelDecoded;
	std::size_t m = _codePtr->m();
	PolarCode * parallelCodePtr = new PolarCode(m, k, ReadSequenceFromFileParallel("polar_sequences/" + std::to_string(n) + ".txt"), _codePtr->CrcPoly());
	ScFanoDecoder parallelDecoder(parallelCodePtr, 0.0, 1);
	//////
#endif // PARALLEL_DECODER

	_decoderPtr->SetSigma(sigma);
	_channelPtr->SetSnr(snr);

	_decoderPtr->ClearOperationsCount();

	while ((tests < _maxTestsCount || _maxTestsCount == -1) && (wrong_dec < _maxRejectionsCount)) {
		tests++;

		std::generate(word.begin(), word.end(), [&]() { return uniform_discrete_dist(randomDevice); });

		codeword = _encoderPtr->Encode(word);
		// Give answer to a decoder for debugging or statistic retreving
		_decoderPtr->SetDecoderAnswer(_encoderPtr->PolarTransform(codeword));

		channelOuput = _channelPtr->Pass(codeword);
		//channelOuput = { 0.446035, 0.997261, 0.0692609, 0.0672047, 0.994517, 0.71172, 0.948817, 0.000515609, 0.145983, 0.00175451, 0.771573, 0.975468, 0.987648, 0.289571, 0.00332316, 0.96251 };
		channelOuput = { 2995.41336494294, -3013.35534387100, -3044.97136616593, -2997.33560358008, 2946.57282761782, 2997.97709825946, 2977.83928017960, 2958.05823010378, 2942.86691749557, -3107.50244143779, -2942.69771009553, -2965.05754714817, 2925.56742136829, 3107.54188636563, -3003.08819043366, -2923.81156604213, -2971.15697300063, 2919.27050225658, -3028.75280661210, 3010.10814628387, -3133.55318745752, -2893.84767014759, -2952.86760652498, 2954.26455989276, 2910.91927110516, -3045.21603406567, 3034.01199506910, -3051.03916987710, 2940.26909546812, 3013.22405912085, -2957.76711880748, 2924.73999739717 };

		decoded = _decoderPtr->Decode(channelOuput);
		//try {
		//
		//}
		//catch (std::exception e) {
		//	std::string filename = "results/Creeper.debug";
		//	DumpInfo(filename, VecToStr<double>(channelOuput));
		//	//break;
		//}
#ifdef PARALLEL_DECODER1
		// HERE parallel decoder to comparision, SC - worse decoder
		/*for (std::size_t i = 0; i < n; i++)
		{
			channelOuput1[i] = LlrToP1(channelOuput[i]);
		}*/
		parallelDecoded = parallelDecoder.Decode(channelOuput);
		if (word != decoded && parallelDecoded == word) {
			std::cout << "Find" << std::endl;
			std::string debugInfo = parallelDecoder.GetPathInfo();
			std::string filename = "results/dump_SCS1.debug";
			DumpInfo(filename, VecToStr<double>(channelOuput));
			DumpInfo(filename, VecToStr<int>(word));
			DumpInfo(filename, VecToStr<int>(decoded));
			DumpInfo(filename, VecToStr<int>(parallelDecoded));
			DumpInfo(filename, debugInfo);
			DumpInfo(filename, "");
		}
		///////////
#endif //PARALLEL_DECODER

		if (tests % 1000 == 0)
			std::cout << tests << std::endl;
		
		if (decoded != word) {
			wrong_dec += 1;

			/*for (std::size_t i = 0; i < n; i++)
			{
				if (decoded[i] != word[i])
					wrong_bits += 1;
			}*/
		}
	}

	//std::cout << "BER: " << (double)wrong_bits / n / tests << std::endl;

	OperationsCount operationsCount = _decoderPtr->GetOperationsCount();

	auto t2 = std::chrono::steady_clock::now();

	result.snr = snr;
	result.ebn0 = GetEbN0(snr, k, n);
	result.sigma = sigma;

	result.fer = (double)wrong_dec / tests;

	result.elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
	result.testsCount = tests;
	result.rejectionsCount = wrong_dec;
	
	result.operationsCount = operationsCount;

	return result;
}
