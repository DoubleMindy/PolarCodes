#include "../include/BpskAwgnChannel.h"

#include <iostream>

int mainasd() {
	BaseChannel * channelPtr = new BpskAwgnChannel();

	std::vector<int> input = { 0 };
	std::vector<double> output;

	std::vector<double> snr_array = { -0.25, 0.25, 0.75};

	std::size_t N = 100000;



	for (std::size_t j = 0; j < snr_array.size(); j++)
	{
		int count = 0;
		for (std::size_t i = 0; i < N; i++)
		{
			channelPtr->SetSnr(snr_array[j]);

			output = channelPtr->Pass(input);

			//std::cout << output[0] << std::endl;

			if (output[0] < 0.0) {
				count++;
				//std::cout << "adsfasdfafds\n";
			}
				
		}

		std::cout << "SNR: " << snr_array[j] << ", p= " << ((double)count )/ N << std::endl;
	}
	

	
	// sqrt(pow(10, -snr / 10) / 2);

	/*double sigma = 0.727755;
	std::normal_distribution<double> normal_dist(0, sigma);
	std::random_device rd{};
	std::mt19937 gen{ rd() };

	int N = 100000;
	int count = 0;
	double x = 0.0;
	for (std::size_t i = 0; i < N; i++)
	{
		x = 1.0 + normal_dist(gen);

		if (x < 0.0) {
			count++;
		}
	}

	std::cout << "P: " << (double)count / N << std::endl;*/

	return 0;
}