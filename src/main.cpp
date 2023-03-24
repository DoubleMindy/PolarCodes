#include <iostream>
#include <vector>
#include <chrono>
#include "../include/Simulate.h"

void PrintUsage() {
#ifdef __linux__ 
	std::cout << "Usage: LDPC config_path" << std::endl;
#elif _WIN32
	std::cout << "Usage: LDPC.exe config_path" << std::endl;
#else

#endif
}

void testRoutine() {

	auto t1 = std::chrono::steady_clock::now();
	
	int N = 1000000000;

	int max_m = 0;
	int max_element;

	for (int i = 0; i < N; i++)
	{
		if (max_m < i) {
			max_m = i;
			max_element = i;
		}
	}
	
	auto t2 = std::chrono::steady_clock::now();

	auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

	std::cout << "Time: " << dt << std::endl;

	return;
}

int main(int argc, char* argv[]) {

	if (argc != 2) {
		PrintUsage();
		return 1;
	}
	
	bool test = false;
	if (test){
		testRoutine();
	}

	bool isMainRun = true;
	if (isMainRun) {
		Simulate(argv[1]);
	}

	return 0;
}