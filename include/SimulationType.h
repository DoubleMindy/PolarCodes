#pragma once

#include <unordered_map>
#include <string>

enum simulatorType { MC, UnknownSimulation };

simulatorType simulatorFromString(std::string str);
std::string simulationTypeToString(simulatorType);
