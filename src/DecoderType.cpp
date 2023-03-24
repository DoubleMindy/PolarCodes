#include "../include//DecoderType.h"

decoderType decoderTypeFromString(std::string str) {
	std::unordered_map<std::string, decoderType> decoderTypeResolver = {
		{"SCRecursive", decoderType::SCRecursive},
		{"SCFano", decoderType::SCFano},
		{"SCFlip", decoderType::SCFlip},
		{"SCFlipProg", decoderType::SCFlipProg},
		{"SC", decoderType::SC},
		{"SCL", decoderType::SCList},
		{"SCListFlipStat", decoderType::SCListFlipStat},
		{"SCListFlipOracleStat", decoderType::SCListFlipOracleStat},
		{"SCFlipFano", decoderType::SCFlipFano},
		{"SCListFano", decoderType::SCListFano},
		{"SCS", decoderType::SCStack},
		{"SCOptimized", decoderType::SCOptimized},
		{"SCSOptimized", decoderType::SCStackOptimized},
		{"SCCreeper", decoderType::SCCreeper},
	};

	if (decoderTypeResolver.count(str) > 0)
		return decoderTypeResolver[str];

	return UnknownDecoder;
}

std::string decoderTypeToString(decoderType type) {

	std::unordered_map<decoderType, std::string> decoderTypeStringResolver = {
		{decoderType::SCRecursive, "SCRecursive"},
		{decoderType::SCFano, "SCFano"},
		{decoderType::SCFlip, "SCFlip"},
		{decoderType::SCFlipProg, "SCFlipProg"},
		{decoderType::SC, "SC"},
		{decoderType::SCList, "SCL"},
		{decoderType::SCListFlipStat, "SCListFlipStat"},
		{decoderType::SCListFlipOracleStat, "SCListFlipOracleStat"},
		{decoderType::SCFlipFano, "SCFlipFano"},
		{decoderType::SCListFano, "SCListFano"},
		{decoderType::SCStack, "SCS"},
		{decoderType::SCOptimized, "SCOptimized"},
		{decoderType::SCStackOptimized, "SCSOptimized"},
		{decoderType::SCCreeper, "SCCreeper"},

	};

	return decoderTypeStringResolver[type];
}