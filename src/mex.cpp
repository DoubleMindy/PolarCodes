#include "PolarCode.h"
#include "BaseDecoder.h"
#include "ScCreeperDecoder.h"
#include "ScDecoder.h"
#include "ScFanoDecoder.h"
#include "ScFlipFanoDecoder.h"
#include "ScStackDecoder.h"
#include "ScStackOptimizedDecoder.h"
#include <mex_cast.h>
#include <mex_commands.h>
#include <stdexcept>
#include <tuple>

void mex(mexbind0x::MXCommands & m) {
	m.on("free code", [](PolarCode * code) { delete code; });
	m.on("free decoder", [](BaseDecoder * ptr) { delete ptr; });
	m.on("code by used", [](int m, std::vector<int> usedBits) {
		return new PolarCode(m, std::move(usedBits));
	});
	m.on("code by reliability",
		[](int m, int k, std::vector<int> reliabilitySequence,
			std::vector<int> crcPoly) {
			return new PolarCode(m, k, std::move(reliabilitySequence),
				std::move(crcPoly));
		});
	m.on("make creeper",
		[](PolarCode * code, double delta) {
			return static_cast<BaseDecoder *>(
				new ScCreeperDecoder(code, delta));
		});
	m.on("make fano",
		[](PolarCode * code, double T, double delta) {
			return static_cast<BaseDecoder *>(
				new ScFanoDecoder(code, T, delta));
		});
	m.on("make flip fano",
		[](PolarCode * code, double T, double delta, double
                   approximationSnr, double L) {
			return static_cast<BaseDecoder *>(
				new ScFlipFanoDecoder(code, T,
					delta,
					approximationSnr,
					L));
		});
	m.on("make stack",
		[](PolarCode * code, int L, int D) {
			return static_cast<BaseDecoder *>(
				new ScStackDecoder(code, L, D));
		});
	m.on("make stack opt",
		[](PolarCode * code, int L, int D) {
			return static_cast<BaseDecoder *>(
				new ScStackOptimizedDecoder(code, L, D));
		});
	m.on("decode", &BaseDecoder::Decode);
	m.on("get op count", [](BaseDecoder * decoder) {
		auto op_count = decoder->GetOperationsCount();
		return std::make_tuple(op_count.Iterations,
                                       op_count.Normilizer,
                                       op_count.Sums,
                                       op_count.Muls,
                                       op_count.Comps,
                                       op_count.Xors);
	});
	m.on("clear op count", [](BaseDecoder * decoder) {
            decoder->ClearOperationsCount();
        });
	if (!m.has_matched())
		throw std::invalid_argument("Command not found");
}

MEX_SIMPLE(mex);
