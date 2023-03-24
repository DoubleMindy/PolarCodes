#pragma once

#include "CRC.h"
#include "ScDecoder.h"

// logically it is abstract class
class ScCrcAidedDecoder : public ScDecoder {

protected:
	
	void DecodeFrom(int position);
	void DecodeFromTo(int position, int endPosition);

public:
	ScCrcAidedDecoder(PolarCode * code);
	
	~ScCrcAidedDecoder() {};
};