#pragma once
#include"consts.h"
class Noise
{
public:
	doublevector getWhiteNoise(double sigma);
	doublevector getColorNoise(const doublevector& white_noise, double betta, double sigma);
};

