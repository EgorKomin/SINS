#include "Noise.h"



doublevector Noise::getWhiteNoise(double sigma)
{
	// белый шум
	double a, b;
	a = -sigma;
	b = sigma;

	auto str = consts::getconsts();
	auto dt = 1 / str.f;

	doublevector white_noise;

	for (int i = 0; i < str.T / dt; i++)
	{
		white_noise.push_back(1. * (b - a) * rand() / RAND_MAX + a);
	}

	return white_noise;
}

doublevector Noise::getColorNoise(const doublevector &white_noise, double betta, double sigma)
{
	//генерация шумов акселерометров или гироскопов

	auto str = consts::getconsts();
	auto dt = 1 / str.f;

	double noise;

	doublevector color_noise;
	color_noise.push_back(0);

	for (int i = 1; i < str.T / dt; i++)
	{
		noise = (1 - betta * dt) * color_noise[i - 1] + sigma * sqrt(2 * betta * dt) * white_noise[i - 1];
		color_noise.push_back(noise);
	}

	return color_noise;
}
