#include "Random.h"



double RandUni()
{
	static thread_local std::mt19937 generator(std::random_device{}());
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	return distribution(generator);
}

int RandUniInt(int iMin, int iMax)
{
	static thread_local std::mt19937 generator(std::random_device{}());
	std::uniform_int_distribution<int> distribution(iMin, iMax);
	return distribution(generator);
}

double RandUni(double iMin, double iMax)
{
	static thread_local std::mt19937 generator(std::random_device{}());
	std::uniform_real_distribution<double> distribution(iMin, iMax);
	return distribution(generator);
}

double RandGaussian(double dMean, double dStdDev)
{
	static thread_local std::mt19937 generator(std::random_device{}());
	std::normal_distribution<double> distribution(dMean, dStdDev);
	return distribution(generator);
}

double RandGaussian()
{
	static thread_local std::mt19937 generator(std::random_device{}());
	std::normal_distribution<double> distribution(0.0, 1.0);
	return distribution(generator);
}

double RandGaussianNN()
{
	static thread_local std::mt19937 generator(std::random_device{}());
	std::normal_distribution<double> distribution(0.0, 1.0);
	double dVal = distribution(generator);
	dVal = (dVal / 2.0) + 0.5;
	if (dVal > 1.0)
		dVal = 1.0;
	if (dVal < 0.0)
		dVal = 0.0;
	return dVal;
}

uint8_t RandUniZeroOne()
{
	return uint8_t(RandUni() < 0.5);
}

int8_t RandUniPlusMinusOne()
{
	return int8_t(RandUni() < 0.5) * 2 - 1;
}


