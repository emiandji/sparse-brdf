

#pragma once


#include <omp.h>
#include <random>


double RandUni();

int RandUniInt(int iMin, int iMax);

double RandUni(double iMin, double iMax);

double RandGaussian(double dMean, double dStdDev);

double RandGaussian();

double RandGaussianNN();

uint8_t RandUniZeroOne();

int8_t RandUniPlusMinusOne();
