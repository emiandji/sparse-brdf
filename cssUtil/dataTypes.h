#pragma once





//define as double or float
#define CCS_INTERNAL_TYPE double

////type used for storing nonzero coefficients
#define CCS_INTERNAL_NZ_VAL_TYPE double
#define CCS_INTERNAL_NZ_LOC_TYPE uint16_t

#define CCS_GPU_INTERNAL_TYPE double	//Changing to float will cause severe inaccuracies in the results

////define as float or half_float::half
//#define CCS_GPU_INTERNAL_TYPE float
////type used for storing location of nonzero coefficients
//#define CCS_GPU_INTERNAL_NZ_LOC_TYPE uint32_t

#if defined(_WIN32) || defined(_WIN64)
#define CCS_PAR_FOR_IDX int
#else
#define CCS_PAR_FOR_IDX size_t
#endif