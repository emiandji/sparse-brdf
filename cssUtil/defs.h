#pragma once 


#include <iostream>
#include <stdint.h>
#include <vector>
#include <string>
#include <assert.h>
#include <cmath>
#include <fstream>
#include <numeric>
#include <limits>
#include <functional>
#include <omp.h>
#include <queue>
#include <map>
#include <iterator>
#include <thread>
#include <boost/dynamic_bitset.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/filesystem.hpp>
#include <boost/random.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/multi_array.hpp>
#include <boost/progress.hpp>
#include <boost/numeric/conversion/bounds.hpp>
#include <Eigen/Core>
#include <Eigen/LU> 
#include <Eigen/SVD>
#include <Eigen/Sparse>
#include "dataTypes.h"







#if defined(_WIN32) || defined(_WIN64)
#define CCS_IS_WINDOWS
#elif defined(__linux__)
#define CCS_IS_LINUX
#elif defined(__APPLE__)
#define CCS_IS_APPLE
#if !(defined(__i386__) || defined(__amd64__))
#define CCS_IS_APPLE_PPC
#else
#define CCS_IS_APPLE_X86
#endif
#elif defined(__OpenBSD__)
#define CCS_IS_OPENBSD
#endif



#if !defined(CCS_IS_APPLE) && !defined(CCS_IS_OPENBSD)
#include <malloc.h> // for _alloca, memalign
#endif
#if !defined(CCS_IS_WINDOWS) && !defined(CCS_IS_APPLE) && !defined(CCS_IS_OPENBSD)
#include <alloca.h>
#endif


#if defined(CCS_IS_WINDOWS)
#define alloca _alloca
#endif
#ifndef CCS_L1_CACHE_LINE_SIZE
#define CCS_L1_CACHE_LINE_SIZE 64
#endif
#ifndef CCS_POINTER_SIZE
#if defined(__amd64__) || defined(_M_X64)
#define CCS_POINTER_SIZE 8
#elif defined(__i386__) || defined(_M_IX86)
#define CCS_POINTER_SIZE 4
#endif
#endif
#ifndef CCS_HAS_64_BIT_ATOMICS
#if (CCS_POINTER_SIZE == 8)
#define CCS_HAS_64_BIT_ATOMICS
#endif
#endif // CCS_HAS_64_BIT_ATOMICS



#if defined(CCS_IS_WINDOWS)
#include <windows.h>
#else
#include <sys/ioctl.h>
#include <unistd.h>
#include <errno.h>
#endif 


#ifdef NDEBUG
#define Assert(expr) ((void)0)
#else
#define Assert(expr) \
    ((expr) ? (void)0 : \
        Severe("Assertion \"%s\" failed in %s, line %d", \
               #expr, __FILE__, __LINE__))
#endif // NDEBUG

//boost::c_storage_order()
#define SAFE_DELETE(x) if(x) {delete x; x = NULL;} 
#define SAFE_DELETE_ARRAY(x) if(x) {delete [] x; x = NULL;}
#define CCS_TENSOR_STORAGE_ORDER (boost::fortran_storage_order())
#define FISHER_NB_BUFFER_SIZE 1024

template <typename T> inline
T absT(const T& v) { return v < 0 ? -v : v; }


enum QMetric   { CCS_QMETRIC_MSE = 0, CCS_QMETRIC_PSNR = 1, CCS_QMETRIC_SSIM = 2 };


inline void YCCToRGB(const double ycc[3], double rgb[3])
{
	rgb[0] = ycc[0] - 0.0000012*ycc[1] + 1.4019995*ycc[2];
	rgb[1] = ycc[0] - 0.3441356*ycc[1] - 0.7141361*ycc[2];
	rgb[2] = ycc[0] + 1.7720000*ycc[1] + 0.0000004*ycc[2];
}

inline void RGBToYCC(const double rgb[3], double ycc[3])
{
	ycc[0] = 0.299000*rgb[0] + 0.587000*rgb[1] + 0.114000*rgb[2];
	ycc[1] = -0.168736*rgb[0] - 0.331264*rgb[1] + 0.500000*rgb[2];
	ycc[2] = 0.500000*rgb[0] - 0.418688*rgb[1] - 0.081312*rgb[2];
}


template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
	if (!v.empty()) {
		out << '[';
		std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
		out << "\b\b]";
	}
	return out;
}
