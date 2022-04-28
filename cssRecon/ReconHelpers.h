#pragma once

#include "../cssUtil/defs.h"





template <typename T = CCS_INTERNAL_TYPE>
void SL0(Eigen::Matrix<T, Eigen::Dynamic, 1>& s,
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
	T sigma_min,
	T sigma_decrease_factor = 0.5,
	T mu_0 = 2.0,
	size_t L = 3);

template <typename T = CCS_INTERNAL_TYPE>
void SL0(Eigen::Matrix<T, Eigen::Dynamic, 1>& s,
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A_pinv,
	T sigma_min,
	T sigma_decrease_factor = 0.5,
	T mu_0 = 2.0,
	size_t L = 3);

template <typename T = CCS_INTERNAL_TYPE>
void OMP(Eigen::Matrix<T, Eigen::Dynamic, 1>& s,
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
	size_t iSparsity);

//This implementation is based on 
//http://www.cs.technion.ac.il/users/wwwb/cgi-bin/tr-get.cgi/2008/CS/CS-2008-08.pdf
//and consumes more memory but it is much faster
template <typename T = CCS_INTERNAL_TYPE>
void BatchOMP(Eigen::Matrix<T, Eigen::Dynamic, 1>& s,
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vA,
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vAtA,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
	size_t iSparsity);



#include "ReconHelpers.inl"