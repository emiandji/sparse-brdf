#pragma once


#include "defs.h"

// template <typename T>
// double MSE(const T* pRef, const T* pRecon, size_t iLength, int64_t iQuantization);
// template <typename T>
// double MSE(const std::vector<const T*>& pRef, const std::vector<const T*>& pRecon, const std::vector<size_t>& vLength, int64_t iQuantization);
// 
// template <typename T>
// double PSNR(const T* pRef, const T* pRecon, size_t iLength, int64_t iQuantization);
// template <typename T>
// double PSNR(const std::vector<const T*>& pRef, const std::vector<const T*>& pRecon, const std::vector<size_t>& vLength, int64_t iQuantization);
// 
// template <typename T>
// double SNR(const T* pRef, const T* pRecon, size_t iLength);
// template <typename T>
// double SNR(const std::vector<const T*>& pRef, const std::vector<const T*>& pRecon, const std::vector<size_t>& vLength);



template <typename T>
double MSE(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mRef,	const Eigen::Matrix<T, Eigen::Dynamic, 1>& mRecon);
template <typename T>
double PSNR(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mRef, const Eigen::Matrix<T, Eigen::Dynamic, 1>& mRecon);
template <typename T>
double SNR(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mRef,	const Eigen::Matrix<T, Eigen::Dynamic, 1>& mRecon);

template <typename T>
double MSE(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& mRef, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& mRecon);
template <typename T>
double PSNR(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& mRef, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& mRecon);
template <typename T>
double SNR(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& mRef, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& mRecon);



template <typename T>
double MSE(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRef, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRecon);
template <typename T>
double PSNR(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRef, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRecon);
template <typename T>
double SNR(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRef, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRecon);

template <typename T>
double MSE(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& mRef,
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& mRecon);
template <typename T>
double PSNR(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& mRef,
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& mRecon);
template <typename T>
double SNR(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& mRef,
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& mRecon);

template <typename T, size_t N>
double MSE(const std::vector<boost::multi_array<T, N> >& mRef, const std::vector<boost::multi_array<T, N> >& mRecon);
template <typename T, size_t N>
double PSNR(const std::vector<boost::multi_array<T, N> >& mRef, const std::vector<boost::multi_array<T, N> >& mRecon);
template <typename T, size_t N>
double SNR(const std::vector<boost::multi_array<T, N> >& mRef,	const std::vector<boost::multi_array<T, N> >& mRecon);


template <typename T>
double PSNR(double mse);


#include "ImgQlty.inl"
