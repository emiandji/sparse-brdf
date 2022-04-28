#pragma once 

#include "defs.h"
#include "../cssData/Data.h"

template <typename T>
int64_t NumOfPatches(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mIm,
	const Eigen::Vector2i& vPatchSize,
	PatchType ePatchType,
	const std::vector<uint32_t>& vSlidingDis);


int64_t NumOfPatches3D(const Eigen::Vector3i& vImSize,
	const Eigen::Vector3i& vPatchSize,
	PatchType ePatchType,
	const std::vector<uint32_t>& vSlidingDis);


//Note that image patches are ordered column-wise
template <typename T>
bool im2col(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mIm,
	const Eigen::Vector2i& vPatchSize,
	PatchType ePatchType,
	const std::vector<uint32_t>& vSlidingDis,
	std::vector<Eigen::Vector2i>& vOutSldPatches,
	T tPadVal,
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mOut,
	size_t& mblocks, size_t& nblocks);


template <typename T>
bool col2im(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mCol,
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mOut,
	const Eigen::Vector2i& vPatchSize,
	const Eigen::Vector2i& vImSize,
	PatchType ePatchType,
	const std::vector<Eigen::Vector2i>& vInSldPatches);



//Note that image patches are ordered column-wise
template <typename T>
bool im2col3D(const boost::multi_array<T, 3>& mIm,
	const Eigen::Vector3i& vPatchSize,
	PatchType ePatchType,
	const std::vector<uint32_t>& vSlidingDis,
	std::vector<Eigen::Vector3i>& vOutSldPatches,
	T tPadVal,
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mOut,
	size_t& mblocks, size_t& nblocks, size_t& pblocks);



template <typename T>
bool col2im3D(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mCol,
	boost::multi_array<T, 3>& mOut,
	const Eigen::Vector3i& vPatchSize,
	const Eigen::Vector3i& vImSize,
	PatchType ePatchType,
	const std::vector<Eigen::Vector3i>& vInSldPatches);




#include "Patchifier.inl"