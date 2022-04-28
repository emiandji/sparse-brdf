#pragma once


#include "defs.h"
#include "Random.h"
#include "ProgressReporter.h"

/*
#ifdef SVD_CLAPACK
// extern "C"
// {
#include <f2c.h>
#include <blaswrap.h>
#include <clapack.h>
//}
#elif SVD_LAPACKE
#include <lapacke.h>
#elif SVD_MKL
#include <mkl.h>
#else
#endif
*/

#ifdef abs
#undef abs
#endif
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif



/* solve L*x = b */
template <typename T>
void backsubst_L(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& L,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& b, Eigen::Matrix<T, Eigen::Dynamic, 1>& x, size_t k);

/* solve L'*x = b */
template <typename T>
void backsubst_Lt(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& L,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& b, Eigen::Matrix<T, Eigen::Dynamic, 1>& x, size_t k);

/* solve U*x = b */
template <typename T>
void backsubst_U(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& U,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& b, Eigen::Matrix<T, Eigen::Dynamic, 1>& x, size_t k);

/* solve U'*x = b */
template <typename T>
void backsubst_Ut(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& U,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& b, Eigen::Matrix<T, Eigen::Dynamic, 1>& x, size_t k);

/* back substitution solver */
template <typename T>
void backsubst(char ul, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& b, Eigen::Matrix<T, Eigen::Dynamic, 1>& x, size_t k);

/* solve equation set using cholesky decomposition */
template <typename T>
void cholsolve(char ul, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& b, Eigen::Matrix<T, Eigen::Dynamic, 1>& x, size_t k);



template <typename T>
T MutualCoherenceSum(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mLeft, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRight);

template <typename T>
T MutualCoherenceAvg(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mLeft, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRight);

template <typename T>
T MutualCoherenceExpSum(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mLeft, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRight);

template <typename T>
T MutualCoherenceExpAvg(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mLeft, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRight);

template <typename T>
T MutualCoherenceExpProd(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mLeft, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRight);

template <typename T>
T MutualCoherenceMax(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mLeft, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRight);

template <typename T>
T MutualCoherenceExpMax(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mLeft, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRight);




//Matrix to vector conversion and vice versa
template <typename Derived1, typename Derived2>
void Mat1Dto2D(const Eigen::DenseBase<Derived1>& mIn, size_t iRows, size_t iCols, Eigen::DenseBase<Derived2> const & mOut);

template <typename Derived1, typename Derived2>
void Mat2Dto1D(const Eigen::DenseBase<Derived1>& mIn, Eigen::DenseBase<Derived2> const & mOut);

//Vector, matrix and tensor distance metrics (only the internal type is allowed)
template <typename Derived1, typename Derived2>
CCS_INTERNAL_TYPE EuclidDist2(const Eigen::MatrixBase<Derived1>& coord1, const Eigen::MatrixBase<Derived2>& coord2);

template <size_t N>
CCS_INTERNAL_TYPE EuclidDist2(const boost::multi_array<CCS_INTERNAL_TYPE, N>& coord1, const boost::multi_array<CCS_INTERNAL_TYPE, N>& coord2);

template <typename Derived1, typename Derived2>
CCS_INTERNAL_TYPE EuclidDist1(const Eigen::MatrixBase<Derived1>& coord1, const Eigen::MatrixBase<Derived2>& coord2);

template <size_t N>
CCS_INTERNAL_TYPE EuclidDist1(const boost::multi_array<CCS_INTERNAL_TYPE, N>& coord1, const boost::multi_array<CCS_INTERNAL_TYPE, N>& coord2);

template <typename Derived1, typename Derived2>
CCS_INTERNAL_TYPE EuclidDistInf(const Eigen::MatrixBase<Derived1>& coord1, const Eigen::MatrixBase<Derived2>& coord2);

template <size_t N>
CCS_INTERNAL_TYPE EuclidDistInf(const boost::multi_array<CCS_INTERNAL_TYPE, N>& coord1, const boost::multi_array<CCS_INTERNAL_TYPE, N>& coord2);



////Singular Value Decomposition of a matrix
//template <typename T>
//void SVD(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mMat,
//	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mU,
//	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mVt,
//	Eigen::Matrix<T, Eigen::Dynamic, 1>& vS,
//	bool bThin);

//Higher order SVD
template <typename T, size_t N>
void HOSVD(const boost::multi_array<T, N>& inputTensor, std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vU, boost::multi_array<T, N>& coreTensor);

template <typename T>
void TensorProdALS(const boost::multi_array<T, 3>& inputTensor, size_t iMode, std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU);
template <typename T>
void TensorProdALS(const boost::multi_array<T, 4>& inputTensor, size_t iMode, std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU);
template <typename T>
void TensorProdALS(const boost::multi_array<T, 5>& inputTensor, size_t iMode, std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU);
template <typename T>
void TensorProdALS(const boost::multi_array<T, 6>& inputTensor, size_t iMode, std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU);
template <typename T, size_t N>
int TensorALSR1(const boost::multi_array<T, N>& inputTensor, size_t iMaxIter, T thresh, std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU);

////Higher order SVD (truncated)
//template <typename T, size_t N>
//void HOSVDtruncated(const boost::multi_array<T, N>& inputTensor,
//	std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vU,
//	boost::multi_array<T, N>& coreTensor);

// template<typename T>
// void PowerIter(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mMat,
// 	double& sigma,
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1>& vU,
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1>& vV,
// 	CCS_INTERNAL_TYPE thresh,
// 	size_t iMaxIter);
// 
// template <typename T>
// void PowerIterFull(uint32_t iNumTerms,
// 	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& F,
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, Eigen::Dynamic>& U,
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, Eigen::Dynamic>& V,
// 	CCS_INTERNAL_TYPE threshold,
// 	uint32_t iMaxIter,
// 	double& dTimeDelta);
// 
// template <typename T>
// void PowerIterFull(uint32_t iNumTerms,
// 	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& F,
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, Eigen::Dynamic>& U,
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, Eigen::Dynamic>& V,
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1>& S,
// 	CCS_INTERNAL_TYPE threshold, uint32_t iMaxIter, double& dTimeDelta);
// 
// template <typename T>
// void PowerIterFullOld(uint32_t iNumTerms,
// 	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& F,
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, Eigen::Dynamic>& U,
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, Eigen::Dynamic>& V,
// 	CCS_INTERNAL_TYPE threshold,
// 	uint32_t iMaxIter,
// 	double& dTimeDelta);


////Pseudo-inverse
//template<typename _Matrix_Type_>
//void PInv(const _Matrix_Type_& mIn, _Matrix_Type_& mOut, double epsilon = boost::numeric::bounds<double>::smallest());


//Kronecker product of two matrices
template <typename Derived1, typename Derived2, typename Derived3>
void Kron(const Eigen::MatrixBase<Derived1>& mLeft,
	const Eigen::MatrixBase<Derived2>& mRight,
	Eigen::MatrixBase<Derived3> const& mResC);

template <typename T>
void KronSelective(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& matrices,
	const std::vector<size_t>& indices, 
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRes);

template <typename T>
double KronSelectiveElem(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& matrices,
	const std::vector<size_t>& indices, size_t iRow, size_t iCol);

template <typename T>
double KronSelectiveElem(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& matrices,
	const std::vector<size_t>& indices, size_t iRow, size_t iCol, const std::vector<size_t>& vRowMult, const std::vector<size_t>& vColMult);

template <typename T, typename Idx>
void MultKronSelBySpMatTran(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& matrices,
	const std::vector<size_t>& indices, Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc, 
	Eigen::Matrix<T, Eigen::Dynamic, 1>& val, size_t iNumColsForResult, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& ret);

template <typename T, typename Idx>
void MultKronSelBySpMatTran(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& matrices,
	const std::vector<size_t>& indices, Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& val, const std::vector<size_t>& vRowMult, const std::vector<size_t>& vColMult, 
	size_t iNumColsForResult, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& ret);

template <typename T>
void KhatriRaoCol(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mLeft,
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRight,
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRes);

template <typename T>
void KhatriRaoColSelective(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& matrices,
	const std::vector<size_t>& indices,
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRes);

template <typename Idx>
void SpMatTrans(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& nzLoc, Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& result);

template <typename Idx>
void SpMatTransInplace(Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& nzLoc);


//Clamp function for 1D data (in-place)
template <typename T>
void Clamp1D(Eigen::Matrix<T, Eigen::Dynamic, 1>& mIn, T minVal, T maxVal);

//Clamp function for 2D data (in-place)
template <typename T>
void Clamp2D(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mIn, T minVal, T maxVal);

//Clamp function for nD data (in-place)
template <typename T, size_t N>
void ClampnD(boost::multi_array<T, N>& mIn, T minVal, T maxVal);



//---------------------------------------------------------------------------------------------------
//Some tensor functions------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

template <typename inType, typename dimType, typename outType>
void Ind2Sub(inType idx, const std::vector<dimType>& vDims, std::vector<outType>& vIdx);

template <typename inType, typename dimType, typename outType>
void Ind2Sub(inType idx, const Eigen::Matrix<dimType, Eigen::Dynamic, 1>& vDims, Eigen::Matrix<outType, Eigen::Dynamic, 1>& vIdx);

template <typename inType, typename dimType, typename outType>
void Ind2Sub(inType idx, const std::vector<dimType>& vDims, Eigen::Matrix<outType, Eigen::Dynamic, 1>& vIdx);

template <typename inType, typename dimType>
size_t Sub2Ind(std::vector<inType>& vIdx, const std::vector<dimType>& vDims);

template <typename inType, typename dimType>
size_t Sub2Ind(const inType* vIdx, const dimType* vDims, size_t iNumDims);

//Assignment-----------------------------------------------------------------------------------------

template <typename T, size_t N>
void TensorSetOne(boost::multi_array<T, N>& tensor);

template <typename T, size_t N>
void TensorSetZero(boost::multi_array<T, N>& tensor);

template <typename T, size_t N>
void TensorSetConst(boost::multi_array<T, N>& tensor, T val);

template <typename T, size_t N>
void TensorSetRandUni(boost::multi_array<T, N>& tensor, T tMin, T tMax);

template <typename T, size_t N>
void TensorSetRandGaussian(boost::multi_array<T, N>& tensor, T mean, T var);

template <typename T, size_t N>
void TensorSetRandNormal(boost::multi_array<T, N>& tensor);

template <typename T, size_t N>
void TensorCopy(const boost::multi_array<T, N>& from, boost::multi_array<T, N>& to);

template <typename T, size_t N, typename Derived2>
void TensorCopy(const boost::multi_array<T, N>& from, Eigen::DenseBase<Derived2> const & to);

template <typename T, size_t N>
void TensorCopy(const boost::multi_array<T, N>& from, std::vector<T>& to);

template <typename T, size_t N, typename Derived2>
void TensorCopy(const Eigen::DenseBase<Derived2>& from, boost::multi_array<T, N>& to);

template <typename T, size_t N>
void TensorCopy(const std::vector<T>& from, boost::multi_array<T, N>& to);

template <typename T, size_t N, typename Derived2>
void TensorCopy(const Eigen::DenseBase<Derived2>& from, boost::multi_array<T, N>& to, const std::vector<typename boost::multi_array<T, N>::size_type>& shape);

template <typename T, size_t N>
void TensorCopy(const std::vector<T>& from, boost::multi_array<T, N>& to, const std::vector<typename boost::multi_array<T, N>::size_type>& shape);



//Tensor-Matrix product------------------------------------------------------------------------------

inline void TensorProdMatSize(size_t* shape, size_t numdim, size_t mode, const Eigen::Vector2i& matDims, size_t* outShape)
{
	std::copy(shape, shape + numdim, outShape);
	outShape[mode] = matDims(0);
}

//dummy
template <typename T, typename Derived1>
void TensorProdMat(const boost::multi_array<T, 1>& tensor,
	const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 1>::index mode,
	boost::multi_array<T, 1>& ret)	{}

//dummy
template <typename T, typename Derived1>
void TensorProdMat(const boost::multi_array<T, 2>& tensor,
	const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 2>::index mode,
	boost::multi_array<T, 2>& ret) {}

template <typename T, typename Derived1>
void TensorProdMat(const boost::multi_array<T, 3>& tensor,
	const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 3>::index mode,
	boost::multi_array<T, 3>& ret);

template <typename T, typename Derived1>
void TensorProdMat(const boost::multi_array<T, 4>& tensor,
	const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 4>::index mode,
	boost::multi_array<T, 4>& ret);

template <typename T, typename Derived1>
void TensorProdMat(const boost::multi_array<T, 5>& tensor,
	const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 5>::index mode,
	boost::multi_array<T, 5>& ret);

template <typename T, typename Derived1>
void TensorProdMat(const boost::multi_array<T, 6>& tensor,
	const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 6>::index mode,
	boost::multi_array<T, 6>& ret);

template <typename T>
void TensorProdMatMulti(const boost::multi_array<T, 3>& tensor,
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vMats,
	boost::multi_array<T, 3>& ret);



//Tensor unfolding functions-------------------------------------------------------------------------

inline Eigen::Vector2i TensorUnfoldSize(const size_t* shape, size_t numdim, size_t mode)
{
	size_t J = 1;
	for (size_t i = 0; i < numdim; ++i)
	{
		if (i == mode)
			continue;
		J *= shape[i];
	}
	return Eigen::Vector2i(shape[mode], J);
}

//dummy
template <typename T, typename Derived1>
void TensorUnfold(const boost::multi_array<T, 1>& tensor,
	typename boost::multi_array<T, 1>::index mode,
	Eigen::MatrixBase<Derived1>& ret)	{}

//dummy
template <typename T, typename Derived1>
void TensorUnfold(const boost::multi_array<T, 2>& tensor,
	typename boost::multi_array<T, 2>::index mode,
	Eigen::MatrixBase<Derived1>& ret)	{}

template <typename T, typename Derived1>
void TensorUnfold(const boost::multi_array<T, 3>& tensor,
	typename boost::multi_array<T, 3>::index mode,
	Eigen::MatrixBase<Derived1>& ret);

template <typename T, typename Derived1>
void TensorUnfold(const boost::multi_array<T, 4>& tensor,
	typename boost::multi_array<T, 4>::index mode,
	Eigen::MatrixBase<Derived1>& ret);

template <typename T, typename Derived1>
void TensorUnfold(const boost::multi_array<T, 5>& tensor,
	typename boost::multi_array<T, 5>::index mode,
	Eigen::MatrixBase<Derived1>& ret);

template <typename T, typename Derived1>
void TensorUnfold(const boost::multi_array<T, 6>& tensor,
	typename boost::multi_array<T, 6>::index mode,
	Eigen::MatrixBase<Derived1>& ret);

//for sparse matrices
template <typename T>
void TensorUnfold(const boost::multi_array<T, 3>& tensor,
	typename boost::multi_array<T, 3>::index mode,
	Eigen::SparseMatrix<T>& ret, size_t iNNZ);

template <typename T>
void TensorUnfold(const boost::multi_array<T, 4>& tensor,
	typename boost::multi_array<T, 4>::index mode,
	Eigen::SparseMatrix<T>& ret, size_t iNNZ);

template <typename T>
void TensorUnfold(const boost::multi_array<T, 5>& tensor,
	typename boost::multi_array<T, 5>::index mode,
	Eigen::SparseMatrix<T>& ret, size_t iNNZ);

template <typename T>
void TensorUnfold(const boost::multi_array<T, 6>& tensor,
	typename boost::multi_array<T, 6>::index mode,
	Eigen::SparseMatrix<T>& ret, size_t iNNZ);

template <typename T, typename Idx>
void SparseTensorUnfold(Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc, Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	const std::vector<size_t>& vShape, size_t iMode, Eigen::SparseMatrix<T>& ret);

template <typename Idx>
void SparseTensorUnfold(Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc, 
	const std::vector<size_t>& vShape, size_t iMode, Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& retLoc);

//dummy
template <typename T, typename Derived1>
void TensorFold(const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 1>::index mode,
	boost::multi_array<T, 1>& tensor) {}

//dummy
template <typename T, typename Derived1>
void TensorFold(const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 2>::index mode,
	boost::multi_array<T, 2>& tensor) {}

template <typename T, typename Derived1>
void TensorFold(const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 3>::index mode,
	boost::multi_array<T, 3>& tensor);

template <typename T, typename Derived1>
void TensorFold(const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 4>::index mode,
	boost::multi_array<T, 4>& tensor);

template <typename T, typename Derived1>
void TensorFold(const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 5>::index mode,
	boost::multi_array<T, 5>& tensor);

template <typename T, typename Derived1>
void TensorFold(const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 6>::index mode,
	boost::multi_array<T, 6>& tensor);

//Tensor math functions-----------------------------------------------------------------------------
template <typename T, size_t N>
void TensorAdd(const boost::multi_array<T, N>& input1, boost::multi_array<T, N>& ret);

template <typename T, size_t N>
void TensorAdd(const boost::multi_array<T, N>& input1, const boost::multi_array<T, N>& input2, boost::multi_array<T, N>& ret);

template <typename T, size_t N>
void TensorAdd(T input1, boost::multi_array<T, N>& ret);

template <typename T, size_t N>
void TensorSubtract(const boost::multi_array<T, N>& input1, boost::multi_array<T, N>& ret);

template <typename T, size_t N>
void TensorSubtract(const boost::multi_array<T, N>& input1, const boost::multi_array<T, N>& input2, boost::multi_array<T, N>& ret);

template <typename T, size_t N>
void TensorSubtract(T input1, boost::multi_array<T, N>& ret);

template <typename T, size_t N>
void TensorProdVal(const boost::multi_array<T, N>& input1, T val, boost::multi_array<T, N>& ret);

template <typename T, size_t N>
void TensorProdVal(T input1, boost::multi_array<T, N>& ret);

template <typename T, size_t N>
void TensorPow(boost::multi_array<T, N>& tns, T pwr);

template <typename T, size_t N>
void TensorPow(boost::multi_array<T, N>& tns, T pwr, boost::multi_array<T, N>& ret);

//Tensor norms (only the internal type is allowed)---------------------------------------------------
template <size_t N>
CCS_INTERNAL_TYPE TensorNorm2(const boost::multi_array<CCS_INTERNAL_TYPE, N>& tensor);

template <size_t N>
CCS_INTERNAL_TYPE TensorNorm2Squared(const boost::multi_array<CCS_INTERNAL_TYPE, N>& tensor);

template <size_t N>
CCS_INTERNAL_TYPE TensorNorm1(const boost::multi_array<CCS_INTERNAL_TYPE, N>& tensor);

template <size_t N>
CCS_INTERNAL_TYPE TensorNormInf(const boost::multi_array<CCS_INTERNAL_TYPE, N>& tensor);


//Convert sparse tensor into a location and value pair
template <typename T, typename Idx, typename Derived1>
size_t DenseToSparse(const Eigen::MatrixBase<Derived1>& mat, Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& val, size_t iNNZ);

template <typename T, typename Idx>
size_t DenseToSparse(const boost::multi_array<T, 3>& tensor, Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& val, size_t iNNZ);

template <typename T, typename Idx>
size_t DenseToSparse(const boost::multi_array<T, 4>& tensor, Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& val, size_t iNNZ);

template <typename T, typename Idx>
size_t DenseToSparse(const boost::multi_array<T, 5>& tensor, Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& val, size_t iNNZ);

template <typename T, typename Idx>
size_t DenseToSparse(const boost::multi_array<T, 6>& tensor, Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& val, size_t iNNZ);

//Convert sparse tensor to dense
template <typename T, typename Idx>
void SparseToDense(const Eigen::Matrix<Idx, Eigen::Dynamic, 1>& loc,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& val, Eigen::Matrix<T, Eigen::Dynamic, 1>& ret);

template <typename T, typename Idx>
void SparseToDense(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& val, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& ret);

template <typename T, typename Idx>
void SparseToDense(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& val, boost::multi_array<T, 3>& ret);

template <typename T, typename Idx>
void SparseToDense(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& val, boost::multi_array<T, 4>& ret);

template <typename T, typename Idx>
void SparseToDense(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& val, boost::multi_array<T, 5>& ret);

template <typename T, typename Idx>
void SparseToDense(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& val, boost::multi_array<T, 6>& ret);

template <typename T, typename Idx>
void SparseMatToEigenSparse(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& val, size_t iRows, size_t iCols, Eigen::SparseMatrix<T>& ret);


//Gets a dictionary, coefficients (as location/value pair) and outputs reconstructed data point (dictionaries can be overcomplete)
template <typename T, typename Idx>
void HOSVDReconSp(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc, const Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vU, Eigen::Matrix<T, Eigen::Dynamic, 1>& coreTensor);

template <typename T, typename Idx>
void HOSVDReconSp(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc, const Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vU, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& coreTensor);

template <typename T, typename Idx>
void HOSVDReconSp(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc, const Eigen::Matrix<T, Eigen::Dynamic, 1>& val, 
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vU, boost::multi_array<T, 3>& coreTensor);

template <typename T, typename Idx>
void HOSVDReconSp(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc, const Eigen::Matrix<T, Eigen::Dynamic, 1>& val, 
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vU, boost::multi_array<T, 4>& coreTensor);

template <typename T, typename Idx>
void HOSVDReconSp(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc, const Eigen::Matrix<T, Eigen::Dynamic, 1>& val, 
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vU, boost::multi_array<T, 5>& coreTensor);

template <typename T, typename Idx>
void HOSVDReconSp(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc, const Eigen::Matrix<T, Eigen::Dynamic, 1>& val, 
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vU, boost::multi_array<T, 6>& coreTensor);

template <typename T>
void HOSVDReconR1(T singVal, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU, boost::multi_array<T, 3>& coreTensor);

template <typename T>
void HOSVDReconR1(T singVal, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU, boost::multi_array<T, 4>& coreTensor);

template <typename T>
void HOSVDReconR1(T singVal, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU, boost::multi_array<T, 5>& coreTensor);

template <typename T>
void HOSVDReconR1(T singVal, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU, boost::multi_array<T, 6>& coreTensor);



// template <typename T>
// void KronSelective_tmp1(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& matrices,
// 	const std::vector<size_t>& indices,
// 	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRes);
// 
// template <typename T>
// void KronSelective_tmp2(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& matrices,
// 	const std::vector<size_t>& indices,
// 	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRes);

#include "LATools.inl"
