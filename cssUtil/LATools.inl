

/* solve L*x = b */
template <typename T>
void backsubst_L(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& L, 
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& b, 
	Eigen::Matrix<T, Eigen::Dynamic, 1>& x, size_t k)
{
	size_t n = L.rows();
	if (k > L.rows() || k > L.cols())
	{
		x.resize(0);
		return;
	}
	x.resize(k);
	T rhs;

	for (size_t i = 0; i < k; ++i)
	{
		rhs = b[i];
		for (size_t j = 0; j < i; ++j)
			rhs -= L.data()[j * n + i] * x(j);
		x(i) = rhs / L.data()[i * n + i];
	}
}


/* solve L'*x = b */
template <typename T>
void backsubst_Lt(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& L,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& x, size_t k)
{
	size_t n = L.rows();
	if (k > L.rows() || k > L.cols())
	{
		x.resize(0);
		return;
	}
	x.resize(k);
	T rhs;

	for (size_t i = k; i >= 1; --i) 
	{
		rhs = b(i - 1);
		for (size_t j = i; j < k; ++j) 
			rhs -= L.data()[(i - 1) * n + j] * x(j);
		x(i - 1) = rhs / L.data()[(i - 1) * n + i - 1];
	}
}


/* solve U*x = b */
template <typename T>
void backsubst_U(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& U, 
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& x, size_t k)
{
	size_t n = U.rows();
	if (k > U.rows() || k > U.cols())
	{
		x.resize(0);
		return;
	}
	x.resize(k);
	T rhs;

	for (size_t i = k; i >= 1; --i) 
	{
		rhs = b(i - 1);
		for (size_t j = i; j < k; ++j) 
			rhs -= U.data()[j * n + i - 1] * x(j);
		x(i - 1) = rhs / U.data()[(i - 1) * n + i - 1];
	}
}


/* solve U'*x = b */
template <typename T>
void backsubst_Ut(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& U,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& x, size_t k)
{
	size_t n = U.rows();
	if (k > U.rows() || k > U.cols())
	{
		x.resize(0);
		return;
	}
	x.resize(k);
	T rhs;

	for (size_t i = 0; i < k; ++i) 
	{
		rhs = b(i);
		for (size_t j = 0; j < i; ++j) 
			rhs -= U.data()[i * n + j] * x(j);
		x(i) = rhs / U.data()[i * n + i];
	}
}


/* back substitution solver */
template <typename T>
void backsubst(char ul, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, 
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& b, 
	Eigen::Matrix<T, Eigen::Dynamic, 1>& x, size_t k)
{
	if (tolower(ul) == 'u') 
		backsubst_U(A, b, x, k);
	else if (tolower(ul) == 'l') 
		backsubst_L(A, b, x, k);
	else 
		std::cerr << "Invalid triangular matrix type for backsubst(): must be ''U'' or ''L''" << std::endl;
}


/* solve equation set using cholesky decomposition */
template <typename T>
void cholsolve(char ul, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& x, size_t k)
{
	Eigen::Matrix<T, Eigen::Dynamic, 1> tmp(k);

	if (tolower(ul) == 'l') 
	{
		backsubst_L(A, b, tmp, k);
		backsubst_Lt(A, tmp, x, k);
	}
	else if (tolower(ul) == 'u') 
	{
		backsubst_Ut(A, b, tmp, k);
		backsubst_U(A, tmp, x, k);
	}
	else 
	{
		std::cerr << "Invalid triangular matrix type for cholsolve(): must be ''U'' or ''L''" << std::endl;
	}
}



template <typename T>
T MutualCoherenceSum(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mLeft, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRight)
{
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mGramMat = mLeft.transpose() * mRight;
	return mGramMat.cwiseAbs().sum();
}


template <typename T>
T MutualCoherenceAvg(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mLeft, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRight)
{
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mGramMat = mLeft.transpose() * mRight;
	return mGramMat.cwiseAbs().sum() / (mLeft.cols() * (mLeft.cols() - 1));
}


template <typename T>
T MutualCoherenceExpSum(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mLeft, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRight)
{
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mGramMat = mLeft.transpose() * mRight;
	return exp(mGramMat.cwiseAbs().sum());
}


template <typename T>
T MutualCoherenceExpAvg(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mLeft, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRight)
{
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mGramMat = mLeft.transpose() * mRight;
	return exp(mGramMat.cwiseAbs().sum() / (mLeft.cols() * (mLeft.cols() - 1)));
}


template <typename T>
T MutualCoherenceExpProd(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mLeft, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRight)
{
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mGramMat = mLeft.transpose() * mRight;
	return exp(mGramMat.cwiseAbs().prod());
}


template <typename T>
T MutualCoherenceMax(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mLeft, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRight)
{
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mGramMat = mLeft.transpose() * mRight;
	for (size_t i = 0; i < mGramMat.rows(); ++i)
		mGramMat(i, i) = 0;
	return mGramMat.cwiseAbs().maxCoeff();
}


template <typename T>
T MutualCoherenceExpMax(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mLeft, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRight)
{
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mGramMat = mLeft.transpose() * mRight;
	for (size_t i = 0; i < mGramMat.rows(); ++i)
		mGramMat(i, i) = 0;
	return exp(mGramMat.cwiseAbs().maxCoeff());
}



template <typename Derived1, typename Derived2>
void Mat1Dto2D(const Eigen::DenseBase<Derived1>& mIn, size_t iRows, size_t iCols, Eigen::DenseBase<Derived2> const & mOut)
{
	if (mIn.rows() != (iRows*iCols))
	{
		Eigen::DenseBase<Derived2>& mOut1 = const_cast<Eigen::DenseBase<Derived2>&>(mOut);
		mOut1.derived().resize(0, 0);
		return;
	}
	Eigen::DenseBase<Derived2>& mOut1 = const_cast<Eigen::DenseBase<Derived2>&>(mOut);
	mOut1.derived().resize(iRows, iCols);
	for (size_t m = 0; m < iCols; ++m)
		for (size_t n = 0; n < iRows; ++n)
			mOut1(n, m) = mIn(n + m*iRows);
}


template <typename Derived1, typename Derived2>
void Mat2Dto1D(const Eigen::DenseBase<Derived1>& mIn, Eigen::DenseBase<Derived2> const & mOut)
{
	Eigen::DenseBase<Derived2>& mOut1 = const_cast<Eigen::DenseBase<Derived2>&>(mOut);
	mOut1.derived().resize(mIn.rows()*mIn.cols());
	for (size_t m = 0; m < mIn.cols(); ++m)
		for (size_t n = 0; n < mIn.rows(); ++n)
			mOut1(n + m*mIn.rows()) = mIn(n, m);
}


template <typename Derived1, typename Derived2>
CCS_INTERNAL_TYPE EuclidDist2(const Eigen::MatrixBase<Derived1>& coord1, const Eigen::MatrixBase<Derived2>& coord2)
{
	return (coord1 - coord2).squaredNorm();
}

template <size_t N>
CCS_INTERNAL_TYPE EuclidDist2(const boost::multi_array<CCS_INTERNAL_TYPE, N>& coord1, const boost::multi_array<CCS_INTERNAL_TYPE, N>& coord2)
{
	typedef typename boost::multi_array<CCS_INTERNAL_TYPE, N>::size_type szType;
	typedef typename Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1> Vec;
	uint32_t iLength = std::accumulate(coord1.shape(), coord1.shape() + N, szType(1), std::multiplies<szType>());
	return (Vec::Map(coord1.data(), iLength) - Vec::Map(coord2.data(), iLength)).squaredNorm();
}


template <typename Derived1, typename Derived2>
CCS_INTERNAL_TYPE EuclidDist1(const Eigen::MatrixBase<Derived1>& coord1, const Eigen::MatrixBase<Derived2>& coord2)
{
	return (coord1 - coord2).cwiseAbs().sum();
}

template <size_t N>
CCS_INTERNAL_TYPE EuclidDist1(const boost::multi_array<CCS_INTERNAL_TYPE, N>& coord1, const boost::multi_array<CCS_INTERNAL_TYPE, N>& coord2)
{
	typedef typename boost::multi_array<CCS_INTERNAL_TYPE, N>::size_type szType;
	typedef typename Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1> Vec;
	uint32_t iLength = std::accumulate(coord1.shape(), coord1.shape() + N, szType(1), std::multiplies<szType>());
	return (Vec::Map(coord1.data(), iLength) - Vec::Map(coord2.data(), iLength)).cwiseAbs().sum();
}


template <typename Derived1, typename Derived2>
CCS_INTERNAL_TYPE EuclidDistInf(const Eigen::MatrixBase<Derived1>& coord1, const Eigen::MatrixBase<Derived2>& coord2)
{
	return (coord1 - coord2).cwiseAbs().maxCoeff();
}

template <size_t N>
CCS_INTERNAL_TYPE EuclidDistInf(const boost::multi_array<CCS_INTERNAL_TYPE, N>& coord1, const boost::multi_array<CCS_INTERNAL_TYPE, N>& coord2)
{
	typedef typename boost::multi_array<CCS_INTERNAL_TYPE, N>::size_type szType;
	typedef typename Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1> Vec;
	uint32_t iLength = std::accumulate(coord1.shape(), coord1.shape() + N, szType(1), std::multiplies<szType>());
	return (Vec::Map(coord1.data(), iLength) - Vec::Map(coord2.data(), iLength)).cwiseAbs().maxCoeff();
}




/*
template <>
inline void SVD<float>(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& mMat,
	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& mU,
	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& mVt,
	Eigen::Matrix<float, Eigen::Dynamic, 1>& vS, 
	bool bThin)
{
	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> mat = mMat;
#ifdef SVD_CLAPACK
	integer m = mMat.rows();
	integer n = mMat.cols();
	integer lda = m;
	integer ldu = m;
	integer ldvt;
	integer INFO;
	char JOBU, JOBVT;
	if (bThin)
	{
		ldvt = std::min(m, n);
		JOBU = 'S';
		JOBVT = 'S';
		mU.resize(m, std::min(m, n));
		mVt.resize(std::min(m, n), n);
	}
	else
	{
		ldvt = n;
		JOBU = 'A';
		JOBVT = 'A';
		mU.resize(m, m);
		mVt.resize(n, n);
	}
	vS.resize(std::min(m, n));
	integer LWORK_Z1;
	float tmp;
	float* WORK_Z1 = NULL;
	LWORK_Z1 = -1;
	sgesvd_(&JOBU, &JOBVT, &m, &n, mat.data(), &lda, vS.data(), mU.data(), &ldu, mVt.data(), &ldvt, &tmp, &LWORK_Z1, &INFO);
	LWORK_Z1 = integer(tmp);
	WORK_Z1 = new float[LWORK_Z1];
	sgesvd_(&JOBU, &JOBVT, &m, &n, mat.data(), &lda, vS.data(), mU.data(), &ldu, mVt.data(), &ldvt, WORK_Z1, &LWORK_Z1, &INFO);
	SAFE_DELETE_ARRAY(WORK_Z1);
#elif SVD_LAPACKE
	lapack_int m = mMat.rows();
	lapack_int n = mMat.cols();
	lapack_int lda = m;
	lapack_int ldu = m;
	lapack_int ldvt;
	lapack_int info;
	char jobu, jobvt;
	if (bThin)
	{
		ldvt = std::min(m, n);
		jobu = 'S';
		jobvt = 'S';
		mU.resize(m, std::min(m, n));
		mVt.resize(std::min(m, n), n);
		vS.resize(std::min(m, n));
}
	else
	{
		ldvt = n;
		jobu = 'A';
		jobvt = 'A';
		mU.resize(m, m);
		mVt.resize(n, n);
		vS.resize(std::min(m, n));
	}
	std::vector<float> superb(std::min(m, n) - 1);
	info = LAPACKE_sgesvd(LAPACK_COL_MAJOR, jobu, jobvt, m, n, mat.data(), lda, vS.data(), mU.data(), ldu, mVt.data(), ldvt, &superb[0]);
	if (info > 0)
		printf("SVD failed.\n");
#elif SVD_MKL
	MKL_INT m = mMat.rows();
	MKL_INT n = mMat.cols();
	MKL_INT lda = m;
	MKL_INT ldu = m;
	MKL_INT ldvt;
	MKL_INT info;
	char jobu, jobvt;
	if (bThin)
	{
		ldvt = std::min(m, n);
		jobu = 'S';
		jobvt = 'S';
		mU.resize(m, std::min(m, n));
		mVt.resize(std::min(m, n), n);
		vS.resize(std::min(m, n));
	}
	else
	{
		ldvt = n;
		jobu = 'A';
		jobvt = 'A';
		mU.resize(m, m);
		mVt.resize(n, n);
		vS.resize(std::min(m, n));
	}
	std::vector<float> superb(std::min(m, n) - 1);
	info = LAPACKE_sgesvd(LAPACK_COL_MAJOR, jobu, jobvt, m, n, mat.data(), lda, vS.data(), mU.data(), ldu, mVt.data(), ldvt, &superb[0]);
	if (info > 0)
		printf("SVD failed.\n");
#else
#endif
}

template <>
inline void SVD<double>(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& mMat,
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& mU,
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& mVt,
	Eigen::Matrix<double, Eigen::Dynamic, 1>& vS,
	bool bThin)
{
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat = mMat;
#ifdef SVD_CLAPACK
	integer m = mMat.rows();
	integer n = mMat.cols();
	integer INFO;
	integer lda = m;
	integer ldu = m;
	integer ldvt;
	char JOBU, JOBVT;
	if (bThin)
	{
		ldvt = std::min(m, n);
		JOBU = 'S';
		JOBVT = 'S';
		mU.resize(m, std::min(m, n));
		mVt.resize(std::min(m, n), n);
}
	else
	{
		ldvt = n;
		JOBU = 'A';
		JOBVT = 'A';
		mU.resize(m, m);
		mVt.resize(n, n);
	}
	vS.resize(std::min(m, n));
	integer LWORK_Z1;
	double tmp;
	double* WORK_Z1 = NULL;
	LWORK_Z1 = -1;
	dgesvd_(&JOBU, &JOBVT, &m, &n, mat.data(), &lda, vS.data(), mU.data(), &ldu, mVt.data(), &ldvt, &tmp, &LWORK_Z1, &INFO);
	LWORK_Z1 = integer(tmp);
	WORK_Z1 = new double[LWORK_Z1];
	dgesvd_(&JOBU, &JOBVT, &m, &n, mat.data(), &lda, vS.data(), mU.data(), &ldu, mVt.data(), &ldvt, WORK_Z1, &LWORK_Z1, &INFO);
	SAFE_DELETE_ARRAY(WORK_Z1);
#elif SVD_LAPACKE
	lapack_int m = mMat.rows();
	lapack_int n = mMat.cols();
	lapack_int lda = m;
	lapack_int ldu = m;
	lapack_int ldvt;
	lapack_int info;
	char jobu, jobvt;
	if (bThin)
	{
		ldvt = std::min(m, n);
		jobu = 'S';
		jobvt = 'S';
		mU.resize(m, std::min(m, n));
		mVt.resize(std::min(m, n), n);
	}
	else
	{
		ldvt = n;
		jobu = 'A';
		jobvt = 'A';
		mU.resize(m, m);
		mVt.resize(n, n);
	}
	vS.resize(std::min(m, n));
	std::vector<double> superb(std::min(m, n) - 1);
	info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, jobu, jobvt, m, n, mat.data(), lda, vS.data(), mU.data(), ldu, mVt.data(), ldvt, &superb[0]);
	if (info > 0)
		printf("SVD failed.\n");
#elif SVD_MKL
	MKL_INT m = mMat.rows();
	MKL_INT n = mMat.cols();
	MKL_INT lda = m;
	MKL_INT ldu = m;
	MKL_INT ldvt;
	MKL_INT info;
	char jobu, jobvt;
	if (bThin)
	{
		ldvt = std::min(m, n);
		jobu = 'S';
		jobvt = 'S';
		mU.resize(m, std::min(m, n));
		mVt.resize(std::min(m, n), n);
	}
	else
	{
		ldvt = n;
		jobu = 'A';
		jobvt = 'A';
		mU.resize(m, m);
		mVt.resize(n, n);
	}
	vS.resize(std::min(m, n));
	std::vector<double> superb(std::min(m, n) - 1);
	info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, jobu, jobvt, m, n, mat.data(), lda, vS.data(), mU.data(), ldu, mVt.data(), ldvt, &superb[0]);
	if (info > 0)
		printf("SVD failed.\n");
#else
#endif
}
*/



template<typename T, size_t N>
void HOSVD(const boost::multi_array<T, N>& inputTensor, std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vU, boost::multi_array<T, N>& coreTensor)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	std::vector<szType> vShape;
	vShape.assign(inputTensor.shape(), inputTensor.shape() + inputTensor.num_dimensions());
	szType iLength = std::accumulate(vShape.begin(), vShape.end(), 1, std::multiplies<szType>());

	vU.resize(vShape.size());
	for (size_t i = 0; i < vShape.size(); ++i)
		vU[i].resize(vShape[i], vShape[i]);

	//Compute bases U_i
	for (size_t i = 0; i < vShape.size(); ++i)
	{
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mUnfolded;
		TensorUnfold<T>(inputTensor, i, mUnfolded);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mU, mV;
		Eigen::Matrix<T, Eigen::Dynamic, 1> vS;
		Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > svd(mUnfolded, Eigen::ComputeThinU | Eigen::ComputeThinV);
		mU = svd.matrixU();
		vS = svd.singularValues();
		vU[i] = mU;
	}

	//Compute core tensor (coefficients)
	boost::multi_array<T, N> tmp(vShape, CCS_TENSOR_STORAGE_ORDER);
	TensorCopy<T, N>(inputTensor, tmp);
	for (size_t i = 0; i < vShape.size(); ++i)
	{
		TensorProdMat(tmp, vU[i].transpose(), i, coreTensor);
		TensorCopy(coreTensor, tmp);
	}
}




template <typename T>
void TensorProdALS(const boost::multi_array<T, 3>& inputTensor, size_t iMode, std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU)
{
	if (inputTensor.num_dimensions() != vU.size() || iMode >= vU.size())
		return;
	typedef typename boost::multi_array<T, 3>::size_type szType;
	typedef typename boost::multi_array<T, 3>::index idxType;

	const szType* shape = inputTensor.shape();
	std::vector<size_t> vOtherModes;
	vOtherModes.reserve(2);
	for (szType i = 0; i < inputTensor.num_dimensions(); ++i)
		if (i != iMode)
			vOtherModes.push_back(i);

	for (idxType i = 0; i < shape[iMode]; ++i)
	{
		T accum = 0;
		for (idxType i2 = 0; i2 < shape[vOtherModes[1]]; ++i2)
		{
			for (idxType i1 = 0; i1 < shape[vOtherModes[0]]; ++i1)
			{
				boost::array<idxType, 3> idx;
				idx[vOtherModes[0]] = i1;
				idx[vOtherModes[1]] = i2;
				idx[iMode] = i;
				accum += inputTensor(idx) * (vU[vOtherModes[0]])(i1) * (vU[vOtherModes[1]])(i2);
			}
		}
		(vU[iMode])(i) = accum;
	}
	vU[iMode] /= vU[iMode].norm();
}

template <typename T>
void TensorProdALS(const boost::multi_array<T, 4>& inputTensor, size_t iMode, std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU)
{
	if (inputTensor.num_dimensions() != vU.size() || iMode >= vU.size())
		return;
	typedef typename boost::multi_array<T, 4>::size_type szType;
	typedef typename boost::multi_array<T, 4>::index idxType;

	const szType* shape = inputTensor.shape();
	std::vector<size_t> vOtherModes;
	vOtherModes.reserve(3);
	for (szType i = 0; i < inputTensor.num_dimensions(); ++i)
		if (i != iMode)
			vOtherModes.push_back(i);

	for (idxType i = 0; i < shape[iMode]; ++i)
	{
		T accum = 0;
		for (idxType i3 = 0; i3 < shape[vOtherModes[2]]; ++i3)
		{
			for (idxType i2 = 0; i2 < shape[vOtherModes[1]]; ++i2)
			{
				for (idxType i1 = 0; i1 < shape[vOtherModes[0]]; ++i1)
				{
					boost::array<idxType, 4> idx;
					idx[vOtherModes[0]] = i1;
					idx[vOtherModes[1]] = i2;
					idx[vOtherModes[2]] = i3;
					idx[iMode] = i;
					accum += inputTensor(idx) * (vU[vOtherModes[0]])(i1) * (vU[vOtherModes[1]])(i2) * (vU[vOtherModes[2]])(i3);
				}
			}
		}
		(vU[iMode])(i) = accum;
	}
	vU[iMode] /= vU[iMode].norm();
}

template <typename T>
void TensorProdALS(const boost::multi_array<T, 5>& inputTensor, size_t iMode, std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU)
{
	if (inputTensor.num_dimensions() != vU.size() || iMode >= vU.size())
		return;
	typedef typename boost::multi_array<T, 5>::size_type szType;
	typedef typename boost::multi_array<T, 5>::index idxType;

	const szType* shape = inputTensor.shape();
	std::vector<size_t> vOtherModes;
	vOtherModes.reserve(4);
	for (szType i = 0; i < inputTensor.num_dimensions(); ++i)
		if (i != iMode)
			vOtherModes.push_back(i);

	for (idxType i = 0; i < shape[iMode]; ++i)
	{
		T accum = 0;
		for (idxType i4 = 0; i4 < shape[vOtherModes[3]]; ++i4)
		{
			for (idxType i3 = 0; i3 < shape[vOtherModes[2]]; ++i3)
			{
				for (idxType i2 = 0; i2 < shape[vOtherModes[1]]; ++i2)
				{
					for (idxType i1 = 0; i1 < shape[vOtherModes[0]]; ++i1)
					{
						boost::array<idxType, 5> idx;
						idx[vOtherModes[0]] = i1;
						idx[vOtherModes[1]] = i2;
						idx[vOtherModes[2]] = i3;
						idx[vOtherModes[3]] = i4;
						idx[iMode] = i;
						accum += inputTensor(idx) * (vU[vOtherModes[0]])(i1) * (vU[vOtherModes[1]])(i2) * (vU[vOtherModes[2]])(i3) * (vU[vOtherModes[3]])(i4);
					}
				}
			}
		}
		(vU[iMode])(i) = accum;
	}
	vU[iMode] /= vU[iMode].norm();
}

template <typename T>
void TensorProdALS(const boost::multi_array<T, 6>& inputTensor, size_t iMode, std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU)
{
	if (inputTensor.num_dimensions() != vU.size() || iMode >= vU.size())
		return;
	typedef typename boost::multi_array<T, 6>::size_type szType;
	typedef typename boost::multi_array<T, 6>::index idxType;

	const szType* shape = inputTensor.shape();
	std::vector<size_t> vOtherModes;
	vOtherModes.reserve(5);
	for (szType i = 0; i < inputTensor.num_dimensions(); ++i)
		if (i != iMode)
			vOtherModes.push_back(i);

	for (idxType i = 0; i < shape[iMode]; ++i)
	{
		T accum = 0;
		for (idxType i5 = 0; i5 < shape[vOtherModes[4]]; ++i5)
		{
			for (idxType i4 = 0; i4 < shape[vOtherModes[3]]; ++i4)
			{
				for (idxType i3 = 0; i3 < shape[vOtherModes[2]]; ++i3)
				{
					for (idxType i2 = 0; i2 < shape[vOtherModes[1]]; ++i2)
					{
						for (idxType i1 = 0; i1 < shape[vOtherModes[0]]; ++i1)
						{
							boost::array<idxType, 6> idx;
							idx[vOtherModes[0]] = i1;
							idx[vOtherModes[1]] = i2;
							idx[vOtherModes[2]] = i3;
							idx[vOtherModes[3]] = i4;
							idx[vOtherModes[4]] = i5;
							idx[iMode] = i;
							accum += inputTensor(idx) * (vU[vOtherModes[0]])(i1) * (vU[vOtherModes[1]])(i2) * (vU[vOtherModes[2]])(i3) * (vU[vOtherModes[3]])(i4) * (vU[vOtherModes[4]])(i5);
						}
					}
				}
			}
		}
		(vU[iMode])(i) = accum;
	}
	vU[iMode] /= vU[iMode].norm();
}

template <typename T, size_t N>
int TensorALSR1(const boost::multi_array<T, N>& inputTensor, size_t iMaxIter, T thresh, std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU)
{
	if (inputTensor.num_dimensions() != vU.size())
		return -1;
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;
	std::vector<szType> vShape;
	vShape.assign(inputTensor.shape(), inputTensor.shape() + inputTensor.num_dimensions());

	std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> > vLastU = vU;

	size_t i = 0;
	for (i = 0; i < iMaxIter; ++i)
	{
		for (size_t j = 0; j < vShape.size(); ++j)
			TensorProdALS(inputTensor, j, vU);
		//Calculate error
		T err = 0;
		for (size_t j = 0; j < vU.size(); ++j)
			err += (vLastU[j] - vU[j]).squaredNorm();
		if (err < thresh)
			break;
		for (size_t j = 0; j < vU.size(); ++j)
			vLastU[j] = vU[j];
	}
	return i;
}



// template <typename T, size_t N>
// void HOSVDtruncated(const boost::multi_array<T, N>& inputTensor,
// 	std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vU,
// 	boost::multi_array<T, N>& coreTensor)
// {
// 	typedef typename boost::multi_array<T, N>::size_type szType;
// 	typedef typename boost::multi_array<T, N>::index idxType;
// 
// 	std::vector<szType> vShape;
// 	vShape.assign(inputTensor.shape(), inputTensor.shape() + inputTensor.num_dimensions());
// 	szType iLength = std::accumulate(vShape.begin(), vShape.end(), 1, std::multiplies<szType>());
// 	vU.resize(vShape.size());
// 
// 	//Compute bases U_i
// 	for (size_t i = 0; i < vShape.size(); ++i)
// 	{
// 		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mUnfolded;
// 		TensorUnfold<T>(inputTensor, i, mUnfolded);
//		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mU, mVt;
//		Eigen::Matrix<T, Eigen::Dynamic, 1> vS;
//		SVD(mUnfolded, mU, mVt, vS, true);
// 		vU[i] = mU;
// 	}
// 
// 	//Compute core tensor (coefficients)
// 	boost::multi_array<T, N> tmp;
// 	TensorCopy<T, N>(inputTensor, tmp);
// 	for (size_t i = 0; i < vShape.size(); ++i)
// 	{
// 		TensorProdMat(tmp, vU[i].transpose(), i, coreTensor);
// 		TensorCopy(coreTensor, tmp);
// 	}
// }


// 
// template<typename T>
// void PowerIter(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mMat,
// 	double& sigma,
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1>& vU,
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1>& vV,
// 	CCS_INTERNAL_TYPE thresh,
// 	size_t iMaxIter)
// {
// 	thresh *= std::max(mMat.rows(), mMat.cols());
// 
// 	vU.setZero(mMat.rows());
// 	vV.resize(mMat.cols());
// 	for (size_t i = 0; i < vV.size(); ++i)
// 		vV(i) = RandUni();
// 	vV /= vV.norm();
// 
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1> vU1(vU.size());
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1> vV1(vV.size());
// 	vU1 = vU;
// 	vV1 = vV;
// 	for (size_t i = 0; i < iMaxIter; ++i)
// 	{
// #pragma omp parallel for schedule(guided)
// 		for(size_t j = 0; j < mMat.rows(); ++j)
// 			vU(j) = mMat.row(j).template cast<CCS_INTERNAL_TYPE>() * vV;
// 		//vU = mMat * vV;
// 		vU /= vU.norm();
// #pragma omp parallel for schedule(guided)
// 		for (size_t j = 0; j < mMat.cols(); ++j)
// 			vV(j) = mMat.col(j).template cast<CCS_INTERNAL_TYPE>().dot(vU);
// 		//vV = mMat.transpose() * vU;
// 		vV /= vV.norm();
// 		if ((vV - vV1).norm() < thresh)
// 		{
// 			break;
// 		}
// 		else
// 		{
// 			vV1 = vV;
// 			vU1 = vU;
// 		}
// 	}
// 
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1> tmp(mMat.cols());
// #pragma omp parallel for schedule(guided)
// 	for (size_t j = 0; j < mMat.cols(); ++j)
// 		tmp(j) = vU.dot(mMat.col(j).template cast<CCS_INTERNAL_TYPE>());
// 	sigma = tmp.dot(vV);
// 	//sigma = vU.transpose() * mMat * vV;
// }
// 
// template<typename T>
// void PowerIterFull(uint32_t iNumTerms,
// 	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& F, 
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, Eigen::Dynamic>& U,
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, Eigen::Dynamic>& V,
// 	CCS_INTERNAL_TYPE threshold, 
// 	uint32_t iMaxIter, 
// 	double& dTimeDelta)
// {
// 	double dTime = omp_get_wtime();
// 	int iNodeID = 0;
// 	//MPI_Comm_rank(MPI_COMM_WORLD, &iNodeID);
// 
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, Eigen::Dynamic> Fd = F.template cast<CCS_INTERNAL_TYPE>();
// 	CProgressReporter reporter(iNumTerms, iNodeID);
// 	for (uint32_t iTerm = 0; iTerm < iNumTerms; ++iTerm)
// 	{
// 		Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1> vU;
// 		Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1> vV;
// 		CCS_INTERNAL_TYPE sigma = 0;
// 		PowerIter(Fd, sigma, vU, vV, threshold, iMaxIter);
// #pragma omp parallel for schedule(guided)
// 		for (int i = 0; i < vV.size(); ++i)
// 			Fd.col(i) -= (sigma*vV(i)) * vU;
// //		F -= (sigma*vU) * vV.transpose();
// 		U.col(iTerm) = vU;
// 		V.col(iTerm) = vV;
// 		reporter.Update();
// 	}
// 	reporter.Done();
// 	dTimeDelta = omp_get_wtime() - dTime;
// }
// 
// 
// template <typename T>
// void PowerIterFull(uint32_t iNumTerms,
// 	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& F,
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, Eigen::Dynamic>& U,
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, Eigen::Dynamic>& V,
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1>& S,
// 	CCS_INTERNAL_TYPE threshold, uint32_t iMaxIter, double& dTimeDelta)
// {
// 	double dTime = omp_get_wtime();
// 	int iNodeID = 0;
// 	MPI_Comm_rank(MPI_COMM_WORLD, &iNodeID);
// 
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, Eigen::Dynamic> Fd = F.template cast<CCS_INTERNAL_TYPE>();
// 	CProgressReporter reporter(iNumTerms, iNodeID);
// 	for (uint32_t iTerm = 0; iTerm < iNumTerms; ++iTerm)
// 	{
// 		Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1> vU;
// 		Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1> vV;
// 		CCS_INTERNAL_TYPE sigma = 0;
// 		PowerIter(Fd, sigma, vU, vV, threshold, iMaxIter);
// #pragma omp parallel for schedule(static)
// 		for (size_t i = 0; i < vV.size(); ++i)
// 			Fd.col(i) -= (sigma*vV(i)) * vU;
// 		//F -= (sigma*vU) * vV.transpose();
// 		U.col(iTerm) = vU;
// 		V.col(iTerm) = vV;
// 		S(iTerm) = sigma;
// 		reporter.Update();
// 	}
// 	reporter.Done();
// 	dTimeDelta = omp_get_wtime() - dTime;
// }
// 
// 
// template <typename T>
// void PowerIterFullOld(uint32_t iNumTerms,
// 	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& F,
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, Eigen::Dynamic>& U,
// 	Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, Eigen::Dynamic>& V,
// 	CCS_INTERNAL_TYPE threshold,
// 	uint32_t iMaxIter, 
// 	double& dTimeDelta)
// {
// 	double dLastTime = omp_get_wtime();
// 	uint32_t M = F.rows();//m_pInitInfo->m_iRows * m_pInitInfo->m_iColumns;
// 	uint32_t N = F.cols();//m_pInitInfo->m_iNumSamplesX * m_pInitInfo->m_iNumSamplesY;
// 	uint32_t K = iNumTerms;
// 
// 	int iNodeID = 0;
// 	//MPI_Comm_rank(MPI_COMM_WORLD, &iNodeID);
// 	CProgressReporter reporter(K, iNodeID);
// 
// 	Eigen::Matrix<T, Eigen::Dynamic, 1> UU;
// 	Eigen::Matrix<T, Eigen::Dynamic, 1> VV;
// 	Eigen::Matrix<T, Eigen::Dynamic, 1> VV2;
// 	Eigen::Matrix<T, Eigen::Dynamic, 1> UU2;
// 	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Ap;
// 	T lambda_p, maxu, maxv, alpha, diff, lastDiff;
// 	lastDiff = std::numeric_limits<T>::max();
// 
// 	// 	if(M>=N)
// 	// 	{
// 	UU.resize(M);
// 	VV.resize(N);
// 	VV2.resize(N);
// 	Ap.resize(N, N);
// 
// 	for (uint32_t p = 1; p <= K; ++p)
// 	{
// 		Ap = F.transpose() * F;
// 		for (uint32_t i = 0; i < VV2.rows(); ++i)
// 			VV2(i) = RandUni();
// 		VV2 = VV2.norm() * VV2;
// 
// 		uint32_t iNumItr;
// 		for (iNumItr = 0; iNumItr < iMaxIter; ++iNumItr)
// 		{
// 			VV = Ap * VV2;
// 			lambda_p = VV.norm();
// 			if (lambda_p == 0.0)
// 			{
// 				lambda_p = 1.0;
// 				T fTmp = 1.0 / sqrt((T)(N));
// 				for (uint32_t i = 0; i < N; ++i)
// 					VV(i, 0) = fTmp;
// 				break;
// 			}
// 			else
// 			{
// 				VV = VV / lambda_p;
// 				diff = (VV - VV2).norm();
// 				if ((diff <= threshold) || ((iNumItr >= 1) && (diff > lastDiff)))
// 					break;
// 				lastDiff = diff;
// 				VV2 = VV;
// 			}
// 		}
// 
// 		lambda_p = sqrt(lambda_p);
// 		UU = F * (VV / lambda_p);
// 		maxu = (UU.array().abs()).matrix().maxCoeff();
// 		maxv = (VV.array().abs()).matrix().maxCoeff();
// 		alpha = sqrt(maxu * lambda_p / maxv);
// 		UU = (lambda_p * UU) / alpha;
// 		VV = alpha * VV;
// 		U.col(p - 1) = UU;
// 		V.col(p - 1) = VV;
// #pragma omp parallel for schedule(static)
// 		for (int i = 0; i < VV.size(); ++i)
// 			F.col(i) -= VV(i) * UU;
// 		// #pragma omp parallel for schedule(static)
// 		// 		for (int i = 0; i < F.rows(); ++i)
// 		// 			for (int j = 0; j < F.cols(); ++j)
// 		// 				F(i, j) -= UU(i) * VV(j);// (*F) -= (UU * VV.transpose());
// 		reporter.Update();
// 	}
// 	reporter.Done();
// 	dTimeDelta = omp_get_wtime() - dLastTime;
// 	// 	}
// 	// 	else
// 	// 	{
// 	// 		UU.resize(M);
// 	// 		VV.resize(N);
// 	// 		UU2.resize(M);
// 	// 		Ap.resize(M,M);
// 	// 		for(uint32_t p = 1; p <= K; ++p)
// 	// 		{
// 	// 			Ap =  (*F)*F->transpose();
// 	// 			boost::minstd_rand generator(42u);
// 	// 			boost::uniform_real<> uni_dist(0,1);
// 	// 			boost::variate_generator<boost::minstd_rand&, boost::uniform_real<> > uni(generator, uni_dist);
// 	// 			for(uint32_t i = 0; i < UU2.rows(); ++i)
// 	// 				UU2(i) = uni();
// 	// 			UU2 = UU2.norm() * UU2;
// 	// 			uint32_t iNumItr;
// 	// 			for(iNumItr = 0; iNumItr < SLF_PCA_MAX_ITER; ++iNumItr)
// 	// 			{
// 	// 				UU = Ap * UU2;
// 	// 				lambda_p = UU.norm();
// 	// 				if(lambda_p == 0.0f)
// 	// 				{
// 	// 					lambda_p = 1.0f;
// 	// 					float fTmp = 1.0f / sqrtf((float)(M));
// 	// 					for(uint32_t i = 0; i < M; ++i)
// 	// 						UU(i,0) = fTmp;
// 	// 					break;
// 	// 				}
// 	// 				else
// 	// 				{
// 	// 					UU = UU / lambda_p;
// 	// 					diff = (UU - UU2).norm();
// 	// 					if( (diff <= SLF_PCA_THRESHOLD) || ((iNumItr>=1) && (diff>lastDiff)) )
// 	// 						break;
// 	// 					lastDiff = diff;
// 	// 					UU2 = UU;
// 	// 				}
// 	// 			}
// 	// 
// 	// 			lambda_p = sqrtf(lambda_p);
// 	// 			VV =  ((UU.transpose()*(*F))/ lambda_p).transpose();
// 	// 			maxu = (UU.array().abs()).matrix().maxCoeff();
// 	// 			maxv = (VV.array().abs()).matrix().maxCoeff();
// 	// 			alpha = sqrtf(maxu * lambda_p / maxv);
// 	// 			VV = (lambda_p * VV) / alpha;
// 	// 			UU = alpha * UU;
// 	// 			U->col(p-1) = UU;
// 	// 			V->col(p-1) = VV;
// 	// 			(*F) -= (UU * VV.transpose());
// 	// 			pReporter->Update();
// 	// 		}
// 	// 	}
// }


////Pseudo-inverse
//template<typename _Matrix_Type_>
//void PInv(const _Matrix_Type_& mIn, _Matrix_Type_& mOut, double epsilon)
//{
//	Eigen::JacobiSVD< _Matrix_Type_ > svd(mIn, Eigen::ComputeThinU | Eigen::ComputeThinV);
//	double tolerance = epsilon * std::max(mIn.cols(), mIn.rows()) *svd.singularValues().array().abs()(0);
//	mOut = svd.matrixV() * (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
//}


//Kronecker product of matrices
template <typename Derived1, typename Derived2, typename Derived3>
void Kron(const Eigen::MatrixBase<Derived1>& mLeft,
	const Eigen::MatrixBase<Derived2>& mRight,
	Eigen::MatrixBase<Derived3> const & mResC)
{
	Eigen::MatrixBase<Derived3>& mRes = const_cast<Eigen::MatrixBase<Derived3>&>(mResC);
	if ((mLeft.rows()*mRight.rows() != mResC.rows()) || (mLeft.cols()*mRight.cols() != mResC.cols()))
		mRes.derived().resize(mLeft.rows()*mRight.rows(), mLeft.cols()*mRight.cols());

	for (size_t j = 0; j < mLeft.rows(); j++)
		for (size_t k = 0; k < mLeft.cols(); k++)
			mRes.block(j*mRight.rows(), k*mRight.cols(), mRight.rows(), mRight.cols()) = mLeft(j, k)*mRight;
}


template <typename T>
void KronSelective(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& matrices,
	const std::vector<size_t>& indices,
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRes)
{
	if (matrices.size() < 2)
	{
		mRes.resize(0, 0);
		return;
	}

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp1;
	tmp1 = matrices[indices[0]];

	for (size_t i = 0; (i + 1) < indices.size(); i++)
	{
		Kron(tmp1, matrices[indices[i + 1]], mRes);
		tmp1 = mRes;
	}
}


template <typename T>
double KronSelectiveElem(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& matrices,
	const std::vector<size_t>& indices, size_t iRow, size_t iCol)
{
	std::vector<size_t> vRowMult(indices.size() + 1);
	std::vector<size_t> vColMult(indices.size() + 1);
	for (size_t i = 0; i < indices.size(); ++i)
	{
		vRowMult[i] = matrices[indices[i]].rows();
		vColMult[i] = matrices[indices[i]].cols();
	}
	vRowMult.back() = 1;
	vColMult.back() = 1;
	for (int i = vRowMult.size() - 2; i >= 0; i--)
		vRowMult[i] *= vRowMult[i + 1];
	for (int i = vColMult.size() - 2; i >= 0; i--)
		vColMult[i] *= vColMult[i + 1];

	double res = 1.0;
	for (size_t p = 0; p < indices.size(); ++p)
	{
		size_t i_p = (iRow / vRowMult[p + 1]) % matrices[indices[p]].rows();
		size_t j_p = (iCol / vColMult[p + 1]) % matrices[indices[p]].cols();
		res *= matrices[indices[p]](i_p, j_p);
	}
	return res;
}


template <typename T>
double KronSelectiveElem(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& matrices,
	const std::vector<size_t>& indices, size_t iRow, size_t iCol, const std::vector<size_t>& vRowMult, const std::vector<size_t>& vColMult)
{
	double res = 1.0;
	for (size_t p = 0; p < indices.size(); ++p)
	{
		size_t i_p = (iRow / vRowMult[p + 1]) % matrices[indices[p]].rows();
		size_t j_p = (iCol / vColMult[p + 1]) % matrices[indices[p]].cols();
		res *= matrices[indices[p]](i_p, j_p);
	}
	return res;
}


template <typename T, typename Idx>
void MultKronSelBySpMatTran(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& matrices,
	const std::vector<size_t>& indices,
	Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	size_t iNumColsForResult,
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& ret)
{
	size_t iNumRowsForResult = 1;
	for (size_t i = 0; i < indices.size(); ++i)
		iNumRowsForResult *= matrices[indices[i]].rows();
	if ((ret.rows() != iNumRowsForResult) || (ret.cols() != iNumColsForResult))
		ret.resize(iNumRowsForResult, iNumColsForResult);

	ret.setZero();

	for (size_t i = 0; i < loc.cols(); ++i)
	{
		Idx iRow = loc(1, i);	//Switched because we take the transpose of the sparse matrix
		Idx iCol = loc(0, i);
		for (size_t j = 0; j < iNumRowsForResult; ++j)
			ret(j, iCol) += KronSelectiveElem(matrices, indices, j, iRow) * val(i);
	}
}


template <typename T, typename Idx>
void MultKronSelBySpMatTran(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& matrices,
	const std::vector<size_t>& indices,
	Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	const std::vector<size_t>& vRowMult,
	const std::vector<size_t>& vColMult,
	size_t iNumColsForResult,
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& ret)
{
	size_t iNumRowsForResult = 1;
	for (size_t i = 0; i < indices.size(); ++i)
		iNumRowsForResult *= matrices[indices[i]].rows();
	if ((ret.rows() != iNumRowsForResult) || (ret.cols() != iNumColsForResult))
		ret.resize(iNumRowsForResult, iNumColsForResult);

	ret.setZero();

	for (size_t i = 0; i < loc.cols(); ++i)
	{
		Idx iRow = loc(1, i);	//Switched because we take the transpose of the sparse matrix
		Idx iCol = loc(0, i);
		for (size_t j = 0; j < iNumRowsForResult; ++j)
			ret(j, iCol) += KronSelectiveElem(matrices, indices, j, iRow, vRowMult, vColMult) * val(i);
	}
}


template <typename T>
void KhatriRaoCol(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mLeft,
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRight,
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRes)
{
	if (mLeft.cols() != mRight.cols())
	{
		mRes.resize(0, 0);
		return;
	}
	if (mLeft.rows() * mRight.rows() != mRes.rows())
		mRes.resize(mLeft.rows() * mRight.rows(), mLeft.cols());

	for (size_t i = 0; i < mRes.cols(); ++i)
		for (size_t j = 0; j < mLeft.rows(); ++j)
			mRes.col(i).segment(j*mRight.rows(), mRight.rows()) = mLeft(j, i) * mRight.col(i);
}

template <typename T>
void KhatriRaoColSelective(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& matrices,
	const std::vector<size_t>& indices,
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRes)
{
	if (matrices.size() < 2)
	{
		mRes.resize(0, 0);
		return;
	}

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp1;
	tmp1 = matrices[indices[0]];

	for (size_t i = 0; (i + 1) < indices.size(); i++)
	{
		KhatriRaoCol(tmp1, matrices[indices[i + 1]], mRes);
		tmp1 = mRes;
	}
}


template <typename Idx>
void SpMatTrans(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& nzLoc, Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& result)
{
	if ((nzLoc.rows() != result.rows()) || (nzLoc.cols() != result.cols()))
		result.resize(nzLoc.rows(), nzLoc.cols());
	result.row(0) = nzLoc.row(1);
	result.row(1) = nzLoc.row(0);
}

template <typename Idx>
void SpMatTransInplace(Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& nzLoc)
{
	Eigen::Matrix<Idx, Eigen::Dynamic, 1> tmp = nzLoc.row(0);
	nzLoc.row(0) = nzLoc.row(1);
	nzLoc.row(1) = tmp;
}


//Clamp function for 1D data (in-place)
template <typename T>
void Clamp1D(Eigen::Matrix<T, Eigen::Dynamic, 1>& mIn, T minVal, T maxVal)
{
	for (size_t j = 0; j < mIn.size(); ++j)
	{
		if (mIn(j) > maxVal)
			mIn(j) = maxVal;
		if (mIn(j) < minVal)
			mIn(j) = minVal;
	}
}


//Clamp function for 2D data (in-place)
template <typename T>
void Clamp2D(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mIn, T minVal, T maxVal)
{
	for (size_t j = 0; j < mIn.cols(); ++j)
	{
		for (size_t k = 0; k < mIn.rows(); ++k)
		{
			if (mIn(k, j) > maxVal)
				mIn(k, j) = maxVal;
			if (mIn(k, j) < minVal)
				mIn(k, j) = minVal;
		}
	}
}


template <typename T, size_t N>
void ClampnD(boost::multi_array<T, N>& mIn, T minVal, T maxVal)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;
	std::vector<szType> shape;

	shape.assign(mIn.shape(), mIn.shape() + mIn.num_dimensions());
	szType iLength = std::accumulate(shape.begin(), shape.end(), szType(1), std::multiplies<szType>());
	for (idxType j = 0; j < iLength; ++j)
	{
		if (mIn.data()[j] > maxVal)
			mIn.data()[j] = maxVal;
		if (mIn.data()[j] < minVal)
			mIn.data()[j] = minVal;
	}
}


template <typename inType, typename dimType, typename outType>
void Ind2Sub(inType idx, const std::vector<dimType>& vDims, std::vector<outType>& vIdx)
{
	vIdx.resize(vDims.size());
	std::vector<dimType> vK(vDims.size());
	std::partial_sum(vDims.begin(), vDims.end(), vK.begin(), std::multiplies<dimType>());
	for (size_t i = vDims.size() - 1; i > 0; i--)
	{
		size_t vi = idx % vK[i - 1];
		size_t vj = (idx - vi) / vK[i - 1];
		vIdx[i] = vj;
		idx = vi;
	}
	vIdx[0] = idx % vK[0];
}


template <typename inType, typename dimType, typename outType>
void Ind2Sub(inType idx, const Eigen::Matrix<dimType, Eigen::Dynamic, 1>& vDims, Eigen::Matrix<outType, Eigen::Dynamic, 1>& vIdx)
{
	vIdx.resize(vDims.size());
	Eigen::Matrix<dimType, Eigen::Dynamic, 1> vK(vDims.size());
	std::partial_sum(vDims.data(), vDims.data() + vDims.size(), vK.data(), std::multiplies<dimType>());
	for (size_t i = vDims.size() - 1; i > 0; i--)
	{
		size_t vi = idx % vK(i - 1);
		size_t vj = (idx - vi) / vK(i - 1);
		vIdx(i) = vj;
		idx = vi;
	}
	vIdx(0) = idx % vK(0);
}


template <typename inType, typename dimType, typename outType>
void Ind2Sub(inType idx, const std::vector<dimType>& vDims, Eigen::Matrix<outType, Eigen::Dynamic, 1>& vIdx)
{
	vIdx.resize(vDims.size());

	std::vector<dimType> vK(vDims.size());
	std::partial_sum(vDims.begin(), vDims.end(), vK.begin(), std::multiplies<dimType>());

	for (size_t i = vDims.size() - 1; i > 0; i--)
	{
		size_t vi = idx % vK[i - 1];
		size_t vj = (idx - vi) / vK[i - 1];
		vIdx(i) = vj;
		idx = vi;
	}
	vIdx(0) = idx % vK[0];
}


template <typename inType, typename dimType>
size_t Sub2Ind(const inType* vIdx, const dimType* vDims, size_t iNumDims)
{
	std::vector<dimType> vK(iNumDims);
	std::partial_sum(vDims, vDims + iNumDims, vK.begin(), std::multiplies<dimType>());
	size_t idx = vIdx[0];
	for (size_t i = 1; i < iNumDims; ++i)
		idx += size_t(vIdx[i]) * vK[i - 1];
	return idx;
}


template <typename T, size_t N>
void TensorSetOne(boost::multi_array<T, N>& tensor)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	const szType* shape = tensor.shape();
	szType iLength = std::accumulate(shape, shape + tensor.num_dimensions(), szType(1), std::multiplies<szType>());
	for (idxType i = 0; i < iLength; ++i)
		tensor.data()[i] = T(1);
}


template <typename T, size_t N>
void TensorSetZero(boost::multi_array<T, N>& tensor)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	const szType* shape = tensor.shape();
	szType iLength = std::accumulate(shape, shape + tensor.num_dimensions(), szType(1), std::multiplies<szType>());
	for (idxType i = 0; i < iLength; ++i)
		tensor.data()[i] = T(0);
}


template <typename T, size_t N>
void TensorSetConst(boost::multi_array<T, N>& tensor, T val)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	const szType* shape = tensor.shape();
	szType iLength = std::accumulate(shape, shape + tensor.num_dimensions(), szType(1), std::multiplies<szType>());
	for (idxType i = 0; i < iLength; ++i)
		tensor.data()[i] = val;
}

template<typename T, size_t N>
void TensorSetRandUni(boost::multi_array<T, N>& tensor, T tMin, T tMax)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	const szType* shape = tensor.shape();
	szType iLength = std::accumulate(shape, shape + tensor.num_dimensions(), szType(1), std::multiplies<szType>());
	for (idxType i = 0; i < iLength; ++i)
		tensor.data()[i] = RandUni(tMin, tMax);
}

template <typename T, size_t N>
void TensorSetRandGaussian(boost::multi_array<T, N>& tensor, T mean, T var)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	const szType* shape = tensor.shape();
	szType iLength = std::accumulate(shape, shape + tensor.num_dimensions(), szType(1), std::multiplies<szType>());
	for (idxType i = 0; i < iLength; ++i)
		tensor.data()[i] = RandGaussian(mean, var);
}

template <typename T, size_t N>
void TensorSetRandNormal(boost::multi_array<T, N>& tensor)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	const szType* shape = tensor.shape();
	szType iLength = std::accumulate(shape, shape + tensor.num_dimensions(), szType(1), std::multiplies<szType>());
	for (idxType i = 0; i < iLength; ++i)
		tensor.data()[i] = RandGaussian();
}



template <typename T, size_t N>
void TensorCopy(const boost::multi_array<T, N>& from, boost::multi_array<T, N>& to)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	const szType* shapeFrom = from.shape();
	const szType* shapeTo = to.shape();

	bool bEqual = std::equal(shapeFrom, shapeFrom + from.num_dimensions(), shapeTo);

	if (bEqual)
	{
		szType iLength = std::accumulate(shapeFrom, shapeFrom + from.num_dimensions(), szType(1), std::multiplies<szType>());
		for (idxType i = 0; i < iLength; ++i)
			to.data()[i] = from.data()[i];
	}
	else
	{
		std::vector<szType> vShape;
		vShape.assign(shapeFrom, shapeFrom + from.num_dimensions());
		to.resize(vShape);
		szType iLength = std::accumulate(shapeFrom, shapeFrom + from.num_dimensions(), szType(1), std::multiplies<szType>());
		for (idxType i = 0; i < iLength; ++i)
			to.data()[i] = from.data()[i];
	}
}



template <typename T, size_t N, typename Derived2>
void TensorCopy(const boost::multi_array<T, N>& from, Eigen::DenseBase<Derived2> const & to)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;
	const szType* shapeFrom = from.shape();
	szType iLength = std::accumulate(shapeFrom, shapeFrom + from.num_dimensions(), szType(1), std::multiplies<szType>());

	Eigen::DenseBase<Derived2>& to1 = const_cast<Eigen::DenseBase<Derived2>&>(to);

	if (iLength != to.size())
		to1.derived().resize(iLength);

	for (idxType i = 0; i < iLength; ++i)
		to1(i) = (from.data())[i];
}


template <typename T, size_t N>
void TensorCopy(const boost::multi_array<T, N>& from, std::vector<T>& to)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;
	const szType* shapeFrom = from.shape();
	szType iLength = std::accumulate(shapeFrom, shapeFrom + from.num_dimensions(), szType(1), std::multiplies<szType>());

	if (iLength != to.size())
		to.resize(iLength);

	for (idxType i = 0; i < iLength; ++i)
		to[i] = (from.data())[i];
}


template <typename T, size_t N, typename Derived2>
void TensorCopy(const Eigen::DenseBase<Derived2>& from, boost::multi_array<T, N>& to)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	szType iLength = to.num_elements();
	if (iLength != from.size())
		return;


	for (idxType i = 0; i < iLength; ++i)
		(to.data())[i] = from(i);
}


template <typename T, size_t N>
void TensorCopy(const std::vector<T>& from, boost::multi_array<T, N>& to)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	szType iLength = to.num_elements();
	if (iLength != from.size())
		return;


	for (idxType i = 0; i < iLength; ++i)
		(to.data())[i] = from[i];
}


template <typename T, size_t N, typename Derived2>
void TensorCopy(const Eigen::DenseBase<Derived2>& from, boost::multi_array<T, N>& to, const std::vector<typename boost::multi_array<T, N>::size_type>& shape)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	szType iLength = std::accumulate(shape.begin(), shape.end(), szType(1), std::multiplies<szType>());
	if (iLength != from.size())
		return;
	to.resize(shape);

	for (idxType i = 0; i < iLength; ++i)
		(to.data())[i] = from(i);
}


template <typename T, size_t N>
void TensorCopy(const std::vector<T>& from, boost::multi_array<T, N>& to, const std::vector<typename boost::multi_array<T, N>::size_type>& shape)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	szType iLength = std::accumulate(shape.begin(), shape.end(), szType(1), std::multiplies<szType>());
	if (iLength != from.size())
		return;
	to.resize(shape);

	for (idxType i = 0; i < iLength; ++i)
		(to.data())[i] = from[i];
}


template <typename T, typename Derived1>
void TensorProdMat(const boost::multi_array<T, 3>& tensor,
	const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 3>::index mode,
	boost::multi_array<T, 3>& ret)
{
	typedef typename boost::multi_array<T, 3>::size_type szType;
	typedef typename boost::multi_array<T, 3>::index idxType;
	const szType* inputShape = tensor.shape();

	if (inputShape[mode] != mat.cols())
		return;
	if (mode >= 3)
		return;

	szType J = mat.rows();
	boost::array<idxType, 3> outputShape = { { inputShape[0], inputShape[1], inputShape[2] } };
	outputShape[mode] = J;

	if (!std::equal(outputShape.begin(), outputShape.end(), ret.shape()))
		ret.resize(outputShape);

	for (idxType i1 = 0; i1 < outputShape[0]; ++i1)
	{
		for (idxType i2 = 0; i2 < outputShape[1]; ++i2)
		{
			for (idxType i3 = 0; i3 < outputShape[2]; ++i3)
			{
				boost::array<idxType, 3> idx = { { i1,i2,i3 } };
				boost::array<idxType, 3> idx1 = { { i1,i2,i3 } };
				T acumm = 0;
				for (idxType i_n = 0; i_n < inputShape[mode]; ++i_n)
				{
					idx1[mode] = i_n;
					acumm += tensor(idx1) * mat(idx[mode], i_n);
				}
				ret(idx) = acumm;
			}
		}
	}
}


template <typename T, typename Derived1>
void TensorProdMat(const boost::multi_array<T, 4>& tensor,
	const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 4>::index mode,
	boost::multi_array<T, 4>& ret)
{
	typedef typename boost::multi_array<T, 4>::size_type szType;
	typedef typename boost::multi_array<T, 4>::index idxType;
	const szType* inputShape = tensor.shape();

	if (inputShape[mode] != mat.cols())
		return;
	if (mode >= 4)
		return;

	szType J = mat.rows();
	boost::array<idxType, 4> outputShape = { { inputShape[0], inputShape[1], inputShape[2], inputShape[3] } };
	outputShape[mode] = J;

	if (!std::equal(outputShape.begin(), outputShape.end(), ret.shape()))
		ret.resize(outputShape);

	for (idxType i1 = 0; i1 < outputShape[0]; ++i1)
	{
		for (idxType i2 = 0; i2 < outputShape[1]; ++i2)
		{
			for (idxType i3 = 0; i3 < outputShape[2]; ++i3)
			{
				for (idxType i4 = 0; i4 < outputShape[3]; ++i4)
				{
					boost::array<idxType, 4> idx = { { i1,i2,i3,i4 } };
					boost::array<idxType, 4> idx1 = { { i1,i2,i3,i4 } };
					T acumm = 0;
					for (idxType i_n = 0; i_n < inputShape[mode]; ++i_n)
					{
						idx1[mode] = i_n;
						acumm += tensor(idx1) * mat(idx[mode], i_n);
					}
					ret(idx) = acumm;
				}
			}
		}
	}
}

template <typename T, typename Derived1>
void TensorProdMat(const boost::multi_array<T, 5>& tensor,
	const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 5>::index mode,
	boost::multi_array<T, 5>& ret)
{
	typedef typename boost::multi_array<T, 5>::size_type szType;
	typedef typename boost::multi_array<T, 5>::index idxType;
	const szType* inputShape = tensor.shape();

	if (inputShape[mode] != mat.cols())
		return;
	if (mode >= 5)
		return;

	szType J = mat.rows();
	boost::array<idxType, 5> outputShape = { { inputShape[0], inputShape[1], inputShape[2], inputShape[3], inputShape[4] } };
	outputShape[mode] = J;

	if (!std::equal(outputShape.begin(), outputShape.end(), ret.shape()))
		ret.resize(outputShape);

	for (idxType i1 = 0; i1 < outputShape[0]; ++i1)
	{
		for (idxType i2 = 0; i2 < outputShape[1]; ++i2)
		{
			for (idxType i3 = 0; i3 < outputShape[2]; ++i3)
			{
				for (idxType i4 = 0; i4 < outputShape[3]; ++i4)
				{
					for (idxType i5 = 0; i5 < outputShape[4]; ++i5)
					{
						boost::array<idxType, 5> idx = { { i1,i2,i3,i4,i5 } };
						boost::array<idxType, 5> idx1 = { { i1,i2,i3,i4,i5 } };
						T acumm = 0;
						for (idxType i_n = 0; i_n < inputShape[mode]; ++i_n)
						{
							idx1[mode] = i_n;
							acumm += tensor(idx1) * mat(idx[mode], i_n);
						}
						ret(idx) = acumm;
					}
				}
			}
		}
	}
}



template <typename T, typename Derived1>
void TensorProdMat(const boost::multi_array<T, 6>& tensor,
	const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 6>::index mode,
	boost::multi_array<T, 6>& ret)
{
	typedef typename boost::multi_array<T, 6>::size_type szType;
	typedef typename boost::multi_array<T, 6>::index idxType;
	const szType* inputShape = tensor.shape();

	if (inputShape[mode] != mat.cols())
		return;
	if (mode > 5)
		return;

	szType J = mat.rows();
	boost::array<idxType, 6> outputShape = { { inputShape[0], inputShape[1], inputShape[2], inputShape[3], inputShape[4], inputShape[5] } };
	outputShape[mode] = J;

	if (!std::equal(outputShape.begin(), outputShape.end(), ret.shape()))
		ret.resize(outputShape);

	for (idxType i1 = 0; i1 < outputShape[0]; ++i1)
	{
		for (idxType i2 = 0; i2 < outputShape[1]; ++i2)
		{
			for (idxType i3 = 0; i3 < outputShape[2]; ++i3)
			{
				for (idxType i4 = 0; i4 < outputShape[3]; ++i4)
				{
					for (idxType i5 = 0; i5 < outputShape[4]; ++i5)
					{
						for (idxType i6 = 0; i6 < outputShape[5]; ++i6)
						{
							boost::array<idxType, 6> idx = { { i1,i2,i3,i4,i5,i6 } };
							boost::array<idxType, 6> idx1 = { { i1,i2,i3,i4,i5,i6 } };
							T acumm = 0;
							for (idxType i_n = 0; i_n < inputShape[mode]; ++i_n)
							{
								idx1[mode] = i_n;
								acumm += tensor(idx1) * mat(idx[mode], i_n);
							}
							ret(idx) = acumm;
						}
					}
				}
			}
		}
	}
}

template<typename T>
inline void TensorProdMatMulti(const boost::multi_array<T, 3>& tensor,
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vMats,
	boost::multi_array<T, 3>& ret)
{
	typedef typename boost::multi_array<T, 3>::size_type szType;
	typedef typename boost::multi_array<T, 3>::index idxType;
	const szType* inputShape = tensor.shape();

	if (vMats.size() != 3)
		return;
	for (size_t i = 0; i < vMats.size(); ++i)
		if (inputShape[i] != vMats[i].cols())
			return;

	boost::array<idxType, 3> outputShape = { { vMats[0].rows(), vMats[1].rows(), vMats[2].rows() } };

	if (!std::equal(outputShape.begin(), outputShape.end(), ret.shape()))
		ret.resize(outputShape);

	for (idxType i1 = 0; i1 < outputShape[0]; ++i1)
	{
		for (idxType i2 = 0; i2 < outputShape[1]; ++i2)
		{
			for (idxType i3 = 0; i3 < outputShape[2]; ++i3)
			{
				boost::array<idxType, 3> idx = { { i1,i2,i3 } };
				T acumm = 0;
				for (size_t mode = 0; mode < 3; ++mode)
				{
					boost::array<idxType, 3> idx1 = { { i1,i2,i3 } };
					for (idxType i_n = 0; i_n < inputShape[mode]; ++i_n)
					{
						idx1[mode] = i_n;
						acumm += tensor(idx1) * vMats[mode](idx[mode], i_n);
					}
				}
				ret(idx) = acumm;
			}
		}
	}

}


template <typename T, typename Derived1>
void TensorUnfold(const boost::multi_array<T, 3>& tensor,
	typename boost::multi_array<T, 3>::index mode,
	Eigen::MatrixBase<Derived1>& ret)
{
	typedef typename boost::multi_array<T, 3>::size_type szType;
	typedef typename boost::multi_array<T, 3>::index idxType;
	const szType* inputShape = tensor.shape();

	idxType J = 1;
	for (idxType i = 0; i < 3; ++i)
	{
		if (i == mode)
			continue;
		J *= inputShape[i];
	}

	if (inputShape[mode] != ret.rows() || J != ret.cols())
		ret.derived().resize(inputShape[mode], J);

	for (idxType i3 = 0; i3 < inputShape[2]; ++i3)
	{
		for (idxType i2 = 0; i2 < inputShape[1]; ++i2)
		{
			for (idxType i1 = 0; i1 < inputShape[0]; ++i1)
			{
				idxType i[3] = { i1, i2, i3 };
				idxType sum = 0;
				for (idxType k = 0; k < 3; ++k)
				{
					if (k == mode)
						continue;
					idxType prod = 1;
					for (idxType m = 0; m < k; ++m)
					{
						if (m == mode)
							continue;
						prod *= inputShape[m];
					}
					sum += i[k] * prod;
				}
				ret(i[mode], sum) = tensor[i1][i2][i3];
			}
		}
	}
}


template <typename T, typename Derived1>
void TensorUnfold(const boost::multi_array<T, 4>& tensor,
	typename boost::multi_array<T, 4>::index mode,
	Eigen::MatrixBase<Derived1>& ret)
{
	typedef typename boost::multi_array<T, 4>::size_type szType;
	typedef typename boost::multi_array<T, 4>::index idxType;
	const szType* inputShape = tensor.shape();

	idxType J = 1;
	for (idxType i = 0; i < 4; ++i)
	{
		if (i == mode)
			continue;
		J *= inputShape[i];
	}

	if (inputShape[mode] != ret.rows() || J != ret.cols())
		ret.derived().resize(inputShape[mode], J);

	for (idxType i1 = 0; i1 < inputShape[0]; ++i1)
	{
		for (idxType i2 = 0; i2 < inputShape[1]; ++i2)
		{
			for (idxType i3 = 0; i3 < inputShape[2]; ++i3)
			{
				for (idxType i4 = 0; i4 < inputShape[3]; ++i4)
				{
					idxType i[4] = { i1, i2, i3, i4 };
					idxType sum = 0;
					for (idxType k = 0; k < 4; ++k)
					{
						if (k == mode)
							continue;
						idxType prod = 1;
						for (idxType m = 0; m < k; ++m)
						{
							if (m == mode)
								continue;
							prod *= inputShape[m];
						}
						sum += i[k] * prod;
					}
					ret(i[mode], sum) = tensor[i1][i2][i3][i4];
				}
			}
		}
	}
}

template <typename T, typename Derived1>
void TensorUnfold(const boost::multi_array<T, 5>& tensor,
	typename boost::multi_array<T, 5>::index mode,
	Eigen::MatrixBase<Derived1>& ret)
{
	typedef typename boost::multi_array<T, 5>::size_type szType;
	typedef typename boost::multi_array<T, 5>::index idxType;
	const szType* inputShape = tensor.shape();

	idxType J = 1;
	for (idxType i = 0; i < 5; ++i)
	{
		if (i == mode)
			continue;
		J *= inputShape[i];
	}

	if (inputShape[mode] != ret.rows() || J != ret.cols())
		ret.derived().resize(inputShape[mode], J);

	for (idxType i1 = 0; i1 < inputShape[0]; ++i1)
	{
		for (idxType i2 = 0; i2 < inputShape[1]; ++i2)
		{
			for (idxType i3 = 0; i3 < inputShape[2]; ++i3)
			{
				for (idxType i4 = 0; i4 < inputShape[3]; ++i4)
				{
					for (idxType i5 = 0; i5 < inputShape[4]; ++i5)
					{
						idxType i[5] = { i1, i2, i3, i4, i5 };
						idxType sum = 0;
						for (idxType k = 0; k < 5; ++k)
						{
							if (k == mode)
								continue;
							idxType prod = 1;
							for (idxType m = 0; m < k; ++m)
							{
								if (m == mode)
									continue;
								prod *= inputShape[m];
							}
							sum += i[k] * prod;
						}
						ret(i[mode], sum) = tensor[i1][i2][i3][i4][i5];
					}
				}
			}
		}
	}
}



template <typename T, typename Derived1>
void TensorUnfold(const boost::multi_array<T, 6>& tensor,
	typename boost::multi_array<T, 6>::index mode,
	Eigen::MatrixBase<Derived1>& ret)
{
	typedef typename boost::multi_array<T, 6>::size_type szType;
	typedef typename boost::multi_array<T, 6>::index idxType;
	const szType* inputShape = tensor.shape();

	idxType J = 1;
	for (idxType i = 0; i < 6; ++i)
	{
		if (i == mode)
			continue;
		J *= inputShape[i];
	}

	if (inputShape[mode] != ret.rows() || J != ret.cols())
		ret.derived().resize(inputShape[mode], J);

	for (idxType i1 = 0; i1 < inputShape[0]; ++i1)
	{
		for (idxType i2 = 0; i2 < inputShape[1]; ++i2)
		{
			for (idxType i3 = 0; i3 < inputShape[2]; ++i3)
			{
				for (idxType i4 = 0; i4 < inputShape[3]; ++i4)
				{
					for (idxType i5 = 0; i5 < inputShape[4]; ++i5)
					{
						for (idxType i6 = 0; i6 < inputShape[5]; ++i6)
						{
							idxType i[6] = { i1, i2, i3, i4, i5, i6 };
							idxType sum = 0;
							for (idxType k = 0; k < 6; ++k)
							{
								if (k == mode)
									continue;
								idxType prod = 1;
								for (idxType m = 0; m < k; ++m)
								{
									if (m == mode)
										continue;
									prod *= inputShape[m];
								}
								sum += i[k] * prod;
							}
							ret(i[mode], sum) = tensor[i1][i2][i3][i4][i5][i6];
						}
					}
				}
			}
		}
	}
}

template<typename T>
void TensorUnfold(const boost::multi_array<T, 3>& tensor,
	typename boost::multi_array<T, 3>::index mode,
	Eigen::SparseMatrix<T>& ret,
	size_t iNNZ)
{
	typedef typename boost::multi_array<T, 3>::size_type szType;
	typedef typename boost::multi_array<T, 3>::index idxType;
	const szType* inputShape = tensor.shape();

	idxType J = 1;
	for (idxType i = 0; i < 3; ++i)
	{
		if (i == mode)
			continue;
		J *= inputShape[i];
	}

	std::vector<Eigen::Triplet<T> > vNNZCoeffs;
	ret.resize(inputShape[mode], J);
	if (iNNZ > 0)
		vNNZCoeffs.reserve(iNNZ);


	for (idxType i3 = 0; i3 < inputShape[2]; ++i3)
	{
		for (idxType i2 = 0; i2 < inputShape[1]; ++i2)
		{
			for (idxType i1 = 0; i1 < inputShape[0]; ++i1)
			{
				if (fabs(tensor[i1][i2][i3]) > std::numeric_limits<T>::min())
				{
					idxType i[3] = { i1, i2, i3 };
					idxType sum = 0;
					for (idxType k = 0; k < 3; ++k)
					{
						if (k == mode)
							continue;
						idxType prod = 1;
						for (idxType m = 0; m < k; ++m)
						{
							if (m == mode)
								continue;
							prod *= inputShape[m];
						}
						sum += i[k] * prod;
					}
					//ret.insert(i[mode], sum) = tensor[i1][i2][i3];
					vNNZCoeffs.push_back(Eigen::Triplet<T>(i[mode], sum, tensor[i1][i2][i3]));
				}
			}
		}
	}

	ret.setFromTriplets(vNNZCoeffs.begin(), vNNZCoeffs.end());
}


template<typename T>
void TensorUnfold(const boost::multi_array<T, 4>& tensor,
	typename boost::multi_array<T, 4>::index mode,
	Eigen::SparseMatrix<T>& ret,
	size_t iNNZ)
{
	typedef typename boost::multi_array<T, 4>::size_type szType;
	typedef typename boost::multi_array<T, 4>::index idxType;
	const szType* inputShape = tensor.shape();

	idxType J = 1;
	for (idxType i = 0; i < 4; ++i)
	{
		if (i == mode)
			continue;
		J *= inputShape[i];
	}

	std::vector<Eigen::Triplet<T> > vNNZCoeffs;
	ret.resize(inputShape[mode], J);
	if (iNNZ > 0)
		vNNZCoeffs.reserve(iNNZ);

	for (idxType i1 = 0; i1 < inputShape[0]; ++i1)
	{
		for (idxType i2 = 0; i2 < inputShape[1]; ++i2)
		{
			for (idxType i3 = 0; i3 < inputShape[2]; ++i3)
			{
				for (idxType i4 = 0; i4 < inputShape[3]; ++i4)
				{
					if (fabs(tensor[i1][i2][i3][i4]) > std::numeric_limits<T>::min())
					{
						idxType i[4] = { i1, i2, i3, i4 };
						idxType sum = 0;
						for (idxType k = 0; k < 4; ++k)
						{
							if (k == mode)
								continue;
							idxType prod = 1;
							for (idxType m = 0; m < k; ++m)
							{
								if (m == mode)
									continue;
								prod *= inputShape[m];
							}
							sum += i[k] * prod;
						}
						//ret.insert(i[mode], sum) = tensor[i1][i2][i3][i4];
						vNNZCoeffs.push_back(Eigen::Triplet<T>(i[mode], sum, tensor[i1][i2][i3][i4]));
					}
				}
			}
		}
	}

	ret.setFromTriplets(vNNZCoeffs.begin(), vNNZCoeffs.end());
}


template<typename T>
void TensorUnfold(const boost::multi_array<T, 5>& tensor,
	typename boost::multi_array<T, 5>::index mode,
	Eigen::SparseMatrix<T>& ret,
	size_t iNNZ)
{
	typedef typename boost::multi_array<T, 5>::size_type szType;
	typedef typename boost::multi_array<T, 5>::index idxType;
	const szType* inputShape = tensor.shape();

	idxType J = 1;
	for (idxType i = 0; i < 5; ++i)
	{
		if (i == mode)
			continue;
		J *= inputShape[i];
	}

	std::vector<Eigen::Triplet<T> > vNNZCoeffs;
	ret.resize(inputShape[mode], J);
	if (iNNZ > 0)
		vNNZCoeffs.reserve(iNNZ);

	for (idxType i1 = 0; i1 < inputShape[0]; ++i1)
	{
		for (idxType i2 = 0; i2 < inputShape[1]; ++i2)
		{
			for (idxType i3 = 0; i3 < inputShape[2]; ++i3)
			{
				for (idxType i4 = 0; i4 < inputShape[3]; ++i4)
				{
					for (idxType i5 = 0; i5 < inputShape[4]; ++i5)
					{
						if (fabs(tensor[i1][i2][i3][i4][i5]) > std::numeric_limits<T>::min())
						{
							idxType i[5] = { i1, i2, i3, i4, i5 };
							idxType sum = 0;
							for (idxType k = 0; k < 5; ++k)
							{
								if (k == mode)
									continue;
								idxType prod = 1;
								for (idxType m = 0; m < k; ++m)
								{
									if (m == mode)
										continue;
									prod *= inputShape[m];
								}
								sum += i[k] * prod;
							}
							//ret.insert(i[mode], sum) = tensor[i1][i2][i3][i4][i5];
							vNNZCoeffs.push_back(Eigen::Triplet<T>(i[mode], sum, tensor[i1][i2][i3][i4][i5]));
						}
					}
				}
			}
		}
	}

	ret.setFromTriplets(vNNZCoeffs.begin(), vNNZCoeffs.end());
}


template<typename T>
void TensorUnfold(const boost::multi_array<T, 6>& tensor,
	typename boost::multi_array<T, 6>::index mode,
	Eigen::SparseMatrix<T>& ret,
	size_t iNNZ)
{
	typedef typename boost::multi_array<T, 6>::size_type szType;
	typedef typename boost::multi_array<T, 6>::index idxType;
	const szType* inputShape = tensor.shape();

	idxType J = 1;
	for (idxType i = 0; i < 6; ++i)
	{
		if (i == mode)
			continue;
		J *= inputShape[i];
	}

	std::vector<Eigen::Triplet<T> > vNNZCoeffs;
	ret.resize(inputShape[mode], J);
	if (iNNZ > 0)
		vNNZCoeffs.reserve(iNNZ);

	for (idxType i1 = 0; i1 < inputShape[0]; ++i1)
	{
		for (idxType i2 = 0; i2 < inputShape[1]; ++i2)
		{
			for (idxType i3 = 0; i3 < inputShape[2]; ++i3)
			{
				for (idxType i4 = 0; i4 < inputShape[3]; ++i4)
				{
					for (idxType i5 = 0; i5 < inputShape[4]; ++i5)
					{
						for (idxType i6 = 0; i6 < inputShape[5]; ++i6)
						{
							if (fabs(tensor[i1][i2][i3][i4][i5][i6]) > std::numeric_limits<T>::min())
							{
								idxType i[6] = { i1, i2, i3, i4, i5, i6 };
								idxType sum = 0;
								for (idxType k = 0; k < 6; ++k)
								{
									if (k == mode)
										continue;
									idxType prod = 1;
									for (idxType m = 0; m < k; ++m)
									{
										if (m == mode)
											continue;
										prod *= inputShape[m];
									}
									sum += i[k] * prod;
								}
								vNNZCoeffs.push_back(Eigen::Triplet<T>(i[mode], sum, tensor[i1][i2][i3][i4][i5][i6]));
								//ret.insert(i[mode], sum) = tensor[i1][i2][i3][i4][i5][i6];
							}
						}
					}
				}
			}
		}
	}

	ret.setFromTriplets(vNNZCoeffs.begin(), vNNZCoeffs.end());
}


template <typename T, typename Idx>
void SparseTensorUnfold(Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc, 
	Eigen::Matrix<T, Eigen::Dynamic, 1>& val, 
	const std::vector<size_t>& vShape, 
	size_t iMode,
	Eigen::SparseMatrix<T>& ret)
{
	size_t J = 1;
	for (size_t i = 0; i < loc.rows(); ++i)
	{
		if (i == iMode)
			continue;
		J *= vShape[i];
	}
	ret.resize(vShape[iMode], J);

	std::vector<Eigen::Triplet<T> > vNNZCoeffs;
	for (size_t i = 0; i < loc.cols(); ++i)
	{
		size_t sum = 0;
		for (size_t k = 0; k < loc.rows(); ++k)
		{
			if (k == iMode)
				continue;
			size_t prod = 1;
			for (size_t m = 0; m < k; ++m)
			{
				if (m == iMode)
					continue;
				prod *= vShape[m];
			}
			sum += loc(k, i) * prod;
		}
		vNNZCoeffs.push_back(Eigen::Triplet<T>(loc(iMode, i), sum, val(i)));
	}
	ret.setFromTriplets(vNNZCoeffs.begin(), vNNZCoeffs.end());
}


template <typename Idx>
void SparseTensorUnfold(Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	const std::vector<size_t>& vShape,
	size_t iMode,
	Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& retLoc)
{

	retLoc.resize(2, loc.cols());

	for (size_t i = 0; i < loc.cols(); ++i)
	{
		size_t sum = 0;
		for (size_t k = 0; k < loc.rows(); ++k)
		{
			if (k == iMode)
				continue;
			size_t prod = 1;
			for (size_t m = 0; m < k; ++m)
			{
				if (m == iMode)
					continue;
				prod *= vShape[m];
			}
			sum += loc(k, i) * prod;
		}
		retLoc(0, i) = loc(iMode, i);
		retLoc(1, i) = sum;
	}
}


template <typename T, typename Derived1>
void TensorFold(const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 3>::index mode,
	boost::multi_array<T, 3>& tensor)
{
	typedef typename boost::multi_array<T, 3>::size_type szType;
	typedef typename boost::multi_array<T, 3>::index idxType;
	const szType* inputShape = tensor.shape();

	idxType J = 1;
	for (idxType i = 0; i < 3; ++i)
	{
		if (i == mode)
			continue;
		J *= inputShape[i];
	}

	for (idxType i3 = 0; i3 < inputShape[2]; ++i3)
	{
		for (idxType i2 = 0; i2 < inputShape[1]; ++i2)
		{
			for (idxType i1 = 0; i1 < inputShape[0]; ++i1)
			{
				idxType i[3] = { i1, i2, i3 };
				idxType sum = 0;
				for (idxType k = 0; k < 3; ++k)
				{
					if (k == mode)
						continue;
					idxType prod = 1;
					for (idxType m = 0; m < k; ++m)
					{
						if (m == mode)
							continue;
						prod *= inputShape[m];
					}
					sum += i[k] * prod;
				}
				tensor[i1][i2][i3] = mat(i[mode], sum);
			}
		}
	}
}

template <typename T, typename Derived1>
void TensorFold(const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 4>::index mode,
	boost::multi_array<T, 4>& tensor)
{
	typedef typename boost::multi_array<T, 4>::size_type szType;
	typedef typename boost::multi_array<T, 4>::index idxType;
	const szType* inputShape = tensor.shape();

	idxType J = 1;
	for (idxType i = 0; i < 4; ++i)
	{
		if (i == mode)
			continue;
		J *= inputShape[i];
	}

	for (idxType i1 = 0; i1 < inputShape[0]; ++i1)
	{
		for (idxType i2 = 0; i2 < inputShape[1]; ++i2)
		{
			for (idxType i3 = 0; i3 < inputShape[2]; ++i3)
			{
				for (idxType i4 = 0; i4 < inputShape[3]; ++i4)
				{
					idxType i[4] = { i1, i2, i3, i4 };
					idxType sum = 0;
					for (idxType k = 0; k < 4; ++k)
					{
						if (k == mode)
							continue;
						idxType prod = 1;
						for (idxType m = 0; m < k; ++m)
						{
							if (m == mode)
								continue;
							prod *= inputShape[m];
						}
						sum += i[k] * prod;
					}
					tensor[i1][i2][i3][i4] = mat(i[mode], sum);
				}
			}
		}
	}
}

template <typename T, typename Derived1>
void TensorFold(const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 5>::index mode,
	boost::multi_array<T, 5>& tensor)
{
	typedef typename boost::multi_array<T, 5>::size_type szType;
	typedef typename boost::multi_array<T, 5>::index idxType;
	const szType* inputShape = tensor.shape();

	idxType J = 1;
	for (idxType i = 0; i < 5; ++i)
	{
		if (i == mode)
			continue;
		J *= inputShape[i];
	}

	for (idxType i1 = 0; i1 < inputShape[0]; ++i1)
	{
		for (idxType i2 = 0; i2 < inputShape[1]; ++i2)
		{
			for (idxType i3 = 0; i3 < inputShape[2]; ++i3)
			{
				for (idxType i4 = 0; i4 < inputShape[3]; ++i4)
				{
					for (idxType i5 = 0; i5 < inputShape[4]; ++i5)
					{
						idxType i[5] = { i1, i2, i3, i4, i5 };
						idxType sum = 0;
						for (idxType k = 0; k < 5; ++k)
						{
							if (k == mode)
								continue;
							idxType prod = 1;
							for (idxType m = 0; m < k; ++m)
							{
								if (m == mode)
									continue;
								prod *= inputShape[m];
							}
							sum += i[k] * prod;
						}
						tensor[i1][i2][i3][i4][i5] = mat(i[mode], sum);
					}
				}
			}
		}
	}
}

template <typename T, typename Derived1>
void TensorFold(const Eigen::MatrixBase<Derived1>& mat,
	typename boost::multi_array<T, 6>::index mode,
	boost::multi_array<T, 6>& tensor)
{
	typedef typename boost::multi_array<T, 6>::size_type szType;
	typedef typename boost::multi_array<T, 6>::index idxType;
	const szType* inputShape = tensor.shape();

	idxType J = 1;
	for (idxType i = 0; i < 6; ++i)
	{
		if (i == mode)
			continue;
		J *= inputShape[i];
	}

	for (idxType i1 = 0; i1 < inputShape[0]; ++i1)
	{
		for (idxType i2 = 0; i2 < inputShape[1]; ++i2)
		{
			for (idxType i3 = 0; i3 < inputShape[2]; ++i3)
			{
				for (idxType i4 = 0; i4 < inputShape[3]; ++i4)
				{
					for (idxType i5 = 0; i5 < inputShape[4]; ++i5)
					{
						for (idxType i6 = 0; i6 < inputShape[5]; ++i6)
						{
							idxType i[6] = { i1, i2, i3, i4, i5, i6 };
							idxType sum = 0;
							for (idxType k = 0; k < 6; ++k)
							{
								if (k == mode)
									continue;
								idxType prod = 1;
								for (idxType m = 0; m < k; ++m)
								{
									if (m == mode)
										continue;
									prod *= inputShape[m];
								}
								sum += i[k] * prod;
							}
							tensor[i1][i2][i3][i4][i5][i6] = mat(i[mode], sum);
						}
					}
				}
			}
		}
	}
}

template <typename T, size_t N>
void TensorAdd(const boost::multi_array<T, N>& input1, boost::multi_array<T, N>& ret)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	if (!std::equal(input1.shape(), input1.shape() + input1.num_dimensions(), ret.shape()))
		return;

	szType iLength = input1.num_elements();
	for (idxType i = 0; i < iLength; ++i)
		ret.data()[i] += input1.data()[i];
}

template <typename T, size_t N>
void TensorAdd(const boost::multi_array<T, N>& input1, const boost::multi_array<T, N>& input2, boost::multi_array<T, N>& ret)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	if (!std::equal(input1.shape(), input1.shape() + input1.num_dimensions(), ret.shape()))
	{
		std::vector<szType> shape;
		shape.assign(input1.shape(), input1.shape() + input1.num_dimensions());
		ret.resize(shape);
	}

	szType iLength = std::accumulate(input1.shape(), input1.shape() + input1.num_dimensions(), szType(1), std::multiplies<szType>());
	for (idxType i = 0; i < iLength; ++i)
		ret.data()[i] = input1.data()[i] + input2.data()[i];
}

template <typename T, size_t N>
void TensorAdd(T input1, boost::multi_array<T, N>& ret)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;
	szType iLength = ret.num_elements();
	for (idxType i = 0; i < iLength; ++i)
		ret.data()[i] += input1;
}


template <typename T, size_t N>
void TensorSubtract(const boost::multi_array<T, N>& input1, boost::multi_array<T, N>& ret)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;
	szType iLength = input1.num_elements();
	if (iLength != ret.num_elements())
	{
		std::cout << "error in TensorSubtract()" << std::endl;
		return;
	}
	for (idxType i = 0; i < iLength; ++i)
		ret.data()[i] -= input1.data()[i];
}

template <typename T, size_t N>
void TensorSubtract(const boost::multi_array<T, N>& input1, const boost::multi_array<T, N>& input2, boost::multi_array<T, N>& ret)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	if (!std::equal(input1.shape(), input1.shape() + input1.num_dimensions(), ret.shape()))
	{
		std::vector<szType> shape;
		shape.assign(input1.shape(), input1.shape() + N);
		ret.resize(shape);
	}

	szType iLength = std::accumulate(input1.shape(), input1.shape() + input1.num_dimensions(), szType(1), std::multiplies<szType>());

	for (idxType i = 0; i < iLength; ++i)
		ret.data()[i] = input1.data()[i] - input2.data()[i];
}

template <typename T, size_t N>
void TensorSubtract(T input1, boost::multi_array<T, N>& ret)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;
	szType iLength = ret.num_elements();
	for (idxType i = 0; i < iLength; ++i)
		ret.data()[i] -= input1;
}


template <typename T, size_t N>
void TensorProdVal(const boost::multi_array<T, N>& input1, T input2, boost::multi_array<T, N>& ret)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;
	std::vector<szType> shape;
	shape.assign(input1.shape(), input1.shape() + N);
	ret.resize(shape);
	szType iLength = std::accumulate(shape.begin(), shape.end(), szType(1), std::multiplies<szType>());
	for (idxType i = 0; i < iLength; ++i)
		ret.data()[i] = input2 * input1.data()[i];
}


template <typename T, size_t N>
void TensorProdVal(T input1, boost::multi_array<T, N>& ret)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;
	std::vector<szType> shape;
	shape.assign(ret.shape(), ret.shape() + N);
	ret.resize(shape);
	szType iLength = std::accumulate(shape.begin(), shape.end(), szType(1), std::multiplies<szType>());
	for (idxType i = 0; i < iLength; ++i)
		ret.data()[i] *= input1;
}


template <typename T, size_t N>
void TensorPow(boost::multi_array<T, N>& tns, T pwr)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	for (idxType i = 0; i < tns.num_elements(); ++i)
		tns.data()[i] = pow(double(tns.data()[i]), double(pwr));
}


template <typename T, size_t N>
void TensorPow(boost::multi_array<T, N>& tns, T pwr, boost::multi_array<T, N>& ret)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	for (idxType i = 0; i < tns.num_elements(); ++i)
		ret.data()[i] = pow(double(tns.data()[i]), double(pwr));
}


template <size_t N>
CCS_INTERNAL_TYPE TensorNorm2(const boost::multi_array<CCS_INTERNAL_TYPE, N>& tensor)
{
	typedef typename Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1> Vec;
	return Vec::Map(tensor.data(), tensor.num_elements()).norm();
}


template <size_t N>
CCS_INTERNAL_TYPE TensorNorm2Squared(const boost::multi_array<CCS_INTERNAL_TYPE, N>& tensor)
{
	typedef typename Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1> Vec;
	return Vec::Map(tensor.data(), tensor.num_elements()).squaredNorm();
}


template <size_t N>
CCS_INTERNAL_TYPE TensorNorm1(const boost::multi_array<CCS_INTERNAL_TYPE, N>& tensor)
{
	typedef typename Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1> Vec;
	return Vec::Map(tensor.data(), tensor.num_elements()).template lpNorm<1>();
}

template<size_t N>
CCS_INTERNAL_TYPE TensorNormInf(const boost::multi_array<CCS_INTERNAL_TYPE, N>& tensor)
{
	typedef typename Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1> Vec;
	return Vec::Map(tensor.data(), tensor.num_elements()).cwiseAbs().maxCoeff();
}




template<typename T, typename Idx, typename Derived1>
size_t DenseToSparse(const Eigen::MatrixBase<Derived1>& mat,
	Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	size_t iNNZ)
{
	size_t iCounter = 0;
	size_t iDims = 0;
	if (mat.rows() > 1 && mat.cols() > 1)
		iDims = 2;
	else if ((mat.rows() > 1 && mat.cols() == 1) || (mat.cols() > 1 && mat.rows() == 1))
		iDims = 1;
	else
		return 0;

	if (loc.rows() != iDims || loc.cols() != iNNZ)
		loc.setZero(iDims, iNNZ);
	if (val.size() != iNNZ)
		val.setZero(iNNZ);

	for (size_t i2 = 0; i2 < mat.cols(); ++i2)
	{
		for (size_t i1 = 0; i1 < mat.rows(); ++i1)
		{
			if ((fabs(mat(i1, i2)) > std::numeric_limits<T>::min()) && (iCounter < iNNZ))
			{
				if (iDims == 1)
				{
					loc(0, iCounter) = i1 + i2 * mat.rows();
				}
				else
				{
					loc(0, iCounter) = i1;
					loc(1, iCounter) = i2;
				}
				val(iCounter) = mat(i1, i2);
				iCounter++;
			}
		}
	}

	if (iCounter < iNNZ)
	{
		loc = loc.block(0, 0, iDims, iCounter).eval();
		val = val.segment(0, iCounter).eval();
	}

	return iCounter;
}


template<typename T, typename Idx>
size_t DenseToSparse(const boost::multi_array<T, 3>& tensor,
	Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	size_t iNNZ)
{
	typedef typename boost::multi_array<T, 3>::size_type szType;
	typedef typename boost::multi_array<T, 3>::index idxType;
	const szType* pShape = tensor.shape();

	if (loc.cols() != iNNZ || loc.rows() != 3)
		loc.setZero(3, iNNZ);
	if (val.size() != iNNZ)
		val.setZero(iNNZ);

	size_t iCounter = 0;
	for (idxType i3 = 0; i3 < pShape[2]; ++i3)
	{
		for (idxType i2 = 0; i2 < pShape[1]; ++i2)
		{
			for (idxType i1 = 0; i1 < pShape[0]; ++i1)
			{
				T tmp = tensor[i1][i2][i3];
				if ((fabs(tmp) > std::numeric_limits<T>::min()) && (iCounter < iNNZ))
				{
					loc(0, iCounter) = i1;
					loc(1, iCounter) = i2;
					loc(2, iCounter) = i3;
					val(iCounter) = tmp;
					iCounter++;
				}
			}
		}
	}

	if (iCounter < iNNZ)
	{
		loc = loc.block(0, 0, 3, iCounter).eval();
		val = val.segment(0, iCounter).eval();
	}

	return iCounter;
}

template <typename T, typename Idx>
size_t DenseToSparse(const boost::multi_array<T, 4>& tensor,
	Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& val, size_t iNNZ)
{
	typedef typename boost::multi_array<T, 4>::size_type szType;
	typedef typename boost::multi_array<T, 4>::index idxType;
	const szType* pShape = tensor.shape();

	if (loc.cols() != iNNZ || loc.rows() != 4)
		loc.setZero(4, iNNZ);
	if (val.size() != iNNZ)
		val.setZero(iNNZ);

	size_t iCounter = 0;
	for (idxType i4 = 0; i4 < pShape[3]; ++i4)
	{
		for (idxType i3 = 0; i3 < pShape[2]; ++i3)
		{
			for (idxType i2 = 0; i2 < pShape[1]; ++i2)
			{
				for (idxType i1 = 0; i1 < pShape[0]; ++i1)
				{
					T tmp = tensor[i1][i2][i3][i4];
					if ((fabs(tmp) > std::numeric_limits<T>::min()) && (iCounter < iNNZ))
					{
						loc(0, iCounter) = i1;
						loc(1, iCounter) = i2;
						loc(2, iCounter) = i3;
						loc(3, iCounter) = i4;
						val(iCounter) = tmp;
						iCounter++;
					}
				}
			}
		}
	}

	if (iCounter < iNNZ)
	{
		loc = loc.block(0, 0, 4, iCounter).eval();
		val = val.segment(0, iCounter).eval();
	}

	return iCounter;
}

template <typename T, typename Idx>
size_t DenseToSparse(const boost::multi_array<T, 5>& tensor,
	Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& val, size_t iNNZ)
{
	typedef typename boost::multi_array<T, 5>::size_type szType;
	typedef typename boost::multi_array<T, 5>::index idxType;
	const szType* pShape = tensor.shape();

	if (loc.cols() != iNNZ || loc.rows() != 5)
		loc.setZero(5, iNNZ);
	if (val.size() != iNNZ)
		val.setZero(iNNZ);

	size_t iCounter = 0;
	for (idxType i5 = 0; i5 < pShape[4]; ++i5)
	{
		for (idxType i4 = 0; i4 < pShape[3]; ++i4)
		{
			for (idxType i3 = 0; i3 < pShape[2]; ++i3)
			{
				for (idxType i2 = 0; i2 < pShape[1]; ++i2)
				{
					for (idxType i1 = 0; i1 < pShape[0]; ++i1)
					{
						T tmp = tensor[i1][i2][i3][i4][i5];
						if ((fabs(tmp) > std::numeric_limits<T>::min()) && (iCounter < iNNZ))
						{
							loc(0, iCounter) = i1;
							loc(1, iCounter) = i2;
							loc(2, iCounter) = i3;
							loc(3, iCounter) = i4;
							loc(4, iCounter) = i5;
							val(iCounter) = tmp;
							iCounter++;
						}
					}
				}
			}
		}
	}

	if (iCounter < iNNZ)
	{
		loc = loc.block(0, 0, 5, iCounter).eval();
		val = val.segment(0, iCounter).eval();
	}

	return iCounter;
}

template <typename T, typename Idx>
size_t DenseToSparse(const boost::multi_array<T, 6>& tensor,
	Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& val, size_t iNNZ)
{
	typedef typename boost::multi_array<T, 6>::size_type szType;
	typedef typename boost::multi_array<T, 6>::index idxType;
	const szType* pShape = tensor.shape();

	if (loc.cols() != iNNZ || loc.rows() != 6)
		loc.setZero(6, iNNZ);
	if (val.size() != iNNZ)
		val.setZero(iNNZ);

	size_t iCounter = 0;
	for (idxType i6 = 0; i6 < pShape[5]; ++i6)
	{
		for (idxType i5 = 0; i5 < pShape[4]; ++i5)
		{
			for (idxType i4 = 0; i4 < pShape[3]; ++i4)
			{
				for (idxType i3 = 0; i3 < pShape[2]; ++i3)
				{
					for (idxType i2 = 0; i2 < pShape[1]; ++i2)
					{
						for (idxType i1 = 0; i1 < pShape[0]; ++i1)
						{
							T tmp = tensor[i1][i2][i3][i4][i5][i6];
							if ((fabs(tmp) > std::numeric_limits<T>::min()) && (iCounter < iNNZ))
							{
								loc(0, iCounter) = i1;
								loc(1, iCounter) = i2;
								loc(2, iCounter) = i3;
								loc(3, iCounter) = i4;
								loc(4, iCounter) = i5;
								loc(5, iCounter) = i6;
								val(iCounter) = tmp;
								iCounter++;
							}
						}
					}
				}
			}
		}
	}

	if (iCounter < iNNZ)
	{
		loc = loc.block(0, 0, 6, iCounter).eval();
		val = val.segment(0, iCounter).eval();
	}

	return iCounter;
}


template <typename T, typename Idx>
void SparseToDense(const Eigen::Matrix<Idx, Eigen::Dynamic, 1>& loc,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	Eigen::Matrix<T, Eigen::Dynamic, 1>& ret)
{
	if (loc.size() > ret.size() || loc.size() != val.size())
		return;

	ret.setZero();
	for (size_t i = 0; i < loc.size(); ++i)
		ret(loc(i)) = val(i);
}

template <typename T, typename Idx>
void SparseToDense(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& ret)
{
	if (loc.cols() > ret.size() || loc.cols() != val.rows() || loc.rows() != 2)
		return;

	ret.setZero();
	for (size_t i = 0; i < loc.cols(); ++i)
		ret(loc(0, i), loc(1, i)) = val(i);
}

template <typename T, typename Idx>
void SparseToDense(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	boost::multi_array<T, 3>& ret)
{
	if (loc.cols() > ret.num_elements() || loc.cols() != val.rows() || loc.rows() != 3)
		return;

	TensorSetZero(ret);
	for (size_t i = 0; i < loc.cols(); ++i)
		ret[loc(0, i)][loc(1, i)][loc(2, i)] = val(i);
}

template <typename T, typename Idx>
void SparseToDense(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	boost::multi_array<T, 4>& ret)
{
	if (loc.cols() > ret.num_elements() || loc.cols() != val.rows() || loc.rows() != 4)
		return;

	TensorSetZero(ret);
	for (size_t i = 0; i < loc.cols(); ++i)
		ret[loc(0, i)][loc(1, i)][loc(2, i)][loc(3, i)] = val(i);
}

template <typename T, typename Idx>
void SparseToDense(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	boost::multi_array<T, 5>& ret)
{
	if (loc.cols() > ret.num_elements() || loc.cols() != val.rows() || loc.rows() != 5)
		return;

	TensorSetZero(ret);
	for (size_t i = 0; i < loc.cols(); ++i)
		ret[loc(0, i)][loc(1, i)][loc(2, i)][loc(3, i)][loc(4, i)] = val(i);
}

template <typename T, typename Idx>
void SparseToDense(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	boost::multi_array<T, 6>& ret)
{
	if (loc.cols() > ret.num_elements() || loc.cols() != val.rows() || loc.rows() != 6)
		return;

	TensorSetZero(ret);
	for (size_t i = 0; i < loc.cols(); ++i)
		ret[loc(0, i)][loc(1, i)][loc(2, i)][loc(3, i)][loc(4, i)][loc(5, i)] = val(i);
}


template <typename T, typename Idx>
void SparseMatToEigenSparse(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& val, size_t iRows, size_t iCols, Eigen::SparseMatrix<T>& ret)
{
	std::vector<Eigen::Triplet<T> > vNNZCoeffs;
	ret.resize(iRows, iCols);
	vNNZCoeffs.reserve(val.size());

	for (size_t i = 0; i < val.size(); ++i)
		vNNZCoeffs.push_back(Eigen::Triplet<T>(loc(0, i), loc(1, i), val(i)));
	ret.setFromTriplets(vNNZCoeffs.begin(), vNNZCoeffs.end());
}


template <typename T, typename Idx>
void HOSVDReconSp(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc, const Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vU, Eigen::Matrix<T, Eigen::Dynamic, 1>& coreTensor)
{
	if (loc.rows() != 1 || vU.size() != 1 || loc.cols() != val.size())
		return;

	size_t iSparsity = loc.cols();

	if (coreTensor.size() != vU[0].rows())
		coreTensor.resize(vU[0].rows());

	for (size_t i1 = 0; i1 < coreTensor.size(); ++i1)
	{
		T accum = 0;
		for (size_t s = 0; s < iSparsity; ++s)
			accum += val(s) * vU[0](i1, loc(0, s));
		coreTensor(i1) = accum;
	}
}


template <typename T, typename Idx>
void HOSVDReconSp(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc, const Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vU, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& coreTensor)
{
	if (loc.rows() != 2 || vU.size() != 2 || loc.cols() != val.size())
		return;

	size_t iSparsity = loc.cols();

	if ((coreTensor.rows() != vU[0].rows()) || (coreTensor.cols() != vU[1].rows()))
		coreTensor.resize(vU[0].rows(), vU[1].rows());

	for (size_t i2 = 0; i2 < coreTensor.cols(); ++i2)
	{
		for (size_t i1 = 0; i1 < coreTensor.rows(); ++i1)
		{
			T accum = 0;
			for (size_t s = 0; s < iSparsity; ++s)
				accum += val(s) * vU[0](i1, loc(0, s)) * vU[1](i2, loc(1, s));
			coreTensor(i1, i2) = accum;
		}
	}
}

template <typename T, typename Idx>
void HOSVDReconSp(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc, const Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vU, boost::multi_array<T, 3>& coreTensor)
{
	typedef typename boost::multi_array<T, 3>::size_type szType;
	typedef typename boost::multi_array<T, 3>::index idxType;

	if (loc.rows() != 3 || vU.size() != 3 || loc.cols() != val.size())
		return;

	const szType* pShape = coreTensor.shape();
	size_t iSparsity = loc.cols();

	if ((pShape[0] != vU[0].rows()) || (pShape[1] != vU[1].rows()) || (pShape[2] != vU[2].rows()))
	{
		std::array<szType, 3> newShape = { vU[0].rows(), vU[1].rows(), vU[2].rows() };
		coreTensor.resize(newShape);
	}

	for (idxType i3 = 0; i3 < pShape[2]; ++i3)
	{
		for (idxType i2 = 0; i2 < pShape[1]; ++i2)
		{
			for (idxType i1 = 0; i1 < pShape[0]; ++i1)
			{
				T accum = 0;
				for (size_t s = 0; s < iSparsity; ++s)
					accum += val(s) * vU[0](i1, loc(0, s)) * vU[1](i2, loc(1, s)) * vU[2](i3, loc(2, s));
				coreTensor[i1][i2][i3] = accum;
			}
		}
	}
}

template <typename T, typename Idx>
void HOSVDReconSp(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc, const Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vU, boost::multi_array<T, 4>& coreTensor)
{
	typedef typename boost::multi_array<T, 4>::size_type szType;
	typedef typename boost::multi_array<T, 4>::index idxType;

	if (loc.rows() != 4 || vU.size() != 4 || loc.cols() != val.size())
		return;

	const szType* pShape = coreTensor.shape();
	size_t iSparsity = loc.cols();

	if ((pShape[0] != vU[0].rows()) || (pShape[1] != vU[1].rows()) || (pShape[2] != vU[2].rows()) || (pShape[3] != vU[3].rows()))
	{
		std::array<szType, 4> newShape = { vU[0].rows(), vU[1].rows(), vU[2].rows(), vU[3].rows() };
		coreTensor.resize(newShape);
	}

	for (idxType i4 = 0; i4 < pShape[3]; ++i4)
	{
		for (idxType i3 = 0; i3 < pShape[2]; ++i3)
		{
			for (idxType i2 = 0; i2 < pShape[1]; ++i2)
			{
				for (idxType i1 = 0; i1 < pShape[0]; ++i1)
				{
					T accum = 0;
					for (size_t s = 0; s < iSparsity; ++s)
						accum += val(s) * vU[0](i1, loc(0, s)) * vU[1](i2, loc(1, s)) * vU[2](i3, loc(2, s)) * vU[3](i4, loc(3, s));
					coreTensor[i1][i2][i3][i4] = accum;
				}
			}
		}
	}
}

template <typename T, typename Idx>
void HOSVDReconSp(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc, const Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vU, boost::multi_array<T, 5>& coreTensor)
{
	typedef typename boost::multi_array<T, 5>::size_type szType;
	typedef typename boost::multi_array<T, 5>::index idxType;

	if (loc.rows() != 5 || vU.size() != 5 || loc.cols() != val.size())
		return;

	const szType* pShape = coreTensor.shape();
	size_t iSparsity = loc.cols();

	if ((pShape[0] != vU[0].rows()) || (pShape[1] != vU[1].rows()) || (pShape[2] != vU[2].rows()) || (pShape[3] != vU[3].rows()) || (pShape[4] != vU[4].rows()))
	{
		std::array<szType, 5> newShape = { vU[0].rows(), vU[1].rows(), vU[2].rows(), vU[3].rows(), vU[4].rows() };
		coreTensor.resize(newShape);
	}

	for (idxType i5 = 0; i5 < pShape[4]; ++i5)
	{
		for (idxType i4 = 0; i4 < pShape[3]; ++i4)
		{
			for (idxType i3 = 0; i3 < pShape[2]; ++i3)
			{
				for (idxType i2 = 0; i2 < pShape[1]; ++i2)
				{
					for (idxType i1 = 0; i1 < pShape[0]; ++i1)
					{
						T accum = 0;
						for (size_t s = 0; s < iSparsity; ++s)
							accum += val(s) * vU[0](i1, loc(0, s)) * vU[1](i2, loc(1, s)) * vU[2](i3, loc(2, s)) * vU[3](i4, loc(3, s)) * vU[4](i5, loc(4, s));
						coreTensor[i1][i2][i3][i4][i5] = accum;
					}
				}
			}
		}
	}
}

template <typename T, typename Idx>
void HOSVDReconSp(const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>& loc, const Eigen::Matrix<T, Eigen::Dynamic, 1>& val,
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& vU, boost::multi_array<T, 6>& coreTensor)
{
	typedef typename boost::multi_array<T, 6>::size_type szType;
	typedef typename boost::multi_array<T, 6>::index idxType;

	if (loc.rows() != 6 || vU.size() != 6 || loc.cols() != val.size())
		return;

	const szType* pShape = coreTensor.shape();
	size_t iSparsity = loc.cols();

	if ((pShape[0] != vU[0].rows()) || (pShape[1] != vU[1].rows()) || (pShape[2] != vU[2].rows()) ||
		(pShape[3] != vU[3].rows()) || (pShape[4] != vU[4].rows()) || (pShape[5] != vU[5].rows()))
	{
		std::array<szType, 6> newShape = { vU[0].rows(), vU[1].rows(), vU[2].rows(), vU[3].rows(), vU[4].rows(), vU[5].rows() };
		coreTensor.resize(newShape);
	}

	for (idxType i6 = 0; i6 < pShape[5]; ++i6)
	{
		for (idxType i5 = 0; i5 < pShape[4]; ++i5)
		{
			for (idxType i4 = 0; i4 < pShape[3]; ++i4)
			{
				for (idxType i3 = 0; i3 < pShape[2]; ++i3)
				{
					for (idxType i2 = 0; i2 < pShape[1]; ++i2)
					{
						for (idxType i1 = 0; i1 < pShape[0]; ++i1)
						{
							T accum = 0;
							for (size_t s = 0; s < iSparsity; ++s)
								accum += val(s) * vU[0](i1, loc(0, s)) * vU[1](i2, loc(1, s)) * vU[2](i3, loc(2, s)) * vU[3](i4, loc(3, s)) * vU[4](i5, loc(4, s)) * vU[5](i6, loc(5, s));
							coreTensor[i1][i2][i3][i4][i5][i6] = accum;
						}
					}
				}
			}
		}
	}
}



template <typename T>
void HOSVDReconR1(T singVal, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU, boost::multi_array<T, 3>& coreTensor)
{
	typedef typename boost::multi_array<T, 3>::size_type szType;
	typedef typename boost::multi_array<T, 3>::index idxType;

	if (vU.size() != 3)
		return;

	const szType* pShape = coreTensor.shape();
	if ((pShape[0] != vU[0].size()) || (pShape[1] != vU[1].size()) || (pShape[2] != vU[2].size()))
	{
		std::array<szType, 3> newShape = { vU[0].size(), vU[1].size(), vU[2].size() };
		coreTensor.resize(newShape);
	}

	for (idxType i3 = 0; i3 < pShape[2]; ++i3)
		for (idxType i2 = 0; i2 < pShape[1]; ++i2)
			for (idxType i1 = 0; i1 < pShape[0]; ++i1)
				coreTensor[i1][i2][i3] = singVal * (vU[0])(i1) * (vU[1])(i2) * (vU[2])(i3);
}

template <typename T>
void HOSVDReconR1(T singVal, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU, boost::multi_array<T, 4>& coreTensor)
{
	typedef typename boost::multi_array<T, 4>::size_type szType;
	typedef typename boost::multi_array<T, 4>::index idxType;

	if (vU.size() != 4)
		return;

	const szType* pShape = coreTensor.shape();
	if ((pShape[0] != vU[0].size()) || (pShape[1] != vU[1].size()) || (pShape[2] != vU[2].size()) ||
		(pShape[3] != vU[3].size()))
	{
		std::array<szType, 4> newShape = { vU[0].size(), vU[1].size(), vU[2].size() , vU[3].size() };
		coreTensor.resize(newShape);
	}

	for (idxType i4 = 0; i4 < pShape[3]; ++i4)
		for (idxType i3 = 0; i3 < pShape[2]; ++i3)
			for (idxType i2 = 0; i2 < pShape[1]; ++i2)
				for (idxType i1 = 0; i1 < pShape[0]; ++i1)
					coreTensor[i1][i2][i3][i4] = singVal * (vU[0])(i1) * (vU[1])(i2) * (vU[2])(i3) * (vU[3])(i4);
}


template <typename T>
void HOSVDReconR1(T singVal, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU, boost::multi_array<T, 5>& coreTensor)
{
	typedef typename boost::multi_array<T, 5>::size_type szType;
	typedef typename boost::multi_array<T, 5>::index idxType;

	if (vU.size() != 5)
		return;

	const szType* pShape = coreTensor.shape();
	if ((pShape[0] != vU[0].size()) || (pShape[1] != vU[1].size()) || (pShape[2] != vU[2].size()) ||
		(pShape[3] != vU[3].size()) || (pShape[4] != vU[4].size()))
	{
		std::array<szType, 5> newShape = { vU[0].size(), vU[1].size(), vU[2].size(), vU[3].size(), vU[4].size() };
		coreTensor.resize(newShape);
	}

	for (idxType i5 = 0; i5 < pShape[4]; ++i5)
		for (idxType i4 = 0; i4 < pShape[3]; ++i4)
			for (idxType i3 = 0; i3 < pShape[2]; ++i3)
				for (idxType i2 = 0; i2 < pShape[1]; ++i2)
					for (idxType i1 = 0; i1 < pShape[0]; ++i1)
						coreTensor[i1][i2][i3][i4][i5] = singVal * (vU[0])(i1) * (vU[1])(i2) * (vU[2])(i3) * (vU[3])(i4) * (vU[4])(i5);
}


template <typename T>
void HOSVDReconR1(T singVal, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& vU, boost::multi_array<T, 6>& coreTensor)
{
	typedef typename boost::multi_array<T, 6>::size_type szType;
	typedef typename boost::multi_array<T, 6>::index idxType;

	if (vU.size() != 6)
		return;

	const szType* pShape = coreTensor.shape();
	if ((pShape[0] != vU[0].size()) || (pShape[1] != vU[1].size()) || (pShape[2] != vU[2].size()) ||
		(pShape[3] != vU[3].size()) || (pShape[4] != vU[4].size()) || (pShape[5] != vU[5].size()))
	{
		std::array<szType, 6> newShape = { vU[0].size(), vU[1].size(), vU[2].size(), vU[3].size(), vU[4].size(), vU[5].size() };
		coreTensor.resize(newShape);
	}

	for (idxType i6 = 0; i6 < pShape[5]; ++i6)
		for (idxType i5 = 0; i5 < pShape[4]; ++i5)
			for (idxType i4 = 0; i4 < pShape[3]; ++i4)
				for (idxType i3 = 0; i3 < pShape[2]; ++i3)
					for (idxType i2 = 0; i2 < pShape[1]; ++i2)
						for (idxType i1 = 0; i1 < pShape[0]; ++i1)
							coreTensor[i1][i2][i3][i4][i5][i6] = singVal * (vU[0])(i1) * (vU[1])(i2) * (vU[2])(i3) * (vU[3])(i4) * (vU[4])(i5) * (vU[5])(i6);
}


// ***tmp***
// template <typename T>
// void KronSelective_tmp1(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& matrices,
// 	const std::vector<size_t>& indices,
// 	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRes)
// {
// 	std::vector<size_t> vRowMult(indices.size());
// 	std::vector<size_t> vColMult(indices.size());
// 	for (size_t i = 0; i < indices.size(); ++i)
// 	{
// 		vRowMult[i] = matrices[indices[i]].rows();
// 		vColMult[i] = matrices[indices[i]].cols();
// 	}
// 	for (int i = vRowMult.size() - 2; i >= 0; i--)
// 		vRowMult[i] *= vRowMult[i + 1];
// 	for (int i = vColMult.size() - 2; i >= 0; i--)
// 		vColMult[i] *= vColMult[i + 1];
// 	size_t iNumRows = vRowMult[0];
// 	size_t iNumCols = vColMult[0];
// 
// 	mRes.setOnes(iNumRows, iNumCols);
// 
// 	for (size_t m = 0; m < iNumRows; ++m)
// 	{
// 		for (size_t n = 0; n < iNumCols; ++n)
// 		{
// 			for (size_t p = 0; p < indices.size(); ++p)
// 			{
// 				size_t i_p = (m / (p == (indices.size() - 1) ? 1 : vRowMult[p + 1])) % matrices[indices[p]].rows();
// 				size_t j_p = (n / (p == (indices.size() - 1) ? 1 : vColMult[p + 1])) % matrices[indices[p]].cols();
// 				mRes(m, n) *= matrices[indices[p]](i_p, j_p);
// 			}
// 		}
// 	}
// }
// 
// 
// ***tmp***
// template <typename T>
// void KronSelective_tmp2(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& matrices,
// 	const std::vector<size_t>& indices,
// 	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRes)
// {
// 	std::vector<size_t> vRowMultInv(indices.size());
// 	std::vector<size_t> vColMultInv(indices.size());
// 	std::vector<size_t> vRowMult(indices.size());
// 	std::vector<size_t> vColMult(indices.size());
// 	for (size_t i = 0; i < indices.size(); ++i)
// 	{
// 		vRowMult[i] = matrices[indices[i]].rows();
// 		vColMult[i] = matrices[indices[i]].cols();
// 		vRowMultInv[i] = matrices[indices[i]].rows();
// 		vColMultInv[i] = matrices[indices[i]].cols();
// 	}
// 	for (int i = vRowMultInv.size() - 2; i >= 0; i--)
// 		vRowMultInv[i] *= vRowMultInv[i + 1];
// 	for (int i = vColMultInv.size() - 2; i >= 0; i--)
// 		vColMultInv[i] *= vColMultInv[i + 1];
// 	for (int i = 1; i < vRowMult.size(); i++)
// 		vRowMult[i] *= vRowMult[i - 1];
// 	for (int i = 1; i < vColMult.size(); i++)
// 		vColMult[i] *= vColMult[i - 1];
// 	size_t iNumRows = vRowMultInv[0];
// 	size_t iNumCols = vColMultInv[0];
// 
// 	mRes.setOnes(iNumRows, iNumCols);
// 
// 	for (size_t p = 0; p < indices.size(); ++p)
// 	{
// 		size_t iRowExtent = p == (indices.size() - 1) ? 1 : vRowMultInv[p + 1];
// 		size_t iColExtent = p == (indices.size() - 1) ? 1 : vColMultInv[p + 1];
// 		size_t iRows = matrices[indices[p]].rows();
// 		size_t iCols = matrices[indices[p]].cols();
// 
// 		for (size_t i = 0; i < vRowMult[p]; ++i)
// 		{
// 			size_t i_p = i % iRows;
// 			for (size_t j = 0; j < vColMult[p]; ++j)
// 			{
// 				mRes.block(i*iRowExtent, j*iColExtent, iRowExtent, iColExtent) *= matrices[indices[p]](i_p, j % iCols);
// 			}
// 		}
// 	}
// }