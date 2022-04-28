





template <typename T>
int64_t NumOfPatches(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mIm,
	const Eigen::Vector2i& vPatchSize,
	PatchType ePatchType,
	const std::vector<uint32_t>& vSlidingDis)
{
	if (ePatchType == CCS_PATCHTYPE_NOV)
	{
		size_t m = mIm.rows();
		size_t n = mIm.cols();

		if (m == 0 || n == 0 || mIm.rows() < vPatchSize(0) || mIm.cols() < vPatchSize(1))
			return -1;

		size_t mpad = m % vPatchSize(0);
		if (mpad > 0)
			mpad = vPatchSize(0) - mpad;

		size_t npad = n % vPatchSize(1);
		if (npad > 0)
			npad = vPatchSize(1) - npad;

		return ((m + mpad) / vPatchSize(0)) * ((n + npad) / vPatchSize(1));
	}
	else if (ePatchType == CCS_PATCHTYPE_OV)
	{
		size_t mm = mIm.rows();
		size_t nn = mIm.cols();
		size_t m = vPatchSize(0);
		size_t n = vPatchSize(1);
		if (m <= 0 || n <= 0 || mm < m || nn < n || vSlidingDis.size() < 2 || vSlidingDis[0] == 0 || vSlidingDis[1] == 0 || vSlidingDis[0] > m || vSlidingDis[1] > n)
			return -1;

		size_t iBlocksM = mm - m + 1;
		size_t iBlocksN = nn - n + 1;

		Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mIdx; mIdx.setConstant(iBlocksM, iBlocksN, 0);
		for (size_t i = 0; i < mIdx.rows(); i += vSlidingDis[0])
			for (size_t j = 0; j < mIdx.cols(); j += vSlidingDis[1])
				mIdx(i, j) = 1;
		for (size_t i = 0; i < mIdx.rows(); i += vSlidingDis[0])
			mIdx(i, mIdx.cols() - 1) = 1;
		for (size_t i = 0; i < mIdx.cols(); i += vSlidingDis[1])
			mIdx(mIdx.rows() - 1, i) = 1;
		mIdx(mIdx.rows() - 1, mIdx.cols() - 1) = 1;

		return (mIdx.array() > 0).count();
	}
	else
		return -1;
}



int64_t NumOfPatches3D(const Eigen::Vector3i& vImSize,
	const Eigen::Vector3i& vPatchSize,
	PatchType ePatchType,
	const std::vector<uint32_t>& vSlidingDis)
{
	if (ePatchType == CCS_PATCHTYPE_NOV)
	{
		size_t m = vImSize(0);
		size_t n = vImSize(1);
		size_t p = vImSize(2);

		if (m == 0 || n == 0 || p == 0 || m < vPatchSize(0) || n < vPatchSize(1) || p < vPatchSize(2))
			return -1;

		size_t mpad = m % vPatchSize(0);
		if (mpad > 0)
			mpad = vPatchSize(0) - mpad;

		size_t npad = n % vPatchSize(1);
		if (npad > 0)
			npad = vPatchSize(1) - npad;

		size_t ppad = p % vPatchSize(2);
		if (ppad > 0)
			ppad = vPatchSize(2) - ppad;

		return ((m + mpad) / vPatchSize(0)) * ((n + npad) / vPatchSize(1)) * ((p + ppad) / vPatchSize(2));
	}
	else if (ePatchType == CCS_PATCHTYPE_OV)
	{
		size_t mm = vImSize(0);
		size_t nn = vImSize(1);
		size_t pp = vImSize(2);
		size_t m = vPatchSize(0);
		size_t n = vPatchSize(1);
		size_t p = vPatchSize(2);
		if (m <= 0 || n <= 0 || p <= 0 || mm < m || nn < n || pp < p || vSlidingDis.size() < 3 || vSlidingDis[0] == 0 || vSlidingDis[1] == 0 || vSlidingDis[2] == 0 ||
			vSlidingDis[0] > m || vSlidingDis[1] > n || vSlidingDis[2] > p)
			return -1;

		size_t iBlocksM = mm - m + 1;
		size_t iBlocksN = nn - n + 1;
		size_t iBlocksP = pp - p + 1;

		std::vector<size_t> vTmp = { iBlocksM, iBlocksN, iBlocksP };
		boost::multi_array<int, 3> mIdx(vTmp, CCS_TENSOR_STORAGE_ORDER);
		for (size_t i = 0; i < mIdx.num_elements(); ++i)
			mIdx.data()[i] = 0;
		for (size_t i = 0; i < vTmp[0]; i += vSlidingDis[0])
			for (size_t j = 0; j < vTmp[1]; j += vSlidingDis[1])
				for (size_t k = 0; k < vTmp[2]; k += vSlidingDis[2])
					mIdx[i][j][k] = 1;

		size_t iCount = 0;
		for (size_t i = 0; i < mIdx.num_elements(); ++i)
			if (mIdx.data()[i])
				iCount++;

		for (size_t i = 0; i < vTmp[0]; i += vSlidingDis[0])
			for (size_t j = 0; j < vTmp[1]; j += vSlidingDis[1])
				mIdx[i][j][vTmp[2] - 1] = 1;
		for (size_t i = 0; i < vTmp[0]; i += vSlidingDis[0])
			for (size_t k = 0; k < vTmp[2]; k += vSlidingDis[2])
				mIdx[i][vTmp[1] - 1][k] = 1;
		for (size_t j = 0; j < vTmp[1]; j += vSlidingDis[1])
			for (size_t k = 0; k < vTmp[2]; k += vSlidingDis[2])
				mIdx[vTmp[0] - 1][j][k] = 1;

		for (size_t i = 0; i < vTmp[0]; i += vSlidingDis[0])
			mIdx[i][vTmp[1] - 1][vTmp[2] - 1] = 1;
		for (size_t j = 0; j < vTmp[1]; j += vSlidingDis[1])
			mIdx[vTmp[0] - 1][j][vTmp[2] - 1] = 1;
		for (size_t k = 0; k < vTmp[2]; k += vSlidingDis[2])
			mIdx[vTmp[0] - 1][vTmp[1] - 1][k] = 1;

		mIdx[vTmp[0] - 1][vTmp[1] - 1][vTmp[2] - 1] = 1;

		iCount = 0;
		for (size_t i = 0; i < mIdx.num_elements(); ++i)
			if (mIdx.data()[i])
				iCount++;
		int x = 0;
		return iCount;
	}
	else
		return -1;

}


//Note that image patches are ordered column-wise
template <typename T>
bool im2col(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mIm,
	const Eigen::Vector2i& vPatchSize,
	PatchType ePatchType,
	const std::vector<uint32_t>& vSlidingDis,
	std::vector<Eigen::Vector2i>& vOutSldPatches,
	T tPadVal,
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mOut,
	size_t& mblocks, size_t& nblocks)
{
	if (ePatchType == CCS_PATCHTYPE_NOV)
	{
		size_t m = mIm.rows();
		size_t n = mIm.cols();

		if (m == 0 || n == 0 || mIm.rows() < vPatchSize(0) || mIm.cols() < vPatchSize(1))
		{
			mOut.resize(0, 0);
			return false;
		}

		size_t mpad = m % vPatchSize(0);
		if (mpad > 0)
			mpad = vPatchSize(0) - mpad;

		size_t npad = n % vPatchSize(1);
		if (npad > 0)
			npad = vPatchSize(1) - npad;

		mblocks = (m + mpad) / vPatchSize(0);
		nblocks = (n + npad) / vPatchSize(1);

		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> aa; aa.setConstant(m + mpad, n + npad, tPadVal);
		aa.topLeftCorner(m, n) = mIm;

		mOut.setZero(vPatchSize.prod(), mblocks*nblocks);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mPatch; mPatch.resize(vPatchSize(0), vPatchSize(1));

		for (size_t i = 0; i < mblocks; ++i)
		{
			for (size_t j = 0; j < nblocks; ++j)
			{
				size_t iIdx = i + j*mblocks;
				mPatch = aa.block(i*vPatchSize(0), j*vPatchSize(1), vPatchSize(0), vPatchSize(1));
				for (size_t k = 0; k < mPatch.cols(); k++)
					for (size_t l = 0; l < mPatch.rows(); ++l)
						mOut(l + k*mPatch.rows(), iIdx) = mPatch(l, k);
			}
		}
	}
	else if (ePatchType == CCS_PATCHTYPE_OV)
	{
		size_t mm = mIm.rows();
		size_t nn = mIm.cols();
		size_t m = vPatchSize(0);
		size_t n = vPatchSize(1);

		if (m <= 0 || n <= 0 || mm < m || nn < n || vSlidingDis.size() < 2 || vSlidingDis[0] == 0 || vSlidingDis[1] == 0 || vSlidingDis[0] > m || vSlidingDis[1] > n)
		{
			mOut.resize(0, 0);
			return false;
		}

		mblocks = mm - m + 1;
		nblocks = nn - n + 1;

		Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mIdx; mIdx.setConstant(mblocks, nblocks, 0);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mPatch; mPatch.resize(m, n);

		for (size_t i = 0; i < mIdx.rows(); i += vSlidingDis[0])
			for (size_t j = 0; j < mIdx.cols(); j += vSlidingDis[1])
				mIdx(i, j) = 1;
		for (size_t i = 0; i < mIdx.rows(); i += vSlidingDis[0])
			mIdx(i, mIdx.cols() - 1) = 1;
		for (size_t i = 0; i < mIdx.cols(); i += vSlidingDis[1])
			mIdx(mIdx.rows() - 1, i) = 1;
		mIdx(mIdx.rows() - 1, mIdx.cols() - 1) = 1;

		mblocks = (mIdx.row(0).array() > 0).count();
		nblocks = (mIdx.col(0).array() > 0).count();

		size_t iNumPatches = (mIdx.array() > 0).count();
		mOut.resize(m*n, iNumPatches);
		vOutSldPatches.resize(iNumPatches);
		size_t iCount = 0;
		for (size_t i = 0; i < mIdx.cols(); ++i)
		{
			for (size_t j = 0; j < mIdx.rows(); ++j)
			{
				size_t iIdx = j + i*mIdx.rows();
				if (mIdx(j, i))
				{
					mPatch = mIm.block(j, i, m, n);
					for (size_t k = 0; k < mPatch.cols(); k++)
						for (size_t l = 0; l < mPatch.rows(); ++l)
							mOut(l + k*mPatch.rows(), iCount) = mPatch(l, k);
					vOutSldPatches[iCount] = Eigen::Vector2i(j, i);
					iCount++;
				}
			}
		}
	}
	else
	{
		std::cerr << "ERROR: Wrong patch type in im2col()." << std::endl;
		return false;
	}

	return true;
}


template <typename T>
bool col2im(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mCol, 
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mOut, 
            const Eigen::Vector2i& vPatchSize, 
            const Eigen::Vector2i& vImSize, 
            PatchType ePatchType,
            const std::vector<Eigen::Vector2i>& vInSldPatches)
{
	size_t m = vPatchSize(0);
	size_t n = vPatchSize(1);
	if (vPatchSize.prod() != mCol.rows())
	{
		mOut.resize(0, 0);
		return false;
	}

	if (ePatchType == CCS_PATCHTYPE_NOV)
	{
		size_t mpad = vImSize(0) % m;
		if (mpad > 0)
			mpad = vPatchSize(0) - mpad;

		size_t npad = vImSize(1) % n;
		if (npad > 0)
			npad = vPatchSize(1) - npad;

		mpad = vImSize(0) + mpad;
		npad = vImSize(1) + npad;

		size_t mblocks = mpad / m;
		size_t nblocks = npad / n;
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> aa; aa.resize(mpad, npad);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mPatch; mPatch.resize(m, n);

		for (size_t i = 0; i < mblocks; ++i)
		{
			for (size_t j = 0; j < nblocks; ++j)
			{
				size_t iIdx = i + j*mblocks;
				for (size_t k = 0; k < mPatch.cols(); k++)
					for (size_t l = 0; l < mPatch.rows(); ++l)
						mPatch(l,k) = mCol(l+k*mPatch.rows(), iIdx);
				aa.block(i*m, j*n, m, n) = mPatch;
			}
		}
		mOut.resize(vImSize(0), vImSize(1));
		mOut = aa.topLeftCorner(vImSize(0), vImSize(1));
	}
	else if (ePatchType == CCS_PATCHTYPE_OV)
	{
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mOutDouble;
		mOutDouble.setConstant(vImSize(0), vImSize(1), 0.0);
		mOut.setConstant(vImSize(0), vImSize(1), 0.0);
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mWt; mWt.setConstant(vImSize(0), vImSize(1), 0.0);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mPatch; mPatch.resize(m,n);
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mPatchDouble; mPatchDouble.resize(m, n);
		size_t iBlocksM = vImSize(0) - m + 1;
		size_t iBlocksN = vImSize(1) - n + 1;
		for (size_t i = 0; i < vInSldPatches.size(); ++i)
		{
			for (size_t k = 0; k < mPatch.cols(); k++)
				for (size_t l = 0; l < mPatch.rows(); ++l)
					mPatch(l, k) = mCol(l + k*mPatch.rows(), i);
			CastDataPointTTypeToInternal(mPatch, mPatchDouble);
			mOutDouble.block(vInSldPatches[i](0), vInSldPatches[i](1), m, n) += mPatchDouble;
			mWt.block(vInSldPatches[i](0), vInSldPatches[i](1), m, n) += Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Constant(m, n, 1.0);
		}
		mOutDouble = mOutDouble.array() / mWt.array();
		CastDataPointInternalToTType(mOutDouble, mOut);
	}
	else
	{
		std::cerr << "ERROR: Wrong patch type in col2im()." << std::endl;
		return false;
	}

	return true;
}





//Note that image patches are ordered column-wise
template <typename T>
bool im2col3D(const boost::multi_array<T, 3>& mIm,
	const Eigen::Vector3i& vPatchSize,
	PatchType ePatchType,
	const std::vector<uint32_t>& vSlidingDis,
	std::vector<Eigen::Vector3i>& vOutSldPatches,
	T tPadVal,
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mOut,
	size_t& mblocks, size_t& nblocks, size_t& pblocks)
{
	typedef typename boost::multi_array<T, 3>::size_type szType;
	typedef typename boost::multi_array<T, 3>::index idxType;

	std::vector<int> vvPatchSize = { vPatchSize(0), vPatchSize(1), vPatchSize(2) };
	boost::multi_array<T, 3> mPatch(vvPatchSize, CCS_TENSOR_STORAGE_ORDER);

	if (ePatchType == CCS_PATCHTYPE_NOV)
	{
		const szType* shape = mIm.shape();
		size_t m = shape[0];
		size_t n = shape[1];
		size_t p = shape[2];
		if (m == 0 || n == 0 || p == 0 || m < vPatchSize(0) || n < vPatchSize(1) || p < vPatchSize(2))
		{
			mOut.resize(0, 0);
			return false;
		}

		size_t mpad = m % vPatchSize(0);
		if (mpad > 0)
			mpad = vPatchSize(0) - mpad;

		size_t npad = n % vPatchSize(1);
		if (npad > 0)
			npad = vPatchSize(1) - npad;

		size_t ppad = p % vPatchSize(2);
		if (ppad > 0)
			ppad = vPatchSize(2) - ppad;

		mblocks = (m + mpad) / vPatchSize(0);
		nblocks = (n + npad) / vPatchSize(1);
		pblocks = (p + ppad) / vPatchSize(2);

		std::vector<szType> vTmp = { (m + mpad) , (n + npad) , (p + ppad) };
		boost::multi_array<T, 3> aa(vTmp, CCS_TENSOR_STORAGE_ORDER); 
		for (size_t i = 0; i < aa.num_elements(); ++i)
			aa.data()[i] = tPadVal;
		for (size_t t0 = 0; t0 < m; ++t0)
			for (size_t t1 = 0; t1 < n; ++t1)
				for (size_t t2 = 0; t2 < p; ++t2)
					aa[t0][t1][t2] = mIm[t0][t1][t2];

		mOut.setZero(vPatchSize.prod(), mblocks*nblocks*pblocks);

		for (size_t i = 0; i < mblocks; ++i)
		{
			for (size_t j = 0; j < nblocks; ++j)
			{
				for (size_t k = 0; k < pblocks; ++k)
				{
					size_t iIdx = i + j*mblocks + k*nblocks*mblocks;
					for (size_t t0 = 0; t0 < vPatchSize(0); ++t0)
						for (size_t t1 = 0; t1 < vPatchSize(1); ++t1)
							for (size_t t2 = 0; t2 < vPatchSize(2); ++t2)
								mPatch[t0][t1][t2] = aa[t0 + i*vPatchSize(0)][t1 + j*vPatchSize(1)][t2 + k*vPatchSize(2)];
					for (size_t t = 0; t < vPatchSize.prod(); ++t)
						mOut(t, iIdx) = mPatch.data()[t];
				}
			}
		}
	}
	else if (ePatchType == CCS_PATCHTYPE_OV)
	{
		const szType* shape = mIm.shape();
		size_t mm = shape[0];
		size_t nn = shape[1];
		size_t pp = shape[2];
		size_t m = vPatchSize(0);
		size_t n = vPatchSize(1);
		size_t p = vPatchSize(2);
		if (m <= 0 || n <= 0 || p <= 0 || mm < m || nn < n || pp < p || vSlidingDis.size() < 3 || vSlidingDis[0] == 0 || vSlidingDis[1] == 0 || vSlidingDis[2] == 0 ||
			vSlidingDis[0] > m || vSlidingDis[1] > n || vSlidingDis[2] > p)
		{
			mOut.resize(0, 0);
			return false;
		}
		mblocks = mm - m + 1;
		nblocks = nn - n + 1;
		pblocks = pp - p + 1;

		std::vector<size_t> vTmp = { mblocks, nblocks, pblocks };
		boost::multi_array<int, 3> mIdx(vTmp, CCS_TENSOR_STORAGE_ORDER);
		for (size_t i = 0; i < mIdx.num_elements(); ++i)
			mIdx.data()[i] = 0;
		for (size_t i = 0; i < vTmp[0]; i += vSlidingDis[0])
			for (size_t j = 0; j < vTmp[1]; j += vSlidingDis[1])
				for (size_t k = 0; k < vTmp[2]; k += vSlidingDis[2])
					mIdx[i][j][k] = 1;

		for (size_t i = 0; i < vTmp[0]; i += vSlidingDis[0])
			for (size_t j = 0; j < vTmp[1]; j += vSlidingDis[1])
				mIdx[i][j][vTmp[2] - 1] = 1;
		for (size_t i = 0; i < vTmp[0]; i += vSlidingDis[0])
			for (size_t k = 0; k < vTmp[2]; k += vSlidingDis[2])
				mIdx[i][vTmp[1] - 1][k] = 1;
		for (size_t j = 0; j < vTmp[1]; j += vSlidingDis[1])
			for (size_t k = 0; k < vTmp[2]; k += vSlidingDis[2])
				mIdx[vTmp[0] - 1][j][k] = 1;

		for (size_t i = 0; i < vTmp[0]; i += vSlidingDis[0])
			mIdx[i][vTmp[1] - 1][vTmp[2] - 1] = 1;
		for (size_t j = 0; j < vTmp[1]; j += vSlidingDis[1])
			mIdx[vTmp[0] - 1][j][vTmp[2] - 1] = 1;
		for (size_t k = 0; k < vTmp[2]; k += vSlidingDis[2])
			mIdx[vTmp[0] - 1][vTmp[1] - 1][k] = 1;

		mIdx[vTmp[0] - 1][vTmp[1] - 1][vTmp[2] - 1] = 1;

		mblocks = 0; nblocks = 0; pblocks = 0;
		for (size_t i = 0; i < vTmp[0]; ++i)
			if (mIdx[i][0][0] > 0)
				mblocks++;
		for (size_t j = 0; j < vTmp[1]; ++j)
			if (mIdx[0][j][0] > 0)
				nblocks++;
		for (size_t k = 0; k < vTmp[2]; ++k)
			if (mIdx[0][0][k] > 0)
				pblocks++;

		size_t iNumPatches = 0;
		for (size_t i = 0; i < mIdx.num_elements(); ++i)
			if (mIdx.data()[i])
				iNumPatches++;
		mOut.resize(m*n*p, iNumPatches);
		vOutSldPatches.resize(iNumPatches);
		size_t iCount = 0;
		for (size_t i = 0; i < mIdx.shape()[0]; ++i)
		{
			for (size_t j = 0; j < mIdx.shape()[1]; ++j)
			{
				for (size_t k = 0; k < mIdx.shape()[2]; ++k)
				{
					size_t iIdx = i + j*mIdx.shape()[0] + k*mIdx.shape()[1] * mIdx.shape()[0];
					if (mIdx[i][j][k])
					{
						for (size_t t0 = 0; t0 < m; ++t0)
							for (size_t t1 = 0; t1 < n; ++t1)
								for (size_t t2 = 0; t2 < p; ++t2)
									mPatch[t0][t1][t2] = mIm[t0 + i][t1 + j][t2 + k];
						for (size_t t = 0; t < vPatchSize.prod(); ++t)
							mOut(t, iCount) = mPatch.data()[t];
						vOutSldPatches[iCount] = Eigen::Vector3i(i, j, k);
						iCount++;
					}
				}
			}
		}
	}
	else
	{
		std::cerr << "ERROR: Wrong patch type in im2col()." << std::endl;
		return false;
	}

	return true;
}



template <typename T>
bool col2im3D(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mCol,
	boost::multi_array<T, 3>& mOut,
	const Eigen::Vector3i& vPatchSize,
	const Eigen::Vector3i& vImSize,
	PatchType ePatchType,
	const std::vector<Eigen::Vector3i>& vInSldPatches)
{
	size_t m = vPatchSize(0);
	size_t n = vPatchSize(1);
	size_t p = vPatchSize(2);
	if (vPatchSize.prod() != mCol.rows())
		return false;

	std::vector<int> vvPatchSize = { vPatchSize(0), vPatchSize(1), vPatchSize(2) };
	boost::multi_array<T, 3> mPatch(vvPatchSize, CCS_TENSOR_STORAGE_ORDER);
	std::vector<int> vvImSize = { vImSize(0), vImSize(1), vImSize(2) };
	mOut = boost::multi_array<T, 3>(vvImSize, CCS_TENSOR_STORAGE_ORDER);
	for (size_t i = 0; i < mOut.num_elements(); ++i)
		mOut.data()[i] = 0;
	

	if (ePatchType == CCS_PATCHTYPE_NOV)
	{
		size_t mpad = vImSize(0) % m;
		if (mpad > 0)
			mpad = vPatchSize(0) - mpad;

		size_t npad = vImSize(1) % n;
		if (npad > 0)
			npad = vPatchSize(1) - npad;

		size_t ppad = vImSize(2) % p;
		if (ppad > 0)
			ppad = vPatchSize(2) - ppad;

		mpad = vImSize(0) + mpad;
		npad = vImSize(1) + npad;
		ppad = vImSize(2) + ppad;

		size_t mblocks = mpad / m;
		size_t nblocks = npad / n;
		size_t pblocks = ppad / p;
		std::vector<size_t> aaSize = { mpad, npad, ppad };
		boost::multi_array<T, 3> aa(aaSize, CCS_TENSOR_STORAGE_ORDER);

		for (size_t i = 0; i < mblocks; ++i)
		{
			for (size_t j = 0; j < nblocks; ++j)
			{
				for (size_t k = 0; k < pblocks; ++k)
				{
					size_t iIdx = i + j*mblocks + k*nblocks*mblocks;
					for (size_t t = 0; t < vPatchSize.prod(); ++t)
						mPatch.data()[t] = mCol(t, iIdx);
					for (size_t t0 = 0; t0 < vPatchSize(0); ++t0)
						for (size_t t1 = 0; t1 < vPatchSize(1); ++t1)
							for (size_t t2 = 0; t2 < vPatchSize(2); ++t2)
								aa[t0 + i*m][t1 + j*n][t2 + k*p] = mPatch[t0][t1][t2];
				}
			}
		}
		for (size_t t0 = 0; t0 < mOut.shape()[0]; ++t0)
			for (size_t t1 = 0; t1 < mOut.shape()[1]; ++t1)
				for (size_t t2 = 0; t2 < mOut.shape()[2]; ++t2)
					mOut[t0][t1][t2] = aa[t0][t1][t2];
	}
	else if (ePatchType == CCS_PATCHTYPE_OV)
	{
		boost::multi_array<size_t, 3> mWt(vvImSize, CCS_TENSOR_STORAGE_ORDER);
		for (size_t i = 0; i < mWt.num_elements(); ++i)
			mWt.data()[i] = 0;
		size_t iBlocksM = vImSize(0) - m + 1;
		size_t iBlocksN = vImSize(1) - n + 1;
		size_t iBlocksP = vImSize(2) - p + 1;
		for (size_t i = 0; i < vInSldPatches.size(); ++i)
		{
			for (size_t t0 = 0; t0 < m; ++t0)
			{
				for (size_t t1 = 0; t1 < n; ++t1)
				{
					for (size_t t2 = 0; t2 < p; ++t2)
					{
						//mPatch[t0][t1][t2] = mCol(t0 + t1*m + t2*m*n, i);
						size_t iIdx = t0 + t1*m + t2*m*n;
						mOut[vInSldPatches[i](0) + t0][vInSldPatches[i](1) + t1][vInSldPatches[i](2) + t2] += mCol(iIdx, i);
						mWt[vInSldPatches[i](0) + t0][vInSldPatches[i](1) + t1][vInSldPatches[i](2) + t2] += 1;
					}
				}
			}
		}
		for(size_t t = 0; t < vImSize.prod(); ++t)
			mOut.data()[t] /= T(mWt.data()[t]);
	}
	else
	{
		std::cerr << "ERROR: Wrong patch type in col2im()." << std::endl;
		return false;
	}

	return true;
}