




template <typename T, size_t N>
CDictTestEnsOrthnD<T, N>::CDictTestEnsOrthnD(const std::string& strCoeffsInfoFileName/* = "SInfo"*/,
	const std::string& strVarMemb/* = "M"*/,
	const std::string& strVarCoeffVal/* = "NZ"*/,
	const std::string& strVarCoeffLoc/* = "LNZ"*/,
	const std::string& strVarNNZ/* = "NNZ"*/,
	const std::string& strWeightCos, 
	const std::string& strWeightCosTi) : CDictTest<T, N>(strCoeffsInfoFileName, strVarMemb, strVarCoeffVal, strVarCoeffLoc, strVarNNZ)
{
	this->m_bStoreCoeffs = false;
	m_strWeightCos = strWeightCos;
	m_strWeightCosTi = strWeightCosTi;
}


template <typename T, size_t N>
CDictTestEnsOrthnD<T, N>::~CDictTestEnsOrthnD()
{
	CleanUp();
}


template <typename T, size_t N>
bool CDictTestEnsOrthnD<T, N>::Init(bool bStoreCoeffs)
{
	typedef typename CData<T, N>::InternalArraynD InternalTensor;
	CleanUp();

	this->m_bStoreCoeffs = bStoreCoeffs;

// 	std::vector<size_t> vTmpSize = { 1,1,1 };
// 	m_WeightCos = typename CData<T, N>::InternalArraynD(vTmpSize, CCS_TENSOR_STORAGE_ORDER);
// 	m_WeightCosTi = typename CData<T, N>::InternalArraynD(vTmpSize, CCS_TENSOR_STORAGE_ORDER);

	VclMatio readFile;
	readFile.openForReading(m_strWeightCos);
	readFile.readEigenMatrixndNamed("weight_matrix", m_WeightCos);
	readFile.close();

	readFile.openForReading(m_strWeightCosTi);
	readFile.readEigenMatrixndNamed("weight_matrix", m_WeightCosTi);
	readFile.close();

	return true;
}


template <typename T, size_t N>
bool CDictTestEnsOrthnD<T, N>::Test(const CData<T, N>* pData, const CDictionary* pDict, uint32_t iSparsity, CCS_INTERNAL_TYPE tThreshold, double& dTimeDelta)
{
	typedef typename CData<T, N>::InternalArray1D InternalVec;
	typedef typename CData<T, N>::InternalArray2D InternalMat;
	typedef typename CData<T, N>::InternalArraynD InternalTensor;

	assert(pData);
	assert(pDict);

	if (!pData || !pDict)
		return false;

	int iNumThreads = omp_get_max_threads();

	typedef typename CData<T, N>::ArraynD::size_type szType;
	size_t iNumDataPoints = pData->GetNumDataElems();
	pData->GetDataElemDim(this->m_vShape);
	szType iLength = std::accumulate(this->m_vShape.begin(), this->m_vShape.end(), szType(1), std::multiplies<szType>());
	size_t iNumDims = this->m_vShape.size();
	size_t iNullables = iLength - iSparsity;


	if ((pData->GetDatanD() == NULL) || (pData->GetDatanD()->empty()))
	{
		std::cerr << "ERROR: Invalid data for testing." << std::endl;
		return false;
	}

	//Make sure we have a correct ensemble
	if (pDict->GetNumDict() < 1)
	{
		std::cerr << "ERROR: Invalid dictionary." << std::endl;
		return false;
	}
	for (size_t i = 0; i < pDict->GetNumDict(); ++i)
	{
		if (pDict->GetDict(i).size() != N)
		{
			std::cerr << "ERROR: Invalid dictionary." << std::endl;
			return false;
		}
	}

	for (size_t i = 0; i < pDict->GetNumDict(); ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			if ((pDict->GetDictElem(i, j).rows() != this->m_vShape[j]) || (pDict->GetDictElem(i, j).cols() != this->m_vShape[j]))
			{
				std::cerr << "ERROR: Dictionary not compatible with data." << std::endl;
				return false;
			}
		}
	}

	//Initialize M, S and NNZ
	this->m_mMemb.setZero(iNumDataPoints, 1);
	if (this->m_bStoreCoeffs)
	{
		this->m_vNZValue.resize(iNumDataPoints);
		this->m_vNZLoc.resize(iNumDataPoints);
		this->m_vNNZ.setZero(iNumDataPoints);
	}

	double dLastTime = omp_get_wtime();
	CProgressReporter reporter(iNumDataPoints, 0);
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < iNumDataPoints; ++i)
	{
		int iThreadID = omp_get_thread_num();

		std::vector<size_t>		vIndices;
		Eigen::VectorXi			vNNZ;
		InternalVec				vError;

		vIndices.resize(iLength);
		vNNZ.setZero(pDict->GetNumDict());
		vError.setZero(pDict->GetNumDict());

		//Cast data point
		InternalTensor dataPoint(this->m_vShape, CCS_TENSOR_STORAGE_ORDER);
		CastDataPointTTypeToInternal<T, N>(pData->GetDatanD()->at(i), dataPoint);
		InternalTensor dataPointWeighted(this->m_vShape, CCS_TENSOR_STORAGE_ORDER);

		//Apply weights to data point
		for (size_t j = 0; j < iLength; ++j)
			dataPointWeighted.data()[j] = dataPoint.data()[j] * m_WeightCosTi.data()[j];

		CCS_INTERNAL_TYPE dAvgErrorRef = TensorNorm2Squared(dataPointWeighted) / iLength;

		InternalTensor tmpS1(this->m_vShape, CCS_TENSOR_STORAGE_ORDER);
		InternalTensor tmpS2(this->m_vShape, CCS_TENSOR_STORAGE_ORDER);

		for (uint32_t a = 0; a < pDict->GetNumDict(); ++a)
		{
			TensorCopy(dataPointWeighted, tmpS1);
			for (size_t j = 0; j < iNumDims; ++j)
			{
				TensorProdMat(tmpS1, pDict->GetDictElem(a, j).transpose(), j, tmpS2);
				TensorCopy(tmpS2, tmpS1);
			}

			typename InternalTensor::element* pCoeffs = tmpS1.data();
			std::iota(vIndices.begin(), vIndices.end(), 0);
			std::sort(vIndices.begin(), vIndices.end(), [&pCoeffs](size_t i1, size_t i2) { return fabs(pCoeffs[i1]) > fabs(pCoeffs[i2]); });

			vError(a) = dAvgErrorRef;
			if (dAvgErrorRef < tThreshold) //when data point norm is too small
			{
				vError(a) = tThreshold;
				vNNZ(a) = vNNZ(a) + 1;
			}
			else
			{
				//Add coefficients until the threshold or maximum number of non-zero coefficients are reached
				while ((vError(a) > tThreshold) && (vNNZ(a) < iSparsity))
				{
					vError(a) = vError(a) - powf(pCoeffs[vIndices[vNNZ(a)]], 2.0) / iLength;
					vNNZ(a) = vNNZ(a) + 1;
				}
			}
		}

		Eigen::VectorXi::Index iSparsest;
		int iMinSparsity;
		iMinSparsity = vNNZ.minCoeff(&iSparsest);

		if (iMinSparsity == iSparsity)
		{
			typename InternalVec::Index iLeastError;
			vError.minCoeff(&iLeastError);
			this->m_mMemb(i, 0) = iLeastError;
			iMinSparsity = iSparsity;
		}
		else
		{
			this->m_mMemb(i, 0) = iSparsest;
		}

		if (this->m_bStoreCoeffs)
		{
			
			TensorCopy(dataPoint, tmpS1);
			for (size_t j = 0; j < iNumDims; ++j)
			{
				TensorProdMat(tmpS1, pDict->GetDictElem(this->m_mMemb(i, 0), j).transpose(), j, tmpS2);
				TensorCopy(tmpS2, tmpS1);
			}

			typename InternalTensor::element* pCoeffs = tmpS1.data();
			std::iota(vIndices.begin(), vIndices.end(), 0);
			std::sort(vIndices.begin(), vIndices.end(), [&pCoeffs](size_t i1, size_t i2) { return fabs(pCoeffs[i1]) > fabs(pCoeffs[i2]); });

			CCS_INTERNAL_TYPE dErr = TensorNorm2Squared(dataPoint) / iLength;
			size_t iNNZ = 0;
			InternalTensor Sk(this->m_vShape, CCS_TENSOR_STORAGE_ORDER);
			TensorSetZero(Sk);

			for (size_t j = 0; j < iMinSparsity; ++j)
				(Sk.data())[vIndices[j]] = pCoeffs[vIndices[j]];

			this->m_vNNZ(i) = iMinSparsity;
			this->m_vNZLoc[i].resize(this->m_vShape.size(), iMinSparsity);
			this->m_vNZValue[i].resize(iMinSparsity);
			DenseToSparse(Sk, this->m_vNZLoc[i], this->m_vNZValue[i], iMinSparsity);
		}

		reporter.Update();
	}
	reporter.Done();

	dTimeDelta = omp_get_wtime() - dLastTime;

	return true;
}


template <typename T, size_t N>
bool CDictTestEnsOrthnD<T, N>::Save(const boost::filesystem::path& strPath, const std::string& strFileName, size_t iNumQuantPoints, bool bSaveRaw, double& dTimeDelta)
{
	if (bSaveRaw)
	{
		VclMatio writeCoeffs;
		double dLastTime = omp_get_wtime();

		if (!writeCoeffs.openForWriting((strPath / fs::path(strFileName)).string()))
		{
			assert(0);
			return false;
		}
		if (!writeCoeffs.writeEigenMatrix2dNamed(this->m_strVarMemb, this->m_mMemb))
		{
			assert(0);
			return false;
		}
		if (this->m_bStoreCoeffs)
		{
			Eigen::VectorXi vDims(this->m_vShape.size());
			for (size_t i = 0; i < this->m_vShape.size(); ++i)
				vDims(i) = this->m_vShape[i];
			if (!writeCoeffs.writeEigenVectorNamed(this->m_strVarDim, vDims))
			{
				assert(0);
				return false;
			}
			//Count nonzeros
			size_t iTotalNNZ = std::accumulate(this->m_vNNZ.data(), this->m_vNNZ.data() + this->m_vNNZ.size(), 0);

			//Gather all nonzeros in pairs of location and value.
			typename CDictTest<T, N>::DataPointCoeffLocs mAllCoeffsLoc(this->m_vShape.size(), iTotalNNZ);
			typename CDictTest<T, N>::DataPointCoeffVals mAllCoeffsVal(iTotalNNZ);

			//Write nonzeros
			size_t iCounter = 0;
			for (size_t i = 0; i < this->m_vNZLoc.size(); ++i)
			{
				mAllCoeffsLoc.block(0, iCounter, this->m_vNZLoc[i].rows(), this->m_vNNZ(i)) = this->m_vNZLoc[i];
				mAllCoeffsVal.segment(iCounter, this->m_vNNZ(i)) = this->m_vNZValue[i];
				iCounter += this->m_vNNZ(i);
			}
			if (!writeCoeffs.writeEigenMatrix2dNamed(this->m_strVarCoeffLoc, mAllCoeffsLoc))
			{
				assert(0);
				return false;
			}
			if (!writeCoeffs.writeEigenVectorNamed(this->m_strVarCoeffVal, mAllCoeffsVal))
			{
				assert(0);
				return false;
			}
			if (!writeCoeffs.writeEigenVectorNamed(this->m_strVarNNZ, this->m_vNNZ))
			{
				assert(0);
				return false;
			}
		}

		writeCoeffs.close();
		dTimeDelta = omp_get_wtime() - dLastTime;
	}
	else
	{
		// 		double dLastTime = omp_get_wtime();
		// 
		// 		int iNodeID = 0;
		// 		MPI_Comm_rank(MPI_COMM_WORLD, &iNodeID);
		// 
		// 		CCoeffCompression<T, N> coeffComp;
		// 		if (!coeffComp.Init(iNumQuantPoints, strPath))
		// 		{
		// 			assert(0);
		// 			return false;
		// 		}
		// 
		// 		if (iNodeID == 0)
		// 		{
		// 			if (!coeffComp.EncodeMemb(this->m_mMemb))
		// 			{
		// 				assert(0);
		// 				return false;
		// 			}
		// 		}
		// 
		// 		if (this->m_bStoreCoeffs)
		// 		{
		// 			if (iNodeID == 0)
		// 			{
		// 				//Write MATLAB file containing coefficient info.
		// 				VclMatio outFile;
		// 				if (!outFile.openForWriting((strPath / fs::path(this->m_strCoeffsInfoFileName + std::string(".mat"))).string()))
		// 					return false;
		// 				Eigen::VectorXi vDims;
		// 				vDims.resize(this->m_mCoeffsnD->at(0).num_dimensions() + 1);
		// 				for (size_t i = 0; i < vDims.size() - 1; ++i)
		// 					vDims(i) = (this->m_mCoeffsnD->at(0).shape())[i];
		// 				vDims(vDims.size() - 1) = this->m_mCoeffsnD->size();
		// 				if (!outFile.writeEigenVectorNamed(this->m_strVarDim, vDims))
		// 					return false;
		// 				outFile.close();
		// 
		// 				if (!coeffComp.EncodeNNZ(this->m_vNNZ))
		// 				{
		// 					assert(0);
		// 					return false;
		// 				}
		// 				if (!coeffComp.EncodeCoeffs(*this->m_mCoeffsnD))
		// 				{
		// 					assert(0);
		// 					return false;
		// 				}
		// 			}
		// 		}
		// 
		// 		dTimeDelta = omp_get_wtime() - dLastTime;
	}

	return true;
}


template <typename T, size_t N>
bool CDictTestEnsOrthnD<T, N>::Load(const boost::filesystem::path& strPath, const std::string& strFileName, bool bLoadRaw)
{
	typedef typename CData<T, N>::ArraynD::size_type szType;

	if (bLoadRaw)
	{
		VclMatio readCoeffs;

		if (!readCoeffs.openForReading((strPath / fs::path(strFileName)).string()))
		{
			assert(0);
			return false;
		}
		if (!readCoeffs.readEigenMatrix2dNamed(this->m_strVarMemb, this->m_mMemb))
		{
			assert(0);
			return false;
		}
		if (this->m_bStoreCoeffs)
		{
			Eigen::VectorXi vDims;
			if (!readCoeffs.readEigenVectorNamed(this->m_strVarDim, vDims))
			{
				assert(0);
				return false;
			}
			this->m_vShape.resize(vDims.size());
			for (size_t i = 0; i < this->m_vShape.size(); ++i)
				this->m_vShape[i] = vDims(i);
			if (!readCoeffs.readEigenVectorNamed(this->m_strVarNNZ, this->m_vNNZ))
			{
				assert(0);
				return false;
			}
			if (this->m_vNNZ.size() != this->m_mMemb.rows())
			{
				assert(0);
				return false;
			}

			typename CDictTest<T, N>::DataPointCoeffLocs mAllCoeffsLoc;
			typename CDictTest<T, N>::DataPointCoeffVals mAllCoeffsVal;
			if (!readCoeffs.readEigenMatrix2dNamed(this->m_strVarCoeffLoc, mAllCoeffsLoc))
			{
				assert(0);
				return false;
			}
			if (!readCoeffs.readEigenVectorNamed(this->m_strVarCoeffVal, mAllCoeffsVal))
			{
				assert(0);
				return false;
			}
			if (mAllCoeffsLoc.rows() != vDims.size() || mAllCoeffsLoc.cols() != mAllCoeffsVal.size())
			{
				assert(0);
				return false;
			}
			size_t iTotalNNZ = std::accumulate(this->m_vNNZ.data(), this->m_vNNZ.data() + this->m_vNNZ.size(), 0);
			if (iTotalNNZ != mAllCoeffsLoc.cols())
			{
				assert(0);
				return false;
			}

			size_t iCounter = 0;
			this->m_vNZLoc.resize(this->m_vNNZ.size());
			this->m_vNZValue.resize(this->m_vNNZ.size());
			for (size_t i = 0; i < this->m_vNZLoc.size(); ++i)
			{
				this->m_vNZLoc[i].resize(vDims.size(), this->m_vNNZ(i));
				this->m_vNZLoc[i] = mAllCoeffsLoc.block(0, iCounter, mAllCoeffsLoc.rows(), this->m_vNNZ(i));
				this->m_vNZValue[i].resize(this->m_vNNZ(i));
				this->m_vNZValue[i] = mAllCoeffsVal.segment(iCounter, this->m_vNNZ(i));
				iCounter += this->m_vNNZ(i);
			}
		}

		readCoeffs.close();
	}
	else
	{
		// 		CCoeffCompression<T, N> coeffComp;
		// 		if (!coeffComp.Init(-1, strPath))
		// 		{
		// 			assert(0);
		// 			return false;
		// 		}
		// 
		// 		if (!coeffComp.DecodeMemb(this->m_mMemb))
		// 		{
		// 			assert(0);
		// 			return false;
		// 		}
		// 
		// 		if (this->m_bStoreCoeffs)
		// 		{
		// 			//Load coefficient info. from file
		// 			VclMatio inFile;
		// 			if (!inFile.openForReading((strPath / fs::path(this->m_strCoeffsInfoFileName + std::string(".mat"))).string()))
		// 				return false;
		// 			Eigen::VectorXi vDims;
		// 			if (!inFile.readEigenVectorNamed(this->m_strVarDim, vDims))
		// 				return false;
		// 			inFile.close();
		// 
		// 			if (!coeffComp.DecodeNNZ(this->m_vNNZ))
		// 			{
		// 				assert(0);
		// 				return false;
		// 			}
		// 			std::vector<szType> vShape;
		// 			vShape.resize(vDims.size() - 1);
		// 			for (size_t i = 0; i < vShape.size(); ++i)
		// 				vShape[i] = vDims(i);
		// 			SAFE_DELETE(this->m_mCoeffsnD);
		// 			this->m_mCoeffsnD = new std::vector<typename CData<T, N>::ArraynD>(vDims[vDims.size() - 1], typename CData<T, N>::ArraynD(vShape, CCS_TENSOR_STORAGE_ORDER));
		// 			if (!coeffComp.DecodeCoeffs(*this->m_mCoeffsnD))
		// 			{
		// 				assert(0);
		// 				return false;
		// 			}
		// 		}
	}

	return true;
}

template<typename T, size_t N>
inline bool CDictTestEnsOrthnD<T, N>::Concatenate(const std::vector<CDictTest<T, N>*>& vToConcat, const std::vector<size_t>& vShape)
{
	typedef typename CData<T, N>::InternalArray1D InternalVec;
	typedef typename CData<T, N>::InternalArray2D InternalMat;
	typedef typename CData<T, N>::InternalArraynD InternalTensor;

	for (size_t i = 0; i < vToConcat.size(); ++i)
		if (!vToConcat[i])
			return false;

	this->m_vShape = vShape;

	for (size_t i = 0; i < vToConcat.size(); ++i)
	{
		if (!std::equal(vToConcat[i]->GetShape().begin(), vToConcat[i]->GetShape().end(), vToConcat[0]->GetShape().begin()))
		{
			std::cerr << "ERROR: Testing data to concatenate have incompatible dimensionality." << std::endl;
			return false;
		}
	}

	std::vector<size_t> vNumDataPoints(vToConcat.size());
	for (size_t i = 0; i < vToConcat.size(); ++i)
		vNumDataPoints[i] = vToConcat[i]->GetM().rows();
	std::vector<size_t> vStartPos(vNumDataPoints.size()+1);
	vStartPos[0] = 0;
	std::partial_sum(vNumDataPoints.begin(), vNumDataPoints.end(), vStartPos.begin()+1);
	vStartPos.pop_back();
	size_t iTotalNumDataPoints = std::accumulate(vNumDataPoints.begin(), vNumDataPoints.end(), 0);

	//Copy membership indices
	this->m_mMemb.resize(iTotalNumDataPoints, 1);
	for (size_t i = 0; i < vToConcat.size(); ++i)
		this->m_mMemb.block(vStartPos[i], 0, vNumDataPoints[i], vToConcat[i]->GetM().cols()) = vToConcat[i]->GetM();

	//Copy number of nonzeros (NNZ) indices
	this->m_vNNZ.resize(iTotalNumDataPoints);
	for (size_t i = 0; i < vToConcat.size(); ++i)
		this->m_vNNZ.segment(vStartPos[i], vNumDataPoints[i]) = vToConcat[i]->GetNNZ();

	//Copy coefficient values
	this->m_vNZValue.resize(iTotalNumDataPoints);
	for (size_t i = 0; i < vToConcat.size(); ++i)
		for (size_t j = 0; j < vNumDataPoints[i]; ++j)
			this->m_vNZValue[j + i * vNumDataPoints[i]] = (vToConcat[i]->GetNZValue())[j];

	//Copy coefficient locations
	this->m_vNZLoc.resize(iTotalNumDataPoints);
	for (size_t i = 0; i < vToConcat.size(); ++i)
		for (size_t j = 0; j < vNumDataPoints[i]; ++j)
			this->m_vNZLoc[j + i * vNumDataPoints[i]] = (vToConcat[i]->GetNZLoc())[j];

	return true;
}


template <typename T, size_t N>
void CDictTestEnsOrthnD<T, N>::CleanUp()
{
	this->m_mMemb.resize(0, 0);
	this->m_vNNZ.resize(0);
	this->m_vNZValue.clear();
	this->m_vNZLoc.clear();
	this->m_vShape.clear();
}
