


namespace fs = ::boost::filesystem;




template<typename T>
size_t CDataBRDFMerl3D<T>::GetNumDataElems() const
{
	if (this->m_eDataType == CData<T, 3>::DataType::CCS_DATA_1D)
		return this->m_mData1D->size();
	else
		return this->m_mDatanD->size();
}

template <typename T>
void CDataBRDFMerl3D<T>::GetDataElemDim(std::vector<size_t>& vDims) const
{
	if (this->m_eDataType == CData<T, 3>::DataType::CCS_DATA_1D)
	{
		vDims.clear();
		vDims.resize(1);
		vDims[0] = m_iTheta_h * m_iTheta_d * m_iPhi_d;
	}
	else
	{
		vDims.clear();
		vDims.resize(3);
		vDims[0] = m_iTheta_h;
		vDims[1] = m_iTheta_d;
		vDims[2] = m_iPhi_d;
	}
}


template <typename T>
bool CDataBRDFMerl3D<T>::Init()
{
	CleanUp();

	return true;
}

template <typename T>
CDataBRDFMerl3D<T>::CDataBRDFMerl3D(uint32_t iTheta_h, uint32_t iTheta_d, uint32_t iPhi_d, const std::string& strBRDFVarName) : CData<T, 3>()
{
	this->m_eDataType = CData<T, 3>::DataType::CCS_DATA_2D;
	m_iTheta_h = iTheta_h;
	m_iTheta_d = iTheta_d;
	m_iPhi_d = iPhi_d;
	m_strBRDFVarName = strBRDFVarName;
}


template <typename T>
CDataBRDFMerl3D<T>::~CDataBRDFMerl3D()
{
	CleanUp();
}


template <typename T>
bool CDataBRDFMerl3D<T>::LoadFromDisk(const std::string& strDataFolder)
{
	this->m_strDataFolder = strDataFolder;

	fs::path rootPath(this->m_strDataFolder);
	if (!fs::exists(rootPath) || !fs::is_directory(rootPath))
	{
		std::cerr << "ERROR: The input folder does not exist." << std::endl;
		return false;
	}

	//Find BRDFs in the root folder
	m_vFileNames.clear();
	fs::directory_iterator it_end;
	for (fs::directory_iterator it(rootPath); it != it_end; ++it)
		if (fs::is_regular_file(it->status()))
			m_vFileNames.push_back(it->path().filename());

	//Sort BRDF file names
	std::sort(m_vFileNames.begin(), m_vFileNames.end(), [](const fs::path& p1, const fs::path& p2) {return p1.filename().string() < p2.filename().string(); });

	std::vector<size_t> vDims;
	GetDataElemDim(vDims);
	vDims.push_back(3);
	m_vBRDF.resize(m_vFileNames.size(), typename CData<T, 4>::ArraynD(vDims, CCS_TENSOR_STORAGE_ORDER));

	for (size_t i = 0; i < m_vFileNames.size(); ++i)
	{
		VclMatio matReader;

		if (!matReader.openForReading((rootPath / m_vFileNames[i]).string()))
		{
			std::cerr << "ERROR: Cannot read BRDF file at " << (rootPath / m_vFileNames[i]).string() << std::endl;
			return false;
		}
		if(!matReader.readEigenMatrixndNamed(m_strBRDFVarName, m_vBRDF[i]))
		{
			std::cerr << "ERROR: Cannot read BRDF variable at " << (rootPath / m_vFileNames[i]).string() << std::endl;
			return false;
		}

		matReader.close();
	}

	return true;
}


template <typename T>
bool CDataBRDFMerl3D<T>::PrepareData(const std::vector<uint32_t>& vPatchSize, PatchType ePatchType, const std::vector<uint32_t>& vSlidingDis, bool bFreeLoadedData)
{
	this->m_vPatchSize = vPatchSize;
	this->m_ePatchType = ePatchType;
	this->m_vSlidingDis = vSlidingDis;

	if (m_vBRDF.size() == 0)
	{
		std::cerr << "ERROR: Data not loaded yet. Cannot prepare data." << std::endl;
		return false;
	}

	m_vNumPatches.resize(m_vFileNames.size());

	//Compute the total number of patches
	for (size_t i = 0; i < m_vNumPatches.size(); ++i)
		m_vNumPatches[i] = 3;
	size_t iTotalNumPatches = m_vNumPatches.size() * 3;

	//Allocate space for input data
	std::vector<size_t> vDims;
	GetDataElemDim(vDims);
	size_t iLength = std::accumulate(vDims.begin(), vDims.end(), 1, std::multiplies<size_t>());
	SAFE_DELETE(this->m_mDatanD);
	this->m_mDatanD = new std::vector<typename CData<T, 3>::ArraynD>(iTotalNumPatches, typename CData<T, 3>::ArraynD(vDims, CCS_TENSOR_STORAGE_ORDER));

	for (size_t i = 0; i < m_vNumPatches.size(); ++i)
	{
		for (size_t j = 0; j < iLength; ++j)
		{
			this->m_mDatanD->at(i * 3).data()[j] = m_vBRDF[i].data()[j];
			this->m_mDatanD->at(i * 3 + 1).data()[j] = m_vBRDF[i].data()[iLength + j];
			this->m_mDatanD->at(i * 3 + 2).data()[j] = m_vBRDF[i].data()[2 * iLength + j];
		}
	}

	if (bFreeLoadedData)
		m_vBRDF.clear();

	return true;
}



template<typename T>
void CDataBRDFMerl3D<T>::CreateMask(MaskRandMode eRandMode, float fNonZeroRatio)
{

}


template <typename T>
void CDataBRDFMerl3D<T>::Resize(size_t iNumPoints, const std::vector<size_t> vDims)
{
	if (this->m_eDataType == CData<T, 3>::DataType::CCS_DATA_1D)
	{
		size_t iLength = std::accumulate(vDims.begin(), vDims.end(), 1, std::multiplies<size_t>());
		SAFE_DELETE(this->m_mData1D);
		this->m_mData1D = new std::vector<typename CData<T, 3>::Array1D>(iNumPoints);
		for (size_t i = 0; i < this->m_mData1D->size(); ++i)
			this->m_mData1D->at(i).setZero(iLength);
	}
	else
	{
		SAFE_DELETE(this->m_mDatanD);
		this->m_mDatanD = new std::vector<typename CData<T, 3>::ArraynD>(iNumPoints, typename CData<T, 3>::ArraynD(vDims, CCS_TENSOR_STORAGE_ORDER));
		for (size_t i = 0; i < this->m_mDatanD->size(); ++i)
			TensorSetZero(this->m_mDatanD->at(i));
	}
}


template <typename T>
template <typename U>
void CDataBRDFMerl3D<T>::ResizeLike(const CData<U, 3>* pData)
{
	CleanUp();
	CDataBRDFMerl3D<U>* pDataBRDF = dynamic_cast<CDataBRDFMerl3D<U>*>(const_cast<CData<U, 3>*>(pData));
	if (!pDataBRDF)
		return;

	CopyPropsFrom(pDataBRDF);
	std::vector<size_t> vDims;
	pDataBRDF->GetDataElemDim(vDims);
	if (pDataBRDF->GetDataType() == CData<T, 3>::DataType::CCS_DATA_1D)
	{
		size_t iLength = std::accumulate(vDims.begin(), vDims.end(), 1, std::multiplies<size_t>());
		SAFE_DELETE(this->m_mData1D);
		this->m_mData1D = new std::vector<typename CData<T, 3>::Array1D>(pDataBRDF->GetNumDataElems());
		for (size_t i = 0; i < this->m_mData1D->size(); ++i)
			this->m_mData1D->at(i).setZero(iLength);
	}
	else
	{
		SAFE_DELETE(this->m_mDatanD);
		this->m_mDatanD = new std::vector<typename CData<T, 3>::ArraynD>(pDataBRDF->GetNumDataElems(), typename CData<T, 3>::ArraynD(vDims, CCS_TENSOR_STORAGE_ORDER));
		for (size_t i = 0; i < this->m_mDatanD->size(); ++i)
			TensorSetZero(this->m_mDatanD->at(i));
	}
}


template <typename T>
template <typename U>
void CDataBRDFMerl3D<T>::CopyPropsFrom(const CData<U, 3>* pData)
{
	const CDataBRDFMerl3D<U>* pDataBRDF = dynamic_cast<const CDataBRDFMerl3D<U>*>(pData);
	if (!pDataBRDF)
		return;
	m_vFileNames = pDataBRDF->m_vFileNames;
	m_vNumPatches = pDataBRDF->m_vNumPatches;
	m_iTheta_h = pDataBRDF->m_iTheta_h;
	m_iTheta_d = pDataBRDF->m_iTheta_d;
	m_iPhi_d = pDataBRDF->m_iPhi_d;
	m_strBRDFVarName = pDataBRDF->m_strBRDFVarName;
	this->m_eDataType = pDataBRDF->m_eDataType;
	this->m_vPatchSize = pDataBRDF->m_vPatchSize;
	this->m_ePatchType = pDataBRDF->m_ePatchType;
	this->m_vSlidingDis = pDataBRDF->m_vSlidingDis;
	this->m_iQuantization = pDataBRDF->m_iQuantization;
	this->m_strDataFolder = pDataBRDF->m_strDataFolder;
	this->m_strOutputFolder = pDataBRDF->m_strOutputFolder;
}


template <typename T>
void CDataBRDFMerl3D<T>::ConvertTo1D(bool bDeleteOld)
{
	size_t iNumPoints = GetNumDataElems();
	std::vector<size_t> vDims;
	GetDataElemDim(vDims);
	size_t iLength = std::accumulate(vDims.begin(), vDims.end(), 1, std::multiplies<size_t>());
	this->m_eDataType = CData<T, 3>::DataType::CCS_DATA_1D;
	SAFE_DELETE(this->m_mData1D);
	this->m_mData1D = new std::vector<typename CData<T, 3>::Array1D>(iNumPoints);
	for (size_t i = 0; i < iNumPoints; ++i)
	{
		this->m_mData1D->at(i).resize(iLength);
		for (size_t j = 0; j < iLength; ++j)
			(this->m_mData1D->at(i))(j) = (this->m_mDatanD->at(i).data())[j];
	}
	if (bDeleteOld)
		SAFE_DELETE(this->m_mDatanD);

	if (this->m_mMasknD)
	{
		SAFE_DELETE(this->m_mMask1D);
		this->m_mMask1D = new std::vector<typename CData<T, 3>::Msk1D>(iNumPoints);
		for (size_t i = 0; i < iNumPoints; ++i)
		{
			this->m_mMask1D->at(i).resize(iLength);
			for (size_t j = 0; j < iLength; ++j)
				(this->m_mMask1D->at(i))(j) = (this->m_mMasknD->at(i).data())[j];
		}
		if (bDeleteOld)
			SAFE_DELETE(this->m_mMasknD);
	}
}


template<typename T>
void CDataBRDFMerl3D<T>::ConvertBack(bool bDeleteOld)
{
	size_t iNumPoints = GetNumDataElems();
	this->m_eDataType = CData<T, 3>::DataType::CCS_DATA_nD;
	std::vector<size_t> vDims;
	GetDataElemDim(vDims);
	size_t iLength = std::accumulate(vDims.begin(), vDims.end(), 1, std::multiplies<size_t>());
	SAFE_DELETE(this->m_mDatanD);
	this->m_mDatanD = new std::vector<typename CData<T, 3>::ArraynD>(iNumPoints, typename CData<T, 3>::ArraynD(vDims, CCS_TENSOR_STORAGE_ORDER));
	for (size_t i = 0; i < iNumPoints; ++i)
		for (size_t j = 0; j < iLength; ++j)
			(this->m_mDatanD->at(i).data())[j] = (this->m_mData1D->at(i))(j);
	if (bDeleteOld)
		SAFE_DELETE(this->m_mData1D);

	if (this->m_mMask1D)
	{
		SAFE_DELETE(this->m_mMasknD);
		this->m_mMasknD = new std::vector<typename CData<T, 3>::MsknD>(iNumPoints, typename CData<T, 3>::MsknD(vDims, CCS_TENSOR_STORAGE_ORDER));
		for (size_t i = 0; i < iNumPoints; ++i)
			for (size_t j = 0; j < iLength; ++j)
				(this->m_mMasknD->at(i).data())[j] = (this->m_mMask1D->at(i))(j);
		if (bDeleteOld)
			SAFE_DELETE(this->m_mMask1D);
	}
}


template <typename T>
void CDataBRDFMerl3D<T>::Clamp(T minVal, T maxVal)
{
	if (this->m_eDataType == CData<T, 3>::DataType::CCS_DATA_1D)
	{
#pragma omp parallel for schedule(dynamic,1)
		for (int i = 0; i < GetNumDataElems(); ++i)
			Clamp1D(this->m_mData1D->at(i), minVal, maxVal);
	}
	else
	{
#pragma omp parallel for schedule(dynamic,1)
		for (int i = 0; i < GetNumDataElems(); ++i)
			ClampnD(this->m_mDatanD->at(i), minVal, maxVal);
	}
}


template <typename T>
void CDataBRDFMerl3D<T>::RemoveZeros(CCS_INTERNAL_TYPE thresh)
{
	if (thresh < 0.0)
		return;

	std::vector<size_t> vDims;
	GetDataElemDim(vDims);
	size_t iLength = std::accumulate(vDims.begin(), vDims.end(), 1, std::multiplies<size_t>());
	std::vector<size_t> vIdx;
	vIdx.reserve(GetNumDataElems());
	for (size_t i = 0; i < GetNumDataElems(); ++i)
	{
		typename CData<T, 3>::InternalArraynD dataPoint(vDims, CCS_TENSOR_STORAGE_ORDER);
		CastDataPointTTypeToInternal<T, 3>(this->m_mDatanD->at(i), dataPoint);
		if (TensorNorm2(dataPoint) > thresh)
			vIdx.push_back(i);
	}

	std::vector<typename CData<T, 3>::ArraynD>* mDatanD = new std::vector<typename CData<T, 3>::ArraynD>(vIdx.size(), typename CData<T, 3>::ArraynD(vDims, CCS_TENSOR_STORAGE_ORDER));
	for (size_t i = 0; i < vIdx.size(); ++i)
		for (size_t j = 0; j < iLength; ++j)
			(mDatanD->at(i).data())[j] = (this->m_mDatanD->at(vIdx[i]).data())[j];
	SAFE_DELETE(this->m_mDatanD);
	this->m_mDatanD = mDatanD;
}


template <typename T>
bool CDataBRDFMerl3D<T>::AssembleData()
{
	if (m_vNumPatches.empty() || this->m_mDatanD->empty())
	{
		std::cerr << "ERROR: Data has not been read to disk." << std::endl;
		return false;
	}

	std::vector<size_t> vBRDFdims = { m_iTheta_h , m_iTheta_d , m_iPhi_d , 3 };
	size_t iLength = m_iTheta_h * m_iTheta_d * m_iPhi_d;
	m_vBRDF.resize(m_vNumPatches.size(), typename  CData<T, 4>::ArraynD(vBRDFdims, CCS_TENSOR_STORAGE_ORDER));

	for (size_t i = 0; i < m_vNumPatches.size(); ++i)
	{
		for (size_t j = 0; j < iLength; ++j)
		{
			m_vBRDF[i].data()[j] = this->m_mDatanD->at(i * 3).data()[j];
			m_vBRDF[i].data()[iLength + j] = this->m_mDatanD->at(i * 3 + 1).data()[j];
			m_vBRDF[i].data()[2 * iLength + j] = this->m_mDatanD->at(i * 3 + 2).data()[j];
		}
	}

	return true;
}


template<typename T>
void CDataBRDFMerl3D<T>::CalcQuality(const CData<T, 3>* pOther, ReconQualityMetric eMetric, std::vector<CCS_INTERNAL_TYPE>& vQlty)
{
	const CDataBRDFMerl3D<T>* pDataRGB = dynamic_cast<const CDataBRDFMerl3D<T>*>(pOther);
	if (!pDataRGB)
		return;

	vQlty.resize(m_vBRDF.size(), 0.0f);
	for (size_t i = 0; i < m_vBRDF.size(); ++i)
	{
		vQlty[i] += (CData<T, 3>::Array1D::Map(m_vBRDF[i].data(), m_vBRDF[i].num_elements()).template cast<CCS_INTERNAL_TYPE>() - 
			CData<T, 3>::Array1D::Map(pDataRGB->m_vBRDF[i].data(), pDataRGB->m_vBRDF[i].num_elements()).template cast<CCS_INTERNAL_TYPE>()).squaredNorm();
		vQlty[i] /= m_vBRDF[i].num_elements();
		switch (eMetric)
		{
		case CCS_RECON_QLTY_MSE:
			break;
		case CCS_RECON_QLTY_PSNR:
			vQlty[i] = PSNR<T>(vQlty[i]);
			break;
		}
	}
}



template <typename T>
bool CDataBRDFMerl3D<T>::WriteAssembled(const std::string& strOutputFolder)
{
	if (m_vBRDF.size() != m_vFileNames.size())
	{
		std::cerr << "ERROR: Number of BRDF filenames not equal to the number of BRDFs." << std::endl;
		return false;
	}

	//Create the output folder (if there is none)
	this->m_strOutputFolder = strOutputFolder;
	fs::path outputPath(this->m_strOutputFolder);
	if (!fs::exists(outputPath) || !fs::is_directory(outputPath))
	{
		if (!fs::create_directory(outputPath))
		{
			std::cerr << "ERROR: Unable to create output directory to store results." << std::endl;
			return false;
		}
	}

	//Iterate over each BRDF and write them to disk as MAT files
	for (int i = 0; i < m_vFileNames.size(); ++i)
	{
		VclMatio matWriter;
		if (!matWriter.openForWriting((outputPath / m_vFileNames[i]).string()))
		{
			std::cerr << "ERROR: Cannot write BRDF file at " << (outputPath / m_vFileNames[i]).string() << std::endl;
			return false;
		}
		if (!matWriter.writeEigenMatrixndNamed(m_strBRDFVarName, m_vBRDF[i]))
		{
			std::cerr << "ERROR: Cannot write BRDF variable at " << (outputPath / m_vFileNames[i]).string() << std::endl;
			return false;
		}
		matWriter.close();
	}

	return true;
}


template <typename T>
void CDataBRDFMerl3D<T>::CleanUp()
{
	SAFE_DELETE(this->m_mData1D);
	SAFE_DELETE(this->m_mData2D);
	SAFE_DELETE(this->m_mDatanD);
	SAFE_DELETE(this->m_mMask1D);
	SAFE_DELETE(this->m_mMask2D);
	SAFE_DELETE(this->m_mMasknD);
	this->m_strDataFolder.clear();
	this->m_strOutputFolder.clear();
	this->m_ePatchType = CCS_PATCHTYPE_NOV;
	this->m_vSlidingDis.clear();
	this->m_iQuantization = 0;
	this->m_vPatchSize.clear();

	m_vBRDF.clear();
	m_vFileNames.clear();
	m_vNumPatches.clear();
}
