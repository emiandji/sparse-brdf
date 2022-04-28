#pragma once


#include "../cssUtil/defs.h"



enum PatchType { CCS_PATCHTYPE_NOV = 0, CCS_PATCHTYPE_OV = 1 };
enum ColorMode { CCS_COLOR_RGB = 0, CCS_COLOR_BW = 1, CCS_COLOR_YCC = 2 };
enum MaskRandMode { CCS_MASK_RND_UNIFORM = 0, CCS_MASK_RND_HALTON = 1 };
enum ReconQualityMetric { CCS_RECON_QLTY_MSE = 0, CCS_RECON_QLTY_PSNR = 1 };


template <typename T, size_t N>
class CData
{

public: 

	enum DataType {CCS_DATA_1D = 0, CCS_DATA_2D = 1, CCS_DATA_nD = 2};

	typedef Eigen::Matrix<T, Eigen::Dynamic, 1>					Array1D;
	typedef Eigen::Matrix<T, 1, Eigen::Dynamic>					Array1DRow;
	typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>	Array2D;
	typedef boost::multi_array<T, N>							ArraynD;

	typedef Eigen::Matrix<uint8_t, Eigen::Dynamic, 1>				Msk1D;
	typedef Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>	Msk2D;
	typedef boost::multi_array<uint8_t, N>							MsknD;

	typedef Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, 1>					InternalArray1D;
	typedef Eigen::Matrix<CCS_INTERNAL_TYPE, 1, Eigen::Dynamic>					InternalArray1DRow;
	typedef Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, Eigen::Dynamic>	InternalArray2D;
	typedef boost::multi_array<CCS_INTERNAL_TYPE, N>							InternalArraynD;

	typedef Eigen::Matrix<CCS_GPU_INTERNAL_TYPE, Eigen::Dynamic, 1>					GPUArray1D;
	typedef Eigen::Matrix<CCS_GPU_INTERNAL_TYPE, 1, Eigen::Dynamic>					GPUArray1DRow;
	typedef Eigen::Matrix<CCS_GPU_INTERNAL_TYPE, Eigen::Dynamic, Eigen::Dynamic>	GPUArray2D;
	typedef boost::multi_array<CCS_GPU_INTERNAL_TYPE, N>							GPUArraynD;

public:
    
	CData()
	{
		m_eDataType = CCS_DATA_1D;
		m_mData1D = NULL;
		m_mData2D = NULL;
		m_mDatanD = NULL;
		m_mMask1D = NULL;
		m_mMask2D = NULL;
		m_mMasknD = NULL;
		m_ePatchType = CCS_PATCHTYPE_NOV;
		m_iQuantization = -1;
	}
	virtual ~CData() 
	{
		SAFE_DELETE(m_mData1D);
		SAFE_DELETE(m_mData2D);
		SAFE_DELETE(m_mDatanD);
	}

	virtual const std::vector<Array1D>* GetData1D() const { return this->m_mData1D; }
	virtual const std::vector<Array2D>* GetData2D() const { return this->m_mData2D; }
	virtual const std::vector<ArraynD>* GetDatanD() const { return this->m_mDatanD; }


	virtual const std::vector<Msk1D>* GetMask1D() const { return this->m_mMask1D; }
	virtual const std::vector<Msk2D>* GetMask2D() const { return this->m_mMask2D; }
	virtual const std::vector<MsknD>* GetMasknD() const { return this->m_mMasknD; }

	virtual size_t GetNumDataElems() const = 0;
	virtual void GetDataElemDim(std::vector<size_t>& vDims) const = 0;		//call after LoadFromDisk() and PrepareData()

	virtual DataType GetDataType()					const	{ return m_eDataType; }
	virtual const std::string& GetDataFolder()		const	{ return m_strDataFolder; }
	virtual const std::string& GetOutputFolder()	const	{ return m_strOutputFolder; }
	virtual int64_t GetQuantization()				const	{ return m_iQuantization; }

	virtual void SetDataPoint(size_t iIdx, const InternalArray1D& dataPoint);
	virtual void SetDataPoint(size_t iIdx, const InternalArray2D& dataPoint);
	virtual void SetDataPoint(size_t iIdx, const InternalArraynD& dataPoint);

	virtual void GetScaleNormalizedDataPoint(size_t iIdx, InternalArray1D& dataPoint);
	virtual void GetScaleNormalizedDataPoint(size_t iIdx, InternalArray1DRow& dataPoint);
	virtual void GetScaleNormalizedDataPoint(size_t iIdx, InternalArray2D& dataPoint);
	virtual void GetScaleNormalizedDataPoint(size_t iIdx, InternalArraynD& dataPoint);

	virtual void GetNormalizedDataPoint(size_t iIdx, InternalArray1D& dataPoint);
	virtual void GetNormalizedDataPoint(size_t iIdx, InternalArray1DRow& dataPoint);
	virtual void GetNormalizedDataPoint(size_t iIdx, InternalArray2D& dataPoint);
	virtual void GetNormalizedDataPoint(size_t iIdx, InternalArraynD& dataPoint);

	virtual bool Init() = 0;
	virtual bool LoadFromDisk(const std::string& strDataFolder) = 0;
	virtual bool PrepareData(const std::vector<uint32_t>& vPatchSize, PatchType ePatchType, const std::vector<uint32_t>& vSlidingDis, bool bFreeLoadedData) = 0;
	virtual void CreateMask(MaskRandMode eRandMode, float fNonZeroRatio) = 0;		//This function should be called after PrepareData()
	virtual void Resize(size_t iNumPoints, const std::vector<size_t> vDims) = 0;
// 	template <typename U> virtual void ResizeLike(const CData<U, N>* pData) = 0;
// 	template <typename U> virtual void CopyPropsFrom(const CData<U, N>* pData) = 0;
	virtual void ConvertTo1D(bool bDeleteOld) = 0;
	virtual void ConvertBack(bool bDeleteOld) = 0;
	virtual void Clamp(T minVal, T maxVal) = 0;
	virtual void RemoveZeros(CCS_INTERNAL_TYPE thresh) = 0;			//Remove data points with zero norm. AssembleData() cannot be called after this function.
	virtual bool AssembleData() = 0;
	virtual void CalcQuality(const CData<T, N>* pOther, ReconQualityMetric eMetric, std::vector<CCS_INTERNAL_TYPE>& vQlty) = 0;
	virtual bool WriteAssembled(const std::string& strOutputFolder) = 0;
	virtual void CleanUp() = 0;
    
protected:

	DataType m_eDataType;
    
	std::string m_strDataFolder;
	std::string m_strOutputFolder;

	std::vector<Array1D>* m_mData1D;
	std::vector<Array2D>* m_mData2D;
	std::vector<ArraynD>* m_mDatanD;

	std::vector<Msk1D>* m_mMask1D;
	std::vector<Msk2D>* m_mMask2D;
	std::vector<MsknD>* m_mMasknD;

	std::vector<uint32_t> m_vPatchSize;
	PatchType m_ePatchType;
	std::vector<uint32_t> m_vSlidingDis;
	int64_t m_iQuantization;					//Set to -1 to use floating point values. It is assumed that all the data share one quantization factor
};




//It is assumed that the input is in [-1,1] or [0,1]
template <typename T>
T CastInternalToType(CCS_INTERNAL_TYPE input)
{
	throw std::runtime_error("Type is not supported to cast from InternalType.");
}

template <>
inline uint8_t CastInternalToType<uint8_t>(CCS_INTERNAL_TYPE input)
{
	if (input > 1.0)
		input = 1.0;
	if (input < 0.0)
		input = 0.0;
	return uint8_t(std::round(input * std::numeric_limits<uint8_t>::max()));
}

template <>
inline int8_t CastInternalToType<int8_t>(CCS_INTERNAL_TYPE input)
{
	if (input > 1.0)
		input = 1.0;
	if (input < -1.0)
		input = -1.0;
	return int8_t(std::round(input * std::numeric_limits<int8_t>::max()));
}

template <>
inline uint16_t CastInternalToType<uint16_t>(CCS_INTERNAL_TYPE input)
{
	if (input > 1.0)
		input = 1.0;
	if (input < 0.0)
		input = 0.0;
	return uint16_t(std::round(input * std::numeric_limits<uint16_t>::max()));
}

template <>
inline int16_t CastInternalToType<int16_t>(CCS_INTERNAL_TYPE input)
{
	if (input > 1.0)
		input = 1.0;
	if (input < -1.0)
		input = -1.0;
	return int16_t(std::round(input * std::numeric_limits<int16_t>::max()));
}

template <>
inline uint32_t CastInternalToType<uint32_t>(CCS_INTERNAL_TYPE input)
{
	if (input > 1.0)
		input = 1.0;
	if (input < 0.0)
		input = 0.0;
	return uint32_t(std::round(input * std::numeric_limits<uint32_t>::max()));
}

template <>
inline int32_t CastInternalToType<int32_t>(CCS_INTERNAL_TYPE input)
{
	if (input > 1.0)
		input = 1.0;
	if (input < -1.0)
		input = -1.0;
	return int32_t(std::round(input * std::numeric_limits<int32_t>::max()));
}

template <>
inline uint64_t CastInternalToType<uint64_t>(CCS_INTERNAL_TYPE input)
{
	if (input > 1.0)
		input = 1.0;
	if (input < 0.0)
		input = 0.0;
	return uint64_t(std::round(input * std::numeric_limits<uint64_t>::max()));
}

template <>
inline int64_t CastInternalToType<int64_t>(CCS_INTERNAL_TYPE input)
{
	if (input > 1.0)
		input = 1.0;
	if (input < -1.0)
		input = -1.0;
	return int64_t(std::round(input * std::numeric_limits<int64_t>::max()));
}

template <>
inline double CastInternalToType<double>(CCS_INTERNAL_TYPE input)
{
	return input;
}

template <>
inline float CastInternalToType<float>(CCS_INTERNAL_TYPE input)
{
	return input;
}




//It is assumed that the input is in [-1,1] or [0,1]
template <typename T>
CCS_INTERNAL_TYPE CastTypeToInternal(T input)
{
	throw std::runtime_error("Type is not supported to cast to InternalType.");
}

template <>
inline CCS_INTERNAL_TYPE CastTypeToInternal<int8_t>(int8_t input)
{
	return (CCS_INTERNAL_TYPE)(input) / std::numeric_limits<int8_t>::max();
}

template <>
inline CCS_INTERNAL_TYPE CastTypeToInternal<uint8_t>(uint8_t input)
{
	return (CCS_INTERNAL_TYPE)(input) / std::numeric_limits<uint8_t>::max();
}

template <>
inline CCS_INTERNAL_TYPE CastTypeToInternal<int16_t>(int16_t input)
{
	return (CCS_INTERNAL_TYPE)(input) / std::numeric_limits<int16_t>::max();
}

template <>
inline CCS_INTERNAL_TYPE CastTypeToInternal<uint16_t>(uint16_t input)
{
	return (CCS_INTERNAL_TYPE)(input) / std::numeric_limits<uint16_t>::max();
}

template <>
inline CCS_INTERNAL_TYPE CastTypeToInternal<int32_t>(int32_t input)
{
	return (CCS_INTERNAL_TYPE)(input) / std::numeric_limits<int32_t>::max();
}

template <>
inline CCS_INTERNAL_TYPE CastTypeToInternal<uint32_t>(uint32_t input)
{
	return (CCS_INTERNAL_TYPE)(input) / std::numeric_limits<uint32_t>::max();
}

template <>
inline CCS_INTERNAL_TYPE CastTypeToInternal<int64_t>(int64_t input)
{
	return (CCS_INTERNAL_TYPE)(input) / std::numeric_limits<int64_t>::max();
}

template <>
inline CCS_INTERNAL_TYPE CastTypeToInternal<uint64_t>(uint64_t input)
{
	return (CCS_INTERNAL_TYPE)(input) / std::numeric_limits<uint64_t>::max();
}

template <>
inline CCS_INTERNAL_TYPE CastTypeToInternal<double>(double input)
{
	return input;
}

template <>
inline CCS_INTERNAL_TYPE CastTypeToInternal<float>(float input)
{
	return input;
}


//T -> CCS_INTERNAL_TYPE
template <typename Derived1, typename Derived2>
inline void CastDataPointTTypeToInternal(const Eigen::MatrixBase<Derived1>& src, Eigen::MatrixBase<Derived2> const & dst)
{
	Eigen::MatrixBase<Derived2>& _dst = const_cast<Eigen::MatrixBase<Derived2>&>(dst);

	for (size_t i = 0; i < src.rows(); ++i)
		for (size_t j = 0; j < src.cols(); ++j)
			_dst(i, j) = CastTypeToInternal<typename Derived1::Scalar>(src(i, j));
}

////T -> CCS_INTERNAL_TYPE
//template <typename T>
//inline void CastDataPointTTypeToInternal(const typename CData<T, 0>::Array1D& src, typename CData<T, 0>::InternalArray1D& dst)
//{
//	for (size_t i = 0; i < src.size(); ++i)
//		dst(i) = CastTypeToInternal<T>(src(i));
//}
//
////T -> CCS_INTERNAL_TYPE
//template <typename T>
//inline void CastDataPointTTypeToInternal(const typename CData<T, 0>::Array2D& src, typename CData<T, 0>::InternalArray2D& dst)
//{
//	for (size_t i = 0; i < src.rows(); ++i)
//		for (size_t i = 0; i < src.cols(); ++i)
//			dst(i, j) = CastTypeToInternal<T>(src(i, j));
//}

//T -> CCS_INTERNAL_TYPE
template <typename T, size_t N>
inline void CastDataPointTTypeToInternal(const typename CData<T, N>::ArraynD& src, typename CData<T, N>::InternalArraynD& dst)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	szType iLength = std::accumulate(src.shape(), src.shape() + N, szType(1), std::multiplies<szType>());
	for (size_t i = 0; i < iLength; ++i)
		dst.data()[i] = CastTypeToInternal<T>(src.data()[i]);
}

//CCS_INTERNAL_TYPE -> T
template <typename Derived1, typename Derived2>
inline void CastDataPointInternalToTType(const Eigen::MatrixBase<Derived1>& src, Eigen::MatrixBase<Derived2> const & dst)
{
	Eigen::MatrixBase<Derived2>& _dst = const_cast<Eigen::MatrixBase<Derived2>&>(dst);

	for (size_t i = 0; i < src.rows(); ++i)
		for (size_t j = 0; j < src.cols(); ++j)
			_dst(i, j) = CastInternalToType<typename Derived2::Scalar>(src(i, j));
}

////CCS_INTERNAL_TYPE -> T
//template <typename T>
//inline void CastDataPointInternalToTType(const typename CData<T, 0>::InternalArray1D& src, typename CData<T, 0>::Array1D& dst)
//{
//	for (size_t i = 0; i < src.size(); ++i)
//		dst(i) = CastInternalToType<T>(src(i));
//}
//
////CCS_INTERNAL_TYPE -> T
//template <typename T>
//inline void CastDataPointInternalToTType(const typename CData<T, 0>::InternalArray2D& src, typename CData<T, 0>::Array2D& dst)
//{
//	for (size_t i = 0; i < src.rows(); ++i)
//		for (size_t i = 0; i < src.cols(); ++i)
//			dst(i, j) = CastInternalToType<T>(src(i, j));
//}

//CCS_INTERNAL_TYPE -> T
template <typename T, size_t N>
inline void CastDataPointInternalToTType(const typename CData<T, N>::InternalArraynD& src, typename CData<T, N>::ArraynD& dst)
{
	typedef typename boost::multi_array<CCS_INTERNAL_TYPE, N>::size_type szType;
	szType iLength = std::accumulate(src.shape(), src.shape() + N, szType(1), std::multiplies<szType>());
	for (size_t i = 0; i < iLength; ++i)
		dst.data()[i] = CastInternalToType<T>(src.data()[i]);
}






template <typename T, size_t N>
void CData<T, N>::SetDataPoint(size_t iIdx, const InternalArray1D& dataPoint)
{
	CastDataPointInternalToTType(dataPoint, m_mData1D->at(iIdx));
// 	for (size_t i = 0; i < (*m_mData1D)[iIdx].size(); ++i)
// 		((*m_mData1D)[iIdx])(i) = CastInternalToType<T>(dataPoint(i));
}

template <typename T, size_t N>
void CData<T, N>::SetDataPoint(size_t iIdx, const InternalArray2D& dataPoint)
{
	CastDataPointInternalToTType(dataPoint, m_mData2D->at(iIdx));
// 	for (size_t i = 0; i < (*m_mData2D)[iIdx].rows(); ++i)
// 		for (size_t j = 0; j < (*m_mData2D)[iIdx].cols(); ++j)
// 			((*m_mData2D)[iIdx])(i, j) = CastInternalToType<T>(dataPoint(i, j));
}

template <typename T, size_t N>
void CData<T, N>::SetDataPoint(size_t iIdx, const InternalArraynD& dataPoint)
{
	CastDataPointInternalToTType<T, N>(dataPoint, m_mDatanD->at(iIdx));
// 	typedef typename boost::multi_array<T, N>::size_type szType;
// 	typedef typename boost::multi_array<T, N>::index idxType;
// 	szType iLength = std::accumulate(dataPoint.shape(), dataPoint.shape() + N, szType(1), std::multiplies<szType>());
// 	for (idxType i = 0; i < iLength; ++i)
// 		((*m_mDatanD)[iIdx]).data()[i] = CastInternalToType<T>(dataPoint.data()[i]);
}



template<typename T, size_t N>
void CData<T, N>::GetNormalizedDataPoint(size_t iIdx, InternalArray1D& dataPoint)
{
	CastDataPointTTypeToInternal(m_mData1D->at(iIdx), dataPoint);
	dataPoint /= dataPoint.squaredNorm();
}

template<typename T, size_t N>
void CData<T, N>::GetNormalizedDataPoint(size_t iIdx, InternalArray1DRow& dataPoint)
{
	CastDataPointTTypeToInternal(Array1DRow::Map(m_mData1D->at(iIdx).data(), m_mData1D->at(iIdx).size()), dataPoint);
	dataPoint /= dataPoint.squaredNorm();
}

template<typename T, size_t N>
void CData<T, N>::GetNormalizedDataPoint(size_t iIdx, InternalArray2D& dataPoint)
{
	CastDataPointTTypeToInternal(m_mData2D->at(iIdx), dataPoint);
	dataPoint /= dataPoint.squaredNorm();
}

template<typename T, size_t N>
void CData<T, N>::GetNormalizedDataPoint(size_t iIdx, InternalArraynD& dataPoint)
{
	CastDataPointTTypeToInternal<T, N>(m_mDatanD->at(iIdx), dataPoint);
	size_t iNumElems = dataPoint.num_elements();
	CCS_INTERNAL_TYPE maxVal = InternalArray1D::Map(dataPoint.data(), iNumElems).squaredNorm();
	for (size_t i = 0; i < iNumElems; ++i)
		dataPoint.data()[i] /= maxVal;
}



template<typename T, size_t N>
void CData<T, N>::GetScaleNormalizedDataPoint(size_t iIdx, InternalArray1D& dataPoint)
{
	CastDataPointTTypeToInternal(m_mData1D->at(iIdx), dataPoint);
	dataPoint /= dataPoint.cwiseAbs().maxCoeff();
}

template<typename T, size_t N>
void CData<T, N>::GetScaleNormalizedDataPoint(size_t iIdx, InternalArray1DRow& dataPoint)
{
	CastDataPointTTypeToInternal(Array1DRow::Map(m_mData1D->at(iIdx).data(), m_mData1D->at(iIdx).size()), dataPoint);
	dataPoint /= dataPoint.cwiseAbs().maxCoeff();
}

template<typename T, size_t N>
void CData<T, N>::GetScaleNormalizedDataPoint(size_t iIdx, InternalArray2D& dataPoint)
{
	CastDataPointTTypeToInternal(m_mData2D->at(iIdx), dataPoint);
	dataPoint /= dataPoint.cwiseAbs().maxCoeff();
}

template<typename T, size_t N>
void CData<T, N>::GetScaleNormalizedDataPoint(size_t iIdx, InternalArraynD& dataPoint)
{
	CastDataPointTTypeToInternal<T, N>(m_mDatanD->at(iIdx), dataPoint);
	size_t iNumElems = dataPoint.num_elements();
	CCS_INTERNAL_TYPE maxVal = InternalArray1D::Map(dataPoint.data(), iNumElems).cwiseAbs().maxCoeff();
	for (size_t i = 0; i < iNumElems; ++i)
		dataPoint.data()[i] /= maxVal;
}




template <typename Derived1>
inline void NormalizeDataPoint(Eigen::MatrixBase<Derived1> const & dataPoint)
{
	const_cast< Eigen::MatrixBase<Derived1>& >(dataPoint) /= dataPoint.squaredNorm();
}

template<typename T, size_t N>
inline void NormalizeDataPoint(typename CData<T, N>::InternalArraynD& dataPoint)
{
	size_t iNumElems = dataPoint.num_elements();
	CCS_INTERNAL_TYPE nrm = CData<T, N>::InternalArray1D::Map(dataPoint.data(), iNumElems).squaredNorm();
	for (size_t i = 0; i < iNumElems; ++i)
		dataPoint.data()[i] /= nrm;
}

template <typename Derived1>
inline void ScaleNormalizeDataPoint(Eigen::MatrixBase<Derived1> const & dataPoint)
{
	const_cast<Eigen::MatrixBase<Derived1>&>(dataPoint) /= dataPoint.cwiseAbs().maxCoeff();
}

template<typename T, size_t N>
inline void ScaleNormalizeDataPoint(typename CData<T, N>::InternalArraynD& dataPoint)
{
	size_t iNumElems = dataPoint.num_elements();
	CCS_INTERNAL_TYPE nrm = CData<T, N>::InternalArray1D::Map(dataPoint.data(), iNumElems).cwiseAbs().maxCoeff();
	for (size_t i = 0; i < iNumElems; ++i)
		dataPoint.data()[i] /= nrm;
}

