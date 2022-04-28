#pragma once



#include "Data.h"
#include "../cssUtil/VclMatlab.h"
#include "../cssUtil/Patchifier.h"
#include "../cssUtil/LATools.h"
#include "../cssUtil/ImgQlty.h"


template <typename T>
class CDataBRDFMerl3D : public CData<T, 3>
{

public:

	CDataBRDFMerl3D(uint32_t iTheta_h, uint32_t iTheta_d, uint32_t iPhi_d, const std::string& strBRDFVarName);
	virtual ~CDataBRDFMerl3D();

	virtual const std::vector<boost::filesystem::path>& GetFileNames()	const { return m_vFileNames; }

	virtual size_t GetNumDataElems() const;
	virtual void GetDataElemDim(std::vector<size_t>& vDims) const;

	template <typename U> void ResizeLike(const CData<U, 3>* pData);
	template <typename U> void CopyPropsFrom(const CData<U, 3>* pData);

	virtual bool Init();
	virtual bool LoadFromDisk(const std::string& strDataFolder);
	virtual bool PrepareData(const std::vector<uint32_t>& vPatchSize, PatchType ePatchType, const std::vector<uint32_t>& vSlidingDis, bool bFreeLoadedData);
	virtual void CreateMask(MaskRandMode eRandMode, float fNonZeroRatio);
	virtual void Resize(size_t iNumPoints, const std::vector<size_t> vDims);
	virtual void ConvertTo1D(bool bDeleteOld);
	virtual void ConvertBack(bool bDeleteOld);
	virtual void Clamp(T minVal, T maxVal);
	virtual void RemoveZeros(CCS_INTERNAL_TYPE thresh);
	virtual bool AssembleData();
	virtual void CalcQuality(const CData<T, 3>* pOther, ReconQualityMetric eMetric, std::vector<CCS_INTERNAL_TYPE>& vQlty);
	virtual bool WriteAssembled(const std::string& strOutputFolder);
	virtual void CleanUp();

	const std::vector<typename CData<T, 4>::ArraynD>& GetBRDFs() { return m_vBRDF; }
	const std::vector<boost::filesystem::path>& GetFileNames() { return m_vFileNames; }

protected:

	uint32_t m_iTheta_h;
	uint32_t m_iTheta_d;
	uint32_t m_iPhi_d;
	std::string m_strBRDFVarName;

	//Each BRDF is of size [iTheta_h, iTheta_d, iPhi_d, 3]
	std::vector<typename CData<T, 4>::ArraynD> m_vBRDF;

private:


	//File names for BRDFs
	std::vector<boost::filesystem::path> m_vFileNames;

	//Number of patches per BRDF
	std::vector<size_t> m_vNumPatches;
};



#include "DataBRDFMerl3D.inl"
