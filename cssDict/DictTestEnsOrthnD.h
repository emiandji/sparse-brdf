#pragma once

#include "DictTest.h"


template <typename T, size_t N>
class CDictTestEnsOrthnD : public CDictTest<T, N>
{

public:
    
	CDictTestEnsOrthnD(const std::string& strCoeffsInfoFileName = "SInfo", const std::string& strVarMemb = "M",
		const std::string& strVarCoeffVal = "NZ", const std::string& strVarCoeffLoc = "LNZ", const std::string& strVarNNZ = "NNZ",
		const std::string& strWeightCos = "MERL-cosweight_matrix.mat", const std::string& strWeightCosTi = "MERL-costiweight_matrix.mat");
	virtual ~CDictTestEnsOrthnD();

	virtual bool Init(bool bStoreCoeffs);
	virtual bool Test(const CData<T, N>* pData, const CDictionary* pDict, uint32_t iSparsity, CCS_INTERNAL_TYPE tThreshold, double& dTimeDelta);
	virtual bool Save(const boost::filesystem::path& strPath, const std::string& strFileName, size_t iNumQuantPoints, bool bSaveRaw, double& dTimeDelta);
	virtual bool Load(const boost::filesystem::path& strPath, const std::string& strFileName, bool bLoadRaw);
	virtual bool Concatenate(const std::vector< CDictTest<T, N>* >& vToConcat, const std::vector<size_t>& vShape);
	virtual void CleanUp();

	std::string m_strWeightCos;
	std::string m_strWeightCosTi;

	typename CData<T, N>::InternalArraynD m_WeightCos;
	typename CData<T, N>::InternalArraynD m_WeightCosTi;
};

#include "DictTestEnsOrthnD.inl"
