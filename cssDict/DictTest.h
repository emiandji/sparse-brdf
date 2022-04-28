#pragma once

#include "../cssUtil/ProgressReporter.h"
#include "Dictionary.h"
#include "CoeffCompression.h"



template <typename T, size_t N>
class CDictTest
{

public:

	typedef typename CData<CCS_INTERNAL_NZ_LOC_TYPE, N>::Array2D DataPointCoeffLocs;
	typedef typename CData<CCS_INTERNAL_TYPE, N>::Array1D DataPointCoeffVals;
	typedef std::vector<DataPointCoeffVals> CoeffVals;
	typedef std::vector<DataPointCoeffLocs> CoeffLocs;
    
	CDictTest(const std::string& strCoeffsInfoFileName, const std::string& strVarMemb, const std::string& strVarCoeffVal,
		const std::string& strVarCoeffLoc, const std::string& strVarNNZ)
	{
		m_strCoeffsInfoFileName = strCoeffsInfoFileName;
		m_strVarMemb = strVarMemb;
		m_strVarCoeffVal = strVarCoeffVal;
		m_strVarCoeffLoc = strVarCoeffLoc;
		m_strVarNNZ = strVarNNZ;
		m_strVarDim = "dims";
	}
	virtual ~CDictTest() {}

	//Get functions
	virtual const Eigen::MatrixXi&		GetM()			const	{ return m_mMemb; }
	virtual const Eigen::VectorXi&		GetNNZ()		const	{ return m_vNNZ; }
	virtual const std::vector<size_t>&	GetShape()		const	{ return m_vShape; }
	virtual const CoeffLocs&			GetNZLoc()		const	{ return m_vNZLoc; }
	virtual const CoeffVals&			GetNZValue()	const	{ return m_vNZValue; }
	virtual DataPointCoeffLocs&			GetElemNZLoc(size_t iDataElemIdx) { return m_vNZLoc[iDataElemIdx]; }
	virtual DataPointCoeffVals&			GetElemNZValue(size_t iDataElemIdx) { return m_vNZValue[iDataElemIdx]; }

	//Set functions
	virtual void SetM(const Eigen::MatrixXi& mM)			{ m_mMemb = mM; }
	virtual void SetNNZ(const Eigen::VectorXi& vNNZ)		{ m_vNNZ = vNNZ; }
	virtual void SetNZLoc(const CoeffLocs& vNZLoc)			{ m_vNZLoc = vNZLoc; }
	virtual void SetNZValue(const CoeffVals& vNZValue)		{ m_vNZValue = vNZValue; }
	virtual void SetElemNZLoc(size_t iDataElemIdx, const DataPointCoeffLocs& NZLoc) { m_vNZLoc[iDataElemIdx] = NZLoc; }
	virtual void SetElemNZValue(size_t iDataElemIdx, const DataPointCoeffVals& NZValue) { m_vNZValue[iDataElemIdx] = NZValue; }

	virtual void SetElemNNZ(size_t iDataElemIdx, int iVal) { m_vNNZ(iDataElemIdx) = iVal; }

	virtual bool GetStoreCoeffs() { return m_bStoreCoeffs; }

	virtual bool Init(bool bStoreCoeffs) = 0;
	virtual bool Test(const CData<T, N>* pData, const CDictionary* pDict, uint32_t iSparsity, CCS_INTERNAL_TYPE tThreshold, double& dTimeDelta) = 0;
	virtual bool Save(const boost::filesystem::path& strPath, const std::string& strFileName, size_t iNumQuantPoints, bool bSaveRaw, double& dTimeDelta) = 0;
	virtual bool Load(const boost::filesystem::path& strPath, const std::string& strFileName, bool bLoadRaw) = 0;
	virtual bool Concatenate(const std::vector< CDictTest<T, N>* >& vToConcat, const std::vector<size_t>& vShape) = 0;
	virtual void CleanUp() = 0;
    
protected:

	std::string m_strCoeffsInfoFileName;
	std::string m_strVarMemb;
	std::string m_strVarCoeffVal;
	std::string m_strVarCoeffLoc;
	std::string m_strVarNNZ;
	std::string m_strVarDim;

	//Exemplar membership vector
	Eigen::MatrixXi m_mMemb;
    
	//Number of non-zero coefficients per data point
	Eigen::VectorXi m_vNNZ;
	
    //Nonzero coefficients for each data point stacked in a vector (m_vNNZ can be used to extract coeffs for each data point)
	CoeffVals m_vNZValue;

	//Location of nonzero coefficients for each data point stacked (each matrix is of size (NNZ, numDims))
	CoeffLocs m_vNZLoc;

	//This stores the dimensionality of each data point (because we only store nonzeros and the dimensionality gets lost in the process)
	std::vector<size_t> m_vShape;

	bool m_bStoreCoeffs;
};
