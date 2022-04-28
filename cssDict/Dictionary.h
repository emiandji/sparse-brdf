#pragma once



#include "../cssUtil/VclMatlab.h"
#include "../cssUtil/LATools.h"

class CDictionary
{

public:

	enum DictType { DICT_INVALID = 0, DICT_ORTHO = 1, DICT_ORTHO_ENSEMBLE = 2, DICT_OVERCOMPLETE = 3, DICT_OVERCOMPLETE_ENSEMBLE = 4 };

    typedef Eigen::Matrix<CCS_INTERNAL_TYPE, Eigen::Dynamic, Eigen::Dynamic> DictElem;
	typedef std::vector<DictElem> Dict;
	typedef std::vector<Dict> DictEnsemble;

public:

	CDictionary()    { m_eDictType = DICT_INVALID; }
	CDictionary(const CDictionary& dict);
	virtual ~CDictionary()	{ Clear(); }

	void SetTrainTime(double dTrainTime)		{ m_dTrainTime = dTrainTime; }
	double GetTrainTime()						{ return m_dTrainTime; }

	// Memory Allocation. Must be called before using this class
	void AllocOrtho(std::vector<size_t> vSignalDims);
	void AllocOrthoEnsemble(std::vector<size_t> vSignalDims, size_t iNumDict);
	void AllocOvercomplete(std::vector<size_t> vSignalDims, std::vector<CCS_INTERNAL_TYPE> vOvFactor);									//Works with 1D signals only
	void AllocOvercompleteEnsemble(std::vector<size_t> vSignalDims, std::vector<CCS_INTERNAL_TYPE> vOvFactor, size_t iNumDict);		//Works with 1D signals only

    // Getters
	DictType GetDictType()												const	{ return m_eDictType; }
	size_t GetNumDict()													const	{ return m_Ensemble.size(); }
	size_t GetNumDictElems(size_t iDicIdx)								const	{ return m_Ensemble[iDicIdx].size(); }
	const DictEnsemble& GetEnsemble()									const	{ return m_Ensemble; }
	const Dict& GetDict(size_t iDicIdx)									const	{ return m_Ensemble[iDicIdx]; }
	const DictElem& GetDictElem(size_t iDicIdx, size_t iElemIdx)		const	{ return m_Ensemble[iDicIdx][iElemIdx]; }
	const Eigen::VectorXi& GetEnsHist()									const	{ return m_vEnsHist; }
	

    // Setters
	void SetDictElem(size_t iDicIdx, size_t iElemIdx, const DictElem& dictElem); 
	void SetDict(size_t iDicIdx, const Dict& dict);
	void SetEnsHist(const Eigen::VectorXi& vEnsHist);
	
	// Add/Remove
	int AddDict(const Dict& dict)									{ m_Ensemble.push_back(dict); }												//Returns the index of the added dictionary
	int AddDictElem(size_t iDictIdx, const DictElem& dictElem)		{ m_Ensemble[iDictIdx].push_back(dictElem); }								//Returns the index of the added dictionary
	void RemoveDict(size_t iDictIdx)								{ m_Ensemble.erase(m_Ensemble.begin() + iDictIdx); }						//Returns the index of the added dictionary
	void RemoveDictElem(size_t iDictIdx, size_t iElemIdx)			{ m_Ensemble[iDictIdx].erase(m_Ensemble[iDictIdx].begin() + iElemIdx); }	//Returns the index of the added dictionary

	// Conversion
	void GetDict1D(size_t iDictIdx, DictElem& dictElem);	//Converts and returns. Internal data is not changed
	void GetEnsemble1D(DictEnsemble& vDictEns);				//Converts and returns. Internal data is not changed
	void ConvertTo1D();										//Internal data is changed (only works on 2D ensemble for now)

	// File I/O
	bool Load(const boost::filesystem::path& pathEnsemble, const std::string& strDictName = "Dict", const std::string& strEnsHistName = "EnsHist");
	bool Save(const boost::filesystem::path& pathEnsemble, const std::string& strDictName = "Dict", const std::string& strEnsHistName = "EnsHist", const std::string& strVarTrainTime = "TrainTime");

	// Clean up
	void Clear();

private:

	DictType m_eDictType;

	double m_dTrainTime;			//Time took to train the dictionary

	//A collection of nD dictionaries, m_vDict[i][j] -> i : dictionary index, j : the index for jth matrix (j is the dimensionality of the dictionary)
	DictEnsemble m_Ensemble;

	//Histogram of the ensemble's membership based on the training set (this is computed outside of this class during training)
	Eigen::VectorXi m_vEnsHist;
};

