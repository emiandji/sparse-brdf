#include "Dictionary.h"







void CDictionary::SetDictElem(size_t iDicIdx, size_t iElemIdx, const DictElem& dictElem)
{
	for (size_t i = 0; i < dictElem.rows(); ++i)
		for (size_t j = 0; j < dictElem.cols(); ++j)
			(m_Ensemble[iDicIdx][iElemIdx])(i, j) = dictElem(i, j);
}


void CDictionary::SetDict(size_t iDicIdx, const Dict& dict)
{
	m_Ensemble[iDicIdx].clear();
	m_Ensemble[iDicIdx].resize(dict.size());
	for (size_t i = 0; i < dict.size(); ++i)
		m_Ensemble[iDicIdx][i].resize(dict[i].rows(), dict[i].cols());
	for (size_t i = 0; i < dict.size(); ++i)
		for (size_t j = 0; j < dict[i].rows(); ++j)
			for (size_t k = 0; k < dict[i].cols(); ++k)
				(m_Ensemble[iDicIdx][i])(j, k) = (dict[i])(j, k);
}

void CDictionary::SetEnsHist(const Eigen::VectorXi& vEnsHist)
{
	if (vEnsHist.size() != m_Ensemble.size())
	{
		m_vEnsHist.resize(0);
		return;
	}
	m_vEnsHist = vEnsHist;
}


CDictionary::CDictionary(const CDictionary& dict)
{
	m_eDictType = dict.m_eDictType;
	m_dTrainTime = dict.m_dTrainTime;

	m_Ensemble.resize(dict.m_Ensemble.size());
	for (size_t i = 0; i < dict.m_Ensemble.size(); ++i)
		m_Ensemble[i].resize(dict.m_Ensemble[i].size());
	for (size_t i = 0; i < dict.m_Ensemble.size(); ++i)
		for (size_t j = 0; j < dict.m_Ensemble[i].size(); ++j)
			m_Ensemble[i][j].resize(dict.m_Ensemble[i][j].rows(), dict.m_Ensemble[i][j].cols());

	for (size_t i = 0; i < m_Ensemble.size(); ++i)
		for (size_t j = 0; j < m_Ensemble[i].size(); ++j)
			for (size_t m = 0; m < m_Ensemble[i][j].rows(); ++m)
				for (size_t n = 0; n < m_Ensemble[i][j].cols(); ++n)
					(m_Ensemble[i][j])(m, n) = (dict.m_Ensemble[i][j])(m, n);
}


void CDictionary::AllocOrtho(std::vector<size_t> vSignalDims)
{
	m_eDictType = DICT_ORTHO;
	m_Ensemble.resize(1);
	m_Ensemble[0].resize(vSignalDims.size());
	for (size_t i = 0; i < vSignalDims.size(); ++i)
		m_Ensemble[0][i].setConstant(vSignalDims[i], vSignalDims[i], 0);
}


void CDictionary::AllocOrthoEnsemble(std::vector<size_t> vSignalDims, size_t iNumDict)
{
	m_eDictType = DICT_ORTHO_ENSEMBLE;
	m_Ensemble.resize(iNumDict);

	for (size_t i = 0; i < iNumDict; ++i)
		m_Ensemble[i].resize(vSignalDims.size());

	for (size_t i = 0; i < iNumDict; ++i)
		for (size_t j = 0; j < vSignalDims.size(); ++j)
			m_Ensemble[i][j].setConstant(vSignalDims[j], vSignalDims[j], 0);
}


void CDictionary::AllocOvercomplete(std::vector<size_t> vSignalDims, std::vector<CCS_INTERNAL_TYPE> vOvFactor)
{
	m_eDictType = DICT_OVERCOMPLETE;
	m_Ensemble.resize(1);
	m_Ensemble[0].resize(vSignalDims.size());
	for (size_t i = 0; i < vSignalDims.size(); ++i)
		m_Ensemble[0][i].setConstant(vSignalDims[i], size_t(vSignalDims[i] * vOvFactor[i]), 0);
}


void CDictionary::AllocOvercompleteEnsemble(std::vector<size_t> vSignalDims, std::vector<CCS_INTERNAL_TYPE> vOvFactor, size_t iNumDict)
{
	m_eDictType = DICT_OVERCOMPLETE_ENSEMBLE;
	m_Ensemble.resize(iNumDict);

	for (size_t i = 0; i < iNumDict; ++i)
		m_Ensemble[i].resize(vSignalDims.size());

	for (size_t i = 0; i < iNumDict; ++i)
		for (size_t j = 0; j < vSignalDims.size(); ++j)
			m_Ensemble[i][j].setConstant(vSignalDims[j], size_t(vSignalDims[j] * vOvFactor[j]), 0);
}


void CDictionary::GetDict1D(size_t iDictIdx, DictElem& dictElem)
{
	std::vector<size_t> vIndices;
	vIndices.resize(m_Ensemble[iDictIdx].size());
	std::iota(vIndices.begin(), vIndices.end(), 0);
	std::reverse(vIndices.begin(), vIndices.end());
	KronSelective(m_Ensemble[iDictIdx], vIndices, dictElem);
}


void CDictionary::GetEnsemble1D(DictEnsemble& vDictEns)
{
	vDictEns.resize(m_Ensemble.size());
	for (size_t i = 0; i < vDictEns.size(); ++i)
		vDictEns[i].resize(1);
	for (size_t i = 0; i < vDictEns.size(); ++i)
		GetDict1D(i, vDictEns[i][0]);
}


void CDictionary::ConvertTo1D()
{
	//tmp
#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < m_Ensemble.size(); ++i)
	{
		CDictionary::DictElem tmp;
		GetDict1D(i, tmp);
		m_Ensemble[i].resize(1);
		m_Ensemble[i][0].resize(tmp.rows(), tmp.cols());
		for (size_t m = 0; m < tmp.rows(); ++m)
			for (size_t n = 0; n < tmp.cols(); ++n)
				m_Ensemble[i][0](m, n) = tmp(m, n);
	}
}


bool CDictionary::Save(const boost::filesystem::path& pathDict, const std::string& strDictName, const std::string& strEnsHistName, const std::string& strVarTrainTime)
{
	if (m_Ensemble.empty())
	{
		return false;
	}

	// Open the ensemble file
	VclMatio writeDict;
	if (!writeDict.openForWriting(pathDict.string()))
	{
		std::cout << "Error: Cannot open file for dictionary writing." << std::endl;
		return false;
	}

	// Iterate and save the ensembles of dictionary elements to disk (one ensemble at a time)
	size_t iNumDim = m_Ensemble[0].size();
	for (size_t i = 0; i < iNumDim; ++i)
	{
		Dict ensTmp;
		ensTmp.resize(m_Ensemble.size());
		for (size_t j = 0; j < m_Ensemble.size(); ++j)
		{
			ensTmp[j].resize(m_Ensemble[j][i].rows(), m_Ensemble[j][i].cols());
			ensTmp[j] = m_Ensemble[j][i];
		}
		if (!writeDict.writeEigenMatrix3dNamed(strDictName + std::string("Dim") + boost::lexical_cast<std::string>(i + 1), ensTmp))
		{
			std::cout << "Error: Cannot write the dictionary." << std::endl;
			return false;
		}
	}

	//Save ensemble's histogram
	if (!writeDict.writeEigenVectorNamed(strEnsHistName, m_vEnsHist))
	{
		std::cout << "Error: Cannot write the dictionary." << std::endl;
		return false;
	}

	// Save the dimensionality of the ensemble
	if (!writeDict.writeValueNamed("numDicts", double(m_Ensemble.size())))
	{
		std::cout << "Error: Cannot write the dictionary." << std::endl;
		return false;
	}
	if (!writeDict.writeValueNamed("numDim", double(iNumDim)))
	{
		std::cout << "Error: Cannot write the dictionary." << std::endl;
		return false;
	}

	//Save training time
	if (!writeDict.writeValueNamed(strVarTrainTime, m_dTrainTime))
	{
		std::cout << "Error: Cannot write the dictionary." << std::endl;
		return false;
	}

	writeDict.close();

	return true;
}


bool CDictionary::Load(const boost::filesystem::path& pathDict, const std::string& strDictName, const std::string& strEnsHistName)
{
	// Open the ensemble file
	VclMatio readDict;
	if (!readDict.openForReading(pathDict.string()))
	{
		std::cout << "Error: Cannot open file for dictionary reading." << std::endl;
		return false;
	}

	double tNumDicts = 0;
	if (!readDict.readValueNamed("numDicts", tNumDicts))
	{
		std::cout << "Error: Cannot read the dictionary." << std::endl;
		return false;
	}
	size_t iNumDicts = (size_t)(tNumDicts);

	double tNumDim = 0;
	if (!readDict.readValueNamed("numDim", tNumDim))
	{
		std::cout << "Error: Cannot read the dictionary." << std::endl;
		return false;
	}
	size_t iNumDim = (size_t)(tNumDim);


	m_Ensemble.resize(iNumDicts);
	for (size_t i = 0; i < m_Ensemble.size(); ++i)
		m_Ensemble[i].resize(iNumDim);

	for (size_t i = 0; i < iNumDim; ++i)
	{
		Dict ensTmp;
		ensTmp.resize(iNumDicts);
		if (!readDict.readEigenMatrix3dNamed(strDictName + std::string("Dim") + boost::lexical_cast<std::string>(i + 1), ensTmp))
		{
			std::cout << "Error: Cannot read the dictionary." << std::endl;
			return false;
		}
		for (size_t j = 0; j < iNumDicts; ++j)
		{
			m_Ensemble[j][i].resize(ensTmp[j].rows(), ensTmp[j].cols());
			m_Ensemble[j][i] = ensTmp[j];
		}
	}

	if (!readDict.readEigenVectorNamed(strEnsHistName, m_vEnsHist))
	{
		std::cout << "Error: Cannot read the dictionary." << std::endl;
		return false;
	}

	readDict.close();

	return true;
}


void CDictionary::Clear()
{
	m_Ensemble.clear();
}

