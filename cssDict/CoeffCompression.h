#pragma once

#include "../cssUtil/HuffmanDict.h"
#include "../cssData/Data.h"
#include "../cssUtil/FisherNB.h"
#include "../cssUtil/ProbTools.h"

namespace fs = ::boost::filesystem;


template <typename T, size_t N>
class CCoeffCompression
{
	
public:

	CCoeffCompression();
	virtual ~CCoeffCompression();

	bool Init(size_t iNumQuantPoints, 
		const fs::path& strAddr,
		const fs::path& strCoeffsIntFileName = "CoeffsInt.bin",
		const fs::path& strCoeffsFracFileName = "CoeffsFrac.bin",
		const fs::path& strNNZFileName = "NNZ.bin",
		const fs::path& strMembFileName = "Memb.bin",
		const fs::path& strCoeffsIntDictFileName = "CoeffsIntDict.bin",
		const fs::path& strCoeffsFracDictFileName = "CoeffsFracDict.bin",
		const fs::path& strNNZDictFileName = "NNZDict.bin",
		const fs::path& strMembDictFileName = "MembDict.bin");

	bool EncodeCoeffs(const std::vector< typename CData<T, N>::Array1D >& pCoeffs);
	bool EncodeCoeffs(const std::vector< typename CData<T, N>::Array2D >& pCoeffs);
	bool EncodeCoeffs(const std::vector< typename CData<T, N>::ArraynD >& pCoeffs);

	bool DecodeCoeffs(std::vector< typename CData<T, N>::Array1D >& pCoeffs);
	bool DecodeCoeffs(std::vector< typename CData<T, N>::Array2D >& pCoeffs);
	bool DecodeCoeffs(std::vector< typename CData<T, N>::ArraynD >& pCoeffs);

	bool EncodeNNZ(const Eigen::VectorXi& vNNZ);
	bool EncodeMemb(const Eigen::VectorXi& vMemb);

	bool DecodeNNZ(Eigen::VectorXi& vNNZ);
	bool DecodeMemb(Eigen::VectorXi& vMemb);

private:

	fs::path m_strAddr;
	fs::path m_strCoeffsIntFileName;
	fs::path m_strCoeffsFracFileName;
	fs::path m_strNNZFileName;
	fs::path m_strMembFileName;
	fs::path m_strCoeffsIntDictFileName;
	fs::path m_strCoeffsFracDictFileName;
	fs::path m_strNNZDictFileName;
	fs::path m_strMembDictFileName;
	size_t m_iNumQuantPoints;

	bool EncodeTask(const std::vector<int>& vSignal, const fs::path& pathCode, const fs::path& pathDict);
	bool EncodeTask(const std::vector<T>& vSignal, int iNumQuantPoints, const fs::path& pathCode, const fs::path& pathDict);

	bool DecodeTask(const fs::path& pathCode, const fs::path& pathDict, std::vector<int>& vSignal);
	bool DecodeTask(const fs::path& pathCode, const fs::path& pathDict, std::vector<T>& vSignal);
};



#include "CoeffCompression.inl"