

template<typename T, size_t N>
inline CCoeffCompression<T, N>::CCoeffCompression()
{

}

template<typename T, size_t N>
inline CCoeffCompression<T, N>::~CCoeffCompression()
{

}

template<typename T, size_t N>
inline bool CCoeffCompression<T, N>::Init(size_t iNumQuantPoints,
	const fs::path& strAddr,
	const fs::path& strCoeffsIntFileName,
	const fs::path& strCoeffsFracFileName,
	const fs::path& strNNZFileName,
	const fs::path& strMembFileName,
	const fs::path& strCoeffsIntDictFileName,
	const fs::path& strCoeffsFracDictFileName,
	const fs::path& strNNZDictFileName,
	const fs::path& strMembDictFileName)
{
	m_iNumQuantPoints = iNumQuantPoints;
	m_strAddr = strAddr;


	m_strCoeffsIntFileName = strCoeffsIntFileName;
	m_strCoeffsFracFileName = strCoeffsFracFileName;
	m_strNNZFileName = strNNZFileName;
	m_strMembFileName = strMembFileName;

	m_strCoeffsIntDictFileName = strCoeffsIntDictFileName;
	m_strCoeffsFracDictFileName = strCoeffsFracDictFileName;
	m_strNNZDictFileName = strNNZDictFileName;
	m_strMembDictFileName = strMembDictFileName;

	return true;
}



template <typename T, size_t N>
bool CCoeffCompression<T, N>::EncodeTask(const std::vector<int>& vSignal, const fs::path& pathCode, const fs::path& pathDict)
{
	int iNodeID = 0;

	std::vector<float> vfProb_int;
	std::vector<int> vi_int_symbol;
	std::vector<int> Signal_inttmp(vSignal);
	std::vector<float> vfedges;

	std::sort(Signal_inttmp.begin(), Signal_inttmp.end());
	std::unique_copy(Signal_inttmp.begin(), Signal_inttmp.end(), std::back_inserter(vi_int_symbol));
	Signal_inttmp.clear();

	vfProb_int.resize(vi_int_symbol.size(), 0.0f);

	vfedges.push_back((float)vi_int_symbol[0] - 0.5f);
	for (int i = 0; i < vi_int_symbol.size(); i++)
		vfedges.push_back((float)vi_int_symbol[i] + 0.5f);

	HistCount<int, float>(vSignal, vfedges, vfProb_int);

	CHuffmanDict<int> huff;

	huff.CreateTree(vi_int_symbol, vfProb_int);
	huff.GenerateDict(huff.GetRoot(), boost::dynamic_bitset<unsigned long>());
	if (iNodeID == 0)
	{
		if (!huff.WriteDict(pathDict.string()))
		{
			std::cerr << "ERROR: Something went wrong trying to save the dictionary!" << std::endl;
			return false;
		}
	}

	boost::dynamic_bitset<unsigned long> huffencoded;
	huff.Encode(vSignal, huffencoded);
	if (iNodeID == 0)
	{
		if (!huff.WriteCodedSignal(huffencoded, pathCode.string(), vSignal.size()))
		{
			std::cerr << "ERROR: Something went wrong trying to save the encoded signal!" << std::endl;
			return false;
		}
	}

	return true;
}



template<typename T, size_t N>
bool CCoeffCompression<T, N>::EncodeTask(const std::vector<T>& vSignal, int iNumQuantPoints, const fs::path & pathCode, const fs::path & pathDict)
{
	int iNodeID = 0;

	double dTimeDelta = 0;
	double dDummy;
	CFisherNB<T> fisher;
	std::vector<int> viIndx;
	std::vector<T> vPartitions;
	std::vector<T> vCodebook;
	fisher.Cluster(vSignal, iNumQuantPoints, viIndx, vPartitions, vCodebook, dDummy, dTimeDelta);
	vPartitions.clear();

	/// convert from Eigen to vector
	std::vector<int> viSortedIndx = viIndx;
	
	std::vector<int> vUniqueIndex;
	std::sort(viSortedIndx.begin(), viSortedIndx.end());
	std::unique_copy(viSortedIndx.begin(), viSortedIndx.end(), std::back_inserter(vUniqueIndex));
	viSortedIndx.clear();

	std::vector<float> edges;
	std::vector<float> vfProb;
	edges.push_back((float)vUniqueIndex[0] - 0.5);
	for (size_t i = 0; i < vUniqueIndex.size(); i++)
		edges.push_back((float)vUniqueIndex[i] + 0.5);
	vfProb.resize(vUniqueIndex.size(), 0.0);

	HistCount<int, float>(viIndx, edges, vfProb);

	std::vector<T> vfQuantizSig;
	vfQuantizSig.resize(vSignal.size());

	for (size_t i = 0; i < vfQuantizSig.size(); i++)
		vfQuantizSig[i] = vCodebook[viIndx[i]];

// 	/// convert from Eigen to vector
// 	std::vector<T> vSymbol_frac;
// 	vSymbol_frac.assign(centroids.data(), centroids.data() + centroids.size());

	CHuffmanDict<T> huff;
	huff.CreateTree(vCodebook, vfProb);
	huff.GenerateDict(huff.GetRoot(), boost::dynamic_bitset<unsigned long>());
	if (iNodeID == 0)
	{
		if (!huff.WriteDict(pathDict.string()))
		{
			std::cerr << "ERROR: Something went wrong trying to save the dictionary!" << std::endl;
			return false;
		}
	}

	boost::dynamic_bitset<unsigned long> huffencoded_frac;
	huff.Encode(vfQuantizSig, huffencoded_frac);
	if (iNodeID == 0)
	{
		if (!huff.WriteCodedSignal(huffencoded_frac, pathCode.string(), vSignal.size()))
		{
			std::cerr << "ERROR: Something went wrong trying to save the encoded signal!" << std::endl;
			return false;
		}
	}

	return true;
}



template <typename T, size_t N>
bool CCoeffCompression<T, N>::DecodeTask(const fs::path& pathCode, const fs::path& pathDict, std::vector<int>& vSignal)
{
	boost::dynamic_bitset<unsigned long> huffencoded;
	CHuffmanDict<int> huff;
	if (!huff.ReadDict(pathDict.string()))
	{
		std::cerr << "ERROR: failed to load the dictionary" << std::endl;
		return false;
	}
	size_t iSize = 0;
	if (!huff.ReadCodedSignal(huffencoded, pathCode.string(), iSize))
	{
		std::cerr << "ERROR: failed to load the encoded signal" << std::endl;
		return false;
	}
	huff.Decode(huffencoded, iSize, vSignal);

	return true;
}



template <typename T, size_t N>
bool CCoeffCompression<T, N>::DecodeTask(const fs::path& pathCode, const fs::path& pathDict, std::vector<T>& vSignal)
{
	boost::dynamic_bitset<unsigned long> huffencoded;
	CHuffmanDict<T> huff;
	if (!huff.ReadDict(pathDict.string()))
	{
		std::cerr << "ERROR: failed to load the dictionary" << std::endl;
		return false;
	}
	size_t iSize = 0;
	if (!huff.ReadCodedSignal(huffencoded, pathCode.string(), iSize))
	{
		std::cerr << "ERROR: failed to load the encoded signal" << std::endl;
		return false;
	}
	huff.Decode(huffencoded, iSize, vSignal);

	return true;
}



template<typename T, size_t N>
bool CCoeffCompression<T, N>::EncodeCoeffs(const std::vector<typename CData<T, N>::Array1D>& pCoeffs)
{
	if (pCoeffs.empty())
		return false;

	size_t iN = pCoeffs.size();
	size_t iLength = pCoeffs.at(0).size();

	//Separate integer and fractional part of the signal
	std::vector<T> vSignal;

	vSignal.resize(iN*iLength);
	for (size_t i = 0; i < iN; ++i)
	{
		for (size_t j = 0; j < iLength; ++j)
		{
			size_t iIdx = j + i*iLength;
			vSignal[iIdx] = (pCoeffs[i].data())[j];
		}
	}

	//Encode the fractional part of coefficients
	if (!EncodeTask(vSignal, m_iNumQuantPoints, m_strAddr / m_strCoeffsFracFileName, m_strAddr / m_strCoeffsFracDictFileName))
	{
		std::cerr << "ERROR: failed to encode the fractional part of the signal" << std::endl;
		return false;
	}

	return true;
}




template <typename T, size_t N>
bool CCoeffCompression<T, N>::EncodeCoeffs(const std::vector< typename CData<T, N>::Array2D >& pCoeffs)
{
	if (pCoeffs.empty())
		return false;

	size_t iN = pCoeffs.size();
	size_t iLength = pCoeffs.at(0).size();

	//Separate integer and fractional part of the signal
	std::vector<T> vSignal;

	vSignal.resize(iN*iLength);
	for (size_t i = 0; i < iN; ++i)
	{
		for (size_t j = 0; j < iLength; ++j)
		{
			size_t iIdx = j + i*iLength;
			vSignal[iIdx] = (pCoeffs[i].data())[j];
		}
	}

	//Encode the fractional part of coefficients
	if (!EncodeTask(vSignal, m_iNumQuantPoints, m_strAddr / m_strCoeffsFracFileName, m_strAddr / m_strCoeffsFracDictFileName))
	{
		std::cerr << "ERROR: failed to encode the fractional part of the signal" << std::endl;
		return false;
	}

	return true;
}



template <typename T, size_t N>
bool CCoeffCompression<T, N>::EncodeCoeffs(const std::vector< typename CData<T, N>::ArraynD >& pCoeffs)
{
	if (pCoeffs.empty())
		return false;

	typedef typename CData<T, N>::ArraynD::size_type szType;
	size_t iN = pCoeffs.size();
	size_t iLength = std::accumulate(pCoeffs.at(0).shape() , pCoeffs.at(0).shape() + pCoeffs.at(0).num_dimensions(), 1, std::multiplies<szType>());

	//Separate integer and fractional part of the signal
	std::vector<T> vSignal;

	vSignal.resize(iN*iLength);
	for (size_t i = 0; i < iN; ++i)
	{
		for (size_t j = 0; j < iLength; ++j)
		{
			size_t iIdx = j + i*iLength;
			vSignal[iIdx] = (pCoeffs[i].data())[j];
		}
	}

	//Encode the fractional part of coefficients
	if (!EncodeTask(vSignal, m_iNumQuantPoints, m_strAddr / m_strCoeffsFracFileName, m_strAddr / m_strCoeffsFracDictFileName))
	{
		std::cerr << "ERROR: failed to encode the fractional part of the signal" << std::endl;
		return false;
	}

	return true;
}




template <typename T, size_t N>
bool CCoeffCompression<T, N>::DecodeCoeffs(std::vector< typename CData<T, N>::Array1D >& pCoeffs)
{
	size_t iN = pCoeffs.size();
	size_t iLength = pCoeffs.at(0).size();

	std::vector<T> vDecoded;

	if (!DecodeTask(m_strAddr / m_strCoeffsFracFileName, m_strAddr / m_strCoeffsFracDictFileName, vDecoded))
	{
		std::cerr << "ERROR: failed to decode fractional part of coefficients" << std::endl;
		return false;
	}

	if ((iN*iLength) != vDecoded.size())
	{
		std::cerr << "ERROR: number of decoded coefficients does not match the requested number" << std::endl;
		return false;
	}

	for (size_t i = 0; i < iN; ++i)
	{
		for (size_t j = 0; j < iLength; ++j)
		{
			size_t iIdx = j + i*iLength;
			(pCoeffs[i].data())[j] = vDecoded[iIdx];
		}
	}

	return true;
}



template <typename T, size_t N>
bool CCoeffCompression<T, N>::DecodeCoeffs(std::vector< typename CData<T, N>::Array2D >& pCoeffs)
{
	size_t iN = pCoeffs.size();
	size_t iLength = pCoeffs.at(0).size();

	std::vector<T> vDecoded;

	if (!DecodeTask(m_strAddr / m_strCoeffsFracFileName, m_strAddr / m_strCoeffsFracDictFileName, vDecoded))
	{
		std::cerr << "ERROR: failed to decode fractional part of coefficients" << std::endl;
		return false;
	}

	if ((iN*iLength) != vDecoded.size())
	{
		std::cerr << "ERROR: number of decoded coefficients does not match the requested number" << std::endl;
		return false;
	}

	for (size_t i = 0; i < iN; ++i)
	{
		for (size_t j = 0; j < iLength; ++j)
		{
			size_t iIdx = j + i*iLength;
			(pCoeffs[i].data())[j] = vDecoded[iIdx];
		}
	}

	return true;
}




template <typename T, size_t N>
bool CCoeffCompression<T, N>::DecodeCoeffs(std::vector< typename CData<T, N>::ArraynD >& pCoeffs)
{
	typedef typename CData<T, N>::ArraynD::size_type szType;

	size_t iN = pCoeffs.size();
	size_t iLength = std::accumulate(pCoeffs.at(0).shape(), pCoeffs.at(0).shape() + pCoeffs.at(0).num_dimensions(), 1, std::multiplies<szType>());

	std::vector<T> vDecoded;

	if (!DecodeTask(m_strAddr / m_strCoeffsFracFileName, m_strAddr / m_strCoeffsFracDictFileName, vDecoded))
	{
		std::cerr << "ERROR: failed to decode fractional part of coefficients" << std::endl;
		return false;
	}

	if ((iN*iLength) != vDecoded.size())
	{
		std::cerr << "ERROR: number of decoded coefficients does not match the requested number" << std::endl;
		return false;
	}

	for (size_t i = 0; i < iN; ++i)
	{
		for (size_t j = 0; j < iLength; ++j)
		{
			size_t iIdx = j + i*iLength;
			(pCoeffs[i].data())[j] = vDecoded[iIdx];
		}
	}

	return true;
}




template<typename T, size_t N>
bool CCoeffCompression<T, N>::EncodeNNZ(const Eigen::VectorXi& vNNZ)
{
	std::vector<int> vTmpNNZ;
	vTmpNNZ.assign(vNNZ.data(), vNNZ.data() + vNNZ.size());
	if (!EncodeTask(vTmpNNZ, m_strAddr / m_strNNZFileName, m_strAddr / m_strNNZDictFileName))
	{
		std::cerr << "ERROR: failed to encode NNZ vector" << std::endl;
		return false;
	}
	return true;
}




template<typename T, size_t N>
bool CCoeffCompression<T, N>::EncodeMemb(const Eigen::VectorXi& vMemb)
{
	std::vector<int> vTmpMemb;
	vTmpMemb.assign(vMemb.data(), vMemb.data() + vMemb.size());
	if (!EncodeTask(vTmpMemb, m_strAddr / m_strMembFileName, m_strAddr / m_strMembDictFileName))
	{
		std::cerr << "ERROR: failed to encode the membership vector" << std::endl;
		return false;
	}
	return true;
}




template<typename T, size_t N>
bool CCoeffCompression<T, N>::DecodeNNZ(Eigen::VectorXi& vNNZ)
{
	std::vector<int> vTmpNNZ;
	if (!DecodeTask(m_strAddr / m_strNNZFileName, m_strAddr / m_strNNZDictFileName, vTmpNNZ))
	{
		std::cerr << "ERROR: failed to decode NNZ vector" << std::endl;
		return false;
	}
	vNNZ.resize(vTmpNNZ.size());
	vNNZ = Eigen::VectorXi::Map(&vTmpNNZ[0], vTmpNNZ.size());

// 	boost::dynamic_bitset<unsigned long> huffencoded;
// 	CHuffmanDict<int> huffInt;
// 	if (!huffInt.ReadDict(boost::filesystem::path(m_strAddr + m_strNNZDictFileName)))
// 	{
// 		std::cerr << "ERROR: failed to load the dictionary" << std::endl;
// 		return false;
// 	}
// 	size_t iSize = 0;
// 	if (!huffInt.ReadCodedSignal(huffencoded, pathCode.string(), iSize))
// 	{
// 		std::cerr << "ERROR: failed to load the encoded signal" << std::endl;
// 		return false;
// 	}
// 	std::vector<int> vTmpNNZ;
// 	huffInt.Decode(huffencoded, iSize, vTmpNNZ);
// 	vNNZ.resize(vTmpNNZ.size());
// 	vNNZ = Eigen::VectorXi::Map(&vTmpNNZ[0], vTmpNNZ.size());

	return true;
}




template<typename T, size_t N>
bool CCoeffCompression<T, N>::DecodeMemb(Eigen::VectorXi& vMemb)
{
	std::vector<int> vTmpMemb;
	if (!DecodeTask(m_strAddr / m_strMembFileName, m_strAddr / m_strMembDictFileName, vTmpMemb))
	{
		std::cerr << "ERROR: failed to decode Memb vector" << std::endl;
		return false;
	}
	vMemb.resize(vTmpMemb.size());
	vMemb = Eigen::VectorXi::Map(&vTmpMemb[0], vTmpMemb.size());


// 	boost::dynamic_bitset<unsigned long> huffencoded;
// 	CHuffmanDict<int> huffInt;
// 	if (!huffInt.ReadDict(pathDict.string()))
// 	{
// 		std::cerr << "ERROR: failed to load the dictionary" << std::endl;
// 		return false;
// 	}
// 	size_t iSize = 0;
// 	if (!huffInt.ReadCodedSignal(huffencoded, pathCode.string(), iSize))
// 	{
// 		std::cerr << "ERROR: failed to load the encoded signal" << std::endl;
// 		return false;
// 	}
// 	std::vector<int> vTmpMemb;
// 	huffInt.Decode(huffencoded, iSize, vTmpMemb);
// 	vMemb.resize(vTmpMemb.size());
// 	vMemb = Eigen::VectorXi::Map(&vTmpMemb[0], vTmpMemb.size());

	return true;
}



// template<typename T, size_t N>
// void CCoeffComp<T, N>::HistCount(std::vector<float>& vprob_, const std::vector<int>& _vsig, const std::vector<float>& _vedges)
// {
// 	if (vprob_.empty()) //add exception here
// 		return;
// 	size_t nsizeSig = _vsig.size();
// 	size_t nsizeEdges = _vedges.size();
// 	std::vector<size_t> vn(vprob_.size(), 0);
// 
// 
// 	size_t sum = 0;
// 	for (size_t i = 0; i < nsizeSig; i++)
// 	{
// 		for (size_t j = 0; j < nsizeEdges - 1; j++)
// 		{
// 			if ((float)_vsig[i] > _vedges[j] && (float)_vsig[i] < _vedges[j + 1])
// 			{
// 				vn[j]++;
// 				sum++;
// 			}
// 		}
// 	}
// 
// 	for (size_t i = 0; i < vprob_.size(); i++)
// 	{
// 
// 		vprob_[i] = (float)vn[i] / (float)sum;
// 	}
// }

