



template <typename T>
CHuffmanDict<T>::CHuffmanDict()
{
	m_cRoot = NULL;
}


template <typename T>
CHuffmanDict<T>::CHuffmanDict(typename std::map<T, boost::dynamic_bitset<unsigned long> >& _mHuffmanDictMap)
{
	m_cRoot = NULL;
}


template <typename T>
CHuffmanDict<T>::~CHuffmanDict()
{
	SAFE_DELETE(m_cRoot);
}


template <typename T>
void CHuffmanDict<T>::CreateTree(std::vector<T>& _Sig, std::vector<float>& _fProb)
{
	std::priority_queue<ClNode*, std::vector<ClNode*>, SNodeCompare> Hufftree_queue;
	size_t uniqueSig = _Sig.size();
	for (size_t i = 0; i < uniqueSig; i++)
	{
		if (_fProb[i] > 1e-12)
			Hufftree_queue.push(new CLeafNode<T>(_fProb[i], _Sig[i]));
	}

	while (Hufftree_queue.size() > 1)
	{
		ClNode * childR = Hufftree_queue.top();
		Hufftree_queue.pop();

		ClNode * childL = Hufftree_queue.top();
		Hufftree_queue.pop();

		ClNode * parent = new CInternalNode(childR, childL);

		Hufftree_queue.push(parent);
	}
	m_cRoot = Hufftree_queue.top();

}




template <typename T>
void CHuffmanDict<T>::GenerateDict(const ClNode* _node, const boost::dynamic_bitset<unsigned long>& _prefix)
{
	if (const CLeafNode<T>* lf = dynamic_cast<const CLeafNode<T>*>(_node))
	{
		m_mHuffmanDictMap[lf->m_fC] = _prefix;
	}
	else if (const CInternalNode* in = dynamic_cast<const CInternalNode*>(_node))
	{
		boost::dynamic_bitset<unsigned long> leftPrefix = _prefix;
		leftPrefix.push_back(true);
		GenerateDict(in->m_cLeft, leftPrefix);

		boost::dynamic_bitset<unsigned long> rightPrefix = _prefix;
		rightPrefix.push_back(false);
		GenerateDict(in->m_cRight, rightPrefix);
	}
}

template <typename T>
void CHuffmanDict<T>::Encode(const std::vector<T> &_testSignal, boost::dynamic_bitset<unsigned long>& _Huffencoded)
{

	int idxCode = 1;
	for (int i = 0; i < _testSignal.size(); i++)
	{
		// For each signal value, search sequentially through the dictionary to
		// find the code for the given signal
		boost::dynamic_bitset<unsigned long> tempcode;
		typename std::map<T, boost::dynamic_bitset<unsigned long> >::const_iterator it;

		it = m_mHuffmanDictMap.find(_testSignal[i]);

		if (it != m_mHuffmanDictMap.end())
		{
			tempcode.clear();
			tempcode = it->second;
			injectLoopDyn(_Huffencoded, tempcode);
		}
	}
}

template <typename T>
void CHuffmanDict<T>::Decode(const boost::dynamic_bitset<unsigned long>& _encodedSig, size_t iOrigSignalSize, vector<T>& _decodedSig)
{
	std::map<std::string, T> tmpMap;

	for (typename std::map<T, boost::dynamic_bitset<unsigned long> >::const_iterator it = m_mHuffmanDictMap.begin(); it != m_mHuffmanDictMap.end(); ++it)
	{
		string tmpstr;
		boost::to_string(it->second, tmpstr);
		tmpMap[tmpstr] = it->first;
	}

	size_t i = 0;
	size_t j = 0;
	boost::dynamic_bitset<unsigned long> tempCode;
	_decodedSig.resize(iOrigSignalSize);
	while (i < _encodedSig.size())
	{
		tempCode.push_back(_encodedSig[i]);
		string tmpstr;
		boost::to_string(tempCode, tmpstr);
		typename std::map<std::string, T>::const_iterator it = tmpMap.find(tmpstr);

		if (it != tmpMap.end())
		{
			_decodedSig[j++] = it->second;
			tempCode.clear();
		}
		i++;
	}
}

template <typename T>
bool CHuffmanDict<T>::ReadDict(const std::string _filename)
{
	m_mHuffmanDictMap.clear();
	ifstream inFile(_filename, std::ios::binary);

	if (inFile.is_open())
	{
		while (!inFile.eof())
		{
			T key;
			if (!inFile.read(reinterpret_cast<char*>(&key), sizeof(T)))
				break;

			//read size of coded dictionary
			size_t len;
			inFile.read((char*)&len, sizeof(size_t));

			//read coded dictionary
			unsigned long n;
			inFile.read(reinterpret_cast<char*>(&n), sizeof(n));
			boost::dynamic_bitset<unsigned long> coded(len, n);

			m_mHuffmanDictMap[key] = coded;
		}

		inFile.close();
		return true;

	}

	return false;
}

template <typename T>
bool CHuffmanDict<T>::WriteDict(const std::string _filename)
{
	ofstream output(_filename, std::ios::binary);

	if (output.is_open())
	{
		for (typename std::map<T, boost::dynamic_bitset<unsigned long> >::const_iterator it = m_mHuffmanDictMap.begin(); it != m_mHuffmanDictMap.end(); ++it)
		{
			//write  key value
			T key = it->first;

			output.write(reinterpret_cast<const char *>(&key), sizeof(T));

			//write size of coded dictionary
			size_t len;
			len = it->second.size();
			output.write(reinterpret_cast<const char *>(&len), sizeof(size_t));

			//write the coded dictionary 
			unsigned long n = it->second.to_ulong();
			output.write(reinterpret_cast<const char*>(&n), sizeof(n));
		}

		output.close();
		return true;
	}

	return false;
}


template <typename T>
bool CHuffmanDict<T>::WriteCodedSignal(boost::dynamic_bitset<unsigned long> &_encodedSig, const std::string _filename, size_t iOrigSignalSize)
{
	ofstream output(_filename, ios::binary);
	if (output.is_open())
	{
		boost::archive::binary_oarchive oa(output);

		size_t lenCoded;
		lenCoded = _encodedSig.size();

		oa << iOrigSignalSize;
		oa << lenCoded;
		oa << _encodedSig;

		output.close();
		return true;
	}
	return false;


}


template <typename T>
bool CHuffmanDict<T>::ReadCodedSignal(boost::dynamic_bitset<unsigned long> &_encodedSig, const std::string _filename, size_t& iOrigSignalSize)
{
	_encodedSig.clear();
	ifstream ifs(_filename, ios::binary);
	if (ifs.is_open())
	{
		boost::archive::binary_iarchive ia(ifs);
		size_t lenCoded;

		ia >> iOrigSignalSize;
		ia >> lenCoded;
		_encodedSig.resize(lenCoded);
		ia >> _encodedSig;

		ifs.close();
		return true;
	}
	return false;

}
