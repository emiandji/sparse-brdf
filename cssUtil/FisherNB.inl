






template <typename T>
JenksFisher<T>::JenksFisher(const std::vector<std::pair<T, size_t> >& vcpc, size_t k) : m_M(vcpc.size())
, m_K(k)
, m_BufSize(vcpc.size() - (k - 1))
, m_PrevSSM(m_BufSize)
, m_CurrSSM(m_BufSize)
, m_CB(m_BufSize * (m_K - 1))
, m_CBPtr()
{
	m_CumulValues.reserve(vcpc.size());
	T cwv = 0;
	size_t cw = 0, w;

	for (size_t i = 0; i != m_M; ++i)
	{
		assert(!i || vcpc[i].first > vcpc[i - 1].first); // PRECONDITION: the value sequence must be strictly increasing

		w = vcpc[i].second;
		assert(w > 0); // PRECONDITION: all weights must be positive

		cw += w;
		assert(cw > w || !i); // No overflow? No loss of precision?

		cwv += w * vcpc[i].first;
		m_CumulValues.push_back(std::pair<T, size_t>(cwv, cw));
		if (i < m_BufSize)
			m_PrevSSM[i] = cwv * cwv / cw; // prepare SSM for first class. Last (k-1) values are omitted
	}
}

template <typename T>
T JenksFisher<T>::GetW(size_t b, size_t e)
{
	assert(b);    // First element always belongs to class 0, thus queries should never include it.
	assert(b <= e);
	assert(e < m_M);

	T res = m_CumulValues[e].second;
	res -= m_CumulValues[b - 1].second;
	return res;
}

template <typename T>
T JenksFisher<T>::GetWV(size_t b, size_t e)
{
	assert(b);
	assert(b <= e);
	assert(e < m_M);

	T res = m_CumulValues[e].first;
	res -= m_CumulValues[b - 1].first;
	return res;
}

template <typename T>
T JenksFisher<T>::GetSSM(size_t b, size_t e)
{
	T res = GetWV(b, e);
	return res * res / GetW(b, e);
}

template <typename T>
size_t JenksFisher<T>::FindMaxBreakIndex(size_t i, size_t bp, size_t ep)
{
	assert(bp < ep);
	assert(bp <= i);
	assert(ep <= i + 1);
	assert(i < m_BufSize);
	assert(ep <= m_BufSize);

	T minSSM = m_PrevSSM[bp] + GetSSM(bp + m_NrCompletedRows, i + m_NrCompletedRows);
	size_t foundP = bp;
	while (++bp < ep)
	{
		T currSSM = m_PrevSSM[bp] + GetSSM(bp + m_NrCompletedRows, i + m_NrCompletedRows);
		if (currSSM > minSSM)
		{
			minSSM = currSSM;
			foundP = bp;
		}
	}
	m_CurrSSM[i] = minSSM;
	return foundP;
}

template <typename T>
void JenksFisher<T>::CalcRange(size_t bi, size_t ei, size_t bp, size_t ep)
{
	assert(bi <= ei);

	assert(ep <= ei);
	assert(bp <= bi);

	if (bi == ei)
		return;
	assert(bp < ep);

	size_t mi = (bi + ei) / 2;
	size_t mp = FindMaxBreakIndex(mi, bp, std::min<size_t>(ep, mi + 1));

	assert(bp <= mp);
	assert(mp < ep);
	assert(mp <= mi);

	// solve first half of the sub-problems with lower 'half' of possible outcomes
	CalcRange(bi, mi, bp, std::min<size_t>(mi, mp + 1));

	m_CBPtr[mi] = mp; // store result for the middle element.

					  // solve second half of the sub-problems with upper 'half' of possible outcomes
	CalcRange(mi + 1, ei, mp, ep);
}

template <typename T>
void JenksFisher<T>::CalcAll()
{
	if (m_K >= 2)
	{
		m_CBPtr = m_CB.begin();
		for (m_NrCompletedRows = 1; m_NrCompletedRows < m_K - 1; ++m_NrCompletedRows)
		{
			CalcRange(0, m_BufSize, 0, m_BufSize); // complexity: O(m*log(m))

			m_PrevSSM.swap(m_CurrSSM);
			m_CBPtr += m_BufSize;
		}
	}
}


template <typename T>
size_t GetTotalCount(const std::vector<std::pair<T, size_t> >& vcpc)
{
	size_t sum = 0;
	typename std::vector<std::pair<T, size_t> >::const_iterator
		i = vcpc.begin(),
		e = vcpc.end();
	for (sum = 0; i != e; ++i)
		sum += (*i).second;
	return sum;
}


template <typename T>
void GetCountsDirect(std::vector<std::pair<T, size_t> >& vcpc, const T* values, size_t size)
{
	assert(size <= FISHER_NB_BUFFER_SIZE);
	assert(size > 0);
	assert(vcpc.empty());

	T buffer[FISHER_NB_BUFFER_SIZE];

	std::copy(values, values + size, buffer);
	std::sort(buffer, buffer + size);

	T currValue = buffer[0];
	size_t     currCount = 1;
	for (size_t index = 1; index != size; ++index)
	{
		if (currValue < buffer[index])
		{
			vcpc.push_back(std::pair<T, size_t>(currValue, currCount));
			currValue = buffer[index];
			currCount = 1;
		}
		else
			++currCount;
	}
	vcpc.push_back(std::pair<T, size_t>(currValue, currCount));
}


template <typename T>
void MergeToLeft(std::vector<std::pair<T, size_t> >& vcpcLeft, const std::vector<std::pair<T, size_t> >& vcpcRight, std::vector<std::pair<T, size_t> >& vcpcDummy)
{
	assert(vcpcDummy.empty());
	vcpcDummy.swap(vcpcLeft);
	vcpcLeft.resize(vcpcRight.size() + vcpcDummy.size());

	std::merge(vcpcRight.begin(), vcpcRight.end(), vcpcDummy.begin(), vcpcDummy.end(), vcpcLeft.begin(), 
		[](const std::pair<T, size_t>& lhs, const std::pair<T, size_t>& rhs) {return lhs.first < rhs.first;});

	typename std::vector<std::pair<T, size_t> >::iterator
		currPair = vcpcLeft.begin(),
		lastPair = vcpcLeft.end();


	typename std::vector<std::pair<T, size_t> >::iterator index = currPair + 1;
	while (index != lastPair && currPair->first < index->first)
	{
		currPair = index;
		++index;
	}

	T currValue = currPair->first;
	size_t     currCount = currPair->second;
	for (; index != lastPair; ++index)
	{
		if (currValue < index->first)
		{
			*currPair++ = std::pair<T, size_t>(currValue, currCount);
			currValue = index->first;
			currCount = index->second;
		}
		else
			currCount += index->second;
	}
	*currPair++ = std::pair<T, size_t>(currValue, currCount);
	vcpcLeft.erase(currPair, lastPair);

	vcpcDummy.clear();
}


template <typename T>
void ValueCountPairContainerArray<T>::resize(size_t k)
{
	assert(this->capacity() >= k);
	while (this->size() < k)
	{
		this->push_back(std::vector<std::pair<T, size_t> >());
		this->back().reserve(FISHER_NB_BUFFER_SIZE);
	}
}

template <typename T>
void ValueCountPairContainerArray<T>::GetValueCountPairs(std::vector<std::pair<T, size_t> >& vcpc, const T* values, size_t size, unsigned int nrUsedContainers)
{
	assert(vcpc.empty());
	if (size <= FISHER_NB_BUFFER_SIZE)
		GetCountsDirect(vcpc, values, size);
	else
	{
		resize(nrUsedContainers + 2);

		unsigned int m = size / 2;

		GetValueCountPairs(vcpc, values, m, nrUsedContainers);
		GetValueCountPairs(this->begin()[nrUsedContainers], values + m, size - m, nrUsedContainers + 1);

		MergeToLeft(vcpc, this->begin()[nrUsedContainers], this->begin()[nrUsedContainers + 1]);
		this->begin()[nrUsedContainers].clear();
	}
	assert(GetTotalCount(vcpc) == size);
}



template <typename T>
void GetValueCountPairs(std::vector<std::pair<T, size_t> >& vcpc, const T* values, size_t n)
{
	vcpc.clear();

	if (n)
	{
		ValueCountPairContainerArray<T> vcpca;
		// max nr halving is log2(max cardinality / FISHER_NB_BUFFER_SIZE); max cardinality is size_t(-1)
		vcpca.reserve(3 + 8 * sizeof(size_t) - 10);
		vcpca.GetValueCountPairs(vcpc, values, n, 0);

		assert(vcpc.size());
	}
}



template <typename T>
void ClassifyJenksFisherFromValueCountPairs(std::vector<T>& breaksArray, size_t k, const std::vector<std::pair<T, size_t> >& vcpc)
{
	breaksArray.resize(k);
	size_t m = vcpc.size();

	assert(k <= m); // PRECONDITION

	if (!k)
		return;

	JenksFisher<T> jf(vcpc, k);

	if (k > 1)
	{
		jf.CalcAll();

		size_t lastClassBreakIndex = jf.FindMaxBreakIndex(jf.m_BufSize - 1, 0, jf.m_BufSize);

		while (--k)
		{
			breaksArray[k] = vcpc[lastClassBreakIndex + k].first;
			assert(lastClassBreakIndex < jf.m_BufSize);
			if (k > 1)
			{
				jf.m_CBPtr -= jf.m_BufSize;
				lastClassBreakIndex = jf.m_CBPtr[lastClassBreakIndex];
			}
		}
		assert(jf.m_CBPtr == jf.m_CB.begin());
	}
	assert(k == 0);
	breaksArray[0] = vcpc[0].first; // break for the first class is the minimum of the dataset.
}


template <typename T>
CFisherNB<T>::CFisherNB()
{

}

template<typename T>
CFisherNB<T>::~CFisherNB()
{

}


template<typename T>
void CFisherNB<T>::quantiz(T sig, const std::vector<T>& vPartition,	int& iIdx)
{
	for (size_t i = 0; i < vPartition.size(); ++i)
	{
		if (sig <= vPartition[i])
		{
			iIdx = i;
			break;
		}
	}
	if (sig > vPartition[vPartition.size() - 1])
		iIdx = vPartition.size();
	//dDistortion = pow(sig - vCodebook(iIdx + 1), T(2.0));
}


template<typename T>
void CFisherNB<T>::quantiz(const std::vector<T>& vSig, const std::vector<T>& vPartition, std::vector<T>& vCodebook, std::vector<int>& vIdx)
{
	for (size_t i = 0; i < vSig.size(); ++i)
		quantiz(vSig[i], vPartition, vIdx[i]);

	vCodebook.resize(vPartition.size() + 1, 0);
	std::vector<size_t> vAccum(vCodebook.size(), 0);

	for (size_t i = 0; i < vSig.size(); ++i)
	{
		vCodebook[vIdx[i]] += vSig[i];
		vAccum[vIdx[i]]++;
	}
	for (size_t i = 0; i < vCodebook.size(); ++i)
		vCodebook[i] /= T(vAccum[i]);
}


template<typename T>
void CFisherNB<T>::Cluster(const std::vector<T>& vData, 
	size_t iNumBreaks,
	std::vector<int>& vMemb, 
	std::vector<T>& vPartitions,
	std::vector<T>& vCodebook,
	double& dDistortion, 
	double& dTimeDelta)
{
	double dLastTime = omp_get_wtime();

	std::vector<std::pair<T, size_t> > sortedUniqueValueCounts;
	GetValueCountPairs(sortedUniqueValueCounts, &vData[0], vData.size());

	ClassifyJenksFisherFromValueCountPairs(vPartitions, iNumBreaks, sortedUniqueValueCounts);

	vMemb.resize(vData.size());
	quantiz(vData, vPartitions, vCodebook, vMemb);

	dTimeDelta = omp_get_wtime() - dLastTime;
}
