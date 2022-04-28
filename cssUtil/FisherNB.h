

//Fisher Natural Breaks Classification
//This is a slightly modified version of 
//http://wiki.objectvision.nl/index.php/CalcNaturalBreaksCode


#pragma once



#include "defs.h"





// helper struct JenksFisher
// captures the intermediate data and methods for the calculation of Natural Class Breaks.
template <typename T>
class JenksFisher  
{
public:
	JenksFisher(const std::vector<std::pair<T, size_t> >& vcpc, size_t k);

	// Gets sum of weighs for elements b..e.
	T GetW(size_t b, size_t e);

	// Gets sum of weighed values for elements b..e
	T GetWV(size_t b, size_t e);

	// Gets the Squared Mean for elements b..e, multiplied by weight.
	// Note that n*mean^2 = sum^2/n when mean := sum/n
	T GetSSM(size_t b, size_t e);

	// finds CB[i+m_NrCompletedRows] given that 
	// the result is at least bp+(m_NrCompletedRows-1)
	// and less than          ep+(m_NrCompletedRows-1)
	// Complexity: O(ep-bp) <= O(m)
	size_t FindMaxBreakIndex(size_t i, size_t bp, size_t ep);

	// find CB[i+m_NrCompletedRows]
	// for all i>=bi and i<ei given that
	// the results are at least bp+(m_NrCompletedRows-1)
	// and less than            ep+(m_NrCompletedRows-1)
	// Complexity: O(log(ei-bi)*Max((ei-bi),(ep-bp))) <= O(m*log(m))
	void CalcRange(size_t bi, size_t ei, size_t bp, size_t ep);


	// complexity: O(m*log(m)*k)
	void CalcAll();


	size_t                   m_M, m_K, m_BufSize;
	std::vector<std::pair<T, size_t> > m_CumulValues;
	std::vector<T> m_PrevSSM;
	std::vector<T> m_CurrSSM;
	std::vector<size_t>               m_CB;
	std::vector<size_t>::iterator     m_CBPtr;
	size_t                  m_NrCompletedRows;
};


template <typename T>
size_t GetTotalCount(const  std::vector<std::pair<T, size_t> >& vcpc);

template <typename T>
void GetCountsDirect(std::vector<std::pair<T, size_t> >& vcpc, const T* values, size_t size);

template <typename T>
void MergeToLeft(std::vector<std::pair<T, size_t> >& vcpcLeft, const std::vector<std::pair<T, size_t> >& vcpcRight, std::vector<std::pair<T, size_t> >& vcpcDummy);

template <typename T>
struct ValueCountPairContainerArray : std::vector<typename std::vector<std::pair<T, size_t> > >
{
	void resize(size_t k);
	void GetValueCountPairs(std::vector<std::pair<T, size_t> >& vcpc, const T* values, size_t size, unsigned int nrUsedContainers);
};

template <typename T>
void GetValueCountPairs(std::vector<std::pair<T, size_t> >& vcpc, const T* values, size_t n);

template <typename T>
void ClassifyJenksFisherFromValueCountPairs(std::vector<T>& breaksArray, size_t k, const std::vector<std::pair<T, size_t> >& vcpc);



template <typename T>
class CFisherNB
{
private:

	//finds the closest codeword to the signal in the codebook
	void quantiz(T sig, const std::vector<T>& vPartition, int& iIdx);
	void quantiz(const std::vector<T>& vSig, const std::vector<T>& vPartition, std::vector<T>& vCodebook, std::vector<int>& vIdx);

public:

	CFisherNB();
	virtual ~CFisherNB();

	void Cluster(const std::vector<T>& vData, size_t iNumBreaks, std::vector<int>& vMemb, std::vector<T>& vPartitions, 
		std::vector<T>& vCodebook, double& dDistortion, double& dTimeDelta);
};






#include "FisherNB.inl"