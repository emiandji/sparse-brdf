#pragma once


#include "../cssUtil/ProgressReporter.h"
#include "../cssData/Data.h"
#include "../cssDict/Dictionary.h"
#include "../cssDict/DictTest.h"

template <typename T, size_t N>
class CRecon
{
public:

	CRecon()	{}
	virtual ~CRecon()	{}

	virtual bool Init() = 0;
	virtual bool Reconstruct(const CDictTest<T, N>* pTest, const CDictionary* pDict, uint32_t iSparsity, CCS_INTERNAL_TYPE tThreshold,
		CData<T, N>* pOutput, double& dTimeDelta) = 0;
	virtual void CleanUp() = 0;

private:

};
