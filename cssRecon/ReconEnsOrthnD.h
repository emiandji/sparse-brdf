#pragma once

#include "Recon.h"
#include "ReconHelpers.h"



template <typename T, size_t N>
class CReconEnsOrthnD : public CRecon<T, N>
{
public:

	CReconEnsOrthnD();
	~CReconEnsOrthnD();

	virtual bool Init();
	virtual bool Reconstruct(const CDictTest<T, N>* pTest, const CDictionary* pDict, uint32_t iSparsity, CCS_INTERNAL_TYPE tThreshold,
		CData<T, N>* pOutput, double& dTimeDelta);
	virtual void CleanUp();

private:

};




#include "ReconEnsOrthnD.inl"
