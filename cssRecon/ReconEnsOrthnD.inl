


template <typename T, size_t N>
CReconEnsOrthnD<T, N>::CReconEnsOrthnD()
{

}


template <typename T, size_t N>
CReconEnsOrthnD<T, N>::~CReconEnsOrthnD()
{
	CleanUp();
}


template <typename T, size_t N>
bool CReconEnsOrthnD<T, N>::Init()
{
	CleanUp();

	return true;
}


template <typename T, size_t N>
bool CReconEnsOrthnD<T, N>::Reconstruct(const CDictTest<T, N>* pTest,
	const CDictionary* pDict,
	uint32_t iSparsity, 
	CCS_INTERNAL_TYPE tThreshold, 
	CData<T, N>* pOutput, 
	double& dTimeDelta)
{
	size_t iN = pTest->GetNZValue().size();
	size_t iLength = std::accumulate(pTest->GetShape().begin(), pTest->GetShape().end(), size_t(1), std::multiplies<size_t>());
	if (pTest->GetNZValue().size() != pTest->GetNZLoc().size())
	{
		std::cerr << "ERROR: Size of coefficient values and locations don't match." << std::endl;
		return false;
	}
	std::vector<size_t> vDims = pTest->GetShape();
	std::vector<size_t> vDimsOutput;
	pOutput->GetDataElemDim(vDimsOutput);
	if (!std::equal(vDims.begin(), vDims.end(), vDimsOutput.begin()) || iN != pOutput->GetNumDataElems())
	{
		std::cerr << "ERROR: output data incompatible with coefficients" << std::endl;
		return false;
	}

	for (size_t i = 0; i < pDict->GetNumDict(); ++i)
	{
		for (size_t j = 0; j < pDict->GetNumDictElems(i); ++j)
		{
			if (pDict->GetDictElem(i, j).rows() != vDims[j])
			{
				std::cerr << "ERROR: invalid dictionary for reconstruction" << std::endl;
				return false;
			}
		}
	}

	if ((pTest->GetM().rows() != iN) || (pTest->GetNNZ().size() != iN))
	{
		std::cerr << "ERROR: invalid membership or NNZ vector" << std::endl;
		return false;
	}

	double dLastTime = omp_get_wtime();
	CProgressReporter reporter(iN, 0);
#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < iN; ++i)
	{
		typename CData<T, N>::InternalArraynD dataPoint(vDims, CCS_TENSOR_STORAGE_ORDER);
		HOSVDReconSp(pTest->GetNZLoc().at(i), pTest->GetNZValue().at(i), pDict->GetDict(pTest->GetM()(i, 0)), dataPoint);
		pOutput->SetDataPoint(i, dataPoint);
		reporter.Update();
	}
	reporter.Done();
	dTimeDelta = omp_get_wtime() - dLastTime;

	return true;
}


template <typename T, size_t N>
void CReconEnsOrthnD<T, N>::CleanUp()
{

}
