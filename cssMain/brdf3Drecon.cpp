
#include "cssData/DataBRDFMerl3D.h"
#include "cssDict/DictTestEnsOrthnD.h"
#include "cssRecon/ReconEnsOrthnD.h"


#if defined(_MSC_VER) && _MSC_VER >= 1400 
#pragma warning(push) 
#pragma warning(disable:4996) 
#endif 



int main(int argc, char* argv[])
{
	omp_set_dynamic(0);
	omp_set_nested(1);

	double dTimeDelta = 0.0;

	typedef float DataType;



	std::vector<uint32_t> vPatchSize(2);
	vPatchSize[0] = 0;
	vPatchSize[1] = 0;
	std::vector<uint32_t> vSlidingDist(2);
	vSlidingDist[0] = 0;
	vSlidingDist[1] = 0;

	uint32_t iSparsity = 262;
	CCS_INTERNAL_TYPE threshold = 1e-16;

	CDictTestEnsOrthnD<DataType, 3>* pTest = NULL;
	CDictionary* pDict = NULL;
	CDataBRDFMerl3D<DataType>* pTestData = NULL;
	CReconEnsOrthnD<DataType, 3>* pRecon = NULL;
	CDataBRDFMerl3D<DataType>* pOutputData = NULL;

	std::string strTestSetDir;
	std::string strOutputDir;
	std::string strCoeffsDir;
	std::string strDictAddr;
	std::string strCoeffFileName;


	strOutputDir = "DataBRDF/Output/";
	strCoeffsDir = "DataBRDF/";
	strTestSetDir = strCoeffsDir + std::string("TestSet");
	strDictAddr = strCoeffsDir + std::string("DictEnsOrth3D.mat");
	strCoeffFileName = std::string("Coeffs.mat");



	pTestData = new CDataBRDFMerl3D<DataType>(90, 90, 180, "brdf4d");
	if (!pTestData->Init()) return -1;
	if (!pTestData->LoadFromDisk(strTestSetDir)) return -1;
	if (!pTestData->PrepareData(vPatchSize, CCS_PATCHTYPE_NOV, vSlidingDist, false)) return -1;
	std::cout << "Number of patches: " << pTestData->GetNumDataElems() << std::endl;

	pDict = new CDictionary();
	if (!pDict->Load(strDictAddr))	return -1;

	pTest = new CDictTestEnsOrthnD<DataType, 3>();
	if (!pTest->Init(true)) return -1;
	if (!pTest->Load(strCoeffsDir, strCoeffFileName, true))	return -1;
	pRecon = new CReconEnsOrthnD<DataType, 3>();
	pOutputData = new CDataBRDFMerl3D<DataType>(90, 90, 180, "brdf4d");
	pOutputData->ResizeLike(pTestData);
	if (!pRecon->Init())	return -1;
	std::cout << "Reconstructing" << std::endl;
	if (!pRecon->Reconstruct(pTest, pDict, iSparsity, threshold, pOutputData, dTimeDelta))	return -1;
	std::cout << "Done in " << dTimeDelta << " seconds" << std::endl;
	if (!pOutputData->AssembleData()) return -1;
	if (!pOutputData->WriteAssembled(strOutputDir)) return -1;



	CCS_INTERNAL_TYPE mse = MSE<DataType>(*pTestData->GetDatanD(), *pOutputData->GetDatanD());
	std::cout << "MSE: " << mse << std::endl;


	SAFE_DELETE(pOutputData);
	SAFE_DELETE(pRecon);
	SAFE_DELETE(pDict);
	SAFE_DELETE(pTest);
	SAFE_DELETE(pTestData);


	int x;
	std::cin >> x;

	return 0;
}