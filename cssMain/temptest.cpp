


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

	CDictTestEnsOrthnD<DataType, 3>* pTestB1 = NULL;
	CDictTestEnsOrthnD<DataType, 3>* pTestB2 = NULL;
	CDictionary* pDict = NULL;
	CDataBRDFMerl3D<DataType>* pTestDataB1 = NULL;
	CDataBRDFMerl3D<DataType>* pTestDataB2 = NULL;
	CDataBRDFMerl3D<DataType>* pOutputData = NULL;

	std::string strTestSetDirB1;
	std::string strTestSetDirB2;
	std::string strOutputDir;
	std::string strCoeffsDir;
	std::string strDataSetName;
	std::string strDictAddr;

	strTestSetDirB1 = "DataBRDF/interp/B1/";
	strTestSetDirB2 = "DataBRDF/interp/B2/";
	strCoeffsDir = "DataBRDF/interp/";
	strOutputDir = "DataBRDF/interp/Output/";

	uint32_t iSparsity = 1024;
	CCS_INTERNAL_TYPE threshold = 1e-16;
	std::vector<CCS_INTERNAL_TYPE> vAlpha = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };


	pTestDataB1 = new CDataBRDFMerl3D<DataType>(90, 90, 180, "brdf4d");
	if (!pTestDataB1->Init()) return -1;
	if (!pTestDataB1->LoadFromDisk(strTestSetDirB1)) return -1;
	if (!pTestDataB1->PrepareData(vPatchSize, CCS_PATCHTYPE_NOV, vSlidingDist, false)) return -1;
	std::cout << "Number of patches: " << pTestDataB1->GetNumDataElems() << std::endl;

	pTestDataB2 = new CDataBRDFMerl3D<DataType>(90, 90, 180, "brdf4d");
	if (!pTestDataB2->Init()) return -1;
	if (!pTestDataB2->LoadFromDisk(strTestSetDirB2)) return -1;
	if (!pTestDataB2->PrepareData(vPatchSize, CCS_PATCHTYPE_NOV, vSlidingDist, false)) return -1;
	std::cout << "Number of patches: " << pTestDataB2->GetNumDataElems() << std::endl;

	pDict = new CDictionary();
	if (!pDict->Load(strCoeffsDir))	return -1;

	std::vector<size_t> vSigDim;
	pTestDataB1->GetDataElemDim(vSigDim);
	size_t iNumDims = vSigDim.size();


	//Compress B1
	std::ostringstream strCoeffFileNameB1;
	strCoeffFileNameB1 << "Coeffs_BRDF_3D_B1" << strDataSetName << "_tau" << iSparsity << "_eps" << threshold << ".mat";
	pTestB1 = new CDictTestEnsOrthnD<DataType, 3>();
	if (!pTestB1->Init(true)) return -1;
	std::cout << "Computing coefficients" << std::endl;
	if (!pTestB1->Test(pTestDataB1, pDict, iSparsity, threshold, dTimeDelta)) return -1;
	std::cout << "Done in " << dTimeDelta << " seconds" << std::endl;
	std::cout << "Encoding and saving coefficients" << std::endl;
	if (!pTestB1->Save(strCoeffsDir, strCoeffFileNameB1.str(), 200, true, dTimeDelta))	return false;
	std::cout << "Done in " << dTimeDelta << " seconds" << std::endl;

	//Compress B2
	std::ostringstream strCoeffFileNameB2;
	strCoeffFileNameB2 << "Coeffs_BRDF_3D_B2" << strDataSetName << "_tau" << iSparsity << "_eps" << threshold << ".mat";
	pTestB2 = new CDictTestEnsOrthnD<DataType, 3>();
	if (!pTestB2->Init(true)) return -1;
	std::cout << "Computing coefficients" << std::endl;
	if (!pTestB2->Test(pTestDataB2, pDict, iSparsity, threshold, dTimeDelta)) return -1;
	std::cout << "Done in " << dTimeDelta << " seconds" << std::endl;
	std::cout << "Encoding and saving coefficients" << std::endl;
	if (!pTestB2->Save(strCoeffsDir, strCoeffFileNameB2.str(), 200, true, dTimeDelta))	return false;
	std::cout << "Done in " << dTimeDelta << " seconds" << std::endl;

	size_t iNumColorChannels = 3;
	//Get coefficients for B1 and B2
	std::vector<typename CData<DataType, 3>::InternalArraynD> vS1(3, CData<DataType, 3>::InternalArraynD(vSigDim, CCS_TENSOR_STORAGE_ORDER));
	std::vector<typename CData<DataType, 3>::InternalArraynD> vS2(3, CData<DataType, 3>::InternalArraynD(vSigDim, CCS_TENSOR_STORAGE_ORDER));
	for (size_t i = 0; i < iNumColorChannels; ++i)
	{
		SparseToDense(pTestB1->GetElemNZLoc(i), pTestB1->GetElemNZValue(i), vS1[i]);
		SparseToDense(pTestB2->GetElemNZLoc(i), pTestB2->GetElemNZValue(i), vS2[i]);
	}


	CDictionary R12;
	R12.AllocOrthoEnsemble(vSigDim, 3);		//One dictionary for each color channel
	for (size_t i = 0; i < iNumColorChannels; ++i)
	{
		for (size_t j = 0; j < iNumDims; ++j)
		{
			typename CDictionary::DictElem D1tD2 = pDict->GetDictElem(pTestB1->GetM()(i), j).transpose() * pDict->GetDictElem(pTestB2->GetM()(i), j);
			Eigen::JacobiSVD<CDictionary::DictElem> svd(D1tD2, Eigen::ComputeFullU | Eigen::ComputeFullV);
			typename CDictionary::DictElem tmpRot = svd.matrixU() * svd.matrixV().transpose();
			R12.SetDictElem(i, j, tmpRot);
		}
	}
	R12.Save(strCoeffsDir + "Rot.mat");


	for (size_t i = 0; i < vAlpha.size(); ++i)
	{
		pOutputData = new CDataBRDFMerl3D<DataType>(90, 90, 180, "brdf4d");
		pOutputData->ResizeLike(pTestDataB1);

		for (size_t j = 0; j < iNumColorChannels; ++j)
		{
			//Rotate S2 using R12
			typename CData<DataType, 3>::InternalArraynD tmp1(vSigDim, CCS_TENSOR_STORAGE_ORDER);
			typename CData<DataType, 3>::InternalArraynD S2_rot(vSigDim, CCS_TENSOR_STORAGE_ORDER);
			TensorCopy(vS2[j], tmp1);
			for (size_t k = 0; k < iNumDims; ++k)
			{
				TensorProdMat(tmp1, R12.GetDictElem(j, k), k, S2_rot);
				TensorCopy(S2_rot, tmp1);
			}
			//Interpolate coefficients
			typename CData<DataType, 3>::InternalArraynD S_interped(vSigDim, CCS_TENSOR_STORAGE_ORDER);
			CDictTest<DataType, 3>::DataPointCoeffVals::Map(S_interped.data(), S_interped.num_elements()) =
				(1.0 - vAlpha[i]) * CDictTest<DataType, 3>::DataPointCoeffVals::Map(vS1[j].data(), vS1[j].num_elements()) +
				vAlpha[i] * CDictTest<DataType, 3>::DataPointCoeffVals::Map(S2_rot.data(), S2_rot.num_elements());

			typename CData<DataType, 3>::InternalArraynD interped(vSigDim, CCS_TENSOR_STORAGE_ORDER);
			TensorCopy(S_interped, tmp1);
			for (size_t k = 0; k < iNumDims; ++k)
			{
				TensorProdMat(tmp1, pDict->GetDictElem(pTestB1->GetM()(j), k), k, interped);
				TensorCopy(interped, tmp1);
			}

			pOutputData->SetDataPoint(j, interped);

			//Use pTestB1 as container for interpolated coefficients
			typename CDictTest<DataType, 3>::DataPointCoeffLocs interpCoeffLocs;
			typename CDictTest<DataType, 3>::DataPointCoeffVals interpCoeffVals;

			size_t iNNZ = DenseToSparse(S_interped, interpCoeffLocs, interpCoeffVals, S_interped.num_elements());
			pTestB1->SetElemNNZ(j, iNNZ);
			pTestB1->SetElemNZLoc(j, interpCoeffLocs);
			pTestB1->SetElemNZValue(j, interpCoeffVals);
		}

		std::ostringstream strCoeffFileName;
		strCoeffFileName << "Coeffs_BRDF_3D_Interp" << strDataSetName << "_tau" << iSparsity << "_eps" << threshold << "_alpha" << vAlpha[i] << ".mat";
		if (!pTestB1->Save(strCoeffsDir, strCoeffFileName.str(), 200, true, dTimeDelta))	return false;

		if (!pOutputData->AssembleData()) return -1;
		std::ostringstream strOutputDirStream;
		strOutputDirStream << "alpha" << vAlpha[i] << "/";
		if (!pOutputData->WriteAssembled(strOutputDir + strOutputDirStream.str())) return -1;

		SAFE_DELETE(pOutputData);
	}


	SAFE_DELETE(pTestB1);
	SAFE_DELETE(pTestB2);
	SAFE_DELETE(pDict);
	SAFE_DELETE(pTestDataB1);
	SAFE_DELETE(pTestDataB2);
	SAFE_DELETE(pOutputData);


	int x;
	std::cin >> x;

	return 0;
}



