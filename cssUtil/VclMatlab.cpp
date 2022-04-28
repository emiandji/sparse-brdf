


#include "VclMatlab.h"

//
// NOTE:
// This is not the prettiest code ever!
// Because the Matio library is rather lame, one needs rewind the variable iterator
// pointer each time one accesses a variable.  Therefore, each function is required
// to rewind the pointer when it starts to make sure all is well.
//





VclMatio::VclMatio()
    : m_matfp(NULL)
    , m_filename("")
{ }


VclMatio::VclMatio(const std::string& filename, bool forReading)
    : m_filename(filename)
{
    if (forReading)
        openForReading(m_filename);
    else
        openForWriting(m_filename);
}


VclMatio::~VclMatio()
{
    close();
}


// Reading
bool VclMatio::openForReading(const std::string& filename)
{
    close();
    m_matfp = Mat_Open(filename.c_str(), MAT_ACC_RDONLY);
    if (m_matfp == NULL)
        return false;
    return true;
}


// Writing
bool VclMatio::openForWriting(const std::string& filename)
{
    close();
    m_matfp = Mat_CreateVer(filename.c_str(), NULL, MAT_FT_MAT73);
    if (m_matfp == NULL)
    {
		std::cerr << "Unable to create mat file: " << filename << std::endl;
		return false;
	}
	return true;
}


// Info
std::vector<std::string> VclMatio::variableNames()
{
    std::vector<std::string> vNames;
    
    matvar_t* matvar;
    Mat_Rewind(m_matfp);
    while ((matvar = Mat_VarReadNextInfo(m_matfp)) != NULL) 
    {
        vNames.push_back(matvar->name);
        Mat_VarFree(matvar);
        matvar = NULL;
    }
    return vNames;
}


bool VclMatio::readVariableNamedDim(const std::string& name, std::vector<size_t>& vDim)
{
	if (!matFileOpen(m_matfp))
		return false;

	matvar_t* matvarRead = NULL;

	// Read the data into memory
	Mat_Rewind(m_matfp);
	matvarRead = Mat_VarReadInfo(m_matfp, name.c_str());

	vDim.clear();
	vDim.assign(matvarRead->dims, matvarRead->dims + matvarRead->rank);

	Mat_VarFree(matvarRead);
	return true;
}


bool VclMatio::isValue(const std::string& name)
{
    if (!matFileOpen(m_matfp))
        return false;
    
    // Make sure the variable exists and is the proper rank
    bool good = true;
    Mat_Rewind(m_matfp);
    matvar_t* matvar = Mat_VarReadInfo(m_matfp, name.c_str());
    if (matvar == NULL)
    {
        std::cerr << "Variable named '" << name << "' does not exist." << std::endl;
        good = false;
    }
    if (matvar->rank == 2 && (matvar->dims[0] != 1 || matvar->dims[1] != 1))
        good = false;
    if (matvar->rank >= 3)
        good = false;
    
    Mat_VarFree(matvar);
    return good;
}


bool VclMatio::isRowVec(const std::string& name)
{
    if (!matFileOpen(m_matfp))
        return false;

    // Make sure the variable exists and is the proper rank
    bool good = true;
    Mat_Rewind(m_matfp);
    matvar_t* matvar = Mat_VarReadInfo(m_matfp, name.c_str());
    if (matvar == NULL)
    {
        std::cerr << "Variable named '" << name << "' does not exist." << std::endl;
        good = false;
    }
    if (matvar->rank != 2 || matvar->dims[0] != 1)
        good = false;

    Mat_VarFree(matvar);
    return good;
}


bool VclMatio::isColVec(const std::string& name)
{
    if (!matFileOpen(m_matfp))
        return false;

    // Make sure the variable exists and is the proper rank
    bool good = true;
    Mat_Rewind(m_matfp);
    matvar_t* matvar = Mat_VarReadInfo(m_matfp, name.c_str());
    if (matvar == NULL)
    {
        std::cerr << "Variable named '" << name << "' does not exist." << std::endl;
        good = false;
    }
    if (matvar->rank != 2 || matvar->dims[1] != 1)
        good = false;

    Mat_VarFree(matvar);
    return good;
}


bool VclMatio::isArray2d(const std::string& name)
{
    if (!matFileOpen(m_matfp))
        return false;
    
    // Make sure the variable exists and is the proper rank
    bool good = true;
    Mat_Rewind(m_matfp);
    matvar_t* matvar = Mat_VarReadInfo(m_matfp, name.c_str());
    if (matvar == NULL)
    {
        std::cerr << "Variable named '" << name << "' does not exist." << std::endl;
        good = false;
    }
	if (matvar->rank != 2)
		good = false;
    if (matvar->rank == 2 && matvar->dims[0] == 1 && matvar->dims[1] == 1)
        good = false;
    
    Mat_VarFree(matvar);
    return good;
}


bool VclMatio::isArray3d(const std::string& name)
{
    if (!matFileOpen(m_matfp))
        return false;
    
    // Make sure the variable exists and is the proper rank
    bool good = true;
    Mat_Rewind(m_matfp);
    matvar_t* matvar = Mat_VarReadInfo(m_matfp, name.c_str());
    if (matvar == NULL)
    {
        std::cerr << "Variable named '" << name << "' does not exist." << std::endl;
        good = false;
    }
    if (matvar->rank != 3)
        good = false;
	if (std::accumulate(matvar->dims, matvar->dims + matvar->rank, 0) == matvar->rank)
		good = false;
    
    Mat_VarFree(matvar);
    return good;
}


bool VclMatio::isArraynd(const std::string& name)
{
	if (!matFileOpen(m_matfp))
		return false;

	// Make sure the variable exists and is the proper rank
	bool good = true;
	Mat_Rewind(m_matfp);
	matvar_t* matvar = Mat_VarReadInfo(m_matfp, name.c_str());
	if (matvar == NULL)
	{
		std::cerr << "Variable named '" << name << "' does not exist." << std::endl;
		good = false;
	}
	if (matvar->rank < 3)
		good = false;
	if (std::accumulate(matvar->dims, matvar->dims + matvar->rank, 0) == matvar->rank)
		good = false;

	Mat_VarFree(matvar);
	return good;
}

matio_types VclMatio::GetTypeFromClass(matio_classes cls)
{
	switch (cls) {
	case MAT_C_DOUBLE:
		return MAT_T_DOUBLE;
	case MAT_C_SINGLE:
		return MAT_T_SINGLE;
#ifdef HAVE_MAT_INT64_T
	case MAT_C_INT64:
		return MAT_T_INT64;
#endif
#ifdef HAVE_MAT_UINT64_T
	case MAT_C_UINT64:
		return MAT_T_UINT64;
#endif
	case MAT_C_INT32:
		return MAT_T_INT32;
	case MAT_C_UINT32:
		return MAT_T_UINT32;
	case MAT_C_INT16:
		return MAT_T_INT16;
	case MAT_C_UINT16:
		return MAT_T_UINT16;
	case MAT_C_INT8:
		return MAT_T_INT8;
	case MAT_C_CHAR:
		return MAT_T_UINT8;
	case MAT_C_UINT8:
		return MAT_T_UINT8;
	case MAT_C_CELL:
		return MAT_T_CELL;
	case MAT_C_STRUCT:
		return MAT_T_STRUCT;
	default:
		return MAT_T_UNKNOWN;
	}
}


// Closing
bool VclMatio::close()
{
    Mat_Close(m_matfp);
    m_matfp = NULL;
    return true;
}


// Private functions
bool VclMatio::matFileOpen(const mat_t* mat)
{
    if (mat == NULL)
    {
        std::cerr << "No mat file opened for reading." << std::endl;
        return false;
    }
    return true;
}




