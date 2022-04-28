//
// VCL header
//

#pragma once

#include "defs.h"
#include <matio.h>


template <class matio_type_t>
matio_classes GetMatioClassFromT()
{
	throw std::runtime_error("Intrinsic type not supported by the Matio library data type.");
}

template <>
inline matio_classes GetMatioClassFromT<int8_t>()
{
	return (MAT_C_INT8);
}

template <>
inline matio_classes GetMatioClassFromT<uint8_t>()
{
	return (MAT_C_UINT8);
}

template <>
inline matio_classes GetMatioClassFromT<int16_t>()
{
	return (MAT_C_INT16);
}

template <>
inline matio_classes GetMatioClassFromT<uint16_t>()
{
	return (MAT_C_UINT16);
}

template <>
inline matio_classes GetMatioClassFromT<int32_t>()
{
	return (MAT_C_INT32);
}

template <>
inline matio_classes GetMatioClassFromT<uint32_t>()
{
	return (MAT_C_UINT32);
}

template <>
inline matio_classes GetMatioClassFromT<int64_t>()
{
	return (MAT_C_INT64);
}

template <>
inline matio_classes GetMatioClassFromT<uint64_t>()
{
	return (MAT_C_UINT64);
}

template <>
inline matio_classes GetMatioClassFromT<float>()
{
	return (MAT_C_SINGLE);
}

template <>
inline matio_classes GetMatioClassFromT<double>()
{
	return (MAT_C_DOUBLE);
}








class VclMatio
{
public:
	VclMatio();
	VclMatio(const std::string& filename, bool forReading = true);
	~VclMatio();

	// Reading//////////////////////////////////////////////////////////////////////////////////////////////////////
	bool openForReading(const std::string& filename);

	template <typename T>	bool readValueNamed(const std::string& name, T& value);
	template <typename T>	bool readEigenVectorNamed(const std::string& name, Eigen::Matrix<T, Eigen::Dynamic, 1>& value);
	template <typename T>	bool readEigenVectorNamed(const std::string& name, Eigen::Matrix<T, 1, Eigen::Dynamic>& value);
	template <typename T>	bool readEigenVector2dNamed(const std::string& name, std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& value);
	template <typename T>	bool readEigenVector2dNamed(const std::string& name, std::vector<Eigen::Matrix<T, 1, Eigen::Dynamic> >& value);
	template <typename T>	bool readEigenMatrix2dNamed(const std::string& name, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& value);
	template <typename T>	bool readEigenMatrix3dNamed(const std::string& name, std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& value);
	template <typename T, size_t N> bool readEigenMatrixndNamed(const std::string& name, boost::multi_array<T, N>& value);
	template <typename T, size_t N> bool readEigenMatrixndNamed(const std::string& name, std::vector<boost::multi_array<T, N> >& value);
	

	// Writing//////////////////////////////////////////////////////////////////////////////////////////////////////
	bool openForWriting(const std::string& filename);

	template <typename T>	bool writeValueNamed(const std::string& name, const T& value);
	template <typename T>	bool writeEigenVectorNamed(const std::string& name, const Eigen::Matrix<T, Eigen::Dynamic, 1>& value);
	template <typename T>	bool writeEigenVectorNamed(const std::string& name, const Eigen::Matrix<T, 1, Eigen::Dynamic>& value);
	template <typename T>	bool writeEigenVector2dNamed(const std::string& name, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& value);
	template <typename T>	bool writeEigenVector2dNamed(const std::string& name, const std::vector<Eigen::Matrix<T, 1, Eigen::Dynamic> >& value);
	template <typename T>	bool writeEigenMatrix2dNamed(const std::string& name, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& value);
	template <typename T>	bool writeEigenMatrix3dNamed(const std::string& name, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& value);
	template <typename T, size_t N>	bool writeEigenMatrixndNamed(const std::string& name, const boost::multi_array<T, N>& value);
	template <typename T, size_t N>	bool writeEigenMatrixndNamed(const std::string& name, const std::vector<boost::multi_array<T, N> >& value);


	// Info
	std::vector<std::string> variableNames();
	bool readVariableNamedDim(const std::string& name, std::vector<size_t>& vDim);
	bool isValue(const std::string& name);
	bool isRowVec(const std::string& name);
	bool isColVec(const std::string& name);
	bool isArray2d(const std::string& name);
	bool isArray3d(const std::string& name);
	bool isArraynd(const std::string& name);
	matio_types GetTypeFromClass(matio_classes cls);

	// Closing
	bool close();

private:
	mat_t* m_matfp;
	std::string m_filename;

	bool matFileOpen(const mat_t* mat);

	//Template helper functions

	template <typename T>			void		CopyData(void* src, matio_classes matioClass, T& dst);
	template <typename T>			void		CopyData(void* src, matio_classes matioClass, std::vector<T>& dst);
	template <typename T>			void		CopyData(void* src, matio_classes matioClass, Eigen::Matrix<T, Eigen::Dynamic, 1>& dst);
	template <typename T>			void		CopyData(void* src, matio_classes matioClass, Eigen::Matrix<T, 1, Eigen::Dynamic>& dst);
	template <typename T>			void		CopyData(void* src, matio_classes matioClass, std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& dst);
	template <typename T>			void		CopyData(void* src, matio_classes matioClass, std::vector<Eigen::Matrix<T, 1, Eigen::Dynamic> >& dst);
	template <typename T>			void		CopyData(void* src, matio_classes matioClass, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& dst);
	template <typename T>			void		CopyData(void* src, matio_classes matioClass, std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& dst);
	template <typename T, size_t N>	void		CopyData(void* src, matio_classes matioClass, boost::multi_array<T, N>& dst);
	template <typename T, size_t N>	void		CopyData(void* src, matio_classes matioClass, std::vector<boost::multi_array<T, N> >& dst);
};





template <typename T>
void VclMatio::CopyData(void* src, matio_classes matioClass, T& dst)
{
	switch (matioClass)
	{
	case MAT_C_INT8:
		dst = ((int8_t*)src)[0];
		break;
	case MAT_C_UINT8:
		dst = ((uint8_t*)src)[0];
		break;
	case MAT_C_INT16:
		dst = ((int16_t*)src)[0];
		break;
	case MAT_C_UINT16:
		dst = ((uint16_t*)src)[0];
		break;
	case MAT_C_INT32:
		dst = ((int32_t*)src)[0];
		break;
	case MAT_C_UINT32:
		dst = ((uint32_t*)src)[0];
		break;
	case MAT_C_INT64:
		dst = ((int64_t*)src)[0];
		break;
	case MAT_C_UINT64:
		dst = ((uint64_t*)src)[0];
		break;
	case MAT_C_SINGLE:
		dst = ((float*)src)[0];
		break;
	case MAT_C_DOUBLE:
		dst = ((double*)src)[0];
		break;
	}
}

template <typename T>
void VclMatio::CopyData(void* src, matio_classes matioClass, std::vector<T>& dst)
{
	switch (matioClass)
	{
	case MAT_C_INT8:
		for (size_t i = 0; i < dst.size(); i++)
			dst[i] = ((int8_t*)src)[i];
		break;
	case MAT_C_UINT8:
		for (size_t i = 0; i < dst.size(); i++)
			dst[i] = ((uint8_t*)src)[i];
		break;
	case MAT_C_INT16:
		for (size_t i = 0; i < dst.size(); i++)
			dst[i] = ((int16_t*)src)[i];
		break;
	case MAT_C_UINT16:
		for (size_t i = 0; i < dst.size(); i++)
			dst[i] = ((uint16_t*)src)[i];
		break;
	case MAT_C_INT32:
		for (size_t i = 0; i < dst.size(); i++)
			dst[i] = ((int32_t*)src)[i];
		break;
	case MAT_C_UINT32:
		for (size_t i = 0; i < dst.size(); i++)
			dst[i] = ((uint32_t*)src)[i];
		break;
	case MAT_C_INT64:
		for (size_t i = 0; i < dst.size(); i++)
			dst[i] = ((int64_t*)src)[i];
		break;
	case MAT_C_UINT64:
		for (size_t i = 0; i < dst.size(); i++)
			dst[i] = ((uint64_t*)src)[i];
		break;
	case MAT_C_SINGLE:
		for (size_t i = 0; i < dst.size(); i++)
			dst[i] = ((float*)src)[i];
		break;
	case MAT_C_DOUBLE:
		for (size_t i = 0; i < dst.size(); i++)
			dst[i] = ((double*)src)[i];
		break;
	}
}

template <typename T>
void VclMatio::CopyData(void* src, matio_classes matioClass, Eigen::Matrix<T, Eigen::Dynamic, 1>& dst)
{
	switch (matioClass)
	{
	case MAT_C_INT8:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((int8_t*)src)[i];
		break;
	case MAT_C_UINT8:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((uint8_t*)src)[i];
		break;
	case MAT_C_INT16:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((int16_t*)src)[i];
		break;
	case MAT_C_UINT16:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((uint16_t*)src)[i];
		break;
	case MAT_C_INT32:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((int32_t*)src)[i];
		break;
	case MAT_C_UINT32:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((uint32_t*)src)[i];
		break;
	case MAT_C_INT64:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((int64_t*)src)[i];
		break;
	case MAT_C_UINT64:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((uint64_t*)src)[i];
		break;
	case MAT_C_SINGLE:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((float*)src)[i];
		break;
	case MAT_C_DOUBLE:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((double*)src)[i];
		break;
	}
}

template <typename T>
void VclMatio::CopyData(void* src, matio_classes matioClass, Eigen::Matrix<T, 1, Eigen::Dynamic>& dst)
{
	switch (matioClass)
	{
	case MAT_C_INT8:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((int8_t*)src)[i];
		break;
	case MAT_C_UINT8:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((uint8_t*)src)[i];
		break;
	case MAT_C_INT16:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((int16_t*)src)[i];
		break;
	case MAT_C_UINT16:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((uint16_t*)src)[i];
		break;
	case MAT_C_INT32:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((int32_t*)src)[i];
		break;
	case MAT_C_UINT32:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((uint32_t*)src)[i];
		break;
	case MAT_C_INT64:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((int64_t*)src)[i];
		break;
	case MAT_C_UINT64:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((uint64_t*)src)[i];
		break;
	case MAT_C_SINGLE:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((float*)src)[i];
		break;
	case MAT_C_DOUBLE:
		for (size_t i = 0; i < dst.size(); i++)
			dst(i) = ((double*)src)[i];
		break;
	}
}


template <typename T>
void VclMatio::CopyData(void* src, matio_classes matioClass, std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& dst)
{
	switch (matioClass)
	{
	case MAT_C_INT8:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((int8_t*)src)[j + i*dst[i].size()];
		break;
	case MAT_C_UINT8:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((uint8_t*)src)[j + i*dst[i].size()];
		break;
	case MAT_C_INT16:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((int16_t*)src)[j + i*dst[i].size()];
		break;
	case MAT_C_UINT16:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((uint16_t*)src)[j + i*dst[i].size()];
		break;
	case MAT_C_INT32:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((int32_t*)src)[j + i*dst[i].size()];
		break;
	case MAT_C_UINT32:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((uint32_t*)src)[j + i*dst[i].size()];
		break;
	case MAT_C_INT64:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((int64_t*)src)[j + i*dst[i].size()];
		break;
	case MAT_C_UINT64:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((uint64_t*)src)[j + i*dst[i].size()];
		break;
	case MAT_C_SINGLE:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((float*)src)[j + i*dst[i].size()];
		break;
	case MAT_C_DOUBLE:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((double*)src)[j + i*dst[i].size()];
		break;
	}
}

template <typename T>
void VclMatio::CopyData(void* src, matio_classes matioClass, std::vector<Eigen::Matrix<T, 1, Eigen::Dynamic> >& dst)
{
	switch (matioClass)
	{
	case MAT_C_INT8:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((int8_t*)src)[i + j*dst.size()];
		break;
	case MAT_C_UINT8:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((uint8_t*)src)[i + j*dst.size()];
		break;
	case MAT_C_INT16:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((int16_t*)src)[i + j*dst.size()];
		break;
	case MAT_C_UINT16:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((uint16_t*)src)[i + j*dst.size()];
		break;
	case MAT_C_INT32:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((int32_t*)src)[i + j*dst.size()];
		break;
	case MAT_C_UINT32:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((uint32_t*)src)[i + j*dst.size()];
		break;
	case MAT_C_INT64:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((int64_t*)src)[i + j*dst.size()];
		break;
	case MAT_C_UINT64:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((uint64_t*)src)[i + j*dst.size()];
		break;
	case MAT_C_SINGLE:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((float*)src)[i + j*dst.size()];
		break;
	case MAT_C_DOUBLE:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].size(); j++)
				dst[i](j) = ((double*)src)[i + j*dst.size()];
		break;
	}
}

template <typename T>
void VclMatio::CopyData(void* src, matio_classes matioClass, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& dst)
{
	switch (matioClass)
	{
	case MAT_C_INT8:
		for (size_t i = 0; i < dst.cols(); i++)
			for (size_t j = 0; j < dst.rows(); j++)
				dst(j, i) = ((int8_t*)src)[j + i*dst.rows()];
		break;
	case MAT_C_UINT8:
		for (size_t i = 0; i < dst.cols(); i++)
			for (size_t j = 0; j < dst.rows(); j++)
				dst(j, i) = ((uint8_t*)src)[j + i*dst.rows()];
		break;
	case MAT_C_INT16:
		for (size_t i = 0; i < dst.cols(); i++)
			for (size_t j = 0; j < dst.rows(); j++)
				dst(j, i) = ((int16_t*)src)[j + i*dst.rows()];
		break;
	case MAT_C_UINT16:
		for (size_t i = 0; i < dst.cols(); i++)
			for (size_t j = 0; j < dst.rows(); j++)
				dst(j, i) = ((uint16_t*)src)[j + i*dst.rows()];
		break;
	case MAT_C_INT32:
		for (size_t i = 0; i < dst.cols(); i++)
			for (size_t j = 0; j < dst.rows(); j++)
				dst(j, i) = ((int32_t*)src)[j + i*dst.rows()];
		break;
	case MAT_C_UINT32:
		for (size_t i = 0; i < dst.cols(); i++)
			for (size_t j = 0; j < dst.rows(); j++)
				dst(j, i) = ((uint32_t*)src)[j + i*dst.rows()];
		break;
	case MAT_C_INT64:
		for (size_t i = 0; i < dst.cols(); i++)
			for (size_t j = 0; j < dst.rows(); j++)
				dst(j, i) = ((int64_t*)src)[j + i*dst.rows()];
		break;
	case MAT_C_UINT64:
		for (size_t i = 0; i < dst.cols(); i++)
			for (size_t j = 0; j < dst.rows(); j++)
				dst(j, i) = ((uint64_t*)src)[j + i*dst.rows()];
		break;
	case MAT_C_SINGLE:
		for (size_t i = 0; i < dst.cols(); i++)
			for (size_t j = 0; j < dst.rows(); j++)
				dst(j, i) = ((float*)src)[j + i*dst.rows()];
		break;
	case MAT_C_DOUBLE:
		for (size_t i = 0; i < dst.cols(); i++)
			for (size_t j = 0; j < dst.rows(); j++)
				dst(j, i) = ((double*)src)[j + i*dst.rows()];
		break;
	}
}

template <typename T>
void VclMatio::CopyData(void* src, matio_classes matioClass, std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& dst)
{
	switch (matioClass)
	{
	case MAT_C_INT8:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].cols(); j++)
				for (size_t k = 0; k < dst[i].rows(); k++)
					dst[i](k, j) = ((int8_t*)src)[k + j*dst[i].rows() + i*dst[i].rows()*dst[i].cols()];
		break;
	case MAT_C_UINT8:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].cols(); j++)
				for (size_t k = 0; k < dst[i].rows(); k++)
					dst[i](k, j) = ((uint8_t*)src)[k + j*dst[i].rows() + i*dst[i].rows()*dst[i].cols()];
		break;
	case MAT_C_INT16:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].cols(); j++)
				for (size_t k = 0; k < dst[i].rows(); k++)
					dst[i](k, j) = ((int16_t*)src)[k + j*dst[i].rows() + i*dst[i].rows()*dst[i].cols()];
		break;
	case MAT_C_UINT16:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].cols(); j++)
				for (size_t k = 0; k < dst[i].rows(); k++)
					dst[i](k, j) = ((uint16_t*)src)[k + j*dst[i].rows() + i*dst[i].rows()*dst[i].cols()];
		break;
	case MAT_C_INT32:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].cols(); j++)
				for (size_t k = 0; k < dst[i].rows(); k++)
					dst[i](k, j) = ((int32_t*)src)[k + j*dst[i].rows() + i*dst[i].rows()*dst[i].cols()];
		break;
	case MAT_C_UINT32:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].cols(); j++)
				for (size_t k = 0; k < dst[i].rows(); k++)
					dst[i](k, j) = ((uint32_t*)src)[k + j*dst[i].rows() + i*dst[i].rows()*dst[i].cols()];
		break;
	case MAT_C_INT64:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].cols(); j++)
				for (size_t k = 0; k < dst[i].rows(); k++)
					dst[i](k, j) = ((int64_t*)src)[k + j*dst[i].rows() + i*dst[i].rows()*dst[i].cols()];
		break;
	case MAT_C_UINT64:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].cols(); j++)
				for (size_t k = 0; k < dst[i].rows(); k++)
					dst[i](k, j) = ((uint64_t*)src)[k + j*dst[i].rows() + i*dst[i].rows()*dst[i].cols()];
		break;
	case MAT_C_SINGLE:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].cols(); j++)
				for (size_t k = 0; k < dst[i].rows(); k++)
					dst[i](k, j) = ((float*)src)[k + j*dst[i].rows() + i*dst[i].rows()*dst[i].cols()];
		break;
	case MAT_C_DOUBLE:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < dst[i].cols(); j++)
				for (size_t k = 0; k < dst[i].rows(); k++)
					dst[i](k, j) = ((double*)src)[k + j*dst[i].rows() + i*dst[i].rows()*dst[i].cols()];
		break;
	}
}

template <typename T, size_t N>
void VclMatio::CopyData(void* src, matio_classes matioClass, boost::multi_array<T, N>& dst)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;
	size_t iLength = std::accumulate(dst.shape(), dst.shape() + dst.num_dimensions(), 1, std::multiplies<szType>());

	switch (matioClass)
	{
	case MAT_C_INT8:
		for (size_t i = 0; i < iLength; i++)
			dst.data()[i] = ((int8_t*)src)[i];
		break;
	case MAT_C_UINT8:
		for (size_t i = 0; i < iLength; i++)
			dst.data()[i] = ((uint8_t*)src)[i];
		break;
	case MAT_C_INT16:
		for (size_t i = 0; i < iLength; i++)
			dst.data()[i] = ((int16_t*)src)[i];
		break;
	case MAT_C_UINT16:
		for (size_t i = 0; i < iLength; i++)
			dst.data()[i] = ((uint16_t*)src)[i];
		break;
	case MAT_C_INT32:
		for (size_t i = 0; i < iLength; i++)
			dst.data()[i] = ((int32_t*)src)[i];
		break;
	case MAT_C_UINT32:
		for (size_t i = 0; i < iLength; i++)
			dst.data()[i] = ((uint32_t*)src)[i];
		break;
	case MAT_C_INT64:
		for (size_t i = 0; i < iLength; i++)
			dst.data()[i] = ((int64_t*)src)[i];
		break;
	case MAT_C_UINT64:
		for (size_t i = 0; i < iLength; i++)
			dst.data()[i] = ((uint64_t*)src)[i];
		break;
	case MAT_C_SINGLE:
		for (size_t i = 0; i < iLength; i++)
			dst.data()[i] = ((float*)src)[i];
		break;
	case MAT_C_DOUBLE:
		for (size_t i = 0; i < iLength; i++)
			dst.data()[i] = ((double*)src)[i];
		break;
	}
}


template <typename T, size_t N>
void VclMatio::CopyData(void* src, matio_classes matioClass, std::vector<boost::multi_array<T, N> >& dst)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;
	size_t iLength = std::accumulate(dst[0].shape(), dst[0].shape() + dst[0].num_dimensions(), 1, std::multiplies<szType>());

	switch (matioClass)
	{
	case MAT_C_INT8:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < iLength; j++)
				dst[i].data()[j] = ((int8_t*)src)[j + i*iLength];
		break;
	case MAT_C_UINT8:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < iLength; j++)
				dst[i].data()[j] = ((uint8_t*)src)[j + i*iLength];
		break;
	case MAT_C_INT16:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < iLength; j++)
				dst[i].data()[j] = ((int16_t*)src)[j + i*iLength];
		break;
	case MAT_C_UINT16:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < iLength; j++)
				dst[i].data()[j] = ((uint16_t*)src)[j + i*iLength];
		break;
	case MAT_C_INT32:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < iLength; j++)
				dst[i].data()[j] = ((int32_t*)src)[j + i*iLength];
		break;
	case MAT_C_UINT32:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < iLength; j++)
				dst[i].data()[j] = ((uint32_t*)src)[j + i*iLength];
		break;
	case MAT_C_INT64:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < iLength; j++)
				dst[i].data()[j] = ((int64_t*)src)[j + i*iLength];
		break;
	case MAT_C_UINT64:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < iLength; j++)
				dst[i].data()[j] = ((uint64_t*)src)[j + i*iLength];
		break;
	case MAT_C_SINGLE:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < iLength; j++)
				dst[i].data()[j] = ((float*)src)[j + i*iLength];
		break;
	case MAT_C_DOUBLE:
		for (size_t i = 0; i < dst.size(); i++)
			for (size_t j = 0; j < iLength; j++)
				dst[i].data()[j] = ((double*)src)[j + i*iLength];
		break;
	}
}




template <typename T>
bool VclMatio::readValueNamed(const std::string& name, T& value)
{
	value = 0.0f;
	if (!matFileOpen(m_matfp))
		return false;

	Mat_Rewind(m_matfp);
	matvar_t* matvarRead = Mat_VarRead(m_matfp, name.c_str());
	if (!matvarRead)
		return false;
	CopyData<T>(matvarRead->data, matvarRead->class_type, value);

	Mat_VarFree(matvarRead);
	return true;
}


template <typename T>
bool VclMatio::readEigenVectorNamed(const std::string& name, Eigen::Matrix<T, Eigen::Dynamic, 1>& value)
{
	value.resize(0);
	if (!matFileOpen(m_matfp))
		return false;

	// Read the data into memory
	Mat_Rewind(m_matfp);
	matvar_t* matvarRead = Mat_VarRead(m_matfp, name.c_str());
	if (!matvarRead)
		return false;

	size_t dim = matvarRead->dims[0] * matvarRead->dims[1];

	// Assign the actual value to our return variable
	value.resize(dim);
	CopyData<T>(matvarRead->data, matvarRead->class_type, value);

	Mat_VarFree(matvarRead);
	return true;
}


template <typename T>
bool VclMatio::readEigenVectorNamed(const std::string& name, Eigen::Matrix<T, 1, Eigen::Dynamic>& value)
{
	value.resize(0);
	if (!matFileOpen(m_matfp))
		return false;

	// Read the data into memory
	Mat_Rewind(m_matfp);
	matvar_t* matvarRead = Mat_VarRead(m_matfp, name.c_str());
	if (!matvarRead)
		return false;

	size_t dim = matvarRead->dims[0] * matvarRead->dims[1];

	// Assign the actual value to our return variable
	value.resize(dim);
	CopyData<T>(matvarRead->data, matvarRead->class_type, value);

	Mat_VarFree(matvarRead);
	return true;
}


template <typename T>
bool VclMatio::readEigenVector2dNamed(const std::string& name, std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& value)
{
	value.clear();
	if (!matFileOpen(m_matfp))
		return false;

	// Make sure the variable exists, is floating point, and is a 2d matrix
	if (!isArray2d(name))
		return false;

	// Read the data into memory
	Mat_Rewind(m_matfp);
	matvar_t* matvarRead = Mat_VarRead(m_matfp, name.c_str());

	const size_t rows = matvarRead->dims[0];
	const size_t cols = matvarRead->dims[1];
	// Allocate space
	value.resize(cols);
	for (size_t i = 0; i < cols; ++i)
		value[i].resize(rows);

	// Assign the actual value to our return variable
	CopyData<T>(matvarRead->data, matvarRead->class_type, value);

	Mat_VarFree(matvarRead);
	return true;
}


template <typename T>
bool VclMatio::readEigenVector2dNamed(const std::string& name, std::vector<Eigen::Matrix<T, 1, Eigen::Dynamic> >& value)
{
	value.clear();
	if (!matFileOpen(m_matfp))
		return false;

	// Make sure the variable exists, is floating point, and is a 2d matrix
	if (!isArray2d(name))
		return false;

	// Read the data into memory
	Mat_Rewind(m_matfp);
	matvar_t* matvarRead = Mat_VarRead(m_matfp, name.c_str());

	const size_t rows = matvarRead->dims[0];
	const size_t cols = matvarRead->dims[1];
	// Allocate space
	value.resize(rows);
	for (size_t i = 0; i < rows; ++i)
		value[i].resize(cols);

	// Assign the actual value to our return variable
	CopyData<T>(matvarRead->data, matvarRead->class_type, value);

	Mat_VarFree(matvarRead);
	return true;
}


template <typename T>
bool VclMatio::readEigenMatrix2dNamed(const std::string& name, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& value)
{
	value.resize(0, 0);
	if (!matFileOpen(m_matfp))
		return false;

	// Make sure the variable exists, and is a 2d matrix
	if (!isArray2d(name))
		return false;

	// Read the data into memory
	Mat_Rewind(m_matfp);
	matvar_t* matvarRead = Mat_VarRead(m_matfp, name.c_str());

	// Assign the actual value to our return variable
	const size_t rows = matvarRead->dims[0];
	const size_t cols = matvarRead->dims[1];
	value.resize(rows, cols);

	// Assign the actual value to our return variable
	CopyData<T>(matvarRead->data, matvarRead->class_type, value);

	Mat_VarFree(matvarRead);
	return true;
}


template <typename T>
bool VclMatio::readEigenMatrix3dNamed(const std::string& name, std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& value)
{
	value.clear();
	if (!matFileOpen(m_matfp))
		return false;

	// Make sure the variable existsand is a 3d matrix
	if (!isArray3d(name))
		return false;

	// Read the data into memory
	Mat_Rewind(m_matfp);
	matvar_t* matvarRead = Mat_VarRead(m_matfp, name.c_str());

	const size_t dimX = matvarRead->dims[0];
	const size_t dimY = matvarRead->dims[1];
	const size_t dimZ = matvarRead->dims[2];

	// Allocate space
	value.resize(dimZ);
	for (size_t i = 0; i < dimZ; ++i)
		value[i].resize(dimX, dimY);

	// Assign the actual value to our return variable
	CopyData<T>(matvarRead->data, matvarRead->class_type, value);

	Mat_VarFree(matvarRead);
	return true;
}


template <typename T, size_t N>
bool VclMatio::readEigenMatrixndNamed(const std::string& name, boost::multi_array<T, N>& value)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	if (!matFileOpen(m_matfp))
		return false;

	// Make sure the variable exists and is an nD array
	if (!isArraynd(name))
		return false;

	// Read the data into memory
	Mat_Rewind(m_matfp);
	matvar_t* matvarRead = Mat_VarRead(m_matfp, name.c_str());

	// Allocate space
	std::vector<szType> vShape;
	vShape.assign(matvarRead->dims, matvarRead->dims + matvarRead->rank);
	szType iLength = std::accumulate(vShape.begin(), vShape.end(), 1, std::multiplies<szType>());
	value.resize(vShape);

	// Assign the actual value to our return variable
	CopyData<T>(matvarRead->data, matvarRead->class_type, value);

	Mat_VarFree(matvarRead);
	return true;
}


template <typename T, size_t N> 
bool VclMatio::readEigenMatrixndNamed(const std::string& name, std::vector<boost::multi_array<T, N> >& value)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	value.clear();
	if (!matFileOpen(m_matfp))
		return false;

	// Make sure the variable exists and is an nD array
	if (!isArraynd(name))
		return false;

	// Read the data into memory
	Mat_Rewind(m_matfp);
	matvar_t* matvarRead = Mat_VarRead(m_matfp, name.c_str());

	// Allocate space
	std::vector<szType> vShape;
	vShape.assign(matvarRead->dims, matvarRead->dims + matvarRead->rank);
	size_t iN = vShape.back();
	vShape.pop_back();
	value.resize(iN);
	for (size_t i = 0; i < iN; ++i)
		value[i].resize(vShape);

	// Assign the actual value to our return variable
	CopyData<T>(matvarRead->data, matvarRead->class_type, value);

	Mat_VarFree(matvarRead);
	return true;
}




template <typename T>
bool VclMatio::writeValueNamed(const std::string& name, const T& value)
{
	matvar_t* matvar;
	size_t dims[2] = { 1, 1 };
	matio_classes cls = GetMatioClassFromT<T>();
	matvar = Mat_VarCreate(name.c_str(), cls, GetTypeFromClass(cls), 2, dims, (void*)&value, MAT_F_DONT_COPY_DATA);
	if (matvar == NULL)
		return false;
	Mat_VarWrite(m_matfp, matvar, MAT_COMPRESSION_NONE);
	Mat_VarFree(matvar);
	return true;
}


template <typename T>
bool VclMatio::writeEigenVectorNamed(const std::string& name, const Eigen::Matrix<T, Eigen::Dynamic, 1>& value)
{
	matvar_t* matvar;
	size_t dims[2] = { value.size(), 1 };
	matio_classes cls = GetMatioClassFromT<T>();
	matvar = Mat_VarCreate(name.c_str(), cls, GetTypeFromClass(cls), 2, dims, (void*)value.data(), MAT_F_DONT_COPY_DATA);
	if (matvar == NULL)
		return false;
	Mat_VarWrite(m_matfp, matvar, MAT_COMPRESSION_NONE);
	Mat_VarFree(matvar);
	return true;
}


template <typename T>
bool VclMatio::writeEigenVectorNamed(const std::string& name, const Eigen::Matrix<T, 1, Eigen::Dynamic>& value)
{
	matvar_t* matvar;
	size_t dims[2] = { 1, value.size() };
	matio_classes cls = GetMatioClassFromT<T>();
	matvar = Mat_VarCreate(name.c_str(), cls, GetTypeFromClass(cls), 2, dims, (void*)value.data(), MAT_F_DONT_COPY_DATA);
	if (matvar == NULL)
		return false;
	Mat_VarWrite(m_matfp, matvar, MAT_COMPRESSION_NONE);
	Mat_VarFree(matvar);
	return true;
}


template <typename T>
bool VclMatio::writeEigenVector2dNamed(const std::string& name, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& value)
{
	if (value.empty())
		return false;

	matvar_t* matvar;
	// Copy data to matrix
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp;
	tmp.resize(value[0].size(), value.size());
	for (size_t i = 0; i < tmp.rows(); ++i)
		for (size_t j = 0; j < tmp.cols(); ++j)
			tmp(i, j) = (value[j])(i);
	size_t dims[2] = { tmp.rows(), tmp.cols() };
	matio_classes cls = GetMatioClassFromT<T>();
	matvar = Mat_VarCreate(name.c_str(), cls, GetTypeFromClass(cls), 2, dims, (void*)tmp.data(), MAT_F_DONT_COPY_DATA);
	if (matvar == NULL)
		return false;
	Mat_VarWrite(m_matfp, matvar, MAT_COMPRESSION_NONE);
	Mat_VarFree(matvar);
	return true;
}


template <typename T>
bool VclMatio::writeEigenVector2dNamed(const std::string& name, const std::vector<Eigen::Matrix<T, 1, Eigen::Dynamic> >& value)
{
	if (value.empty())
		return false;

	matvar_t* matvar;
	// Copy data to matrix
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp;
	tmp.resize(value.size(), value[0].size());
	for (size_t i = 0; i < tmp.rows(); ++i)
		for (size_t j = 0; j < tmp.cols(); ++j)
			tmp(i, j) = (value[i])(j);
	size_t dims[2] = { tmp.rows(), tmp.cols() };
	matio_classes cls = GetMatioClassFromT<T>();
	matvar = Mat_VarCreate(name.c_str(), cls, GetTypeFromClass(cls), 2, dims, (void*)tmp.data(), MAT_F_DONT_COPY_DATA);
	if (matvar == NULL)
		return false;
	Mat_VarWrite(m_matfp, matvar, MAT_COMPRESSION_NONE);
	Mat_VarFree(matvar);
	return true;
}


template <typename T>
bool VclMatio::writeEigenMatrix2dNamed(const std::string& name, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& value)
{
	matvar_t* matvar;
	size_t dims[2] = { value.rows(), value.cols() };
	matio_classes cls = GetMatioClassFromT<T>();
	matvar = Mat_VarCreate(name.c_str(), cls, GetTypeFromClass(cls), 2, dims, (void*)value.data(), MAT_F_DONT_COPY_DATA);
	if (matvar == NULL)
		return false;
	Mat_VarWrite(m_matfp, matvar, MAT_COMPRESSION_NONE);
	Mat_VarFree(matvar);
	return true;
}


template <typename T>
bool VclMatio::writeEigenMatrix3dNamed(const std::string& name, const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& value)
{
	if (value.empty())
		return false;

	matvar_t* matvar;
	size_t dims[3] = { value[0].rows(), value[0].cols(), value.size() };

	// Makes a copy of the data before writing it to the file. 
	std::vector<T> allData;
	size_t iCounter = 0;
	for (size_t i = 0; i < value.size(); ++i)
		iCounter += value[i].size();
	allData.resize(iCounter);
	for (size_t i = 0; i < value.size(); i++)
		for (size_t j = 0; j < value[i].cols(); j++)
			for (size_t k = 0; k < value[i].rows(); k++)
				allData[k + j*value[i].rows() + i*value[i].size()] = value[i](k, j);
	matio_classes cls = GetMatioClassFromT<T>();
	matvar = Mat_VarCreate(name.c_str(), cls, GetTypeFromClass(cls), 3, dims, (void*)allData.data(), 0);
	if (matvar == NULL)
		return false;
	Mat_VarWrite(m_matfp, matvar, MAT_COMPRESSION_NONE);
	Mat_VarFree(matvar);
	return true;

	/*
	// Attempts to write the data without copying.  Doesn't work yet...
	// http://groups.inf.ed.ac.uk/vision/MAJECKA/Detector/Resources/libraries/matio-1.3.3/test/test_mat.c
	matvar = Mat_VarCreate(name.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, 3, dims, NULL, 0);
	Mat_VarWriteInfo(m_matfp, matvar);
	for (int i = 0; i < value.size(); i++)
	{
	int start[2] = {0,0};
	int stride[2] = {1,1};
	int edge[2] = {3,5};
	Mat_VarWriteData(m_matfp, matvar, (void*)value[i].data(), start, stride, edge);
	}
	*/
}


template <typename T, size_t N>
bool VclMatio::writeEigenMatrixndNamed(const std::string& name, const boost::multi_array<T, N>& value)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;

	matvar_t* matvar;
	std::vector<size_t> vDims;
	vDims.assign(value.shape(), value.shape() + value.num_dimensions());

	matio_classes cls = GetMatioClassFromT<T>();
	matvar = Mat_VarCreate(name.c_str(), cls, GetTypeFromClass(cls), N, vDims.data(), (void*)value.data(), 0);
	if (matvar == NULL)
		return false;
	Mat_VarWrite(m_matfp, matvar, MAT_COMPRESSION_NONE);
	Mat_VarFree(matvar);
	return true;
}


template <typename T, size_t N>
bool VclMatio::writeEigenMatrixndNamed(const std::string& name, const std::vector<boost::multi_array<T, N> >& value)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename boost::multi_array<T, N>::index idxType;
	if (value.empty())
		return false;
	size_t iN = value.size();
	for (size_t i = 0; i < iN; ++i)
		if (!std::equal(value[0].shape(), value[0].shape() + value[0].num_dimensions(), value[i].shape()))
			return false;

	matvar_t* matvar;
	std::vector<size_t> vDims;
	vDims.assign(value[0].shape(), value[0].shape() + value[0].num_dimensions());
	szType iLength = std::accumulate(vDims.begin(), vDims.end(), 1, std::multiplies<size_t>());
	vDims.push_back(iN);

	std::vector<T> allData(iN * iLength);
	for (size_t i = 0; i < iN; ++i)
		for (size_t j = 0; j < iLength; ++j)
			allData[j + i*iLength] = value[i].data()[j];

	matio_classes cls = GetMatioClassFromT<T>();
	matvar = Mat_VarCreate(name.c_str(), cls, GetTypeFromClass(cls), N + 1, vDims.data(), (void*)allData.data(), 0);
	if (matvar == NULL)
		return false;
	Mat_VarWrite(m_matfp, matvar, MAT_COMPRESSION_NONE);
	Mat_VarFree(matvar);
	return true;
}