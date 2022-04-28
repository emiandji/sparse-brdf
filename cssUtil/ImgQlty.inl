


// 
// template <typename T>
// double MSE(const T* pRef, const T* pRecon, size_t iLength, int64_t iQuantization)
// {
// 	if (!pRef || !pRecon || (iLength == 0) || (iQuantization == 0))
// 		return -1.0;
// 
// 	double tSum = 0;
// 	if (iQuantization == -1)
// 	{
// 		for (size_t i = 0; i < iLength; ++i)
// 			tSum += powf((pRef[i] - pRecon[i]),2.0);
// 		return tSum / iLength;
// 	}
// 	else
// 	{
// 		T tQuantization = iQuantization;
// 		for (size_t i = 0; i < iLength; ++i)
// 			tSum += powf((pRef[i] - pRecon[i]) * tQuantization, 2.0);
// 		return tSum / iLength;
// 	}
// }
// 
// 
// template <typename T>
// double MSE(const std::vector<const T*>& pRef, const std::vector<const T*>& pRecon, const std::vector<size_t>& vLength, int64_t iQuantization)
// {
// 	if ((vLength.size() != pRef.size()) || (pRef.size() != pRecon.size()) || vLength.empty() || (iQuantization == 0))
// 		return -1.0;
// 	for (size_t i = 0; i < pRef.size(); ++i)
// 		if (!pRef[i] || !pRecon[i] || vLength[i] == 0)
// 			return -1.0;
// 
// 	double tSum = 0;
// 	size_t iTotal = std::accumulate(vLength.begin(), vLength.end(), 0);
// 	if (iQuantization == -1)
// 	{
// #pragma omp parallel for schedule(guided) reduction(+:tSum)
// 		for (size_t i = 0; i < pRef.size(); ++i)
// 		{
// 			double tSumLocal = 0;
// 			for (size_t j = 0; j < vLength[i]; ++j)
// 				tSumLocal += powf((pRef[i][j] - pRecon[i][j]), 2.0);
// 			tSum += tSumLocal;
// 		}
// 		return tSum / iTotal;
// 	}
// 	else
// 	{
// 		T tQuantization = iQuantization;
// #pragma omp parallel for schedule(guided) reduction(+:tSum)
// 		for (size_t i = 0; i < pRef.size(); ++i)
// 		{
// 			T tSumLocal = 0;
// 			for (size_t j = 0; j < vLength[i]; ++j)
// 				tSumLocal += powf((pRef[i][j] - pRecon[i][j]) * tQuantization, 2.0);
// 			tSum += tSumLocal;
// 		}
// 		return tSum / iTotal;
// 	}
// }
// 
// 
// template <typename T>
// double PSNR(const T* pRef, const T* pRecon, size_t iLength, int64_t iQuantization)
// {
// 	if (!pRef || !pRecon || (iLength == 0) || (iQuantization == 0))
// 		return -1.0;
// 
// 	if (iQuantization == -1)
// 	{
// 		T refMax = 0;
// 		T reconMax = 0;
// 		for (size_t i = 0; i < iLength; ++i)
// 		{
// 			if (absT<T>(pRef[i]) > refMax)
// 				refMax = absT<T>(pRef[i]);
// 			if (absT<T>(pRecon[i]) > reconMax)
// 				reconMax = absT<T>(pRecon[i]);
// 		}
// 		T tQuantization = std::max<T>(refMax, reconMax);
// 		return 10.0 * log10f(powf(tQuantization, 2.0) / MSE<T>(pRef, pRecon, iLength, -1));
// 	}
// 	return 10.0 * log10f(powf(iQuantization, 2.0) / MSE<T>(pRef, pRecon, iLength, iQuantization));
// }
// 
// 
// template <typename T>
// double PSNR(const std::vector<const T*>& pRef, const std::vector<const T*>& pRecon, const std::vector<size_t>& vLength, int64_t iQuantization)
// {
// 	if ((vLength.size() != pRef.size()) || (pRef.size() != pRecon.size()) || vLength.empty() || (iQuantization == 0))
// 		return -1.0;
// 	for (size_t i = 0; i < pRef.size(); ++i)
// 		if (!pRef[i] || !pRecon[i] || vLength[i] == 0)
// 			return -1.0;
// 
// 	if (iQuantization == -1)
// 	{
// 		T refMax = 0;
// 		T reconMax = 0;
// 		for (size_t i = 0; i < pRef.size(); ++i)
// 		{
// 			for (size_t j = 0; j < vLength[i]; ++j)
// 			{
// 				if (absT<T>(pRef[i][j]) > refMax)
// 					refMax = absT<T>(pRef[i]);
// 				if (absT<T>(pRecon[i][j]) > reconMax)
// 					reconMax = absT<T>(pRecon[i]);
// 			}
// 		}
// 		T tQuantization = std::max<T>(refMax, reconMax);
// 		return 10.0 * log10f(powf(tQuantization, 2.0) / MSE<T>(mRef, mRecon, -1));
// 	}
// 	return 10.0 * log10f(powf(iQuantization, 2.0) / MSE<T>(mRef, mRecon, iQuantization));
// }
// 
// 
// template <typename T>
// double SNR(const T* pRef, const T* pRecon, size_t iLength)
// {
// 	if (!pRef || !pRecon || (iLength == 0))
// 		return -1.0;
// 
// 	double p_signal = 0.0;
// 	for (size_t i = 0; i < iLength; ++i)
// 		p_signal += powf(pRef[i], 2.0);
// 
// 	double p_noise = 0.0;
// 	for (size_t i = 0; i < iLength; ++i)
// 		p_noise += powf(pRef[i] - pRecon[i], 2.0);
// 
// 	return 10.0 * log10f(p_signal / p_noise);
// }
// 
// 
// template <typename T>
// double SNR(const std::vector<const T*>& pRef, const std::vector<const T*>& pRecon, const std::vector<size_t>& vLength)
// {
// 	if ((vLength.size() != pRef.size()) || (pRef.size() != pRecon.size()) || vLength.empty())
// 		return -1.0;
// 	for (size_t i = 0; i < pRef.size(); ++i)
// 		if (!pRef[i] || !pRecon[i] || vLength[i] == 0)
// 			return -1.0;
// 
// 	double p_signal = 0.0;
// #pragma omp parallel for schedule(guided) reduction(+:p_signal)
// 	for (size_t i = 0; i < pRef.size(); ++i)
// 	{
// 		double p_signal_local = 0.0;
// 		for (size_t j = 0; j < vLength[i]; ++j)
// 			p_signal_local += pow(pRef[i][j], 2.0);
// 		p_signal += p_signal_local;
// 	}
// 
// 	double p_noise = 0.0;
// #pragma omp parallel for schedule(guided) reduction(+:p_noise)
// 	for (size_t i = 0; i < pRef.size(); ++i)
// 	{
// 		double p_noise_local = 0.0;
// 		for (size_t j = 0; j < vLength[i]; ++j)
// 			p_noise_local += pow(pRef[i][j] - pRecon[i][j], 2.0);
// 		p_noise += p_noise_local;
// 	}
// 
// 	return 10.0 * log10f(p_signal / p_noise);
// }




template <typename T>
double MSE(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mRef,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& mRecon)
{
	if (mRef.size() != mRecon.size())
		return -1.0;

	double tSum = (mRef.template cast<double>() - mRecon.template cast<double>()).squaredNorm();
	return tSum / mRef.size();
}

template <typename T>
double PSNR(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mRef,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& mRecon)
{
	if (mRef.size() != mRecon.size())
		return -1.0;

	return 10.0 * log10f(powf(std::numeric_limits<T>::max(), 2.0) / MSE<T>(mRef, mRecon));
}

template <typename T>
double SNR(const Eigen::Matrix<T, Eigen::Dynamic, 1>& mRef,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& mRecon)
{
	if (mRef.size() != mRecon.size())
		return -1.0;

	double p_signal = mRef.template cast<double>().squaredNorm();
	double p_noise = (mRef.template cast<double>() - mRecon.template cast<double>()).squaredNorm();
	return 10.0 * log10f(p_signal / p_noise);
}

template <typename T>
double MSE(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& mRef, 
      const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& mRecon)
{
	if (mRef.size() != mRecon.size())
		return -1.0;
	for (size_t i = 0; i < mRef.size(); ++i)
		if (mRef[i].size() != mRecon[i].size())
			return -1.0;

	double tSum = 0;
#pragma omp parallel for schedule(guided) reduction(+:tSum)
	for (size_t i = 0; i < mRef.size(); ++i)
		tSum = tSum + (mRef[i].template cast<double>() - mRecon[i].template cast<double>()).squaredNorm();

	return tSum / (mRef.size()*mRef[0].rows());
}

template <typename T>
double PSNR(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& mRef, 
       const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& mRecon)
{
	if (mRef.size() != mRecon.size())
		return -1.0;
	for (size_t i = 0; i < mRef.size(); ++i)
		if (mRef[i].size() != mRecon[i].size())
			return -1.0;

	return 10.0 * log10f(powf(std::numeric_limits<T>::max(), 2.0) / MSE<T>(mRef, mRecon));
}

template <typename T>
double SNR(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& mRef,
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1> >& mRecon)
{
	if (mRef.size() != mRecon.size())
		return -1.0;
	for (size_t i = 0; i < mRef.size(); ++i)
		if (mRef[i].size() != mRecon[i].size())
			return -1.0;

	double p_signal = 0.0;
#pragma omp parallel for schedule(guided) reduction(+:p_signal)
	for (size_t i = 0; i < mRef.size(); ++i)
		p_signal = p_signal + mRef[i].template cast<double>().squaredNorm();

	double p_noise = 0.0;
#pragma omp parallel for schedule(guided) reduction(+:p_noise)
	for (size_t i = 0; i < mRef.size(); ++i)
		p_noise = p_noise + (mRef[i].template cast<double>() - mRecon[i].template cast<double>()).squaredNorm();

	return 10.0 * log10f(p_signal / p_noise);
}






template <typename T>
double MSE(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRef, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRecon)
{
	if (mRef.size() != mRecon.size())
		return -1.0;

	double tSum = (mRef.template cast<double>() - mRecon.template cast<double>()).squaredNorm();
	return tSum / (mRef.size());
}

template <typename T>
double PSNR(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRef, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRecon)
{
	if (mRef.size() != mRecon.size())
		return -1.0;

	return 10.0 * log10f(powf(std::numeric_limits<T>::max(), 2.0) / MSE<T>(mRef, mRecon));
}

template <typename T>
double SNR(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRef, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mRecon)
{
	if (mRef.size() != mRecon.size())
		return -1.0;

	double p_signal = mRef.template cast<double>().squaredNorm();
	double p_noise = (mRef.template cast<double>() - mRecon.template cast<double>()).squaredNorm();
	return 10.0 * log10f(p_signal / p_noise);
}

template <typename T>
double MSE(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& mRef, 
      const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& mRecon)
{
	if (mRef.size() != mRecon.size())
		return -1.0;
	for(size_t i = 0; i < mRef.size(); ++i)
		if(mRef[i].size() != mRecon[i].size())
			return -1.0;

	double tSum = 0;
#pragma omp parallel for schedule(guided) reduction(+:tSum)
	for (size_t i = 0; i < mRef.size(); ++i)
		tSum = tSum + (mRef[i].template cast<double>() - mRecon[i].template cast<double>()).squaredNorm();

	return tSum / (mRef.size()*mRef[0].size());
}

template <typename T>
double PSNR(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& mRef, 
       const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& mRecon)
{
	if (mRef.size() != mRecon.size())
		return -1.0;
	for (size_t i = 0; i < mRef.size(); ++i)
		if (mRef[i].size() != mRecon[i].size())
			return -1.0;

	return 10.0 * log10f(powf(std::numeric_limits<T>::max(), 2.0) / MSE<T>(mRef, mRecon));
}

template <typename T>
double SNR(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& mRef,
	const std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& mRecon)
{
	if (mRef.size() != mRecon.size())
		return -1.0;
	for (size_t i = 0; i < mRef.size(); ++i)
		if (mRef[i].size() != mRecon[i].size())
			return -1.0;

	double p_signal = 0.0;
#pragma omp parallel for schedule(guided) reduction(+:p_signal)
	for (size_t i = 0; i < mRef.size(); ++i)
		p_signal = p_signal + mRef[i].template cast<double>().squaredNorm();

	double p_noise = 0.0;
#pragma omp parallel for schedule(guided) reduction(+:p_noise)
	for (size_t i = 0; i < mRef.size(); ++i)
		p_noise = p_noise + (mRef[i].template cast<double>() - mRecon[i].template cast<double>()).squaredNorm();

	return 10.0 * log10f(p_signal / p_noise);
}







template <typename T, size_t N>
double MSE(const boost::multi_array<T, N>& mRef, const boost::multi_array<T, N>& mRecon)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	typedef typename Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
	szType iLengthRef = mRef.num_elements();
	szType iLengthRecon = mRecon.num_elements();
	if (iLengthRef != iLengthRecon)
		return -1;

	Vec x = Vec::Map(mRef.data(), iLengthRef);
	Vec y = Vec::Map(mRecon.data(), iLengthRecon);
	double tSum = (x.template cast<double>() - y.template cast<double>()).squaredNorm();
	return tSum / (mRef.size()*iLengthRef);
}


template <typename T, size_t N>
double PSNR(const boost::multi_array<T, N>& mRef, const boost::multi_array<T, N>& mRecon)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	szType iLengthRef = mRef.num_elements();
	szType iLengthRecon = mRecon.num_elements();
	if (iLengthRef != iLengthRecon)
		return -1;

	return 10.0 * log10f(powf(std::numeric_limits<T>::max(), 2.0) / MSE<T>(mRef, mRecon));
}


template <typename T, size_t N>
double SNR(const boost::multi_array<T, N>& mRef, const boost::multi_array<T, N>& mRecon)
{
	typedef typename boost::multi_array<T, N>::size_type szType;
	szType iLengthRef = mRef.num_elements();
	szType iLengthRecon = mRecon.num_elements();
	typedef typename Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
	if (iLengthRef != iLengthRecon)
		return -1;

	Vec x = Vec::Map(mRef.data(), iLengthRef);
	Vec y = Vec::Map(mRecon.data(), iLengthRecon);
	double p_signal = x.template cast<double>().squaredNorm();
	double p_noise = (x.template cast<double>() - y.template cast<double>()).squaredNorm();
	return 10.0 * log10f(p_signal / p_noise);
}

template <typename T, size_t N>
double MSE(const std::vector<boost::multi_array<T, N> >& mRef, const std::vector<boost::multi_array<T, N> >& mRecon)
{
	if (mRef.size() != mRecon.size())
		return -1;
	for (size_t i = 0; i < mRef.size(); ++i)
		if (mRef[i].num_elements() != mRecon[i].num_elements())
			return -1;

	typedef typename boost::multi_array<T, N>::size_type szType;
	szType iLengthRef = mRef[0].num_elements();
	szType iLengthRecon = mRecon[0].num_elements();
	typedef typename Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;

	double tSum = 0;
#pragma omp parallel for schedule(guided) reduction(+:tSum)
	for (int i = 0; i < mRef.size(); ++i)
	{
		Vec x = Vec::Map(mRef[i].data(), iLengthRef);
		Vec y = Vec::Map(mRecon[i].data(), iLengthRecon);
		tSum = tSum + (x.template cast<double>() - y.template cast<double>()).squaredNorm();
	}

	return tSum / (mRef.size()*iLengthRef);
}


template <typename T, size_t N>
double PSNR(const std::vector<boost::multi_array<T, N> >& mRef, const std::vector<boost::multi_array<T, N> >& mRecon)
{
	if (mRef.size() != mRecon.size())
		return -1;
	for (size_t i = 0; i < mRef.size(); ++i)
		if (mRef[i].num_elements() != mRecon[i].num_elements())
			return -1;

	typedef typename boost::multi_array<T, N>::size_type szType;
	szType iLengthRef = mRef[0].num_elements();
	szType iLengthRecon = mRecon[0].num_elements();

	return 10.0 * log10f(powf(std::numeric_limits<T>::max(), 2.0) / MSE<T>(mRef, mRecon));
}


template <typename T, size_t N>
double SNR(const std::vector<boost::multi_array<T, N> >& mRef, const std::vector<boost::multi_array<T, N> >& mRecon)
{
	if (mRef.size() != mRecon.size())
		return -1;
	for (size_t i = 0; i < mRef.size(); ++i)
		if (mRef[i].num_elements() != mRecon[i].num_elements())
			return -1;

	typedef typename boost::multi_array<T, N>::size_type szType;
	szType iLengthRef = mRef[0].num_elements();
	szType iLengthRecon = mRecon[0].num_elements();
	typedef typename Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;

	double p_signal = 0.0;
#pragma omp parallel for schedule(guided) reduction(+:p_signal)
	for (int i = 0; i < mRef.size(); ++i)
	{
		Vec x = Vec::Map(mRef[i].data(), iLengthRef);
		p_signal = p_signal + (x.template cast<double>()).squaredNorm();
	}

	double p_noise = 0.0;
#pragma omp parallel for schedule(guided) reduction(+:p_noise)
	for (int i = 0; i < mRef.size(); ++i)
	{
		Vec x = Vec::Map(mRef[i].data(), iLengthRef);
		Vec y = Vec::Map(mRecon[i].data(), iLengthRecon);
		p_noise = p_noise + (x.template cast<double>() - y.template cast<double>()).squaredNorm();
	}

	return 10.0 * log10f(p_signal / p_noise);
}

template<typename T>
inline double PSNR(double mse)
{
	return 10.0 * log10f(powf(std::numeric_limits<T>::max(), 2.0) / mse);
}

