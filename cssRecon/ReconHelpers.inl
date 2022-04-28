



template <typename T>
void SL0(Eigen::Matrix<T, Eigen::Dynamic, 1>& s,
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
	T sigma_min,
	T sigma_decrease_factor,
	T mu_0,
	size_t L)
{
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A_pinv = A.transpose() * (A * A.transpose()).inverse();

// 	Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
// 	double tolerance = std::numeric_limits<double>::epsilon() * std::max(A.rows(), A.cols()) * fabs(svd.singularValues()(0));
// 	A_pinv = svd.matrixV() * (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().transpose();
	
	s = A_pinv * x;
	T sigma = s.cwiseAbs().maxCoeff();
	Eigen::Matrix<T, Eigen::Dynamic, 1> delta(s.size());

	while (sigma > sigma_min)
	{
		for (size_t i = 0; i < L; ++i)
		{
			delta = (s.array() * (-1.0*s.array().pow(2.0) / (2.0 * sigma*sigma)).exp()).matrix();
			s = s - mu_0*delta;
			s = s - A_pinv*(A*s - x);
		}
		sigma = sigma * sigma_decrease_factor;
	}
}

template <typename T>
void SL0(Eigen::Matrix<T, Eigen::Dynamic, 1>& s,
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A_pinv,
	T sigma_min,
	T sigma_decrease_factor,
	T mu_0,
	size_t L)
{
	s = A_pinv * x;
	T sigma = s.cwiseAbs().maxCoeff();
	Eigen::Matrix<T, Eigen::Dynamic, 1> delta(s.size());

	while (sigma > sigma_min)
	{
		for (size_t i = 0; i < L; ++i)
		{
			delta = (s.array() * (-1.0*s.array().pow(2.0) / (2.0 * sigma*sigma)).exp()).matrix();
			s = s - mu_0*delta;
			s = s - A_pinv*(A*s - x);
		}
		sigma = sigma * sigma_decrease_factor;
	}
}



template <typename T>
void OMP(Eigen::Matrix<T, Eigen::Dynamic, 1>& s,
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
	size_t iSparsity)
{
	s.setZero(A.cols());
	Eigen::Matrix<T, Eigen::Dynamic, 1> r(x.size());
	r = x;

	std::vector<size_t> vLambda;
	vLambda.reserve(iSparsity);
	Eigen::Matrix<T, Eigen::Dynamic, 1> s_nz;

	for (size_t k = 0; k < iSparsity; ++k)
	{
		size_t lambda_k;
		(A.transpose() * r).array().abs().maxCoeff(&lambda_k);
		vLambda.push_back(lambda_k);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A_lambda(A.rows(), vLambda.size());
		for (size_t i = 0; i < vLambda.size(); ++i)
			A_lambda.col(i) = A.col(vLambda[i]);
//		s_nz = A_lambda.householderQr().solve(x);
		Eigen::JacobiSVD<typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > svd(A_lambda, Eigen::ComputeThinU | Eigen::ComputeThinV);
		s_nz = svd.solve(x);
		r = x - A_lambda * s_nz;
	}
	for (size_t k = 0; k < iSparsity; ++k)
		s(vLambda[k]) = s_nz(k);
}


//This implementation is based on 
//http://www.cs.technion.ac.il/users/wwwb/cgi-bin/tr-get.cgi/2008/CS/CS-2008-08.pdf
//and consumes more memory but it is much faster
template <typename T>
void BatchOMP(Eigen::Matrix<T, Eigen::Dynamic, 1>& s,
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vA,
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vAtA,
	const Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
	size_t iSparsity)
{
	s.setZero(vA.cols());

	std::vector<size_t> vLambda;
	vLambda.reserve(iSparsity);
	Eigen::Matrix<T, Eigen::Dynamic, 1> vAlpha0 = vA.transpose() * x;
	Eigen::Matrix<T, Eigen::Dynamic, 1> vAlpha = vAlpha0;
	Eigen::Matrix<T, Eigen::Dynamic, 1> vAlpha0_lambda(iSparsity);
	vAlpha0_lambda.setZero();

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> AtA_lambda(vA.cols(), iSparsity);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> L(iSparsity, iSparsity);
	L.setZero();
	Eigen::Matrix<T, Eigen::Dynamic, 1> s_nz;

	size_t i = 0; 

	while (i < iSparsity)
	{
		//Atom selection
		size_t lambda_k;
		vAlpha.array().abs().maxCoeff(&lambda_k);
		vLambda.push_back(lambda_k);

		//Update intermediate variables based on selected atom
		AtA_lambda.col(i) = vAtA.col(vLambda[i]);
		vAlpha0_lambda(i) = vAlpha0(vLambda[i]);

		//Incremental Cholesky update stage
		if (i == 0)
		{
			L(0, 0) = 1.0;
		}
		else
		{
			Eigen::Matrix<T, Eigen::Dynamic, 1> nu(vLambda.size());
			Eigen::Matrix<T, Eigen::Dynamic, 1> omega;
			for (size_t j = 0; j < vLambda.size(); ++j)
				nu(j) = AtA_lambda(vLambda[j], i);
			backsubst('L', L, nu, omega, i);
			for (size_t j = 0; j < omega.size(); ++j)
				L(i, j) = omega(j);
			T omega_nrm = omega.squaredNorm();
			if ((1.0 - omega_nrm) < 1e-14)
				break;
			L(i, i) = sqrt(1.0 - omega_nrm);
		}

		i++;

		//Residual update
		Eigen::Matrix<T, Eigen::Dynamic, 1> vTmp = vAlpha0_lambda.segment(0, vLambda.size());
		cholsolve('L', L, vTmp, s_nz, i);
		vAlpha = vAlpha0 - AtA_lambda.block(0, 0, AtA_lambda.rows(), i) * s_nz;
	}

	//Write nonzero coefficients
	for (size_t k = 0; k < i; ++k)
		s(vLambda[k]) = s_nz(k);
}

