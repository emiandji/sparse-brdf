

#pragma once

#include <assert.h>
#include <limits.h>



template <typename T>
class Distribution1D 
{
public:
	// Distribution1D Public Methods
	Distribution1D(const T *f, int n) {
		count = n;
		func = new T[n];
		memcpy(func, f, n * sizeof(T));
		cdf = new T[n + 1];
		// Compute integral of step function at $x_i$
		cdf[0] = 0.;
		for (int i = 1; i < count + 1; ++i)
			cdf[i] = cdf[i - 1] + func[i - 1] / n;

		// Transform step function integral into CDF
		funcInt = cdf[count];
		if (funcInt == 0.f) {
			for (int i = 1; i < n + 1; ++i)
				cdf[i] = T(i) / T(n);
		}
		else {
			for (int i = 1; i < n + 1; ++i)
				cdf[i] /= funcInt;
		}
	}
	~Distribution1D() {
		delete[] func;
		delete[] cdf;
	}

	T SampleContinuous(T u, T *pdf, int *off = NULL) const 
	{
		// Find surrounding CDF segments and _offset_
		T *ptr = std::upper_bound(cdf, cdf + count + 1, u);
		int offset = std::max<int>(0, int(ptr - cdf - 1));
		if (off) *off = offset;
		assert(offset < count);
		assert(u >= cdf[offset] && u < cdf[offset + 1]);

		// Compute offset along CDF segment
		T du = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
		assert(!isnan(du));

		// Compute PDF for sampled offset
		if (pdf) *pdf = func[offset] / funcInt;

		// Return $x\in{}[0,1)$ corresponding to sample
		return (offset + du) / count;
	}

	int SampleDiscrete(T u, T *pdf) const 
	{
		// Find surrounding CDF segments and _offset_
		T *ptr = std::upper_bound(cdf, cdf + count + 1, u);
		int offset = std::max<int>(0, int(ptr - cdf - 1));
		assert(offset < count);
		assert(u >= cdf[offset] && u < cdf[offset + 1]);
		if (pdf) *pdf = func[offset] / (funcInt * count);
		return offset;
	}

private:
	friend struct Distribution2D;
	// Distribution1D Private Data
	T *func, *cdf;
	T funcInt;
	int count;
};


template <typename TypeSig, typename TypeProb>
void HistCount(const std::vector<TypeSig>& vSig, const std::vector<TypeProb>& vEdges, std::vector<TypeProb>& vProb);

#include "ProbTools.inl"