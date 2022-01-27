#pragma once

#include <complex>

typedef double fp;

template <class TComplexContainer, class TRealContainer>
inline void CA2DATA(const TComplexContainer &container, TRealContainer &re_buffer, TRealContainer &im_buffer)
{
	#pragma omp parallel for
	for (size_t idv = 0; idv < container.size(); idv++)
	{
		re_buffer[idv] = container[idv].real();
		im_buffer[idv] = container[idv].imag();
	};
}

