#pragma once

#include <cmath>
#include <complex>
#include <functional>

////Simple test:
//#include <iostream>
//#include <iomanip>
//#include "TFilonIntegrator.h"
//int main(int argc, char **argv)
//{
//	const double a = -1.;
//	const double b = 1.;
//	const double N = 11.;
//	const double w = 10.;
//	auto F = [](double x){return std::complex<double>(pow(x, 2.), 0.);};
//	TFilonIntegrator <double> filon_integrator;
//	std::complex<double> result = filon_integrator.Integrate(F, w, a, b, N);
//	std::cout<<std::fixed<<std::setprecision(10) << result << std::endl;		//-0.1401909989
//};

///Вычисляем интегралы вида
///b
///I f(x) * exp(iwx) dx
///a
///методом Филона:
///Михалкович С.С. Интегралы от быстро осциллирующих функций. Многомерные Интегралы. Методические указания к выполнению индивидуальных заданий на ЭВМ для студентов 2 курса физического факультета. Ростов-на-Дону:РГУ.-2000.-20с.
template <class fp>
class TFilonIntegrator {
protected:
	std::complex<fp> D(const fp &p, const short int j) const
	{
		if (p == 0.)
			return 0.;
		std::complex<fp> result;
		const fp	pow_p_3 = pow(p, fp(-3)),
							cos_p = cos(p),
							pow_p_2 = pow(p, fp(2)),
							sin_p = sin(p);

		switch ( j )
		{
			case 1:
				result.real( pow_p_3 * (fp(2) * p * cos_p - (fp(2) - pow_p_2) * sin_p) );
				result.imag( pow_p_3 * (pow_p_2 * cos_p - p * sin_p) );
				break;
			case 2:
				result.real( pow_p_3 * (fp(4) * sin_p - fp(4) * p * cos_p) );
				result.imag( 0 );
				break;
			case 3:
				result.real( pow_p_3 * (fp(2) * p * cos_p + (pow_p_2 - fp(2)) * sin_p) );
				result.imag( pow_p_3 * (p * sin_p - pow_p_2 * cos_p) );
		};
		return result;
	}
public:
	std::complex<fp> Integrate(
		std::function<std::complex<fp>(const fp&)> f,	///std::complex<fp> f( const fp &)
		const fp &w,								///множитель в экспоненте перед переменной интегрирования
		const fp &a,								///нижний предел интегрирования
		const fp &b,								///верхний предел интегрирования
		const int N									///Количество отрезков разбиения
	) const {
		std::complex<fp> result = 0.;

		const fp h = (b - a) / N;

		for (int i = 1; i <= N; i++)		//i = 1
		{
			fp	ai = a + h * (i-1),
					bi = a + h * i;
			fp arg = w * (ai + bi)/fp(2.);
			fp mid = (bi - ai)/fp(2.);
			std::complex<fp> exponent = {cos(arg), sin(arg)};
			for (int j = 1; j <= 3; j++)
			{
				//Ключевой момент здесь автоматически раскрываются скобки и перемножаются все комплексные числа
				result += exponent * D(w * mid, j) * mid * f(ai + (j - 1) * mid);
			};
		};

		return result;
	}			//Integrate
};
