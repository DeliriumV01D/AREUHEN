#include "TFFT.h"
#include "TThread.h"

TFFT :: TFFT(const TFFTProperties &fft_properties)
{
	SetFFTProperties(fft_properties);
}

void TFFT :: SetFFTProperties(const TFFTProperties &fft_properties) noexcept
{
	FFTProperties = fft_properties;
}

TFFTProperties TFFT :: GetFFTProperties() const noexcept
{
	return FFTProperties;
}

void TFFT :: operator () (const Double_t * re_in, const Double_t * im_in, Double_t * re_out, Double_t * im_out)
{
	if (FFTProperties.FFTOperation != TFFTOperation::C2C_FORWARD && FFTProperties.FFTOperation != TFFTOperation::C2C_BACKWARD)
		throw std::runtime_error("TFFT :: operator () error: incorrect FFTOperation");
	
	//Проблема в том, что TVirtualFFT использует глобальный gFFT
	//поэтому он не годится для многопоточного использования
	//на данный момент приходится использовать глобальный мьютекс Lock();Unlock();
	std::unique_ptr<TVirtualFFT> fft;
	TThread::Lock();
	try {		
		if (this->FFTProperties.FFTOperation == TFFTOperation::C2C_FORWARD)
		{
			fft = std::unique_ptr<TVirtualFFT>(TVirtualFFT::FFT(1, &FFTProperties.N, "C2C ES K"));
      fft->SetPointsComplex(re_in, im_in);
			fft->Transform();
			fft->GetPointsComplex(re_out, im_out);
    } else {
			fft = std::unique_ptr<TVirtualFFT>(TVirtualFFT::FFT(1, &FFTProperties.N, "C2R ES K"));
      fft->SetPointsComplex(re_in, im_in);
			fft->Transform();
      fft->GetPointsComplex(re_out, im_out);
			if (FFTProperties.FFTNormalization == TFFTNormalization::BACKWARD)
			{
				//Нужно *1/Nbins - тоже самое, что 1/2pi при непрерывном преобразовании
				#pragma omp parallel for
				for (int i = 0; i < FFTProperties.N; i++)
				{
					re_out[i] = re_out[i]/FFTProperties.N;
					im_out[i] = im_out[i]/FFTProperties.N;
				}
			}
		}
		if (FFTProperties.FFTNormalization == TFFTNormalization::TWO_SIDED)
		{
			//Нужно *1/sqrt(Nbins) - тоже самое, что 1/sqrt(2pi) при непрерывном преобразовании
			const Double_t sqrtn = std::sqrt(static_cast<Double_t>(FFTProperties.N));
			#pragma omp parallel for
			for (int i = 0; i < FFTProperties.N; i++)
			{
				re_out[i] = re_out[i]/sqrtn;
				im_out[i] = im_out[i]/sqrtn;
			}
		}
  } catch (std::exception &e) {
		TThread::UnLock();
		throw e;
	}
	TThread::UnLock();
}


void TFFT :: operator () (Double_t * re_c, Double_t * im_c, Double_t * r)
{
	if (FFTProperties.FFTOperation != TFFTOperation::R2C_FORWARD && FFTProperties.FFTOperation != TFFTOperation::C2R_BACKWARD)
    throw std::runtime_error("TFFT :: operator () error: incorrect FFTOperation");
	
	//Проблема в том, что TVirtualFFT использует глобальный gFFT
	//поэтому он не годится для многопоточного использования
	//на данный момент приходится использовать глобальный мьютекс Lock();Unlock();
	std::unique_ptr<TVirtualFFT> fft;
	TThread::Lock();
	try {		
		//const Double_t SQRTN = sqrt(Double_t(NPnt));
		if (FFTProperties.FFTOperation == TFFTOperation::R2C_FORWARD)
		{
			fft = std::unique_ptr<TVirtualFFT>(TVirtualFFT::FFT(1, &FFTProperties.N, "R2C ES K"));
      fft->SetPoints(r);
			fft->Transform();
      fft->GetPointsComplex(re_c, im_c);
			if (FFTProperties.FFTNormalization == TFFTNormalization::TWO_SIDED)
			{
				const Double_t sqrtn = std::sqrt(static_cast<Double_t>(FFTProperties.N));
				//Нужно *1/sqrt(Nbins) - тоже самое, что 1/sqrt(2pi) при непрерывном преобразовании
				#pragma omp parallel for
				for (int i = 0; i < FFTProperties.N; i++)
				{
					re_c[i] = re_c[i]/sqrtn;
					im_c[i] = im_c[i]/sqrtn;
				}
			}
    } else {
			fft = std::unique_ptr<TVirtualFFT>(TVirtualFFT::FFT(1, &FFTProperties.N, "C2R ES K"));
      fft->SetPointsComplex(re_c, im_c);
			fft->Transform();
      fft->GetPoints(r);
			if (FFTProperties.FFTNormalization == TFFTNormalization::TWO_SIDED)
			{
				//Нужно *1/sqrt(Nbins) - тоже самое, что 1/sqrt(2pi) при непрерывном преобразовании
				const Double_t sn = std::sqrt(static_cast<Double_t>(FFTProperties.N));
				#pragma omp parallel for
				for (int i = 0; i < FFTProperties.N; i++)
					r[i] = r[i]/sn;
			}
			if (FFTProperties.FFTNormalization == TFFTNormalization::BACKWARD)
			{
				//Нужно *1/Nbins - тоже самое, что 1/2pi при непрерывном преобразовании
				const Double_t sn = static_cast<Double_t>(FFTProperties.N);
				#pragma omp parallel for
				for (int i = 0; i < FFTProperties.N; i++)
					r[i] = r[i]/sn;
			}
		}
  } catch (std::exception &e) {
		TThread::UnLock();
		throw e;
	}
	TThread::UnLock();
}

//void TFFT :: operator () (const Double_t * r_in, Double_t * r_out)
//{
//	if (FFTOperation != TFFTOperation::R2R_FORWARD && FFTOperation != TFFTOperation::R2R_BACKWARD)
//		throw std::runtime_error("TFFT :: operator () error: incorrect FFTOperation");
//}
