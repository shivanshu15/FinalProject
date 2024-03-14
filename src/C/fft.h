// #ifndef FFT_H
// #define FFT_H

// #if __cplusplus
// #   include <complex>
//     typedef std::complex< double > cplx;
// #else
// #   include <complex.h>
// #if defined(__GNUC__) || defined(__GNUG__)
// typedef double complex cplx;
// #elif defined(_MSC_VER)
// typedef _Dcomplex cplx;
// #endif
// #endif

// #include <stdlib.h>
// #include <string.h>
// #ifndef CMPLX
// #define CMPLX(x, y) ((cplx)((double)(x) + _Complex_I * (double)(y)))
// #endif
// extern void twiddles(cplx a[], int size);
// // extern void _fft(cplx a[], cplx out[], int size, int step, cplx tw[]);
// extern void fft(cplx a[], int size, cplx tw[]);
// extern void ifft(cplx a[], int size, cplx tw[]);
// #endif



#ifndef FFT_H
#define FFT_H

#if __cplusplus
#   include <complex>
    typedef std::complex<double> cplx;
#else
#   include <complex.h>
#if defined(__GNUC__) || defined(__GNUG__)
typedef double complex cplx;
#elif defined(_MSC_VER)
typedef _Dcomplex cplx;
#endif
#endif

#include <stdlib.h>
#include <string.h>

#ifndef CMPLX
#define CMPLX(x, y) ((cplx)((double)(x) + _Complex_I * (double)(y)))
#endif

typedef struct {
    cplx *F;
    int nFFT;
    cplx *tw;
} fft_parallel_arg;

// Prototype for the twiddles function
extern void twiddles(cplx a[], int size);

extern void *fft_parallel(void *args);

// Prototype for the parallel FFT function
extern void fft(cplx a[], int size);

#endif  // FFT_H

































// #pragma once

// #include <complex.h>
// #include <stdbool.h>
// #include <stddef.h>


// #ifdef __cplusplus
// extern "C" {
// #define restrict
// #endif


// /* 
//  * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
//  * The vector can have any length. This is a wrapper function. The inverse transform does not perform scaling,
//  * so it is not a true inverse. Returns true if successful, false otherwise (out of memory).
//  */
// bool Fft_transform(double complex vec[], size_t n, bool inverse);


// /* 
//  * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
//  * The vector's length must be a power of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
//  * Returns true if successful, false otherwise (n is not a power of 2, or out of memory).
//  */
// bool Fft_transformRadix2(double complex vec[], size_t n, bool inverse);


// /* 
//  * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
//  * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 FFT function.
//  * Uses Bluestein's chirp z-transform algorithm. Returns true if successful, false otherwise (out of memory).
//  */
// bool Fft_transformBluestein(double complex vec[], size_t n, bool inverse);


// /* 
//  * Computes the circular convolution of the given complex vectors. Each vector's length must be the same.
//  * Returns true if successful, false otherwise (out of memory).
//  */
// bool Fft_convolve(const double complex xvec[restrict], const double complex yvec[restrict], double complex outvec[restrict], size_t n);


// #ifdef __cplusplus
// #undef restrict
// }
// #endif