#if __cplusplus
#   include <complex>
typedef std::complex< double > cplx;
#else
#   include <complex.h>
#if defined(__GNUC__) || defined(__GNUG__)
typedef double complex cplx;
#elif defined(_MSC_VER)
typedef _Dcomplex cplx;
#endif
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>

#include "stats.h"
#include "fft.h"
#include "histcounts.h"

#include "helper_functions.h"

// #ifndef CMPLX
// #define CMPLX(x, y) ((cplx)((double)(x) + _Imaginary_I * (double)(y)))
// #endif

#ifndef CMPLX
#define CMPLX(x, y) ((cplx)((double)(x) + (double)(y) * I))
#endif

#define pow2(x) (1 << x)

int nextpow2(int n)
{
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
}

/*
static void apply_conj(cplx a[], int size, int normalize)
{
    switch(normalize) {
        case(1):
            for (int i = 0; i < size; i++) {
                a[i] = conj(a[i]) / size;
            }
            break;
        default:
            for (int i = 0; i < size; i++) {
                a[i] = conj(a[i]);
            }
            break;
    }
}
 */

void dot_multiply(cplx a[], cplx b[], int size)
{
    for (int i = 0; i < size; i++) {
        a[i] = _Cmulcc(a[i], conj(b[i]));
    }
}

// double * CO_AutoCorr(const double y[], const int size, const int tau[], const int tau_size)
// {
//     double m, nFFT;
//     m = mean(y, size);
//     nFFT = nextpow2(size) << 1;

//     cplx * F = malloc(nFFT * sizeof *F);
//     cplx * tw = malloc(nFFT * sizeof *tw);
//     for (int i = 0; i < size; i++) {

//         #if defined(__GNUC__) || defined(__GNUG__)
//                 F[i] = CMPLX(y[i] - m, 0.0);
//         #elif defined(_MSC_VER)
//                 cplx tmp = { y[i] - m, 0.0 };
//                 F[i] = tmp;
//         #endif

//     }
//     for (int i = size; i < nFFT; i++) {
//         #if defined(__GNUC__) || defined(__GNUG__)
//                 F[i] = CMPLX(0.0, 0.0);
//         #elif defined(_MSC_VER)
//                 cplx tmp = { 0.0, 0.0 };
//                 F[i] = tmp; // CMPLX(0.0, 0.0);
//         #endif

//     }
//     //size = nFFT;

//     // twiddles(tw, nFFT);
//     // //fft(F, nFFT, tw);
//     //     fft_parallel_arg arg = {F, nFFT, tw};
//     // fft_parallel(&arg);
//     // dot_multiply(F, F, nFFT);
//     // //fft(F, nFFT, tw);
//     //     fft_parallel_arg arg = {F, nFFT, tw};
//     // fft_parallel(&arg);

//    // printf("%s\n", "Calculating twiddles");
//     twiddles(tw, nFFT);
    
//     //printf("%s\n", "Calculating fft");
    
//     fft(F, nFFT); // Use the fft function for parallel processing
//     dot_multiply(F, F, nFFT);
//     fft(F, nFFT); // Use the fft function again after dot_multiply

//     cplx divisor = F[0];
//     for (int i = 0; i < nFFT; i++) {
//         //F[i] = F[i] / divisor;
//         F[i] = _Cdivcc(F[i], divisor);
//     }

//     // double * out = malloc(tau_size * sizeof(out));
//     // for (int i = 0; i < tau_size; i++) {
//     //     out[i] = creal(F[tau[i]]);

//     double *out = malloc(tau_size * sizeof(double)); // Correct allocation size
//     for (int i = 0; i < tau_size; i++) {
//         out[i] = creal(F[tau[i]]);
//     }
    
//     free(F);
//     free(tw);
//     return out;
// }




double * CO_AutoCorr(const double y[], const int size, const int tau[], const int tau_size)
{
    double m = mean(y, size);
    int nFFT = nextpow2(size) << 1;

    cplx *F = (cplx *)malloc(nFFT * sizeof(cplx));
    for (int i = 0; i < size; i++) {
        F[i] = CMPLX(y[i] - m, 0.0);
    }
    for (int i = size; i < nFFT; i++) {
        F[i] = CMPLX(0.0, 0.0);
    }

    fft(F, nFFT); // Use FFTW-based fft
    dot_multiply(F, F, nFFT);
    ifft(F, nFFT); // Use FFTW-based ifft

    // Normalize the result of IFFT by nFFT
    for (int i = 0; i < nFFT; i++) {
        F[i] /= nFFT;
    }

    double *out = (double *)malloc(tau_size * sizeof(double));
    for (int i = 0; i < tau_size; i++) {
        out[i] = creal(F[tau[i]]);
    }

    free(F);
    return out;
}


// double * co_autocorrs(const double y[], const int size)
// {
//     //printf("%s\n", "Calculating co_autocorr");
//     double m, nFFT;
//     m = mean(y, size);
//     nFFT = nextpow2(size) << 1;

//     cplx * F = malloc(nFFT * 2 * sizeof *F);
//     cplx * tw = malloc(nFFT * 2 * sizeof *tw);
//     for (int i = 0; i < size; i++) {

//         #if defined(__GNUC__) || defined(__GNUG__)
//                 F[i] = CMPLX(y[i] - m, 0.0);
//         #elif defined(_MSC_VER)
//                 cplx tmp = { y[i] - m, 0.0 };
//                 F[i] = tmp;
//         #endif
//     }
//     for (int i = size; i < nFFT; i++) {

//         #if defined(__GNUC__) || defined(__GNUG__)
//             F[i] = CMPLX(0.0, 0.0);
//         #elif defined(_MSC_VER)
//             cplx tmp = { 0.0, 0.0 };
//             F[i] = tmp;
//         #endif
//     }
//     //size = nFFT;

//     // twiddles(tw, nFFT);
//     // //fft(F, nFFT, tw);
//     // fft_parallel_arg arg = {F, nFFT, tw};
//     // fft_parallel(&arg);
//     // //parallel_fft(F, nFFT, tw);

//     // dot_multiply(F, F, nFFT);
//     // //fft(F, nFFT, tw);
//     // fft_parallel_arg arg = {F, nFFT, tw};
//     // fft_parallel(&arg);

//     twiddles(tw, nFFT);
//     fft(F, nFFT); // Use the fft function for parallel processing
//     for (int i = 0; i < nFFT; i++) {
//     //printf("After FFT, F[%d] = %f + %fi\n", i, creal(F[i]), cimag(F[i]));
// }
//     dot_multiply(F, F, nFFT);
//     fft(F, nFFT); // Use the fft function again after dot_multiply

//     cplx divisor = F[0];
//     for (int i = 0; i < nFFT; i++) {
//         F[i] = _Cdivcc(F[i], divisor); // F[i] / divisor;
//     }

//     // double * out = malloc(nFFT * 2 * sizeof(out));
//     // for (int i = 0; i < nFFT; i++) {  
//     //     out[i] = creal(F[i]);
//     // }

//         double *out = malloc(nFFT * sizeof(double)); // Correct allocation size
//     for (int i = 0; i < nFFT; i++) {
//         out[i] = creal(F[i]);
//     }
    
//     free(F);
//     free(tw);
//     return out;
// }





double * co_autocorrs(const double y[], const int size) {
    double m = mean(y, size); 
    int nFFT = nextpow2(size) << 1; 

    cplx *F = (cplx *)malloc(nFFT * sizeof(cplx));

    for (int i = 0; i < size; i++) {
        F[i] = CMPLX(y[i] - m, 0.0); 
    }
    for (int i = size; i < nFFT; i++) {
        F[i] = CMPLX(0.0, 0.0);
    }

    fft(F, nFFT); //FFTW3

    for (int i = 0; i < nFFT; i++) {
        F[i] *= conj(F[i]);
    }

    ifft(F, nFFT);

    double *out = (double *)malloc(size * sizeof(double));

    for (int i = 0; i < size; i++) {
        out[i] = creal(F[i]) / size;
    }

    free(F);

    return out;
}


//int co_firstzero(const double y[], const int size, const int maxtau, const double * autocorrs)
int co_firstzero(const double y[], const int size, const int maxtau,const double * autocorrs)
{

    // double * autocorrs = malloc(size * sizeof * autocorrs);
    // autocorrs = co_autocorrs(y, size);

    //double * autocorrs = autocorrs;
     //double * autocorrs = co_autocorrs(y, size);
    int zerocrossind = 0;
    while(autocorrs[zerocrossind] > 0 && zerocrossind < maxtau)
    {
        zerocrossind += 1;
    }

    //free(autocorrs);
    return zerocrossind;

}


#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NUM_THREADS 8 // You can adjust the number of threads

typedef struct {
    const double* y;
    const double* autocorrs;
    int size;
    double thresh;
    int startIndex;
    int endIndex;
    double result;
} ThreadData;

void* threadFunction(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    double localResult = (double)data->size;
    for (int i = data->startIndex; i <= data->endIndex && i < data->size - 2; i++) {
        if (data->autocorrs[i+1] < data->thresh) {
            double m = data->autocorrs[i+1] - data->autocorrs[i];
            double dy = data->thresh - data->autocorrs[i];
            double dx = dy / m;
            localResult = ((double)i) + dx;
            break; // Found the threshold breach, exit loop
        }
    }
    data->result = localResult;
    pthread_exit(NULL);
}

double CO_f1ecac(const double y[], const int size, const double * autocorrs) {
    double thresh = 1.0 / exp(1);
    pthread_t threads[NUM_THREADS];
    ThreadData threadData[NUM_THREADS];
    int chunkSize = (size - 2) / NUM_THREADS;
    double globalResult = (double)size;

    for (int i = 0; i < NUM_THREADS; i++) {
        threadData[i].y = y;
        threadData[i].autocorrs = autocorrs;
        threadData[i].size = size;
        threadData[i].thresh = thresh;
        threadData[i].startIndex = i * chunkSize;
        threadData[i].endIndex = (i + 1) * chunkSize - 1;
        if (i == NUM_THREADS - 1) threadData[i].endIndex = size - 3; // Last thread takes any remaining elements
        int rc = pthread_create(&threads[i], NULL, threadFunction, (void*)&threadData[i]);
        if (rc) {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
        if (threadData[i].result < globalResult) {
            globalResult = threadData[i].result;
        }
    }

    return globalResult;
}

    
   






// double CO_f1ecac(const double y[], const int size, const double * autocorrs)
// {

//     // NaN check
//     for(int i = 0; i < size; i++)
//     {
//         if(isnan(y[i]))
//         {
//             return 0;
//         }
//     }

//     // compute autocorrelations
//     //double * autocorrs = co_autocorrs(y, size);
//     //double * autocorrs = autocorrs;
//     // threshold to cross
//     double thresh = 1.0/exp(1);

//     double out = (double)size;
//     for(int i = 0; i < size-2; i++){
//         // printf("i=%d autocorrs_i=%1.3f\n", i, autocorrs[i]);
//         if ( autocorrs[i+1] < thresh ){
//             double m = autocorrs[i+1] - autocorrs[i];
//             double dy = thresh - autocorrs[i];
//             double dx = dy/m;
//             out = ((double)i) + dx;
//             // printf("thresh=%1.3f AC(i)=%1.3f AC(i-1)=%1.3f m=%1.3f dy=%1.3f dx=%1.3f out=%1.3f\n", thresh, autocorrs[i], autocorrs[i-1], m, dy, dx, out);
//             //free(autocorrs);
//             return out;
//         }
//     }

//     //free(autocorrs);

//     return out;

// }

double CO_Embed2_Basic_tau_incircle(const double y[], const int size, const double radius, const int tau,const double * autocorrs)
{
    int tauIntern = 0;

    if(tau < 0)
    {
        tauIntern = co_firstzero(y, size, size,autocorrs);
    }
    else{
        tauIntern = tau;
    }

    double insidecount = 0;
    for(int i = 0; i < size-tauIntern; i++)
    {
        if(y[i]*y[i] + y[i+tauIntern]*y[i+tauIntern] < radius)
        {
            insidecount += 1;
        }
    }

    return insidecount/(size-tauIntern);
}

// double CO_Embed2_Dist_tau_d_expfit_meandiff(const double y[], const int size,const double * autocorrs)
// {

//     // NaN check
//     for(int i = 0; i < size; i++)
//     {
//         if(isnan(y[i]))
//         {
//             return NAN;
//         }
//     }

//     int tau = co_firstzero(y, size, size,autocorrs);

//     //printf("co_firstzero ran\n");

//     if (tau > (double)size/10){
//         tau = floor((double)size/10);
//     }
//     //printf("tau = %i\n", tau);

//     double * d = malloc((size-tau) * sizeof(double));
//     for(int i = 0; i < size-tau-1; i++)
//     {

//         d[i] = sqrt((y[i+1]-y[i])*(y[i+1]-y[i]) + (y[i+tau]-y[i+tau+1])*(y[i+tau]-y[i+tau+1]));

//         //printf("d[%i]: %1.3f\n", i, d[i]);
//         if (isnan(d[i])){
//             free(d);
//             return NAN;
//         }

//     }


//     // mean for exponential fit
//     double l = mean(d, size-tau-1);


//     int nBins = num_bins_auto(d, size-tau-1);
//     if (nBins == 0){
//         return 0;
//     }
//     int * histCounts = malloc(nBins * sizeof(double));
//     double * binEdges = malloc((nBins + 1) * sizeof(double));
//     histcounts_preallocated(d, size-tau-1, nBins, histCounts, binEdges);


//     // normalise to probability
//     double * histCountsNorm = malloc(nBins * sizeof(double));
//     for(int i = 0; i < nBins; i++){
//         //printf("histCounts %i: %i\n", i, histCounts[i]);
//         histCountsNorm[i] = (double)histCounts[i]/(double)(size-tau-1);
//         //printf("histCounts norm %i: %1.3f\n", i, histCountsNorm[i]);
//     }



//     double * d_expfit_diff = malloc(nBins * sizeof(double));
//     for(int i = 0; i < nBins; i++){
//         double expf = exp(-(binEdges[i] + binEdges[i+1])*0.5/l)/l;
//         if (expf < 0){
//             expf = 0;
//         }
//         d_expfit_diff[i] = fabs(histCountsNorm[i]-expf);
//         //printf("d_expfit_diff %i: %1.3f\n", i, d_expfit_diff[i]);
//     }

//     double out = mean(d_expfit_diff, nBins);


//     // arrays created dynamically in function histcounts
//     free(d);
//     free(d_expfit_diff);
//     free(binEdges);
//     free(histCountsNorm);
//     free(histCounts);

//     return out;

// }


#include <pthread.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

typedef struct {
    const double* y;
    double* d;
    int start;
    int end;
    int tau;
} DistanceCalcArgs;

typedef struct {
    const int* histCounts;
    double* histCountsNorm;
    int start;
    int end;
    int size_tau_1;
} NormalizeArgs;

void* calcDistances(void* args) {
    DistanceCalcArgs* data = (DistanceCalcArgs*)args;
    for (int i = data->start; i < data->end; i++) {
        data->d[i] = sqrt((data->y[i+1]-data->y[i])*(data->y[i+1]-data->y[i]) + (data->y[i+data->tau]-data->y[i+data->tau+1])*(data->y[i+data->tau]-data->y[i+data->tau+1]));
    }
    pthread_exit(NULL);
}

void* normalizeHistCounts(void* args) {
    NormalizeArgs* data = (NormalizeArgs*)args;
    for (int i = data->start; i < data->end; i++) {
        data->histCountsNorm[i] = (double)data->histCounts[i] / (double)data->size_tau_1;
    }
    pthread_exit(NULL);
}

double CO_Embed2_Dist_tau_d_expfit_meandiff(const double y[], const int size, const double* autocorrs) {
    


    for(int i = 0; i < size; i++)
    {
        if(isnan(y[i]))
        {
            return NAN;
        }
    }

    int tau = co_firstzero(y, size, size, autocorrs);
    if (tau > (double)size / 10) {
        tau = floor((double)size / 10);
    }

    double* d = malloc((size - tau) * sizeof(double));
    if (d == NULL) {
        
        return NAN;
    }

    int numThreads = 4; 
    pthread_t threads[numThreads];
    DistanceCalcArgs threadArgs[numThreads];
    int chunkSize = (size - tau - 1) / numThreads;

    for (int i = 0; i < numThreads; i++) {
        threadArgs[i].y = y;
        threadArgs[i].d = d;
        threadArgs[i].start = i * chunkSize;
        threadArgs[i].end = (i == numThreads - 1) ? size - tau - 1 : (i + 1) * chunkSize;
        threadArgs[i].tau = tau;
        pthread_create(&threads[i], NULL, calcDistances, (void*)&threadArgs[i]);
    }

    for (int i = 0; i < numThreads; i++) {
        pthread_join(threads[i], NULL);
    }

   
    double l = mean(d, size-tau-1);
  
    int nBins = num_bins_auto(d, size - tau - 1);
    double* histCountsNorm = malloc(nBins * sizeof(double));
    int* histCounts = malloc(nBins * sizeof(int));
    double* binEdges = malloc((nBins + 1) * sizeof(double));
    histcounts_preallocated(d, size-tau-1, nBins, histCounts, binEdges);
    if (histCountsNorm == NULL) {
        // Handle malloc failure
        free(d);
        return NAN;
    }

    NormalizeArgs normArgs[numThreads];
    chunkSize = nBins / numThreads;

    for (int i = 0; i < numThreads; i++) {
        normArgs[i].histCounts = histCounts;
        normArgs[i].histCountsNorm = histCountsNorm;
        normArgs[i].start = i * chunkSize;
        normArgs[i].end = (i == numThreads - 1) ? nBins : (i + 1) * chunkSize;
        normArgs[i].size_tau_1 = size - tau - 1;
        pthread_create(&threads[i], NULL, normalizeHistCounts, (void*)&normArgs[i]);
    }

    for (int i = 0; i < numThreads; i++) {
        pthread_join(threads[i], NULL);
    }

   

        double * d_expfit_diff = malloc(nBins * sizeof(double));
    for(int i = 0; i < nBins; i++){
        double expf = exp(-(binEdges[i] + binEdges[i+1])*0.5/l)/l;
        if (expf < 0){
            expf = 0;
        }
        d_expfit_diff[i] = fabs(histCountsNorm[i]-expf);
        //printf("d_expfit_diff %i: %1.3f\n", i, d_expfit_diff[i]);
    }

    double out = mean(d_expfit_diff, nBins);


    free(d);
    free(d_expfit_diff);
    free(binEdges);
    free(histCountsNorm);
    free(histCounts);
    return out;
}







#include <pthread.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#define NUM_THREADS 8 

typedef struct {
    const double* autocorrs;
    int start;
    int end;
    int size;
    int thread_min_index;
} ThreadData0;

void* find_min(void* arg) {
    ThreadData0* data = (ThreadData0*)arg;
    int minInd = data->size;
    for (int i = data->start; i < data->end - 1 && i < data->size - 1; i++) {
        if (data->autocorrs[i] < data->autocorrs[i - 1] && data->autocorrs[i] < data->autocorrs[i + 1]) {
            minInd = i;
            break;
        }
    }
    data->thread_min_index = minInd;
    pthread_exit(NULL);
}

int CO_FirstMin_ac(const double y[], const int size, const double* autocorrs) {
    pthread_t threads[NUM_THREADS];
    ThreadData0 thread_data[NUM_THREADS];
    int segment_size = size / NUM_THREADS;

    for (int i = 0; i < NUM_THREADS; i++) {
        thread_data[i].autocorrs = autocorrs;
        thread_data[i].size = size;
        thread_data[i].start = i * segment_size;
        thread_data[i].end = (i + 1) * segment_size;
        thread_data[i].thread_min_index = size;

        int rc = pthread_create(&threads[i], NULL, find_min, (void*)&thread_data[i]);
        if (rc) {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    int global_min_index = size;
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
        if (thread_data[i].thread_min_index < global_min_index) {
            global_min_index = thread_data[i].thread_min_index;
        }
    }

    return global_min_index;
}




// int CO_FirstMin_ac(const double y[], const int size, const double * autocorrs)
// {

//     // NaN check
//     for(int i = 0; i < size; i++)
//     {
//         if(isnan(y[i]))
//         {
//             return 0;
//         }
//     }

//     //double * autocorrs = co_autocorrs(y, size);
//     //double * autocorrs = autocorrs;

//     int minInd = size;
//     for(int i = 1; i < size-1; i++)
//     {
//         if(autocorrs[i] < autocorrs[i-1] && autocorrs[i] < autocorrs[i+1])
//         {
//             minInd = i;
//             break;
//         }
//     }

//     //free(autocorrs);

//     return minInd;

// }


#include <math.h>
#include <stdbool.h>


#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    const double* y;
    int start;
    int end;
    double sumCubes;
} ThreadData1;

void* threadFunction1(void* arg) {
    ThreadData1* data = (ThreadData1*)arg;
    data->sumCubes = 0;
    for(int i = data->start; i < data->end; i++) {
        double diff = data->y[i+1] - data->y[i];
        data->sumCubes += diff * diff * diff;
    }
    pthread_exit(NULL);
}

double CO_trev_1_num(const double y[], const int size) {
    int numThreads = 8;
    pthread_t threads[numThreads];
    ThreadData1 threadData[numThreads];
    int chunkSize = (size - 1) / numThreads;
    double totalSumCubes = 0;

    // Create threads
    for (int i = 0; i < numThreads; ++i) {
        threadData[i].y = y;
        threadData[i].start = i * chunkSize;
        threadData[i].end = (i == numThreads - 1) ? (size - 1) : (threadData[i].start + chunkSize);
        if(pthread_create(&threads[i], NULL, threadFunction1, (void*)&threadData[i])) {
            fprintf(stderr, "Error creating thread\n");
            return 1;
        }
    }

    // Join threads
    for (int i = 0; i < numThreads; ++i) {
        pthread_join(threads[i], NULL);
        totalSumCubes += threadData[i].sumCubes;
    }

    // Calculate mean directly
    double mean = totalSumCubes / (size - 1);
    return mean;
}

// Don't forget to link with -lpthread when you compile this program.


// double CO_trev_1_num(const double y[], const int size)
// {
//     // Check for NaN in the input array
//     for(int i = 0; i < size; i++)
//     {
//         if(isnan(y[i]))
//         {
//             return NAN;
//         }
//     }

//     // NaN check passed, compute sum of cubes
//     double sumCubes = 0;
//     for(int i = 0; i < size - 1; i++)
//     {
//         double diff = y[i+1] - y[i];
//         sumCubes += diff * diff * diff; // Replace pow with direct multiplication
//     }

//     // Calculate mean directly
//     double mean = sumCubes / (size - 1);
//     return mean;
// }



// double CO_trev_1_num(const double y[], const int size)
// {

//     // NaN check
//     for(int i = 0; i < size; i++)
//     {
//         if(isnan(y[i]))
//         {
//             return NAN;
//         }
//     }

//     int tau = 1;

//     double * diffTemp = malloc((size-1) * sizeof * diffTemp);

//     for(int i = 0; i < size-tau; i++)
//     {
//         diffTemp[i] = pow(y[i+1] - y[i],3);
//     }

//     double out;

//     out = mean(diffTemp, size-tau);

//     free(diffTemp);

//     return out;
// }

#define tau 2
#define numBins 5

double CO_HistogramAMI_even_2_5(const double y[], const int size)
{
    // NaN check
    for (int i = 0; i < size; i++) {
        if (isnan(y[i])) {
            return NAN;
        }
    }

    // const int tau = 2;
    // const int numBins = 5;

    // Allocate memory for y1 and y2
    double *y1 = malloc((size - tau) * sizeof(double));
    double *y2 = malloc((size - tau) * sizeof(double));

    // Populate y1 and y2
    for (int i = 0; i < size - tau; i++) {
        y1[i] = y[i];
        y2[i] = y[i + tau];
    }

    // Calculate min and max values of y
    double minValue = y[0];
    double maxValue = y[0];
    for (int i = 1; i < size; i++) {
        if (y[i] < minValue) {
            minValue = y[i];
        }
        if (y[i] > maxValue) {
            maxValue = y[i];
        }
    }

    double binStep = (maxValue - minValue + 0.2) / numBins;
    double binEdges[numBins + 1];
    for (int i = 0; i < numBins + 1; i++) {
        binEdges[i] = minValue + binStep * i - 0.1;
    }

    int *bins1 = histbinassign(y1, size - tau, binEdges, numBins + 1);
    int *bins2 = histbinassign(y2, size - tau, binEdges, numBins + 1);

    int jointHistLinear[(numBins + 1) * (numBins + 1)] = {0};

    for (int i = 0; i < size - tau; i++) {
        int binIndex1 = bins1[i];
        int binIndex2 = bins2[i];
        jointHistLinear[binIndex1 * (numBins + 1) + binIndex2]++;
    }

    double pij[numBins][numBins];
    int sumBins = 0;
    for (int i = 0; i < numBins; i++) {
        for (int j = 0; j < numBins; j++) {
            pij[j][i] = jointHistLinear[i * (numBins + 1) + j];
            sumBins += pij[j][i];
        }
    }

    for (int i = 0; i < numBins; i++) {
        for (int j = 0; j < numBins; j++) {
            pij[j][i] /= sumBins;
        }
    }

    double pi[numBins] = {0};
    double pj[numBins] = {0};
    for (int i = 0; i < numBins; i++) {
        for (int j = 0; j < numBins; j++) {
            pi[i] += pij[i][j];
            pj[j] += pij[i][j];
        }
    }

    double ami = 0;
    for (int i = 0; i < numBins; i++) {
        for (int j = 0; j < numBins; j++) {
            if (pij[i][j] > 0) {
                ami += pij[i][j] * log(pij[i][j] / (pj[j] * pi[i]));
            }
        }
    }

    free(bins1);
    free(bins2);

    free(y1);
    free(y2);

    return ami;
}

// double CO_HistogramAMI_even_2_5(const double y[], const int size)
// {

//     // NaN check
//     for(int i = 0; i < size; i++)
//     {
//         if(isnan(y[i]))
//         {
//             return NAN;
//         }
//     }

//     //const int tau = 2;
//     //const int numBins = 5;

//     double * y1 = malloc((size-tau) * sizeof(double));
//     double * y2 = malloc((size-tau) * sizeof(double));

//     for(int i = 0; i < size-tau; i++){
//         y1[i] = y[i];
//         y2[i] = y[i+tau];
//     }

//     // set bin edges
//     const double maxValue = max_(y, size);
//     const double minValue = min_(y, size);

//     double binStep = (maxValue - minValue + 0.2)/5;
//     //double binEdges[numBins+1] = {0};
// 	double binEdges[5+1] = {0};
//     for(int i = 0; i < numBins+1; i++){
//         binEdges[i] = minValue + binStep*i - 0.1;
//         // printf("binEdges[%i] = %1.3f\n", i, binEdges[i]);
//     }


//     // count histogram bin contents
//     int * bins1;
//     bins1 = histbinassign(y1, size-tau, binEdges, numBins+1);

//     int * bins2;
//     bins2 = histbinassign(y2, size-tau, binEdges, numBins+1);

//     /*
//     // debug
//     for(int i = 0; i < size-tau; i++){
//         printf("bins1[%i] = %i, bins2[%i] = %i\n", i, bins1[i], i, bins2[i]);
//     }
//     */

//     // joint
//     double * bins12 = malloc((size-tau) * sizeof(double));
//     //double binEdges12[(numBins + 1) * (numBins + 1)] = {0};
// 	double binEdges12[(5 + 1) * (5 + 1)] = {0};

//     for(int i = 0; i < size-tau; i++){
//         bins12[i] = (bins1[i]-1)*(numBins+1) + bins2[i];
//         // printf("bins12[%i] = %1.3f\n", i, bins12[i]);
//     }

//     for(int i = 0; i < (numBins+1)*(numBins+1); i++){
//         binEdges12[i] = i+1;
//         // printf("binEdges12[%i] = %1.3f\n", i, binEdges12[i]);
//     }

//     // fancy solution for joint histogram here
//     int * jointHistLinear;
//     jointHistLinear = histcount_edges(bins12, size-tau, binEdges12, (numBins + 1) * (numBins + 1));

//     /*
//     // debug
//     for(int i = 0; i < (numBins+1)*(numBins+1); i++){
//         printf("jointHistLinear[%i] = %i\n", i, jointHistLinear[i]);
//     }
//     */

//     // transfer to 2D histogram (no last bin, as in original implementation)
//     double pij[numBins][numBins];
//     int sumBins = 0;
//     for(int i = 0; i < numBins; i++){
//         for(int j = 0; j < numBins; j++){
//             pij[j][i] = jointHistLinear[i*(numBins+1)+j];

//             // printf("pij[%i][%i]=%1.3f\n", i, j, pij[i][j]);

//             sumBins += pij[j][i];
//         }
//     }

//     // normalise
//     for(int i = 0; i < numBins; i++){
//         for(int j = 0; j < numBins; j++){
//             pij[j][i] /= sumBins;
//         }
//     }

//     // marginals
//     //double pi[numBins] = {0};
// 	double pi[5] = {0};
//     //double pj[numBins] = {0};
// 	double pj[5] = {0};
//     for(int i = 0; i < numBins; i++){
//         for(int j = 0; j < numBins; j++){
//             pi[i] += pij[i][j];
//             pj[j] += pij[i][j];
//             // printf("pij[%i][%i]=%1.3f, pi[%i]=%1.3f, pj[%i]=%1.3f\n", i, j, pij[i][j], i, pi[i], j, pj[j]);
//         }
//     }

//     /*
//     // debug
//     for(int i = 0; i < numBins; i++){
//         printf("pi[%i]=%1.3f, pj[%i]=%1.3f\n", i, pi[i], i, pj[i]);
//     }
//     */

//     // mutual information
//     double ami = 0;
//     for(int i = 0; i < numBins; i++){
//         for(int j = 0; j < numBins; j++){
//             if(pij[i][j] > 0){
//                 //printf("pij[%i][%i]=%1.3f, pi[%i]=%1.3f, pj[%i]=%1.3f, logarg=, %1.3f, log(...)=%1.3f\n",
//                 //       i, j, pij[i][j], i, pi[i], j, pj[j], pij[i][j]/(pi[i]*pj[j]), log(pij[i][j]/(pi[i]*pj[j])));
//                 ami += pij[i][j] * log(pij[i][j]/(pj[j]*pi[i]));
//             }
//         }
//     }

//     free(bins1);
//     free(bins2);
//     free(jointHistLinear);

//     free(y1);
//     free(y2);
//     free(bins12);

//     return ami;
// }
