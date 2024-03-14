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

    for (int i = 0; i < numThreads; ++i) {
        threadData[i].y = y;
        threadData[i].start = i * chunkSize;
        threadData[i].end = (i == numThreads - 1) ? (size - 1) : (threadData[i].start + chunkSize);
        if(pthread_create(&threads[i], NULL, threadFunction1, (void*)&threadData[i])) {
            fprintf(stderr, "Error creating thread\n");
            return 1;
        }
    }

    for (int i = 0; i < numThreads; ++i) {
        pthread_join(threads[i], NULL);
        totalSumCubes += threadData[i].sumCubes;
    }

    double mean = totalSumCubes / (size - 1);
    return mean;
}



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

    double *y1 = malloc((size - tau) * sizeof(double));
    double *y2 = malloc((size - tau) * sizeof(double));

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

