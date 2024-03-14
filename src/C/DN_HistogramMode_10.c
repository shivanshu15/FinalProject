#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "stats.h"
#include "histcounts.h"

// double DN_HistogramMode_10(const double y[], const int size)
// {
//     // NaN check
//     for(int i = 0; i < size; i++)
//     {
//         if(isnan(y[i]))
//         {
//             return NAN;
//         }
//     }
    
//     const int nBins = 10;
    
//     int * histCounts;
//     double * binEdges;
    
//     histcounts(y, size, nBins, &histCounts, &binEdges);
    
//     double maxCount = 0;
//     int numMaxs = 1;
//     double out = 0;;
//     for(int i = 0; i < nBins; i++)
//     {
//         // printf("binInd=%i, binCount=%i, binEdge=%1.3f \n", i, histCounts[i], binEdges[i]);
        
//         if (histCounts[i] > maxCount)
//         {
//             maxCount = histCounts[i];
//             numMaxs = 1;
//             out = (binEdges[i] + binEdges[i+1])*0.5;
//         }
//         else if (histCounts[i] == maxCount){
            
//             numMaxs += 1;
//             out += (binEdges[i] + binEdges[i+1])*0.5;
//         }
//     }
//     out = out/numMaxs;
    
//     // arrays created dynamically in function histcounts
//     free(histCounts);
//     free(binEdges);
    
//     return out;
// }





#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>

#define NUM_THREADS 8 



typedef struct {
    const double* y;
    int start;
    int end;
    double result;
} ThreadData;


void* checkNaN(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    for(int i = data->start; i < data->end; i++) {
        if(isnan(data->y[i])) {
            data->result = NAN;
            return NULL;
        }
    }
    data->result = 0; 
    return NULL;
}


double DN_HistogramMode_10(const double y[], const int size) {
    pthread_t threads[NUM_THREADS];
    ThreadData threadData[NUM_THREADS];
    int segmentLength = size / NUM_THREADS;
    
    for(int i = 0; i < NUM_THREADS; i++) {
        threadData[i].y = y;
        threadData[i].start = i * segmentLength;
        threadData[i].end = (i == NUM_THREADS - 1) ? size : (i + 1) * segmentLength;
        pthread_create(&threads[i], NULL, checkNaN, (void*)&threadData[i]);
    }
    
    for(int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
        if(threadData[i].result == NAN) {
            return NAN;
        }
    }
    
    const int nBins = 10;
    int *histCounts;
    double *binEdges;
    histcounts(y, size, nBins, &histCounts, &binEdges);
    
    double maxCount = 0;
    int numMaxs = 1;
    double out = 0;
    for(int i = 0; i < nBins; i++) {
        if (histCounts[i] > maxCount) {
            maxCount = histCounts[i];
            numMaxs = 1;
            out = (binEdges[i] + binEdges[i+1]) * 0.5;
        } else if (histCounts[i] == maxCount) {
            numMaxs += 1;
            out += (binEdges[i] + binEdges[i+1]) * 0.5;
        }
    }
    out = out / numMaxs;
    
    free(histCounts);
    free(binEdges);
    
    return out;
}









/*
 double DN_HistogramMode_10(double y[], int size)
 {
 
 double min = DBL_MAX, max=-DBL_MAX;
 for(int i = 0; i < size; i++)
 {
 if (y[i] < min)
 {
 min = y[i];
 }
 if (y[i] > max)
 {
 max = y[i];
 }
 }
 
 double binStep = (max - min)/10;
 
 // fprintf(stdout, "min=%f, max=%f, binStep=%f \n", min, max, binStep);
 
 int histCounts[10] = {0};
 for(int i = 0; i < size; i++)
 {
 int binsLeft = 10;
 int lowerInd = 0, upperInd = 10;
 while(binsLeft > 1)
 {
 int limitInd = (upperInd - lowerInd)/2 + lowerInd;
 double limit = limitInd * binStep + min;
 
 if (y[i] < limit)
 {
 upperInd = limitInd;
 }
 else
 {
 lowerInd = limitInd;
 }
 binsLeft = upperInd - lowerInd;
 }
 histCounts[lowerInd] += 1;
 }
 
 double maxCount = 0;
 int maxCountInd = 0;
 for(int i = 0; i < 10; i++)
 {
 // fprintf(stdout, "binInd=%i, binCount=%i \n", i, histCounts[i]);
 
 if (histCounts[i] > maxCount)
 {
 maxCountInd = i;
 maxCount = histCounts[i];
 }
 }
 return binStep*(maxCountInd+0.5) + min;
 }
 */
