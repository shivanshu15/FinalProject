
#include <pthread.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>

#define MAX_THREADS 16 

typedef struct {
    const double *y;
    int start;
    int end;
    int result; 
} ThreadArg;

// Thread function
void *nan_check_hist(void *arg) {
    ThreadArg *targ = (ThreadArg *)arg;
    for(int i = targ->start; i < targ->end; ++i) {
        if(isnan(targ->y[i])) {
            targ->result = 1;
            return NULL;
        }
    }
    targ->result = 0;
    return NULL;
}

double DN_HistogramMode_5(const double y[], const int size) {
    pthread_t threads[MAX_THREADS];
    ThreadArg args[MAX_THREADS];
    int segment_length = size / MAX_THREADS;

    for(int i = 0; i < MAX_THREADS; ++i) {
        args[i].y = y;
        args[i].start = i * segment_length;
        args[i].end = (i == MAX_THREADS - 1) ? size : args[i].start + segment_length;
        if(pthread_create(&threads[i], NULL, nan_check_hist, &args[i])) {
            fprintf(stderr, "Error creating thread\n");
            return NAN;
        }
    }

    for(int i = 0; i < MAX_THREADS; ++i) {
        pthread_join(threads[i], NULL);
        if(args[i].result) return NAN; 
    }


    const int nBins = 5;
    
    int * histCounts;
    double * binEdges;
    
    histcounts(y, size, nBins, &histCounts, &binEdges);
    
    /*
    for(int i = 0; i < nBins; i++){
        printf("histCounts[%i] = %i\n", i, histCounts[i]);
    }
    for(int i = 0; i < nBins+1; i++){
        printf("binEdges[%i] = %1.3f\n", i, binEdges[i]);
    }
     */
    
    double maxCount = 0;
    int numMaxs = 1;
    double out = 0;;
    for(int i = 0; i < nBins; i++)
    {
        // printf("binInd=%i, binCount=%i, binEdge=%1.3f \n", i, histCounts[i], binEdges[i]);
        
        if (histCounts[i] > maxCount)
        {
            maxCount = histCounts[i];
            numMaxs = 1;
            out = (binEdges[i] + binEdges[i+1])*0.5;
        }
        else if (histCounts[i] == maxCount){
            
            numMaxs += 1;
            out += (binEdges[i] + binEdges[i+1])*0.5;
        }
    }
    out = out/numMaxs;
    
    free(histCounts);
    free(binEdges);
    
    return out;
}

