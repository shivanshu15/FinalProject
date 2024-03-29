

#include <math.h>
#include <stdlib.h>
#include <pthread.h>

#include "IN_AutoMutualInfoStats.h"
#include "CO_AutoCorr.h"
#include "stats.h"

typedef struct {
    const double *y;
    int size;
    int start;
    int end;
    double *ami;
} ThreadData;

void *compute_ami(void *args) {
    ThreadData *data = (ThreadData *)args;
    for (int i = data->start; i < data->end; i++) {
        double ac = autocorr_lag(data->y, data->size, i + 1);
        data->ami[i] = -0.5 * log(1 - ac * ac);
    }
    pthread_exit(NULL);
}

double IN_AutoMutualInfoStats_40_gaussian_fmmi(const double y[], const int size) {
    int tau = 40;
    if (tau > ceil((double)size / 2)) {
        tau = ceil((double)size / 2);
    }

    double *ami = malloc(size * sizeof(double));
    const int numThreads = 8;
    pthread_t threads[numThreads];
    ThreadData threadData[numThreads];
    int chunkSize = tau / numThreads;

    for (int i = 0; i < numThreads; i++) {
        threadData[i].y = y;
        threadData[i].size = size;
        threadData[i].ami = ami;
        threadData[i].start = i * chunkSize;
        threadData[i].end = (i == numThreads - 1) ? tau : (i + 1) * chunkSize;
        pthread_create(&threads[i], NULL, compute_ami, (void *)&threadData[i]);
    }

    for (int i = 0; i < numThreads; i++) {
        pthread_join(threads[i], NULL);
    }

    double fmmi = tau;
    for (int i = 1; i < tau - 1; i++) {
        if (ami[i] < ami[i - 1] && ami[i] < ami[i + 1]) {
            fmmi = i;
            break;
        }
    }

    free(ami);
    return fmmi;
}





//
//  IN_AutoMutualInfoStats.c
//  C_polished
//
//  Created by Carl Henning Lubba on 22/09/2018.
//  Copyright © 2018 Carl Henning Lubba. All rights reserved.
//
// #include <math.h>

// #include "IN_AutoMutualInfoStats.h"
// #include "CO_AutoCorr.h"
// #include "stats.h"

// double IN_AutoMutualInfoStats_40_gaussian_fmmi(const double y[], const int size)
// {
//     // NaN check
//     for(int i = 0; i < size; i++)
//     {
//         if(isnan(y[i]))
//         {
//             return NAN;
//         }
//     }
    
//     // maximum time delay
//     int tau = 40;
    
//     // don't go above half the signal length
//     if(tau > ceil((double)size/2)){
//         tau = ceil((double)size/2);
//     }
    
//     // compute autocorrelations and compute automutual information
//     double * ami = malloc(size * sizeof(double));
//     for(int i = 0; i < tau; i++){
//         double ac = autocorr_lag(y,size, i+1);
//         ami[i] = -0.5 * log(1 - ac*ac);
//         // printf("ami[%i]=%1.7f\n", i, ami[i]);
//     }
    
//     // find first minimum of automutual information
//     double fmmi = tau;
//     for(int i = 1; i < tau-1; i++){
//         if(ami[i] < ami[i-1] & ami[i] < ami[i+1]){
//             fmmi = i;
//             // printf("found minimum at %i\n", i);
//             break;
//         }
//     }
    
//     free(ami);
    
//     return fmmi;
// }