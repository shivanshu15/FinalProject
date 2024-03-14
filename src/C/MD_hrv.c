//
//  MD_hrv.c
//  C_polished
//
//  Created by Carl Henning Lubba on 22/09/2018.
//  Copyright © 2018 Carl Henning Lubba. All rights reserved.
//

// #include "MD_hrv.h"
// #include "stats.h"

// double MD_hrv_classic_pnn40(const double y[], const int size){
    
//     // NaN check
//     for(int i = 0; i < size; i++)
//     {
//         if(isnan(y[i]))
//         {
//             return NAN;
//         }
//     }
    
//     const int pNNx = 40;
    
//     // compute diff
//     double * Dy = malloc((size-1) * sizeof(double));
//     diff(y, size, Dy);
    
//     double pnn40 = 0;
//     for(int i = 0; i < size-1; i++){
//         if(fabs(Dy[i])*1000 > pNNx){
//             pnn40 += 1;
//         }
//     }
    
//     free(Dy);
    
//     return pnn40/(size-1);
// }




#include <pthread.h>
#include <stdlib.h>
#include <math.h>
#include "MD_hrv.h"
#include "stats.h"

typedef struct {
    const double* y;
    int start;
    int end;
    double pnn40Partial;
} ThreadArg;

void* calculatePartialPnn40(void* arg) {
    ThreadArg* data = (ThreadArg*) arg;
    const double pNNx = 40;
    double pnn40 = 0;
    for(int i = data->start; i < data->end; i++) {
        if(fabs(data->y[i])*1000 > pNNx){
            pnn40 += 1;
        }
    }
    data->pnn40Partial = pnn40;
    pthread_exit(NULL);
}


#include <stdio.h>
#include <pthread.h>

#define NUM_THREADS 8 // Adjust based on your needs and system capabilities

double MD_hrv_classic_pnn40(const double y[], const int size) {
    double *Dy = malloc((size - 1) * sizeof(double));
    diff(y, size, Dy);

    pthread_t threads[NUM_THREADS];
    ThreadArg args[NUM_THREADS];
    int segmentSize = (size - 1) / NUM_THREADS;

    for (int i = 0; i < NUM_THREADS; i++) {
        args[i].y = Dy;
        args[i].start = i * segmentSize;
        args[i].end = (i == NUM_THREADS - 1) ? (size - 1) : (i + 1) * segmentSize;
        pthread_create(&threads[i], NULL, calculatePartialPnn40, (void*)&args[i]);
    }

    double totalPnn40 = 0;
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
        totalPnn40 += args[i].pnn40Partial;
    }

    free(Dy);
    
    return totalPnn40 / (size - 1);
}
