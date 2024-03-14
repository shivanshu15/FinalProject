//
//  SB_BinaryStats.c
//  C_polished
//
//  Created by Carl Henning Lubba on 22/09/2018.
//  Copyright © 2018 Carl Henning Lubba. All rights reserved.
//

#include "SB_BinaryStats.h"
#include "stats.h"
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>




typedef struct {
    const double* y;
    int* yBin;
    int start;
    int end;
} ThreadArgs;


void* threadFunc(void* arg) {
    ThreadArgs* args = (ThreadArgs*)arg;
    for (int i = args->start; i < args->end; i++) {
        
        // Binarize, skipping the last element as it does not have a pair
        if (i < args->end - 1) {
            double diffTemp = args->y[i + 1] - args->y[i];
            args->yBin[i] = diffTemp < 0 ? 0 : 1;
        }
    }
    pthread_exit(NULL);
}


#define NUM_THREADS 8

double SB_BinaryStats_diff_longstretch0(const double y[], const int size) {
    pthread_t threads[NUM_THREADS];
    ThreadArgs args[NUM_THREADS];
    int segmentSize = (size - 1) / NUM_THREADS; 

    int* yBin = malloc((size - 1) * sizeof(int));
    if (yBin == NULL) {
        fprintf(stderr, "Failed to allocate memory.\n");
        return -1; 
    }

    for (int i = 0; i < NUM_THREADS; i++) {
        args[i].y = y;
        args[i].yBin = yBin;
        args[i].start = i * segmentSize;
        args[i].end = (i == NUM_THREADS - 1) ? size : args[i].start + segmentSize;
        int rc = pthread_create(&threads[i], NULL, threadFunc, (void*)&args[i]);
        if (rc) {
            fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

   for (int i = 0; i < NUM_THREADS; i++) {
    void* status;
    pthread_join(threads[i], &status);
    if (status != NULL) { 
        free(yBin);
        return *((double*)status);
    }
}

    int maxstretch0 = 0;
    int last1 = 0;
    for (int i = 0; i < size - 1; i++) {
        if (yBin[i] == 1 || i == size - 2) {
            double stretch0 = i - last1;
            if (stretch0 > maxstretch0) {
                maxstretch0 = stretch0;
            }
            last1 = i;
        }
    }

    free(yBin);

    return maxstretch0;
}




#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

typedef struct {
    const double* y;
    double yMean;
    int start;
    int end;
    int maxstretch1;
} ThreadData;

double mean(const double* y, int size);

void* processStretch(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    int* yBin = malloc((data->end - data->start) * sizeof(int));

    for (int i = data->start; i < data->end; i++) {
        yBin[i - data->start] = (data->y[i] - data->yMean <= 0) ? 0 : 1;
    }

    int maxstretch1 = 0;
    int last1 = -1; 
    for (int i = 0; i < data->end - data->start; i++) {
        if (yBin[i] == 0 || i == data->end - data->start - 1) {
            int stretch1 = i - last1;
            if (stretch1 > maxstretch1) {
                maxstretch1 = stretch1;
            }
            last1 = i;
        }
    }

    data->maxstretch1 = maxstretch1;
    free(yBin);

    pthread_exit(NULL);
}

double SB_BinaryStats_mean_longstretch1(const double y[], const int size) {
    double yMean = mean(y, size);

    int numThreads = 8;
    pthread_t threads[numThreads];
    ThreadData data[numThreads];
    int segmentLength = (size - 1) / numThreads;

    for (int i = 0; i < numThreads; i++) {
        data[i].y = y;
        data[i].yMean = yMean;
        data[i].start = i * segmentLength;
        data[i].end = (i == numThreads - 1) ? size - 1 : (i + 1) * segmentLength;
        pthread_create(&threads[i], NULL, processStretch, (void*)&data[i]);
    }

    int maxstretch1 = 0;
    for (int i = 0; i < numThreads; i++) {
        pthread_join(threads[i], NULL);
        if (data[i].maxstretch1 > maxstretch1) {
            maxstretch1 = data[i].maxstretch1;
        }
    }

    return maxstretch1;
}



