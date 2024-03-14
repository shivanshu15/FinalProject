// #include <stdio.h>
// #include "stats.h"

// double DN_Spread_Std(const double a[], const int size)
// {
//     double m = mean(a, size);
//     double sd = 0.0;
//     for (int i = 0; i < size; i++) {
//         sd += pow(a[i] - m, 2);
//     }
//     sd = sqrt(sd / (size - 1));
//     return sd;
// }


#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include "stats.h"

#define NUM_THREADS 4 // Adjust based on your requirements

typedef struct {
    const double* array;
    double mean;
    int start;
    int end;
    double sum;
} ThreadData;

void* partial_sum_std(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    data->sum = 0.0;
    for (int i = data->start; i < data->end; ++i) {
        data->sum += pow(data->array[i] - data->mean, 2);
    }
    pthread_exit(NULL);
}

double DN_Spread_Std(const double a[], const int size) {
    pthread_t threads[NUM_THREADS];
    ThreadData thread_data[NUM_THREADS];
    int segment_size = size / NUM_THREADS;
    double total_sum = 0.0;
    double m = mean(a, size);

    // Creating threads for standard deviation computation
    for (int i = 0; i < NUM_THREADS; ++i) {
        thread_data[i].array = a;
        thread_data[i].mean = m;
        thread_data[i].start = i * segment_size;
        thread_data[i].end = (i == NUM_THREADS - 1) ? size : thread_data[i].start + segment_size;
        pthread_create(&threads[i], NULL, partial_sum_std, (void*)&thread_data[i]);
    }

    // Joining threads and accumulating results
    for (int i = 0; i < NUM_THREADS; ++i) {
        pthread_join(threads[i], NULL);
        total_sum += thread_data[i].sum;
    }

    return sqrt(total_sum / (size - 1));
}

