// #include <stdio.h>

// double DN_Mean(const double a[], const int size)
// {
//     double m = 0.0;
//     for (int i = 0; i < size; i++) {
//         m += a[i];
//     }
//     m /= size;
//     return m;
// }


#include <stdio.h>
#include <pthread.h>

#define NUM_THREADS 4 

typedef struct {
    const double* array;
    int start;
    int end;
    double sum;
} ThreadData;

void* partial_sum(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    data->sum = 0.0;
    for (int i = data->start; i < data->end; ++i) {
        data->sum += data->array[i];
    }
    pthread_exit(NULL);
}

double DN_Mean(const double a[], const int size) {
    pthread_t threads[NUM_THREADS];
    ThreadData thread_data[NUM_THREADS];
    int segment_size = size / NUM_THREADS;
    double total_sum = 0.0;

    // Creating threads
    for (int i = 0; i < NUM_THREADS; ++i) {
        thread_data[i].array = a;
        thread_data[i].start = i * segment_size;
        thread_data[i].end = (i == NUM_THREADS - 1) ? size : thread_data[i].start + segment_size;
        pthread_create(&threads[i], NULL, partial_sum, (void*)&thread_data[i]);
    }

    for (int i = 0; i < NUM_THREADS; ++i) {
        pthread_join(threads[i], NULL);
        total_sum += thread_data[i].sum;
    }

    return total_sum / size;
}


