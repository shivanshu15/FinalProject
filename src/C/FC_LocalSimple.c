#include <math.h>
#include <string.h>
#include "stats.h"
#include "CO_AutoCorr.h"

static void abs_diff(const double a[], const int size, double b[])
{
    for (int i = 1; i < size; i++) {
        b[i - 1] = fabs(a[i] - a[i - 1]);
    }
}

double fc_local_simple(const double y[], const int size, const int train_length)
{
    double * y1 = malloc((size - 1) * sizeof *y1);
    abs_diff(y, size, y1);
    double m = mean(y1, size - 1);
    free(y1);
    return m;
}






#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct {
    const double *y;
    double *res;
    int startIdx;
    int endIdx;
    int train_length;
} ThreadArg0;

void* processPortion(void *arg) {
    ThreadArg0 *targ = (ThreadArg0*)arg;
    for (int i = targ->startIdx; i < targ->endIdx; i++) {
        double yest = 0;
        for (int j = 0; j < targ->train_length; j++) {
            yest += targ->y[i + j];
        }
        yest /= targ->train_length;
        targ->res[i] = targ->y[i + targ->train_length] - yest;
    }
    return NULL;
}

double FC_LocalSimple_mean_tauresrat(const double y[], const int size, const int train_length, const double *autocorrs) {
    int num_threads = 8; 
    pthread_t threads[num_threads];
    ThreadArg0 args[num_threads];

    double *res = malloc((size - train_length) * sizeof(*res));
    int chunk_size = (size - train_length) / num_threads;

    for (int i = 0; i < num_threads; i++) {
        args[i].y = y;
        args[i].res = res;
        args[i].startIdx = i * chunk_size;
        args[i].endIdx = (i == num_threads - 1) ? (size - train_length) : (i + 1) * chunk_size;
        args[i].train_length = train_length;
        pthread_create(&threads[i], NULL, processPortion, &args[i]);
    }

    // Join threads
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    double resAC1stZ = co_firstzero(res, size - train_length, size - train_length, autocorrs);
    double yAC1stZ = co_firstzero(y, size, size, autocorrs);
    double output = resAC1stZ / yAC1stZ;

    free(res);
    return output;
}








#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct {
    const double* y;
    double* res;
    int start;
    int end;
    int train_length;
} ThreadArg;



void* thread_func(void* arg) {
    ThreadArg* t_arg = (ThreadArg*)arg;
    const double* y = t_arg->y;
    double* res = t_arg->res;
    int start = t_arg->start;
    int end = t_arg->end;
    int train_length = t_arg->train_length;
    
    for (int i = start; i < end; i++) {
        double yest = 0;
        for (int j = 0; j < train_length; j++) {
            yest += y[i+j];
        }
        yest /= train_length;
        res[i] = y[i+train_length] - yest;
    }
    return NULL;
}

double FC_LocalSimple_mean_stderr(const double y[], const int size, const int train_length) {
    int num_threads = 8; 
    pthread_t threads[num_threads];
    ThreadArg args[num_threads];
    
    double *res = malloc((size - train_length) * sizeof(*res));
    if (res == NULL) {
        fprintf(stderr, "Failed to allocate memory\n");
        return -1; 
    }
    
    int chunk_size = (size - train_length) / num_threads;
    
    for (int i = 0; i < num_threads; i++) {
        args[i].y = y;
        args[i].res = res;
        args[i].start = i * chunk_size;
        args[i].end = (i == num_threads - 1) ? (size - train_length) : (i + 1) * chunk_size;
        args[i].train_length = train_length;
        
        if (pthread_create(&threads[i], NULL, thread_func, &args[i])) {
            fprintf(stderr, "Error creating thread\n");
            free(res);
            return -1; 
        }
    }
    
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }
    
    double output = stddev(res, size - train_length);
    
    free(res);
    return output;
}









double FC_LocalSimple_mean3_stderr(const double y[], const int size)
{
    return FC_LocalSimple_mean_stderr(y, size, 3);
}

double FC_LocalSimple_mean1_tauresrat(const double y[], const int size,const double * autocorrs){
    return FC_LocalSimple_mean_tauresrat(y, size, 1, autocorrs);
}

double FC_LocalSimple_mean_taures(const double y[], const int size, const int train_length,const double * autocorrs)
{
    double * res = malloc((size - train_length) * sizeof *res);
    
    // first z-score
    // no, assume ts is z-scored!!
    //zscore_norm(y, size);
    
    for (int i = 0; i < size - train_length; i++)
    {
        double yest = 0;
        for (int j = 0; j < train_length; j++)
        {
            yest += y[i+j];
            
        }
        yest /= train_length;
        
        res[i] = y[i+train_length] - yest;
    }
    
    int output = co_firstzero(res, size - train_length, size - train_length, autocorrs);
    
    free(res);
    return output;
    
}

double FC_LocalSimple_lfit_taures(const double y[], const int size, const double * autocorrs)
{
    // set tau from first AC zero crossing
    int train_length = co_firstzero(y, size, size,autocorrs);
    
    double * xReg = malloc(train_length * sizeof * xReg);
    // double * yReg = malloc(train_length * sizeof * yReg);
    for(int i = 1; i < train_length+1; i++)
    {
        xReg[i-1] = i;
    }
    
    double * res = malloc((size - train_length) * sizeof *res);
    
    double m = 0.0, b = 0.0;
    
    for (int i = 0; i < size - train_length; i++)
    {
        linreg(train_length, xReg, y+i, &m, &b);
        
        // fprintf(stdout, "i=%i, m=%f, b=%f\n", i, m, b);
        
        res[i] = y[i+train_length] - (m * (train_length+1) + b);
    }
    
    int output = co_firstzero(res, size - train_length, size - train_length,autocorrs);
    
    free(res);
    free(xReg);
    // free(yReg);
    
    return output;
    
}


