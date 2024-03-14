#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include "stats.h"

double DN_OutlierInclude_np_001_mdrmd(const double y[], const int size, const int sign)
{
    
    // NaN check
    for(int i = 0; i < size; i++)
    {
        if(isnan(y[i]))
        {
            return NAN;
        }
    }
    
    double inc = 0.01;
    int tot = 0;
    double * yWork = malloc(size * sizeof(double));
    
    // apply sign and check constant time series
    int constantFlag = 1;
    for(int i = 0; i < size; i++)
    {
        if(y[i] != y[0])
        {
            constantFlag = 0;
        }
        
        // apply sign, save in new variable
        yWork[i] = sign*y[i];
        
        // count pos/ negs
        if(yWork[i] >= 0){
            tot += 1;
        }
        
    }
    if(constantFlag) return 0; // if constant, return 0
    
    // find maximum (or minimum, depending on sign)
    double maxVal = max_(yWork, size);
    
    // maximum value too small? return 0
    if(maxVal < inc){
        return 0;
    }
    
    int nThresh = maxVal/inc + 1;
    
    // save the indices where y > threshold
    double * r = malloc(size * sizeof * r);
    
    // save the median over indices with absolute value > threshold
    double * msDti1 = malloc(nThresh * sizeof(double));
    double * msDti3 = malloc(nThresh * sizeof(double));
    double * msDti4 = malloc(nThresh * sizeof(double));
    
    for(int j = 0; j < nThresh; j++)
    {
        //printf("j=%i, thr=%1.3f\n", j, j*inc);
        
        int highSize = 0;
        
        for(int i = 0; i < size; i++)
        {
            if(yWork[i] >= j*inc)
            {
                r[highSize] = i+1;
                //printf("r[%i]=%1.f \n", highSize, r[highSize]);
                highSize += 1;
            }
        }
        
        // intervals between high-values
        double * Dt_exc = malloc(highSize * sizeof(double));
        
        for(int i = 0; i < highSize-1; i++)
        {
            //printf("i=%i, r[i+1]=%1.f, r[i]=%1.f \n", i, r[i+1], r[i]);
            Dt_exc[i] = r[i+1] - r[i];
        }

        /*
        // median
        double medianOut;
        medianOut = median(r, highSize);
        */
         
        msDti1[j] = mean(Dt_exc, highSize-1);
        msDti3[j] = (highSize-1)*100.0/tot;
        msDti4[j] = median(r, highSize) / ((double)size/2) - 1;
        
        //printf("msDti1[%i] = %1.3f, msDti13[%i] = %1.3f, msDti4[%i] = %1.3f\n",
        //       j, msDti1[j], j, msDti3[j], j, msDti4[j]);
        
        free(Dt_exc);
        
    }
    
    int trimthr = 2;
    int mj = 0;
    int fbi = nThresh-1;
    for(int i = 0; i < nThresh; i ++)
    {
        if (msDti3[i] > trimthr)
        {
            mj = i;
        }
        if (isnan(msDti1[nThresh-1-i]))
        {
            fbi = nThresh-1-i;
        }
    }
    
    double outputScalar;
    int trimLimit = mj < fbi ? mj : fbi;
    outputScalar = median(msDti4, trimLimit+1);
    
    free(r);
    free(yWork);
    free(msDti1);
    free(msDti3);
    free(msDti4);
    
    return outputScalar;
}



// #include <math.h>
// #include <pthread.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include "stats.h"

// // Assuming definitions of max_, mean, median functions are available

// typedef struct {
//     const double *yWork;
//     int size;
//     double inc;
//     int startIndex;
//     int endIndex;
//     double *msDti1;
//     double *msDti3;
//     double *msDti4;
//     int tot;
// } ThreadData;

// void *calculateMetrics(void *arg) {
//     ThreadData *data = (ThreadData *)arg;
//     for (int j = data->startIndex; j < data->endIndex; j++) {
//         int highSize = 0;
//         double *r = (double *)malloc(data->size * sizeof(double));
//         for (int i = 0; i < data->size; i++) {
//             if (data->yWork[i] >= j * data->inc) {
//                 r[highSize++] = i + 1;
//             }
//         }
//         if (highSize > 1) { // Only proceed if there are elements to calculate
//             double *Dt_exc = (double *)malloc((highSize - 1) * sizeof(double));
//             for (int i = 0; i < highSize - 1; i++) {
//                 Dt_exc[i] = r[i + 1] - r[i];
//             }
//             data->msDti1[j] = mean(Dt_exc, highSize - 1);
//             data->msDti4[j] = median(r, highSize) / ((double)data->size / 2) - 1;
//             free(Dt_exc);
//         }
//         data->msDti3[j] = (highSize - 1) * 100.0 / data->tot;
//         free(r);
//     }
//     pthread_exit(NULL);
// }

// double DN_OutlierInclude_np_001_mdrmd(const double y[], const int size, const int sign) {
//     // Your existing code for initial checks and setup before parallelization

//     double inc = 0.01;
//     double maxAbs = 0;
//     double * yAbs = malloc(size * sizeof * yAbs);
    
//     for(int i = 0; i < size; i++)
//     {
//         // yAbs[i] = (y[i] > 0) ? y[i] : -y[i];
//         yAbs[i] = (y[i] > 0) ? y[i] : -y[i];
        
//         if(yAbs[i] > maxAbs)
//         {
//             maxAbs = yAbs[i];
//         }
//     }

//     double *yWork = (double *)malloc(size * sizeof(double));
//     int *highInds = (int *)malloc(size * sizeof(int)); // Declaration and allocation for highInds
//     int tot = 0;
//     for (int i = 0; i < size; i++) {
//         yWork[i] = sign * y[i];
//         if (yWork[i] >= 0) tot++;
//     }

//     double maxVal = max_(yWork, size);
//     if (maxVal < 0.01) {
//         free(yWork);
//         free(highInds); // Don't forget to free highInds if returning early
//         return 0;
//     }


//     int nThresh = maxVal / inc + 1;

//     double *msDti1 = malloc(nThresh * sizeof(double));
//     double *msDti3 = malloc(nThresh * sizeof(double));
//     double *msDti4 = malloc(nThresh * sizeof(double));

//     // Thread creation
//     int NUM_THREADS = 4; // Or dynamically adjust based on your system capabilities
//     pthread_t threads[NUM_THREADS];
//     ThreadData threadData[NUM_THREADS];
//     int range = nThresh / NUM_THREADS;

//     for (int i = 0; i < NUM_THREADS; i++) {
//         threadData[i] = (ThreadData){
//             .yWork = yWork,
//             .size = size,
//             .inc = inc,
//             .startIndex = i * range,
//             .endIndex = (i + 1) * range,
//             .msDti1 = msDti1,
//             .msDti3 = msDti3,
//             .msDti4 = msDti4,
//             .tot = tot
//         };

//         if (i == NUM_THREADS - 1) {
//             // Ensure the last thread covers any remaining part due to integer division
//             threadData[i].endIndex = nThresh;
//         }

//         pthread_create(&threads[i], NULL, calculateMetrics, (void *)&threadData[i]);
//     }

//     for (int i = 0; i < NUM_THREADS; i++) {
//         pthread_join(threads[i], NULL);
//     }

//     // Continue with your original code that follows after the parallelized part
//     // This includes finding the final output value after all threads have completed their calculations


//         for(int j = 0; j < nThresh; j++)
//     {
//         int highSize = 0;
        
//         for(int i = 0; i < size; i++)
//         {
//             if(yAbs[i] >= j*inc)
//             {
//                 // fprintf(stdout, "%i, ", i);
                
//                 highInds[highSize] = i;
//                 highSize += 1;
//             }
//         }
        
//         // median
//         double medianOut;
//         medianOut = median(highInds, highSize);
        
//         msDti3[j] = (highSize-1)*100.0/size;
//         msDti4[j] = medianOut / (size/2) - 1;
        
//     }
    
//     int trimthr = 2;
//     int mj = 0;
//     for(int i = 0; i < nThresh; i ++)
//     {
//         if (msDti3[i] > trimthr)
//         {
//             mj = i;
//         }
//     }
    
//     double outputScalar;
//     outputScalar = median(msDti4, mj);

//     // Cleanup
//     free(msDti1);
//     free(msDti3);
//     free(msDti4);
//     // Add any additional cleanup here

//     return outputScalar; // Make sure you compute this as needed after the threads are joined
// }



double DN_OutlierInclude_p_001_mdrmd(const double y[], const int size)
{
    return DN_OutlierInclude_np_001_mdrmd(y, size, 1.0);
}

double DN_OutlierInclude_n_001_mdrmd(const double y[], const int size)
{
    return DN_OutlierInclude_np_001_mdrmd(y, size, -1.0);
}

double DN_OutlierInclude_abs_001(const double y[], const int size)
{
    double inc = 0.01;
    double maxAbs = 0;
    double * yAbs = malloc(size * sizeof * yAbs);
    
    for(int i = 0; i < size; i++)
    {
        // yAbs[i] = (y[i] > 0) ? y[i] : -y[i];
        yAbs[i] = (y[i] > 0) ? y[i] : -y[i];
        
        if(yAbs[i] > maxAbs)
        {
            maxAbs = yAbs[i];
        }
    }
    
    int nThresh = maxAbs/inc + 1;
    
    printf("nThresh = %i\n", nThresh);
    
    double * highInds = malloc(size * sizeof * highInds);
    
    double * msDti3 = malloc(nThresh * sizeof * msDti3);
    double * msDti4 = malloc(nThresh * sizeof * msDti4);

    for(int j = 0; j < nThresh; j++)
    {
        int highSize = 0;
        
        for(int i = 0; i < size; i++)
        {
            if(yAbs[i] >= j*inc)
            {
                // fprintf(stdout, "%i, ", i);
                
                highInds[highSize] = i;
                highSize += 1;
            }
        }
        
        // median
        double medianOut;
        medianOut = median(highInds, highSize);
        
        msDti3[j] = (highSize-1)*100.0/size;
        msDti4[j] = medianOut / (size/2) - 1;
        
    }
    
    int trimthr = 2;
    int mj = 0;
    for(int i = 0; i < nThresh; i ++)
    {
        if (msDti3[i] > trimthr)
        {
            mj = i;
        }
    }
    
    double outputScalar;
    outputScalar = median(msDti4, mj);

    free(highInds);
    free(yAbs);
    free(msDti4);
    
    return outputScalar;
}
