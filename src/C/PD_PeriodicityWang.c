

#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include "PD_PeriodicityWang.h"
#include "splinefit.h"
#include "stats.h"

typedef struct {
    const double *y;
    int start;
    int end;
    int result; // 1 for no NaN found, 0 for NaN found
} NanCheckArgs;

typedef struct {
    const double *ySub;
    double *acf;
    int size;
    int startLag;
    int endLag;
} AutocorrArgs;



void* autocorrelation(void* arguments) {
    AutocorrArgs* args = (AutocorrArgs*)arguments;
    for (int tau = args->startLag; tau <= args->endLag; tau++) {
        args->acf[tau-1] = autocov_lag(args->ySub, args->size, tau);
    }
    pthread_exit(NULL);
}

int PD_PeriodicityWang_th0_01(const double * y, const int size) {
    int num_threads = 8;
    pthread_t threads[num_threads];
    NanCheckArgs nan_args[num_threads];
    int chunk_size = size / num_threads;


    const double th = 0.01;
    
    double * ySpline = malloc(size * sizeof(double));
    
    splinefit(y, size, ySpline);
    
    //printf("spline fit complete.\n");
    
    // subtract spline from data to remove trend
    double * ySub = malloc(size * sizeof(double));
    for(int i = 0; i < size; i++){
        ySub[i] = y[i] - ySpline[i];
        //printf("ySub[%i] = %1.5f\n", i, ySub[i]);
    }
    
    int acmax = (int)ceil((double)size/3);
    
    double * acf = malloc(acmax*sizeof(double));
    for(int tau = 1; tau <= acmax; tau++){
        acf[tau-1] = autocov_lag(ySub, size, tau);
        //printf("acf[%i] = %1.9f\n", tau-1, acf[tau-1]);
    }
    
    //printf("ACF computed.\n");
    
    // find troughts and peaks
   


    AutocorrArgs autocorr_args[num_threads];
    //double *acf = malloc((size / 3) * sizeof(double)); // Example, adjust as needed

    for (int i = 0; i < num_threads; i++) {
        autocorr_args[i].ySub = y; // Assuming ySub is defined after spline fitting and subtraction
        autocorr_args[i].acf = acf;
        autocorr_args[i].size = size;
        autocorr_args[i].startLag = i * (size / 3) / num_threads + 1;
        autocorr_args[i].endLag = (i + 1) * (size / 3) / num_threads;
        pthread_create(&threads[i], NULL, autocorrelation, (void*)&autocorr_args[i]);
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }



     double * troughs = malloc(acmax * sizeof(double));
    double * peaks = malloc(acmax * sizeof(double));
    int nTroughs = 0;
    int nPeaks = 0;
    double slopeIn = 0;
    double slopeOut = 0;
    for(int i = 1; i < acmax-1; i ++){
        slopeIn = acf[i] - acf[i-1];
        slopeOut = acf[i+1] - acf[i];
        
        if(slopeIn < 0 & slopeOut > 0)
        {
            // printf("trough at %i\n", i);
            troughs[nTroughs] = i;
            nTroughs += 1;
        }
        else if(slopeIn > 0 & slopeOut < 0)
        {
            // printf("peak at %i\n", i);
            peaks[nPeaks] = i;
            nPeaks += 1;
        }
    }
    
    //printf("%i troughs and %i peaks found.\n", nTroughs, nPeaks);
    
    
    // search through all peaks for one that meets the conditions:
    // (a) a trough before it
    // (b) difference between peak and trough is at least 0.01
    // (c) peak corresponds to positive correlation
    int iPeak = 0;
    double thePeak = 0;
    int iTrough = 0;
    double theTrough = 0;
    
    int out = 0;
    
    for(int i = 0; i < nPeaks; i++){
        iPeak = peaks[i];
        thePeak = acf[iPeak];
        
        //printf("i=%i/%i, iPeak=%i, thePeak=%1.3f\n", i, nPeaks-1, iPeak, thePeak);
        
        // find trough before this peak
        int j = -1;
        while(troughs[j+1] < iPeak && j+1 < nTroughs){
            // printf("j=%i/%i, iTrough=%i, theTrough=%1.3f\n", j+1, nTroughs-1, (int)troughs[j+1], acf[(int)troughs[j+1]]);
            j++;
        }
        if(j == -1)
            continue;
        
        iTrough = troughs[j];
        theTrough = acf[iTrough];
        
        // (a) should be implicit
        
        // (b) different between peak and trough it as least 0.01
        if(thePeak - theTrough < th)
            continue;
        
        if(thePeak < 0)
            continue;
        
        out = iPeak;
        break;
    }
    

    free(acf);

    free(ySpline);
    free(ySub);
    //free(acf);
    free(troughs);
    free(peaks);
    
    return out;

}

