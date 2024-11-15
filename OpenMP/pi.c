#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <time.h>
#include <stdlib.h>

double pi(double N){
    double sum = 0;
    for(int i=0; i<=N; i++){
            sum += 4/(1+(i/N)*(i/N));
    }
    sum /= N;
    return sum;
}

double pi_OpenMP1(double N){
    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for(int i=0; i <= (int) N; i++){
            sum += 4/(1+(i/N)*(i/N));
    }
    sum /= N;
    return sum;
}

double pi_OpenMP2(double N){
    
    double* values = (double*) malloc(((int)N+1)*sizeof(double));
    #pragma omp parallel for 
    for(int i=0; i <= (int) N; i++){
            values[i] = 4/(1+(i/N)*(i/N));
    }
    
    double sum = 0;
    for(int i=0; i <= (int) N; i++){
            sum += values[i]; 
    }
    sum /= N;
    free(values);
    return sum;
}

double pi_OpenMP3(double N){
    
    double total_Sum = 0;

    #pragma omp parallel 
    {
        double partial_Sum = 0;

        #pragma omp for
        for(int i = 0; i <= (int) N; i++){
            partial_Sum += 4/(1+(i/N)*(i/N));
        }

        //Create thread safe region.
        #pragma omp critical
        {
            //add each threads partial sum to the total sum
            total_Sum += partial_Sum;
        }
    }
    return total_Sum/N;
}

double pi_OpenMP4(double N){
    
    double total_Sum = 0;

    #pragma omp parallel 
    {
        double partial_Sum = 0;

        #pragma omp for
        for(int i = 0; i <= (int) N; i++){
            partial_Sum += 4/(1+(i/N)*(i/N));
        }

        //Create thread safe region.
        #pragma omp atomic update
        total_Sum += partial_Sum;
    }
    return total_Sum/N;
}


int main(){

    const double N = 1000000000.;

    clock_t start, end, time;

    double start_openmp, end_openmp, time_openmp;

    printf("# OpenMP 1 # \n");
    start_openmp = omp_get_wtime();
    double p = pi_OpenMP1(N);
    end_openmp = omp_get_wtime();
    time_openmp = end_openmp - start_openmp;
    printf("- result: %lf\n",p);
    printf("- time: %lf\n",time_openmp);

    printf("# OpenMP 2 # \n");
    start_openmp = omp_get_wtime();
    p = pi_OpenMP2(N);
    end_openmp = omp_get_wtime();
    time_openmp = end_openmp - start_openmp;
    printf("- result: %lf\n",p);
    printf("- time: %lf\n",time_openmp);

    printf("# OpenMP 3 # \n");
    start_openmp = omp_get_wtime();
    p = pi_OpenMP3(N);
    end_openmp = omp_get_wtime();
    time_openmp = end_openmp - start_openmp;
    printf("- result: %lf\n",p);
    printf("- time: %lf\n",time_openmp);

    printf("# OpenMP 4 # \n");
    start_openmp = omp_get_wtime();
    p = pi_OpenMP4(N);
    end_openmp = omp_get_wtime();
    time_openmp = end_openmp - start_openmp;
    printf("- result: %lf\n",p);
    printf("- time: %lf\n",time_openmp);

    printf("# SEQUENTIAL # \n");
    start = clock();
    p = pi(N);
    end = clock();
    time = end - start;
    printf("- result: %lf\n",p);
    printf("- time: %lf\n",((double) (end - start)) / CLOCKS_PER_SEC);



}