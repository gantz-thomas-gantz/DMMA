#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

void cs(void** A, int N){
    void** B = (void**) malloc(N*sizeof(void*));
    int next = 0;
    for(int i=0; i<N; i++){
        if(A[i]!=NULL){
            B[next] = A[i];
            next++;
        }        
    }
    free(B);
}

void cs_OpenMP1(void** A, int N){
    
    void** B = (void**) malloc(N*sizeof(void*));
    int next = 0;

    #pragma omp parallel 
    {   
        int this;

        #pragma omp for
        for(int i=0; i<N; i++){
            if(A[i]!=NULL){
                #pragma omp critical
                {
                    B[next] = A[i];
                    next++;
                }
                
                
            }        
        }
    }
    free(B);
}

void cs_OpenMP2(void** A, int N){
    
    void** B = (void**) malloc(N*sizeof(void*));
    int next = 0;

    #pragma omp parallel 
    {   
        int this;

        #pragma omp for
        for(int i=0; i<N; i++){
            if(A[i]!=NULL){
                #pragma omp atomic capture
                {
                    this = next; 
                    next++;
                }
                
                B[this] = A[i];
            }        
        }
    }
    free(B);
}

// keep the order
void cs_OpenMP3(void** A, int N){
    
    void** B = (void**) malloc(N*sizeof(void*));
    int threads;
    int* nexts = (int*) calloc(100,sizeof(int)); // Assumption: not more than 100 threads

    #pragma omp parallel 
    {   
        
        int rank = omp_get_thread_num();
        threads = omp_get_num_threads();

        #pragma omp for schedule(static)
        for(int i=0; i<N; i++){
            if(A[i]!=NULL){
                B[rank*(N/threads)+nexts[rank]] = A[i];
                nexts[rank]++;
            }        
        }
    }
 
    void** C = (void**) malloc(N*sizeof(void*));
    int counter = 0;
    for(int rank=0; rank<threads;rank++){
        for(int k=rank*(N/threads); k<nexts[rank]; k++){
            C[counter] = B[k];
            counter++;
        }
    }
    free(B);
    free(C);
    free(nexts);
}


int main(){

    int x = 3;
    const int N = 1e8; 
    void** A = (void**) malloc(N*sizeof(void*));
    for(int i=0; i<N; i++){
        if((i%10)==0){
            A[i] = NULL;
        }  
        else{
            A[i] = (void*) &x;
        }      
    }

    clock_t start, end, time;

    double start_openmp, end_openmp, time_openmp;

    printf("# OpenMP 1 # \n");
    start_openmp = omp_get_wtime();
    cs_OpenMP1(A,N);
    end_openmp = omp_get_wtime();
    time_openmp = end_openmp - start_openmp;
    printf("- time: %lf\n",time_openmp);

    printf("# OpenMP 2 # \n");
    start_openmp = omp_get_wtime();
    cs_OpenMP2(A,N);
    end_openmp = omp_get_wtime();
    time_openmp = end_openmp - start_openmp;
    printf("- time: %lf\n",time_openmp);

    printf("# OpenMP 3 # \n");
    start_openmp = omp_get_wtime();
    cs_OpenMP3(A,N);
    end_openmp = omp_get_wtime();
    time_openmp = end_openmp - start_openmp;
    printf("- time: %lf\n",time_openmp);

    printf("# SEQUENTIAL # \n");
    start_openmp = omp_get_wtime();
    cs(A,N);
    end_openmp = omp_get_wtime();
    time_openmp = end_openmp - start_openmp;
    printf("- time: %lf\n",time_openmp);

}