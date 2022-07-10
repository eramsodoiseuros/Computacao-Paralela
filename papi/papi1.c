#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "papi.h"

//PAPI events to monitor
#define NUM_EVENTS 4
// PAPI counters' values
long long values[NUM_EVENTS], min_values[NUM_EVENTS];
int retval, EventSet=PAPI_NULL;
int Events[NUM_EVENTS] = { PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_L1_DCM, PAPI_L2_DCM  };


//Number of times the function is executed and measured
#define NUM_RUNS 1

//Estrutura do Bucket
typedef struct bucket {
    int n_elem;
    int start;
    int index;
}*Bucket;

void valid(int* array, int array_size){
    int valid = 1, i;
    for (i=0; i < array_size-1; i++) {
        if (array[i] > array[i+1]){
            valid = 0;
            printf("\n\tfalhou na posicao: %d, <%d> <!%d!> <%d>", i, array[i-1], array[i], array[i+1]);
            break;
        }
    }

    if (valid) printf("\tThe array was correctly sorted.\n");
    else printf("\tThe array was not correctly sorted.\n");
}

int cmpfunc (const void * a, const void * b){
    return ( *(int*)a - *(int*)b );
}

void swap(int *a, int *b) {
  int t = *a;
  *a = *b;
  *b = t;
}

int partition(int array[], int low, int high) {
  int pivot = array[high];
  int i = (low - 1);
  int j;
  for (j = low; j < high; j++) {
    if (array[j] <= pivot) {
      i++;
      swap(&array[i], &array[j]);
    }
  }

  swap(&array[i + 1], &array[high]);

  return (i + 1);
}

void quickSort(int array[], int low, int high) {
  if (low < high) {
    int pi = partition(array, low, high);
    quickSort(array, low, pi - 1);
    quickSort(array, pi + 1, high);
  }
}

void bucketSortParalelo(int array[], int array_baldes[], int array_size, int width, int n_baldes, struct bucket* baldes) {
    int num_threads = n_baldes;
    omp_set_num_threads(num_threads);

    int i, b_index;
    int global_n_elem[n_baldes];
    int global_starting_position[n_baldes];
    memset(global_n_elem, 0, sizeof(int) * n_baldes);
    memset(global_starting_position, 0, sizeof(int) * n_baldes);

    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
        int local, real;
        int my_id = omp_get_thread_num();


        #pragma omp for private(i, local)
        for (i = 0; i < array_size; i++){
            local = array[i]/width;
            if (local > n_baldes-1)
                local = n_baldes-1;

            real = local + my_id * n_baldes;
            baldes[real].n_elem++;
        }

        int local_sum = 0;
        for (i = my_id; i < n_baldes*num_threads; i = i + num_threads){
            local_sum += baldes[i].n_elem;
        }

        global_n_elem[my_id] = local_sum;

        #pragma omp barrier

        #pragma omp single
        {
            for (i = 1; i < n_baldes; i++){
                global_starting_position[i] = global_starting_position[i-1] + global_n_elem[i-1];
                baldes[i].index = baldes[i-1].index + global_n_elem[i-1];
                baldes[i].start = baldes[i-1].start + global_n_elem[i-1];
            }
        }

        for (i = my_id + n_baldes; i < n_baldes * num_threads; i = i + num_threads){
            int previous_index = i - n_baldes;
            baldes[i].start = baldes[previous_index].start + baldes[previous_index].n_elem;
            baldes[i].index = baldes[previous_index].index + baldes[previous_index].n_elem;  
        }

        #pragma omp barrier

        #pragma omp for private(i, b_index)
        for (i = 0; i < array_size; i++){
            local = array[i]/width;
            if (local > n_baldes - 1)
                local = n_baldes - 1;

            real = local + my_id * n_baldes;
            b_index = baldes[real].index++;
            array_baldes[b_index] = array[i];
        }


        #pragma omp for private(i)
        for (i = 0; i < n_baldes; i++){
            qsort(array_baldes+global_starting_position[i], global_n_elem[i], sizeof(int), cmpfunc);
        }

    }
}

void papi_paralelo(int tamanho, int buckets){
    int max = 99999;
    int *array, * array_baldes;
    int width = (int)max / buckets;
    array = (int*)malloc(sizeof(int) * tamanho);
    array_baldes = (int*)malloc(sizeof(int) * tamanho);
    struct bucket* baldes = (struct bucket *) calloc(buckets * buckets, sizeof(struct bucket));

    long long start_usec, end_usec, elapsed_usec, min_usec = 0L;
    int i, run;
    int num_hwcntrs = 0;

    for(i=0; i < tamanho; i++) {
        array[i] = rand() % max;
    }

    fprintf(stdout, "\nSetting up PAPI...");
    // Initialize PAPI
    retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT) {
        fprintf(stderr, "PAPI library init error!\n");
        return;
    }

    /* create event set */
    if (PAPI_create_eventset(&EventSet) != PAPI_OK) {
        fprintf(stderr, "PAPI create event set error\n");
        return;
    }

    /* Get the number of hardware counters available */
    if ((num_hwcntrs = PAPI_num_hwctrs()) <= PAPI_OK) {
        fprintf(stderr, "PAPI error getting number of available hardware counters!\n");
        return;
    }
    fprintf(stdout, "done!\nThis system has %d available counters.\n\n", num_hwcntrs);

    // We will be using at most NUM_EVENTS counters
    if (num_hwcntrs >= NUM_EVENTS) {
        num_hwcntrs = NUM_EVENTS;
    }
    else {
        fprintf(stderr, "Error: there aren't enough counters to monitor %d events!\n", NUM_EVENTS);
        return;
    }

    if (PAPI_add_events(EventSet, Events, NUM_EVENTS) != PAPI_OK) {
        fprintf(stderr, "PAPI library add events error!\n");
        return;
    }

    for (run = 0; run < NUM_RUNS; run++) {
        fprintf(stdout, "run=%d...", run);

        start_usec = PAPI_get_real_usec();

        /* Start counting events */
        if (PAPI_start(EventSet) != PAPI_OK) {
            fprintf(stderr, "PAPI error starting counters!\n");
            return 0;
        }

        bucketSortParalelo(array, array_baldes, tamanho, width, buckets, baldes);

        free(array);
        array = array_baldes;

        /* Stop counting events */
        if (PAPI_stop(EventSet, values) != PAPI_OK) {
            fprintf(stderr, "PAPI error stoping counters!\n");
            return;
        }

        end_usec = PAPI_get_real_usec();
        fprintf(stdout, "done!\n");

        elapsed_usec = end_usec - start_usec;

        if ((run == 0) || (elapsed_usec < min_usec)) {
            min_usec = elapsed_usec;
            for (i = 0; i < NUM_EVENTS; i++) min_values[i] = values[i];
        }

        valid(array, tamanho);
        free(baldes);
        free(array_baldes);

    } // end runs
    fprintf(stdout, "\nWall clock time: %lld usecs\n", min_usec);

    // output PAPI counters' values
    for (i = 0; i < NUM_EVENTS; i++) {
        char EventCodeStr[PAPI_MAX_STR_LEN];

        if (PAPI_event_code_to_name(Events[i], EventCodeStr) == PAPI_OK) {
            fprintf(stdout, "%s = %lld\n", EventCodeStr, min_values[i]);
        }
        else {
            fprintf(stdout, "PAPI UNKNOWN EVENT = %lld\n", min_values[i]);
        }
    }

    #if NUM_EVENTS >1
        // evaluate CPI and Texec here
        if ((Events[0] == PAPI_TOT_CYC) && (Events[1] == PAPI_TOT_INS)) {
            float CPI = ((float)min_values[0]) / ((float)min_values[1]);
            fprintf(stdout, "CPI = %.2f\n", CPI);
        }
    #endif

    fprintf(stdout, "\nThat's all\n");
}

void bucketSortSequencial(int array[], int arrayb[], int n, int w, int n_b, struct bucket* baldes) {
    int i, j, b_index, low, high;

    for (i = 0; i < n; i++) {
        j = array[i] / w;
        if (j > n_b - 1)
            j = n_b - 1;
        baldes[j].n_elem++;
    }

    baldes[0].index = 0;
    baldes[0].start = 0;
    for (i = 1; i < n_b; i++) {
        baldes[i].index = baldes[i - 1].index + baldes[i - 1].n_elem;
        baldes[i].start = baldes[i - 1].start + baldes[i - 1].n_elem;
    }

    for (i = 0; i < n; i++) {
        j = array[i] / w;
        if (j > n_b - 1)
            j = n_b - 1;
        b_index = baldes[j].index++;
        arrayb[b_index] = array[i];
    }

    if (n < 20000000) {
        for (i = 0; i < n_b; i++) {
            low = baldes[i].start;
            high = baldes[i].start + baldes[i].n_elem - 1;
            quickSort(arrayb, low, high);
        }
    }
    else {
        for (i = 0; i < n_b; i++) {
            qsort(arrayb + baldes[i].start, baldes[i].n_elem, sizeof(int), cmpfunc);
        }
    }
}

void papi_sequencial(int tamanho, int buckets){
    
    int max = 99999;
    int *array, *array_baldes;

    long long start_usec, end_usec, elapsed_usec, min_usec = 0L;
    int i, run;
    int num_hwcntrs = 0;

    int width = (int)max / buckets;
    array = (int*)malloc(sizeof(int) * tamanho);
    array_baldes = (int*)malloc(sizeof(int) * tamanho);
    struct bucket* baldes = (struct bucket *) calloc(buckets, sizeof(struct bucket));

    for(i=0; i < tamanho; i++) {
        array[i] = rand() % max;
    }

    fprintf(stdout, "\nSetting up PAPI...");
    // Initialize PAPI
    retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT) {
        fprintf(stderr, "PAPI library init error!\n");
        return;
    }

    /* create event set */
    if (PAPI_create_eventset(&EventSet) != PAPI_OK) {
        fprintf(stderr, "PAPI create event set error\n");
        return;
    }

    /* Get the number of hardware counters available */
    if ((num_hwcntrs = PAPI_num_hwctrs()) <= PAPI_OK) {
        fprintf(stderr, "PAPI error getting number of available hardware counters!\n");
        return;
    }
    fprintf(stdout, "done!\nThis system has %d available counters.\n\n", num_hwcntrs);

    // We will be using at most NUM_EVENTS counters
    if (num_hwcntrs >= NUM_EVENTS) {
        num_hwcntrs = NUM_EVENTS;
    }
    else {
        fprintf(stderr, "Error: there aren't enough counters to monitor %d events!\n", NUM_EVENTS);
        return;
    }

    if (PAPI_add_events(EventSet, Events, NUM_EVENTS) != PAPI_OK) {
        fprintf(stderr, "PAPI library add events error!\n");
        return;
    }

    for (run = 0; run < NUM_RUNS; run++) {
        fprintf(stdout, "run=%d...", run);

        start_usec = PAPI_get_real_usec();

        /* Start counting events */
        if (PAPI_start(EventSet) != PAPI_OK) {
            fprintf(stderr, "PAPI error starting counters!\n");
            return;
        }

        bucketSortSequencial(array, array_baldes, tamanho, width, buckets, baldes);

        free(array);
        array = array_baldes;

        /* Stop counting events */
        if (PAPI_stop(EventSet, values) != PAPI_OK) {
            fprintf(stderr, "PAPI error stoping counters!\n");
            return;
        }

        end_usec = PAPI_get_real_usec();
        fprintf(stdout, "done!\n");

        elapsed_usec = end_usec - start_usec;

        if ((run == 0) || (elapsed_usec < min_usec)) {
            min_usec = elapsed_usec;
            for (i = 0; i < NUM_EVENTS; i++) min_values[i] = values[i];
        }

        valid(array,tamanho);
        free(baldes);
        free(array_baldes);

    } // end runs
    fprintf(stdout, "\nWall clock time: %lld usecs\n", min_usec);

    // output PAPI counters' values
    for (i = 0; i < NUM_EVENTS; i++) {
        char EventCodeStr[PAPI_MAX_STR_LEN];

        if (PAPI_event_code_to_name(Events[i], EventCodeStr) == PAPI_OK) {
            fprintf(stdout, "%s = %lld\n", EventCodeStr, min_values[i]);
        }
        else {
            fprintf(stdout, "PAPI UNKNOWN EVENT = %lld\n", min_values[i]);
        }
    }

    #if NUM_EVENTS >1
        // evaluate CPI and Texec here
        if ((Events[0] == PAPI_TOT_CYC) && (Events[1] == PAPI_TOT_INS)) {
            float CPI = ((float)min_values[0]) / ((float)min_values[1]);
            fprintf(stdout, "CPI = %.2f\n", CPI);
        }
    #endif

    fprintf(stdout, "\nThat's all\n");
}

int main() {
    int array_size = 50000000;
    int n_baldes =  10;

    printf("\n\tsequencial\n");
    papi_sequencial(array_size, n_baldes);


    printf("\n\tparalelo\n");
    papi_paralelo(array_size, n_baldes);

    return 0;
}