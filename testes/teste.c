#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

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

int cmpfunc (const void * a, const void * b){
 return ( *(int*)a - *(int*)b );
}

typedef struct bucket {
	int n_elem;
	int start;
	int index;
}*Bucket;


void insertionSort(int arr[], int n) {
    int i, key, j;
    for (i = 1; i < n; i++) {
        key = arr[i];
        j = i - 1;

        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }

        arr[j + 1] = key;
    }
}

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


void bucketSortSequencialv2(int array[], int arrayb[], int n, int w, int n_b, struct bucket* baldes) {
    int i, j, b_index, low, high;

    for (i = 0; i < n; i++) {
        j = array[i]/w;
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

    
    for (i = 0; i < n_b; i++){
        insertionSort(arrayb+baldes[i].start, baldes[i].n_elem);
    }
}

void bucketSortSequencial(int array[], int arrayb[], int n, int w, int n_b, struct bucket* baldes) {
    int i, j, b_index, low, high;

    for (i = 0; i < n; i++) {
        j = array[i]/w;
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

    if(n_b < 5000 && n < 3000000){
    	for (i = 0; i < n_b; i++){
        	insertionSort(arrayb+baldes[i].start, baldes[i].n_elem);
    	}
    } else {
	    if(n < 20000000){
	        for (i = 0; i < n_b; i++) {
	            low = baldes[i].start;
	            high = baldes[i].start + baldes[i].n_elem - 1;
	            quickSort(arrayb, low, high);
	        }
	    } else {
	        for (i = 0; i < n_b; i++){
	            qsort(arrayb+baldes[i].start, baldes[i].n_elem, sizeof(int), cmpfunc);
	        }
	    }
    }
}

float testar_sequencial(int tamanho, int buckets){
	int max = 99999, i;
	int *array, * array_baldes;
	
	int width = (int)max / buckets;
	array = (int*)malloc(sizeof(int) * tamanho);
	array_baldes = (int*)malloc(sizeof(int) * tamanho);
	struct bucket* baldes = (struct bucket *) calloc(buckets, sizeof(struct bucket));

	for(i=0; i < tamanho; i++) {
		array[i] = rand() % max;
	}

	double t = omp_get_wtime();

	bucketSortSequencial(array, array_baldes, tamanho, width, buckets, baldes);

	float t_exec = omp_get_wtime() - t;

	free(array);
	array = array_baldes;

	//printf("\n(q) Sorting %d elements using %d buckets took %f seconds.\n", tamanho, buckets, t_exec);
	//valid(array,tamanho);

	free(baldes);
	free(array_baldes);

	return t_exec;
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


	    if(array_size < 2000000){
			#pragma omp for private(i)
			for (i = 0; i < n_baldes; i++){
				insertionSort(array_baldes+global_starting_position[i], global_n_elem[i]);
			}
	    } else {
			#pragma omp for private(i)
			for (int i = 0; i < n_baldes; i++){
				qsort(array_baldes+global_starting_position[i], global_n_elem[i], sizeof(int), cmpfunc);
			}
		}
	}
}

void bucketSortParalelov1(int array[], int array_baldes[], int array_size, int width, int n_baldes, struct bucket* baldes) {
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

void bucketSortParalelov2(int array[], int array_baldes[], int array_size, int width, int n_baldes, struct bucket* baldes) {
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
			insertionSort(array_baldes+global_starting_position[i], global_n_elem[i]);
		}

	}
}


float testar_paralela(int tamanho, int buckets){
	int max = 99999, i;
	int *array, * array_baldes;
	int width = (int)max / buckets;
	array = (int*)malloc(sizeof(int) * tamanho);
	array_baldes = (int*)malloc(sizeof(int) * tamanho);
	struct bucket* baldes = (struct bucket *) calloc(buckets * buckets, sizeof(struct bucket));
	
	for(i=0; i < tamanho; i++) {
		array[i] = rand() % max;
	}

	double t = omp_get_wtime();

	bucketSortParalelo(array, array_baldes, tamanho, width, buckets, baldes);

	float t_exec = omp_get_wtime() - t;
	//printf("\n(q) Sorting %d elements using %d buckets took %f seconds.\n", tamanho, buckets, t_exec);
	

	free(array);
	array = array_baldes;

	//valid(array,tamanho);

	free(baldes);
	free(array_baldes);

	return t_exec;
}

float testar_insert_paralela(int tamanho, int buckets){
	int max = 99999, i;
	int *array, * array_baldes;
	int width = (int)max / buckets;
	array = (int*)malloc(sizeof(int) * tamanho);
	array_baldes = (int*)malloc(sizeof(int) * tamanho);
	struct bucket* baldes = (struct bucket *) calloc(buckets * buckets, sizeof(struct bucket));
	
	for(i=0; i < tamanho; i++) {
		array[i] = rand() % max;
	}

	double t = omp_get_wtime();

	bucketSortParalelov2(array, array_baldes, tamanho, width, buckets, baldes);

	float t_exec = omp_get_wtime() - t;
	//printf("\n(q) Sorting %d elements using %d buckets took %f seconds.\n", tamanho, buckets, t_exec);
	

	free(array);
	array = array_baldes;

	//valid(array,tamanho);

	free(baldes);
	free(array_baldes);

	return t_exec;
}


float testar_insert_sequencial(int tamanho, int buckets){
	int max = 99999, i;
	int *array, * array_baldes;
	
	int width = (int)max / buckets;
	array = (int*)malloc(sizeof(int) * tamanho);
	array_baldes = (int*)malloc(sizeof(int) * tamanho);
	struct bucket* baldes = (struct bucket *) calloc(buckets, sizeof(struct bucket));

	for(i=0; i < tamanho; i++) {
		array[i] = rand() % max;
	}

	double t = omp_get_wtime();

	bucketSortSequencialv2(array, array_baldes, tamanho, width, buckets, baldes);

	float t_exec = omp_get_wtime() - t;

	free(array);
	array = array_baldes;

	// printf("\n(i) Sorting %d elements using %d buckets took %f seconds.\n", tamanho, buckets, t_exec);
	//valid(array,tamanho);

	free(baldes);
	free(array_baldes);

	return t_exec;
}


int comparar_QS_IS(int tamanho, int max_b){

	float qsort, isort;

	int ciclos = 0, i;

	float qsort_avg = 0;
	float isort_avg = 0;

	int qsort_wins = 0;
	int insertsort_wins = 0;

	for(i = 5; i <= max_b; i+=5){

		qsort = testar_sequencial(tamanho, i);

		qsort_avg = qsort_avg * ciclos;
		qsort_avg += qsort;
		qsort_avg = qsort_avg / (ciclos+1);

		isort = testar_insert_sequencial(tamanho, i);

		isort_avg = isort_avg * ciclos;
		isort_avg += isort;
		isort_avg = isort_avg / (ciclos+1);

		if(qsort > isort){
			insertsort_wins++;
		} else {
			qsort_wins++;
		}

		ciclos++;
		//printf("%d,%d,%.08f,%.08f,%d,%d,%.08f,%.08f\n", tamanho, i, isort, qsort, insertsort_wins, qsort_wins, isort_avg, qsort_avg);
	}

	printf("\nsize: %d (buckets 1 até %d) |  (i): %d , (q): %d.\n", tamanho, max_b, insertsort_wins, qsort_wins);
	printf("\n\tQuick Sort avg = <%.08f>", qsort_avg);
	// printf("\n\tInsertion Sort avg = <%.08f>", isort_avg);
	// printf("\nNúmero de testes: <%d>\n", ciclos);

	return ciclos;
}

int comparar_QS_ISv2(int tamanho, int max_b){

	float qsort, isort;

	int ciclos = 0, i;

	float qsort_avg = 0;
	float isort_avg = 0;

	int qsort_wins = 0;
	int insertsort_wins = 0;

	for(i = max_b/3; i <= max_b; i+=max_b/3){

		qsort = testar_sequencial(tamanho, i);

		qsort_avg = qsort_avg * ciclos;
		qsort_avg += qsort;
		qsort_avg = qsort_avg / (ciclos+1);

		isort = testar_insert_sequencial(tamanho, i);

		isort_avg = isort_avg * ciclos;
		isort_avg += isort;
		isort_avg = isort_avg / (ciclos+1);

		if(qsort > isort){
			insertsort_wins++;
		} else {
			qsort_wins++;
		}

		ciclos++;
		//printf("%d,%d,%.08f,%.08f,%d,%d,%.08f,%.08f\n", tamanho, i, isort, qsort, insertsort_wins, qsort_wins, isort_avg, qsort_avg);
	}

	//printf("\nsize: %d (buckets 1 até %d) |  (i): %d , (q): %d.\n", tamanho, max_b, insertsort_wins, qsort_wins);
	//printf("\n\tQuick Sort avg = <%.08f>", qsort_avg);
	// printf("\n\tInsertion Sort avg = <%.08f>", isort_avg);
	// printf("\nNúmero de testes: <%d>\n", ciclos);

	return insertsort_wins;
}


int comparar_QS_IS_paralelo(int tamanho, int max_b, int incr){

	float qsort, isort;

	int ciclos = 0, i;

	float qsort_avg = 0;
	float isort_avg = 0;

	int qsort_wins = 0;
	int insertsort_wins = 0;

	for(i = 100; i <= max_b; i+= incr){

		qsort = testar_paralela(tamanho, i);

		qsort_avg = qsort_avg * ciclos;
		qsort_avg += qsort;
		qsort_avg = qsort_avg / (ciclos+1);

		isort = testar_insert_paralela(tamanho, i);

		isort_avg = isort_avg * ciclos;
		isort_avg += isort;
		isort_avg = isort_avg / (ciclos+1);

		if(qsort > isort){
			insertsort_wins++;
		} else {
			qsort_wins++;
		}

		ciclos++;
		printf("%d,%d,%.08f,%.08f,%d,%d,%.08f,%.08f\n", tamanho, i, isort, qsort, insertsort_wins, qsort_wins, isort_avg, qsort_avg);
	}

	//printf("\nsize: %d (buckets 1 até %d) |  (i): %d , (q): %d.\n", tamanho, max_b, insertsort_wins, qsort_wins);
	//printf("\n\tQuick Sort avg = <%.08f>\n", qsort_avg);
	//printf("\n\tInsertion Sort avg = <%.08f>\n", isort_avg);
	// printf("\nNúmero de testes: <%d>\n", ciclos);

	return insertsort_wins;
}


// verificar se compensa usar insertion sort
void teste1(){
	int tamanho, ciclos = 0;

	#pragma omp parallel
	{
		#pragma omp for
		for(tamanho = 500; tamanho <= 40000; tamanho+=500){
			ciclos += comparar_QS_IS(tamanho, 10);
		}
		#pragma omp barrier


		#pragma omp for
		for(tamanho = 500; tamanho <= 40000; tamanho+=500){
			ciclos += comparar_QS_IS(tamanho, 100);
		}
		#pragma omp barrier

		#pragma omp for
		for(tamanho = 500; tamanho <= 40000; tamanho+=500){
			ciclos += comparar_QS_IS(tamanho, 1000);
		}
		#pragma omp barrier

		#pragma omp for
		for(tamanho = 500; tamanho <= 40000; tamanho+=500){
			ciclos += comparar_QS_IS(tamanho, 10000);
		}
	}
}

void teste2(){
	int tamanho, ciclos = 0, total = (5000000-40000)/10000;

	omp_set_num_threads(8);
	#pragma omp parallel
	{
		#pragma omp for
		for(tamanho = 40000; tamanho <= 5000000; tamanho+=10000){
			ciclos += comparar_QS_ISv2(tamanho, 10000);
		}
	}

	printf("ganhou %d out of %d", ciclos, total);
}

void teste3(){

	comparar_QS_IS_paralelo(5000000, 9100, 500);

	comparar_QS_IS_paralelo(4000000, 9100, 500);

	comparar_QS_IS_paralelo(3000000, 9100, 500);

	comparar_QS_IS_paralelo(2000000, 9100, 500);

	comparar_QS_IS_paralelo(1000000, 9100, 500);

	comparar_QS_IS_paralelo(500000, 9100, 500);

	comparar_QS_IS_paralelo(200000, 9100, 500);

	comparar_QS_IS_paralelo(100000, 9100, 500);

	comparar_QS_IS_paralelo(50000000, 9100, 500);

	comparar_QS_IS_paralelo(30000000, 9100, 500);
}

void teste4(){
	int tamanho, buckets = 10;
	float p, s;
	for(tamanho = 100000; tamanho <= 1000000; tamanho+=100000){
		s = testar_sequencial(tamanho, buckets);
		p = testar_paralela(tamanho, buckets);
	
		printf("%d,%d,%.08f,%.08f\n", tamanho, buckets, s, p);
	}

	buckets = 100;

	for(tamanho = 100000; tamanho <= 2000000; tamanho+=100000){
		s = testar_sequencial(tamanho, buckets);
		p = testar_paralela(tamanho, buckets);

		printf("%d,%d,%.08f,%.08f\n", tamanho, buckets, s, p);
	}

	buckets = 1000;

	for(tamanho = 100000; tamanho <= 2000000; tamanho+=100000){
		s = testar_sequencial(tamanho, buckets);
		p = testar_paralela(tamanho, buckets);

		printf("%d,%d,%.08f,%.08f\n", tamanho, buckets, s, p);
	}

	buckets = 2000;

	for(tamanho = 100000; tamanho <= 2000000; tamanho+=100000){
		s = testar_sequencial(tamanho, buckets);
		p = testar_paralela(tamanho, buckets);

		printf("%d,%d,%.08f,%.08f\n", tamanho, buckets, s, p);
	}

	buckets = 3000;

	for(tamanho = 100000; tamanho <= 2000000; tamanho+=100000){
		s = testar_sequencial(tamanho, buckets);
		p = testar_paralela(tamanho, buckets);

		printf("%d,%d,%.08f,%.08f\n", tamanho, buckets, s, p);
	}
}

void teste5(){

	int tamanho, buckets = 10;
	float p, s;
	buckets = 2000;

	for(tamanho = 100000; tamanho <= 50100000; tamanho+=1000000){
		s = testar_sequencial(tamanho, buckets);
		p = testar_paralela(tamanho, buckets);

		printf("%d,%d,%.08f,%.08f\n", tamanho, buckets, s, p);
	}
}

void teste7(){
	int tamanho, buckets = 10;
	float p;
	for(tamanho = 100000; tamanho <= 10000000; tamanho+=100000){
		p = testar_paralela(tamanho, buckets);
	
		printf("%d,%d,%.08f\n", tamanho, buckets, p);
	}

	buckets = 16;
	for(tamanho = 100000; tamanho <= 10000000; tamanho+=100000){
		p = testar_paralela(tamanho, buckets);
	
		printf("%d,%d,%.08f\n", tamanho, buckets, p);
	}

	buckets = 24;
	for(tamanho = 100000; tamanho <= 10000000; tamanho+=100000){
		p = testar_paralela(tamanho, buckets);
	
		printf("%d,%d,%.08f\n", tamanho, buckets, p);
	}

	buckets = 32;
	for(tamanho = 100000; tamanho <= 10000000; tamanho+=100000){
		p = testar_paralela(tamanho, buckets);
	
		printf("%d,%d,%.08f\n", tamanho, buckets, p);
	}

	buckets = 64;
	for(tamanho = 100000; tamanho <= 10000000; tamanho+=100000){
		p = testar_paralela(tamanho, buckets);
	
		printf("%d,%d,%.08f\n", tamanho, buckets, p);
	}
}

int main(){

	teste5();

	return 0;
}