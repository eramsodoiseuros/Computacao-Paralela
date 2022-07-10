#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

/*
	Desenvolvimento da versão sequencial do algoritmo, em C
	para ordenar um vetor de inteiros, incluindo
	a análise do impacto de possíveis otimizações à versão sequencial.

*/

typedef struct bucket {
	int n_elem;
	int start;
	int index;
}*Bucket;


int cmpfunc (const void * a, const void * b){
 return ( *(int*)a - *(int*)b );
}

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

void swap(int *a, int *b) {
  int t = *a;
  *a = *b;
  *b = t;
}

int partition(int array[], int low, int high) {
  int pivot = array[high];
  int i = (low - 1);

  for (int j = low; j < high; j++) {
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

void main(){
	
	int *array, *array_baldes;
	int array_size, n_baldes, max, width, i, j, low, high, valid = 1;
	double t;
	float t_exec;
	struct bucket *baldes;


	printf("Insert the length of array to sort \n");
	if (scanf("%d", &array_size) != 1){
		printf("\n\tNo valid input was given. The program will not resume.\n");
		return;
	}
	
	printf("Insert the max number for random generated array members.\n");
	if (scanf("%d", &max) != 1){
		printf("\n\tNo valid input was given, default value was set to n=1000.\n");
		max = 1000;
	}

	printf("Insert the number of buckets \n");
	if (scanf("%d", &n_baldes) != 1){
		printf("\n\tNo valid input was given, default value was set to n=10.\n");
		n_baldes = 10;
	}

	width = (int)max/n_baldes;
	array = (int *) malloc(sizeof(int)*array_size);
	array_baldes = (int *) malloc(sizeof(int)*array_size);
	baldes = (struct bucket *) calloc(n_baldes, sizeof(struct bucket));

	for(i=0; i < array_size; i++) {
		array[i] = random() % max;
	}

	printf("The input RNG Array was: [");
	for(i=0; i < array_size; i++) {
			printf("<%d> ",array[i]);
		}
	printf("]\n");

	t = omp_get_wtime();

	for (i=0; i < array_size; i++){
		j = array[i]/width;
		if (j > n_baldes-1)
				j = n_baldes-1;
		baldes[j].n_elem++;
	}

	baldes[0].index = 0;
	baldes[0].start = 0;
	for (i=1; i < n_baldes;i++){
		baldes[i].index = baldes[i-1].index + baldes[i-1].n_elem;
		baldes[i].start = baldes[i-1].start + baldes[i-1].n_elem;
	}

	int b_index;
	for (i=0; i < array_size; i++){
		j = array[i]/width;
		if (j > n_baldes -1)
				j = n_baldes-1;
		b_index = baldes[j].index++;
		array_baldes[b_index] = array[i];
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

	t_exec = omp_get_wtime() - t;

	free(array);
	array = array_baldes;

	printf("\nSorting %d elements took %f seconds.\n", array_size, t_exec);

	for(i=0; i < array_size-1; i++) {
		if (array[i] > array[i+1]){
			valid = 0;
			break;
		}
  }

  if(valid)
		printf("\n\tThe array was correctly sorted.\n");
	else printf("\n\tThe array was not correctly sorted.\n");

	for(i=0; i < array_size; i++) {
			printf("<%d> ",array[i]);
		}
	printf("]\n");

	free(array_baldes);
	free(baldes);
}