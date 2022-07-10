#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>


void print(int* a, int n){
  printf("\t\nSize: %d\n", n);
  for (int i = 0; i < n; ++i)
  {
    printf(" <%d> ", a[i]);
  }
  printf("\n\n");
}

int cmpfunc (const void * a, const void * b)
{
 return ( *(int*)a - *(int*)b );
}

typedef struct bucket {
	int n_elem;
	int start;
	int index;
}*Bucket;

int main(){
	
	int *array, *array_baldes;
  int array_size, n_baldes, num_threads, max, width, i, b_index, valid = 1;
	double t;
	float t_exec;
	struct bucket *baldes;


	printf("Insert the length of array to sort \n");
	if (scanf("%d", &array_size) != 1){
		printf("\n\tNo valid input was given. The program will not resume.\n");
		return -1;
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

  num_threads = n_baldes;
  omp_set_num_threads(num_threads);

  width = (int)max/n_baldes;
  array = (int *) malloc(sizeof(int) * array_size);
  array_baldes = (int *) malloc(sizeof(int) * array_size);
  baldes = (struct bucket *) calloc(n_baldes * num_threads, sizeof(struct bucket));

  int global_n_elem[n_baldes];
  int global_starting_position[n_baldes];
  memset(global_n_elem, 0, sizeof(int) * n_baldes);
  memset(global_starting_position, 0, sizeof(int) * n_baldes);

  for(i=0; i < array_size; i++) {
    array[i] = rand() % max;
  }

  t = omp_get_wtime();

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

  t_exec = omp_get_wtime() - t;

  free(array);
  array = array_baldes;

  printf("\nSorting %d elements took %f seconds.\n", array_size, t_exec);


  for (i=0; i < array_size-1; i++) {
    if (array[i] > array[i+1]){
      valid = 0;
      printf("\n\tfalhou na posicao: %d, <%d> <!%d!> <%d>", i, array[i-1], array[i], array[i+1]);
      break;
    }
  }

  if (valid) printf("\n\tThe array was correctly sorted.\n");
  else printf("\n\tThe array was not correctly sorted.\n");
  free(baldes);
  free(array_baldes);
  
  return 0;
}