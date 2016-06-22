/*
* Purpose:  Implement parallel odd-even sort of an array of
*           nonegative ints
* Input:
*    A:     elements of array (optional)
* Output:
*    A:     elements of A after sorting
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#define COUNT 33554432
#define MPI_SWITCH 
//#define SORT_SREIAL 
const int RMAX = 10000;

/* Local functions */
void Usage(char* program);
void Print_list(double local_A[], int local_n, int rank);
void Merge_low(double local_A[], double temp_B[], double temp_C[],
	int local_n);
void Merge_high(double local_A[], double temp_B[], double temp_C[],
	int local_n);
void Generate_list(double local_A[], int local_n, int my_rank);
int  Compare(const void* a_p, const void* b_p);

/* Functions involving communication */
void Get_args(int argc, char* argv[], int* global_n_p, int* local_n_p,
	char* gi_p, int my_rank, int p, MPI_Comm comm);
void Sort(double local_A[], int local_n, int my_rank,
	int p, MPI_Comm comm);
void Odd_even_iter(double local_A[], double temp_B[], double temp_C[],
	int local_n, int phase, int even_partner, int odd_partner,
	int my_rank, int p, MPI_Comm comm);
void Print_local_lists(double local_A[], int local_n,
	int my_rank, int p, MPI_Comm comm);
void Print_global_list(double local_A[], int local_n, int my_rank,
	int p, MPI_Comm comm);
void Read_list(double local_A[], int local_n, int my_rank, int p,
	MPI_Comm comm);
void Print_seri_list(double local_A[], int local_n);


/*-------------------------------------------------------------------*/
int main(int argc, char* argv[]) {
	printf("start：");
	int my_rank, p;
	char g_i;
	double *local_A;
	double *local_B;
	int global_n;
	int local_n;
	int k;
	clock_t start, end; //记录运行和结束时间
	double time[p];  //记录每个进程运行时间
	double times;
	local_B = (double*)malloc(COUNT*sizeof(double));
#ifdef MPI_SWITCH
	MPI_Comm comm;
	

	MPI_Init(&argc, &argv);
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &p);
	MPI_Comm_rank(comm, &my_rank);

	Get_args(argc, argv, &global_n, &local_n, &g_i, my_rank, p, comm);
	local_A = (double*)malloc(local_n*sizeof(double));
	
	if (g_i == 'g') {
		Generate_list(local_A, local_n, my_rank);
		//Print_local_lists(local_A, local_n, my_rank, p, comm);	   
	}
	else {
		Read_list(local_A, local_n, my_rank, p, comm);
#     ifdef DEBUG
		Print_local_lists(local_A, local_n, my_rank, p, comm);
#     endif
	}

#  ifdef DEBUG
	printf("Proc %d > Before Sort\n", my_rank);
	fflush(stdout);
#  endif
	
	//int x = 0;
	//for (k = my_rank*local_n; k < (my_rank + 1)*local_n; k++){
	//	local_B[k] = local_A[x];
	//	/*printf("%f ", local_B[k]);
	//	printf(" ");*/
	//	x++;
	//}
	//Print_global_list(local_B, local_n, my_rank, p, comm);

	start = clock();
	Sort(local_A, local_n, my_rank, p, comm);
	end = clock();
	times = (double)(end - start) / CLOCKS_PER_SEC;
	printf("进程%d并行奇偶排序所花时间为：%f%c\n", my_rank, times, 's');

	if (my_rank != 0){
		/* Create message */
		//sprintf(timechar, ecvt(times,ndig,&dec,&sign));
		/* Send message to process 0 */
		MPI_Send(&times, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	else{
		/* Print my message */
		//printf("Greetings from process %d of %d!\n", my_rank, p);
		time[0] = times;
		int q;
		for (q = 1; q < p; q++) {
			/* Receive message from process q */
			MPI_Recv(&times, 1, MPI_DOUBLE, q,
				0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			/* Print message from process q */
			time[q] = times;
		}
	}
#  ifdef DEBUG
	Print_local_lists(local_A, local_n, my_rank, p, comm);
	fflush(stdout);
#  endif

	//Print_global_list(local_A, local_n, my_rank, p, comm);

	free(local_A);
	double pallTime ;
	
	if (my_rank == 0){
		pallTime = time[0];
		int i;
		for (i = 0; i<p; i++){
			if (pallTime <= time[i])
				pallTime = time[i];
		}

		Generate_list(local_B, COUNT, 0);
		//print_seri_list(local_b, count);
		clock_t st = clock();
		qsort(local_B, COUNT, sizeof(double), Compare);
		double seriTime = (double)(clock() - st) / CLOCKS_PER_SEC;
		printf("串行排序所花时间为：%f%c\n", seriTime, 's');
		//print_seri_list(local_B, count);
		free(local_B);
		printf("并行排序所花时间为：%f%c\n", pallTime, 's');
		printf("加速比为：%f\n", (double)(seriTime / pallTime));

	}
	MPI_Finalize();
#endif

#ifdef SORT_SREIAL	
	Generate_list(local_B, COUNT, 0);
	//print_seri_list(local_b, count);
	clock_t st = clock();
	qsort(local_B, COUNT, sizeof(double), Compare);
	printf("串行排序所花时间为：%f%c\n", (double)(clock() - st)/CLOCKS_PER_SEC, 's');
	//print_seri_list(local_B, count);
	
	free(local_B);
#endif

	return 0;
}  /* main */


/*-------------------------------------------------------------------
* Function:   Generate_list
* Purpose:    Fill list with random ints
* Input Args: local_n, my_rank
* Output Arg: local_A
*/
void Generate_list(double local_A[], int local_n, int my_rank) {
	int i;

	srandom(my_rank + 1);
	for (i = 0; i < local_n; i++)
		local_A[i] = random() % RMAX + 0.001;

}  /* Generate_list */


/*-------------------------------------------------------------------
* Function:  Usage
* Purpose:   Print command line to start program
* In arg:    program:  name of executable
* Note:      Purely local, run only by process 0;
*/
void Usage(char* program) {
	fprintf(stderr, "usage:  mpirun -np <p> %s <g|i> <global_n>\n",
		program);
	fprintf(stderr, "   - p: the number of processes \n");
	fprintf(stderr, "   - g: generate random, distributed list\n");
	fprintf(stderr, "   - i: user will input list on process 0\n");
	fprintf(stderr, "   - global_n: number of elements in global list");
	fprintf(stderr, " (must be evenly divisible by p)\n");
	fflush(stderr);
}  /* Usage */


/*-------------------------------------------------------------------
* Function:    Get_args
* Purpose:     Get and check command line arguments
* Input args:  argc, argv, my_rank, p, comm
* Output args: global_n_p, local_n_p, gi_p
*/
void Get_args(int argc, char* argv[], int* global_n_p, int* local_n_p,
	char* gi_p, int my_rank, int p, MPI_Comm comm) {

	if (my_rank == 0) {
		if (argc != 3) {
			Usage(argv[0]);
			*global_n_p = -1;  /* Bad args, quit */
		}
		else {
			*gi_p = argv[1][0];
			if (*gi_p != 'g' && *gi_p != 'i') {
				Usage(argv[0]);
				*global_n_p = -1;  /* Bad args, quit */
			}
			else {
				*global_n_p = atoi(argv[2]);
				if (*global_n_p % p != 0) {
					Usage(argv[0]);
					*global_n_p = -1;
				}
			}
		}
	}  /* my_rank == 0 */

	MPI_Bcast(gi_p, 1, MPI_CHAR, 0, comm);
	MPI_Bcast(global_n_p, 1, MPI_INT, 0, comm);

	if (*global_n_p <= 0) {
		MPI_Finalize();
		exit(-1);
	}

	*local_n_p = *global_n_p / p;
#  ifdef DEBUG
	printf("Proc %d > gi = %c, global_n = %d, local_n = %d\n",
		my_rank, *gi_p, *global_n_p, *local_n_p);
	fflush(stdout);
#  endif

}  /* Get_args */


/*-------------------------------------------------------------------
* Function:   Read_list
* Purpose:    process 0 reads the list from stdin and scatters it
*             to the other processes.
* In args:    local_n, my_rank, p, comm
* Out arg:    local_A
*/
void Read_list(double local_A[], int local_n, int my_rank, int p,
	MPI_Comm comm) {
	int i;
	double *temp;

	if (my_rank == 0) {
		temp = (double*)malloc(p*local_n*sizeof(double));
		printf("Enter the elements of the list\n");
		for (i = 0; i < p*local_n; i++)
			scanf("%f", &temp[i]);
	}

	MPI_Scatter(temp, local_n, MPI_DOUBLE, local_A, local_n, MPI_DOUBLE,
		0, comm);

	if (my_rank == 0)
		free(temp);
}  /* Read_list */


/*-------------------------------------------------------------------
* Function:   Print_global_list
* Purpose:    Print the contents of the global list A
* Input args:
*    n, the number of elements
*    A, the list
* Note:       Purely local, called only by process 0
*/
void Print_global_list(double local_A[], int local_n, int my_rank, int p,
	MPI_Comm comm) {
	double* A;
	int i, n;

	if (my_rank == 0) {
		n = p*local_n;
		A = (double*)malloc(n*sizeof(double));
		MPI_Gather(local_A, local_n, MPI_DOUBLE, A, local_n, MPI_DOUBLE, 0,
			comm);
		printf("Global list:\n");
		for (i = 0; i < n; i++)
			printf("%f ", A[i]);
		printf("\n\n");
		free(A);
	}
	else {
		MPI_Gather(local_A, local_n, MPI_DOUBLE, A, local_n, MPI_DOUBLE, 0,
			comm);
	}

}  /* Print_global_list */

/*-------------------------------------------------------------------
* Function:    Compare
* Purpose:     Compare 2 ints, return -1, 0, or 1, respectively, when
*              the first int is less than, equal, or greater than
*              the second.  Used by qsort.
*/
int Compare(const void* a_p, const void* b_p) {
	double a = *((double*)a_p);
	double b = *((double*)b_p);

	if (a < b)
		return -1;
	else if (a == b)
		return 0;
	else /* a > b */
		return 1;
}  /* Compare */

/*-------------------------------------------------------------------
* Function:    Sort
* Purpose:     Sort local list, use odd-even sort to sort
*              global list.
* Input args:  local_n, my_rank, p, comm
* In/out args: local_A
*/
void Sort(double local_A[], int local_n, int my_rank,
	int p, MPI_Comm comm) {
	int phase;
	double *temp_B, *temp_C;
	int even_partner;  /* phase is even or left-looking */
	int odd_partner;   /* phase is odd or right-looking */

	/* Temporary storage used in merge-split */
	temp_B = (double*)malloc(local_n*sizeof(double));
	temp_C = (double*)malloc(local_n*sizeof(double));

	/* Find partners:  negative rank => do nothing during phase */
	if (my_rank % 2 != 0) {
		even_partner = my_rank - 1;
		odd_partner = my_rank + 1;
		if (odd_partner == p) odd_partner = MPI_PROC_NULL;  // Idle during odd phase
	}
	else {
		even_partner = my_rank + 1;
		if (even_partner == p) even_partner = MPI_PROC_NULL;  // Idle during even phase
		odd_partner = my_rank - 1;
	}

	/* Sort local list using built-in quick sort */
	qsort(local_A, local_n, sizeof(double), Compare);

#  ifdef DEBUG
	printf("Proc %d > before loop in sort\n", my_rank);
	fflush(stdout);
#  endif

	for (phase = 0; phase < p; phase++)
		Odd_even_iter(local_A, temp_B, temp_C, local_n, phase,
		even_partner, odd_partner, my_rank, p, comm);

	free(temp_B);
	free(temp_C);
}  /* Sort */


/*-------------------------------------------------------------------
* Function:    Odd_even_iter
* Purpose:     One iteration of Odd-even transposition sort
* In args:     local_n, phase, my_rank, p, comm
* In/out args: local_A
* Scratch:     temp_B, temp_C
*/
void Odd_even_iter(double local_A[], double temp_B[], double temp_C[],
	int local_n, int phase, int even_partner, int odd_partner,
	int my_rank, int p, MPI_Comm comm) {
	MPI_Status status;

	if (phase % 2 == 0) {
		if (even_partner >= 0) {
			MPI_Sendrecv(local_A, local_n, MPI_DOUBLE, even_partner, 0,
				temp_B, local_n, MPI_DOUBLE, even_partner, 0, comm,
				&status);
			if (my_rank % 2 != 0)
				Merge_high(local_A, temp_B, temp_C, local_n);
			else
				Merge_low(local_A, temp_B, temp_C, local_n);
		}
	}
	else { /* odd phase */
		if (odd_partner >= 0) {
			MPI_Sendrecv(local_A, local_n, MPI_DOUBLE, odd_partner, 0,
				temp_B, local_n, MPI_DOUBLE, odd_partner, 0, comm,
				&status);
			if (my_rank % 2 != 0)
				Merge_low(local_A, temp_B, temp_C, local_n);
			else
				Merge_high(local_A, temp_B, temp_C, local_n);
		}
	}
}  /* Odd_even_iter */


/*-------------------------------------------------------------------
* Function:    Merge_low
* Purpose:     Merge the smallest local_n elements in my_keys
*              and recv_keys into temp_keys.  Then copy temp_keys
*              back into my_keys.
* In args:     local_n, recv_keys
* In/out args: my_keys
* Scratch:     temp_keys
*/
void Merge_low(
	double  my_keys[],     /* in/out    */
	double  recv_keys[],   /* in        */
	double  temp_keys[],   /* scratch   */
	int  local_n        /* = n/p, in */) {
	int m_i, r_i, t_i;

	m_i = r_i = t_i = 0;
	while (t_i < local_n) {
		if (my_keys[m_i] <= recv_keys[r_i]) {
			temp_keys[t_i] = my_keys[m_i];
			t_i++; m_i++;
		}
		else {
			temp_keys[t_i] = recv_keys[r_i];
			t_i++; r_i++;
		}
	}

	memcpy(my_keys, temp_keys, local_n*sizeof(double));
}  /* Merge_low */

/*-------------------------------------------------------------------
* Function:    Merge_high
* Purpose:     Merge the largest local_n elements in local_A
*              and temp_B into temp_C.  Then copy temp_C
*              back into local_A.
* In args:     local_n, temp_B
* In/out args: local_A
* Scratch:     temp_C
*/
void Merge_high(double local_A[], double temp_B[], double temp_C[],
	int local_n) {
	int ai, bi, ci;

	ai = local_n - 1;
	bi = local_n - 1;
	ci = local_n - 1;
	while (ci >= 0) {
		if (local_A[ai] >= temp_B[bi]) {
			temp_C[ci] = local_A[ai];
			ci--; ai--;
		}
		else {
			temp_C[ci] = temp_B[bi];
			ci--; bi--;
		}
	}

	memcpy(local_A, temp_C, local_n*sizeof(double));
}  /* Merge_high */


/*-------------------------------------------------------------------
* Only called by process 0
*/
void Print_list(double local_A[], int local_n, int rank) {
	int i;
	printf("%d: ", rank);
	for (i = 0; i < local_n; i++)
		printf("%f ", local_A[i]);
	printf("\n");
}  /* Print_list */

/*-------------------------------------------------------------------
* Function:   Print_local_lists
* Purpose:    Print each process' current list contents
* Input args: all
* Notes:
* 1.  Assumes all participating processes are contributing local_n
*     elements
*/
void Print_local_lists(double local_A[], int local_n,
	int my_rank, int p, MPI_Comm comm) {
	double*       A;
	int        q;
	MPI_Status status;

	if (my_rank == 0) {
		A = (double*)malloc(local_n*sizeof(double));
		Print_list(local_A, local_n, my_rank);
		for (q = 1; q < p; q++) {
			MPI_Recv(A, local_n, MPI_DOUBLE, q, 0, comm, &status);
			Print_list(A, local_n, q);
		}
		free(A);
	}
	else {
		MPI_Send(local_A, local_n, MPI_DOUBLE, 0, 0, comm);
	}
}  /* Print_local_lists */

/*-------------------------------------------------------------------
* Function:   Print_seri_lists
* Purpose:    Print  serial list contents
* Input args: all
* 
*/
void Print_seri_list(double local_A[], int local_n){
	local_n = COUNT;
	int i;
	for(i = 0;i<COUNT;i++){
		printf("%f",local_A[i]);
		printf(" ");
	}
	printf("\n");
}
