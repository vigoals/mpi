#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>
#include <string.h>

#define MATRIX1_TAG 1
#define MATRIX2_TAG 2
#define ANS_TAG 3
#define BEGIN_NUM_TAG 4
#define ANS_SIZE_TAG 5
#define GENERATE_MAX 100

int n0, m, n1;
int task_id;
int task_max;
MPI_Comm comm;
int * matrix2;
int * matrix1;
int * ans_mat;
int my_begin_num;
int my_ans_size;
int * begin_nums;
int * ans_size;

void parallel_mm();
void random_init();
void fill_matrix_by_random(int * a, int count);
void printf_matrix(int *a, int n, int m);
void compute_task_begin_num();

int main(int argc, char* argv[]) {
	int i, j, k;

	MPI_Init(NULL, NULL);
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &task_max);
	MPI_Comm_rank(comm, &task_id);
	//printf("%d %d\n", argc, task_id);

	if (argc != 4) {
		if (task_id == 0) {
			printf("usage:mpirun -n <task_num> ./main <n0> <m> <n1>");
		}
		//printf("exit err %d\n", task_id);
		MPI_Finalize();
		exit(-1);
	}
	n0 = atoi(argv[1]);
	m = atoi(argv[2]);
	n1 = atoi(argv[3]);

	matrix1 = (int *)malloc(n0*m*sizeof(int));
	matrix2 = (int *)malloc(m*n1*sizeof(int));
	if (task_id == 0) {
		random_init();
		printf("Generate and send matrix2 %d*%d.\n", m, n1);
		fill_matrix_by_random(matrix1, n0*m);
		fill_matrix_by_random(matrix2, m*n1);
		//printf_matrix(matrix1, n0, m);
		//printf_matrix(matrix2, m, n1);
		compute_task_begin_num();

		my_begin_num = begin_nums[0];
		my_ans_size = ans_size[0];
		for (i=1; i<task_max; i++) {
			MPI_Send(matrix1, n0*m, MPI_INT, i, MATRIX1_TAG, comm);
			MPI_Send(matrix2, m*n1, MPI_INT, i, MATRIX2_TAG, comm);
			MPI_Send(&begin_nums[i], 1, MPI_INT, i, BEGIN_NUM_TAG, comm);
			MPI_Send(&ans_size[i], 1, MPI_INT, i, ANS_SIZE_TAG, comm);
		}
	} else {
		MPI_Recv(matrix1, n0*m, MPI_INT, 0, MATRIX1_TAG, comm, MPI_STATUS_IGNORE);
		MPI_Recv(matrix2, m*n1, MPI_INT, 0, MATRIX2_TAG, comm, MPI_STATUS_IGNORE);
		MPI_Recv(&my_begin_num, 1, MPI_INT, 0, BEGIN_NUM_TAG, comm, MPI_STATUS_IGNORE);
		MPI_Recv(&my_ans_size, 1, MPI_INT, 0, ANS_SIZE_TAG, comm, MPI_STATUS_IGNORE);
	}

	parallel_mm();

	if (task_id == 0) {
		//printf_matrix(ans_mat, n0, n1);
	}

	MPI_Finalize();
	return 0;
}

void random_init() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	srandom((unsigned int)tv.tv_usec);
}

void fill_matrix_by_random(int *a, int count) {
	int i;
	for (i=0; i<count; i++)
		a[i] = random() % GENERATE_MAX;
}


void compute_task_begin_num() {
	int i, k, t;
	k = 0;
	begin_nums = (int *)malloc(task_max*sizeof(int));
	ans_size = (int *)malloc(task_max*sizeof(int));
	for (i=0; i<task_max; i++) {
		t = n0*n1/task_max;
		if (i < (n0*n1)%task_max)
			t++;
		begin_nums[i] = k;
		ans_size[i] = t;
		k += t;
	}
}

//未考虑并行
void printf_matrix(int *a, int n, int m) {
	int i,j;
	for (i=0; i<n; i++) {
		for (j=0; j<m; j++)
			printf("%5d", a[i*m + j]);
		printf("\n");
	}
	printf("[Size of matrix is %d*%d]\n", n, m);
}

void parallel_mm() {
	int i, j, k, an, l, c;
	int * ans = (int *)malloc(my_ans_size*sizeof(int));

	for (i=my_begin_num; i<my_begin_num + my_ans_size; i++) {
		l = i/n1;
		c = i%n1;

		ans[i - my_begin_num] = 0;
		for (j=0; j<m; j++)
			ans[i - my_begin_num] += matrix1[l*m + j]*matrix2[j*n1 + c];
		//printf("Task %d process ans %d(%d, %d), ans is %d.\n", task_id, i, l, c, ans[i - my_begin_num]);
	}

	if (task_id == 0) {
		ans_mat = (int *)malloc(n0*n1*sizeof(int));
		memcpy(ans_mat, ans, my_ans_size*sizeof(int));

		for (i=1; i<task_max; i++) {
			MPI_Recv(ans_mat + begin_nums[i], ans_size[i], MPI_INT, i, ANS_TAG, comm, MPI_STATUS_IGNORE);
		}
	} else {
		MPI_Send(ans, my_ans_size, MPI_INT, 0, ANS_TAG, comm);
	}

	/*
	int i, j, k, ii, l;
	int now = task_id;
	int * matrix1 = (int *)malloc(m*sizeof(int));
	int * ans;
	int ans_size = n0/task_max;
	if (task_id < n0%task_max)
		ans_size++;
	printf("Ans size of task %d is %d.\n", task_id, ans_size);
	ans = (int *)malloc(n1*sizeof(int)*ans_size);

	//srandom((unsigned int)time(0));
	for (i=0; i<ans_size; i++) {
		for (j=0; j<m; j++)
			matrix1[j] = random() % GENERATE_MAX;

		//l = i*task_max + task_id;
		//printf("Task %d process line %d.\n", task_id, l);

		for (ii=0; ii<n1; ii++) {
			k = 0;
			for (j=0; j<m; j++)
				k += matrix1[j]*matrix2[j*n1 + ii];
			ans[i*n1 + ii] = k;
		}
	}

	if (task_id != 0) {
		MPI_Send(ans, n1*ans_size, MPI_INT, 0, ANS_TAG, comm);
	} else { //recv all answer
		int * ans_recv;
		int as = n0/task_max;
		int * ans_all;

		for (i=1; i<task_max; i++) {
			int as_temp = as;
			if (i < n0%task_max)
				as_temp++;
			ans_recv = (int *)malloc(n1*as_temp*sizeof(int));
			MPI_Recv(ans, n1*as_temp, MPI_INT, i, ANS_TAG, comm, MPI_STATUS_IGNORE);
			printf("Recv %d answers from %d.\n", as_temp, i);
			free(ans_recv);
		}
	}
	free(ans);
	*/
}
