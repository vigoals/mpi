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
#define GET_TIME(now) { \
   struct timeval t; \
   gettimeofday(&t, NULL); \
   now = t.tv_sec + t.tv_usec/1000000.0; \
}

typedef double mtyp;

int n0, m, n1;
int task_id;
int task_max;
MPI_Comm comm;
mtyp * matrix2;
mtyp * matrix1;
mtyp * ans_mat;
mtyp * ans_mat_t;
int my_begin_num;
int my_ans_size;
int * begin_nums;
int * ans_size;

void parallel_mm();
void mm();
void random_init();
void fill_matrix_by_random(mtyp * a, int count);
void printf_matrix(mtyp *a, int n, int m);
void compute_task_begin_num();

int main(int argc, char* argv[]) {
	int i, j, k;
	double start_time, end_time;

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

	matrix1 = (mtyp *)malloc(n0*m*sizeof(mtyp));
	matrix2 = (mtyp *)malloc(m*n1*sizeof(mtyp));
	if (task_id == 0) {
		random_init();
		printf("Generate matrix1(%d*%d) and matrix2(%d*%d).\n", n0, m, m, n1);
		fill_matrix_by_random(matrix1, n0*m);
		fill_matrix_by_random(matrix2, m*n1);
		//printf("Matrix1:\n");
		//printf_matrix(matrix1, n0, m);
		//printf("Matrix2:\n");
		//printf_matrix(matrix2, m, n1);
		
		printf("Send matrix1 and matrix2.\n");
		GET_TIME(start_time);
		compute_task_begin_num();
		my_begin_num = begin_nums[0];
		my_ans_size = ans_size[0];
		for (i=1; i<task_max; i++) {
			MPI_Send(matrix1, n0*m, MPI_DOUBLE, i, MATRIX1_TAG, comm);
			MPI_Send(matrix2, m*n1, MPI_DOUBLE, i, MATRIX2_TAG, comm);
			MPI_Send(&begin_nums[i], 1, MPI_INT, i, BEGIN_NUM_TAG, comm);
			MPI_Send(&ans_size[i], 1, MPI_INT, i, ANS_SIZE_TAG, comm);
		}
	} else {
		MPI_Recv(matrix1, n0*m, MPI_DOUBLE, 0, MATRIX1_TAG, comm, MPI_STATUS_IGNORE);
		MPI_Recv(matrix2, m*n1, MPI_DOUBLE, 0, MATRIX2_TAG, comm, MPI_STATUS_IGNORE);
		MPI_Recv(&my_begin_num, 1, MPI_INT, 0, BEGIN_NUM_TAG, comm, MPI_STATUS_IGNORE);
		MPI_Recv(&my_ans_size, 1, MPI_INT, 0, ANS_SIZE_TAG, comm, MPI_STATUS_IGNORE);
	}

	//parallel_mm();

	if (task_id == 0) {
		double t;
		GET_TIME(end_time);
		t = end_time - start_time;
		printf("Time used by parallel mm is %.9fs.\n", t);
		//printf("Ans matrix by parallel mm:\n");
		//printf_matrix(ans_mat, n0, n1);
		mm();
	}

	MPI_Finalize();
	return 0;
}

void random_init() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	srand((unsigned int)tv.tv_usec);
}

void fill_matrix_by_random(mtyp *a, int count) {
	int i;
	for (i=0; i<count; i++)
		a[i] = (rand()/((double)RAND_MAX)*2 - 1)*GENERATE_MAX;
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
void printf_matrix(mtyp *a, int n, int m) {
	int i,j;
	for (i=0; i<n; i++) {
		for (j=0; j<m; j++) {
			printf("%.5f", a[i*m + j]);
			if (j != m - 1)
				printf(", ");
		}
		printf("\n");
	}
	printf("[Size of matrix is %d*%d]\n", n, m);
}

void parallel_mm() {
	int i, j, k, an, l, c;
	mtyp * ans = (mtyp *)malloc(my_ans_size*sizeof(mtyp));

	if (task_id == 0)
		printf("Begin Parallel mm.\n");

	for (i=my_begin_num; i<my_begin_num + my_ans_size; i++) {
		l = i/n1;
		c = i%n1;

		ans[i - my_begin_num] = 0;
		for (j=0; j<m; j++)
			ans[i - my_begin_num] += matrix1[l*m + j]*matrix2[j*n1 + c];
		//printf("Task %d process ans %d(%d, %d), ans is %d.\n", task_id, i, l, c, ans[i - my_begin_num]);

		/*if (task_id ==0) {
			printf("\r[Run %d of %d.]", i - my_begin_num, my_ans_size);
		}*/
	}

	if (task_id == 0) {
		printf("\n");
		ans_mat = (mtyp *)malloc(n0*n1*sizeof(mtyp));
		memcpy(ans_mat, ans, my_ans_size*sizeof(mtyp));

		for (i=1; i<task_max; i++) {
			MPI_Recv(ans_mat + begin_nums[i], ans_size[i], MPI_DOUBLE, i, ANS_TAG, comm, MPI_STATUS_IGNORE);
		}
	} else {
		MPI_Send(ans, my_ans_size, MPI_DOUBLE, 0, ANS_TAG, comm);
	}

}

void mm() {
	int i, j, k;
	double start, end;
	int p_all, p_go;
	ans_mat_t = (mtyp *)malloc(n0*n1*sizeof(mtyp));

	printf("Begin mm.\n");
	GET_TIME(start);
	p_all = n0*n1;
	p_go = 0;
	for (i=0; i<n0; i++) {
		for (j=0; j<n1; j++) {
			ans_mat_t[i*n1 + j] = 0;
			for (k=0; k<m; k++)
				ans_mat_t[i*n1 + j] += matrix1[i*m + k]*matrix2[k*n1 + j];
		}
		p_go += n1;
		//printf("\r[%2.2f%%]", p_go/((double)p_all)*100);
	}
	GET_TIME(end);
	printf("\nTime used by mm: %.9f.\n", end - start);
	//printf("Ans matrix by mm:\n");
	//printf_matrix(ans_mat_t, n0, n1);
}

