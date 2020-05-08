#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define RAND01 ((double)rand() / (double)RAND_MAX)

#define BLOCK_LOW(id,p,n)      ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n)     (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n)     (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(index,p,n) (((p)*((index)+1)-1)/(n))

void createMatrix(double*** _mz, int num_rows, int num_columns);

/*	
	max_iterations
	num_fs
	num_l
	num_c
	non_zero_entries
*/
int initial_info[5]={0,0,0,0,0};

double alpha = 0;

int* non_zeros_pos = NULL;
double* non_zeros_values = NULL;
double** mz_l = NULL;
double** mz_r = NULL;
double* b_non_zeros = NULL;
double** mz_l_sum = NULL;
double** mz_r_sum = NULL;
int* r_indices = NULL;
int* cols_copy = NULL;

int main(int argc, char* argv[]){
	
    int proc_id, num_procs, name_procs_len, num_threads_proc, thread_id, provided_mpi_support, my_num_nz, ini_line, fin_line, count_copy_cols;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	double secs, lixo;
	
    MPI_Init_thread (&argc, &argv,MPI_THREAD_FUNNELED,&provided_mpi_support);
	if(provided_mpi_support != MPI_THREAD_FUNNELED){//aborted due to failed support for the omp implementation
		MPI_Finalize();
		exit (-3);
	}
	
    if(argc != 3){
		if (!proc_id){
			printf ("Command line: %s <num threads per processor> <instance>.in \n", argv[0]);
		}
		MPI_Finalize();
		exit (-2);
    }
	num_threads_proc = atoi(argv[1]);

	MPI_Comm_rank (MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
	MPI_Get_processor_name(processor_name, &name_procs_len);

	MPI_Barrier (MPI_COMM_WORLD);
    secs = - MPI_Wtime();
	srand(0);
	// printf("Proc: %d here \n",proc_id);//------------------------------------------------
	// fflush(stdout);
	if(!proc_id){
		FILE* fp = NULL;
		fp = fopen(argv[2], "r");
		if(fp == NULL){
			fprintf(stderr, "Erro abrir ficheiro");
			MPI_Finalize();
			exit(-1);
		}
		
		fscanf(fp, "%d\n", &initial_info[0]);
		fscanf(fp, "%lf\n", &alpha);
		fscanf(fp, "%d\n", &initial_info[1]);
		fscanf(fp, "%d %d %d", &initial_info[2], &initial_info[3], &initial_info[4]);
	
		MPI_Bcast (initial_info, 5, MPI_INT, 0, MPI_COMM_WORLD);//share initial information, so that each PC might build their L and R transpose
		
		int* local_buf_non_zeros_pos = (int*) malloc ((BLOCK_SIZE(num_procs-1,num_procs,initial_info[4]) + initial_info[2]) * 2 * sizeof(int));
		double* local_buf_non_zeros_values = (double*) malloc ((BLOCK_SIZE(num_procs-1,num_procs,initial_info[4]) + initial_info[2]) * sizeof(double));
		int lin = 0;
		int col = 0;
		double value = 0;
		int buf_write_pos = 0;
		int cur_line_owner = num_procs-1;
		int prev_line_owner = -1;
		fscanf(fp, "%d %d %lf" , &local_buf_non_zeros_pos[2 * buf_write_pos], &local_buf_non_zeros_pos[2 * buf_write_pos + 1], &local_buf_non_zeros_values[buf_write_pos]);
		for(int i = initial_info[4]-2; i >= 0; i--){
			fscanf(fp, "%d %d %lf" , &lin, &col, &value);
			if(local_buf_non_zeros_pos[2 * (buf_write_pos)] != lin){//line change
				prev_line_owner = cur_line_owner;
				cur_line_owner = BLOCK_OWNER(i,num_procs,initial_info[4]);
				if(prev_line_owner != cur_line_owner){//non zeros per process were reached, the extra nz in the same line are also sent, leaving less work for the master that has to coordenate things
					buf_write_pos++;
					MPI_Send(&buf_write_pos, 1, MPI_INT, prev_line_owner, prev_line_owner, MPI_COMM_WORLD);
					local_buf_non_zeros_values[buf_write_pos] = alpha;
					MPI_Send(local_buf_non_zeros_pos, 2 * buf_write_pos, MPI_INT, prev_line_owner, prev_line_owner, MPI_COMM_WORLD);
					MPI_Send(local_buf_non_zeros_values, buf_write_pos + 1, MPI_DOUBLE, prev_line_owner, prev_line_owner, MPI_COMM_WORLD);
					buf_write_pos = -1;
				}
			}
			buf_write_pos++;
			local_buf_non_zeros_pos[2 * buf_write_pos] = lin;
			local_buf_non_zeros_pos[2 * buf_write_pos + 1] = col;
			local_buf_non_zeros_values[buf_write_pos] = value;
		}
		my_num_nz = buf_write_pos + 1;
		non_zeros_pos = local_buf_non_zeros_pos;//what remains are the nz for the master
		non_zeros_values = local_buf_non_zeros_values;//what remains are the nz for the master	
	}
	else{
		MPI_Status status;
		MPI_Bcast (initial_info, 5, MPI_INT, 0, MPI_COMM_WORLD);//receive initial information
		MPI_Recv(&my_num_nz, 1, MPI_INT, 0, proc_id, MPI_COMM_WORLD, &status);
		non_zeros_pos = (int*) malloc (2 * my_num_nz * sizeof(int));
		non_zeros_values = (double*) malloc ((my_num_nz + 1) * sizeof(double));//+1 is for alpha, that being a double ill be sent along with the values
		MPI_Recv(non_zeros_pos, 2 * my_num_nz, MPI_INT, 0, proc_id, MPI_COMM_WORLD, &status);
		MPI_Recv(non_zeros_values, my_num_nz + 1, MPI_DOUBLE, 0, proc_id, MPI_COMM_WORLD, &status);
		alpha = non_zeros_values[my_num_nz];
	}
	
	//L and L_sum size
	
	ini_line = non_zeros_pos[0];
	fin_line = non_zeros_pos[2 * (my_num_nz - 1)];
	
	//R_mini size
	
	r_indices = (int*) malloc (my_num_nz * sizeof(int));
	cols_copy = (int*) malloc (my_num_nz * sizeof(int));
	cols_copy[0] = non_zeros_pos[1];
	r_indices[0] = 0;
	count_copy_cols = 1;
	int j = 0;
	for (int i = 1; i < my_num_nz; i++){
		for (j = 0; j < count_copy_cols; j++){//check if it is a repeted col
			if(non_zeros_pos[2 * i + 1] == cols_copy[j]){
				r_indices[i] = j;
				j++;
				break;
			}			
		}
		if(j == count_copy_cols){//add new element
			cols_copy[count_copy_cols] = non_zeros_pos[2 * i + 1];
			r_indices[i] = count_copy_cols;
			count_copy_cols++;
		}
	}	
	
	
	
	createMatrix(&mz_l, (fin_line - ini_line + 1), initial_info[1]);//Create matrix L
    createMatrix(&mz_r, initial_info[3], initial_info[1]);//Create matrix R transpose
    b_non_zeros = (double*) malloc (my_num_nz * sizeof(double));//create B non zeros(correspondent to A's non_zeros)
    createMatrix(&mz_l_sum, (fin_line - ini_line + 1), initial_info[1]);//Create matrix L_sum (aka prev)
	createMatrix(&mz_r_sum, count_copy_cols, initial_info[1]);//Create matrix R_sum (aka mini) without ordered cols
	
	//randomFillLR
	//each one creates their L with the needed values	
	for(int i = 0; i < ini_line; i++){
        for(int j = 0; j < initial_info[1]; j++){
            lixo = RAND01 / (double) initial_info[1];
		}
    }
    for(int i = 0; i <= (fin_line-ini_line); i++){
        for(int j = 0; j < initial_info[1]; j++){
            mz_l[i][j] = RAND01 / (double) initial_info[1];
			// mz_l_sum[i][j] = mz_l[i][j];
		}
    } 
	for(int i = fin_line + 1; i < initial_info[2]; i++){
        for(int j = 0; j < initial_info[1]; j++){
            lixo = RAND01 / (double) initial_info[1];
		}
    }

	for(int i = 0; i < initial_info[1]; i++){//every one has a full R for the reduce
        for(int j = 0; j < initial_info[3]; j++){
            mz_r[j][i] = RAND01 / (double) initial_info[1];
		}
    }

// #pragma omp parallel for num_threads(2)
	// for(int i = 0; i < count_copy_cols; i++){
        // for(int j = 0; j < initial_info[1]; j++){
            // mz_r_sum[i][j] = mz_r[cols_copy[i]][j];
			// // num_threads_proc = omp_get_num_threads();
			// // thread_id = omp_get_thread_num();
			// // printf("Hello from thread %d out of %d from process %d out of %d on %s\n",thread_id, num_threads_proc, proc_id, num_procs, processor_name);
        // }
    // }

	double sum = 0;
#pragma omp parallel for firstprivate(sum) num_threads(num_threads_proc)
	for(int i = 0; i < my_num_nz; i++){
		for (int k = 0; k < initial_info[1]; k++) {
			mz_l_sum[non_zeros_pos[2 * i] - ini_line][k] = mz_l[non_zeros_pos[2 * i] - ini_line][k];
			mz_r_sum[r_indices[i]][k] = mz_r[cols_copy[r_indices[i]]][k];
			sum += (mz_l[non_zeros_pos[2 * i] - ini_line][k])*(mz_r[cols_copy[r_indices[i]]][k]);
		}
		b_non_zeros[i] = sum;
		sum = 0;
		// num_threads_proc = omp_get_num_threads();
		// thread_id = omp_get_thread_num();
		// printf("Hello from thread %d out of %d from process %d out of %d on %s\n",thread_id, num_threads_proc, proc_id, num_procs, processor_name);
	}
	
	MPI_Barrier (MPI_COMM_WORLD);
    secs += MPI_Wtime();
	// printf("Proc: %d :Ini stuff: %d %d %d %d %d \n",proc_id,initial_info[0],initial_info[1],initial_info[2],initial_info[3],initial_info[4]);//------------------------------------------------
	// fflush(stdout);
	// printf("Proc: %d :my_nz = %d \n",proc_id,my_num_nz);//------------------------------------------------
	// fflush(stdout);
	// for (int i = 0;i<my_num_nz;i++){
		// printf("Proc: %d :A_nz: %d %d %f \n",proc_id,non_zeros_pos[2*i],non_zeros_pos[2*i+1],non_zeros_values[i]);//------------------------------------------------
		// fflush(stdout);
	// }
	// printf("Proc: %d :alpha %f \n",proc_id,alpha);//------------------------------------------------
	// fflush(stdout);
	if(!proc_id){
		for (int i = 0;i<my_num_nz;i++){
			printf("Proc: %d :A_nz: %d %d %f B_nz: %f\n",proc_id,non_zeros_pos[2*i],non_zeros_pos[2*i+1],non_zeros_values[i],b_non_zeros[i]);//------------------------------------------------
			fflush(stdout);
		}
		
		for (int i = 0;i<my_num_nz;i++){
			printf("Proc: %d R_idx: %d \n",proc_id,r_indices[i]);//------------------------------------------------
			fflush(stdout);
		}
		
		
		printf("Proc: %d : counts_copy = %d fin_lin-ini_lin = %d \n",proc_id,count_copy_cols,fin_line - ini_line);//------------------------------------------------
		fflush(stdout);
			
		for (int i = 0;i<count_copy_cols;i++){
			printf("Proc: %d cols_copy: %d \n",proc_id,cols_copy[i]);//------------------------------------------------
			fflush(stdout);
		}
		
		
		
		for(int i = 0; i < (fin_line - ini_line + 1); i++){
			for(int j = 0; j < initial_info[1]; j++){
				printf("Proc: %d : L_ele: %d %d -> %f \n", proc_id,i,j,mz_l[i][j]);
				fflush(stdout);
			}
		}

		for(int i = 0; i < initial_info[3]; i++){
			for(int j = 0; j < initial_info[1]; j++){
				printf("Proc: %d : R_ele: %d %d -> %f \n", proc_id,i,j,mz_r[i][j]);
				fflush(stdout);
			}
		}
		
		for(int i = 0; i < (fin_line - ini_line + 1); i++){
			for(int j = 0; j < initial_info[1]; j++){
				printf("Proc: %d : L_sum_ele: %d %d -> %f \n", proc_id,i,j,mz_l_sum[i][j]);
				fflush(stdout);
			}
		}

		for(int i = 0; i < count_copy_cols; i++){
			for(int j = 0; j < initial_info[1]; j++){
				printf("Proc: %d : R_sum_ele: %d %d -> %f \n", proc_id,i,j,mz_r_sum[i][j]);
				fflush(stdout);
			}
		}
	}
	// if(!proc_id){
		// printf("GG nos  Time = %12.6f sec  \n",secs);//--------------------------------------------------------------------------------------------------------------------------
		// fflush(stdout);
	// }
	MPI_Finalize();
    return 0;
}

void createMatrix(double*** _mz, int num_rows, int num_columns){
    *_mz = (double**)malloc(num_rows * sizeof(double*));
    if(*_mz == NULL){
        fprintf(stderr, "Erro alocar memoria - linhas\n");
        exit(-1);
    }
    for(int a = 0; a < num_rows; a++){
        (*_mz)[a] = (double*)malloc(num_columns * sizeof(double));
        if((*_mz)[a] == NULL){
            fprintf(stderr, "Erro alocar memoria - colunas\n");
            exit(-1);
        }
    }   
}