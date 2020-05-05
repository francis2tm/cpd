#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <mpi.h>
#include <omp.h>

#define RAND01 ((double)rand() / (double)RAND_MAX)

#define BLOCK_LOW(id,p,n)      ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n)     (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n)     (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(index,p,n) (((p)*((index)+1)-1)/(n))

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

int main(int argc, char* argv[]){
    
	MPI_Status status;
    int proc_id, num_procs, name_procs_len, num_threads_proc, thread_id, provided_mpi_support, my_num_nz;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	double secs;
	
    MPI_Init_thread (&argc, &argv,MPI_THREAD_FUNNELED,&provided_mpi_support);
	if(provided_mpi_support != MPI_THREAD_FUNNELED){//aborted due to failed support for the omp implementation
		MPI_Finalize();
		exit (-3);
	}
    MPI_Comm_rank (MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
	MPI_Get_processor_name(processor_name, &name_procs_len);
	
    if(argc != 3){
		if (!proc_id){
			printf ("Command line: %s <nº threads per processor> <instance>.in \n", argv[0]);
		}
		MPI_Finalize();
		exit (-2);
    }
	
    MPI_Barrier (MPI_COMM_WORLD);
    secs = - MPI_Wtime();
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
	

	
	MPI_Barrier (MPI_COMM_WORLD);
    secs += MPI_Wtime();
	printf("Proc: %d :Ini stuff: %d %d %d %d %d \n",proc_id,initial_info[0],initial_info[1],initial_info[2],initial_info[3],initial_info[4]);//------------------------------------------------
	fflush(stdout);
	printf("Proc: %d :my_nz = %d \n",proc_id,my_num_nz);//------------------------------------------------
	fflush(stdout);
	for (int i = 0;i<my_num_nz;i++){
		printf("Proc: %d :A_nz: %d %d %f \n",proc_id,non_zeros_pos[2*i],non_zeros_pos[2*i+1],non_zeros_values[i]);//------------------------------------------------
		fflush(stdout);
	}
	printf("Proc: %d :alpha %f \n",proc_id,alpha);//------------------------------------------------
	fflush(stdout);
	// if(!proc_id){
		// printf("GG nos  Time = %12.6f sec  \n",secs);//--------------------------------------------------------------------------------------------------------------------------
		// fflush(stdout);
	// }
	MPI_Finalize();
    return 0;
}