#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define RAND01 ((double)rand() / (double)RAND_MAX)

#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(index,p,n) (((p)*((index)+1)-1)/(n))


void createMatrix(double*** _mz, int num_rows, int num_columns);
void randomFillLR(int nU, int nI, int nF);
void multMatrices_final();
void factorization();
void result();

int max_iterations = 0;
double alpha = 0;
int num_fs = 0;
int num_l = 0;
int num_c = 0;
int non_zero_entries = 0;
double** mz_l = NULL;
double** mz_r = NULL;
double** mz_b = NULL;
double** mz_l_sum = NULL;
double** mz_r_sum = NULL;

typedef struct Non_zeros{ 
    int x, y;
    double val;
} Vect;

Vect* mz_a2 = NULL;

int main(int argc, char* argv[]){
    
	MPI_Status status;
    int id, p,
	i, rounds;
    double secs;

    MPI_Init (&argc, &argv);

    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    if(argc != 2){
	if (!id) printf ("Command line: %s <instance>.in \n", argv[0]);
	MPI_Finalize();
	exit (-1);
    }

    MPI_Barrier (MPI_COMM_WORLD);
    secs = - MPI_Wtime();
	
	FILE* fp = NULL;
    int lin = 0;
    int col = 0;
    double value = 0;
    int count = 0;
    
	fp = fopen(argv[1], "r");
	if(fp == NULL){
		fprintf(stderr, "Erro abrir ficheiro");
		exit(-1);
	}

	fscanf(fp, "%d\n", &max_iterations);
	fscanf(fp, "%lf\n", &alpha);
	fscanf(fp, "%d\n", &num_fs);
	fscanf(fp, "%d %d %d", &num_l, &num_c, &non_zero_entries);

	mz_a2 = (Vect*) malloc(non_zero_entries * sizeof(Vect));//Create matrix A_non_zeros
	if(mz_a2 == NULL){
		fprintf(stderr, "Erro alocar nova matrix\n");
		exit(-1);
	}
    
    createMatrix(&mz_l, num_l, num_fs);//Create matrix L
    createMatrix(&mz_r, num_c, num_fs);//Create matrix R transpose
    createMatrix(&mz_b, num_l, num_c);//Create matrix B
    createMatrix(&mz_l_sum, num_l, num_fs);//Create matrix L_sum (aka prev)
    createMatrix(&mz_r_sum, num_c, num_fs);//Create matrix R_sum transpose (aka prev)
    
    // while(!feof(fp) && fscanf(fp, "%d %d %lf" , &lin, &col, &value)){
        // mz_a2[count].x = lin;
        // mz_a2[count].y = col;
        // mz_a2[count].val = value;
        // count++;
    // }
    randomFillLR(num_l, num_c, num_fs);//L and R_transpose ramdom fill in and copy to L_sum and R_sum.
    
	if(!id){
	    
		fprintf(stdout, "ME:%d  Before ->  L\n",id);
		for(int i = 0; i < num_l; i++){
			for(int j = 0; j < num_fs; j++){
				fprintf(stdout, "%f ", mz_l[i][j]);
			}
			fprintf(stdout, "\n");
		}
		fprintf(stdout, "\n");
		
		MPI_Recv(buf,2 , MPI_DOUBLE, p-1, 999, MPI_COMM_WORLD, &status);
		
		
		
		fprintf(stdout, "ME:%d  Before ->  L\n",id);
		for(int i = 0; i < num_l; i++){
			for(int j = 0; j < num_fs; j++){
				fprintf(stdout, "%f ", mz_l[i][j]);
			}
			fprintf(stdout, "\n");
		}
		fprintf(stdout, "\n");
		
		// fprintf(stdout, "ME:%d  ->  R\n",id);
		// for(int i = 0; i < num_c; i++){
			// for(int j = 0; j < num_fs; j++){
				// fprintf(stdout, "%f ", mz_r[i][j]);
			// }
			// fprintf(stdout, "\n");
		// }
		// fprintf(stdout, "\n");
		
	}
	else{
	    MPI_Send(&i, 1, MPI_INT, (id+1)%p, 999, MPI_COMM_WORLD);
	}
	
	
	
	
	// multMatrices_final();//B initial
    // factorization();
    // multMatrices_final();//B final    
    // result();
	
	
	MPI_Barrier (MPI_COMM_WORLD);
    secs += MPI_Wtime();
	
	fprintf(stderr, "GG nos  Time = %12.6f sec  \n",secs);
	
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

void randomFillLR(int nU, int nI, int nF){
    srand(0);
    for(int i = 0; i < nU; i++){
        for(int j = 0; j < nF; j++)
            mz_l[i][j] = RAND01 / (double) nF;
    }  
    for(int i = 0; i < nF; i++){
        for(int j = 0; j < nI; j++)
            mz_r[j][i] = RAND01 / (double) nF;
    } 
    // for(int i = 0; i < nU; i++){
        // for(int j = 0; j < nF; j++){
            // mz_l_sum[i][j] = mz_l[i][j];
        // }
    // }
    // for(int i = 0; i < nI; i++){
        // for(int j = 0; j < nF; j++){
            // mz_r_sum[i][j] = mz_r[i][j];
        // }
    // }
}


void multMatrices_final(){
    double sum = 0;
    for (int e = 0; e < num_l; e++) {
        for (int d = 0; d < num_c; d++) {
            for (int k = 0; k < num_fs; k++) {
                sum += (mz_l[e][k])*(mz_r[d][k]);
            }
            mz_b[e][d] = sum;
            sum = 0;
        }
    }  
}

void factorization(){
    double aux = 0;
    double sum = 0;  
    for(int count = 0; count < max_iterations; count++){
        for(int i = 0; i < non_zero_entries; i++){
            aux = -1 * alpha * 2 * (mz_a2[i].val - mz_b[mz_a2[i].x][mz_a2[i].y]);
            for(int k = 0; k < num_fs; k++){
                mz_l[mz_a2[i].x][k] += aux * (-1 * mz_r_sum[mz_a2[i].y][k]);
                mz_r[mz_a2[i].y][k] += aux * (-1 * mz_l_sum[mz_a2[i].x][k]);
            }
        }     
        for(int i = 0; i < non_zero_entries; i++){
            for (int k = 0; k < num_fs; k++) {
                mz_l_sum[mz_a2[i].x][k] = mz_l[mz_a2[i].x][k];
                mz_r_sum[mz_a2[i].y][k] = mz_r[mz_a2[i].y][k];
                sum += (mz_l[mz_a2[i].x][k])*(mz_r[mz_a2[i].y][k]);
            }
            mz_b[mz_a2[i].x][mz_a2[i].y] = sum;
            sum = 0;
        }
    }  
}

void result(){
    double max_row = 0;
    int* items = (int*)malloc(num_l * sizeof(int));
    if(items == NULL){
        fprintf(stderr, "Erro alocar vetor result\n");
        exit(-1);
    }
    for(int i = 0; i < non_zero_entries; i++){
        mz_b[mz_a2[i].x][mz_a2[i].y] = 0;
    } 
    for(int i = 0; i < num_l; i++){
        for(int j = 0; j < num_c; j++){
            if(mz_b[i][j] > max_row){
                max_row = mz_b[i][j];
                items[i]=j;
            }
        }
        max_row = 0;
    }
    for(int i = 0; i < num_l; i++){
        fprintf(stdout, "%d\n", items[i]);
    }
}