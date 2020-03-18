#include <stdlib.h>
#include <stdio.h>
#include "io.h"

#define RAND01 ((double)rand() / (double)RAND_MAX)

void createMatrix(double*** _mz, int num_rows, int num_columns);
void randomFillLR(int nU, int nI, int nF);
void copyMatrix();
void multMatrices();
void factorization();
void swapMatrices();

int max_iterations = 0;
double alpha = 0;
int num_fs = 0;
int num_l = 0;
int num_c = 0;
int non_zero_entries = 0;

double** mz_a = NULL;
double** mz_l = NULL;
double** mz_r = NULL;
double** mz_b = NULL;
double** mz_l_prev = NULL;
double** mz_r_prev = NULL;

int main(int argc, char* argv[]){
    FILE* fp = NULL;
    int lin = 0;
    int col = 0;
    double value = 0;

    if(argc != 2){
        fprintf(stderr, "Erro nos argumentos");
        exit(-1);
    }
    
    fp = openFile(argv[1]);
    
    fscanf(fp, "%d\n", &max_iterations);
    fscanf(fp, "%lf\n", &alpha);
    fscanf(fp, "%d\n", &num_fs);
    fscanf(fp, "%d %d %d", &num_l, &num_c, &non_zero_entries);
       
    createMatrix(&mz_a, num_l, num_c);//Create matrix A
    createMatrix(&mz_l, num_l, num_fs);//Create matrix L
    createMatrix(&mz_r, num_fs, num_c);//Create matrix R
    createMatrix(&mz_b, num_l, num_c);//Create matrix B
    createMatrix(&mz_l_prev, num_l, num_fs);//Create matrix L_p
    createMatrix(&mz_r_prev, num_fs, num_c);//Create matrix R_p
    
    while(!feof(fp) && fscanf(fp, "%d %d %lf" , &lin, &col, &value)){
        mz_a[lin][col] = value;
    }
    randomFillLR(num_l, num_c, num_fs);//L and R ramdom fill in. Copy to L_p e R_p
    multMatrices();//B initial
    //------------------------------------------------
            // fprintf(stdout, "A\n");
            // for(int i = 0; i < num_l; i++){
                // for(int j = 0; j < num_c; j++){
                    // fprintf(stdout, "%.3f ", mz_a[i][j]);
                // }
                // fprintf(stdout, "\n");
            // }
            // fprintf(stdout, "\n");
            
            
    for(int a = 0; a < max_iterations; a++){  
        factorization();
    }
    fprintf(stdout, "B\n");
            for(int i = 0; i < num_l; i++){
                for(int j = 0; j < num_c; j++){
                    fprintf(stdout, "%f ", mz_b[i][j]);
                }
                fprintf(stdout, "\n");
            }
            fprintf(stdout, "\n");
    fclose(fp);
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
    for(int i = 0; i < num_rows; i++){
        for(int j = 0; j < num_columns; j++){
            (*_mz)[i][j] = 0;
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
            mz_r[i][j] = RAND01 / (double) nF;
    }
    copyMatrix();
}

void copyMatrix(){
    for(int i = 0; i < num_l; i++){
        for(int j = 0; j < num_fs; j++){
            mz_l_prev[i][j] = mz_l[i][j];
        }
    }
    for(int i = 0; i < num_fs; i++){
        for(int j = 0; j < num_c; j++){
            mz_r_prev[i][j] = mz_r[i][j];
        }
    }
}

void multMatrices(){
    double sum = 0;
    for (int e = 0; e < num_l; e++) {
        for (int d = 0; d < num_c; d++) {
            for (int k = 0; k < num_fs; k++) {
                sum = sum + (mz_l[e][k])*(mz_r[k][d]);
            }
            mz_b[e][d] = sum;
            sum = 0;
        }
    }  
}

void swapMatrices(){
    double** backup_l = mz_l_prev;
    double** backup_r = mz_r_prev;
    mz_l_prev=mz_l;
    mz_r_prev=mz_r;
    mz_l=backup_l;
    mz_r=backup_r; 
}

void factorization(){
    double sum_l=0;
    double sum_r=0;
    for(int i = 0; i < num_l; i++){
        for(int j = 0; j < num_c; j++){
            if (mz_a[i][j] != 0){
                for(int k = 0; k < num_fs; k++){
                    for(int n = 0; n < num_c; n++){
                        if (mz_a[i][n] != 0){
                            sum_l = sum_l + 2 * (mz_a[i][n] - mz_b[i][n]) * (-1 * mz_r_prev[k][n]);
                        }
                    }
                    for(int n = 0; n < num_l; n++){
                        if (mz_a[n][j] != 0){
                            sum_r = sum_r + 2 * (mz_a[n][j] - mz_b[n][j]) * (-1 * mz_l_prev[n][k]);
                        }                            
                    }
                    mz_l[i][k] = mz_l_prev[i][k] - alpha * sum_l;
                    sum_l=0;
                    mz_r[k][j] = mz_r_prev[k][j] - alpha * sum_r;
                    sum_r=0;
                }
            }
        }
    }
    multMatrices();
    //---------------------------------------------------------------------------------------------------
    // fprintf(stdout, "L\n");
    // for(int i = 0; i < num_l; i++){
    // for(int j = 0; j < num_fs; j++){
            // fprintf(stdout, "%f ", mz_l[i][j]);
        // }
        // fprintf(stdout, "\n");
    // }
    // fprintf(stdout, "\n");
    // fprintf(stdout, "R\n");
    // for(int i = 0; i < num_fs; i++){
        // for(int j = 0; j < num_c; j++){
            // fprintf(stdout, "%f ", mz_r[i][j]);
        // }
        // fprintf(stdout, "\n");
    // }
    // fprintf(stdout, "\n");
    // fprintf(stdout, "B\n");
    // for(int i = 0; i < num_l; i++){
        // for(int j = 0; j < num_c; j++){
            // fprintf(stdout, "%f ", mz_b[i][j]);
        // }
        // fprintf(stdout, "\n");
    // }
    // fprintf(stdout, "\n");
    // fprintf(stdout, "L_p\n");
    // for(int i = 0; i < num_l; i++){
        // for(int j = 0; j < num_fs; j++){
            // fprintf(stdout, "%f ", mz_l_prev[i][j]);
        // }
        // fprintf(stdout, "\n");
    // }
    // fprintf(stdout, "\n");
    // fprintf(stdout, "R_p\n");
    // for(int i = 0; i < num_fs; i++){
        // for(int j = 0; j < num_c; j++){
            // fprintf(stdout, "%f ", mz_r_prev[i][j]);
        // }
        // fprintf(stdout, "\n");
    // }
    // fprintf(stdout, "\n");
    //-----------------------------------------------------------------------------------------------------------    
        
    swapMatrices();
}

