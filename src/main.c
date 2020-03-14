#include <stdlib.h>
#include <stdio.h>
#include "io.h"

#define RAND01 ((double)rand() / (double)RAND_MAX)

typedef struct Matrix_Struct{
    double** mz;
    int num_l;
    int num_c;
}Matrix;

void createMatrix(Matrix* _mz, int num_rows, int num_columns);
void randomFillLR(int nU, int nI, int nF);
void copyMatrix(Matrix* mz_src, Matrix* mz_dest);
void swapMatrices(Matrix** a, Matrix** b);
void multMatrices(Matrix* c);
void factorization(Matrix* l,Matrix* r);

int max_iterations = 0;
double alpha = 0;
int num_fs = 0;
int num_l = 0;
int num_c = 0;
int non_zero_entries = 0;

Matrix mz_a = {.mz = NULL};
Matrix mz_l =  {.mz = NULL};
Matrix mz_r =  {.mz = NULL};
Matrix mz_b =  {.mz = NULL};
Matrix mz_l_prev =  {.mz = NULL};
Matrix mz_r_prev =  {.mz = NULL};
Matrix *mz_l_ptr = &mz_l;
Matrix *mz_r_ptr = &mz_r;
Matrix *mz_l_prev_ptr = &mz_l_prev;
Matrix *mz_r_prev_ptr = &mz_r_prev;

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
    
    //fprintf(stderr, "%d\n %f\n %d\n %d\n %d\n %d\n",max_iterations,alpha,num_fs,num_l,num_c,non_zero_entries);---------------------------------------------------------------
    
    createMatrix(&mz_a, num_l, num_c);
    createMatrix(&mz_l, num_l, num_fs);
    createMatrix(&mz_r, num_fs, num_c);
    createMatrix(&mz_l_prev, num_l, num_fs);
    createMatrix(&mz_r_prev, num_fs, num_c);
    createMatrix(&mz_b, num_l, num_c);
    
    while(!feof(fp) && fscanf(fp, "%d %d %lf" , &lin, &col, &value)){
        mz_a.mz[lin][col] = value;
    }
    
    //Sรณ imprime---------------------------------------------------------------------------------
    /*
    for(int i = 0; i < num_l; i++){
        for(int j = 0; j < num_c; j++){
            fprintf(stdout, "%.3lf ", mz_a.mz[i][j]);
        }
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
    */
    
    randomFillLR(num_l, num_c, num_fs);
    copyMatrix(&mz_l, &mz_l_prev);
    copyMatrix(&mz_r, &mz_r_prev);
    multMatrices(&mz_b);
    factorization(&mz_l,&mz_r);
    swapMatrices(&mz_l_ptr,&mz_l_prev_ptr);
    swapMatrices(&mz_r_ptr,&mz_r_prev_ptr);
    for(int a = 0; a < 2; a++){//max_iterations---------------------------------------------------------------------------------------------
        multMatrices(&mz_b);
        factorization(&mz_l,&mz_r);
        swapMatrices(&mz_l_ptr,&mz_l_prev_ptr);
        swapMatrices(&mz_r_ptr,&mz_r_prev_ptr);
        //--------------------------------------------------------------------------------------------------
        fprintf(stdout, "A\n");
        for(int i = 0; i < num_l; i++){
            for(int j = 0; j < num_c; j++){
                fprintf(stdout, "%lf ", mz_a.mz[i][j]);
            }
            fprintf(stdout, "\n");
        }
        fprintf(stdout, "\n\n");
        fprintf(stdout, "L\n");
        for(int i = 0; i < num_l; i++){
            for(int j = 0; j < num_fs; j++){
                fprintf(stdout, "%lf ", mz_l.mz[i][j]);
            }
            fprintf(stdout, "\n");
        }
        fprintf(stdout, "\n\n");
        fprintf(stdout, "R\n");    
        for(int i = 0; i < num_fs; i++){
            for(int j = 0; j < num_c; j++){
                fprintf(stdout, "%lf ", mz_r.mz[i][j]);
            }
            fprintf(stdout, "\n");
        }
        fprintf(stdout, "\n\n");
        fprintf(stdout, "B\n");
        for(int i = 0; i < num_l; i++){
            for(int j = 0; j < num_c; j++){
                fprintf(stdout, "%lf ", mz_b.mz[i][j]);
            }
            fprintf(stdout, "\n");
        }
        fprintf(stdout, "\n\n\n");
        fprintf(stdout, "L_p\n");
        for(int i = 0; i < num_l; i++){
            for(int j = 0; j < num_fs; j++){
                fprintf(stdout, "%lf ", mz_l_prev.mz[i][j]);
            }
            fprintf(stdout, "\n");
        }
        fprintf(stdout, "\n\n");
        fprintf(stdout, "R_p\n");
        for(int i = 0; i < num_fs; i++){
            for(int j = 0; j < num_c; j++){
                fprintf(stdout, "%lf ", mz_r_prev.mz[i][j]);
            }
            fprintf(stdout, "\n");
        }
        //--------------------------------------------------------------------------------------------------
    }
    fclose(fp);
    return 0;
}

void createMatrix(Matrix* _mz, int num_rows, int num_columns){

    _mz->num_l = num_rows;
    _mz->num_c = num_columns;

    _mz->mz = (double**)malloc(num_rows * sizeof(double*));
    if(_mz->mz == NULL){
        fprintf(stderr, "Erro alocar memoria");
        exit(-1);
    }
    for(int a = 0; a < num_rows; a++){
        _mz->mz[a] = (double*)malloc(num_columns * sizeof(double));
        if(_mz->mz[a] == NULL){
            fprintf(stderr, "Erro alocar memoria");
            exit(-1);
        }
    }

    for(int i = 0; i < num_rows; i++){
        for(int j = 0; j < num_columns; j++){
            _mz->mz[i][j] = 0;
        }
    }
}

void randomFillLR(int nU, int nI, int nF){
    srand(0);

    for(int i = 0; i < nU; i++){
        for(int j = 0; j < nF; j++)
            mz_l.mz[i][j] = RAND01 / (double) nF;
    }
        
    for(int i = 0; i < nF; i++){
        for(int j = 0; j < nI; j++)
            mz_r.mz[i][j] = RAND01 / (double) nF;
    }
}

void copyMatrix(Matrix* mz_src, Matrix* mz_dest){
    for(int i = 0; i < mz_src->num_l; i++){
        for(int j = 0; j < mz_src->num_c; j++){
            mz_dest->mz[i][j] = mz_src->mz[i][j];
        }
    }
}

void swapMatrices(Matrix** a, Matrix** b){
    Matrix* backup = *a;
    *a = *b;
    *b = backup;
}

void multMatrices(Matrix* c){
    double sum = 0;
    int n=mz_l.num_c;
    int m=mz_l.num_l;
    int q=mz_r.num_c;
    int p=mz_r.num_l;
    if (n != p){
        fprintf(stderr, "Erro multiplacar matrizes L e R");
        exit(-1);
    }else{    
        for (int e = 0; e < m; e++) {
          for (int d = 0; d < q; d++) {
            for (int k = 0; k < p; k++) {
                sum = sum + (mz_l.mz[e][k])*(mz_r.mz[k][d]);
            }
            c->mz[e][d] = sum;
            sum = 0;
          }
        }  
    }
}

void factorization(Matrix* l,Matrix* r){
    double sum_l=0;
    double sum_r=0;
    for(int i = 0; i < num_l; i++){
        for(int j = 0; j < num_c; j++){
            if (mz_a.mz[i][j] != 0){
                for(int k = 0; k < num_fs; k++){
                    for(int n = 0; n < num_c; n++){
                        sum_l = sum_l + 2 * (mz_a.mz[i][n] - mz_b.mz[i][n]) * (-1 * mz_r_prev.mz[k][n]);
                    }
                    for(int n = 0; n < num_l; n++){
                        sum_r = sum_r + 2 * (mz_a.mz[n][j] - mz_b.mz[n][j]) * (-1 * mz_l_prev.mz[n][k]); 
                    }
                    l->mz[i][k] = mz_l_prev.mz[i][k] - alpha * sum_l;
                    sum_l=0;
                    r->mz[k][j] = mz_r_prev.mz[k][j] - alpha * sum_r;
                    sum_r=0;
                }
            }
        }
    }
}