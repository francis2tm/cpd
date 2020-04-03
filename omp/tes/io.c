#include <stdio.h>
#include <stdlib.h>
#include "io.h"

/*****************************************************************************************************
 * openFile ()
 * Arguments: name: nome do ficheiro a ser aberto
 *            mode: modo de abertura
 * Returns: fp. ponteiro para o ficheiro aberto
 * Description: Abre o ficheiro que tem como nome "name"
 ****************************************************************************************************/

FILE* openFile(char* name){
    FILE *fp = NULL;
    fp = fopen(name, "r");

    if(fp == NULL){
        fprintf(stderr, "Erro abrir ficheiro");
        exit(-1);
    }

    return fp;
}