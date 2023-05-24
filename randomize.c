#include <stdio.h>
#include <stdlib.h>
#include <time.h>

const int size = 100;

int main()
{
    srand(time(0));
    printf("%d", size);
    //dynamically allocating 2d matrix
    int **matrixA = (int **)malloc(size * sizeof(int *));
    int **matrixB = (int **)malloc(size * sizeof(int *));
    int **matrixC = (int **)malloc(size * sizeof(int *));
    for (int i = 0; i < size; i++) 
    {
        matrixA[i] = (int *)malloc(size * sizeof(int));
        matrixB[i] = (int *)malloc(size * sizeof(int));
        matrixC[i] = (int *)malloc(size * sizeof(int));
    }
    
    //randomly generating matrix A and B
    for (int i=0; i < size; i++)
    {
        for (int j=0; j< size ; j++)
        {
            matrixA[i][j]=rand() % 100;
            matrixB[i][j]=rand() % 100;
        }
    }

    //multiplying the two matrices
    for (int i = 0; i < size; i++) 
    {
      for (int j = 0; j < size; j++) 
      {
         for (int k = 0; k < size; k++) 
         {
            matrixC[i][j] += matrixA[i][k] * matrixB[k][j];
         }
      }
   }





    //Now writing the randomly generated matrix to file
    FILE* fp = fopen("matrixA.txt", "w");
    //first printing the size on the first line
    fprintf(fp, "%d\n", size);
    //Now writing the matrix
    for (int i = 0; i < size; i++) 
    {
        for (int j = 0; j < size; j++) 
        {
            fprintf(fp, "%d ", matrixA[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    fp = fopen("matrixB.txt", "w");
    //first printing the size on the first line
    fprintf(fp, "%d\n", size);
    //Now writing the matrix
    for (int i = 0; i < size; i++) 
    {
        for (int j = 0; j < size; j++) 
        {
            fprintf(fp, "%d ", matrixB[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    fp = fopen("matrixC.txt", "w");
    //Now writing the matrix
    for (int i = 0; i < size; i++) 
    {
        for (int j = 0; j < size; j++) 
        {
            fprintf(fp, "%d ", matrixC[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);



    //freeing the memory
    for (int i = 0; i < size; i++) {
        free(matrixA[i]);
        free(matrixB[i]);
        free(matrixC[i]);
    }
    free(matrixA);
    free(matrixB);
    return 0;
}
