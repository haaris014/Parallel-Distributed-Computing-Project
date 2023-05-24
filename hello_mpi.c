#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


struct keyVal {   
  int iIndex;
  int kIndex;         
  int value;
}; 


  
struct sortedKeyVal {
    int iIndex;
    int kIndex;
    int * arr;
};


struct keyVal* mapper(int matrixSize, int localSize, int * matrix, int *, int);
struct sortedKeyVal* shuffler(struct keyVal* keyValArray, int keyValSize, int );
int* reducer(int * iArray, int * kArray, int *valArray, int size, int);
bool compareMatrix(int * valArray, int * serialArray, int);



int main(int argc, char** argv) 
{
    // Initialize MPI environment
    MPI_Init(&argc, &argv);
    

    // Get rank and size of the current process
    int rank, numProc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    int matrixSize=0;
    int * matrix1, *matrix2, *localMatrix1, *localMatrix2;
    int localMatrixSize;
    int *sendcounts;
    int *displs;
    struct keyVal *keyValArray;
    struct sortedKeyVal *sortedKeyValArray;
    int structSize; //size of structure array
    int *iArray, *kArray, *valArray;
    char processorName[MPI_MAX_PROCESSOR_NAME];
    int nameLen;
    MPI_Get_processor_name(processorName, &nameLen);
    //printf("Rank %d is running on %s\n", rank, processorName);

    //checking if filenames were provided as parameter
    
    if (argc < 3) 
    {
        if (rank == 0) 
        {
            printf("Error! No input file provided\n");
        }
        MPI_Finalize();
        
        return 0;
    }

    if(rank == 0)
    {
        printf("Master with process_id %d running on %s\n", rank, processorName);
        const char* filename1 = argv[1];
        const char* filename2 = argv[2];
        //reading the matrix A from the file using the name provided
        FILE* fp = fopen(filename1, "r");
        if (!fp) 
        {
            printf("Error! Could not open file '%s'\n", filename1);
            MPI_Finalize();
        }
        //reading the number of rows/columns from the first line
        fscanf(fp, "%d", &matrixSize);
        matrix1 = (int *)malloc(matrixSize * matrixSize * sizeof(int ));
        
        for (int i = 0 ; i < matrixSize * matrixSize; i++)
        {
           fscanf(fp, "%d", &matrix1[i]); 
        }
        fclose(fp);

        //reading the matrix B from the file using the name provided
        fp = fopen(filename2, "r");
        if (!fp) 
        {
            printf("Error! Could not open file '%s'\n", filename2);
            MPI_Finalize();
        }
        //reading the number of rows/columns from the first line
        fscanf(fp, "%d", &matrixSize);
        matrix2 = (int *)malloc(matrixSize * matrixSize * sizeof(int ));
        
        for (int i = 0 ; i < matrixSize * matrixSize; i++)
        {
           fscanf(fp, "%d", &matrix2[i]); 
        }
        fclose(fp);
        localMatrix2=matrix2;

    }

    //now broadcasting the matrix size to all the threads
    MPI_Bcast(&matrixSize, 1, MPI_INT, 0, MPI_COMM_WORLD);


    //now the matrix would be distributed to the mappers
    sendcounts = (int *)malloc(numProc * sizeof(int));
    displs = (int *)malloc(numProc * sizeof(int));
    sendcounts[0]=0;    //so that master doesn't receive any portion
    displs[0]=0;

    for (int i = 1; i < numProc; i++) 
    {
        sendcounts[i] = (matrixSize / (numProc - 1)) * matrixSize;
        displs[i] = (i - 1) * sendcounts[i];
    }
    //to catter remaining values
    int remaining=matrixSize * matrixSize  - (numProc - 1 ) * sendcounts[numProc - 1];
    sendcounts[numProc - 1] = sendcounts[numProc - 1 ] + remaining;

    
    localMatrixSize = sendcounts[rank];
    localMatrix1 = (int *)malloc(localMatrixSize * sizeof(int));
    if (rank != 0)
        localMatrix2 = (int *)malloc(matrixSize * matrixSize * sizeof(int ));
    
   
    
    //scattering first matrix to the processes
    if(rank != 0)
    {
        printf("Task Map assigned to process %d\n", rank);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Scatterv(&matrix1[0], sendcounts, displs, MPI_INT, localMatrix1, localMatrixSize, MPI_INT, 0, MPI_COMM_WORLD);

    //brodcasting the second matrix as the whole matrix has to be sent
    MPI_Bcast(localMatrix2, matrixSize * matrixSize , MPI_INT, 0, MPI_COMM_WORLD);
    if(rank != 0)
    {
        printf("Process %d received task Map on %s\n", rank, processorName);
    }
     MPI_Barrier(MPI_COMM_WORLD);

    if (rank != 0)
    {
        //key-value array will be recieved by each non-root process
        keyValArray= mapper(matrixSize, localMatrixSize, localMatrix1, localMatrix2, rank);
        printf("Process %d has completed task Map\n", rank);
        structSize=(long long int)localMatrixSize/matrixSize * matrixSize * matrixSize;
    }
    //waiting for all mappers to finish the work
    MPI_Barrier(MPI_COMM_WORLD);

    free(localMatrix1);
    free(localMatrix2);

    //issue unidentified (why can't i free matrix 1 and 2)
    // if(rank == 0)
    // {
    // free(matrix1);
    // free(matrix2);
    // }
    
    if (rank == 0)
    {
        structSize = (long long int)matrixSize * matrixSize * matrixSize;
        keyValArray = (struct keyVal*)malloc(structSize * sizeof(struct keyVal));
    }
   
    
    //modiftying the sendCount and displ according to the expected number of key-value pairs
    for (int i = 1; i < numProc; i++) 
    {
        
        sendcounts[i] = sendcounts[i]/matrixSize * matrixSize * matrixSize;
        displs[i] = (i - 1) * sendcounts[i - 1];
    }
    
    MPI_Datatype mystruct_type;
    MPI_Type_contiguous(3, MPI_INT, &mystruct_type);
    MPI_Type_commit(&mystruct_type);

    //recieving the key-value pairs from all the mappers
    MPI_Gatherv(keyValArray, structSize, mystruct_type, keyValArray, sendcounts, displs, mystruct_type, 0, MPI_COMM_WORLD);
    
    
    free(sendcounts);
    free(displs);
    
  
    //now node manager will perform shuffling
    if(rank == 0)
    {
        //sorted key Val Array has the sorted keys
        sortedKeyValArray=shuffler(keyValArray, structSize, matrixSize);
        //dividing the recieved structure array into individual int arrays
        iArray=(int *)malloc (matrixSize * matrixSize *sizeof(int));
        kArray=(int *)malloc (matrixSize * matrixSize *sizeof(int));
        //values for each key will be combined in one large array
        valArray=(int *)malloc (matrixSize * matrixSize * matrixSize *sizeof(int));
        int index=0;
        for (int i=0; i< matrixSize * matrixSize; i++)
        {
            iArray[i] = sortedKeyValArray[i].iIndex;
            kArray[i] = sortedKeyValArray[i].kIndex;

            for (int j=0; j< matrixSize; j++)
            {
              valArray[index++]=  sortedKeyValArray[i].arr[j];
            }
        }
        //now freeing sorted key val array as no longer needed
        free(sortedKeyValArray);
    }

    //now freeing key val array as no longer needed
    free(keyValArray);

    

    //now the matrix would be distributed to the mappers
    //numProc has the number of reducers
    int acutalNumProc=numProc;
    numProc /=2;

    sendcounts = (int *)malloc(acutalNumProc * sizeof(int));
    displs = (int *)malloc(acutalNumProc * sizeof(int));
    sendcounts[0]=0;    //so that master doesn't receive any portion
    displs[0]=0;

    for (int i = 1; i < numProc; i++) 
    {
        sendcounts[i] = (int)(matrixSize * matrixSize)  / (numProc - 1);//10
        displs[i] = (i - 1) * sendcounts[i];
    }
    remaining=matrixSize * matrixSize - (numProc - 1 ) * sendcounts[numProc - 1];
    sendcounts[numProc - 1] = sendcounts[numProc - 1 ] + remaining;

    for (int i = numProc; i < acutalNumProc; i++) 
    {
        sendcounts[i] = 0;
        displs[i] = 0;
    }

    localMatrixSize = sendcounts[rank];
    if(rank != 0 && rank <= numProc)
    {
        //this will be used to get values from node manager
        iArray=(int*)malloc(localMatrixSize * sizeof(int));
        kArray=(int*)malloc(localMatrixSize * sizeof(int));
        valArray=(int*)malloc(localMatrixSize * matrixSize * sizeof(int));
    }
    
     if(rank != 0 && rank <= numProc)
    {
        printf("Task Reduce assigned to process %d\n", rank);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //scattering the key array
    MPI_Scatterv(iArray, sendcounts, displs, MPI_INT, iArray, localMatrixSize, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(kArray, sendcounts, displs, MPI_INT, kArray, localMatrixSize, MPI_INT, 0, MPI_COMM_WORLD);
    
    
    //modifying sendcount and displ according to requirement
    for (int i = 1; i < numProc; i++) 
    {
        sendcounts[i] *= matrixSize;
        displs[i] = displs[i - 1] + sendcounts[i - 1];
    }

    //scattering the value array
    MPI_Scatterv(valArray, sendcounts, displs, MPI_INT, valArray, sendcounts[rank], MPI_INT, 0, MPI_COMM_WORLD);
    if(rank != 0 && rank <= numProc)
    {
        printf("Process %d received task Reduce on %s\n", rank, processorName);
    }
     MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0)
        free(valArray);


    //now sending the individual arrays to the reducer function:
    if(rank != 0 && rank <= numProc)
    {
        valArray=reducer(iArray, kArray, valArray, localMatrixSize, matrixSize);
         printf("Process %d has completed task Reduce\n", rank);
    }

    //waiting for all reducers to finish their work
    MPI_Barrier(MPI_COMM_WORLD);

    //now gathering result back to node manager

    //modifying sendcount and displ according to the valArray
    for (int i = 1; i < numProc; i++) 
    {
        sendcounts[i] /= matrixSize;
        displs[i] = displs[i - 1] + sendcounts[i - 1];
    }

    //initlaizing array for the node manager
    if(rank == 0)
    {
        valArray=(int*)malloc(matrixSize * matrixSize * sizeof(int));
    }

    MPI_Gatherv(valArray, localMatrixSize, MPI_INT, valArray, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);


    
    if(rank == 0)
    {
        //writing values to the file
        FILE *fp = fopen("Output.txt", "w");
        for (int i = 0; i < matrixSize; i++) 
        {
            for (int j = 0; j < matrixSize; j++) 
            {
                fprintf(fp, "%d ", valArray[i*matrixSize + j]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);


        //master will verify that the serial and parallel version's outputs are same
        //reading values from the file which has output of the serial version
        int *serialArray=(int*)malloc(matrixSize * matrixSize * sizeof(int));
        fp = fopen("matrixC.txt", "r");
        if (!fp) 
        {
            printf("Error! Could not open file\n");
            MPI_Finalize();
        }
        
        for (int i = 0 ; i < matrixSize * matrixSize; i++)
        {
           fscanf(fp, "%d", &serialArray[i]); 
        }
        fclose(fp);

        //now comapring the values
        bool equal=compareMatrix(valArray, serialArray, matrixSize);
        printf("Matrix Comparison function returned: ");

        if(equal)
        {
            printf("True\n");
        }
        else
        {
            printf("False\n");
        }

    }


    // // Finalize MPI environment
    MPI_Finalize();

    return 0;
}


bool compareMatrix(int * valArray, int * serialArray, int matrixSize)
{
    bool equal=1;
    for(int i = 0; i< matrixSize * matrixSize; i++)
    {
        if(valArray[i] != serialArray[i])
        {
            equal = 0;
            printf("valI: %d\tvalS: %d\n", valArray[i], serialArray[i]);
        }
    }
    return equal;
}

int* reducer(int * iArray, int * kArray, int *valArray, int size, int matrixSize)
{
    int index=0;
    int *sumArray=(int*)malloc(size * sizeof(int));
    //intializing each index to 0
    for (int i = 0; i < size; i++)
    {
        sumArray[i]= 0;
    }

    for (int i=0; i< size; i++)
    {
        for (int j=0; j< matrixSize; j++)
        {
            sumArray[i] += valArray[index];
            index++;
        }
    }
   
    return sumArray;
}




struct sortedKeyVal* shuffler(struct keyVal* keyValArray, int keyValSize, int matrixSize)
{
    //declaring struct array according to the no. of postitions on output matrix
    struct sortedKeyVal *SortedkeyValArray= ( struct sortedKeyVal *)malloc(matrixSize * matrixSize * sizeof(struct sortedKeyVal));
    for (int i = 0; i < matrixSize * matrixSize; i++)
    {
        int temp=matrixSize * sizeof(int);
        SortedkeyValArray[i].arr= malloc(temp);
    }
    int tempI = 1;
    int tempK=1;
    int index= 0;   //used for struct array
    int innerIndex=0;   //used for array in the struct array
    for (int i = 0; i < keyValSize; i++)
    {
        if(tempI != keyValArray[i].iIndex || tempK != keyValArray[i].kIndex)
        {
            innerIndex= 0;
            tempI = keyValArray[i].iIndex;
            tempK=keyValArray[i].kIndex;
            SortedkeyValArray[index].iIndex=keyValArray[i].iIndex;
            SortedkeyValArray[index].kIndex=keyValArray[i].kIndex;
            index++;
        }
        SortedkeyValArray[index - 1].arr[innerIndex]=keyValArray[i].value;
        innerIndex++;
    }

    return SortedkeyValArray;
}



struct keyVal* mapper(int matrixSize, int localMatrixSize, int * localMatrix1, int *localMatrix2, int rank)
{
    int localMatrixRowSize=localMatrixSize/matrixSize;
    //calculating assigned row number
    int numProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    int startingRow=(rank- 1) * (matrixSize/ (numProc - 1)) + 1;
    int endingRow = startingRow + (localMatrixSize/ matrixSize) - 1;
    struct keyVal *keyValArray= ( struct keyVal *)malloc(localMatrixRowSize * matrixSize * matrixSize * sizeof(struct keyVal));
    int temp=startingRow;
    long long int index=0;
    if (rank)
    {
        temp=startingRow;
        
        for (int i = 0; i< localMatrixRowSize; i++)
        {
            for (int j = 0; j< matrixSize; j++)
            {
                for(int k=0; k < matrixSize; k++)
                {
                    keyValArray[index].iIndex=temp - 1;
                    keyValArray[index].kIndex=j;
                    keyValArray[index].value=localMatrix1[i * matrixSize +k] * localMatrix2[k * matrixSize + j];
                    index++;
                }
                
            }
            temp++;
        }
    }

    return keyValArray;
}    
   


