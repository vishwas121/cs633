#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>

#define FLOAT_MAX 3.40282e+038

FILE* file;
FILE* out;

// Data format :: Blocks of 10000(block_size) row capacity, rows variable == rows in last block, each row of 43 floats.

int checkEOF(){    		// return 1 if EOF
	int c;
	c = fgetc(file);
	if (c == EOF) return 1;
   	ungetc(c, file);
   	return 0;
}

int get_column_size(){
	int c;
	int column_size = 1;
	c = fgetc(file);
	while (c != '\n'){
		if (c==',') column_size++;
		c = fgetc(file);
	}
	return column_size;
}
int main(int argc, char* argv[]){

	MPI_Init(NULL,NULL);

	file = fopen(argv[1], "r");
	out = fopen("output.txt", "w");
	int num_of_blocks = 100000, block_size = 500, column_size = 0;    // change this accordingly
	int rows = 0, block_count = 0;

	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	float** data = (float**)malloc(sizeof(float*)*num_of_blocks);
	int rows_left = 1;          						// flag to check if any input is left 
	int i = 0, j = 0, k = 0;

	if (rank == 0){
		column_size = get_column_size();
		while (rows_left){
			block_count++;
			rows = 0;
			data[i] = (float*)malloc(sizeof(float)*block_size*column_size);                          // May use data format 2
			// for (j=0; j<block_size; j++) data[i][j] = (float*)malloc(sizeof(float)*column_size);

			int is_end = 0;
			j = 0;
			while (!is_end && j<block_size){
				rows++;
				for (k=0; k<(column_size-1); k++){
					fscanf (file, "%f,", &data[i][j*column_size+k]);
				}
				fscanf (file, "%f\n", &data[i][j*column_size + k]);
				if (checkEOF()){
					is_end = 1;
					rows_left = 0;
				}
				else j++;
			}
			if (j<block_size) block_count--;
			i++;
		}
	}
	int initials[4];
		if (!rank){
			initials[0] = block_size;
			initials[1] = column_size;
			initials[2] = block_count;
			initials[3] = rows;
		}
		double x =  MPI_Wtime();
		MPI_Bcast(initials, 4, MPI_INT, 0, MPI_COMM_WORLD);
		x = MPI_Wtime()-x;
		printf("%lf\n",x);
	double s_time = MPI_Wtime();


	// For one process.
	if (size == 1){				// sequential
		float min_val[column_size-2];
		for (i=0; i<column_size-2; i++) min_val[i] = FLOAT_MAX;
		float global_min = FLOAT_MAX;

		for (i=0; i<(block_count-1); i++){
			for (j=0; j<block_size; j++){
				for (k=0; k<(column_size-2); k++){
					if (min_val[k] > data[i][j*column_size + k+2]) min_val[k] = data[i][j*column_size +k+2];
				}
			}	
		}
		for (j=0; j<rows; j++){
			for (k=0; k<(column_size-2); k++){
				if (min_val[k] > data[i][j*column_size + k+2]) min_val[k] = data[i][j*column_size +k+2];
			}
		}
		for (k=0; k<(column_size-2); k++){
			if(global_min > min_val[k]) global_min = min_val[k];
		}

		double e_time = MPI_Wtime() - s_time;

		fprintf(out, "%f", min_val[0]);
		for (k=1; k<(column_size-2); k++){
			fprintf (out, ",%f", min_val[k]);
		}
		fprintf (out, "\n%f\n%f\n", global_min, e_time);
		printf ("%f\n", e_time);
	}

	/* 	Data distribution
		Assumption block_size is optimal i.e. for data distribution
				try different blocksizes to test with differnet column sizes
				same amount of data with single large send or multiple small sends.
		Strategy :: divide block_count/processes each
					for remaining blocks + rows ==> (blocks+rows)/processes 
		block_count_per_process are used. Last block is partially filled.

	*/
	if (size > 1){
		if (rank){
			block_size = initials[0];
			column_size = initials[1];
			block_count = initials[2];
			rows = initials[3];
		}

		int tot_rows = block_count*block_size + rows;
		int block_count_per_process = block_count/size;
		int remain_data = (block_count%size)*block_size + rows;

		//Data format 1::
		// float*** data_per_process = (float***)malloc(sizeof(float**)*(block_count_per_process));
		// for(i=0;i<block_count_per_process;i++){
		// 	data_per_process[i] = (float**)malloc(sizeof(float*)*block_size);
		// 	for(j=0;j<block_size;j++){
		// 		data_per_process[i][j] = (float*)malloc(sizeof(float)*column_size);
		// 	}
		// }

		//Data format 2::
		float** data_per_process = (float**)malloc(sizeof(float*)*(block_count_per_process));
		for (i=0; i<block_count_per_process; i++){
			data_per_process[i] = (float*)malloc(sizeof(float)*block_size*column_size);
		}

		//Data format 3::
		// float** data_per_process = (float**)malloc(sizeof(float*)*(block_count_per_process*block_size));
		// for(i=0;i<(block_count_per_process*block_size);i++){
		// 	data_per_process[i] = (float*)malloc(sizeof(float)*column_size);
		// }

		// Data distribution code::
		// Method 1 :: Scatter blocks 1,2,3.. to ranks 1,2,3.. and then 4,5,6.. to ranks 1,2,3..
		// Seems better than sending blocks 1,2,3 to rank 1 and blocks 4,5,6 to rank 2 and so on. 


		// Assuming data format 2 ::  
		/*  For vectored scatter.
		int sendcounts[size],sdispls[size], recvcounts[size], rdispls[size];
    	for (i=0; i<size; i++){
    		sendcounts[i] = block_size*column_size;
        	sdispls[i] = i*block_size*column_size;
        	recvcounts[i] = block_size*column_size;
        	rdispls[i] = i*block_size*column_size;
    	}
    	*/
		// int sendcounts,recvcounts;
		// sendcounts = block_size*column_size;
		// recvcounts = sendcounts;

		// for(i=0;i<(block_count_per_process-1);i++){
		// 	for(j=0;j<size;j++){
		// 		MPI_Iscatterv(data[i*size],sendcounts[j],sdispls[j],MPI_FLOAT,data_per_process[i],recvcounts[j],MPI_FLOAT,0,MPI_COMM_WORLD);
		// 		MPI_Scatter(data[i*size],sendcounts,MPI_FLOAT,data_per_process[i],recvcounts,MPI_FLOAT,0,MPI_COMM_WORLD);
		// 	}
		// 	// Computation of minimum on received data.
		// }

		//Method 2 :: Individual sendrecvs for each block and processinng in between.
		// Seems better than scatterv as computation on some nodes starts concurrently while root sends to others.
		// If Computation time of blocksize > communication time of blocksize ==> send seems better. otherwise scatter. Have to test with data.
		// Root will bottleneck on sequential sends as it will have data to compute after sends. Have to distribute more data to others.
		// make root compute==0? Not good for lesser number of processes. Have to compute some data on root also.
		// sending equal data.
		double x_time = MPI_Wtime();
		if (!rank){
			MPI_Request req[block_count_per_process][size-1];
			float min_val[column_size-2];
			for (i=0; i<column_size-2; i++) min_val[i] = FLOAT_MAX;
			for (i=0; i<block_count_per_process; i++){
				for(j=1; j<size; j++){
					MPI_Isend (data[i*size+j], block_size*column_size, MPI_FLOAT, j, i*size+j, MPI_COMM_WORLD, &req[i][j-1]);
				}
				for (j=0; j<block_size; j++){
					for (k=0; k<column_size-2; k++){
						if (min_val[k] > data[i*size][j*column_size + k+2]) min_val[k] = data[i*size][j*column_size +k+2];
					}
				}
				MPI_Waitall(size-1, req[i], MPI_STATUS_IGNORE);
				// Computation on own data. Good performance if computation of the first rank is root finish at same time.
				// Have to make computation of root + communication to other processes == computation on rank 1. Have to experiment with different block_sizes.
			}
			float final[column_size-2];
			MPI_Reduce(min_val, final, column_size-2, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
			float global_min = FLOAT_MAX;
			for (k=0; k<column_size-2; k++){
				if(global_min > final[k]) global_min = final[k];
			}

			double e_time = MPI_Wtime() - s_time;
			double max_time;
			MPI_Reduce(&e_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

			fprintf(out, "%f", final[0]);
			for (k=1; k<column_size-2; k++){
				fprintf (out, ",%f", final[k]);
			}
			fprintf (out, "\n%f\n%f\n", global_min, max_time);
			printf ("rank %d : %f\n", rank, e_time);
		}
		if (rank){
			MPI_Recv(data_per_process[0], block_size*column_size, MPI_FLOAT, 0, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Request req;
			float min_val[block_count_per_process][column_size-2];
			for(i=1; i<block_count_per_process; i++){
				MPI_Irecv(data_per_process[i], block_size*column_size, MPI_FLOAT, 0, i*size+rank, MPI_COMM_WORLD, &req);
				//Computation of minimum.
				for (j=0; j<column_size-2; j++) min_val[i-1][j] = FLOAT_MAX;
				for (j=0; j<block_size; j++){
					for (k=0; k<column_size-2; k++){
						if (min_val[i-1][k] > data_per_process[i-1][j*column_size + k + 2]) min_val[i-1][k] = data_per_process[i-1][j*column_size + k + 2];
					}	
				}
				MPI_Wait(&req, MPI_STATUS_IGNORE);
			}
			//Computation of minimum.
			for (j=0; j<column_size-2; j++) min_val[block_count_per_process-1][j] = FLOAT_MAX;
			for (j=0; j<block_size; j++){
				for (k=0; k<column_size-2; k++){
					if (min_val[block_count_per_process-1][k] > data_per_process[block_count_per_process-1][j*column_size + k + 2]) min_val[block_count_per_process-1][k] = data_per_process[block_count_per_process-1][j*column_size + k + 2];
				}	
			}
			for (j=0; j<column_size-2; j++){
				for (i=0; i<block_count_per_process; i++){
					min_val[0][j] = min_val[0][j] > min_val[i][j] ? min_val[i][j] : min_val[0][j];
				}
			}
			double final[column_size-2];
			MPI_Reduce(min_val[0], final, column_size-2, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
			double e_time = MPI_Wtime() - s_time;
			double max_time;
			MPI_Reduce(&e_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			printf ("rank %d : %f\n", rank, e_time);
		}
	}

	MPI_Finalize();
	return 0;
}

// handle remaining data
// get rid of 2 cols, dont communicate them
// segfault checkss_time = MPI_Wtime();
		