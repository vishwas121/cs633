#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>

int N, P, time_steps, myrank;
int row, col; // row and col of process in process matrix
int up, down, left, right; // bool values to check boundary conditions of a process in process matrix

double** init_mat(){
    double** mat = malloc(sizeof(double*)*N);
    for (int i=0; i<N; i++){
        mat[i] = malloc(sizeof(double)*N);
    }
    return mat;
}

void random_init(double** mat){
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            mat[i][j] = drand48();
        }
    }
}

void free_mat(double** mat){
    for (int i=0; i<N; i++){
        free(mat[i]);
    }
    free(mat);
}

void compute_edge_points(double** calc_mat, double** mat, double* recvup, double* recvdown, double* recvleft, double* recvright){
    for (int i=1; i<N-1; i++){
        if (up) calc_mat[0][i] = (mat[0][i+1]+mat[0][i-1]+mat[1][i]+recvup[i])/4;
        else calc_mat[0][i] = (mat[0][i+1]+mat[0][i-1]+mat[1][i])/3;
        if (down) calc_mat[N-1][i] = (mat[N-1][i+1]+mat[N-1][i-1]+mat[N-2][i]+recvdown[i])/4;
        else calc_mat[N-1][i] = (mat[N-1][i+1]+mat[N-1][i-1]+mat[N-2][i])/3;
        if (left) calc_mat[i][0] = (mat[i-1][0]+mat[i+1][0]+mat[i][1]+recvleft[i])/4;
        else calc_mat[i][0] = (mat[i-1][0]+mat[i+1][0]+mat[i][1])/3;
        if (right) calc_mat[i][N-1] = (mat[i-1][N-1]+mat[i+1][N-1]+mat[i][N-2]+recvright[i])/4;
        else calc_mat[i][N-1] = (mat[i-1][N-1]+mat[i+1][N-1]+mat[i][N-2])/3;
    }
}

void compute_inner_points(double** calc_mat, double** mat){
    for (int i=1; i<N-1; i++){
        for (int j=1; j<N-1; j++){
            calc_mat[i][j] = (mat[i-1][j]+mat[i+1][j]+mat[i][j-1]+mat[i][j+1])/4;
        }
    }
}

void compute_corner_points(double** calc_mat, double** mat, double* recvup, double* recvdown, double* recvleft, double* recvright){
    if (up){
        if (left) calc_mat[0][0] = (mat[1][0]+mat[0][1]+recvup[0]+recvleft[0])/4;
        else calc_mat[0][0] = (mat[1][0]+mat[0][1]+recvup[0])/3;
        if (right) calc_mat[0][N-1] = (mat[1][N-1]+mat[0][N-2]+recvup[N-1]+recvright[0])/4;
        else calc_mat[0][N-1] = (mat[1][N-1]+mat[0][N-2]+recvup[N-1])/3;
    }
    else {
        if (left) calc_mat[0][0] = (mat[1][0]+mat[0][1]+recvleft[0])/3;
        else calc_mat[0][0] = (mat[1][0]+mat[0][1])/2;
        if (right) calc_mat[0][N-1] = (mat[1][N-1]+mat[0][N-2]+recvright[0])/3;
        else calc_mat[0][N-1] = (mat[1][N-1]+mat[0][N-2])/2;
    }
    if (down){
        if (left) calc_mat[N-1][0] = (mat[N-2][0]+mat[N-1][1]+recvdown[0]+recvleft[N-1])/4;
        else calc_mat[N-1][0] = (mat[N-2][0]+mat[N-1][1]+recvdown[0])/3;
        if (right) calc_mat[N-1][N-1] = (mat[N-2][N-1]+mat[N-1][N-2]+recvdown[N-1]+recvright[N-1])/4;
        else calc_mat[N-1][N-1] = (mat[N-2][N-1]+mat[N-1][N-2]+recvdown[N-1])/3;
    }
    else {
        if (left) calc_mat[N-1][0] = (mat[N-2][0]+mat[N-1][1]+recvleft[N-1])/3;
        else calc_mat[N-1][0] = (mat[N-2][0]+mat[N-1][1])/2;
        if (right) calc_mat[N-1][N-1] = (mat[N-2][N-1]+mat[N-1][N-2]+recvright[N-1])/3;
        else calc_mat[N-1][N-1] = (mat[N-2][N-1]+mat[N-1][N-2])/2;
    }
}

double method1(double** mat, double** calc_mat){
    random_init(calc_mat);
    double recvup[N], recvdown[N], recvleft[N], recvright[N];
    double s_time = MPI_Wtime();
    for (int t=0; t<time_steps; t++){
        // UPDATE ORIGINAL MATRIX
        for (int i=0; i<N; i++){
            for (int j=0; j<N; j++){
                mat[i][j] = calc_mat[i][j];
            }
        }
        // COMMUNICATION STEP
        MPI_Status status;
        // send up and recv from up
        if (up){
            for (int i=0; i<N; i++){
                MPI_Send(&mat[0][i], 1, MPI_DOUBLE, P*(row-1)+col, t, MPI_COMM_WORLD);
                MPI_Recv(&recvup[i], 1, MPI_DOUBLE, P*(row-1)+col, t, MPI_COMM_WORLD, &status);
            }
        }
        // recv from down and send down
        if (down){
            for (int i=0; i<N; i++){
                MPI_Recv(&recvdown[i], 1, MPI_DOUBLE, P*(row+1)+col, t, MPI_COMM_WORLD, &status);
                MPI_Send(&mat[N-1][i], 1, MPI_DOUBLE, P*(row+1)+col, t, MPI_COMM_WORLD);
            }
        }
        // send left and recv from left
        if (left){
            for (int i=0; i<N; i++){
                MPI_Send(&mat[i][0], 1, MPI_DOUBLE, P*row+col-1, t, MPI_COMM_WORLD);
                MPI_Recv(&recvleft[i], 1, MPI_DOUBLE, P*row+col-1, t, MPI_COMM_WORLD, &status);
            }
        }
        // recv from right and send right
        if (right){
            for (int i=0; i<N; i++){
                MPI_Recv(&recvright[i], 1, MPI_DOUBLE, P*row+col+1, t, MPI_COMM_WORLD, &status);
                MPI_Send(&mat[i][N-1], 1, MPI_DOUBLE, P*row+col+1, t, MPI_COMM_WORLD);
            }
        }
        // COMPUTE EDGE POINTS
        compute_edge_points(calc_mat, mat, recvup, recvdown, recvleft, recvright);
        // COMPUTE INNER POINTS
        compute_inner_points(calc_mat, mat);
        // COMPUTE CORNER POINTS
        compute_corner_points(calc_mat, mat, recvup, recvdown, recvleft, recvright);
    }
    return MPI_Wtime()-s_time;
}

double method2(double** mat, double** calc_mat){
    random_init(calc_mat);
    double recvup[N], recvdown[N], recvleft[N], recvright[N];
    double buffup[N], buffdown[N], buffleft[N], buffright[N];
    double sendup[N], senddown[N], sendleft[N], sendright[N];
    int pos = 0;
    double s_time = MPI_Wtime();
    for (int t=0; t<time_steps; t++){
        // UPDATE ORIGINAL MATRIX
        for (int i=0; i<N; i++){
            for (int j=0; j<N; j++){
                mat[i][j] = calc_mat[i][j];
            }
        }
        // COMMUNICATION STEP
        MPI_Status status;
        // send up and recv from up
        if (up){
            for (int i=0; i<N; i++){
                MPI_Pack(&mat[0][i], 1, MPI_DOUBLE, sendup, sizeof(double)*N, &pos, MPI_COMM_WORLD);
            }
            MPI_Send(sendup, pos, MPI_PACKED, P*(row-1)+col, t, MPI_COMM_WORLD);
            pos = 0;
            MPI_Recv(buffup, sizeof(double)*N, MPI_PACKED, P*(row-1)+col, t, MPI_COMM_WORLD, &status);
            for (int i=0; i<N; i++){
                MPI_Unpack(buffup, sizeof(double)*N, &pos, &recvup[i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
            }
            pos = 0;
        }
        // recv from down and send down
        if (down){
            MPI_Recv(buffdown, sizeof(double)*N, MPI_PACKED, P*(row+1)+col, t, MPI_COMM_WORLD, &status);
            for (int i=0; i<N; i++){
                MPI_Unpack(buffdown, sizeof(double)*N, &pos, &recvdown[i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
            }
            pos = 0;
            for (int i=0; i<N; i++){
                MPI_Pack(&mat[N-1][i], 1, MPI_DOUBLE, senddown, sizeof(double)*N, &pos, MPI_COMM_WORLD);
            }
            MPI_Send(senddown, pos, MPI_PACKED, P*(row+1)+col, t, MPI_COMM_WORLD);
            pos = 0;
        }
        // send left and recv from left
        if (left){
            for (int i=0; i<N; i++){
                MPI_Pack(&mat[i][0], 1, MPI_DOUBLE, sendleft, sizeof(double)*N, &pos, MPI_COMM_WORLD);
            }
            MPI_Send(sendleft, pos, MPI_PACKED, P*row+col-1, t, MPI_COMM_WORLD);
            pos = 0;
            MPI_Recv(buffleft, sizeof(double)*N, MPI_PACKED, P*row+col-1, t, MPI_COMM_WORLD, &status);
            for (int i=0; i<N; i++){
                MPI_Unpack(buffleft, sizeof(double)*N, &pos, &recvleft[i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
            }
            pos = 0;
        }
        // recv from right and send right
        if (right){
            MPI_Recv(buffright, sizeof(double)*N, MPI_PACKED, P*row+col+1, t, MPI_COMM_WORLD, &status);
            for (int i=0; i<N; i++){
                MPI_Unpack(buffright, sizeof(double)*N, &pos, &recvright[i], 1, MPI_DOUBLE, MPI_COMM_WORLD);
            }
            pos = 0;
            for (int i=0; i<N; i++){
                MPI_Pack(&mat[i][N-1], 1, MPI_DOUBLE, sendright, sizeof(double)*N, &pos, MPI_COMM_WORLD);
            }
            MPI_Send(sendright, pos, MPI_PACKED, P*row+col+1, t, MPI_COMM_WORLD);
            pos = 0;
        }
        // COMPUTE EDGE POINTS
        compute_edge_points(calc_mat, mat, recvup, recvdown, recvleft, recvright);
        // COMPUTE INNER POINTS
        compute_inner_points(calc_mat, mat);
        // COMPUTE CORNER POINTS
        compute_corner_points(calc_mat, mat, recvup, recvdown, recvleft, recvright);
    }
    return MPI_Wtime()-s_time;
}

double method3(double** mat, double** calc_mat){
    random_init(calc_mat);
    double recvup[N], recvdown[N], recvleft[N], recvright[N];
    MPI_Datatype t_row, t_col_left, t_col_right;
    MPI_Type_vector(N, 1, 1, MPI_DOUBLE, &t_row); MPI_Type_commit(&t_row);
    int blocklens[N];
    for (int i=0; i<N; i++) blocklens[i] = 1;
    int offsets_left[N];
    for (int i=0; i<N; i++) offsets_left[i] = mat[i] - mat[0];
    MPI_Type_indexed(N, blocklens, offsets_left, MPI_DOUBLE, &t_col_left); MPI_Type_commit(&t_col_left);
    int offsets_right[N];
    for (int i=0; i<N; i++) offsets_right[i] = mat[i] - mat[0] + N - 1;
    MPI_Type_indexed(N, blocklens, offsets_right, MPI_DOUBLE, &t_col_right); MPI_Type_commit(&t_col_right);
    int pos = 0;
    double s_time = MPI_Wtime();
    for (int t=0; t<time_steps; t++){
        // UPDATE ORIGINAL MATRIX
        for (int i=0; i<N; i++){
            for (int j=0; j<N; j++){
                mat[i][j] = calc_mat[i][j];
            }
        }
        // COMMUNICATION STEP
        MPI_Status status;
        // send up and recv from up
        if (up){
            MPI_Send(&mat[0][0], 1, t_row, P*(row-1)+col, t, MPI_COMM_WORLD);
            MPI_Recv(recvup, N, MPI_DOUBLE, P*(row-1)+col, t, MPI_COMM_WORLD, &status);
        }
        // recv from down and send down
        if (down){
            MPI_Recv(recvdown, N, MPI_DOUBLE, P*(row+1)+col, t, MPI_COMM_WORLD, &status);
            MPI_Send(&mat[N-1][0], 1, t_row, P*(row+1)+col, t, MPI_COMM_WORLD);
        }
        // send left and recv from left
        if (left){
            MPI_Send(&mat[0][0], 1, t_col_left, P*row+col-1, t, MPI_COMM_WORLD);
            MPI_Recv(recvleft, N, MPI_DOUBLE, P*row+col-1, t, MPI_COMM_WORLD, &status);
        }
        // recv from right and send right
        if (right){
            MPI_Recv(recvright, N, MPI_DOUBLE, P*row+col+1, t, MPI_COMM_WORLD, &status);
            MPI_Send(&mat[0][0], 1, t_col_right, P*row+col+1, t, MPI_COMM_WORLD);
        }
        // COMPUTE EDGE POINTS
        compute_edge_points(calc_mat, mat, recvup, recvdown, recvleft, recvright);
        // COMPUTE INNER POINTS
        compute_inner_points(calc_mat, mat);
        // COMPUTE CORNER POINTS
        compute_corner_points(calc_mat, mat, recvup, recvdown, recvleft, recvright);
    }
    MPI_Type_free(&t_row); MPI_Type_free(&t_col_left); MPI_Type_free(&t_col_right);
    return MPI_Wtime()-s_time;
}

int main(int argc, char* argv[]){
    int N_2 = atoi(argv[1]);
    time_steps = atoi(argv[2]);

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int P_2;
    MPI_Comm_size(MPI_COMM_WORLD, &P_2);
    P = (int)(sqrt((double)(P_2)));
    N = (int)(sqrt((double)(N_2)));

    // process matrix coordinates
    row = myrank/P;
    col = myrank%P;
    // possible sendrecv directions of data matrix
    up = row>0 ? 1:0;
    down = row<P-1 ? 1:0;
    left = col>0 ? 1:0;
    right = col<P-1 ? 1:0;

    double** mat = init_mat();
    double** calc_mat = init_mat();

    double times[3];
    times[0] = method1(mat, calc_mat);
    times[1] = method2(mat, calc_mat);
    times[2] = method3(mat, calc_mat);

    double max_times[3];
    MPI_Reduce(times, max_times, 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (!myrank){
        printf("%f\n%f\n%f\n", max_times[0], max_times[1], max_times[2]);
    }
    
    free_mat(mat); free_mat(calc_mat);
    MPI_Finalize();
    return 0;
}
