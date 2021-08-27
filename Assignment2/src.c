#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>

#define KB 1024
int D, size, rank, len, ppn, groups;

// randomly initialises a double array of length len
double* random_init(int len){
    srand48(rank);
    double* buf = (double*)malloc(len*sizeof(double));
    for (int i=0; i<len; i++){
        buf[i] = drand48();
    }
    return buf;
}

double avg(double* max_time){
    double ret = 0;
    for (int i=0; i<5; i++) ret += max_time[i];
    return ret/5;
}

void swap_data(double* buf, int off1, int off2){
    double temp;
    for (int i=0; i<len; i++){
        temp = buf[off1+i];
        buf[off1+i] = buf[off2+i];
        buf[off2+i] = temp;
    }
}

// takes a rank and finds closest rank to it topologically
int find_closest_to(int rank, int* mapped){
    // this simple logic works because with the rank ordering we already have :
    // 1. a rank cannot have an unmapped rank behind it
    // 2. the first unmapped rank ahead will be one of the closest ranks
    for (int i=rank+1; i<size; i++){
        if (!mapped[i]){
            mapped[i] = 1;
            return i;
        }
    }
    return -1;
}

// this function is taken from the cited paper, returns new rank mapping for binomial gather
void binomial_gather_rank_reorder(int* reordered_ranks){
    reordered_ranks[0] = 0;
    for (int i=1; i<size; i++) reordered_ranks[i] = -1;

    int mapped[size];
    mapped[0] = 1; 
    for (int i=1; i<size; i++) mapped[i] = 0;

    int v[size];
    int v_len = 0;
    v[v_len++] = 0;

    int i = size/2;
    while (i > 0){
        for (int j=0; j<v_len; j++){
            int ref_rank = v[j];
            if (ref_rank+i < size){
                if (reordered_ranks[ref_rank+i] != -1) continue;
                int new_rank = ref_rank + i;
                int target_core = find_closest_to(reordered_ranks[ref_rank], mapped);
                if (target_core == -1) return;
                reordered_ranks[new_rank] = target_core;
                v[v_len++] = new_rank;
            }
        }
        i = i/2;
    }
}

double standard_bcast(){
    double* buf = random_init(len); // randomly initialise buf
    double time[5];
    int itr = 5;
    while (itr--){
        MPI_Barrier(MPI_COMM_WORLD); // for fair time measurement
        time[itr] = MPI_Wtime();

        MPI_Bcast(buf, len, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        time[itr] = MPI_Wtime()-time[itr];
    }
    free(buf);
    double max_time[5];
    MPI_Reduce(time, max_time, 5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return avg(max_time);
}

double optimized_bcast(int* nodes_per_group){
    double* buf = random_init(len); // randomly initialise buf
    double time[5];
    int itr = 5;
    while (itr--){
        MPI_Barrier(MPI_COMM_WORLD); // for fair measurement
        time[itr] = MPI_Wtime();
        
        MPI_Request req1, req2, req3;
        MPI_Group inter_group_group, inter_node_group, intra_node_group; 
        MPI_Comm inter_group, inter_node, intra_node; // multilevel communicators
        int flag1 = 0; int flag2 = 0; 
        int group_id = 0;

        MPI_Group world_group;
        MPI_Comm_group(MPI_COMM_WORLD, &world_group);

        if (groups > 1){
            int inter_group_ranks[groups];
            int exp_rank = 0;
            for (int i=0; i<groups; i++){
                if (rank == exp_rank) flag1 = 1; // flag1 = 1 for master ranks at inter group level
                if (rank >= exp_rank && rank < exp_rank+nodes_per_group[i]*ppn) group_id = i; // find group id to which the process belongs
                inter_group_ranks[i] = exp_rank;
                exp_rank += nodes_per_group[i]*ppn;
            }
            MPI_Group_incl(world_group, groups, inter_group_ranks, &inter_group_group);
            MPI_Comm_create_group(MPI_COMM_WORLD, inter_group_group, 0, &inter_group); // communicator for inter group master ranks
            if (flag1) MPI_Ibcast(buf, len, MPI_DOUBLE, 0, inter_group, &req1); // 1st level bcast
        }

        if (nodes_per_group[group_id] > 1){
            flag2 = rank%ppn == 0 ? 1 : 0; // flag2 = 1 for master ranks at inter node level (within a group)
            int inter_node_ranks[nodes_per_group[group_id]];
            int offset = 0;
            for (int i=0; i<group_id; i++){
                offset += nodes_per_group[i];
            }
            for (int i=0; i<nodes_per_group[group_id]; i++){
                inter_node_ranks[i] = ppn*i + offset*ppn; 
            }
            MPI_Group_incl(world_group, nodes_per_group[group_id], inter_node_ranks, &inter_node_group);
            MPI_Comm_create_group(MPI_COMM_WORLD, inter_node_group, 0, &inter_node); // communicator for inter node master ranks
            if (flag2) MPI_Ibcast(buf, len, MPI_DOUBLE, 0, inter_node, &req2); // 2nd level bcast
        }

        if (ppn > 1){
            int intra_node_ranks[ppn];
            int base_rank = (rank/ppn)*ppn;
            for (int i=0; i<ppn; i++){ // find all intra node ranks
                intra_node_ranks[i] = base_rank+i;
            }
            MPI_Group_incl(world_group, ppn, intra_node_ranks, &intra_node_group);
            MPI_Comm_create_group(MPI_COMM_WORLD, intra_node_group, 0, &intra_node); // communicator for intra node communication
            MPI_Ibcast(buf, len, MPI_DOUBLE, 0, intra_node, &req3); // 3rd level bcast
        }

        if (groups > 1 && flag1){
            MPI_Comm_free(&inter_group);
            MPI_Wait(&req1, MPI_STATUS_IGNORE);
        }
        if (nodes_per_group[group_id] > 1 && flag2){
            MPI_Comm_free(&inter_node);
            MPI_Wait(&req2, MPI_STATUS_IGNORE);
        }
        if (ppn > 1){
            MPI_Comm_free(&intra_node);
            MPI_Wait(&req3, MPI_STATUS_IGNORE);
        }
    
        time[itr] = MPI_Wtime()-time[itr];
    }
    free(buf);
    double max_time[5];
    MPI_Reduce(time, max_time, 5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return avg(max_time);
}

double standard_reduce(){
    double* sendbuf = random_init(len); // randomly initialise sendbuf
    double* recvbuf;
    if (!rank) recvbuf = (double*)malloc(len*sizeof(double)); // initialise recvbuf for only root process
    double time[5];
    int itr = 5;
    while (itr--){
        MPI_Barrier(MPI_COMM_WORLD); // for fair time measurement
        time[itr] = MPI_Wtime();

        MPI_Reduce(sendbuf, recvbuf, len, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        time[itr] = MPI_Wtime()-time[itr];
    }
    free(sendbuf); 
    if (!rank) free(recvbuf);
    double max_time[5];
    MPI_Reduce(time, max_time, 5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return avg(max_time);
}

double optimized_reduce(int* nodes_per_group){
    double* sendbuf = random_init(len); // randomly initialise sendbuf
    double* recvbuf;
    if (!rank) recvbuf = (double*)malloc(len*sizeof(double)); // initialise recvbuf for only root process
    double time[5];
    int itr = 5;
    while (itr--){
        MPI_Barrier(MPI_COMM_WORLD); // for fair measurement
        time[itr] = MPI_Wtime();

        MPI_Request req1, req2, req3;
        MPI_Group inter_group_group, inter_node_group, intra_node_group;
        MPI_Comm inter_group, inter_node, intra_node;
        int flag1 = 0; int flag2 = 0;
        int group_id;
        double* buf1 = NULL; // level 1 intermediary buf (inter node master ranks may have data in this)
        double* buf2 = NULL; // level 2 intermediary buf (inter group master ranks may have data in this)

        MPI_Group world_group;
        MPI_Comm_group(MPI_COMM_WORLD, &world_group);

        int low = 0;
        for (int i=0; i<groups; i++){ // need to find group id first as order is reverse as that of bcast
            if (rank >= low && rank < low+nodes_per_group[i]*ppn) group_id = i;
            low += nodes_per_group[i]*ppn;
        }

        if (ppn > 1){
            int intra_node_ranks[ppn];
            int base_rank = (rank/ppn)*ppn;
            for (int i=0; i<ppn; i++){ // find all intra node ranks
                intra_node_ranks[i] = base_rank+i;
            }
            MPI_Group_incl(world_group, ppn, intra_node_ranks, &intra_node_group);
            MPI_Comm_create_group(MPI_COMM_WORLD, intra_node_group, 0, &intra_node); // communicator for intra node communication
            if (nodes_per_group[group_id] > 1){
                if (rank == base_rank) buf1 = (double*)malloc(len*sizeof(double));
                MPI_Ireduce(sendbuf, buf1, len, MPI_DOUBLE, MPI_MAX, 0, intra_node, &req3); // multiple nodes in this group, send reduced data 1 level up to buf1
            }
            else {
                if (groups > 1){ // only 1 group in this node, send reduced data 2 levels up to buf2
                    if (rank == base_rank) buf2 = (double*)malloc(len*sizeof(double));
                    MPI_Ireduce(sendbuf, buf2, len, MPI_DOUBLE, MPI_MAX, 0, intra_node, &req3);
                }
                else MPI_Ireduce(sendbuf, recvbuf, len, MPI_DOUBLE, MPI_MAX, 0, intra_node, &req3); // single group and single node in this group, we are done
            }
        }

        if (nodes_per_group[group_id] > 1){
            flag2 = rank%ppn == 0 ? 1 : 0; // flag2 = 1 for master ranks at inter node level (within a group)
            int inter_node_ranks[nodes_per_group[group_id]];
            int offset = 0;
            for (int i=0; i<group_id; i++){
                offset += nodes_per_group[i];
            }
            for (int i=0; i<nodes_per_group[group_id]; i++){
                inter_node_ranks[i] = ppn*i + offset*ppn; 
            }
            MPI_Group_incl(world_group, nodes_per_group[group_id], inter_node_ranks, &inter_node_group);
            MPI_Comm_create_group(MPI_COMM_WORLD, inter_node_group, 0, &inter_node); // communicator for inter node master ranks
            if (flag2){
                if (ppn > 1){
                    if (groups > 1){
                        if (rank == offset*ppn) buf2 = (double*)malloc(len*sizeof(double));
                        MPI_Ireduce(buf1, buf2, len, MPI_DOUBLE, MPI_MAX, 0, inter_node, &req2); // all three levels exist, send reduced data 1 level up to buf2
                    }  
                    else MPI_Ireduce(buf1, recvbuf, len, MPI_DOUBLE, MPI_MAX, 0, inter_node, &req2);  // single group, we are done
                } 
                else {
                    if (groups > 1){
                        if (rank == offset*ppn) buf2 = (double*)malloc(len*sizeof(double));
                        MPI_Ireduce(sendbuf, buf2, len, MPI_DOUBLE, MPI_MAX, 0, inter_node, &req2); // multiple groups, send reduced data 1 level up to buf2
                    }  
                    else MPI_Ireduce(sendbuf, recvbuf, len, MPI_DOUBLE, MPI_MAX, 0, inter_node, &req2); // single group and ppn = 1, we are done
                }
            }
        }

        if (groups > 1){
            int inter_group_ranks[groups];
            int exp_rank = 0;
            for (int i=0; i<groups; i++){
                if (rank == exp_rank) flag1 = 1; // flag1 = 1 for master ranks at inter group level
                inter_group_ranks[i] = exp_rank;
                exp_rank += nodes_per_group[i]*ppn;
            }
            MPI_Group_incl(world_group, groups, inter_group_ranks, &inter_group_group);
            MPI_Comm_create_group(MPI_COMM_WORLD, inter_group_group, 0, &inter_group); // communicator for inter group master ranks
            if (flag1){
                if (!(ppn > 1 || nodes_per_group[group_id] > 1)) MPI_Ireduce(sendbuf, recvbuf, len, MPI_DOUBLE, MPI_MAX, 0, inter_group, &req1); // ppn = 1, single node in this group, we are done
                else MPI_Ireduce(buf2, recvbuf, len, MPI_DOUBLE, MPI_MAX, 0, inter_group, &req1); // inter group master ranks have data in buf2, send reduced data to root node
            } 
        }

        if (ppn > 1){
            MPI_Comm_free(&intra_node);
            MPI_Wait(&req3, MPI_STATUS_IGNORE);
        }
        if (nodes_per_group[group_id] > 1 && flag2){
            MPI_Comm_free(&inter_node);
            MPI_Wait(&req2, MPI_STATUS_IGNORE);
        }
        if (groups > 1 && flag1){
            MPI_Comm_free(&inter_group);
            MPI_Wait(&req1, MPI_STATUS_IGNORE);
        }

        if (buf1) free(buf1);
        if (buf2) free(buf2);

        time[itr] = MPI_Wtime()-time[itr];
    }
    free(sendbuf); 
    if (!rank) free(recvbuf);
    double max_time[5];
    MPI_Reduce(time, max_time, 5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return avg(max_time);
}

double standard_gather(){
    double* sendbuf = random_init(len); // randomly initialise sendbuf
    double* recvbuf;
    if (!rank) recvbuf = random_init(len*size); // initialise recvbuf for only root process
    double time[5];
    int itr = 5;
    while (itr--){
        MPI_Barrier(MPI_COMM_WORLD); // for fair measurement
        time[itr] = MPI_Wtime();

        MPI_Gather(sendbuf, len, MPI_DOUBLE, recvbuf, len, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        time[itr] = MPI_Wtime()-time[itr];
    }
    free(sendbuf); 
    if (!rank) free(recvbuf);
    double max_time[5];
    MPI_Reduce(time, max_time, 5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return avg(max_time);
}

double optimized_gather(int* nodes_per_group){
    double* sendbuf = random_init(len); // randomly initialise sendbuf
    double* recvbuf;
    if (!rank) recvbuf = random_init(len*size); // initialise recvbuf for only root process
    double time[5];
    int itr = 5;
    while (itr--){
        MPI_Barrier(MPI_COMM_WORLD); // for fair measurement
        time[itr] = MPI_Wtime();

        int reordered_ranks[size];
        binomial_gather_rank_reorder(reordered_ranks); // get reordered ranks

        MPI_Group world_group;
        MPI_Comm_group(MPI_COMM_WORLD, &world_group);
        MPI_Group gather_group;
        MPI_Group_incl(world_group, size, reordered_ranks, &gather_group);
        MPI_Comm gather_comm;
        MPI_Comm_create_group(MPI_COMM_WORLD, gather_group, 0, &gather_comm); // create new global communicator with reordered ranks
    
        MPI_Request req;
        MPI_Igather(sendbuf, len, MPI_DOUBLE, recvbuf, len, MPI_DOUBLE, 0, gather_comm, &req);
        MPI_Wait(&req, MPI_STATUS_IGNORE);

        MPI_Comm_free(&gather_comm);
        
        if (!rank){ // correct the order of data in recvbuf (data received in jumbled order due to reordered ranks)
            int check[size];
            for (int i=0; i<size; i++){
                if (reordered_ranks[i] == i) check[i] = 1;
                else check[i] = 0;
            }
            for (int i=0; i<size; i++){
                if (!check[i]){
                    check[i] = 1;
                    check[reordered_ranks[i]] = 1;
                    swap_data(recvbuf, i*len, reordered_ranks[i]*len);
                }
            }
        }
        
        time[itr] = MPI_Wtime()-time[itr];
    }
    free(sendbuf); 
    if (!rank) free(recvbuf);
    double max_time[5];
    MPI_Reduce(time, max_time, 5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return avg(max_time);
}

double standard_alltoallv(int* sendcounts){
    double* sendbuf = random_init(len*size); // randomly initialise sendbuf
    double* recvbuf = (double*)malloc(len*size*sizeof(double)); // allocate recvbuf
    int sdispls[size], recvcounts[size], rdispls[size];
    for (int i=0; i<size; i++){
        sdispls[i] = i*len;
        recvcounts[i] = len;
        rdispls[i] = i*len;
    }
    double time[5];
    int itr = 5;
    while (itr--){
        MPI_Barrier(MPI_COMM_WORLD); // for fair measurement
        time[itr] = MPI_Wtime();

        MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_DOUBLE, recvbuf, recvcounts, rdispls, MPI_DOUBLE, MPI_COMM_WORLD);

        time[itr] = MPI_Wtime()-time[itr];
    }
    free(sendbuf); free(recvbuf);
    double max_time[5];
    MPI_Reduce(time, max_time, 5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return avg(max_time);
}

double optimized_alltoallv(int* sendcounts){
    double* sendbuf = random_init(len*size); // randomly initialise sendbuf
    double* recvbuf = (double*)malloc(len*size*sizeof(double)); // allocate recvbuf
    int sdispls[size], recvcounts[size], rdispls[size];
    for (int i=0; i<size; i++){
        sdispls[i] = i*len;
        recvcounts[i] = len;
        rdispls[i] = i*len;
    }
    double time[5];
    int itr = 5;
    while (itr--){
        MPI_Barrier(MPI_COMM_WORLD); // for fair measurement
        time[itr] = MPI_Wtime();
        
        if (D < 1*KB){ // for small data, use scatterv
            MPI_Request req[size];
            for (int i=0; i<size; i++){
                MPI_Iscatterv(sendbuf, sendcounts, sdispls, MPI_DOUBLE, recvbuf + rdispls[i], recvcounts[i], MPI_DOUBLE, i, MPI_COMM_WORLD, &req[i]);
            }
            MPI_Waitall(size, req, MPI_STATUS_IGNORE);
        }
        else { // else standard call
            MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_DOUBLE, recvbuf, recvcounts, rdispls, MPI_DOUBLE, MPI_COMM_WORLD);
        }

        time[itr] = MPI_Wtime()-time[itr];
    }
    free(sendbuf); free(recvbuf);
    double max_time[5];
    MPI_Reduce(time, max_time, 5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return avg(max_time);
}

int main(int argc, char* argv[]){
    D = atoi(argv[1]); // data size in KB
    ppn = atoi(argv[2]); // ppn
    groups = argc-3; // number of groups (depends on hostfile)
    int nodes_per_group[groups];
    for (int i=3; i<argc; i++) nodes_per_group[i-3] = atoi(argv[i]); // nodes in every group

    len = (D*KB)/sizeof(double); // length of a double array corresponding to D KB

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    double time_bcast = standard_bcast();
    double time_bcast_opt = optimized_bcast(nodes_per_group);

    double time_reduce = standard_reduce();
    double time_reduce_opt = optimized_reduce(nodes_per_group);

    double time_gather = standard_gather();
    double time_gather_opt = optimized_gather(nodes_per_group);

    // reduce data size to 1024 for alltoallv
    if (D >= 2048){
        D = 1024;
        len = (D*KB)/sizeof(double);
    }
    int sendcounts[size];
    for (int i=0; i<size; i++) sendcounts[i] = rand()%(len+1); // randomize sendcounts and pass to both the calls for fair comparison

    double time_alltoallv = standard_alltoallv(sendcounts);
    double time_alltoallv_opt = optimized_alltoallv(sendcounts);

    // print the results
    if (!rank){
        printf ("%f %f\n", time_bcast, time_bcast_opt);
        printf ("%f %f\n", time_reduce, time_reduce_opt);
        printf ("%f %f\n", time_gather, time_gather_opt);
        printf ("%f %f", time_alltoallv, time_alltoallv_opt);
    }
    
    MPI_Finalize();
    return 0;
}