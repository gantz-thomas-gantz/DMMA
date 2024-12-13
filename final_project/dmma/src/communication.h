#include <mpi.h>
#include <stdlib.h>

#include "utilities.h"

int exchange(u64 **Z_recv, struct u64_darray *Z, int p, int my_rank) {
	/* STEP 1:
	 * Send z dynamic array values to respective processes.
	 * a) send sizes first.
	 * b) send actual zs from dynamic arrays.
	 * Immediate sends to avoid deadlocks.*/

	MPI_Request send_request[p];
	for (int dest_rank = 0; dest_rank < p && dest_rank != my_rank;
	     dest_rank++) {
		int sender_packet[2];
		sender_packet[0] = my_rank;
		sender_packet[1] = Z[dest_rank]->size;
		MPI_Isend(sender_packet, 2, MPI_INT, dest_rank, 0,
			  MPI_COMM_WORLD);
		MPI_Isend(Z[dest_rank]->data, Z[dest_rank]->size, MPI_UINT64_T,
			  dest_rank, 1, MPI_COMM_WORLD,
			  send_request[dest_rank]);
	}

	/* STEP 2:
	 * Receive zs to look up in own dictionary.
	 * a) receive sizes of things to receive.
	 * b) allocate static array to contain all received zs.
	 * c) receive actual zs and store them into static array. */

	MPI_Request request;
	int recv_sizes[p];
	int recv_total_size = 0;
	for (int sender = 0; sender < p - 1; sender++) {
		// sender_packet[0]: source, sender_packet[1]: msg size
		int sender_packet[2];
		MPI_Irecv(sender_packet, 2, MPI_INT, MPI_ANY_SOURCE, 0,
			  MPI_COMM_WORLD, request);
		MPI_Wait(&request, MPI_STATUS_IGNORE);
		recv_sizes[sender_packet[0]] = sender_packet[1];
		recv_total_size += sender_packet[1];
	}
	int curr_size = 0;
	Z_recv_size = recv_total_size + Z[my_rank].size;
	*Z_recv = (u64 *)malloc(sizeof(u64) * Z_recv_size);
	// Receive from others.
	MPI_Status recv_requests[p];
	for (int sender = 0; sender < p && sender != my_rank; sender++) {
		MPI_Irecv(Z_recv + curr_size, recv_sizes[sender], MPI_UINT64_T,
			  sender, 0, MPI_COMM_WORLD, recv_requests[sender]);
		curr_size += recv_sizes[sender];
	}

	/* STEP 3:
	 * Free dynamic arrays of sent zs only when finished sending.*/

	for (int dest_rank = 0; dest_rank < p && dest_rank != my_rank;
	     dest_rank++) {
		MPI_Wait(send_request[dest_rank], MPI_STATUS_IGNORE);
		free_u64_darray(Z[dest_rank]);
	}
	memcpy(Z_recv + curr_size, Z[my_rank].data,
	       sizeof(u64) * Z[my_rank].size);
	free_u64_darray(Z[my_rank]);
	free(Z);

	/* STEP 4:
	 * Wait until all data is received.*/

	for (int sender = 0; sender < p && sender != my_rank; sender++)
		MPI_Wait(recv_requests[sender], MPI_STATUS_IGNORE);

	return Z_recv_size;
}
