#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "utilities.h"

// TODO: Send in smaller batches to occupy less memory.
int exchange(u64 **recv, struct u64_darray *send, int p, int my_rank) {
	/* STEP 1:
	 * Send z dynamic array values to respective processes.
	 * a) send sizes first.
	 * b) send actual zs from dynamic arrays.
	 * Immediate sends to avoid deadlocks.*/

	MPI_Request send_request[p];
	MPI_Request size_request[p];
	int sender_packet[2 * (p - 1)];
	int i = 0;
	for (int dest_rank = 0; dest_rank < p; dest_rank++) {
		if (dest_rank == my_rank) continue;
		sender_packet[i + 0] = my_rank;
		sender_packet[i + 1] = send[dest_rank].size;
		MPI_Isend(sender_packet + i, 2, MPI_INT, dest_rank, 0,
			  MPI_COMM_WORLD, &(size_request[dest_rank]));
		MPI_Isend(send[dest_rank].data, send[dest_rank].size,
			  MPI_UINT64_T, dest_rank, 1, MPI_COMM_WORLD,
			  &(send_request[dest_rank]));
		i += 2;
	}

	/* STEP 2:
	 * Receive zs to look up in own dictionary.
	 * a) receive sizes of things to receive.
	 * b) allocate static array to contain all received zs.
	 * c) receive actual zs and store them into static array. */

	MPI_Request recv_size_request;
	int recv_sizes[p];
	int recv_total_size = 0;
	for (int sender = 0; sender < p; sender++) {
		// sender_packet[0]: source, sender_packet[1]: msg size
		if (sender == my_rank) continue;
		int recv_packet[2];
		MPI_Irecv(recv_packet, 2, MPI_INT, MPI_ANY_SOURCE, 0,
			  MPI_COMM_WORLD, &recv_size_request);
		MPI_Wait(&recv_size_request, MPI_STATUS_IGNORE);
		recv_sizes[recv_packet[0]] = recv_packet[1];
		recv_total_size += recv_packet[1];
	}
	int curr_size = 0;
	int recv_size = recv_total_size + send[my_rank].size;
	*recv = (u64 *)malloc(sizeof(u64) * recv_size);
	// Receive from others.
	MPI_Request recv_requests[p];
	for (int sender = 0; sender < p; sender++) {
		if (sender == my_rank) continue;
		MPI_Irecv(recv + curr_size, recv_sizes[sender], MPI_UINT64_T,
			  sender, 1, MPI_COMM_WORLD, &(recv_requests[sender]));
		curr_size += recv_sizes[sender];
	}

	/* STEP 3:
	 * Free dynamic arrays of sent zs only when finished sending.*/

	for (int dest_rank = 0; dest_rank < p; dest_rank++) {
		if (dest_rank == my_rank) continue;
		MPI_Wait(&(size_request[dest_rank]), MPI_STATUS_IGNORE);
		MPI_Wait(&(send_request[dest_rank]), MPI_STATUS_IGNORE);
		free_u64_darray(&(send[dest_rank]));
	}
	memcpy(*recv + curr_size, send[my_rank].data,
	       sizeof(u64) * send[my_rank].size);
	free_u64_darray(&(send[my_rank]));
	/* STEP 4:
	 * Wait until all data is received.*/

	for (int sender = 0; sender < p; sender++) {
		if (sender == my_rank) continue;
		MPI_Wait(&(recv_requests[sender]), MPI_STATUS_IGNORE);
	}
	return recv_size;
}
