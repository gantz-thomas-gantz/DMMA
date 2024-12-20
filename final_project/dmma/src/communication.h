#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "utilities.h"

// TODO: Send in smaller batches to occupy less memory.

int exchange(u64 **recv, struct u64_darray *send, int p, int my_rank) {
	/* STEP 1:
	 * a) Send sizes of data to be sent.
	 * b) Compute contiguous send displacements.
	 */

	MPI_Request size_request[p];
	int send_sizes[p];
	int send_total_size = 0;
	int sender_packet[2 * (p - 1)];
	int i = 0;
	for (int dest_rank = 0; dest_rank < p; dest_rank++) {
		if (dest_rank == my_rank) continue;
		sender_packet[i + 0] = my_rank;
		sender_packet[i + 1] = send[dest_rank].size;
		MPI_Isend(sender_packet + i, 2, MPI_INT, dest_rank, 0,
			  MPI_COMM_WORLD, &(size_request[dest_rank]));
		send_sizes[dest_rank] = send[dest_rank].size;
		send_total_size += send_sizes[dest_rank];
		i += 2;
	}
	send_sizes[my_rank] = send[my_rank].size;
	send_total_size += send_sizes[my_rank];
	int send_displacements[p];
	send_displacements[0] = 0;
	for (int dest_rank = 1; dest_rank < p; dest_rank++)
		send_displacements[dest_rank] =
		    send_displacements[dest_rank - 1] +
		    send_sizes[dest_rank - 1];

	/* STEP 2:
	 * a) receive sizes of things to receive.
	 * b) compute contiguous receive displacements.
	 * b) allocate static array to store all received zs.
	 */

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
	recv_sizes[my_rank] = send[my_rank].size;
	recv_total_size += recv_sizes[my_rank];
	int recv_displacements[p];
	recv_displacements[0] = 0;
	for (int sender = 1; sender < p; sender++)
		recv_displacements[sender] =
		    recv_displacements[sender - 1] + recv_sizes[sender - 1];

	/* STEP 3:
	 * a) Store dynamic array send values contiguously into a single static
	 * array.
	 * b) Send and receive data.
	 * */
	u64 *contiguous_send = malloc(send_total_size * sizeof(u64));
#pragma omp parallel for
	for (int sender = 0; sender < p; sender++) {
		memcpy(contiguous_send + send_displacements[sender],
		       send[sender].data, send[sender].size * sizeof(u64));
		free_u64_darray(&(send[sender]));
	}
	// This commented version of stepwise allocation of the contiguous send
	// array is better for storage but is slower timewise (because of
	// realloc)
	/*u64 *contiguous_send = malloc(send[0].size * sizeof(u64));
	for (int sender = 0; sender < p; sender++) {
		memcpy(contiguous_send + send_displacements[sender],
		       send[sender].data, send[sender].size * sizeof(u64));
		free_u64_darray(&(send[sender]));
		if (sender != p - 1)
			contiguous_send = realloc(
			    contiguous_send, (send_displacements[sender + 1] +
					      send[sender + 1].size) *
						 sizeof(u64));
	}*/

	// Malloc using alignment for subsequent vectorization
	*recv = aligned_alloc(32, sizeof(u64) * recv_total_size);
	MPI_Alltoallv(contiguous_send, send_sizes, send_displacements,
		      MPI_UINT64_T, *recv, recv_sizes, recv_displacements,
		      MPI_UINT64_T, MPI_COMM_WORLD);
	free(contiguous_send);
	return recv_total_size;
}
