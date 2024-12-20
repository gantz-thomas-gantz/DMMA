#include <assert.h>
#include <err.h>
#include <getopt.h>
#include <inttypes.h>
#include <mpi.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "communication.h"
#include "utilities.h"

typedef uint64_t u64; /* portable 64-bit integer */
typedef uint32_t u32; /* portable 32-bit integer */
struct __attribute__((packed)) entry {
	u32 k;
	u64 v;
}; /* hash table entry */

/***************************** global variables ******************************/

u64 n = 0; /* block size (in bits) */
u64 mask;  /* this is 2**n - 1 */

u64 dict_size;	 /* number of slots in the hash table */
struct entry *A; /* the hash table */

/* (P, C) : two plaintext-ciphertext pairs */
u32 P[2][2] = {{0, 0}, {0xffffffff, 0xffffffff}};
u32 C[2][2];

/************************ tools and utility functions *************************/

double wtime() {
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + ts.tv_usec / 1E6;
}

// murmur64 hash functions, tailorized for 64-bit ints / Cf. Daniel Lemire
u64 murmur64(u64 x) {
	x ^= x >> 33;
	x *= 0xff51afd7ed558ccdull;
	x ^= x >> 33;
	x *= 0xc4ceb9fe1a85ec53ull;
	x ^= x >> 33;
	return x;
}

/* represent n in 4 bytes */
void human_format(u64 n, char *target) {
	if (n < 1000) {
		sprintf(target, "%" PRId64, n);
		return;
	}
	if (n < 1000000) {
		sprintf(target, "%.1fK", n / 1e3);
		return;
	}
	if (n < 1000000000) {
		sprintf(target, "%.1fM", n / 1e6);
		return;
	}
	if (n < 1000000000000ll) {
		sprintf(target, "%.1fG", n / 1e9);
		return;
	}
	if (n < 1000000000000000ll) {
		sprintf(target, "%.1fT", n / 1e12);
		return;
	}
}

/******************************** SPECK block cipher **************************/

#define ROTL32(x, r) (((x) << (r)) | (x >> (32 - (r))))
#define ROTR32(x, r) (((x) >> (r)) | ((x) << (32 - (r))))

#define ER32(x, y, k) \
	(x = ROTR32(x, 8), x += y, x ^= k, y = ROTL32(y, 3), y ^= x)
#define DR32(x, y, k) \
	(y ^= x, y = ROTR32(y, 3), x ^= k, x -= y, x = ROTL32(x, 8))

void Speck64128KeySchedule(const u32 K[], u32 rk[]) {
	u32 i, D = K[3], C = K[2], B = K[1], A = K[0];
	for (i = 0; i < 27;) {
		rk[i] = A;
		ER32(B, A, i++);
		rk[i] = A;
		ER32(C, A, i++);
		rk[i] = A;
		ER32(D, A, i++);
	}
}

void Speck64128Encrypt(const u32 Pt[], u32 Ct[], const u32 rk[]) {
	u32 i;
	Ct[0] = Pt[0];
	Ct[1] = Pt[1];
	for (i = 0; i < 27;) ER32(Ct[1], Ct[0], rk[i++]);
}

void Speck64128Decrypt(u32 Pt[], const u32 Ct[], u32 const rk[]) {
	int i;
	Pt[0] = Ct[0];
	Pt[1] = Ct[1];
	for (i = 26; i >= 0;) DR32(Pt[1], Pt[0], rk[i--]);
}

/******************************** dictionary ********************************/

/*
 * "classic" hash table for 64-bit key-value pairs, with linear probing.
 * It operates under the assumption that the keys are somewhat random 64-bit
 * integers. The keys are only stored modulo 2**32 - 5 (a prime number), and
 * this can lead to some false positives.
 */
static const u32 EMPTY = 0xffffffff;  // Biggest 32 bit number
static const u64 PRIME = 0xfffffffb;

/* allocate a hash table with `size` slots (12*size bytes) */
void dict_setup(u64 size) {
	dict_size = size;
	char hdsize[8];
	human_format(dict_size * sizeof(*A), hdsize);
	printf("Dictionary size: %sB\n", hdsize);

	A = malloc(sizeof(*A) * dict_size);
	if (A == NULL) err(1, "impossible to allocate the dictionnary");
	for (u64 i = 0; i < dict_size; i++) A[i].k = EMPTY;
}

/* Insert the binding key |----> value in the dictionnary */
void dict_insert(u64 key, u64 value) {
	u64 h = murmur64(key) % dict_size;
	// printf("key: %zu; val: %zu \n", key, value);
	// printf("hash_val: %zu\n", murmur64(key));
	// printf("h: %zu\n", h);
	for (;;) {
		if (A[h].k == EMPTY) break;
		h += 1;
		if (h == dict_size) h = 0;
	}
	assert(A[h].k == EMPTY);
	A[h].k = key % PRIME;
	A[h].v = value;
}

/* Query the dictionnary with this `key`.  Write values (potentially)
 *  matching the key in `values` and return their number. The `values`
 *  array must be preallocated of size (at least) `maxval`.
 *  The function returns -1 if there are more than `maxval` results.
 */
int dict_probe(u64 key, int maxval, u64 values[]) {
	u32 k = key % PRIME;
	u64 h = murmur64(key) % dict_size;
	int nval = 0;
	for (;;) {
		if (A[h].k == EMPTY) return nval;
		if (A[h].k == k) {
			if (nval == maxval) return -1;
			values[nval] = A[h].v;
			nval += 1;
		}
		h += 1;
		if (h == dict_size) h = 0;
	}
}

/***************************** MITM problem ***********************************/

/* f : {0, 1}^n --> {0, 1}^n.  Speck64-128 encryption of P[0], using k */
u64 f(u64 k) {
	assert((k & mask) == k);
	u32 K[4] = {k & 0xffffffff, k >> 32, 0, 0};
	u32 rk[27];
	Speck64128KeySchedule(K, rk);
	u32 Ct[2];
	Speck64128Encrypt(P[0], Ct, rk);
	return ((u64)Ct[0] ^ ((u64)Ct[1] << 32)) & mask;
}

/* g : {0, 1}^n --> {0, 1}^n.  speck64-128 decryption of C[0], using k */
u64 g(u64 k) {
	assert((k & mask) == k);
	u32 K[4] = {k & 0xffffffff, k >> 32, 0, 0};
	u32 rk[27];
	Speck64128KeySchedule(K, rk);
	u32 Pt[2];
	Speck64128Decrypt(Pt, C[0], rk);
	return ((u64)Pt[0] ^ ((u64)Pt[1] << 32)) & mask;
}

bool is_good_pair(u64 k1, u64 k2) {
	u32 Ka[4] = {k1 & 0xffffffff, k1 >> 32, 0, 0};
	u32 Kb[4] = {k2 & 0xffffffff, k2 >> 32, 0, 0};
	u32 rka[27];
	u32 rkb[27];
	Speck64128KeySchedule(Ka, rka);
	Speck64128KeySchedule(Kb, rkb);
	u32 mid[2];
	u32 Ct[2];
	Speck64128Encrypt(P[1], mid, rka);
	Speck64128Encrypt(mid, Ct, rkb);
	return (Ct[0] == C[1][0]) && (Ct[1] == C[1][1]);
}

/******************************************************************************/

/* search the "golden collision" */
int golden_claw_search(int maxres, u64 **K1, u64 **K2, int my_rank, int p) {
	double start = wtime();
	u64 N = 1ull << n;

	/* STEP 1: Build partitioned dictionary.
	 * Logical step: Process i holds all keys respective to their rank
	 * (key % i == 0).*/

	struct u64_darray *dict_zx = malloc(p * sizeof(struct u64_darray));

#pragma omp parallel for
	for (int rank = 0; rank < p; rank++)
		initialize_u64_darray(
		    &(dict_zx[rank]),
		    2 * N /
			(p * p));  // Expected amount of z and x for each rank

	omp_lock_t locks[p];
	for (int i = 0; i < p; i++) {
		omp_init_lock(&locks[i]);
	}
#pragma omp parallel
	{
#pragma omp for
		for (u64 x = (my_rank * N) / p; x < ((my_rank + 1) * N) / p;
		     x++) {
			// TODO: AVX
			u64 z = f(x);
			int dest_rank = z % p;
			omp_set_lock(&locks[dest_rank]);
			append(&(dict_zx[dest_rank]), z);
			append(&(dict_zx[dest_rank]), x);
			omp_unset_lock(&locks[dest_rank]);
		}
	}

	u64 *dict_zx_recv;
	int dict_zx_recv_size = exchange(&dict_zx_recv, dict_zx, p, my_rank);
	free(dict_zx);

	for (int i = 0; i < dict_zx_recv_size - 1; i += 2)
		dict_insert(dict_zx_recv[i], dict_zx_recv[i + 1]);
	free(dict_zx_recv);

	double mid = wtime();
	printf("Fill: %.1fs\n", mid - start);

	/* STEP 2:
	 * Attribute zs to their repective processes.
	 * Logical step: Look for the key in the dictionary of the process who
	 * could potentially have it (dest_rank = z % p).
	 * Since we don't know exactly how many there are per process, we use
	 * dynammic arrays and initialize them to their expected capacity
	 * (N/(p*p)). */

	struct u64_darray *dict_zy = malloc(p * sizeof(struct u64_darray));
#pragma omp parallel for
	for (int rank = 0; rank < p; rank++)
		initialize_u64_darray(
		    &(dict_zy[rank]),
		    2 * N / (p * p));  // Expected amount of z for each rank

#pragma omp parallel
	{
#pragma omp for
		for (u64 y = (my_rank * N) / p; y < ((my_rank + 1) * N) / p;
		     y++) {
			// TODO: AVX
			u64 z = g(y);
			int dest_rank = z % p;
			// Lock specific to dest_rank
			omp_set_lock(&locks[dest_rank]);
			append(&(dict_zy[dest_rank]), z);
			append(&(dict_zy[dest_rank]), y);
			omp_unset_lock(&locks[dest_rank]);
		}
	}
	for (int i = 0; i < p; i++) {
		omp_destroy_lock(&locks[i]);
	}

	u64 *dict_zy_recv;

	int dict_zy_recv_size = exchange(&dict_zy_recv, dict_zy, p, my_rank);
	free(dict_zy);
	/* STEP 7:
	 * Now look up my zs in my local dictionaries.*/

	int nres = 0;
	u64 ncandidates = 0;
	u64 k1[16], k2[16];
	int flag = 0;
#pragma omp parallel shared(flag, nres, ncandidates, k1, k2)
	{
		printf("n threads: %d\n", omp_get_num_threads());
		u64 x[256 / omp_get_num_threads()];
#pragma omp for
		for (int i = 0; i < dict_zy_recv_size; i += 2) {
			// printf("Rank: %d. Starting probe.\n", my_rank);
			u64 z = dict_zy_recv[i];
			u64 y = dict_zy_recv[i + 1];
			int nx = dict_probe(z, 256 / omp_get_num_threads(), x);
			assert(nx >= 0);
#pragma omp atomic
			ncandidates += nx;
			for (int i = 0; i < nx; i++) {
				if (is_good_pair(x[i], y)) {
					if (flag) continue;
					if (nres < maxres) {
						printf(
						    "SOLUTION FOUND! "
						    "by thread %d\n",
						    omp_get_thread_num());

#pragma omp critical
						{
							k1[nres] = x[i];
							k2[nres] = y;
							nres++;
						}
					} else {
						flag = 1;  // Set flag to
							   // stop further
							   // iterations
					}
				}
			}
		}
	}
	if (flag) return -1;
	free(dict_zy_recv);

	/* STEP 8:
	 * Gather all locally found solutions into root.*/

	int *global_nres =
	    NULL;	     // Array to store nres values from each process
	int *displs = NULL;  // Array to store displacements for gathered data
	// Root process allocates space for recv_counts and displs
	if (my_rank == 0) {
		global_nres = malloc(p * sizeof(int));
		displs = malloc(p * sizeof(int));
	}
	// Gather nres from all processes to the root
	MPI_Gather(&nres, 1, MPI_INT, global_nres, 1, MPI_INT, 0,
		   MPI_COMM_WORLD);
	// Calculate displacements and total size of gathered arrays on
	// the root process
	int total_nres = 0;
	if (my_rank == 0) {
		displs[0] = 0;
		for (int i = 0; i < p; i++) {
			total_nres += global_nres[i];
			if (i > 0)
				displs[i] = displs[i - 1] + global_nres[i - 1];
		}
		// Allocate global arrays on the root process
		*K1 = malloc(total_nres * sizeof(u64));
		*K2 = malloc(total_nres * sizeof(u64));
	}
	// Gather k1 and k2 from all processes into K1 and K2 on the
	// root
	MPI_Gatherv(k1, nres, MPI_UNSIGNED_LONG_LONG, *K1, global_nres, displs,
		    MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Gatherv(k2, nres, MPI_UNSIGNED_LONG_LONG, *K2, global_nres, displs,
		    MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	// Root process can now use the gathered results
	if (my_rank == 0) {
		printf("Total results gathered: %d\n", total_nres);
		for (int i = 0; i < total_nres; i++) {
			printf("K1[%d] = %lu, K2[%d] = %lu\n", i, *K1[i], i,
			       *K2[i]);
		}

		/*
		printf("Probe: %.1fs. %" PRId64 " candidate pairs
		tested\n", wtime() - mid, ncandidates);
		*/
	}
	return total_nres;
}
/************************** command-line options
 * ****************************/

void usage(char **argv) {
	printf("%s [OPTIONS]\n\n", argv[0]);
	printf("Options:\n");
	printf("--n N                       block size [default 24]\n");
	printf("--C0 N                      1st ciphertext (in hex)\n");
	printf("--C1 N                      2nd ciphertext (in hex)\n");
	printf("\n");
	printf("All arguments are required\n");
	exit(0);
}

void process_command_line_options(int argc, char **argv) {
	struct option longopts[4] = {{"n", required_argument, NULL, 'n'},
				     {"C0", required_argument, NULL, '0'},
				     {"C1", required_argument, NULL, '1'},
				     {NULL, 0, NULL, 0}};
	char ch;
	int set = 0;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
			case 'n':
				n = atoi(optarg);
				mask = (1ull << n) - 1;
				break;
			case '0':
				set |= 1;
				u64 c0 = strtoull(optarg, NULL, 16);
				C[0][0] = c0 & 0xffffffff;
				C[0][1] = c0 >> 32;
				break;
			case '1':
				set |= 2;
				u64 c1 = strtoull(optarg, NULL, 16);
				C[1][0] = c1 & 0xffffffff;
				C[1][1] = c1 >> 32;
				break;
			default:
				errx(1, "Unknown option\n");
		}
	}
	if (n == 0 || set != 3) {
		usage(argv);
		exit(1);
	}
}
/******************************************************************************/

int main(int argc, char **argv) {
	MPI_Init(NULL, NULL);
	process_command_line_options(argc, argv);
	printf("Running with n=%d, C0=(%08x, %08x) and C1=(%08x, %08x)\n",
	       (int)n, C[0][0], C[0][1], C[1][0], C[1][1]);
	// Set number of OpenMP threads
	// omp_set_num_threads(3);
	int my_rank;
	int p;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	printf("p=%d\n", p);
	u64 s = (u64)(n < 10 ? 1.125 * (1ull << n) : 1.125 * (1ull << n) / p);
	dict_setup(s);	// If not divisable no problem because hash table
			// has more space allocated than necessary.
	/* search */
	u64 *K1, *K2;
	double start_time, end_time;

	start_time = MPI_Wtime();
	int nkey = golden_claw_search(16, &K1, &K2, my_rank, p);
	end_time = MPI_Wtime();
	printf("Rank %d: Time taken = %f seconds\n", my_rank,
	       end_time - start_time);

	if (my_rank == 0) {
		assert(nkey > 0);
		/* validation */
		for (int i = 0; i < nkey; i++) {
			printf("Validation step %d:\n", i);
			printf("fK1i=%lu; gK2i=%lu\n", f(K1[i]), g(K2[i]));
			assert(f(K1[i]) == g(K2[i]));
			printf("f(K1[i]) == g(K2[i]\n");
			assert(is_good_pair(K1[i], K2[i]));
			printf("Is good pair.\n");
			printf("Solution found: (%" PRIx64 ", %" PRIx64
			       ") [checked OK]\n",
			       K1[i], K2[i]);
		}
	}
	MPI_Finalize();
}
