#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

/* AUTHOR : Charles Bouillaguet <charles.bouillaguet@lip6.fr>
   USAGE  : compile with -lm (and why not -O3)
	    redirect the standard output to a text file
	    gcc heatsink.c -O3 -lm -o heatsink
	    ./heatsink > steady_state.txt
	    then run the indicated python script for graphical rendering

   DISCLAIMER : this code does not claim to an absolute realism.
		this code could be obviously improved, but has been written so
   as to make as clear as possible the physics principle of the simulation.
*/

/* one can change the matter of the heatsink, its size, the power of the CPU,
 * etc. */
#define ALUMINIUM
#define FAST /* MEDIUM is faster, and FAST is even faster (for debugging) */
#define DUMP_STEADY_STATE

const double L = 0.15;	/* length (x) of the heatsink (m) */
const double l = 0.12;	/* height (y) of the heatsink (m) */
const double E = 0.008; /* width (z) of the heatsink (m) */
const double watercooling_T =
    20; /* temperature of the fluid for water-cooling, (°C) */
const double CPU_TDP = 280; /* power dissipated by the CPU (W) */

/* dl: "spatial step" for simulation (m) */
/* dt: "time step" for simulation (s) */
#ifdef FAST
double dl = 0.004;
double dt = 0.004;
#endif

#ifdef MEDIUM
double dl = 0.002;
double dt = 0.002;
#endif

#ifdef NORMAL
double dl = 0.001;
double dt = 0.001;
#endif

#ifdef CHALLENGE
double dl = 0.0001;
double dt = 0.00001;
#endif

/* sink_heat_capacity: specific heat capacity of the heatsink (J / kg / K) */
/* sink_density: density of the heatsink (kg / m^3) */
/* sink_conductivity: thermal conductivity of the heatsink (W / m / K) */
/* euros_per_kg: price of the matter by kilogram */
#ifdef ALUMINIUM
double sink_heat_capacity = 897;
double sink_density = 2710;
double sink_conductivity = 237;
double euros_per_kg = 1.594;
#endif

#ifdef COPPER
double sink_heat_capacity = 385;
double sink_density = 8960;
double sink_conductivity = 390;
double euros_per_kg = 5.469;
#endif

#ifdef GOLD
double sink_heat_capacity = 128;
double sink_density = 19300;
double sink_conductivity = 317;
double euros_per_kg = 47000;
#endif

#ifdef IRON
double sink_heat_capacity = 444;
double sink_density = 7860;
double sink_conductivity = 80;
double euros_per_kg = 0.083;
#endif

const double Stefan_Boltzmann =
    5.6703e-8; /* (W / m^2 / K^4), radiation of black body */
const double heat_transfer_coefficient =
    10; /* coefficient of thermal convection (W / m^2 / K) */
double CPU_surface;

/*
 * Return True if the CPU is in contact with the heatsink at the point (x,y).
 * This describes an AMD EPYC "Rome".
 */
static inline bool CPU_shape(double x, double y) {
	x -= (L - 0.0754) / 2;
	y -= (l - 0.0585) / 2;
	bool small_y_ok =
	    (y > 0.015 && y < 0.025) || (y > 0.0337 && y < 0.0437);
	bool small_x_ok =
	    (x > 0.0113 && x < 0.0186) || (x > 0.0193 && x < 0.0266) ||
	    (x > 0.0485 && x < 0.0558) || (x > 0.0566 && x < 0.0639);
	bool big_ok = (x > 0.03 && x < 0.045 && y > 0.0155 && y < 0.0435);
	return big_ok || (small_x_ok && small_y_ok);
}

/* returns the total area of the surface of contact between CPU and heatsink (in
 * m^2) */
double CPU_contact_surface() {
	double S = 0;
	for (double x = dl / 2; x < L; x += dl)
		for (double y = dl / 2; y < l; y += dl)
			if (CPU_shape(x, y)) S += dl * dl;
	return S;
}

/* Returns the new temperature of the cell (i, j, k). For this, there is an
 * access to neighbor cells (left, right, top, bottom, front, back), except if
 * (i, j, k) is on the external surface. */
static inline double update_temperature(const double *T, int u, int n, int m,
					int o, int i, int j, int k) {
	/* quantity of thermal energy that must be brought to a cell to make it
	 * heat up by 1°C */
	const double cell_heat_capacity =
	    sink_heat_capacity * sink_density * dl * dl * dl; /* J.K */
	const double dl2 = dl * dl;
	double thermal_flux = 0;

	if (i > 0)
		thermal_flux += (T[u - 1] - T[u]) * sink_conductivity *
				dl; /* neighbor x-1 */
	else {
		thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
		thermal_flux -=
		    heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
	}

	if (i < n - 1)
		thermal_flux += (T[u + 1] - T[u]) * sink_conductivity *
				dl; /* neighbor x+1 */
	else {
		thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
		thermal_flux -=
		    heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
	}

	if (j > 0)
		thermal_flux += (T[u - n] - T[u]) * sink_conductivity *
				dl; /* neighbor y-1 */
	else {
		/* Bottom cell: does it receive it from the CPU ? */
		if (CPU_shape(i * dl, k * dl))
			thermal_flux += CPU_TDP / CPU_surface * dl2;
		else {
			thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
			thermal_flux -= heat_transfer_coefficient * dl2 *
					(T[u] - watercooling_T);
		}
	}

	if (j < m - 1)
		thermal_flux += (T[u + n] - T[u]) * sink_conductivity *
				dl; /* neighbor y+1 */
	else {
		thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
		thermal_flux -=
		    heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
	}

	if (k > 0)
		thermal_flux += (T[u - n * m] - T[u]) * sink_conductivity *
				dl; /* neighbor z-1 */
	else {
		thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
		thermal_flux -=
		    heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
	}

	if (k < o - 1)
		thermal_flux += (T[u + n * m] - T[u]) * sink_conductivity *
				dl; /* neighbor z+1 */
	else {
		thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
		thermal_flux -=
		    heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
	}

	/* adjust temperature depending on the heat flux */
	return T[u] + thermal_flux * dt / cell_heat_capacity;
}

/* Run the simulation on the k-th xy plane.
 * v is the index of the start of the k-th xy plane in the arrays T and R. */
static inline void do_xy_plane(const double *T, double *R, int v, int n, int m,
			       int o, int k) {
	if (k == 0)
		// we do not modify the z = 0 plane: it is maintained at
		// constant temperature via water-cooling
		return;

	for (int j = 0; j < m; j++) {	       // y
		for (int i = 0; i < n; i++) {  // x
			int u = v + j * n + i;
			R[u] = update_temperature(T, u, n, m, o, i, j, k);
		}
	}
}

static inline void compute_cube_dims(const int p, const int n, const int m,
				     const int o, int *i, int *j, int *k) {
	double rootp;
	// Step 1: Calculate the cube root of N
	rootp = cbrt(p);  // Cube root of p
	*i = (int)rootp;
	// Step 2: Adjust i until N is divisible by i
	while (p % *i != 0) {
		*i--;
	}
	// Step 3: Find j such that N is divisible by i * j
	*j = p / *i;
	while (p % *j != 0) {
		*j--;
	}
	// Step 4: Find k such that N is divisible by i * j * k
	*k = p / ((*i) * (*j));
	// Step 5: Scaling
	*i *= n;
	*j *= m;
	*k *= o;
	return
}

struct ijk_index {
	int i;
	int j;
	int k;
};

struct ijk_index *init_rank_root_map(int o, int m, int n, int block_o,
				     int block_m, int block_n, int p) {
	struct ijk_index *rank_root_map = malloc(p * sizeof(ijk_index));
	int rank = 0;
	for (int k = 0; k < o; k = k + block_o)
		for (int j = 0; j < m; j = j + block_m)
			for (int i = 0; i < n; i = i + block_n) {
				struct ijk_index *root =
				    malloc(sizeof(struct ijk_index)) root->i =
					i;
				root->j = j;
				root->k = k;
				rank_root_map[rank] = root;
				rank++;
				if (rank + 1 > p) printf("Whoopsie.")
			}
	return rank_root_map;
}

int get_rank(struct ijk_index *rank_root_map, int p, struct ijk_index root) {
	for (int l = 0; l < p; i++) {
		if (((rank_root_map[l])->i == root->i) &&
		    ((rank_root_map[l])->j == root->j) &&
		    ((rank_root_map[l])->k == root->k))
			return l;
	}
	return -1;
}

struct ijk_index get_root(struct ijk_index *rank_root_map, int p, int rank) {
	if (rank >= p) print("Rank requested out of bounds.");
	return rank_root_map[rank];
}

struct ijk_index add_ijk(struct ijk_index a, struct ijk_index b) {
	struct ijk_index c;
	c->i = (a->i) + (b->i);
	c->j = (a->j) + (b->j);
	c->i = (a->k) + (b->k);
	return c;
}

int *init_neighbors(struct ijk_index *rank_root_map, int p, int my_rank,
		    int block_o, int block_m, int block_n, int o, int m,
		    int n) {
	int *neighbors =
	    malloc(6 * sizeof(int));  // 0 - right, 1 - left, 2 - up, 3 - down,
				      // 4 - front, 5 - back
	struct ijk_index my_root = get_root(rank_root_map, p, my_rank);
	// Get right neighbor root and rank
	struct ijk_index right = {block_n, 0, 0};
	struct ijk_index right_root = add(my_root, right);
	int right_rank = get_rank(rank_root_map, p, right_root);
	neighbors[0] = right_rank;
	// Get left neighbor root and rank
	struct ijk_index left = {-block_n, 0, 0};
	struct ijk_index left_root = add(my_root, left);
	int left_rank = get_rank(rank_root_map, p, left_root);
	neighbors[1] = left_rank;
	// Get up neighbor root and rank
	struct ijk_index up = {0, block_m, 0};
	struct ijk_index up_root = add(my_root, up);
	int up_rank = get_rank(rank_root_map, p, up_root);
	neighbors[2] = up_rank;
	// Get down neighbor root and rank
	struct ijk_index down = {0, -block_m, 0};
	struct ijk_index down_root = add(my_root, down);
	int down_rank = get_rank(rank_root_map, p, down_root);
	neighbors[3] = down_rank;
	// Get front neighbor root and rank
	struct ijk_index front = {0, 0, block_o};
	struct ijk_index front_root = add(my_root, front);
	int front_rank = get_rank(rank_root_map, p, front_root);
	neighbors[4] = front_rank;
	// Get back neighbor root and rank
	struct ijk_index back = {0, 0, -block_o};
	struct ijk_index back_root = add(my_root, back);
	int back_rank = get_rank(rank_root_map, p, back_root);
	neighbors[5] = back_rank;
}

int main() {
	MPI_Init(NULL, NULL);
	CPU_surface = CPU_contact_surface();
	double V = L * l * E;
	int n = ceil(L / dl);
	int m = ceil(E / dl);
	int o = ceil(l / dl);

	fprintf(stderr, "HEATSINK\n");
	fprintf(stderr, "\tDimension (cm) [x,y,z] = %.1f x %.1f x %.1f\n",
		100 * L, 100 * E, 100 * l);
	fprintf(stderr, "\tVolume = %.1f cm^3\n", V * 1e6);
	fprintf(stderr, "\tWeight = %.2f kg\n", V * sink_density);
	fprintf(stderr, "\tPrice = %.2f €\n", V * sink_density * euros_per_kg);
	fprintf(stderr, "SIMULATION\n");
	fprintf(stderr, "\tGrid (x,y,z) = %d x %d x %d (%.1fMo)\n", n, m, o,
		7.6293e-06 * n * m * o);
	fprintf(stderr, "\tdt = %gs\n", dt);
	fprintf(stderr, "CPU\n");
	fprintf(stderr, "\tPower = %.0fW\n", CPU_TDP);
	fprintf(stderr, "\tArea = %.1f cm^2\n", CPU_surface * 10000);

	/* Initialize MPI parameters */
	int my_rank;
	int p;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	/* temperature of each cell, in degree Kelvin. */
	if (my_rank == 0) {
		double *T = malloc(n * m * o * sizeof(*T));
		double *R = malloc(n * m * o * sizeof(*R));
		if (T == NULL || R == NULL) {
			perror("T or R could not be allocated");
			exit(1);
		}

		/* initially the heatsink is at the temperature of the
		 * water-cooling fluid */
		for (int u = 0; u < n * m * o; u++)
			R[u] = T[u] = watercooling_T + 273.15;
	}
	/* let's go! we switch the CPU on and launch the simulation until it
	 * reaches a stationary state. */
	double t = 0;
	int n_steps = 0;
	int convergence = 0;

	int *block_n;
	int *block_m;
	int *block_o;
	compute_cube_dims(p, n, m, o, block_n, block_m, block_o);
	// TODO: The domain can be extended to include the surfaces of the
	// neighboring cubes. Bear in mind the special cases for the boundaries
	// of the global domain.
	double *block_T =
	    malloc(block_n * block_m * block_o * sizeof(*block_T));
	double *block_R =
	    malloc(block_n * block_m * block_o * sizeof(*block_R));

	// Initialize map from rank to their local cube root indices
	struct ijk_index *rank_root_map =
	    init_rank_root_map(o, m, n, block_o, block_m, block_n, p);

	// Initialize neighbors
	int *neighbors = init_neighbors(rank_root_map, p, my_rank, block_o,
					block_m, block_n, o, m, n);

	/* simulating time steps */
	while (convergence == 0) {
		// Update all cells.//

		// TODO: HERE exchange halos. Bear in mind the send-receive and
		// the pattern to avoid deadlocks.
		// TODO: This loop needs to be adapted so that only the inner
		// values of the cube are cokmputed (and not the surface values
		// of the neighbors).
		for (int k = 0; k < block_o; k++) {  // z
			if (k == 0)
				// we do not modify the z = 0 plane: it is
				// maintained a t constant temperature via
				// water-cooling
				return;
			for (int j = 0; j < block_m; j++) {	     // y
				for (int i = 0; i < block_n; i++) {  // x
					int u = k * block_o * block_m +
						j * block_n + i;
					block_R[u] = update_temperature(
					    block_T, u, block_n, block_m,
					    block_o, i, j, k);
				}
			}
		}

		/* each second, we test the convergence, and print a short
		 * progress report */
		if (n_steps % ((int)(1 / dt)) == 0) {
			double block_delta_T[1] = {0};
			double block_T_max = -INFINITY;
			for (int u = 0; u < block_n * block_m * block_o; u++) {
				block_delta_T[0] += (block_R[u] - block_T[u]) *
						    (block_R[u] - block_T[u]);
				if (block_R[u] > block_T_max)
					block_T_max = block_R[u];
			}

			/* Compute global delta_T from local contributions and
			 * reduce the result to all processes so that everyone
			 * knows when to stop. */
			double delta_T[1];
			MPI_Allreduce(block_delta_T, delta_T, 1, MPI_DOUBLE,
				      MPI_SUM, MPI_COMM_WORLD);
			delta_T[0] = sqrt(delta_T[0]) / dt;
			// Compute global T_max from local max block_T_max.
			if (my_rank == 0) {
				double T_max[1];
			}
			MPI_Reduce(block_T_max, T_max, 1, MPI_DOUBLE, MPI_MAX,
				   0, MPI_COMM_WORLD);

			fprintf(
			    stderr,
			    "t = %.1fs ; T_max = %.1f°C ; convergence = %g\n",
			    t, max - 273.15, delta_T);
			if (delta_T / delta_t < 0.1) convergence = 1;
		}

		/* the new temperatures are in R */
		double *tmp = block_R;
		block_R = block_T;
		block_T = tmp;

		t += dt;
		n_steps += 1;
	}
	// TODO: Here gather the local cube temperatures block_T into the global
	// temperatures T. Bear in mind that block_T includes the extensions,
	// and this can be taken care of by using the partition of unity.

#ifdef DUMP_STEADY_STATE
	printf("###### STEADY STATE; t = %.1f\n", t);
	for (int k = 0; k < o; k++) {  // z
		printf("# z = %g\n", k * dl);
		for (int j = 0; j < m; j++) {	       // y
			for (int i = 0; i < n; i++) {  // x
				printf("%.1f ",
				       T[k * n * m + j * n + i] - 273.15);
			}
			printf("\n");
		}
	}
	printf("\n");
	fprintf(stderr,
		"For graphical rendering: python3 rendu_picture_steady.py "
		"[filename.txt] %d %d %d\n",
		n, m, o);
#endif
	MPI_Finalize();
	exit(EXIT_SUCCESS);
}
