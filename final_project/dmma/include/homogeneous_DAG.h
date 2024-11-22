#include <stdbool.h>
#include <stddef.h>  // for size_t

struct Graph {
	int nv;
	int ne;
	int *counts;
	int *offsets;
	int *edges;
};

void read_graph(struct Graph *g, char *filename);

int *get_indegree(struct Graph *g);

int *kahn(struct Graph *g);
