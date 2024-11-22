#include <stdio.h>
#include <stdlib.h>

#include "homogeneous_DAG.h"

int main(int argc, char* argv[]) {
	struct Graph g;
	char* filename = "test_graph.txt";
	read_graph(&g, filename);
	printf("Number of vertices: %d\n", g.nv);
	printf("Number of edges: %d\n", g.ne);
	printf("Counts: ");
	for (int i = 0; i < g.nv; i++) {
		printf("%d ", g.counts[i]);
	}
	printf("\n");
	printf("Offsets: ");
	for (int i = 0; i < g.nv; i++) {
		printf("%d ", g.offsets[i]);
	}
	printf("\n");
	printf("Edges: ");
	for (int i = 0; i < g.ne; i++) {
		printf("%d ", g.edges[i]);
	}
	printf("\n");
	int* indegree = get_indegree(&g);
	printf("Indegree: ");
	for (int i = 0; i < g.nv; i++) {
		printf("%d ", indegree[i]);
	}
	printf("\n");
	int* result = kahn(&g);
	printf("Topological Sort: ");
	for (int i = 0; i < g.nv; i++) {
		printf("%d ", result[i]);
	}
	printf("\n");

	free(g.counts);
	free(g.offsets);
	free(g.edges);
	free(indegree);
	free(result);
	return 0;
}
