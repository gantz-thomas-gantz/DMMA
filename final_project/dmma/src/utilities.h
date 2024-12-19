#pragma once

#include <assert.h>
#include <stdlib.h>

typedef uint64_t u64; /* portable 64-bit integer */

struct u64_darray {
	u64 *data;
	int size;
	int capacity;
};

struct u64_darray *create_u64_darray(int capacity) {
	struct u64_darray *a = malloc(sizeof(struct u64_darray));
	a->capacity = capacity;
	a->size = 0;
	a->data = (u64 *)malloc(a->capacity * sizeof(u64));
	return a;
}

void initialize_u64_darray(struct u64_darray *a, int capacity) {
	a->capacity = capacity;
	a->size = 0;
	a->data = (u64 *)malloc(a->capacity * sizeof(u64));
}

u64 read_u64_darray(struct u64_darray *a, int idx) {
	assert(idx < a->size);
	return a->data[idx];
}

void append(struct u64_darray *a, u64 val) {
	if (a->size >= a->capacity) {
		a->data = realloc(a->data, 2 * a->capacity * sizeof(u64));
		a->capacity *= 2;
	}
	a->data[a->size++] = val;
}

void clear_u64_darray(struct u64_darray *a) { a->size = 0; }

void free_u64_darray(struct u64_darray *a) { free(a->data); }
