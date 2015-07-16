#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static inline void print_index(int ndim, int* index) {
	printf("[");
	for (int i = 0; i <  ndim -1; i++) {
		printf("%d,", index[i]);
	}
	printf("%d]", index[ndim-1]);
}

// OK
static inline int matrix_offset(int ndim, int* index, int* strides) {
	int offset = 0;
	for (int i = 0; i < ndim; i++) {
		offset += index[i] * strides[i];
	}
	return offset;
}

static inline int add_indices(int ndim, int* index1, int* index2, int* result) {
	for (int i = 0; i < ndim; i++) {
		result[i] = index1[i] + index2[i];
	}
}

static inline double squared(double x) {
	return x * x;
}

static inline double index_size_squared(int ndim, int* index) {
	double sum = 0;
	for (int i = 0; i < ndim; i++) {
		sum += squared(index[i]);
	}
}

static inline int inbounds(int ndim, int* index, int* shape) {
	for (int i = 0; i < ndim; i++) {
		if (index[i] < 0) return 0;
		if (index[i] >= shape[i]) return 0;
	}
	return 1;
}

static inline int size(int ndim, int* index) {
	int product = 1;
	for (int i = 0; i < ndim; i++) {
		product *= index[i];
	}
	return product;
}

static inline double gamma_index_point_static(int ndim, int* shape, int* strides, int* first_index, double* first_matrix, double* second_matrix) {
	printf("Calculating ");
	print_index(ndim, first_index);
	printf("...\n");

	double dose1 = first_matrix[matrix_offset(ndim, first_index, strides)];
	double min = abs(dose1 - second_matrix[matrix_offset(ndim, first_index, strides)]);
	int max_d = (int)min;
	double min_gamma = squared(min);
	double range_d = 2 * max_d + 1;	
	int i = 0;
	int j = 0;
	int tot_attempts = 1;
	int* running_index;   // [-d..+d, ...]
	int* second_index;    // [i1-d..., ... ,...iN+d]
	int* running_strides; // [(2*d+1)^(N-1), ..., 1]
	
	// Calculate running strides
	running_strides = (int*)malloc(sizeof(int) * ndim);
	running_strides[ndim - 1] = 1;
	for (i = ndim - 2; i >= 0; i--) {
		running_strides[i] = running_strides[i+1] * range_d;
		tot_attempts *= range_d;
	}
	tot_attempts *= range_d;

	// Running index to start with (all to -max_d)
	running_index = (int*)malloc(sizeof(int) * ndim);
	for (i = 0; i < ndim; i++) {
		running_index[i] = -max_d;
	}
	
	second_index = (int*)malloc(sizeof(int) * ndim);
	for (j = 0; j < tot_attempts; j++) {
		add_indices(ndim, first_index, running_index, second_index);
		print_index(ndim, second_index);
		if (!inbounds(ndim, second_index, shape)) {
			printf("- out\n");
		} else {
			double dose2 = second_matrix[matrix_offset(ndim, second_index, strides)];
			double diff2 = squared(dose1 - dose2);
			double gamma = diff2 + index_size_squared(ndim, running_index);

			printf("= %f\n", gamma);
			if (gamma < min_gamma) {
				min_gamma = gamma;
			}
		}

		running_index[ndim - 1]++;
	 	i = ndim - 1;
		for (i = ndim - 1; i > 0; i--) {
			if (running_index[i] == max_d) {
				running_index[i] = -max_d;
				running_index[i - 1]++;
			} else {
				break;
			}
		}
	}
	
	free(running_index);
	free(running_strides);
	return min_gamma;
}

void gamma_index(int ndim, int* shape, double* first_matrix, double* second_matrix, double* result, double dd, double dta) {
	double rel_dd = dd / dta;
	double* first_matrix_rel;
	double* second_matrix_rel;
	int* strides;
	int* running_index;
	int total_size = 1;
	int i;

	// Calculate strides & total size
	strides = (int*)malloc(sizeof(int) * ndim);
	strides[ndim - 1] = 1;
	for (i = ndim - 1; i > 0; i--) {
		strides[i-1] = strides[i] * shape[i];
		total_size *= shape[i];
	}
	total_size *= shape[0];

	first_matrix_rel = malloc(sizeof(double) * total_size);
	second_matrix_rel = malloc(sizeof(double) * total_size);

	// Scale matrices
	for (i = 0; i < total_size; i++) {
		first_matrix_rel[i] = first_matrix[i] / rel_dd;
		second_matrix_rel[i] = second_matrix[i] / rel_dd;
	}

	// Running index to start with (all to -max_d)
	running_index = (int*)malloc(sizeof(int) * ndim);
	for (i = 0; i < ndim; i++) {
		running_index[i] = 0;
	}	

	// printf("%d\n", total_size);

	// Calculate gammas
	for (int j = 0; j < total_size; j++) {		
		result[j] = sqrt(gamma_index_point_static(ndim, shape, strides, running_index, first_matrix_rel, second_matrix_rel)) / dta;

		running_index[ndim - 1]++;
	 	i = ndim - 1;
		for (i = ndim - 1; i > 0; i--) {
			if (running_index[i] == shape[i]) {
				running_index[i] = 0;
				running_index[i - 1]++;
			} else {
				break;
			}
		}
	}

	free(strides);
	free(running_index);
	free(first_matrix_rel);
	free(second_matrix_rel);
}

int main() {
	double mat1[] = { 1, 2, 6, 4};
	double mat2[] = { 1, 6, 4.2, 7.5};
	double* mat3 = (double*)malloc(sizeof(double) * 4);
	int shape[] = { 2, 2};
	gamma_index(2, shape, mat1, mat2, mat3, 1, 1);

	for (int i = 0; i < 4; i++) {
		printf("%f\n", mat3[i]);
	}
}
