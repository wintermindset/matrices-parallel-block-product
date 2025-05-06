#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MATRIX_TYPE double
#define ACCURACY 1e-05

const int MATRIX_SIZE = 2880;
const int BLOCK_SIZES[] = { 1, 4, 8, 12, 16, 20, 32, 40, 48, 60, 80, 120, 160 };

double random_double_in(double min, double max) {
	double range = (max - min);
	double div = RAND_MAX / range;
	return min + (rand() / div);
}

int to_array_index_if_blocks_in_columns(int alpha, int beta, int i, int j, int MAXIMUM_BLOCKS_IN_COLUMN, int COUNT_ELEMS_IN_BLOCKS, int BLOCK_SIZE) {
	return (MAXIMUM_BLOCKS_IN_COLUMN + MAXIMUM_BLOCKS_IN_COLUMN - beta + 1) * beta / 2 * COUNT_ELEMS_IN_BLOCKS + (alpha - beta) * COUNT_ELEMS_IN_BLOCKS + i * BLOCK_SIZE + j;
}

int to_array_index_if_blocks_in_lines(int alpha, int beta, int i, int j, int COUNT_ELEMS_IN_BLOCKS, int BLOCK_SIZE) {
	return (alpha + alpha * alpha) / 2 * COUNT_ELEMS_IN_BLOCKS + beta * COUNT_ELEMS_IN_BLOCKS + i * BLOCK_SIZE + j;
}

void zero_init(MATRIX_TYPE** A) {
	for (int i = 0; i < MATRIX_SIZE; ++i) {
		for (int j = 0; j < MATRIX_SIZE; ++j) {
			A[i][j] = 0;
		}
	}
}

void random_lower_triangular_matrix(MATRIX_TYPE** A) {
	zero_init(A);
	for (int i = 0; i < MATRIX_SIZE; ++i) {
		for (int j = 0; j <= i; ++j) {
			A[i][j] = random_double_in(-1000, 1000);
		}
	}
}

int main() {
	// Setting random seed
	srand(time(NULL));

	// Defining things for time comparison
	clock_t begin, end, prog_begin, prog_end;

	prog_begin = clock();

	double Time[14];
	for (int i = 0; i < 14; ++i) {
		Time[i] = 0;
	}

	for (int test = 0; test < 40; ++test) {
		// Matrix (as two-dimensional arrays) memory-initializing
		MATRIX_TYPE** A = malloc(sizeof(MATRIX_TYPE*) * MATRIX_SIZE);
		MATRIX_TYPE** B = malloc(sizeof(MATRIX_TYPE*) * MATRIX_SIZE);
		MATRIX_TYPE** C = malloc(sizeof(MATRIX_TYPE*) * MATRIX_SIZE);

		for (int i = 0; i < MATRIX_SIZE; ++i) {
			A[i] = (MATRIX_TYPE*)malloc(sizeof(MATRIX_TYPE) * MATRIX_SIZE);
		}
		for (int i = 0; i < MATRIX_SIZE; ++i) {
			B[i] = (MATRIX_TYPE*)malloc(sizeof(MATRIX_TYPE) * MATRIX_SIZE);
		}
		for (int i = 0; i < MATRIX_SIZE; ++i) {
			C[i] = (MATRIX_TYPE*)malloc(sizeof(MATRIX_TYPE) * MATRIX_SIZE);
		}
		
		// Matrix (as two-dimensional arrays) value-initializing
		random_lower_triangular_matrix(A);
		random_lower_triangular_matrix(B);
		zero_init(C);

		// Matrix multiplication (as two-dimensional arrays) without thinking
		begin = clock();
		for (int i = 0; i < MATRIX_SIZE; ++i) {
			for (int j = 0; j < MATRIX_SIZE; ++j) {
				for (int k = 0; k < MATRIX_SIZE; ++k) {
					C[i][j] += A[i][k] * B[k][j];
				}
			}
		}
		end = clock();

		Time[13] += (double)(end - begin) / CLOCKS_PER_SEC;

		// Matrix multiplication (as one-dimensional arrays) with thinking for different BLOCK_SIZE_s
		for (int index = 0; index < 13; ++index) {
			const int BLOCK_SIZE = BLOCK_SIZES[index];
			const int MAXIMUM_BLOCKS_IN_LINE = MATRIX_SIZE / BLOCK_SIZE;
			const int MAXIMUM_BLOCKS_IN_COLUMN = MAXIMUM_BLOCKS_IN_LINE;
			const int COUNT_NON_ZERO_BLOCKS = (MAXIMUM_BLOCKS_IN_LINE * MAXIMUM_BLOCKS_IN_LINE + MAXIMUM_BLOCKS_IN_LINE) / 2;
			const int COUNT_ELEMS_IN_BLOCKS = BLOCK_SIZE * BLOCK_SIZE;
			const int COUNT_ELEMS_IN_NON_ZERO_BLOCKS = COUNT_NON_ZERO_BLOCKS * COUNT_ELEMS_IN_BLOCKS;

			// Matrix memory-initializing
			MATRIX_TYPE* A_as_array = malloc(sizeof(MATRIX_TYPE) * COUNT_ELEMS_IN_NON_ZERO_BLOCKS);
			MATRIX_TYPE* B_as_array = malloc(sizeof(MATRIX_TYPE) * COUNT_ELEMS_IN_NON_ZERO_BLOCKS);
			MATRIX_TYPE* C_as_array = malloc(sizeof(MATRIX_TYPE) * COUNT_ELEMS_IN_NON_ZERO_BLOCKS);

			// Matrix value-initializing
			int local_k = 0;
			for (int alpha = 0; alpha < MAXIMUM_BLOCKS_IN_LINE; ++alpha) {
				for (int beta = 0; beta <= alpha; ++beta) {
					for (int i = 0; i < BLOCK_SIZE; ++i) {
						for (int j = 0; j < BLOCK_SIZE; ++j) {
							A_as_array[to_array_index_if_blocks_in_columns(alpha, beta, i, j, MAXIMUM_BLOCKS_IN_COLUMN, COUNT_ELEMS_IN_BLOCKS, BLOCK_SIZE)] = A[alpha * BLOCK_SIZE + i][beta * BLOCK_SIZE + j];
							B_as_array[local_k] = B[alpha * BLOCK_SIZE + j][beta * BLOCK_SIZE + i];
							C_as_array[local_k] = 0;
							++local_k;
						}
					}
				}
			}

			// Matrix multiplication (as one-dimensional arrays) with thinking for BLOCK_SIZE = BLOCK_SIZES[index]
			begin = clock();
			int temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7;
			int alpha, alpha_alpha, beta, theta, i, j, k;
			for (alpha = 0; alpha < MAXIMUM_BLOCKS_IN_LINE; ++alpha) {
			// alpha & beta & theta - indexes of current block; i & j & k - local coords of elem in current block
				alpha_alpha = alpha * alpha;
				#pragma omp parallel for
				for (beta = 0; beta <= alpha; ++beta) {
					temp_1 = ((alpha + alpha_alpha) / 2 + beta) * COUNT_ELEMS_IN_BLOCKS;
					for (theta = beta; theta <= alpha; ++theta) {
						temp_2 = ((MAXIMUM_BLOCKS_IN_COLUMN + MAXIMUM_BLOCKS_IN_COLUMN - theta + 1) * theta / 2 + (alpha - theta)) * COUNT_ELEMS_IN_BLOCKS;
						temp_3 = ((theta + theta * theta) / 2 + beta) * COUNT_ELEMS_IN_BLOCKS;
						for (i = 0; i < BLOCK_SIZE; ++i) {
							temp_4 = i * BLOCK_SIZE;
							temp_5 = temp_1 + temp_4;
							temp_6 = temp_2 + temp_4;
							for (j = 0; j < BLOCK_SIZE; ++j) {
								temp_7 = j * BLOCK_SIZE + temp_3;
								for (k = 0; k < BLOCK_SIZE; ++k) {
									C_as_array[temp_5 + j] += A_as_array[temp_6 + k] * B_as_array[temp_7 + k];
								}
							}
						}
					}
				}
			}
			end = clock();

			// Checking result
			int flag = 1;
			for (int i = 0; i < MATRIX_SIZE; ++i) {
				for (int j = 0; j <= i; ++j) {
					if (abs(C[i][j] - C_as_array[to_array_index_if_blocks_in_lines(i / BLOCK_SIZE, j / BLOCK_SIZE,
						i % BLOCK_SIZE, j % BLOCK_SIZE, COUNT_ELEMS_IN_BLOCKS, BLOCK_SIZE)]) >= ACCURACY) {
						flag = 0;
					}
				}
			}
			
			if (flag == 0) {
				printf("ERROR, test = %d", index);
				goto end;
			}

			Time[index] += (double)(end - begin) / CLOCKS_PER_SEC;

			// Memory releasing
			free(A_as_array);
			free(B_as_array);
			free(C_as_array);
		}

		for (int i = 0; i < MATRIX_SIZE; ++i) {
			free(A[i]);
			free(B[i]);
			free(C[i]);
		}
		free(A);
		free(B);
		free(C);
	}

	prog_end = clock();

	// printf("40 tests end successfully!\nTime spent = %f\nMatrix size = %d\nResults:\n{", (double)(prog_end - prog_begin) / CLOCKS_PER_SEC, MATRIX_SIZE);
	// printf("A[i][k] * B[k][j]\t\t %f\n", Time[13] / 40.0);
	// for (int index = 0; index < 13; ++index) {
	//     printf("BLOCK_SIZE = %d\t\t%f\n", BLOCK_SIZES[index], Time[index] / 40.0);
	// }
	// printf("}");

	return 0;

	end:
	return 1;
}
