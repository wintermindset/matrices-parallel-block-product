#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MATRIX_TYPE double

const int MATRIX_SIZE = 2880;
const int BLOCK_SIZE = 48;

const int MAXIMUM_BLOCKS_IN_LINE = 60; // MATRIX_SIZE / BLOCK_SIZE;
const int MAXIMUM_BLOCKS_IN_COLUMN = 60; // MAXIMUM_BLOCKS_IN_LINE;
const int COUNT_NON_ZERO_BLOCKS = (60 * 60 + 60) / 2; // (MAXIMUM_BLOCKS_IN_LINE * MAXIMUM_BLOCKS_IN_LINE + MAXIMUM_BLOCKS_IN_LINE) / 2;
const int COUNT_ELEMS_IN_BLOCKS = 48 * 48; // BLOCK_SIZE * BLOCK_SIZE;
const int COUNT_ELEMS_IN_NON_ZERO_BLOCKS = 1830 * 2304; // COUNT_NON_ZERO_BLOCKS * COUNT_ELEMS_IN_BLOCKS;
const int MAXIMUM_BLOCKS_IN_COLUMN2 = 60 + 60;// MAXIMUM_BLOCKS_IN_COLUMN + MAXIMUM_BLOCKS_IN_COLUMN;

double random_double_in(double min, double max) {
	double range = (max - min);
	double div = RAND_MAX / range;
	return min + (rand() / div);
}

int main() {
   	MATRIX_TYPE* A = malloc(sizeof(MATRIX_TYPE) * COUNT_ELEMS_IN_NON_ZERO_BLOCKS);
	MATRIX_TYPE* B = malloc(sizeof(MATRIX_TYPE) * COUNT_ELEMS_IN_NON_ZERO_BLOCKS);
	MATRIX_TYPE* C = malloc(sizeof(MATRIX_TYPE) * COUNT_ELEMS_IN_NON_ZERO_BLOCKS);

	// Random initializing
	srand(time(NULL));
	int local_k = 0;
	for (int alpha = 0; alpha < MAXIMUM_BLOCKS_IN_LINE; ++alpha) {
		for (int beta = 0; beta <= alpha; ++beta) {
			for (int i = 0; i < BLOCK_SIZE; ++i) {
				for (int j = 0; j < BLOCK_SIZE; ++j) {
					A[local_k] = random_double_in(-1000, 1000);
					B[local_k] = random_double_in(-1000, 1000);
					C[local_k] = 0;
					++local_k;
				}
			}
		}
	}

	int temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7;
	int alpha, alpha_alpha, beta, theta, i, j, k;
	for (alpha = 0; alpha < MAXIMUM_BLOCKS_IN_LINE; ++alpha) {
		// alpha & beta & theta - indexes of current block; i & j & k - local coords of elem in current block
		alpha_alpha = alpha * alpha;
		#pragma omp parallel for
		for (beta = 0; beta <= alpha; ++beta) {
			temp_1 = ((alpha + alpha_alpha) / 2 + beta) * COUNT_ELEMS_IN_BLOCKS;
			for (theta = beta; theta <= alpha; ++theta) {
				temp_2 = ((MAXIMUM_BLOCKS_IN_COLUMN2 - theta + 1) * theta / 2 + (alpha - theta)) * COUNT_ELEMS_IN_BLOCKS;
				temp_3 = ((theta + theta * theta) / 2 + beta) * COUNT_ELEMS_IN_BLOCKS;
				for (i = 0; i < BLOCK_SIZE; ++i) {
					temp_4 = i * BLOCK_SIZE;
					temp_5 = temp_1 + temp_4;
					temp_6 = temp_2 + temp_4;
					for (j = 0; j < BLOCK_SIZE; ++j) {
						temp_7 = j * BLOCK_SIZE + temp_3;
						for (k = 0; k < BLOCK_SIZE; ++k) {
							C[temp_5 + j] += A[temp_6 + k] * B[temp_7 + k];
						}
					}
				}
			}
		}
	}

	// Memory releasing	
	free(A);
	free(B);
	free(C);

	return 0;
}