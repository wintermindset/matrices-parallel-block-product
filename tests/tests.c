#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "../src/matrix_ops.h"

#define MATRIX_TYPE double
#define ACCURACY 1e-5

const int MATRIX_SIZE = 2880;
const int BLOCK_SIZES[] = { 1, 4, 8, 12, 16, 20, 32, 40, 48, 60, 80, 120, 160 };

void zero_init(MATRIX_TYPE** M) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            M[i][j] = 0;
        }
    }
}

void random_lower_triangular_matrix(MATRIX_TYPE** M) {
    zero_init(M);
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j <= i; ++j) {
            double range = 2000.0;
            M[i][j] = ((double)rand() / RAND_MAX) * range - 1000.0;
        }
    }
}

int main() {
    srand(time(NULL));
    clock_t begin, end;

    double Time[14] = {0};

    for (int test = 0; test < 40; ++test) {
        MATRIX_TYPE** A = malloc(sizeof(MATRIX_TYPE*) * MATRIX_SIZE);
        MATRIX_TYPE** B = malloc(sizeof(MATRIX_TYPE*) * MATRIX_SIZE);
        MATRIX_TYPE** C = malloc(sizeof(MATRIX_TYPE*) * MATRIX_SIZE);
        for (int i = 0; i < MATRIX_SIZE; ++i) {
            A[i] = malloc(sizeof(MATRIX_TYPE) * MATRIX_SIZE);
            B[i] = malloc(sizeof(MATRIX_TYPE) * MATRIX_SIZE);
            C[i] = malloc(sizeof(MATRIX_TYPE) * MATRIX_SIZE);
        }

        random_lower_triangular_matrix(A);
        random_lower_triangular_matrix(B);
        zero_init(C);

        // Naive 2D multiplication
        begin = clock();
        for (int i = 0; i < MATRIX_SIZE; ++i) {
            for (int j = 0; j <= i; ++j) {
                for (int k = 0; k <= i; ++k) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        end = clock();
        Time[13] += (double)(end - begin) / CLOCKS_PER_SEC;

        // Optimized test loop
        for (int index = 0; index < 13; ++index) {
            int BLOCK_SIZE = BLOCK_SIZES[index];
            int MAX_BLOCKS = MATRIX_SIZE / BLOCK_SIZE;
            int BLOCK_ELEMS = BLOCK_SIZE * BLOCK_SIZE;
            int NON_ZERO_BLOCKS = (MAX_BLOCKS * MAX_BLOCKS + MAX_BLOCKS) / 2;
            int TOTAL_ELEMS = NON_ZERO_BLOCKS * BLOCK_ELEMS;

            MATRIX_TYPE* A_arr = malloc(sizeof(MATRIX_TYPE) * TOTAL_ELEMS);
            MATRIX_TYPE* B_arr = malloc(sizeof(MATRIX_TYPE) * TOTAL_ELEMS);
            MATRIX_TYPE* C_arr = calloc(TOTAL_ELEMS, sizeof(MATRIX_TYPE)); // initialized to 0

            transform_to_block_format(A, A_arr, BLOCK_SIZE, MATRIX_SIZE, 1); // column-based index
            transform_to_block_format(B, B_arr, BLOCK_SIZE, MATRIX_SIZE, 0); // simple ordering

            begin = clock();
            block_multiply(A_arr, B_arr, C_arr, BLOCK_SIZE, MATRIX_SIZE);
            end = clock();
            Time[index] += (double)(end - begin) / CLOCKS_PER_SEC;

            // Validate result
            int valid = 1;
            for (int i = 0; i < MATRIX_SIZE && valid; ++i) {
                for (int j = 0; j <= i && valid; ++j) {
                    int idx = to_array_index_blocks_lines(i / BLOCK_SIZE, j / BLOCK_SIZE,
                                                          i % BLOCK_SIZE, j % BLOCK_SIZE,
                                                          BLOCK_ELEMS, BLOCK_SIZE);
                    if (fabs(C[i][j] - C_arr[idx]) >= ACCURACY) {
                        printf("Mismatch at test %d, BLOCK_SIZE = %d, index (%d, %d)\n", test, BLOCK_SIZE, i, j);
                        valid = 0;
                    }
                }
            }

            if (!valid) {
                printf("❌ Test %d failed for BLOCK_SIZE = %d\n", test, BLOCK_SIZE);
                exit(1);
            }

            free(A_arr);
            free(B_arr);
            free(C_arr);
        }

        for (int i = 0; i < MATRIX_SIZE; ++i) {
            free(A[i]); free(B[i]); free(C[i]);
        }
        free(A); free(B); free(C);
    }
	printf("✅ All 40 tests passed.\n");
    printf("Timing results (seconds):\n");
    printf("Naive:\t\t\t%.5f\n", Time[13] / 40.0);
    for (int i = 0; i < 13; ++i) {
        printf("BLOCK_SIZE = %4d:\t%.5f\n", BLOCK_SIZES[i], Time[i] / 40.0);
    }

    return 0;
}