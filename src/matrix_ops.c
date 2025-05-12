#include <omp.h>
#include "matrix_ops.h"

int to_array_index_blocks_columns(int alpha, int beta, int i, int j, int MAX_BLOCKS, int BLOCK_ELEMS, int BLOCK_SIZE) {
    return (MAX_BLOCKS + MAX_BLOCKS - beta + 1) * beta / 2 * BLOCK_ELEMS + (alpha - beta) * BLOCK_ELEMS + i * BLOCK_SIZE + j;
}

int to_array_index_blocks_lines(int alpha, int beta, int i, int j, int BLOCK_ELEMS, int BLOCK_SIZE) {
    return (alpha + alpha * alpha) / 2 * BLOCK_ELEMS + beta * BLOCK_ELEMS + i * BLOCK_SIZE + j;
}

void block_multiply(
    MATRIX_TYPE* A, MATRIX_TYPE* B, MATRIX_TYPE* C,
    int BLOCK_SIZE, int MATRIX_SIZE
) {
    int MAX_BLOCKS = MATRIX_SIZE / BLOCK_SIZE;
    int BLOCK_ELEMS = BLOCK_SIZE * BLOCK_SIZE;

    int alpha, beta, theta;
    #pragma omp parallel for private(beta, theta)
    for (alpha = 0; alpha < MAX_BLOCKS; ++alpha) {
        int temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7;
        int i, j, k;
        int alpha_sq = alpha * alpha;
        for (beta = 0; beta <= alpha; ++beta) {
            temp_1 = ((alpha + alpha_sq) / 2 + beta) * BLOCK_ELEMS;
            for (theta = beta; theta <= alpha; ++theta) {
                temp_2 = ((2 * MAX_BLOCKS - theta + 1) * theta / 2 + (alpha - theta)) * BLOCK_ELEMS;
                temp_3 = ((theta + theta * theta) / 2 + beta) * BLOCK_ELEMS;
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
}

void transform_to_block_format(
    MATRIX_TYPE** M, MATRIX_TYPE* M_array,
    int BLOCK_SIZE, int MATRIX_SIZE, int by_columns
) {
    int MAX_BLOCKS = MATRIX_SIZE / BLOCK_SIZE;
    int BLOCK_ELEMS = BLOCK_SIZE * BLOCK_SIZE;
    int k = 0;
    for (int alpha = 0; alpha < MAX_BLOCKS; ++alpha) {
        for (int beta = 0; beta <= alpha; ++beta) {
            for (int i = 0; i < BLOCK_SIZE; ++i) {
                for (int j = 0; j < BLOCK_SIZE; ++j) {
                    if (by_columns) {
                        int idx = to_array_index_blocks_columns(alpha, beta, i, j, MAX_BLOCKS, BLOCK_ELEMS, BLOCK_SIZE);
                        M_array[idx] = M[alpha * BLOCK_SIZE + i][beta * BLOCK_SIZE + j];
                    } else {
                        M_array[k++] = M[alpha * BLOCK_SIZE + j][beta * BLOCK_SIZE + i];
                    }
                }
            }
        }
    }
}