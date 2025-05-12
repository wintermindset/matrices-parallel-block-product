#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

typedef double MATRIX_TYPE;

/*
const int MATRIX_SIZE = 2880;
const int BLOCK_SIZE = 48;
const int MAXIMUM_BLOCKS_IN_LINE = 60; // MATRIX_SIZE / BLOCK_SIZE;
const int MAXIMUM_BLOCKS_IN_COLUMN = 60; // MAXIMUM_BLOCKS_IN_LINE;
const int COUNT_NON_ZERO_BLOCKS = (60 * 60 + 60) / 2; // (MAXIMUM_BLOCKS_IN_LINE * MAXIMUM_BLOCKS_IN_LINE + MAXIMUM_BLOCKS_IN_LINE) / 2;
const int COUNT_ELEMS_IN_BLOCKS = 48 * 48; // BLOCK_SIZE * BLOCK_SIZE;
const int COUNT_ELEMS_IN_NON_ZERO_BLOCKS = 1830 * 2304; // COUNT_NON_ZERO_BLOCKS * COUNT_ELEMS_IN_BLOCKS;
const int MAXIMUM_BLOCKS_IN_COLUMN2 = 60 + 60;// MAXIMUM_BLOCKS_IN_COLUMN + MAXIMUM_BLOCKS_IN_COLUMN;
 */

void block_multiply(
    MATRIX_TYPE* A, MATRIX_TYPE* B, MATRIX_TYPE* C,
    int BLOCK_SIZE, int MATRIX_SIZE
);

void transform_to_block_format(
    MATRIX_TYPE** M, MATRIX_TYPE* M_array,
    int BLOCK_SIZE, int MATRIX_SIZE, int by_columns
);

int to_array_index_blocks_columns(int alpha, int beta, int i, int j, int MAX_BLOCKS, int BLOCK_ELEMS, int BLOCK_SIZE);
int to_array_index_blocks_lines(int alpha, int beta, int i, int j, int BLOCK_ELEMS, int BLOCK_SIZE);

#endif