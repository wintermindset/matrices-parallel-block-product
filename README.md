# Triangular Matrices Block Product Project

## Basic concepts

The main idea of this project was to make a product of two triangular matrices 2880x2880 **as quickly as possible** on my computer using the **C language**. For this purpose following things are used:

- `#pragma omp parallel for` for paralleling the for-loop using **OpenMP**;
- creating multiple temporary variables manually to store a result that doesn't change within the loop;
- complete exclusion of functions to skip a «`goto func` assembler part»;
- using a block product;
- storing only non-zero blocks of triangular matrices;
- storing all matrices as a one-dimensional array by block rows/columns;
- storing all the blocks of second matrix by columns for better use of cache lines;
- using **gcc** with `-Ofast` compilation flag.

The computer's characteristics also were taken into account: the best size of block for product was found during the tests (at least it depends on the size of CPU's cache lines).

40 different triples of triangular matrices (A, B, C) were randomly generated, where AB = C, A and B are lower triangular matrices (and C is also as a consequence).
Then, based on matrices A and B, their one-dimensional variants A' and B' were initialized, block multiplication A'B' = C' was performed, and then a check for C' = C was made (for every block size).
If there was a mismatch in even one element, then a `goto end` statement ends the program with an error. Otherwise there is a statistics with the average time of product for every situation (e.g. for every block size).

## Results

| Situation        | Average time (sec) |
|------------------|--------------------|
| Standard         | 118.430725         |
| BLOCK_SIZE = 1   | 41.237025          |
| BLOCK_SIZE = 4   | 2.764975           |
| BLOCK_SIZE = 8   | 1.686800           |
| BLOCK_SIZE = 12  | 1.585125           |
| BLOCK_SIZE = 16  | 1.417400           |
| BLOCK_SIZE = 20  | 1.308775           |
| BLOCK_SIZE = 32  | 1.213100           |
| BLOCK_SIZE = 40  | 1.213350           |
| BLOCK_SIZE = 48  | 1.203475           |
| BLOCK_SIZE = 60  | 1.337650           |
| BLOCK_SIZE = 80  | 1.437950           |
| BLOCK_SIZE = 120 | 1.558475           |
| BLOCK_SIZE = 160 | 1.704000           |

## Discussion

There is an extra test: **MSVC** v. 17.10 (compiler option: `/O2`) with block size = 48. It completed the task for 2.570000 sec, which is slower than **GCC** v. 13.2.0 (compiler option: `-Ofast`).
