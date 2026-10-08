/*
 * gemm_transformer.c
 *
 * The memory access workload labeled GEMM-AI in the paper.
 *
 * WHAT THIS IS. It is the dominant dense matrix multiplication kernel of
 * transformer inference, expressed in plain C. Within one transformer layer,
 * six large matrix multiplications account for essentially all of the memory
 * traffic: the four projections of self-attention, namely query, key, value,
 * and output, followed by the two matrix multiplications of the feed-forward
 * network, which expand the hidden dimension and project it back.
 *
 * WHAT THIS IS NOT. It is not a transformer model. It performs no attention
 * score computation, no softmax, no normalization, and no nonlinearity. It
 * carries no trained weights and produces no meaningful output. Those omitted
 * operations contribute negligibly to memory traffic relative to the matrix
 * multiplications, which is precisely why the kernel is a defensible stand-in
 * for the memory behavior of inference and is not a stand-in for inference
 * itself. The paper describes it in exactly these terms, and any restatement
 * of this work should preserve that description.
 *
 * DIMENSIONS. Model dimension 512, feed-forward inner dimension 2048, sequence
 * length 128, one layer. These sizes were chosen so that the instrumented run
 * emits on the order of eight hundred million memory access records, producing
 * a raw trace of roughly twenty-five to thirty-two gigabytes. That places it in
 * the same envelope as the molecular dynamics and bandwidth benchmark traces,
 * which is what makes the four workloads comparable in volume.
 *
 * COMPILATION. This must be built as follows and no other way:
 *
 *     gcc -O1 -fno-tree-vectorize -o gemm_transformer gemm_transformer.c -lm
 *
 * The optimization level is not a matter of preference. At higher levels the
 * compiler vectorizes the inner loop, individual scalar memory accesses cease
 * to appear in the instrumentation output, and the resulting trace ceases to be
 * comparable with the other three workloads, all of which were traced from
 * unvectorized code.
 */

#include <stdio.h>
#include <stdlib.h>

#define D_MODEL 512   /* model, or hidden, dimension                */
#define D_FF    2048  /* feed-forward inner dimension               */
#define SEQ     128   /* sequence length, in tokens                 */

/*
 * Allocate a matrix of the given shape and fill it deterministically.
 *
 * The fill is deterministic rather than random so that two runs of this program
 * produce identical traces. The seed parameter distinguishes matrices from one
 * another so that the working set does not collapse onto one shared array,
 * which would make the access pattern artificially cache-friendly and would
 * misrepresent the memory behavior being studied.
 */
static float *make_matrix(int rows, int columns, int seed)
{
    size_t count = (size_t) rows * (size_t) columns;
    float *matrix = (float *) malloc(count * sizeof(float));
    size_t index;

    if (matrix == NULL) {
        fprintf(stderr, "allocation of %dx%d matrix failed\n", rows, columns);
        exit(1);
    }
    for (index = 0; index < count; index++) {
        matrix[index] = (float) (((index + (size_t) seed) % 97) + 1) / 97.0f;
    }
    return matrix;
}

/*
 * A deliberately naive matrix multiplication: C, of shape m by n, receives the
 * product of A, of shape m by k, and B, of shape k by n.
 *
 * The loop order is the textbook i, j, k ordering with no blocking and no
 * transposition. This is intentional. A tiled or transposed implementation
 * would convert much of the traffic into cache hits and would therefore hide
 * the very memory accesses the study exists to observe. The inner loop streams
 * one row of A and one strided column of B for each output element, and it is
 * that pattern which reaches the memory system.
 */
static void gemm(const float *A, const float *B, float *C, int m, int k, int n)
{
    int i, j, p;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            float accumulator = 0.0f;
            for (p = 0; p < k; p++) {
                accumulator += A[(size_t) i * k + p] * B[(size_t) p * n + j];
            }
            C[(size_t) i * n + j] = accumulator;
        }
    }
}

int main(void)
{
    double checksum = 0.0;
    size_t index;

    /* The activation tensor entering the layer: one sequence of token vectors. */
    float *x = make_matrix(SEQ, D_MODEL, 1);

    /* ---- Self-attention projections: query, key, value, and output ---- */
    /* Each applies a weight matrix of shape D_MODEL by D_MODEL to x. */
    float *Wq = make_matrix(D_MODEL, D_MODEL, 11);
    float *Wk = make_matrix(D_MODEL, D_MODEL, 22);
    float *Wv = make_matrix(D_MODEL, D_MODEL, 33);
    float *Wo = make_matrix(D_MODEL, D_MODEL, 44);

    float *Q = make_matrix(SEQ, D_MODEL, 0);
    float *K = make_matrix(SEQ, D_MODEL, 0);
    float *V = make_matrix(SEQ, D_MODEL, 0);
    float *O = make_matrix(SEQ, D_MODEL, 0);

    gemm(x, Wq, Q, SEQ, D_MODEL, D_MODEL);
    gemm(x, Wk, K, SEQ, D_MODEL, D_MODEL);
    gemm(x, Wv, V, SEQ, D_MODEL, D_MODEL);
    gemm(V, Wo, O, SEQ, D_MODEL, D_MODEL);   /* output projection */

    free(Wq); free(Wk); free(Wv); free(Wo);
    free(Q);  free(K);  free(V);

    /* ---- Feed-forward network: expand to D_FF, then project back ---- */
    {
        float *W1 = make_matrix(D_MODEL, D_FF, 55);
        float *W2 = make_matrix(D_FF, D_MODEL, 66);
        float *H  = make_matrix(SEQ, D_FF, 0);      /* hidden expansion */
        float *Y  = make_matrix(SEQ, D_MODEL, 0);   /* layer output     */

        gemm(O, W1, H, SEQ, D_MODEL, D_FF);   /* expand      */
        gemm(H, W2, Y, SEQ, D_FF, D_MODEL);   /* project back */

        free(W1); free(W2); free(H); free(O);
        free(x);
        x = Y;
    }

    /*
     * Accumulate and print a checksum. This exists solely so that the compiler
     * cannot prove the results unused and eliminate the multiplications, which
     * would empty the trace.
     */
    for (index = 0; index < (size_t) SEQ * D_MODEL; index += 512) {
        checksum += (double) x[index];
    }
    printf("checksum = %f\n", checksum);

    free(x);
    return 0;
}
