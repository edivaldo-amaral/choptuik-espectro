/* Nucleo EXATO da Fase 2 (docs/COLUMNWISE_EXTERIOR.md).
 *
 * Para cada coluna (valores diadicos inteiros da coluna de Q ja arredondados
 * pelo driver Python) calcula EXATAMENTE
 *   T2 = sum_i eta_i sum_{M,N} (|re|+|im|) P1[|M|] P2[N] (2 se N>0, 1 se N=0),
 * P1[e]=p1^e q1^(E1-e), E1=m+40;  P2[N]=p2^N q2^(E2-N), E2=n+100,
 * onde (re,im) e a saida inteira (escala 2^-(72+P)) de
 *   2 S Xi Gamma2(RefA, v) [ou Gamma2 sharp] + mu R-termos,
 * na imagem tau=+ (N>0 desdobrado) e imagens sigma dadas pelo driver.
 * Aritmetica: __int128 (cotas verificadas em tempo de execucao) e GMP.
 * Entrada/saida em texto pelo stdin/stdout; o driver reconstroi a fracao.
 */
#include <gmp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef __int128 i128;
typedef unsigned __int128 u128;

typedef struct { int m1, n1; i128 a, b; } Entry;
typedef struct { int count; Entry *e; } Comp;
typedef struct { Comp comp[4]; } Key;

static int ncomp, nkeys;
static long p1, q1, p2, q2;
static long eta[4];
static i128 mu72;
static Key *keys;

static i128 parse128(const char *s) {
    int neg = 0; i128 v = 0;
    if (*s == '-') { neg = 1; s++; }
    while (*s >= '0' && *s <= '9') { v = v*10 + (*s - '0'); s++; }
    return neg ? -v : v;
}

static void mpz_set_u128(mpz_t z, u128 v) {
    mpz_set_ui(z, (unsigned long)(v >> 64));
    mpz_mul_2exp(z, z, 64);
    mpz_add_ui(z, z, (unsigned long)(v & 0xFFFFFFFFFFFFFFFFULL));
}

static i128 abs128(i128 x) { return x < 0 ? -x : x; }
#define LIMIT ((i128)1 << 120)
static int overflow = 0;
static inline i128 chk(i128 x) { if (x > LIMIT || x < -LIMIT) overflow = 1; return x; }

typedef struct {
    long id; int d, m, n, j0, sig, ncol, full;
    int keyidx[2];            /* por paridade de j */
    int nr; int rcomp[8], rpar[8]; long rmult[8];
    int *j; i128 *xr, *xi;
} Job;

static char *line = NULL; static size_t cap = 0;
static char *rd(void) { if (getline(&line, &cap, stdin) < 0) { fprintf(stderr, "EOF\n"); exit(2);} return line; }

int main(void) {
    char buf[64], sa[64], sb[64];
    /* cabecalho */
    sscanf(rd(), "SYS %d %ld %ld %ld %ld %63s", &ncomp, &p1, &q1, &p2, &q2, sa);
    mu72 = parse128(sa);
    { char *s = rd(); int off = 0, k; sscanf(s, "%63s%n", buf, &off);
      for (int i = 0; i < ncomp; i++) { sscanf(s+off, "%ld%n", &eta[i], &k); off += k; } }
    sscanf(rd(), "NKEYS %d", &nkeys);
    keys = calloc(nkeys, sizeof(Key));
    for (int kk = 0; kk < nkeys; kk++) {
        for (int i = 0; i < ncomp; i++) {
            int cnt; sscanf(rd(), "COMP %d", &cnt);
            keys[kk].comp[i].count = cnt;
            keys[kk].comp[i].e = malloc(sizeof(Entry)*(cnt ? cnt : 1));
            for (int t = 0; t < cnt; t++) {
                Entry *e = &keys[kk].comp[i].e[t];
                sscanf(rd(), "%d %d %63s %63s", &e->m1, &e->n1, sa, sb);
                e->a = parse128(sa); e->b = parse128(sb);
                if (abs128(e->a) >= ((i128)1 << 72) || abs128(e->b) >= ((i128)1 << 72)) { fprintf(stderr, "campo grande\n"); return 3; }
            }
        }
    }
    int njobs; sscanf(rd(), "NJOBS %d", &njobs);
    Job *jobs = calloc(njobs, sizeof(Job));
    for (int q = 0; q < njobs; q++) {
        Job *J = &jobs[q];
        char *s = rd(); int off = 0, k;
        sscanf(s, "JOB %ld %d %d %d %d %d %d %d %d %d%n", &J->id, &J->d, &J->m, &J->n, &J->j0, &J->sig,
               &J->full, &J->keyidx[0], &J->keyidx[1], &J->nr, &off);
        for (int r = 0; r < J->nr; r++) { sscanf(s+off, "%d %d %ld%n", &J->rcomp[r], &J->rpar[r], &J->rmult[r], &k); off += k; }
        sscanf(rd(), "NCOL %d", &J->ncol);
        J->j = malloc(sizeof(int)*J->ncol); J->xr = malloc(sizeof(i128)*J->ncol); J->xi = malloc(sizeof(i128)*J->ncol);
        for (int t = 0; t < J->ncol; t++) {
            sscanf(rd(), "%d %63s %63s", &J->j[t], sa, sb);
            J->xr[t] = parse128(sa); J->xi[t] = parse128(sb);
            if (abs128(J->xr[t]) >= ((i128)1 << 40) || abs128(J->xi[t]) >= ((i128)1 << 40)) { fprintf(stderr, "coluna grande\n"); return 4; }
        }
    }
    char **result = calloc(njobs, sizeof(char*));
    int bad = 0;
    #pragma omp parallel for schedule(dynamic,1) reduction(|:bad)
    for (int q = 0; q < njobs; q++) {
        Job *J = &jobs[q];
        int m = J->m, n = J->n;
        int rowlo = (J->sig == 2) ? -(m+40) : (m-40);
        int rows = (J->sig == 2) ? 2*m+81 : 81;
        int nofs = J->full ? n+100 : 0;           /* modo completo: N em [-(n+100), n+100] */
        int NN = J->full ? 2*n+201 : n+101;
        i128 *re = malloc(sizeof(i128)*(size_t)rows*NN), *im = malloc(sizeof(i128)*(size_t)rows*NN);
        mpz_t T, row, tmp, P1v, P2tab_e;
        mpz_inits(T, row, tmp, P1v, P2tab_e, NULL);
        int E1 = m+40, E2 = n+100;
        mpz_t *P1 = malloc(sizeof(mpz_t)*(E1+1)), *P2 = malloc(sizeof(mpz_t)*(E2+1));
        for (int e = 0; e <= E1; e++) { mpz_init(P1[e]); mpz_ui_pow_ui(tmp, p1, e); mpz_ui_pow_ui(P1[e], q1, E1-e); mpz_mul(P1[e], P1[e], tmp); }
        for (int e = 0; e <= E2; e++) { mpz_init(P2[e]); mpz_ui_pow_ui(tmp, p2, e); mpz_ui_pow_ui(P2[e], q2, E2-e); mpz_mul(P2[e], P2[e], tmp); }
        int local_over = 0;
        for (int i = 0; i < ncomp; i++) {
            memset(re, 0, sizeof(i128)*(size_t)rows*NN); memset(im, 0, sizeof(i128)*(size_t)rows*NN);
            for (int sgi = 0; sgi < J->sig; sgi++) {
                int sg = (sgi == 0) ? 1 : -1;
                for (int t = 0; t < J->ncol; t++) {
                    int j = J->j[t];
                    Comp *C = &keys[J->keyidx[j & 1]].comp[i];
                    i128 vx = J->xr[t], vy = (sg == 1) ? J->xi[t] : -J->xi[t];
                    int ntau = (J->full && j > 0) ? 2 : 1;
                    for (int ti = 0; ti < ntau; ti++) {
                        int tg = (ti == 0) ? 1 : -1;
                        for (int u = 0; u < C->count; u++) {
                            Entry *e = &C->e[u];
                            int M = sg*m + e->m1, N = tg*j + e->n1;
                            if (!J->full && N < 1) { local_over = 1; continue; }
                            if (N+nofs < 0 || N+nofs >= NN) { local_over = 1; continue; }
                            size_t idx = (size_t)(M - rowlo)*NN + (N+nofs);
                            re[idx] += e->a*vx - e->b*vy;
                            im[idx] += e->a*vy + e->b*vx;
                        }
                    }
                }
                /* termo mu R: S(N) = c_{N+1} - S(N+2), por paridade com multiplicador */
                for (int r = 0; r < J->nr; r++) {
                    if (J->rcomp[r] != i) continue;
                    int par = J->rpar[r]; i128 mult = J->rmult[r];
                    i128 *Sx = calloc(n+3, sizeof(i128)), *Sy = calloc(n+3, sizeof(i128));
                    i128 *cx = calloc(n+2, sizeof(i128)), *cy = calloc(n+2, sizeof(i128));
                    for (int t = 0; t < J->ncol; t++) if ((J->j[t] & 1) == par) {
                        cx[J->j[t]] = J->xr[t]*mult; cy[J->j[t]] = ((sg == 1) ? J->xi[t] : -J->xi[t])*mult; }
                    for (int N = n; N >= 0; N--) {
                        i128 c_x = (N+1 <= n) ? cx[N+1] : 0, c_y = (N+1 <= n) ? cy[N+1] : 0;
                        Sx[N] = c_x - Sx[N+2]; Sy[N] = c_y - Sy[N+2];
                        size_t idx = (size_t)(sg*m - rowlo)*NN + (N+nofs);
                        re[idx] += chk(2*mu72*Sx[N]); im[idx] += chk(2*mu72*Sy[N]);
                        if (J->full && N > 0) {
                            size_t idx2 = (size_t)(sg*m - rowlo)*NN + (nofs-N);
                            re[idx2] += chk(2*mu72*Sx[N]); im[idx2] += chk(2*mu72*Sy[N]);
                        }
                    }
                    free(Sx); free(Sy); free(cx); free(cy);
                }
            }
            for (int rr = 0; rr < rows; rr++) {
                int M = rr + rowlo; int aM = M < 0 ? -M : M;
                mpz_set_ui(row, 0);
                for (int Ni = 0; Ni < NN; Ni++) {
                    size_t idx = (size_t)rr*NN + Ni;
                    i128 x = re[idx], y = im[idx];
                    if (!x && !y) continue;
                    if (x > LIMIT || x < -LIMIT || y > LIMIT || y < -LIMIT) local_over = 1;
                    int N = Ni - nofs, aN = N < 0 ? -N : N;
                    u128 v = (u128)abs128(x) + (u128)abs128(y);
                    mpz_set_u128(tmp, v);
                    if (!J->full && N > 0) mpz_mul_2exp(tmp, tmp, 1);
                    mpz_addmul(row, tmp, P2[aN]);
                }
                mpz_mul(row, row, P1[aM]);
                mpz_mul_ui(row, row, (unsigned long)eta[i]);
                mpz_add(T, T, row);
            }
        }
        result[q] = mpz_get_str(NULL, 16, T);
        if (local_over || overflow) bad = 1;
        for (int e = 0; e <= E1; e++) mpz_clear(P1[e]);
        for (int e = 0; e <= E2; e++) mpz_clear(P2[e]);
        free(P1); free(P2); free(re); free(im);
        mpz_clears(T, row, tmp, P1v, P2tab_e, NULL);
    }
    if (bad) { fprintf(stderr, "estouro ou indice fora da janela\n"); return 5; }
    for (int q = 0; q < njobs; q++) printf("%ld %s\n", jobs[q].id, result[q]);
    return 0;
}
