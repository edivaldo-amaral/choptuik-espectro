/* Bareiss inteiro exato; a camada Python trata os racionais com Fraction. */
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
char *det_bareiss(int n, const char **entradas) {
    mpz_t *a = malloc((size_t)n*n*sizeof(mpz_t));
    if (!a) return NULL;
    for (int i=0; i<n*n; ++i) { mpz_init(a[i]); if(mpz_set_str(a[i],entradas[i],10)) abort(); }
    mpz_t anterior, numerador;
    mpz_init_set_ui(anterior,1); mpz_init(numerador);
    int sinal=1;
    time_t inicio=time(NULL);
    for(int k=0;k<n-1;++k) {
        if(mpz_sgn(a[k*n+k])==0) {
            int r=k+1; while(r<n && mpz_sgn(a[r*n+k])==0) ++r;
            if(r==n) { mpz_set_ui(a[(n-1)*n+n-1],0); break; }
            for(int j=0;j<n;++j) mpz_swap(a[k*n+j],a[r*n+j]);
            sinal=-sinal;
        }
        for(int i=k+1;i<n;++i) {
            for(int j=k+1;j<n;++j) {
                mpz_mul(numerador,a[k*n+k],a[i*n+j]);
                mpz_submul(numerador,a[i*n+k],a[k*n+j]);
                if(!mpz_divisible_p(numerador,anterior)) { fprintf(stderr,"Divisão não exata\n"); abort(); }
                mpz_divexact(a[i*n+j],numerador,anterior);
            }
            mpz_set_ui(a[i*n+k],0);
        }
        mpz_set(anterior,a[k*n+k]);
        if((k+1)%20==0) { printf("pivô %d/%d (%lds)\n",k+1,n-1,(long)(time(NULL)-inicio)); fflush(stdout); }
    }
    if(sinal<0) mpz_neg(a[(n-1)*n+n-1],a[(n-1)*n+n-1]);
    char *saida=malloc(mpz_sizeinbase(a[(n-1)*n+n-1],10)+3);
    if(saida) mpz_get_str(saida,10,a[(n-1)*n+n-1]);
    for(int i=0;i<n*n;++i) mpz_clear(a[i]);
    free(a); mpz_clear(anterior); mpz_clear(numerador);
    return saida;
}
void liberar_determinante(char *p) { free(p); }
