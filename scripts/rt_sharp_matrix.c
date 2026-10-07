/* Exploratory dyadic Galerkin matrix of K(s)=K(0)+sI, not a PDE certificate. */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <gmp.h>
#include "GI.h"
#include "DyadicQ.h"
#include "Sector.h"
#include "Field.h"
#include "MultiField.h"
#include "IndexDOF.h"
#include "MultiField_approximate.h"

static unsigned long cutoff(const char *text, unsigned long maximum) {
    char *end;
    errno=0;
    unsigned long value=strtoul(text,&end,10);
    if(errno || *end || !value || value>maximum) {
        fprintf(stderr,"invalid cutoff: %s\n",text); exit(2);
    }
    return value;
}

int main(int argc, char **argv) {
    if(argc!=5) { fprintf(stderr,"usage: sharp-matrix RefA.dat M N output\n"); return 2; }
    Sector sec=Sector_new(0,0,cutoff(argv[2],250),cutoff(argv[3],750));
    FILE *input=fopen(argv[1],"r");
    if(!input) { perror(argv[1]); return 2; }
    mu_MultiField background;
    mu_MultiField_from_file_alloc(&background,input,'(',',',')');
    fclose(input);
    /* Same background truncation convention as the selected exporter. */
    Sector original=MultiField_sec(&background.mf);
    Sector clipped=Sector_new(0,0,
        2*sec.num_m<original.num_m ? 2*sec.num_m : original.num_m,
        2*sec.num_n<original.num_n ? 2*sec.num_n : original.num_n);
    MultiField_shrink(&background.mf,clipped);
    FILE *output=fopen(argv[4],"wx");
    if(!output) { perror(argv[4]); return 2; }
    IndexDOF index;
    IndexDOF_build_VBold_SHARP_alloc(&index,sec);
    unsigned long dim=IndexDOF_num(&index);
    fprintf(output,"CHOPTUIK_RT_SHARP_MATRIX_V1\ndimension %lu\n",dim);
    fprintf(output,"sector 0 0 %lu %lu\nbackground RefA.dat\nBEGIN_ADDRESSES\n",sec.num_m,sec.num_n);
    for(unsigned long i=0;i<dim;i++) {
        address a=index.__forward[i];
        fprintf(output,"%lu %lu %lu %lu %lu\n",i,a.d,a.k,a.m,a.n);
    }
    fprintf(output,"BEGIN_ENTRIES\n");
    for(unsigned long col=0;col<dim;col++) {
        MultiField basis,out;
        IndexDOF_index_to_MultiField_alloc(&basis,col,&index);
        Multifield_enlarge_attach_zeros(&basis,sec);
        MultiField_copy_alloc(&out,&basis);
        MultiField_NoENL_JBold_SHARP(&out,&background.mu);
        MultiField_NoENL_c_XiDiv_GAMMA1_SHARP_add_to(&out,&basis,&background.mu);
        /* Sharp has factor ONE, not the factor TWO of the selected pencil. */
        MultiField_NoENL_c_GAMMA2_SHARP_add_to(&out,&background.mf,&basis,1);
        for(unsigned long row=0;row<dim;row++) {
            address a=index.__forward[row];
            Field *field=out.comp+a.d;
            if(!Sector_contains_point(field->sec,a.m,a.n)) continue;
            GI *value=&field->data[a.m-field->sec.off_m][a.n-field->sec.off_n];
            mpz_ptr coefficient=a.k==_Re_ ? value->Re : value->Im;
            if(!mpz_sgn(coefficient)) continue;
            fprintf(output,"%lu %lu ",row,col);
            mpz_out_str(output,10,coefficient);
            fprintf(output," %ld\n",field->TwoExp);
        }
        MultiField_free(&out); MultiField_free(&basis);
    }
    fclose(output);
    IndexDOF_free(&index);
    mu_MultiField_free(&background);
    fprintf(stderr,"sharp matrix: %lu DOFs\n",dim);
    return 0;
}
