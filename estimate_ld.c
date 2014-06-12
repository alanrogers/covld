#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "estimate_ld.h"

/*
 * Excoffier-Slatkin EM algorithm for estimating haplotype frequencies.
 * This code implements the special case of the algorithm for two
 * biallelic loci.  With two biallelic loci, there are 4 types of
 * gamete, which I represent as follows:
 *
 *    AB  Ab  aB  ab
 *     0   1   2   3
 *
 * Here A and a are the two alleles at the first locus and B and
 * b are the alleles at the other locus.  The numbers below the
 * gamete symbols are used below as indexes into arrays.
 *
 * Phenotypes at the 1st locus: AA, Aa, and aa are numbered 0, 1, and 2.
 *
 * Phenotypes at the 2nd locus: BB, Bb, and bb are numbered 0, 1, and 2.
 *
 * Input:
 *
 * hh is a vector of 4 floats, into which results will be put
 *
 * h is a vector of 4 haplotype frequencies, indexed as shown above.
 *
 * x is a 3X3 matrix of phenotype counts.  The phenotypes are coded
 * as explained above.  Thus, x[1][2] is the number of copies of
 * the phenotype Aa/BB.
 *
 * n is the sample size and should equal the sum of x.
 *
 * Function returns a revised estimate of h after a single EM step.
 */
void esem_step(double *hh, double *h, int x[][3]) {

    /*
     * g is a 4X4 matrix of genotype frequencies. g[0][3]
     * is the frequency of the genotype that combines gamete 0 (AB)
     * with gamete 3 (ab).
     *
     * p is a 3X3 matrix of phenotype frequencies, recoded as
     *  described for the input matrix x.
     */
    double g[4][4], p[3][3], n;
    int i, j;

    /* Set only lower triangle of g.  Upper triangle unused. */
    for(i=0; i<4; ++i) {
        g[i][i] = h[i]*h[i];
        for(j=0; j<i; ++j) {
            g[i][j] = 2*h[i]*h[j];
        }
    }
        
    p[0][0] = g[0][0];
    p[0][1] = g[1][0];
    p[0][2] = g[1][1];

    p[1][0] = g[2][0];
    p[1][1] = g[3][0]+g[2][1];
    p[1][2] = g[3][1];
    
    p[2][0] = g[2][2];
    p[2][1] = g[3][2];
    p[2][2] = g[3][3];

    hh[0] = 2*x[0][0] + x[0][1] + x[1][0] + x[1][1]*g[3][0]/p[1][1]; 
    hh[1] = x[0][1] + 2*x[0][2] + x[1][1]*g[2][1]/p[1][1] + x[1][2];
    hh[2] = x[1][0] + x[1][1]*g[2][1]/p[1][1] + 2*x[2][0] + x[2][1];
    hh[3] = x[1][1]*g[3][0]/p[1][1] + x[1][2] + x[2][1] + 2*x[2][2];

    /* convert gamete counts counts to relative frequencies*/
    n = hh[0]+hh[1]+hh[2]+hh[3];
    for(i=0; i<4; ++i)
        hh[i] /= n;

    return;
}

/*
 * On return x[a][b] is the number of copies of the genotype for
 * which Y[i]=a and Z[i] = b.
 */
void count_genotypes(int x[][3], int len, const int *Y, const int *Z) {

    register int i;

    memset(&(x[0][0]), 0, 9*sizeof(x[0][0]));

    for(i=0; i<len; ++i) {
#ifndef NDEBUG
        if(!(Y[i]>=0 && Y[i]<=2))
            printf("BAD Y[%d]: %d len=%d\n", i, Y[i], len);
        if(!(Z[i]>=0 && Z[i]<=2))
            printf("BAD Z[%d]: %d len=%d\n", i, Z[i], len);
        
#endif
        assert(Y[i]>=0 && Y[i]<=2);
        assert(Z[i]>=0 && Z[i]<=2);
        x[Y[i]][Z[i]] += 1;
    }
}

/*
 * On return x[a][b] is the number of copies of the gamete for
 * which y[i]=a and z[i] = b.
 */
void count_gametes(int x[][2], int len, const int *y, const int *z) {

    register int i;

    memset(&(x[0][0]), 0, 4*sizeof(x[0][0]));

    for(i=0; i<len; ++i) {
        assert(y[i]==0 || y[i]==1);
        assert(z[i]==0 || z[i]==1);
        x[y[i]][z[i]] += 1;
    }
}

/*
 * Excoffier-Slatkin EM algorithm for estimating haplotype frequencies.
 * Input:
 *
 * Y is vector of genotype values at 1st locus, coded as 0, 1, and 2
 * to represent genotypes aa, aA, and AA.
 *
 * Z is the corresponding vector for 2nd locus.
 *
 * h holds haplotype frequencies.  Initialize these as you please.
 *   On return, h will contain the frequencies estimated by EM.
 *
 * tol controls convergence.  Algorithm stops when sum of absolute
 * differences between new and old haplotype frequencies is <= tol.
 *
 * Function returns 0 on success; -1 on failure.
 */
int esem(int len, const int *Y, const int *Z, double *h, const double tol,
         const int max_itr) {

    int i, j;
    int x[3][3];
    double dh=1.0, hh[4];

    count_genotypes(x, len, Y, Z);

    for(i=0; i<max_itr; ++i) {
        esem_step(hh, h, x);
        dh = 0.0;
        for(j=0; j<4; ++j) {
            dh += fabs(h[j]-hh[j]);
            h[j] = hh[j];
        }
        if(dh <= tol)
            break;
    }

#ifndef NDEBUG
    if(!__finite(dh)) {
        printf("NOT FINITE in esem: dh=%g\n", dh);
        for(i=0; i<3; i++)
            printf("%4d %4d %4d\n", x[i][0], x[i][1], x[i][2]);
    }
#endif  

    if( dh > tol)
        return -1;
    return 0;
}

/*
 * Use Excoffier-Slatkin EM algorithm to estimate r.
 * Input:
 *
 * Y is vector of genotype values at 1st locus, coded as 0, 1, and 2
 * to represent genotypes aa, aA, and AA.
 *
 * Z is the corresponding vector for 2nd locus.
 *
 * Function sets *r; returns 0 on success and -1 on failure.
 *
 */
int esem_r(double *result, int len, const int *Y, const int *Z) {

    /*
     * tol controls convergence.  Algorithm stops when sum of absolute
     * differences between new and old haplotype frequencies is <= tol.
     *
     * max_itr is the maximum number of EM steps.
     *
     * h is the vector of haplotype frequencies
     *
     * hh holds the haplotype frequencies estimated by esem.
     */
    const double tol = 1e-7;
    const int max_itr = 1000;
    int rval;
    
    double h[4], pA, pB, qA, qB, D, denom;

    h[0]=h[1]=h[2]=h[3]=0.25;

    rval = esem(len, Y, Z, h, tol, max_itr);

    pA = h[0]+h[1];
    pB = h[0]+h[2];
    qA = 1.0 - pA;
    qB = 1.0 - pB;
    D = h[0]*h[3] - h[1]*h[2];

    denom = pA*qA*pB*qB;
    if(denom != 0)
        *result = D/sqrt(denom);
    else {
        assert(D == 0);
        *result = 0.0;
    }

#ifndef NDEBUG
    if(!__finite(*result)) {
        printf("esem_r returning bad val\n");
        printf("h=[%g,%g,%g,%g]\n", h[0], h[1], h[2], h[3]);
    }
#endif  

    return rval;
}

/* Trick the optimizer */
void do_nothing(double x) {
}

/*
 * Calculate r from two vectors of indicator variables.
 * On input, xdata[i] = 1 if the ith gamete carries A a locus 1.  It
 * equals 0 otherwise.  ydata[i] is defined similarly for allele B at
 * the 2nd locus. 
 */
double get_r_gamete(int len, int *xdata, int *ydata) {

    register int i, iy, iz;
    double h[] = {0.0, 0.0, 0.0, 0.0};
    double pA, pB, D, denom;
    /*
     * Count gamete types.
     * h[0] counts AB
     * h[1] counts Ab
     * h[2] counts aB
     * h[3] counts ab
     */
    for(i=0; i < len; ++i) {
        iy = xdata[i];
        iz = ydata[i];
        ++h[3 - 2*iy -iz];
    }
    D = h[0] + h[1] + h[2] + h[3];
    h[0] /= D;
    h[1] /= D;
    h[2] /= D;
    h[3] /= D;
    pA = h[0] + h[1];
    pB = h[0] + h[2];
    D = h[0]*h[3] - h[1]*h[2];
    denom = pA*(1-pA)*pB*(1-pB);
    assert(denom >= 0.0);
    if(denom == 0)
        return 0.0;
    return D/sqrt(denom);
}

void print_gamete_table(FILE *f, int len, int *xdata, int *ydata) {
    int i, x[2][2];

    count_gametes(x, len, xdata, ydata);
    fprintf(f, "%3s: %4d %4d\n", "x\\y", 0, 1);
    for(i=0; i<2; ++i)
        fprintf(f, "%3i: %4d %4d\n", i, x[i][0], x[i][1]);
    return;
}

void print_genotype_table(FILE *f, int len, int *xdata, int *ydata) {
    int i, x[3][3];

    count_genotypes(x, len, xdata, ydata);
    fprintf(f, "%3s: %4d %4d %4d\n", "x\\y", 0, 1, 2);
    for(i=0; i<3; ++i)
        fprintf(f, "%3i: %4d %4d %4d\n", i, x[i][0], x[i][1], x[i][2]);
    return;
}

/*
 * Calculate r from two data vectors.
 */
double get_r_corr(int len, int *xdata, int *ydata) {
    double vx, vy, vxvy, cov, r, nsqr;
    int sumx, sumy, sumxx, sumyy, sumxy;
    register int i, n=0, x, y;

    sumx = sumy = sumxx = sumyy = sumxy = 0;
    for(i=0; i<len; ++i) {
        x = xdata[i];
        y = ydata[i];
#if 0
        if(x==MISSING || y==MISSING) 
            continue;
#endif
        n += 1;
        sumx += x;
        sumy += y;
        sumxx += x*x;
        sumyy += y*y;
        sumxy += x*y;
    }
    /* Calculate numerators in integer arithmetic to
     * avoid roundoff.
     */
    nsqr = (double) n*(n-1);
    cov = (n*sumxy - sumx*sumy)/nsqr;
    vx = (n*sumxx - sumx*sumx)/nsqr;
    vy = (n*sumyy - sumy*sumy)/nsqr;

    vxvy = vx*vy;
    if(vxvy != 0)
        r = cov/sqrt(vxvy);
    else {
        assert(cov == 0);
        r = 0.0;
    }

    return r;
}

/*
 * Calculate r from two data vectors.
 */
double get_r_corr_genotype(int len, int *xdata, int *ydata) {
    int count[3][3];
    double vx, vy, cov, r, nsqr, vxvy;
    int sumx, sumy, sumxx, sumyy, sumxy;
    register int i, j;

    if(len<2)
        return 0;

    count_genotypes(count, len, xdata, ydata);

    sumx = sumy = sumxx = sumyy = sumxy = 0;
    for(i=0; i<3; ++i) {
        for(j=0; j<3; ++j) {
            sumx += count[i][j]*i;
            sumy += count[i][j]*j;
            sumxx += count[i][j]*i*i;
            sumyy += count[i][j]*j*j;
            sumxy += count[i][j]*i*j;
        }
    }
    /* Calculate numerators in integer arithmetic to
     * avoid roundoff.
     */
    nsqr = (double) len*(len-1);
    cov = (len*sumxy - sumx*sumy)/nsqr;
    vx = (len*sumxx - sumx*sumx)/nsqr;
    vy = (len*sumyy - sumy*sumy)/nsqr;
    vxvy = vx*vy;
    if(vxvy != 0)
        r = cov/sqrt(vxvy);
    else {
        assert(cov == 0);
        r = 0.0;
    }

    return r;
}

/*
 * Calculate D from r, pA, and pB.
 */
double r_to_D(double r, double pA, double pB) {
    return r*sqrt(pA*(1.0-pA)*pB*(1.0-pB));
}
        
/*
 * Find Dprime from D, pA, and pB.
 */
double D_to_Dprime(double D, double pA, double pB) {
    double a;
    const double qA=1.0-pA, qB=1.0-pB;

    if (D < 0) {
        /*a = fmax(-pA*pB, -qA*qB); */  /* usual definition */
        a = -fmax(-pA*pB, -qA*qB);     /* retain sign of D */
    }else{
        a = fmin(pA*qB, qA*pB);
    }
    return D/a;
}

/*
 * Find Dprime from r, pA, and pB.
 */
double r_to_Dprime(double r, double pA, double pB) {
    return D_to_Dprime(r_to_D(r, pA, pB), pA, pB);
}


/*
 * Find D from Dprime, pA, and pB.
 * Modified to allow the sign of Dprime to equal that of D.
 */
double Dprime_to_D(double Dp, double pA, double pB) {
    const double qA = 1.0 - pA;
    const double qB = 1.0 - pB;
    double a;

    if (Dp < 0) {
#if 0
        a = fmax(-pA*pB, -qA*qB);  /* for usual definition of Dp */
#else
        a = -fmax(-pA*pB, -qA*qB);  /* allow for negative Dp */
#endif
    }else
        a = fmin(pA*qB, qA*pB);

    return Dp*a;
}

