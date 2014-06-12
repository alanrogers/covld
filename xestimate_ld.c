#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "estimate_ld.h"

int main(void) {
    int len1, len2, len3;
    int rval, i, id;
    int X1[] = {2,0,1,1,2,0,1,1,2,0,1,1,2,0,1,1,2,0,1,1};
    int Y1[] = {2,0,1,1,2,0,1,1,2,0,1,1,2,0,1,1,2,0,1,1};
    const int X2[] = {2, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 2, 2, 1,
              2, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 2, 1, 1, 1, 1, 1,
              0, 0, 0, 1, 0, 2, 0, 1, 1, 0, 1, 1, 0, 0};
    const int Y2[] = {2, 1, 2, 2, 2, 2, 1, 2, 2, 2, 0, 2, 1, 2, 2, 2, 2, 1,
              2, 1, 2, 1, 2, 2, 1, 1, 1, 2, 2, 1, 2, 1, 1, 2, 2, 2,
              2, 2, 1, 1, 2, 2, 1, 2, 1, 2, 2, 2, 1, 1}; 
    double h[] = {0.25, 0.25, 0.25, 0.25};
    double r;
    unsigned x[3][3];
    const double tol = 1e-7;
    const unsigned max_itr = 1000;


    /* loci at positions 90 and 711 */
    /* gametic values */
    int x3[] = {1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int y3[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    /* genotypic values */
    int X3[] = {2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2,
        2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    int Y3[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

    len1 = sizeof(X1)/sizeof(X1[0]);
    len2 = sizeof(X2)/sizeof(X2[0]);
    len3 = sizeof(X3)/sizeof(X3[0]);

#if 1
    /*** test with X1, Y1 ***/
    printf("TEST 1\n");
    if( esem(len1, X1, Y1, h, tol, max_itr) < 0) {
    printf("NO CONVERGENCE with X1, Y1\n");
    exit(1);
    }
    printf("h should equal [0.5, 0, 0, 0.5]\n");
    printf("h: [%f, %f, %f, %f]\n", h[0], h[1], h[2], h[3]);

    rval = esem_r(&r, len1, X1, Y1);
    assert(rval == 0);
    printf("r (should equal 1): %f\n", r);
#endif

    /*** test with X2, Y2 ***/
    printf("\nTEST 2\n");
    count_genotypes(x, len2, X2, Y2);
    id = 0;
    for(i=0; i<len2; ++i) {
    id += 1;
    printf("Id%02d", id);
    switch(X2[i]) {
    case 0:
        printf(" %2s", "TT");
        break;
    case 1:
        printf(" %2s", "AT");
        break;
    case 2:
        printf(" %2s", "AA");
        break;
    }
    putchar(' ');
    switch(Y2[i]) {
    case 0:
        printf(" %2s", "TT");
        break;
    case 1:
        printf(" %2s", "AT");
        break;
    case 2:
        printf(" %2s", "AA");
        break;
    }
    putchar('\n');
    }
    if( esem(len2, X2, Y2, h, tol, max_itr) < 0) {
    printf("NO CONVERGENCE with X1, Y1\n");
    exit(1);
    }
    printf("h: [%f, %f, %f, %f]\n", h[0], h[1], h[2], h[3]);

    rval = esem_r(&r, len2, X2, Y2);
    assert(rval == 0);
    printf("r (should equal 0.374999919073): %f\n", r);

    /*** test with X3, Y3 ***/
    printf("\nTEST 3\n");
    r = get_r_gamete(len3, x3, y3);
    printf("r_true: %f\n", r);
    rval = esem_r(&r, len3, X3, Y3);
    assert(rval == 0);
    printf("r_ES  : %f\n", r);
    r = get_r_corr(len3, X3, Y3);
    printf("r_RH  : %f\n", r);

    /* r_true=0.040346 r_RH=0.084215 r_ES=0.040346 */

    return 0;
}
