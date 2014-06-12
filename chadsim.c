#include <math.h>
#include <time.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#define NSUBJ 50

const int missing = -1;
const int linesize=500;
const long nLoci=2979351;/* number of loci in data */
const long maxdist = 500000;
const int winsize = 1600;
const int nbins = 50;
const float hourspersecond = 1.0/3600;

struct locus {
    bool empty;
    int sim;
    int pos;
    int gtype[NSUBJ];
};

typedef struct locus Locus;

/* On entry, ptr is an array of pointers to char, whose length is
   ptrlen, and str is a character string that we wish to parse.
   On return, ptr[i] points to the beginning of the i'th field
   within str, and each field is terminated with '\0'.
*/
int splitstr(int ptrlen, char **ptr, char *str) {
    int i=0, j;
    for(j=0; j < ptrlen; ++j) {
        /* find next nonwhite character */
        while(str[i]!='\0' && isspace(str[i]))
            ++i;

        if(str[i] == '\0')
            return -1;
        ptr[j] = str+i;

        if(j == ptrlen-1)
            break;
        /* find next white character */
        while(str[i]!='\0' && !isspace(str[i]))
            ++i;
        if(i == '\0')
            return -1;
        str[i] = '\0';
        ++i;
    }
    return 0;
}

void Locus_init(Locus *locus, FILE *f) {
    const int ptrlen = NSUBJ+2;
    char line[linesize];
    char *ptr[ptrlen];
    char *rval;
    int i;

    rval = fgets(line, sizeof(line), f);
    if(rval == NULL){
        locus->empty = true;
        return;
    }
    locus->empty = false;

    if (splitstr(ptrlen, ptr, line) == -1) {
        fprintf(stderr,"splitstr failed");
        exit(1);
    }
    locus->sim = strtol(ptr[0], NULL, 10);
    locus->pos = strtol(ptr[1], NULL, 10);
    for(i=0; i < NSUBJ; ++i)
        locus->gtype[i] = strtol(ptr[i+2], NULL, 10);
    return;
}

/*
 * Calculate rsq from two data vectors.
 * Skip individuals with missing data in either vector.
 */
double get_rsq(int len, int *xdata, int *ydata) {
    double mx, my, vx, vy, cov, ndbl, rsq;
    int sumx, sumy, sumxx, sumyy, sumxy;
    register int i, n=0, x, y;

    sumx = sumy = sumxx = sumyy = sumxy = 0;
    for(i=0; i<len; ++i) {
        x = xdata[i];
        y = ydata[i];
        if(x==missing || y==missing) 
            continue;
        n += 1;
        sumx += x;
        sumy += y;
        sumxx += x*x;
        sumyy += y*y;
        sumxy += x*y;
    }
    if (n>0) {
        ndbl = (double) n;
        mx = sumx/ndbl;
        my = sumy/ndbl;
        vx = sumxx/ndbl;
        vy = sumyy/ndbl;
        cov = sumxy/ndbl;
        cov -= mx * my;
        vx -= mx * mx;
        vy -= my * my;
        rsq = cov*cov/(vx*vy);
    }else
        rsq = 0.0;
    return rsq;
}
        
int main() {

    double rsq;
    int i,j, dist, bin;
    long count, ndx;
    int countvec[nbins];
    double rsqmean[nbins];
    FILE *f;
    char buff[linesize], *rval;
    Locus data[winsize];
    clock_t c0;
    double elapsed, perlocus, timeleft, fracdone;

    memset(countvec, 0, nbins*sizeof(countvec[0]));
    memset(rsqmean, 0, nbins*sizeof(rsqmean[0]));
    f = fopen("alans_data.txt", "r");
    rval = fgets(buff, sizeof(buff), f);
    if(rval==NULL) {
        fprintf(stderr,"No data");
        exit(1);
    }
    for(i=0; i < winsize; ++i) {
        Locus_init(data+i, f);
        if(data[i].empty) {
            fprintf(stderr,"Not enough data for array");
            exit(1);
        }
    }
    count = ndx = 0;
    fputc('\n', stderr);
    c0 = clock();
    while( true ) {
        if(data[ndx].empty)
            break;
        for(j=0; j < winsize; ++j) {
            /*Don't compare locus with itself*/
            if(j == ndx)
                continue;

            /* These arise when we are running out of data*/
            if(data[j].empty)
                continue;

            /* skip pairs from different simulations */
            if(data[ndx].sim != data[j].sim) 
                continue;

            dist = abs(data[j].pos - data[ndx].pos);

            /* skip pairs that are too far apart */
            if( dist > maxdist )
                continue;

            rsq = get_rsq(NSUBJ, data[ndx].gtype, data[j].gtype);
            bin = ((int) floor(nbins*dist/((float)maxdist)));
            if( bin == nbins) {
                bin = nbins-1;
            }
            countvec[bin] += 1;
            rsqmean[bin] += rsq;
        }

        Locus_init(data+ndx, f);
        ndx += 1;
        if (ndx == winsize)
            ndx = 0;
        ++count;
        if(count%1000 == 0) {
            fracdone = count/((float) nLoci);
            elapsed = ((double)(clock() - c0))/CLOCKS_PER_SEC;
            perlocus = elapsed/count;
            timeleft = (nLoci - count)*perlocus*hourspersecond;
            fprintf(stderr, "\r%0.5f%%; s/kloc=%0.5lf; %0.3f hours to go    ",
                    100*fracdone, 1000*perlocus, timeleft);
        }
    }
    fputc('\n', stderr);

    elapsed = ((double)(clock() - c0))/CLOCKS_PER_SEC;
    printf("\nElapsed time: %0.4lf sec\n", elapsed);

    printf("%8s %6s\n", "distance", "rsq");
    for (i=0; i < nbins; ++i) {
        if( countvec[i] == 0)
            continue;
        rsqmean[i] /= (float)countvec[i];
        printf("%8.0f %6.4f\n", (i+0.5)*maxdist/((float) nbins), 
               rsqmean[i]);
    }
    return 0;
}
