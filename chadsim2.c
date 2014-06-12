#include <math.h>
#include <time.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "estimate_ld.h"

#define NSUBJ 50
#define LINESIZE 500
#define NBINS 20
#define PTRLEN (2*NSUBJ+2)
#define WINSIZE 250
#define BADDAT (EOF-1)
#define MAXKB 70  /* max dist btw pairs of SNPs in kb */
#define PER_KB 0.001
#define SHUFFLE 1

#if 0                
const long nLoci=5000;   /* number of loci in dummy data */
#define FNAME "chaddat-short.txt"
#else
const long nLoci=2979350;/* number of loci in data */ 
#define FNAME "chaddat.txt"
#endif
const float hourspersecond = 1.0/3600;

typedef struct {
    bool empty;
    int sim, pos;
    double p;
    int gamete[2*NSUBJ];
    int gtype[NSUBJ];
} Locus;

typedef struct {
    unsigned long long ntot;
    double minx, maxx, range;
    unsigned long long count[NBINS];
    double y[NBINS];
    double truey[NBINS];
    double mse[NBINS];
    const char *xname;
    const char *yname;
} ResTbl;

int Locus_init(Locus *locus, FILE *f, int *rndx);
void Locus_print(Locus *locus, FILE *f);
ResTbl *ResTbl_new(double minx, double maxx, const char *xname_str,
                   const char *yname_str);
void ResTbl_tabulate(ResTbl *tbl, double x, double y, double truval);
void ResTbl_print(ResTbl *tbl, FILE *f);
int skipline(FILE *f);

/* allocate and initialize a table */
ResTbl *ResTbl_new(double minx, double maxx, const char *xname_str,
                   const char *yname_str) {
    assert(minx < maxx);
    ResTbl *tbl = malloc(sizeof(ResTbl));
    if(tbl==NULL)
        return NULL;
    tbl->ntot = 0;
    tbl->minx = minx;
    tbl->maxx = maxx;
    tbl->range = maxx - minx;
    memset(tbl->count, 0, sizeof(tbl->count));
    memset(tbl->y, 0, sizeof(tbl->y));
    memset(tbl->truey, 0, sizeof(tbl->truey));
    memset(tbl->mse, 0, sizeof(tbl->mse));
    tbl->xname = strdup(xname_str);
    tbl->yname = strdup(yname_str);
    return tbl;
}

/* add a value to the table */
void ResTbl_tabulate(ResTbl *tbl, double x, double y, double truval) {
    double err;
    int bin;

    /* find the bin corresponding to x */
    bin = ((int) floor(NBINS*(x-tbl->minx)/tbl->range));
    if( bin >= NBINS ) {
        bin = NBINS-1;
    }
    if( bin < 0)
        bin = 0;

    tbl->ntot += 1;
    tbl->count[bin] += 1;
    tbl->y[bin] += y;
    tbl->truey[bin] += truval;
    err = y - truval;
    tbl->mse[bin] += err*err;
    return;
}

/* Print a table of results */
void ResTbl_print(ResTbl *tbl, FILE *f) {

    int j;
    double y, truey, standard_error, n;
    double my=0.0, mtruey=0.0, mse=0.0;

    fprintf(f, "%13s %10s %10s %6s %8s\n",
            tbl->xname, tbl->yname, "true", "stderr", "N");
    for (j=0; j < NBINS; ++j) {
        if( tbl->count[j] == 0)
            continue;
        n = tbl->count[j];
        y = tbl->y[j] / n;
        truey = tbl->truey[j]/n;
        standard_error = sqrt(tbl->mse[j]/n);
        fprintf(f,"%13.4f %10.7f %10.7f %6.4f %8llu\n",
                tbl->minx + (j+0.5)*tbl->range/((double) NBINS),
                y, truey, standard_error, tbl->count[j]);
        my += tbl->y[j];
        mtruey += tbl->truey[j];
        mse += tbl->mse[j];
    }
    my /= tbl->ntot;
    mtruey /= tbl->ntot;
    mse /= tbl->ntot;
    fprintf(f,"%13s %10.7f %10.7f %6.4f %8llu\n", "Mean", 
            my, mtruey, sqrt(mse), tbl->ntot);
    return;
}


/*
 * Read from stream until either a newline or EOF is reached.
 * Return 0 in the former case and EOF in the latter.
 */
int skipline(FILE *f) {
    int i;

    do{
        i = getc(f);
    }while(i != '\n' && i != EOF);
    if(i==EOF)
        return EOF;
    return 0;
}

/*
 * Initialize an object of type Locus.
 * On entry:
 *
 * locus   points to the object to be initialized
 * f       points to the FILE from which the next line
 *         of data should be read.  The locus is initialized
 *         based on the values in this line of data.
 * rndx    is an array of ints [0,1,...,2*NSUBJ]
 *         in random order.  The genic values in the data
 *         are rearranged in the order specified by rndx.
 *         Then, adjacent genes are combined to form genotypes.
 *         This simulates random mating.  ALL loci are rearranged
 *         in parallel in order to maintain the integrity of
 *         chromosomes.
 *
 * On return:
 *
 *   On success, all components of Locus are set, Locus.empty=false,
 *   and function returns 0.  There are two kinds of failure:
 *
 *   On EOF, Locus.empty=true and function returns EOF.
 *   If bad data are detected, Locus.empty=true and function returns
 *   BADDAT.  It also eats the rest of the offending line of input.
 */
int Locus_init(Locus *locus, FILE *f, int *rndx) {
    int i, rval;
    double p;

    rval = fscanf(f, "%d%d", &(locus->sim), &(locus->pos));
    if(rval!= 2) {
        locus->empty = true;
        if(skipline(f) == EOF)
            return EOF;
        return BADDAT;
    }
    for(i=0; i < 2*NSUBJ; ++i) {
        rval = fscanf(f, "%d", locus->gamete+rndx[i]);
        if(rval!= 1) {
            locus->empty = true;
            if(skipline(f) == EOF)
                return EOF;
            return BADDAT;
        }
    }

    locus->empty = false;

    p = 0.0;
    for(i=0; i < NSUBJ; ++i) {
        locus->gtype[i] = locus->gamete[2*i];
        locus->gtype[i] += locus->gamete[2*i+1];
        p += locus->gtype[i];
    }
    locus->p = p/(2*NSUBJ);

    return 0;
}

void Locus_print(Locus *locus, FILE *f) {
    int i;
    fprintf(f, "Locus.empty : %s\n", (locus->empty?"true":"false"));
    fprintf(f, "Locus.sim   : %d\n", locus->sim);
    fprintf(f, "Locus.pos   : %d\n", locus->pos);
    fprintf(f, "Locus.gamete:");
    for(i=0; i<2*NSUBJ; ++i)
        fprintf(f," %d", locus->gamete[i]);
    putc('\n', f);
    fprintf(f, "Locus.gtype :");
    for(i=0; i<NSUBJ; ++i)
        fprintf(f," %d", locus->gtype[i]);
    putc('\n', f);
    return;
}


int main(void) {
    const bool verbose = false;
    const double maxkb = MAXKB;
    int i,j, esem_rval;
    long count, ndx;
    double dist;
    FILE *f;
    Locus data[WINSIZE];
    time_t  t0, t1;
    double elapsed, perlocus;
    double r, rsq_RH, rsq_ES, rsq_true, Dp_true, r_av;
    int rndx[2*NSUBJ];
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
    ResTbl *tbl_dst_rh = ResTbl_new(0.0, maxkb, "distance",
                                    "rsq_RH");
    ResTbl *tbl_dst_es = ResTbl_new(0.0, maxkb, "distance",
                                    "rsq_ES"); 
    ResTbl *tbl_dst_rh_nocnv = ResTbl_new(0.0, maxkb,
										  "distance", "rsq_RH_nocnv");
	ResTbl *tbl_dst_avrsq = ResTbl_new(0.0, maxkb, "distance", "RHEMavrsq");
	ResTbl *tbl_dst_avr = ResTbl_new(0.0, maxkb, "distance", "RHEMavr");
    ResTbl *tbl_Dp_rh = ResTbl_new(-1.0, 1.0, "Dp", "rsq_RH");
    ResTbl *tbl_Dp_es = ResTbl_new(-1.0, 1.0, "Dp", "rsq_ES");

	t0 = time(NULL);
	printf("Date and time: %s", ctime(&t0));
	fflush(stdout);

    gsl_rng_set(rng, (unsigned) t0);

    /*
     * rndx is a shuffled list of consecutive integers.
     * It is used to randomize the order of genes w/i each locus.
     */
    for(i=0; i<2*NSUBJ; ++i)
        rndx[i] = i;
    if(SHUFFLE) {
		printf("Shuffling gametes\n");
		gsl_ran_shuffle(rng, rndx, 2*NSUBJ, sizeof(rndx[0]));
    }else 
		printf("Did NOT shuffle gametes\n");
    
    f = fopen(FNAME, "r");
    if(f==NULL) {
        fprintf(stderr,"Can't open file '%s'.\n", FNAME);
        exit(1);
    }
    printf("Reading data in file %s\n", FNAME);

    /*
     * Populate data array before starting main loop.
     */
    for(i=0; i < WINSIZE; ++i) {
        do{
            j=Locus_init(data+i, f, rndx);
        }while(j == BADDAT);
        if(j==EOF) {
            fprintf(stderr,"Can't initialize data array");
            exit(1);
        }
    }
    fputc('\n', stderr);
    count = ndx = 0;
    t0 = time(NULL);

    /* Loop over focal loci. */
    while( true ) {
        if(data[ndx].empty)
            break;

        /* Compare focal locus with other loci in window. */
        for(j=0; j < WINSIZE; ++j) {
            /*Don't compare locus with itself*/
            if(j == ndx)
                continue;

            /* These arise when we are running out of data*/
            if(data[j].empty)
                continue;

            /* skip pairs from different simulations */
            if(data[ndx].sim != data[j].sim) 
                continue;

            dist = fabs((double)(data[j].pos - data[ndx].pos))*PER_KB;

            /* skip pairs that are too far apart */
            if( dist > maxkb )
                continue;

  
            /* true values */
			r = get_r_gamete(2*NSUBJ, data[ndx].gamete, data[j].gamete);
			rsq_true = r*r;
            Dp_true = r_to_Dprime(r, data[ndx].p, data[j].p);
			assert(fabs(Dp_true) <= 1.0001);

            /* RH method */
			r_av = r = get_r_corr(NSUBJ, data[ndx].gtype, data[j].gtype);
			/*r = get_r_corr_genotype(NSUBJ, data[ndx].gtype,
			  data[j].gtype);*/
			rsq_RH = r*r;
            ResTbl_tabulate(tbl_dst_rh, dist, rsq_RH, rsq_true);
            ResTbl_tabulate(tbl_Dp_rh, Dp_true, rsq_RH, rsq_true);

            /* ES method */
			esem_rval = esem_r(&r, NSUBJ, data[ndx].gtype, data[j].gtype);
			rsq_ES = r*r;
            if( esem_rval ==0) {
                /* tabulate only if EM converges */
                ResTbl_tabulate(tbl_dst_es, dist, rsq_ES, rsq_true);
                ResTbl_tabulate(tbl_Dp_es, Dp_true, rsq_ES, rsq_true);
				/* Average of rsq and of r */
				r_av = 0.5*(r + r_av);
				ResTbl_tabulate(tbl_dst_avrsq, dist, 0.5*(rsq_RH+rsq_ES),
								rsq_true);
				ResTbl_tabulate(tbl_dst_avr, dist, r_av*r_av,
								rsq_true);
            }else
                ResTbl_tabulate(tbl_dst_rh_nocnv, dist, rsq_RH, rsq_true);


			if(verbose) {
				printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
				printf("dist: %f kb\n", dist);
				printf("Locus %5d gametes:", data[ndx].pos);
				for(i=0; i<2*NSUBJ; ++i)
					printf(" %d", data[ndx].gamete[i]);
				putchar('\n');
				printf("Locus %5d gametes:", data[j].pos);
				for(i=0; i<2*NSUBJ; ++i)
					printf(" %d", data[j].gamete[i]);
				putchar('\n');

				printf("Locus %5d genotypes:", data[ndx].pos);
				for(i=0; i<NSUBJ; ++i)
					printf(" %d", data[ndx].gtype[i]);
				putchar('\n');
				printf("Locus %5d genotypes:", data[j].pos);
				for(i=0; i<NSUBJ; ++i)
					printf(" %d", data[j].gtype[i]);
				putchar('\n');
				printf("Gamete table:\n");
				print_gamete_table(stdout, 2*NSUBJ, data[ndx].gamete,
								   data[j].gamete);
				printf("Genotype table:\n");
				print_genotype_table(stdout, NSUBJ, data[ndx].gtype,
									 data[j].gtype);

				if(esem_rval == 0)
					printf("esem_r converged\n");
				else
					printf("esem_r did not converge\n");
				printf("rsq_true=%f rsq_RH=%f rsq_ES=%f\n",
					   rsq_true, rsq_RH, rsq_ES);
				printf("<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
			}
        }

        /* read another locus */
        while( Locus_init(data+ndx, f, rndx) == BADDAT)
			;
        ++count;
        ndx += 1;
        if (ndx == WINSIZE)
            ndx = 0;
        if(count%5000 == 0) {
            t1 = time(NULL);
            elapsed = t1 - t0;
            perlocus = elapsed/count;
            fprintf(stderr,
                    "count=%ld; s/kloc=%0.5f; sim=%d; pos=%d; eof=%d err=%d\n",
                    count, 1000*perlocus, data[ndx].sim, data[ndx].pos,
                    feof(f), ferror(f));
        }
    }
    fputc('\n', stderr);

    t1 = time(NULL);
    printf("\nElapsed time: %0.4f hours\n", (t1-t0)*hourspersecond);

    printf("Distance versus RH:\n");
    ResTbl_print(tbl_dst_rh, stdout);
    putchar('\n');

    printf("Distance versus ES:\n");
    ResTbl_print(tbl_dst_es, stdout);
    putchar('\n');

    printf("Distance versus average of rsq_RH and rsq_ES:\n");
    ResTbl_print(tbl_dst_avrsq, stdout);
    putchar('\n');

    printf("Distance versus squared average of r_RH and r_ES:\n");
    ResTbl_print(tbl_dst_avr, stdout);
    putchar('\n');

    printf("Distance versus RH when ES didn't converge:\n");
    ResTbl_print(tbl_dst_rh_nocnv, stdout);
    putchar('\n');

    printf("Dp versus RH:\n");
    ResTbl_print(tbl_Dp_rh, stdout);
    putchar('\n');

    printf("Dp versus ES:\n");
    ResTbl_print(tbl_Dp_es, stdout);
    putchar('\n');

    return 0;
}
