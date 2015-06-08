/* -*- mode: C; c-basic-offset: 2; indent-tabs-mode: nil -*- */ 
/*
 * This program uses computer simulation to evaluate the statistical
 * properties of the statistical method described in:
 *
 * Rogers, Alan R.  and Huff, Chad. 2008. Linkage Disequilibrium in
 * Loci with Unknown Phase.
 *
 * Translated from tst_r.py
 */

#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "estimate_ld.h"

/* number of repetitions */
#define NREPS 100000

/* diploid sample size */
#define NGTYPES 45

/* haploid sample size */
#define TWO_N (2*NGTYPES)

#define VECLEN 100

#define ZERO(foo) memset((foo), 0, sizeof(foo));

void nonrecombinant_gamete(int yz[2], double pa, double pb, double D,
                           gsl_rng *rng);
double sim_step(double *r_covld, double *r_em,
                int ngtypes, double epa, double epb,
                double eD, double ef, double c, gsl_rng *rng);
void output(int len, const double x[], const double y[], const char *lbl);

/*
 * Generate pairs from following distribution:
 *
 * y  z  Prob
 * -------------------------
 * 1  1  pa * pb         + D
 * 1  0  pa * (1-pb)     - D
 * 0  1  (1-pa) * pb     - D
 * 0  0  (1-pa) * (1-pb) + D
 * Algorithm:
 * 1. Generate y = Bernoulli(pa) and z = Bernoulli(pb)
 * 2. Change y and z as follows:
 *   If D > 0:
 *     (1,0) --> (1,1) with prob D/(pa*(1-pb))
 *     (0,1) --> (0,0) with prob D/((1-pa)*pb)
 *   If D < 0:
 *     (1,1) --> (1,0) with prob -D/(pa*pb)
 *     (0,0) --> (0,1) with prob -D/((1-pa)*(1-pb))
 * 
 * 3. On return, yz = (y, z)
 */
void nonrecombinant_gamete(int yz[2], double pa, double pb, double D,
                           gsl_rng *rng) {
  int y = (int) gsl_ran_bernoulli(rng, pa);
  int z = (int) gsl_ran_bernoulli(rng, pb);
  if (D > 0.0) {
    if (y==1 && z==0) {
      if (gsl_rng_uniform(rng) < D/(pa*(1.0-pb))) {
        z = 1;
      }
    }else if (y==0 && z==1) {
      if (gsl_rng_uniform(rng) < D/((1.0-pa)*pb)) {
        z = 0;
      }
    }
  }else if (D < 0.0) {
    if (y==1 && z==1) {
      if (gsl_rng_uniform(rng) < -D/(pa*pb)) {
        z = 0;
      }
    }else if(y==0 && z==0) {
      if( gsl_rng_uniform(rng) < -D/((1.0-pa)*(1.0-pb))) {
        z = 1;
      }
    }
  }
  yz[0] = y;
  yz[1] = z;
  return;
}


/* One simulated data set.  Returns value of r. */
double sim_step(double *r_covld, double *r_em,
                int ngtypes, double epa, double epb,
                double eD, double ef, double c, gsl_rng *rng) {
  bool mom_recombinant, dad_recombinant;
  double r_gamete;
  int i, j, sumy, sumz, trial, max_trials = 1000;
  int y[TWO_N], z[TWO_N], yz[2], Y[NGTYPES], Z[NGTYPES];
#if 0
  int gamtab[2][2], gtypetab[3][3];
#endif
        
  assert (ef >= 0);

  /* While loop continues until we get a data set with variance at
   * both loci.
   */
  for(trial = 0; trial < max_trials; ++trial) {
    /* The i'th gamete has value (y[i], z[i]) */
    sumy = sumz = 0;
    for(i=0; i<NGTYPES; ++i) {
      mom_recombinant = (gsl_rng_uniform(rng) <= c);
      dad_recombinant = (gsl_rng_uniform(rng) <= c);
      /*
       * Make sure that if there is at least one non-recombinant,
       * we label it "mom".
       */
      if(dad_recombinant==false && mom_recombinant==true) {
        mom_recombinant = false;
        dad_recombinant = true;
      }
      if (mom_recombinant) {
        assert(dad_recombinant);
        /* Generate mom's gamete */
        yz[0] = (int) gsl_ran_bernoulli(rng, epa);
        yz[1] = (int) gsl_ran_bernoulli(rng, epb);
        y[2*i] = yz[0];
        z[2*i] = yz[1];
        sumy += yz[0];
        sumz += yz[1];
        /* Generate dad's gamete */
        if (gsl_rng_uniform(rng) >= ef) {
          /* Dad's copy of A/a is independent. */
          yz[0] = (int) gsl_ran_bernoulli(rng, epa);
        }
        if (gsl_rng_uniform(rng) >= ef) {
          /* Dad's copy of B/b is independent. */
          yz[1] = (int) gsl_ran_bernoulli(rng, epb);
        }
        y[2*i+1] = yz[0];
        z[2*i+1] = yz[1];
        sumy += yz[0];
        sumz += yz[1];
      }else{
        /* Mom non-recombinant */
        /* Generate mom's gamete */
        nonrecombinant_gamete(yz, epa, epb, eD, rng);
        y[2*i] = yz[0];
        z[2*i] = yz[1];
        sumy += yz[0];
        sumz += yz[1];
        if (dad_recombinant) {
          /*
           * Dad is recombinant but Mom isn't.  Mom's gamete may be
           * (a) IBD with Dad's A/a but not with his B/b, (b) vice
           * versa, or (c) IBD with neither.  The probabilities of
           * these outcomes are f, f, and 1-2*f.  Think of it this
           * way: Dad's gamete has two independent genes, so Mom's has
           * to lottery tickets--there are two genes in the population
           * with which Mom's might be IBD.
           */
          if (gsl_rng_uniform(rng) < ef) {
            /* Mom is IBD with Dad's A/a but not his B/b */
            yz[1] = (int) gsl_ran_bernoulli(rng, epb);
          }else if (gsl_rng_uniform(rng) < ef) {
            /* Mom is IBD with Dad's B/b but not his A/a */
            yz[0] = (int) gsl_ran_bernoulli(rng, epa);
          }else {
            /* Mom is not IBD with Dad */
            yz[0] = (int) gsl_ran_bernoulli(rng, epa);
            yz[1] = (int) gsl_ran_bernoulli(rng, epb);
          }
        }else{
          /* Neither gamete is recombinant */
          /* Generate dad's gamete */
          if (gsl_rng_uniform(rng) >= ef) {
            /* 2nd gamete is independent.  Generate a new one. */
            nonrecombinant_gamete(yz, epa, epb, eD, rng);
          }
        }
        y[2*i+1] = yz[0];
        z[2*i+1] = yz[1];
        sumy += yz[0];
        sumz += yz[1];
      }


    }
    if( sumy>0 && sumy<TWO_N && sumz>0 && sumz<TWO_N )
      break;
  }
  if(trial == max_trials) {
    fprintf(stderr,"Failed to find a polymorphic pair in %d trials\n",
            max_trials);
    exit(1);
  }
#if 0
  printf("-----------------------\n");
  printf("%4s %4s\n", "y", "z");
  for(i=0; i<TWO_N; ++i)
    printf("%4d %4d\n", y[i], z[i]);
  printf("sumy=%d sumz=%d TWO_N=%d\n", sumy, sumz, TWO_N);
  count_gametes(gamtab, TWO_N, y, z);
  printf("gamtab:\n");
  for(i=0; i<2; i++)
    printf("%4d %4d\n", gamtab[i][0], gamtab[i][1]);
#endif
        
  r_gamete = get_r_gamete(TWO_N, y, z);
  /*assert(isfinite(r_gamete));*/
  assert(__finite(r_gamete));

  /*
   * The i'th genotype has value (Y[i], Z[i]).  These vectors
   * will lack information about gametic phase.
   */
  for(i=0; i < NGTYPES; ++i) {
    j = 2*i;
    Y[i] = y[j] + y[j+1];
    Z[i] = z[j] + z[j+1];
  }

#if 0
  count_genotypes(gtypetab, NGTYPES, Y, Z);
  printf("gtypetab:\n");
  for(i=0; i<3; i++)
    printf("%4d %4d %4d\n", gtypetab[i][0], gtypetab[i][1],
           gtypetab[i][2]);
#endif   

  /* Estimate r from Y an Z, using both estimators */
  *r_covld = get_r_corr_genotype(NGTYPES, Y, Z);
  assert(__finite(*r_covld));
  if( esem_r(r_em, NGTYPES, Y, Z) ) {
    *r_em = strtod("NAN()", (char **) NULL); /* generate NAN */
  }

  return r_gamete;
}

void output(int len, const double x[], const double y[], const char *lbl) {
  int i;
  double logx, logzero = -4.0;
  printf("%% %s\n", lbl);
  printf("%% Plot of Y against log10(X), w/ log10(0) set to %f\n", logzero);
  fputs("\\plot\n", stdout);
  for(i=0; i<len; ++i) {
    logx = (x[i]==0.0 ? logzero : log10(x[i]));
    printf("  %6.3f %f\n", logx, y[i]);
  }
  fputs("/\n", stdout);
  fputs("\\multiput {$*$} at\n", stdout);
  for(i=0; i<len; ++i) {
    logx = (x[i]==0.0 ? logzero : log10(x[i]));
    printf("  %6.3f %f\n", logx, y[i]);
  }
  printf("/\n");
  return;
}

/*
 * If eps_in, epb_in, eDp_in, or ef_in are set to None, a random
 * value will be chosen for each iteration.
 */
int main(void) {

  const double c[] = {0.0, 0.001, 0.01, 0.1, 0.5};   /* recombination rate */
  double se_covld_vec[VECLEN], se_em_vec[VECLEN];
  double rsq_gam_vec[VECLEN], rsq_covld_vec[VECLEN], rsq_em_vec[VECLEN];
  double epa = 0.5;
  double epb = 0.5;
  double eDp = 0.0;
  double eD;
  double eD_children = Dprime_to_D(eDp, epa, epb);
  /*double eDp;*/
#if 1
  const double ef[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
#else
  const double ef[] = {0.0};
#endif
  double r_gamete, r_covld, r_em;
  double rsq_gamete, rsq_covld, rsq_em;
  double mean_rsq_gamete, mean_rsq_covld, mean_rsq_em;
  double var_covld, var_em, err;
  int n_em, ief, ic, curr_rep=0;
  unsigned long seed;
  int nef = sizeof(ef)/sizeof(ef[0]);
  int nc = sizeof(c)/sizeof(c[0]);
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus); /* Tausworthe generator */

  assert(nc <= VECLEN);

  /* set seed of rng */
  seed = (unsigned long) time((time_t *) 0);
  gsl_rng_set(rng, seed);
  printf("gsl: seed=%lu\n", seed);
  printf("NGTYPES: %d\n", NGTYPES);
  printf("TWO_N: %d\n", TWO_N);
  printf("NREPS: %d\n", NREPS);
  printf("allele frequencies: pA=%f pB=%f\n", epa, epb);
  printf("LD: D=%f Dp=%f\n", eD_children, eDp);

  mean_rsq_gamete = mean_rsq_covld = mean_rsq_em = 0.0;

  for(ief=0; ief < nef; ++ief) {
    ZERO(se_covld_vec);
    ZERO(se_em_vec);
    ZERO(rsq_gam_vec);
    ZERO(rsq_covld_vec);
    ZERO(rsq_em_vec);
    for(ic=0; ic < nc; ++ic) {
      /* convert D so that it refers to the value among parent */
      if(c[ic] < 1.0)
        eD = eD_children / (1.0-c[ic]);
      else
        eD = eD_children;

      var_covld = var_em = 0.0;
      n_em = 0;
      for(curr_rep=0; curr_rep < NREPS; ++curr_rep) {

        r_gamete = sim_step(&r_covld, &r_em, NGTYPES, epa, epb, eD,
                            ef[ief], c[ic], rng);
        rsq_gamete = r_gamete*r_gamete;
        rsq_covld = r_covld * r_covld;
        rsq_em = r_em * r_em;
        mean_rsq_gamete += rsq_gamete;

        err = rsq_covld - rsq_gamete;
        mean_rsq_covld += rsq_covld;
        var_covld += err*err;

        if( __finite(rsq_em) ) {
          err = rsq_em - rsq_gamete;
          var_em += err*err;
          mean_rsq_em += rsq_em;
          n_em += 1;
        }
      }
      assert( curr_rep == NREPS );
      mean_rsq_gamete /= NREPS;
      mean_rsq_covld /= NREPS;
      mean_rsq_em /= n_em;

      var_covld /= NREPS;
      var_em /= n_em;
      se_covld_vec[ic] = sqrt(var_covld);
      se_em_vec[ic] = sqrt(var_em);
      rsq_gam_vec[ic] = mean_rsq_gamete;
      rsq_covld_vec[ic] = mean_rsq_covld;
      rsq_em_vec[ic] = mean_rsq_em;
    }
    fputs("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n", stdout);
    printf("%% f=%f\n", ef[ief]);
#if 0
    output(nc, c, rsq_gam_vec, "X=c, Y=rsq_gam");
    output(nc, c, rsq_covld_vec, "X=c, Y=rsq_covld");
    output(nc, c, rsq_em_vec, "X=c, Y=rsq_em");
#endif
    output(nc, c, se_covld_vec, "X=c, Y=se_covld");
#if 0
    output(nc, c, se_em_vec, "X=c, Y=se_em");
#endif
  }
  printf("nreps=%d\n", NREPS);

  return 0;
}
