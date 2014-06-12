#ifndef ARR_ESTIMATE_LD_H
#define ARR_ESTIMATE_LD_H
#define MISSING (-1)
void print_genotype_table(FILE *f, int len, int *xdata, int *ydata);
void print_gamete_table(FILE *f, int len, int *xdata, int *ydata);
void count_genotypes(int x[][3], int len, const int *Y, const int *Z);
void count_gametes(int x[][2], int len, const int *y, const int *z);
void esem_step(double *hh, double *h, int x[][3]);
int esem(int len, const int *Y, const int *Z, double *h, const double tol,
		 const int max_itr);
int esem_r(double *result, int len, const int *Y, const int *Z);
void do_nothing(double x);
double get_r_corr(int len, int *xdata, int *ydata);
double get_r_corr_genotype(int len, int *xdata, int *ydata);
double get_r_gamete(int len, int *xdata, int *ydata);
double r_to_Dprime(double r, double pA, double pB);
double Dprime_to_D(double Dp, double pA, double pB);
double D_to_Dprime(double D, double pA, double pB);
double r_to_D(double r, double pA, double pB);
#endif
