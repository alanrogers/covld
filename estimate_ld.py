# This module defines several estimators of linkage disequilibrium
# statistics D and r.  All estimators work with data that lack gametic
# phase.  Each estimator requires two vectors, Y and Z, as input.  The
# i'th entry in Y is 2, 1, or 0, for genotypes AA, Aa, and aa, at the
# first locus.  The entries of Z are defined similarly for genotypes of
# the second locus.
#
# The Rogers-Huff method is described in:
#
# Rogers, Alan R.  and Huff, Chad. 2008. Linkage Disequilibrium in
# Loci with Unknown Phase.
#
# I hereby place this computer program into the public domain.  Alan
# R. Rogers

from math import sqrt
import traceback
from myexcept import *

# tol controls convergence.
tol = 1e-7

def bivmom(vec0, vec1):
    """
    Calculate means, variances, the covariance, from two data vectors.
    On entry, vec0 and vec1 should be vectors of numeric values and
    should have the same length.  Function returns m0, v0, m1, v1,
    cov, where m0 and m1 are the means of vec0 and vec1, v0 and v1 are
    the variances, and cov is the covariance.
    """
    m0 = m1 = v0 = v1 = cov = 0
    for x, y in zip(vec0, vec1):
        m0 += x
        m1 += y
        v0 += x*x
        v1 += y*y
        cov += x*y
    n = len(vec0)
    assert n == len(vec1)
    n = float(n)
    m0 /= n
    m1 /= n
    v0 /= n
    v1 /= n
    cov /= n

    cov -= m0 * m1
    v0 -= m0 * m0
    v1 -= m1 * m1
    return m0, v0, m1, v1, cov

def get_covD(Y,Z):
    """
    get_covD estimates pA, pB, and D w/o info on gametic phase.
    Uses the method of Rogers and Huff 2008.
    """
    pA, v0, pB, v1, cov = bivmom(Y,Z)
    pA = 0.5*pA
    pB = 0.5*pB
    qA=1-pA
    qB=1-pB
    two_1pf = sqrt((v0*v1)/(pA*qA*pB*qB)) # estimates 2(1+f)
    D = cov/two_1pf
    return pA, pB, D

def get_r(Y,Z):
    """
    Estimates r w/o info on gametic phase.  Also works with gametic
    data, in which case Y and Z should be vectors of 0/1 indicator
    variables.
    Uses the method of Rogers and Huff 2008.
    """
    mY, vY, mZ, vZ, cov = bivmom(Y,Z)
    return cov/sqrt(vY*vZ)

def get_r_corr_genotype(Y, Z):
    print "Enter get_r_corr_genotype"
    count = [[0,0,0],[0,0,0],[0,0,0]]

    print "A"
    print "type(Y):", type(Y)
    print "len(Y):", len(Y)

    for i in range(len(Y)):
        print i
        print Y[i]
        print Z[i]
        count[Y[i]][Z[i]] += 1

    print "B",

    sumx = sumy = sumxx = sumyy = sumxy = 0
    for i in range(3):
        for j in range(3):
			sumx += count[i][j]*i
			sumy += count[i][j]*j
			sumxx += count[i][j]*i*i
			sumyy += count[i][j]*j*j
			sumxy += count[i][j]*i*j

    print "sumx=%d sumy=%d sumxx=%d sumyy=%d sumxy=%d" \
        % (sumx, sumy, sumxx, sumyy, sumxy)
    # Calculate numerators in integer arithmetic to
    # avoid roundoff.
    len = len(Y)
    nsqr = float(len*(len-1))
    cov = (len*sumxy - sumx*sumy)/nsqr
    vx = (len*sumxx - sumx*sumx)/nsqr
    vy = (len*sumyy - sumy*sumy)/nsqr
    r = cov/sqrt(vx*vy)
    print "get_r_corr_genotype returning:", r
    return r

# Excoffier-Slatkin EM algorithm for estimating haplotype frequencies.
# This code implements the special case of the algorithm for two
# biallelic loci.  With two biallelic loci, there are 4 types of
# gamete, which I represent as follows:
#
#    AB  Ab  aB  ab
#     0   1   2   3
#
# Here A and a are the two alleles at the first locus and B and
# b are the alleles at the other locus.  The numbers below the
# gamete symbols are used below as indexes into arrays.
#
# Phenotypes at the 1st locus: AA, Aa, and aa are numbered 0, 1, and 2.
#
# Phenotypes at the 2nd locus: BB, Bb, and bb are numbered 0, 1, and 2.
#
# Input:
#
# h is a vector of 4 haplotype frequencies, indexed as shown above.
#
# x is a 3X3 matrix of phenotype counts.  The phenotypes are coded
# as explained above.  Thus, x[1][2] is the number of copies of
# the phenotype Aa/BB.
#
# n is the sample size and should equal the sum of x.
#
# Function returns a revised estimate of h after a single EM step.
def esem_step(h, x):

    # g is a 4X4 matrix of genotype frequencies. g[0][3]
    # is the frequency of the genotype that combines gamete 0 (AB)
    # with gamete 3 (ab).
    g = [[None,None,None,None],[None,None,None,None],
         [None,None,None,None],[None,None,None,None]]
    for i in range(4):
        g[i][i] = h[i]*h[i]
        for j in range(i):
            g[i][j] = 2*h[i]*h[j]

    # p is a 3X3 matrix of phenotype frequencies, recoded as
    # described for the input matrix x.
    p = [[None,None,None],[None,None,None],[None,None,None]]

    p[0][0] = g[0][0]
    p[0][1] = g[1][0]
    p[0][2] = g[1][1]

    p[1][0] = g[2][0]
    p[1][1] = g[3][0]+g[2][1]
    p[1][2] = g[3][1]

    p[2][0] = g[2][2]
    p[2][1] = g[3][2]
    p[2][2] = g[3][3]

    hh = [None, None, None, None]
    hh[0] = 2*x[0][0] + x[0][1] + x[1][0] + x[1][1]*g[3][0]/p[1][1]
    hh[1] = x[0][1] + 2*x[0][2] + x[1][1]*g[2][1]/p[1][1] + x[1][2]
    hh[2] = x[1][0] + x[1][1]*g[2][1]/p[1][1] + 2*x[2][0] + x[2][1]
    hh[3] = x[1][1]*g[3][0]/p[1][1] + x[1][2] + x[2][1] + 2*x[2][2]

    # haploid sample size
    n = float(sum(hh))

    # convert gamete counts counts to relative frequencies
    for i in range(4):
        hh[i] /= n

    return hh

# This is exactly like esem_step.  It's here so that I can count
# calls.
def rhesem_step(h, x):

    # g is a 4X4 matrix of genotype frequencies. g[0][3]
    # is the frequency of the genotype that combines gamete 0 (AB)
    # with gamete 3 (ab).
    g = [[None,None,None,None],[None,None,None,None],
         [None,None,None,None],[None,None,None,None]]
    for i in range(4):
        g[i][i] = h[i]*h[i]
        for j in range(i):
            g[i][j] = 2*h[i]*h[j]

    # p is a 3X3 matrix of phenotype frequencies, recoded as
    # described for the input matrix x.
    p = [[None,None,None],[None,None,None],[None,None,None]]

    p[0][0] = g[0][0]
    p[0][1] = g[1][0]
    p[0][2] = g[1][1]

    p[1][0] = g[2][0]
    p[1][1] = g[3][0]+g[2][1]
    p[1][2] = g[3][1]

    p[2][0] = g[2][2]
    p[2][1] = g[3][2]
    p[2][2] = g[3][3]

    hh = [None, None, None, None]
    hh[0] = 2*x[0][0] + x[0][1] + x[1][0] + x[1][1]*g[3][0]/p[1][1]
    hh[1] = x[0][1] + 2*x[0][2] + x[1][1]*g[2][1]/p[1][1] + x[1][2]
    hh[2] = x[1][0] + x[1][1]*g[2][1]/p[1][1] + 2*x[2][0] + x[2][1]
    hh[3] = x[1][1]*g[3][0]/p[1][1] + x[1][2] + x[2][1] + 2*x[2][2]

    # haploid sample size
    n = float(sum(hh))

    # convert gamete counts counts to relative frequencies
    for i in range(4):
        hh[i] /= n

    return hh


# Excoffier-Slatkin EM algorithm for estimating haplotype frequencies.
# Input:
#
# Y is vector of genotype values at 1st locus, coded as 0, 1, and 2
# to represent genotypes aa, aA, and AA.
#
# Z is the corresponding vector for 2nd locus.
#
# h is the initial vector of haplotype frequencies
#
# Function returns h, a vector of 4 haplotype frequencies.
def esem(Y,Z, h=[0.25,0.25,0.25,0.25], max_itr=1000):
    global tol
    x = [[0,0,0],[0,0,0],[0,0,0]]
    for y,z in zip(Y,Z):
        x[y][z] += 1
    for itr in xrange(max_itr):
        hh = esem_step(h, x)
        dh = 0.0
        for u,v in zip(h, hh):
            dh += abs(u-v)
        if dh <= tol:
            break
        h = hh
    if dh > tol:
        raise ConvergenceError
    return hh

# Exactly like esem.  Here to count calls
def rhesem(Y,Z, h=[0.25,0.25,0.25,0.25], max_itr=1000):
    global tol
    x = [[0,0,0],[0,0,0],[0,0,0]]
    for y,z in zip(Y,Z):
        x[y][z] += 1
    for itr in xrange(max_itr):
        hh = rhesem_step(h, x)
        dh = 0.0
        for u,v in zip(h, hh):
            dh += abs(u-v)
        if dh <= tol:
            break
        h = hh
    if dh > tol:
        raise ConvergenceError
    return hh

# Use Rogers-Huff method to obtain initial values for the
# Excoffier-Slatkin algorithm.  Return r.
#
# Input:
#
# Y is vector of genotype values at 1st locus, coded as 0, 1, and 2
# to represent genotypes aa, aA, and AA.
#
# Z is the corresponding vector for 2nd locus.
def rhesem_r(Y,Z, max_itr=1000):
    global tol
    # RH step
    pA, pB, D = get_covD(Y,Z)
    qA=1.0-pA
    qB=1.0-pB
    h = [pA*pB + D, pA*qB-D, qA*pB-D, qA*qB+D]

    # ES step
    h = rhesem(Y,Z, h, max_itr)
    pA = h[0]+h[1]
    pB = h[0]+h[2]
    qA = 1.0 - pA
    qB = 1.0 - pB
    r = h[0]*h[3] - h[1]*h[2]
    r /= sqrt(pA*qA*pB*qB)
    return r

# Use Excoffier-Slatkin EM algorithm to estimate r.
# Input:
#
# Y is vector of genotype values at 1st locus, coded as 0, 1, and 2
# to represent genotypes aa, aA, and AA.
#
# Z is the corresponding vector for 2nd locus.
#
# h is the initial vector of haplotype frequencies
#
# Function returns r.
def esem_r(Y,Z, h=[0.25,0.25,0.25,0.25], max_itr=1000):
    global tol
    h = esem(Y,Z, h, max_itr)
    pA = h[0]+h[1]
    pB = h[0]+h[2]
    qA = 1.0 - pA
    qB = 1.0 - pB
    r = h[0]*h[3] - h[1]*h[2]
    r /= sqrt(pA*qA*pB*qB)
    return r

# Algorithm of Hill 1974
def HillD(Y, Z, initD=None):
    global tol
    max_itr = 1000
    sampsize = len(Y)
    assert(sampsize > 0)
    n00 = n01 = n02 = 0.0
    n10 = n11 = n12 = 0.0
    n20 = n21 = n22 = 0.0
    pA = pB = 0.0
    for i in xrange(sampsize):
        pA += Y[i]
        pB += Z[i]
        score = 10*Y[i] + Z[i]
        if score==00:
            n00 +=1
        elif score==01:
            n01 += 1
        elif score==10:
            n10 +=1
        elif score==11:
            n11 += 1
        else:
            pass
    pA /= 2.0*sampsize
    pB /= 2.0*sampsize

    if initD:
        x = pA*pB + initD
    else:
        x = (2.0*n00 + n01 + n10)/(2.0*(sampsize - n11))
    itr = 0
    while 1:
        y = n11*x*(1.0-pA-pB-x)
        y /= x*(1-pA-pB-x) + (pA-x)*(pB-x)
        y += 2*n00 + n01 + n10
        y /= 2*sampsize
        if abs(x-y) <= tol  or itr > 1000:
            break
        x = y
        itr += 1
    if itr >= max_itr:
        raise ConvergenceError
    D = x - pA*pB
    return pA, pB, D

def Hill_r(Y, Z):
    pA, pB, D = HillD(Y,Z)
    qA = 1-pA
    qB = 1-pB
    return D/sqrt(pA*qA*pB*qB)

class Estimator:
    """
    This class defines a generic estimator.  Initialize like this:

    e = Estimator('mylabel', myestimator)

    where 'mylabel' is the name of your estimator and myestimator is
    a function that takes two data vectors returns some value (presumably
    an estimate of something).

    Thereafter, e.lbl returns the label, e.estimate(Y,X,truval) returns the
    estimate of your statistic obtained from data vectors Y and X, and stores
    the err and MSE.  e.bias() returns the mean error, and stderr returns
    the standard error (root mean squared error).
    """
    def __init__(self, lbl, estimator):
        self.lbl = lbl
        self.estimator = estimator
        self.clear()
        return

    def clear(self):
        """
        Set all numeric values to zero.
        """
        self.n = 0
        self.sum_err = 0.0
        self.sum_mse = 0.0
        return

    def estimate(self, Y, Z, truval=None):
        """
        Estimate r and record error and squared error.

        On entry: Y and Z are data vectors of genotypic values,
        coded as 0, 1, and 2, where 1 is the heterozygote and 0,2
        are the two homozygotes.  truval is the true value of
        r, obtained in some other way.

        On return: the state of the object has been modified to
        reflect another observation of error and squared error, and the
        estimated value is returned.

        If the estimator raises a ConvergenceError or a
        ZeroDivisionError, the state of the object is unchanged and
        the function returns None.
        """
        try:
            val = self.estimator(Y,Z)
        except (ConvergenceError, ZeroDivisionError):
            # These exceptions mean that the method
            # failed with this data set.  This will
            # show up in self.n
            pass
        except Exception, detail:
            print '-'*60
            traceback.print_exc()
            print '-'*60
            exit(1)
        else:
            if truval:
                err = val - truval
                self.sum_err += err
                self.sum_mse += err*err
            self.n += 1
            return val
        return None

    def bias(self):
        if self.n == 0:
            return None
        return self.sum_err/float(self.n)

    def stderr(self):
        if self.n == 0:
            return None
        mse = self.sum_mse/float(self.n)
        return sqrt(mse)

if __name__ == '__main__':

    estimators = [Estimator("Rogers-Huff", get_r), \
                  Estimator("RH-tabulate", get_r_corr_genotype), \
                  Estimator("Excoffier-Slatkin", esem_r), \
                  Estimator("Hill", Hill_r),
                  Estimator("RHES", rhesem_r)]

    if 1:
        print "Test 1"
        Y = [2,0,1,1,2,0,1,1,2,0,1,1,2,0,1,1,2,0,1,1]
        Z = [2,0,1,1,2,0,1,1,2,0,1,1,2,0,1,1,2,0,1,1]
        r_expect = 1.0
        print "All estimates below should equal %f" % r_expect
        for e in estimators:
            r = e.estimate(Y, Z, 1)
            print "Estimator %s returned %f" % (e.lbl, r)

        print "\nTest 2"
    Y = [2, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 2, 2, 1,
         2, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 2, 1, 1, 1, 1, 1,
         0, 0, 0, 1, 0, 2, 0, 1, 1, 0, 1, 1, 0, 0]

    Z = [2, 1, 2, 2, 2, 2, 1, 2, 2, 2, 0, 2, 1, 2, 2, 2, 2, 1,
         2, 1, 2, 1, 2, 2, 1, 1, 1, 2, 2, 1, 2, 1, 1, 2, 2, 2,
         2, 2, 1, 1, 2, 2, 1, 2, 1, 2, 2, 2, 1, 1]

    r_expect = 0.374999919073

    r = esem_r(Y, Z)
    print "esem_r returned %f; should be %f" % (r, r_expect)
