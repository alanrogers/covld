# This program uses computer simulation to evaluate the statistical
# properties of the statistical method described in:
#
# Rogers, Alan R.  and Huff, Chad. 2008. Linkage Disequilibrium in
# Loci with Unknown Phase.
#
# I hereby place this computer program into the public domain.  Alan
# R. Rogers

from random import random, uniform
from datetime import datetime
from estimate_ld import *

#nreps = 10000  # number of repetitions
nreps = 10  # number of repetitions
ngtypes = 45 # diploid sample size

print "Time: %s" % datetime.now()

def D_to_Dprime(D, pA, pB):
    """
    Find Dprime from D, pA, and pB.  As usually defined,
    Dprime ranges btw 0 and 1.  It is modified here to range
    between -1 and 1, retaining the sign of D.
    """
    qA = 1.0 - pA
    qB = 1.0 - pB
    if D < 0:
        #a = max(-pA*pB, -qA*qB)  # usual definition
        a = -max(-pA*pB, -qA*qB)  # retain sign of D
    else:
        a = min(pA*qB, qA*pB)
    return D/a

def Dprime_to_D(Dp, pA, pB):
    """
    Find D from Dprime, pA, and pB.
    Modified to allow the sign of Dprime to equal that of D.
    """
    qA = 1.0 - pA
    qB = 1.0 - pB
    if Dp < 0:
#        a = max(-pA*pB, -qA*qB)  # for usual definition of Dp
        a = -max(-pA*pB, -qA*qB)  # allow for negative Dp
    else:
        a = min(pA*qB, qA*pB)
    return Dp*a

def bernoulli(p):
    """
    Return 1 with probability p, 0 with probability 1-p.
    """
    if random() < p:
        return 1
    return 0

def generate_gamete(pa, pb, D):
    """
    Generate pairs from following distribution:

    y  z  Prob
    -------------------------
    1  1  pa * pb         + D
    1  0  pa * (1-pb)     - D
    0  1  (1-pa) * pb     - D
    0  0  (1-pa) * (1-pb) + D

    Algorithm:

    1. Generate y = Bernoulli(pa) and z = Bernoulli(pb)
    2. Change y and z as follows:

     If D > 0:
      (1,0) --> (1,1) with prob D/(pa*(1-pb))
      (0,1) --> (0,0) with prob D/((1-pa)*pb)

     If D < 0:
      (1,1) --> (1,0) with prob -D/(pa*pb)
      (0,0) --> (0,1) with prob -D/((1-pa)*(1-pb))
    
    3. Function returns (y, z)
    """
    y = bernoulli(pa)
    z = bernoulli(pb)
    if D > 0.0:
        if y==1 and z==0:
            if random() < D/(pa*(1.0-pb)):
                z = 1
        elif y==0 and z==1:
            if random() < D/((1.0-pa)*pb):
                z = 0
    elif D < 0.0:
        if y==1 and z==1:
            if random() < -D/(pa*pb):
                z = 0
        elif y==0 and z==0:
            if random() < -D/((1.0-pa)*(1.0-pb)):
                z = 1
    return (y, z)

def header(leading_lbl, estimators):
    print

    out = "%7s" % ""
    if leading_lbl != None:
        out += " %7s" % ""
    for e in estimators:
        out += "|%21s" % e.lbl
    print out

    out = "%7s" % "N"
    if leading_lbl != None:
        out += " %7s" % leading_lbl
    for e in estimators:
        out += "|%7s %7s %5s" % ("bias", "stderr", "nconv")
    print out


# One simulated data set
def sim_step(ngtypes, epa, epb, eD, ef, estimators):
    assert ef >= 0
    v = 0.0
    while v == 0.0:
        # The i'th gamete has value (y[i], z[i])
        y = []
        z = []
        # Loop continues until we get a data set with variance
        # at both loci.
        for i in range(ngtypes):
            # Generate 1st gamete in genotype
            yval, zval = generate_gamete(epa, epb, eD)
            y.append(yval)
            z.append(zval)

            # Generate 2nd gamete in genotype
            if random() < ef:
                # 2nd gamete is identical by descent
                y.append(yval)
                z.append(zval)
            else:
                # 2nd gamete is independent
                yval, zval = generate_gamete(epa, epb, eD)
                y.append(yval)
                z.append(zval)

        assert len(y) == 2*ngtypes
        assert len(z) == 2*ngtypes

        pA, vA, pB, vB, cov = bivmom(y,z)
        v = vA*vB
    r = cov/sqrt(v) # estimated from gamete frequencies

    # The i'th genotype has value (Y[i], Z[i]).  These vectors
    # will lack information about gametic phase.
    Y = []
    Z = []
    for i in range(ngtypes):
        j = 2*i
        Y.append(y[j] + y[j+1])
        Z.append(z[j] + z[j+1])

    # Estimate r, err, and stderr from Y an Z, using all estimators
    for e in estimators:
        r_current = e.estimate(Y, Z, r)

        if e.lbl == "Excoffier-Slatkin":  #DEBUG
            print "Y:", Y
            print "Z:", Z
            print "r:", r_current

    return r

# If eps_in, epb_in, eDp_in, or ef_in are set to None, a random
# value will be chosen for each iteration.
def simulate(nreps, ngtypes, epa_in, epb_in, eDp_in, ef_in, estimators):
    # Output format
    fmt1 = "%5.2f"

    epa = epa_in
    epb = epb_in
    eDp = eDp_in
    ef = ef_in

    mean_r = 0.0
    curr_rep = 0

    for e in estimators:
        e.clear()

    while curr_rep < nreps:
        # Choose random values for unspecified parameters
        if epa_in == None:
            epa = uniform(0.05, 0.95)
        if epb_in == None:
            epb = uniform(0.05, 0.95)
        if eDp_in == None:
            eDp = uniform(-0.95, 0.95)
        if ef_in == None:
            ef = uniform(0.0, 0.9)
        eD = Dprime_to_D(eDp, epa, epb)
        
        try:
            mean_r += sim_step(ngtypes, epa, epb, eD, ef, estimators)
        except ZeroDivisionError:
            # We get here if either locus is monorphic, or if
            # there are no copies of the "11" homozygote at locus A.
            # The "continue" statement says to skip such loci.
            continue
        curr_rep += 1

    assert curr_rep == nreps
    mean_r /= float(nreps)

    out = "%7d" % ngtypes
    for val in (epa_in, epb_in, eDp_in, ef_in):
        if val != None:
            out += " %7.4f" % val
    
    for e in estimators:
        out += "|%7.4f %7.4f %5d" % (e.bias(), e.stderr(), e.n)

    print out
    return

two_n = 2*ngtypes  # haploid sample size

# Define vector of estimators
#estimators = [Estimator("Rogers-Huff", get_r), \
#              Estimator("Excoffier-Slatkin", esem_r), \
#              Estimator("Hill", Hill_r)]

#estimators = [Estimator("Rogers-Huff", get_r), \
#              Estimator("Excoffier-Slatkin", esem_r),
#              Estimator("RHES", rhesem_r)]
estimators = [Estimator("Rogers-Huff", get_r), \
              Estimator("Excoffier-Slatkin", esem_r)]
#estimators = [Estimator("Rogers-Huff", get_r)]

print "tolerance in EM algorithm:", tol
print "Replicates per simulation: %d, nconv is number that converged" \
      % nreps

if 0:
    # Loop over epa
    header("pA", estimators)
    for epa in [0.1, 0.33, 0.5,  0.67, 0.9]:
        simulate(nreps, ngtypes, epa, None, None, None, estimators)

    # No point in varying epb, since its effect is exactly the same
    # as that of epa.
    
if 0:
    # Loop over Dp
    header("Dp", estimators)
    for eDp in [-0.99, -0.9, -0.7, -0.5, 0.0, 0.5, 0.7, 0.9, 0.99]:
        simulate(nreps, ngtypes, None, None, eDp, None, estimators)

if 0:
    # Loop over f.
    header("f", estimators)
    for ef in [0.0, 0.25, 0.5, 0.75, 0.99]:
        simulate(nreps, ngtypes, None, None, None, ef, estimators)


    # Loop over ngtypes
    nvec = [25, 50, 100, 200, 400]
    header(None, estimators)
    for ngtypes in nvec:
        simulate(nreps, ngtypes, None, None, None, None, estimators)
    
if 1:
    # Mimic assumptions of coalescent simulations
    ngtypes = 50 # diploid sample size
    print "Mimic assumptions of coalescent simulations"
    header("Dp", estimators)
#    for eDp in [-0.99, -0.9, -0.7, -0.5, 0.0, 0.5, 0.7, 0.9, 0.99]:
    for eDp in [0.5, 0.7]:
        simulate(nreps, ngtypes, None, None, eDp, 0.0, estimators)
