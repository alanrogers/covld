# This program scans the HapMap to evaluate the properties of the
# statistical method described in:
#
# Rogers, Alan R.  and Huff, Chad. 2008. Linkage Disequilibrium in
# Loci with Unknown Phase.
#
# The program examines a large number (nLoci) of focal SNPs.  From
# each focal locus, it examines all neighboring SNP within a given
# distance (window) in one direction.  All comparisons are between the
# focal SNP and a neighboring SNP.  Using each pair of SNPs, a
# randomization test is used to test the hypothesis of linkage
# equilibrium.  The program reports:
#
# nsig: the number of comparisons that were significant with at least
#       one method.
#
# RHwins: the number of significant comparisons for which p was smaller
#         for RH than for ES.
#
# ESwins: the number of significant comparisons for which p was smaller
#         for ES than for RH.
#
# A comparison was significant if fewer than 5% of randomized rsq
# values were as large as the observed one.
#
# I hereby place this computer program into the public domain.  Alan
# R. Rogers

from pgen import hapmap_dataset, hapmap_fname, strmat
from random import shuffle, randrange
from datetime import datetime
from estimate_ld import *

nreps = 500  # number of repetitions in randomization test
nLoci = 500   # number of random focal loci
chromosome = 10
pop = 'CHB'
window = 20

print "Time: %s" % datetime.now()
print "Chromosome:", chromosome
print "Pop:", pop
print "Examining %d focal loci" % nLoci
print "Replicates in randomization test: %d" % nreps

# Define vector of estimators
#estimators = [Estimator("Rogers-Huff", get_r), \
#              Estimator("Excoffier-Slatkin", esem_r), \
#              Estimator("Hill", Hill_r)]
estimators = [Estimator("RH", get_r), \
              Estimator("ES", esem_r)]
#estimators = [Estimator("Rogers-Huff", get_r)]

hds = hapmap_dataset(hapmap_fname(chromosome,pop))

rsqobs = len(estimators)*[0.0]

print "%5s" % ("dist"),
for e in estimators:
    print "%7s" % e.lbl,
print "%7s %7s %7s" % ("RHwins", "ESwins", "nsig")
nsig = RHwins = ESwins = 0
for locus in xrange(nLoci):
    i = randrange(len(hds)) # index of random snp
    # scan right
    for j in xrange(i+1, len(hds)):
        kilobases = abs(hds[j].position - hds[i].position)*0.001
        if kilobases > window:
            break
        # Estimate r using all estimators
        for k, e in enumerate(estimators):
            e.clear()
            r = e.estimate(hds[i].gtype, hds[j].gtype)
            rsqobs[k] = r*r

        # Randomization tests
        y = hds[i].gtype
        tail = len(estimators)*[0]
        for rep in xrange(nreps):
            shuffle(y)
            for k, e in enumerate(estimators):
                e.clear()
                r = e.estimate(y, hds[j].gtype)
                rsqsim = r*r
                if rsqsim >= rsqobs[k]:
                    tail[k] += 1
        pval = [x/float(nreps) for x in tail]
        if min(pval) > 0.05:
            continue
        nsig += 1
        if pval[0] == pval[1]:
            continue
        if pval[0] < pval[1]:
            RHwins += 1
        else:
            assert pval[0] > pval[1]
            ESwins += 1
        print "%5.1f" % kilobases,
        for k in range(len(tail)):
            tail[k] = tail[k]/float(nreps)
            print "%7.4f" % (tail[k]),
        print "%7d %7d %7d" % (RHwins, ESwins, nsig)
print "RHwins:", RHwins
print "ESwins:", ESwins
print "nsig:", nsig

