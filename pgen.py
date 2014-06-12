#--*--python--*--
# Python stuff for Anth/Biol 5221: Human Evolutionary Genetics
from random import random, randrange
from math import exp, log, tan, sqrt, ceil, floor, pi, log10
import bisect
from bisect import bisect_left  
import re
import os

hapmap_dir = None

# Figure out where the HapMap directory is on the current machine.
def get_dirpath():
    dirnames = ["/home/rogers/hapmap/hapmap-r23",
                "/Volumes/anthropology/hapmap-r23",
                "."]
    for dn in dirnames:
        if os.path.exists(dn):
            hapmap_dir = dn
            break
    return dn

class Tabulator:
    """
    Tabulate numerical data into bins.

    Alan Rogers 2007
    """
    def __init__(self, low=0.0, high=1.0, nBins=5):
        """
        Tabulator(...) -> new Tabulator

        low and high specify the range of the data to be tabulated.

        nBins specifies number of bins
        """
        self.posinf = high+1.0  # Python 2.4 lacks support for inf
        self.neginf = low-1.0   # Python 2.4 lacks support for inf
        self.low = float(low)
        self.high = float(high)
        self.binsize = (self.high - self.low)/nBins
        self.nBins = nBins+2     # for values outside range
        self.counts = self.nBins*[0]

        # lbound = [-inf, ..., high] size = nBins
        # ubound = [low, ..., inf] size=nBins
        self.ubound = [self.low + i*self.binsize \
                       for i in range(0, self.nBins)]
        self.ubound[self.nBins-2] = self.high
        self.ubound[self.nBins-1] = self.posinf

    def __len__(self):
        "Return the number of bins"
        return len(self.counts)

    def clear(self):
        """
        Set all counts to zero.
        """
        self.counts = self.nBins*[0]

    def merge(self, other):
        """
        Within each bin, add the counts in "other" to those in "self".
        On return, the counts in "self" are the sums of those
        originally in "self" and "other".  Both Tabulators must have
        the same values of low, high, and nBins.
        """
        if self.low != other.low \
                or self.high != other.high \
                or self.nBins != other.nBins:
            raise ValueError, "can't add dissimilar Tabulators"
        for i in range(self.nBins):
            self.counts[i] += other.counts[i]
        return

    def bin(self, x):
        """
        return index of bin containing x
    
        bin 0           : values < low
        bins 1..nBins-2 : values in [low,high]
        bin nBins-1     : values > high
        """

        # if statements avoid a bug in Python 2.4
        if x > self.high:
            return self.nBins-1
        if x < self.low:
            return 0
        if x == self.low:
            return 1

        i = bisect.bisect_left(self.ubound, x)
        assert i >= 0
        assert i < self.nBins
        return i

    def __iadd__(self, x):
        '"self += x" increments count in bin containing x'

        i = self.bin(x)
        self.counts[i] += 1
        return self

    def __getitem__(self, i):
        'self[i] returns the number of items in bin i'
        assert i >= 0
        assert i < self.nBins
        return self.counts[i]

    def upperBound(self, i):
        "Upper bound of i'th bin"
        return self.ubound[i]

    def lowerBound(self, i):
        "Lower bound of i'th bin"
        assert i >= 0
        assert i < self.nBins
        if i == 0:
            return self.neginf
        return self.ubound[i-1]

    def sampleSize(self):
        'Sample size'
        return sum(self.counts)

    def __str__(self):
        'Convert to string'
        sampsize = float(sum(self.counts))
        s = "%-24s %6s %8s\n" % ("Range", "Count", "Freq")
        s += "[%9s ...%9.4f] %6d %8.4f\n" % ("-inf",
                                         self.upperBound(0),
                                         self.counts[0],
                                         self.counts[0]/sampsize)        
        for i in range(1, self.nBins-1):
            s += "[%9.4f ...%9.4f] %6d %8.4f\n" % (self.lowerBound(i),
                                          self.upperBound(i),
                                          self.counts[i],
                                          self.counts[i]/sampsize)
        s += "[%9.4f ...%9s] %6d %8.4f\n" % (self.lowerBound(self.nBins-1),
                                       "inf",
                                       self.counts[self.nBins-1],
                                       self.counts[self.nBins-1]/sampsize)
        s += "%-24s %6d %8.4f" % ("Total", sampsize, 1)
        return s

    def clear(self):
        'Set all counts to zero'
        self.counts = self.nBins*[0]
    
nold = (-1)
pold = (-1.0)
pc = plog = pclog = en = oldg = 0.0
cof = (76.18009172947146,  -86.50532032941677,  24.01409824083091, \
       -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5)

def gammln(xx):
        """
        Return natural log of Gamma function.
        """

        y = x = xx
        tmp = x + 5.5
        tmp -= (x+0.5)*log(tmp)
        ser = 1.000000000190015
        for j in range(6):
                y += 1
                ser += cof[j]/y
        return -tmp + log(2.5066282746310005*ser/x)

factln_array_size = 101
factln_array = factln_array_size*[0]
def factln(n):
    """
    Return log of factorial of n.
    """
    if n < 0:
        raise ValueError, "argument to factln must be positive"
    if n <= 1:
        return 0
    if n < factln_array_size:
        if factln_array[n] > 0:
            return factln_array[n]
        else:
            factln_array[n] = gammln(n+1.0)
            return factln_array[n]
    else:
        return gammln(n+1.0)


def bnldev(n, pp):
    """
    Return a random deviate drawn from the Binomial distribution
    with parameters n (a positive integer) and pp (a probability).
    """

    # This code was translated from Numerical Recipes in C

    global nold, pold, pc, plog, pclog, en, oldg

    if pp <= 0.5:
        p = pp
    else:
        p = 1.0 - pp

    am = n*p

    if n < 25:
        # Direct method
        bnl = 0.0
        for j in xrange(n):
            if random() < p:
                bnl += 1.0

    elif am < 1.0:
        # Poisson method
        g = exp(-am)
        t = 1.0
        for j in xrange(n):
            t *= random()
            if t < g:
                break
        bnl = j

    else:
        # rejection method
        if n != nold:
            en = n
            oldg = gammln(en + 1.0)
            nold = n;

        if p != pold:
            pc = 1.0 - p
            plog = log(p)
            pclog = log(pc)
            pold = p

        sq = sqrt(1.0*am*pc)

        while True:
            while True:
                angle = pi*random();
                y = tan(angle)
                em = sq*y + am
                if em >= 0 and em < (en+1.0):
                    break
            em = floor(em)
            t = 1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0) \
                                     - gammln(en-em+1.0) \
                                     + em*plog + (en-em)*pclog)
            if random() < t:
                break
        bnl = em

    if p != pp:
        bnl = n - bnl

    return int(round(bnl))

def mnldev(n, p):
    """
    Return a deviate from a Multinomial random variable.
    On entry, n should be a positive integer and p a vector
    of values that are proportional to probabilities.
    """
    yvec = []
    i = 0
    while i < len(p)-1:
        y = bnldev(n, p[i]/float(sum(p[i:])))
        y = int(y)
        n -= y
        yvec.append( y )
        i += 1
    yvec.append(n)
    return yvec

class axis:
    def __init__(self, v, output_size, nticks=5):
        self.low = min(v)
        self.high = max(v)
        self.outsize = int(output_size)

        # Get tick values, the number of character positions needed to
        # print them, and the format to be used.
        self.tickval, self.lblwid, self.tickfmt = \
                      lblaxis(self.low,self.high,nticks)

        # Reset high and low if necessary to accomodate ticks
        if self.tickval[-1] > self.high:
            self.high = self.tickval[-1]
        if self.tickval[0] < self.low:
            self.low = self.tickval[0]

        # Set scale parameter
        if self.high == self.low:
            self.scale = 1.0
        else:
            self.scale = (self.outsize-1.0)/(self.high-self.low)

        # Get the indices of the tick values
        self.tickndx = []
        for tval in self.tickval:
            self.tickndx.append(self.ndx(tval))


    # Return index of val
    def ndx(self, val):
        assert val >= self.low
        assert val <= self.high
        i = int(round(self.scale*(val - self.low)))
        assert i >= 0
        assert i < self.outsize
        return i

maxxcess = 0.6
roundval = [1.0, 1.5, 2.0, 2.5, 5.0]

def lblaxis(asklow, askhigh, askntic):
    """
    Generate labels for an axis.  Input parameters:

    asklow : desired low value of axis
    askhigh: desired high value
    askntic: number of tick marks desired

    Function returns (tic, lng, fmt), where

    tic is a list of floats, the values at which tick marks should be
        placed, and
    lng is the number of character positions needed for printing
        tick labels.
    fmt is a string that specifies the format to be used in printing
        tick labels.

    For example, the tick marks can be printed like this:

      for value in tic:
          print fmt % value

    The range of tick mark values and the number of tick marks should
    be close to the values requested, but will probably not be identical.
    The function tries to generate a list of values that will be pretty.
    """

    high = askhigh
    low = asklow
    assert high > low, "lblaxis: low > high"
    n = askntic

    for i in range(5):
        oldhigh = high
        oldlow = low
        w = (high - low) / (n - 1.0)
        x = log10(w)
        exp10 = floor(x)
        y = 10.0**(x - exp10)
        y = 0.5 * floor(2.0 * y + 0.5)
        # y should now be close to one of the values in roundval
        i0 = 0
        for i in range(1, len(roundval)):
            if abs(y - roundval[i]) < abs(y - roundval[i0]):
                i0 = i

        # width of interval is roundval[i0] * 10^exp10
        w = roundval[i0] * pow(10.0, exp10)

        # calculate new low
        low = ceil(oldlow / float(w)) * w
        if low - asklow > maxxcess * w:
            low = w*(ceil(oldlow / float(w)) - 1.0) # new low

        # calculate new high
        high = w*floor(oldhigh / float(w))
        if askhigh - high > maxxcess * w:
            high = w*(floor(oldhigh / float(w)) + 1.0)
        if oldhigh == high or oldlow == low:
            break

    askntic = n = int((high - low) / float(w) + 1.5) # rounds to nearest int

    eps = (high - low)*0.05

    # fill tick-mark vector and printing format
    tic = []
    lng=frac=0
    for i in range(n):
        this_tic = low + i * w
        l, f, s = prwid(this_tic, eps)
        lng = max(lng, l)
        frac = max(frac, f)
        tic.append(this_tic)

    # format string
    fmt = "%%%d.%df" % (lng, frac)

    return (tic, lng, fmt)


def prwid(w, precision):

    # Step 1: round w to a value determined by the precision
    # needed.
    s = "%f" % w
    s = s.strip('0')
    for frac in range(0,len(s)):
        r = round(w, frac)
        if abs(r-w) <= precision:
            w = r
            break
    fmt = "%%1.%df" % frac
    s = fmt % w
    if '.' in s:
        s = s.rstrip('0')
        s = s.rstrip('.')

    # lng is length of whole string
    lng = len(s)
    if lng == 0:
        lng = 1

    # format string
    fmt = "%%%d.%df" % (lng, frac)
    return (lng, frac, fmt)

# Plots on terminals w/o graphics capabilities.
# On entry: 
#       x is a float vector of x-axis values
#       y is a float vector of y-axis values
# Function writes plot on standard output
def charplot(x, y, nticks=5, output_height=24, output_width=78):

    yax = axis(y, output_height)
    left_pad = yax.lblwid+1
    output_width -= left_pad
    xax = axis(x, output_width)

    blank_row = (output_width)*[' ']
    chrmat = [blank_row[:] for i in range(output_height)]


    # create plot matrix
    for i in range(1, output_height):     # draw side lines
        chrmat[i][0] = '|'
        chrmat[i][output_width - 1] = '|'
    for i in range(output_width):         # draw top & bottom lines
        chrmat[0][i] = '-'
        chrmat[output_height-1][i] = '-'

    for ndx in xax.tickndx:     # tics on x axis 
        chrmat[0][ndx] = '+'
        chrmat[output_height-1][ndx] = '+'

    for ndx in yax.tickndx:     # tics on y axis
        chrmat[ndx][0] = '+'
        chrmat[ndx][output_width-1] = '+'

    for xval, yval in zip(x,y):     # data points 
        chrmat[yax.ndx(yval)][xax.ndx(xval)] = '*'

    # create list of y axis labels
    ylab = []
    for i in range(yax.outsize):
        ylab.append(yax.lblwid*' ')
    for ndx, val in zip(yax.tickndx, yax.tickval):
        ylab[ndx] = yax.tickfmt % val

    assert(yax.outsize == len(chrmat))

    # print chrmat
    for i in range(len(chrmat)):
        ndx = len(chrmat) - i - 1
        row = chrmat[ndx]
        line = ''.join(row)
        print "%s %s" % (ylab[ndx], line)
    s = ''
    curr_pos = -left_pad
    for ndx, val in zip(xax.tickndx, xax.tickval):
        lbl = xax.tickfmt % val
        lbl = lbl.strip()
        if '.' in lbl:
            lbl.rstrip('0')
            lbl.rstrip('.')
        wid = len(lbl)
        halfwid = wid/2
        skip = ndx - halfwid - curr_pos
        s += skip*' '
        s += lbl
        curr_pos += skip + wid
    print s

# Smooth data using an exponentially weighted moving average.  x is
# the data vector to be smoothed, and b is a parameter that controls
# the amount of smoothing.  Choose a value that satisfies 0 <= b < 1.
# If b=0, no smoothing is done.  Larger values of b produce more
# smoothing.
def expmovav(x, b=0.7):
    n = len(x)
    y0 = [(1.0-b)*x[0]]
    for i in range(1, n):
        val = (1.0-b)*x[i] + b*y0[i-1]
        y0.append( val  )

    y1 = n*[None]
    y1[n-1] = (1.0-b)*x[n-1]
    for i in range(n-2, -1, -1):
        val = (1.0-b)*x[i] + b*y1[i+1]
        y1[i] = val
    return [0.5*(yy0+yy1) for yy0,yy1 in zip(y0,y1)]

# Convert a matrix into a printable string
def strmat(m):
    """
    Convert a matrix (i.e. a list of lists) into a printable string.
    
    The items in the matrix can be of any type that Python can
    turn into a string.  Function returns a string in which the
    items of the matrix are right-justified in columns that are
    as narrow as possible.

    Example:

    from pgen import strmat
    m = [[1,'a',3.3333], ["", None, 1.23]]   # define a matrix
    print strmat(m)                          # print it
    """

    # find appropriate width for each column
    colwid = []
    rowmax = 0
    for row in m:
        if len(row) > rowmax:
            rowmax = len(row)
        for i in range(len(row)):
            wid = len(str(row[i]))
            if i >= len(colwid):
                colwid.append(wid)
            else:
                colwid[i] = max(wid, colwid[i], len(str(i)))

    # create string
    s = ""
    for i in range(len(m)-1, -1, -1):
        row = m[i]
        s += "%2d [" % i
        for i in range(len(row)):
            fmt = " %%%ds" % colwid[i]
            s += fmt % str(row[i])
        while i < rowmax-1:
            fmt = " %%%ds" % colwid[i]
            s += fmt % ""
            i += 1
        s += "]\n"
    s += "%2s  " % ""
    for i in range(rowmax):
        fmt = " %%%dd" % colwid[i]
        s += fmt % i
    return s

def scatsmooth(x, y, n=10, lo_x=None, hi_x=None):
    """
    Smooth a set of (x,y) data by dividing the range of
    x values into n equally-spaced bins, and calculating
    the average y-value within each bin.

    Function returns (x_mid, y_mean)
    where x_mid is a vector of n midpoints of bins,
    and y_mean is the vector of means within the bins.
    """

    if lo_x == None:
        lo_x = min(x)

    if hi_x == None:
        hi_x = max(x)

    assert hi_x > lo_x
    scale = (n-1.0)/(hi_x-lo_x)

    ym = n*[0.0]  # vector of means
    yn = n*[0]    # vector of sample sizes
    step = (hi_x-lo_x)/float(n)

    for xval, yval in zip(x,y):
        ndx = int(round(scale*(xval-lo_x)))
        ym[ndx] += yval
        yn[ndx] += 1

    x_midpoint = []
    for ndx in range(n):
        if yn[ndx] > 0:
            ym[ndx] /= float(yn[ndx])
        x_midpoint.append(lo_x + (ndx+0.5)*step)

    return (x_midpoint, ym, yn)

def expran(h):
    """
    Return a random variate drawn from the exponential distribution
    with hazard h and mean 1/h.
    """
    return( -log(1.0-random())/h ) 

# poisson distribution

# Variables whose values are saved between calls.
sq=lnLam=g=oldLam=0.0

def poidev(lam):
    """
    Return a random variate drawn from the Poisson distribution
    with mean lam.
    """
    global sq, lnLam, g, oldLam
    
    pt9 = 0.9
    zero = 0.0
    one = 1.0
    two = 2.0
    twelve = 12.0

    rval = t = y = 0.0

    if (lam < twelve):
        if (lam != oldLam):
            oldLam = lam
            g = exp(-lam)

        rval = -one
        t = one
        while True:
            rval += 1
            t *= random()
            if t <= g:
                break

    else:
        if lam != oldLam:
            oldLam = lam
            sq = sqrt(two * lam)
            lnLam = log(lam)
            g = lam * lnLam - gammln(lam + one)

        while True:
            while True:
                y = tan(pi * random())
                rval = sq * y + lam
                if rval >= zero:
                    break

            rval = floor(rval)
            t = pt9 * (one + y * y) * exp(rval * lnLam - \
                                          gammln(rval + one) - g)
            if random() <= t:
                break

    return int(round(rval))

def bico(n, k):
    """
    Return binomial coefficient: n choose k
    From Numerical Recipes.
    """
    return floor(0.5+exp(factln(n) - factln(k) - factln(n-k)))

def binom_prob(x, n, p):
    """
    Return probability that binomial(n,p) is equal to x.
    """
    if p == 0:
        if x == 0:
            return 1
        else:
            return 0
    if p == 1:
        if x == n:
            return 1
        else:
            return 0
    
    rval = factln(n) - factln(x) - factln(n-x)
    rval += x*log(p)
    rval += (n-x)*log(1.0-p)
    return exp(rval)

class hapmap_snp:
    """
    Represents data from a single HapMap SNP.  Constructor takes
    a single argument, a line of data in HapMap format.

    Instance variables:

    alleles : list of the alleles present at locus

    chromosome: an integer representing the chromosome

    position: an integer; the position of the SNP on the chromosome

    trios: boolean; True if these are trio data (THIS HAS BEEN DISABLED)

    gtype: a list of genotype data.  Each item in the list is an integer,
    the number of copies of alleles[0] in an individual genotype.  The
    Python "None" value is used for genotypes with missing data.

    sampleSize: number of values in gtype, excluding missing values.

    mean: the mean of gtype

    variance: the variance of gtype
    """
    def __init__(self, line):

        # for debug
        self.input = line[:]
        
        # Get rid of white space at beginning and end
        line = line.strip()

        # Skip comments
        if line[0] == '#':
            raise ValueError, "Can't make a snp from a comment"

        # Turn line into a list of fields.
        line = line.split()

        # Skip header
        if line[0] == "rs#":
            raise ValueError, "Can't make a snp from HapMap header"

        self.id = line[0]
        self.alleles = line[1].split('/')
        self.chromosome = line[2][3:]
        self.position = int(line[3])

        #trios_string = line[9].split(':')[4].split('-')[2]
        #if trios_string == 'trios':
        #    self.trios = True
        #else:
        #    self.trios = False

        # Data begins with field 11
        self.gtype = line[11:]

        self.mean = self.variance = 0.0
        self.sampleSize = 0
        # Recode as numeric by counting copies of allele[0]
        for i, x in enumerate(self.gtype):
            if 'N' in x:
                self.gtype[i] = None
            else:
                self.sampleSize += 1
                self.gtype[i] = x.count(self.alleles[0])
                self.mean += self.gtype[i]
                self.variance += self.gtype[i] * self.gtype[i]

        if self.sampleSize == 0:
            return

        self.mean /= float(self.sampleSize)
        self.variance /= float(self.sampleSize)
        self.variance -= self.mean*self.mean

    def __str__(self):
        #if self.trios:
        #    trios_string = "trios"
        #else:
        #    trios_string = "notrios"
        
        return ':'.join(['/'.join(self.alleles), \
                 str(self.chromosome), \
                 str(self.position), \
                 #trios_string, \
                 str(self.gtype)])

    def __getitem__(self, i):
        return self.gtype[i]

    def __len__(self):
        return len(self.gtype)

    def __delitem__(self, i):
        del self.gtype[i]

    def __contains__(self, i):
        return (i in self.gtype)

def gtype_monomorphic(gtype_list):
    """
    Return True if list of genotypes is monomorphic, False if polymorphic.
    Genotypes must be coded as integers in [0,1,2].
    """
    got0 = False
    got2 = False
    for gtype in gtype_list:
        if gtype == None:
            continue
        elif gtype==1:
            # heterozygote: locus is polymorphic
            return False
        elif gtype == 0:
            # 00 homozygote
            got0 = True
            if got2:
                return False # locus is polymorphic
        else:
            assert gtype == 2
            got2 = True
            if got0:
                return False # locus is polymorphic
    return True # locus is monomorphic
            

class hapmap_dataset:
    """
    hapmap_dataset is a class representing the data from a single file.

    Instance variables:

    filename: the name of the original HapMap data file

    snps: a list of objects of type hapmap_snp, each containing data
    from a single snp.

    Methods:

    find_position(pos) Return index of first SNP whose position on
    chromosome is >= pos.  If no such SNP exists, then return len(self)-1. 
    """
    def  __init__(self, infile, dirpath=hapmap_dir):
        if dirpath==None:
            dirpath = get_dirpath()
        self.filename = dirpath + '/' + infile
        try:
            linelist = open(self.filename)
        except:
            raise IOError, "hapmap_dataset.__init__:Can't open '%s' for input" \
                % infile
        self.snps = []
        nmiss = 0
        for line in linelist:
            line = line.strip()
            if line[0] == '#' or line[:3] == "rs#":
                continue
            snp = hapmap_snp(line)
            if None in snp.gtype:
                nmiss += 1
            elif gtype_monomorphic(snp.gtype):
                nmiss += 1
            else:
                self.snps.append(snp)
        print "hapmap_dataset:" + \
            "Skipped %d SNPs w/ missing or monomorhic data; kept %d" % \
            (nmiss, len(self.snps))
        return

    def __str__(self):
        return ':'.join([self.filename,str(len(self.snps))])

    def __getitem__(self, i):
        return self.snps[i]

    def __len__(self):
        return len(self.snps)

    def __delitem__(self, i):
        del self.snps[i]

    def find_position(self, pos):
        """
        Return index of SNP whose position on chromosome
        is closest to pos.
        """
        class poslist:
            def __init__(self, hf):
                self.hf = hf

            def __len__(self):
                return len(self.hf)

            def __getitem__(self, i):
                return self.hf[i].position
        pl = poslist(self)
        ndx = bisect_left(pl, pos)
        if ndx == len(self.snps):
            ndx = len(self.snps)-1
        elif ndx > 0:
            if abs(pos - pl[ndx-1]) < abs(pos - pl[ndx]):
                ndx = ndx-1
        return ndx
            
    def __contains__(self, pos):
        i = self.find_position(pos)
        return (self.snps[i].position == pos)

def hapmap_fname(chromosome, pop, dirpath=hapmap_dir):
    """
    Return a HapMap file name corresponding to the chromosome, population,
    and directory that are specified in the list of arguments.
    """
    pat = "^genotypes_chr%s_%s_.*_fwd\.txt$" % (str(chromosome), pop)
    repat = re.compile(pat) 
    fnames = []
    if dirpath==None:
        dirpath = get_dirpath()
    for fname in os.listdir(dirpath):
        m = repat.match(fname)
        if m:
            fnames.append(fname)
    if len(fnames) == 0:
        raise IOError, "hapmap_fname(%s,%s,%s):Can't find input file" \
                % (chromosome,pop,dirpath)
    elif len(fnames) > 1:
        raise IOError, "can't identify unique file: %s" % fnames
    else:   
        return fnames[0]

def test_monomorphic():
    assert gtype_monomorphic([0,0,0,0]) == True
    assert gtype_monomorphic([2,2,2,2]) == True
    assert gtype_monomorphic([1,1,1,1]) == False
    assert gtype_monomorphic([2,1,2,2]) == False

def test_hapmap_fname(verbose=False):
    for i in range(1,23):
        fname = hapmap_fname(i,'CEU')
        if verbose:
            print 'for chromosome %d CEU: %s' % (i, fname)
    fname = hapmap_fname('X','CEU')
    if verbose:
        print 'for chromosome %d CEU: %s' % (i, fname)
    fname = hapmap_fname('Y','CEU')
    if verbose:
        print 'for chromosome %s CEU: %s' % ('Y', fname)
    print "hapmap_fname OK"

def test_hapmap_dataset(verbose=False):
    hf = hapmap_dataset(hapmap_fname(21, 'CEU'))
    if verbose:
        print hf
        print hf.snps[0]
    i = randrange(len(hf))    # index of random snp
    pos = hf[i].position
    assert(pos in hf)
    assert(not (pos+1 in hf))
    print "hapmap_dataset OK"

def test_mnldev():
    nreps = 1000
    n = 100
    p = [0.1, 0.3, 0.1]
    ex = [n*p[i]/float(sum(p)) for i in range(len(p))]
    for i in xrange(10):
        x = mnldev(n, p)
        print x
    
def test_poidev():
    lam = 3
    nreps = 1000
    maxtab = 10
    tab = (maxtab+1)*[0]
    for i in xrange(nreps):
        x = poidev(lam)
        if x <= maxtab:
            tab[x] += 1
    print tab

    chisq = 0
    for x, nx in enumerate(tab):
        lnpx = x*log(lam) - lam - gammln(x+1)
        px = exp(lnpx)
        ex = nreps*px
        dx = nx - ex
        chisq += dx*dx/ex
        
    print "test_poidev: ChiSq[%d] = %f" % (maxtab, chisq)

    # Test against 0.05 critical value for 10 df
    if chisq > 18.307:
        print "poidev failed"
    else:
        print "poidev OK"

def test_expran():
    nreps = 1000
    h = 0.5
    m1 = m2 = 0.0
    for i in xrange(nreps):
        r = expran(h)
        m1 += r
        m2 += r*r
    m1 /= float(nreps)
    m2 /= float(nreps)
    Em1 = 1.0/h     # expected value of m1
    Em2 = 2.0/(h*h) # expected value of m2

    print "m1=%f Em1=%f m2=%f Em2=%f" % (m1, Em1, m2, Em2)
        
def test_scatsmooth():
    from random import random
    n = 100
    x = [i for i in range(n)]
    y = [i + random() for i in range(n)]

    xs, ys, yn = scatsmooth(x,y)
    for i in range(len(xs)):
        print xs[i], ys[i], yn[i]

def test_charplot():
    yv = [1 - i/20.0 + i*i/400.0 for i in range(20)]
    xv = [i for i in range(20)]
    charplot(xv, yv)

def do_bnldev_test(n, p, nreps):
    counts = (n+1)*[0]
    for i in xrange(nreps):
        x = bnldev(n, p)
        counts[x] += 1

    sumchisq = 0.0
    for x, count in enumerate(counts):
        pr = binom_prob(x, n, p)
        exp_count = nreps*pr
        dev = count - exp_count
        chisq = dev*dev/float(exp_count)
        sumchisq += chisq
        #print "x=%d count=%d exp_count=%f" % (x, count, exp_count)
    return sumchisq, len(counts)

def test_bnldev():
    "Test function bnldev"
    nreps = 10000

    chisq, df = do_bnldev_test(10, 0.2, nreps)  # direct method
    print "test_bnldev direct method: chisq[%d]=%f" % (df, chisq)

    chisq, df = do_bnldev_test(100, 0.005, nreps) # poisson method
    print "test_bnldev poisson method: chisq[%d]=%f" % (df, chisq)

    chisq, df = do_bnldev_test(100, 0.2, nreps) # rejection method
    print "test_bnldev rejection method: chisq[%d]=%f" % (df, chisq)

    chisq, df = do_bnldev_test(100, 0.9, nreps) # rejection method again
    print "test_bnldev rejection method again: chisq[%d]=%f" % (df, chisq)

def test_Tabulator():
    'Test class Tabulator'

    tab = Tabulator(1.0, 2.0, 4)
    assert tab[tab.bin(1.1)] == 0
    tab += 1.1
    assert tab[tab.bin(1.1)] == 1
    tab += 1.1
    assert tab[tab.bin(1.1)] == 2
    tab += 1.222
    assert tab[tab.bin(1.1)] == 3
    tab += 1.3
    assert tab[tab.bin(1.1)] == 3
    assert tab[tab.bin(1.3)] == 1
    assert tab[tab.bin(1.8)] == 0
    assert tab.sampleSize() == 4
    tab += 99.9
    assert tab[tab.bin(99.9)] == 1
    assert tab.sampleSize() == 5
    tab += -9999.9
    assert(tab.bin(-9999.9) == 0)
    assert tab[tab.bin(9999.9)] == 1
    assert tab.sampleSize() == 6
    print tab

    tab.clear()
    assert tab.sampleSize() == 0
    assert tab[tab.bin(1.3)] == 0

    # Test Tabulator.merge
    x = Tabulator()
    y = Tabulator()
    x += 0.1
    y += 0.3
    x += 0.5
    y += 0.7
    x += 0.9
    assert x.sampleSize() == 3
    assert y.sampleSize() == 2
    x.merge(y)
    assert x.sampleSize() == 5
    assert x[x.bin(0.1)] == 1
    assert x[x.bin(0.2)] == 1
    assert x[x.bin(0.3)] == 1
    assert x[x.bin(0.4)] == 1
    assert x[x.bin(0.5)] == 1

    x.clear()
    assert x.sampleSize() == 0
    print "Tabulator OK"

def test_strmat():
    print strmat([[0,1,2],[3,4,5,6],[6,"aaa",8]])

if __name__ == '__main__':
    test_monomorphic()
    test_Tabulator()
    test_bnldev()
    test_charplot()
    test_scatsmooth()
    test_expran()
    test_poidev()
    test_mnldev()
    test_hapmap_fname()
    test_hapmap_dataset()
    test_strmat()

