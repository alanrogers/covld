from estimate_ld import *
import sys

def usage():
    print "Usage: covld [options] [inputfile]"
    print
    print "where options may include:"
    print " -h"
    print " --help  Print this message"
    print " -es     Use Excoffier-Slatkin EM algorithm"
    print
    print "If input file is omitted, reads from standard input."
    exit(1)

estimators = [Estimator("Rogers-Huff", get_r)]
infile = sys.stdin
do_es = False

# Process command line args
for i in range(len(sys.argv)):
    if sys.argv[i] == "-h" or sys.argv[i]=="--help":
        usage()
    elif sys.argv[i] == "-es":
        do_es = True
    else:
        try:
          infile = open(sys.argv[i])
        except:
            traceback.print_exc()

y = []
z = []

for line in infile:
    line = line.strip()
    if len(line)==0:
        continue
    if line[0] == '#':
        continue
    line = line.split()
    yval = int(line[0])
    zval = int(line[1])
    y.append(yval)
    z.append(zval)

print "%-10s %-8s" % ("Estimator", "Estimate")

rRH = get_r(y,z)
print "%-10s %-8.4f" % ("r_RH", rRH)

if do_es:
    rES = esem_r(y,z)
    print "%-10s %-8.4f" % ("r_ES", rES)
