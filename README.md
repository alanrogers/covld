covld
=====

Estimate linkage disequilibrium between unphased loci

The software in this directory implements a method for estimating linkage
disequilibrium, which is described in the following publication.

    @Article{Rogers:G-182-839,
      author =       {Rogers, Alan R. and Chad Huff},
      title =        {Linkage Disequilibrium between Loci of Unknown Phase},
      journal =      {Genetics},
      year =         2009,
      volume =   182,
      number =   3,
      pages =    {839--844},
      month =    {July},
    }

    FILE                 DESCRIPTION
    
    estimate_ld.py      Implements various estimators of LD.  If you write your
                        own Python program, import this as a package.
    
    myexcept.py         Defines exceptions used by the other programs.
    
    tst_r.py            Tests the functions in estimate_ld.py.
    
    covld.py            Reads a data file and estimates LD using the Rogers-Huff
                        method and (optionally) the Excoffier-Slatkin method.
                        For usage info, type "python covld.py -h".
    
    covld-dat.txt       Sample input file for covld.py.  The input format is
                        described in the header of this file.
    
This software is all in the public domain.

Alan R. Rogers
