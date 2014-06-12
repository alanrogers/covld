vpath %.h ..
.c.o:
	$(CC) $(CFLAGS) $(INCL) -c -o ${@F}  $<
FLAGS1 = -Wall -Wshadow -Wpointer-arith \
  -Wcast-qual -Wcast-align -Wwrite-strings  \
  -Waggregate-return -Wstrict-prototypes -Wmissing-prototypes \
  -Wmissing-declarations -Wnested-externs  -pedantic
FLAGS2 = -Wconversion -std=c99 -D_GNU_SOURCE
FLAGS3 = -DNDEBUG -O3
PROFILE = -pg
CFLAGS = -g $(FLAGS1) $(FLAGS2) $(FLAGS3)
#CFLAGS = -g $(FLAGS1) $(FLAGS2)
#CFLAGS = -g $(FLAGS1) $(FLAGS2) $(FLAGS3) $(PROFILE)
CC = gcc
LIB = -pg -lm -lgsl -lgslcblas
EXE = chadsim2

TST_R_OBJ = tst_r.o estimate_ld.o
tst_r : $(TST_R_OBJ)
	$(CC) -o $@ $(TST_R_OBJ) $(LIB)


CHADSIM2_OBJ = chadsim2.o estimate_ld.o
chadsim2 : $(CHADSIM2_OBJ)
	$(CC) -o $@ $(CHADSIM2_OBJ) $(LIB)

XESTIMATE_LD_OBJ = xestimate_ld.o estimate_ld.o
xestimate_ld : $(XESTIMATE_LD_OBJ)
	$(CC) -o $@ $(XESTIMATE_LD_OBJ) $(LIB)

TRY_OBJ = try.o
try : $(TRY_OBJ)
	$(CC) -o $@ $(TRY_OBJ) $(LIB)

clean :
	rm -f *.o core *~ *.tmp chadsim2 chadsim

.PHONY: clean

