CC = mpicc
FC = mpif90
CFLAGS = -O3 -g -Wall -fopenmp -DAdd_ #validate Add_=1 in the code
FFLAGS = -O3 -g -Wall -cpp -ffree-line-length-none -fopenmp

#===========================================================
# MUMPS paths
mumps = $(HOME)/Install/MUMPS_5.7.3
IMUMPS = -I$(mumps)/include
LMUMPS = -L$(mumps)/lib -lzmumps -lmumps_common

pord = $(mumps)/PORD
IPORD = -I$(pord)/include
LPORD = -L$(pord)/lib -lpord

IMETIS    = -I/usr/include/parmetis
LMETIS = -L/usr/lib  -lparmetis -lmetis #parallel version
#IMETIS    = -I/usr/include/metis #sequential version

ISCOTCH   = -I/usr/include/scotch
LSCOTCH   = -L/usr/lib -lptesmumps -lptscotch -lptscotcherr #parallel version
#LSCOTCH   = -L/usr/lib -lesmumps -lscotch -lscotcherr #sequential version

#===========================================================
BIN = .

# Order: MUMPS -> PORD -> SCOTCH -> ScaLAPACK -> LAPACK -> BLAS -> system libraries
LIB = $(LMUMPS) $(LPORD) $(LMETIS) $(LSCOTCH) -lscalapack-openmpi -llapack -lblas -lpthread -lgfortran -lm -lmpi -lmpi_mpifh -lfftw3 
INC = -I. $(IMUMPS) $(IPORD) $(LMETIS) $(ISCOTCH)
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)


all: clean main

main:	$(OBJ)
	$(CC) $(CFLAGS) -o $(BIN)/main $(OBJ) $(LIB)

%.o: %.c 
	$(CC) $(CFLAGS) -c $^ -o $@ $(INC) $(LIB)

clean:
	find . -name "*.o"   -exec rm {} \;
	find . -name "*.c%"  -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
	find . -name "*~"    -exec rm {} \;
	find . -name "\#*"   -exec rm {} \;
	rm -f $(OBJ) main *.mod


