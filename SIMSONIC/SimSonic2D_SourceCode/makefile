#---------------------------------------
# Resources section - Symbol definitions
#---------------------------------------

CC = cc -O3 -fopenmp 
LIB = SimSonic2D.a
AR  = ar -r
DIR = ./

#----------------------
# Makefile head section
#----------------------

all: 	simsonic2D

clean:
	rm -f *.a *.o simsonic2D

build:
	make clean all

simsonic2D: $(LIB)
	$(CC) -o simsonic2D $(LIB) -lm -static -static-libgcc

#----------------------------------
# Dependencies for the source files
#----------------------------------

SIMSONIC2D_C_DEP = \
	$(DIR)/SimSonic2D.h \
	$(DIR)/PMLBordVitesses.c \
	$(DIR)/PMLBordContraintes.c \
	$(DIR)/PMLCoin2DVitesses.c \
	$(DIR)/PMLCoin2DContraintes.c \

#--------------------------------
# Compilation of the source files
#--------------------------------

simsonic2D.o : SimSonic2D.c $(SIMSONIC2D_C_DEP)
	$(CC) -c SimSonic2D.c -o simsonic2D.o

#----------------------------------
# Library dependencies and creation
#----------------------------------

LIB_DEP = \
	simsonic2D.o \

$(LIB) : $(LIB_DEP)
	$(AR) $(LIB) $(LIB_DEP)

