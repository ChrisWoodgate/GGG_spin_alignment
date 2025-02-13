########################################
# Makefile for Spin Calculator         #
########################################

# C compiler and flags
CC=gcc
CFLAGS=-Wall -Wextra -O3 -std=c99 -fopenmp -lm

# Command to use for linking and executable
LD=gcc
LDFLAGS =
EXE=GGG_spins.x

OBJECTS=mt19937ar.o data.o spin_functions.o spin_calculator.o

all: spin_calculator

# Default build target
spin_calculator : $(OBJECTS) 
	$(CC) $(CFLAGS) -o $(EXE) $(OBJECTS) $(LDFLAGS)

# Purge build files and executable
clean :
	rm -rf *.o *.mod $(EXE)

%.o: %.c 
	$(CC) $(CFLAGS) -c -o $@ $<

