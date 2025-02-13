########################################
# Makefile for Spin Calculator         #
########################################

# C compiler and flags
CC=gcc
CFLAGS=-Wall -Wextra -std=c99 -fopenmp -lm

# Command to use for linking and executable
LD=gcc
LDFLAGS =
EXE=pertubation.x

OBJECTS=mt19937ar.o data.o spin_functions.o pertubation.o

all: pertubation

# Default build target
pertubation : $(OBJECTS) 
	$(CC) $(CFLAGS) -o $(EXE) $(OBJECTS) $(LDFLAGS)

# Purge build files and executable
clean :
	rm -rf *.o *.mod $(EXE)

%.o: %.c 
	$(CC) $(CFLAGS) -c -o $@ $<

