# Compiler 
CC = gcc
CFLAGS = -Wall -O2

# Include
INCLUDE += -I/home/x_fedlu/gsl/include

# Libraries
LIB = -static -L/home/x_fedlu/gsl/lib

# Files
SRC = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(SRC))

# Name of the executable
EXECUTABLE = hsmc

# Compile
all: $(EXECUTABLE)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $<

# Link
$(EXECUTABLE): $(OBJS)
	$(CC) $(LIB) $^ -o $@ -lgsl -lgslcblas -lm
	rm *.o
