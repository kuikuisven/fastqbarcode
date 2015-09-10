SHELL = /bin/sh

CC              = clang 

default: fastqbarcode 

fastqbarcode: main.c
	${CC}  -c main.c 
	${CC} -o fastqbarcode  main.o
	rm *.o