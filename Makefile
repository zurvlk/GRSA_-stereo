#Makefile

# CXX = icc
# CXXFLAGS = -Wall -parallel -lm -fPIC

CC=gcc
all:grsa

grsa:graph.o bmp.o ford_fulkerson.o queue_stack.o grsa_drv.o grsa.o -lm
	$(CC) graph.o bmp.o ford_fulkerson.o queue_stack.o grsa_drv.o grsa.o -o grsa -lm -g3
.c.o:
	$(CC) -c $< -Wall -l -g3
remove:
	rm -f *.o
clean:
	rm -f *.o *~ grsa
