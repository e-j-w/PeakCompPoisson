CFLAGS   = -O -g -Wall -fPIC -ansi -lm
ROOT = $(shell $(ROOTSYS)/bin/root-config --glibs) $(shell $(ROOTSYS)/bin/root-config --cflags)

all: peak_comp_poisson_root

peak_comp_poisson_root: peak_comp_poisson_root.c peak_comp_poisson_root.h read_config.c 
	g++ peak_comp_poisson_root.c $(CFLAGS) $(ROOT) -o peak_comp_poisson_root
clean:
	rm -rf *~ *.o peak_comp_poisson_root *tmpdatafile*
