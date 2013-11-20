# Linux flags
CC=g++
CFLAGS= -O3 -Wall -W -g
#Debugging flags
#CFLAGS= -g  -Wall -W 
LFLAGS= -lm
#Linux flags for intel compiler
#CC=icc
#CFLAGS= -fast -Wall
#LFLAGS= -lm

# SGI 
#CC=cc
#CFLAGS = -O2 -mips4
#LFLAGS= -lm

# Alpha
#CC=cc
#CFLAGS = -fast -arch host -tune host
#LFLAGS = -fast -non_shared -om -WL,-om_no_inst_sched -lm

SRCS=ensemble.cpp  gmxmdpio.cpp  gmxpullpot.cpp  gwham_gromacs463.cpp  gwham_gromacs463_umb.cpp mc.cpp ran2.c
OBJSGMX463=ensemble.o  gmxmdpio.o  gmxpullpot.o  gwham_gromacs463.o  
OBJSGMX463UMB=ensemble.o  gmxmdpio.o  gmxpullpot.o  gwham_gromacs463_umb.o 
OBJSMC=ensemble.o mc.o  
HEADER_FILES=densityofstate.hpp  ensemble.hpp  fileio.hpp  fileio_utils.hpp  gmxmdpio.hpp  gmxpullpot.hpp  gmxxvgio.hpp  gnarray.hpp  gwham.hpp  hamiltonian.hpp  typedefs.hpp mc.hpp
all: gwham_gromacs463 gwham_gromacs463_umb mc

ran2.o: ran2.c
	$(CC) $(CFLAGS) -c ran2.c

gwham_gromacs463_umb.o:  gwham.hpp densityofstate.hpp gwham_gromacs463_umb.cpp hamiltonian.hpp gmxxvgio.hpp fileio.hpp gnarray.hpp typedefs.hpp ensemble.o gmxmdpio.o 
	$(CC) $(CFLAGS) -c gwham_gromacs463_umb.cpp

gwham_gromacs463.o: gwham.hpp densityofstate.hpp gwham_gromacs463_umb.cpp hamiltonian.hpp gmxxvgio.hpp fileio.hpp gnarray.hpp typedefs.hpp ensemble.o gmxmdpio.o 
	$(CC) $(CFLAGS) -c gwham_gromacs463.cpp

mc.o: gwham.hpp densityofstate.hpp hamiltonian.hpp gnarray.hpp typedefs.hpp mc.hpp mc.cpp ensemble.o ran2.o
	$(CC) $(CFLAGS) -c mc.cpp

ensemble.o:  ensemble.cpp ensemble.hpp typedefs.hpp 
	$(CC) $(CFLAGS) -c ensemble.cpp

gmxpullpot.o:  gmxpullpot.cpp gmxpullpot.hpp typedefs.hpp 
	$(CC) $(CFLAGS) -c gmxpullpot.cpp

gmxmdpio.o:  gmxmdpio.cpp gmxmdpio.hpp typedefs.hpp fileio_utils.hpp 
	$(CC) $(CFLAGS) -c gmxmdpio.cpp

gwham_gromacs463 : $(OBJSGMX463) 
	$(CC) $(CFLAGS) -o gwham_gromacs463 $(OBJSGMX463) $(LFLAGS) 

gwham_gromacs463_umb : $(OBJSGMX463UMB) 
	$(CC) $(CFLAGS) -o gwham_gromacs463_umb $(OBJSGMX463UMB) $(LFLAGS) 

mc : $(OBJSMC) 
	$(CC) $(CFLAGS) -o mc $(OBJSMC) $(LFLAGS)

#driver: bootstrap.o driver.o $(NR)/ran2.o $(NR)/locate.o
#	$(CC) $(CFLAGS) -o driver bootstrap.o driver.o ran2.o locate.o $(LFLAGS)

depend:
	makedepend -- $(CFLAGS) -- $(SRCS)

clean:
	rm -f gwham_gromacs463 gwham_gromacs463_umb mc *.o .*.swp *~

remake:
	make clean
	make depend
	make all

tags: 
	ctags *

# DO NOT DELETE
