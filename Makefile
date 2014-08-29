# Linux flags
CC=g++
#CFLAGS= -Wall -W -O3 
#Debugging flags
CFLAGS= -g -Wall -W -std=c++11
LFLAGS= -lm -std=c++11
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

SRCS=ensemble.cpp gmxmdp.hpp gmxmdpio.cpp fileio.cpp  gmxpullpot.cpp  gwham_gromacs463_umb.cpp mc.cpp ran2.c gmbar.cpp gmbar_gromacs463_umb.cpp
OBJSGMX463UMB=ensemble.o  gmxmdp.o gmxmdpio.o fileio_utils.o fileio.o  gmxpullpot.o  gwham_gromacs463_umb.o 
OBJSGMX463UMBFB=ensemble.o  gmxmdp.o gmxmdpio.o fileio_utils.o fileio.o  gmxpullpot.o  gwham_gromacs463_umbfb.o 
OBJSMBARGMX463UMB=ensemble.o  gmxmdp.o gmxmdpio.o fileio_utils.o fileio.o  gmxpullpot.o  gmbar_gromacs463_umb.o gmbar.o
OBJSUMB=ensemble.o fileio_utils.o fileio.o gwham_umb.o 
OBJSMC=ensemble.o mc.o fileio_utils.o
OBJSGMBAR=gmbar.o 
HEADER_FILES=densityofstate.hpp  ensemble.hpp  fileio.hpp  fileio_utils.hpp mdp.hpp gmxmdp.hpp gmxmdpio.hpp  gmxpullpot.hpp  gmxxvgio.hpp  gnarray.hpp  gwham.hpp  hamiltonian.hpp  typedefs.hpp mc.hpp gmbar.hpp 
all: gwham_gromacs463_umb gwham_gromacs463_umbfb mc gmbar_gromacs463_umb gwham_umb 

ran2.o: ran2.c
	$(CC) $(CFLAGS) -c ran2.c

hamiltonian.o : hamiltonian.cpp hamiltonian.hpp functor.hpp 
	$(CC) $(CFLAGS) -c hamiltonian.cpp

gwham_umb.o:  gwham.hpp densityofstate.hpp gwham_umb.cpp hamiltonian.hpp fileio.hpp fileio_utils.hpp gnarray.hpp typedefs.hpp ensemble.o
	$(CC) $(CFLAGS) -c gwham_umb.cpp

gwham_gromacs463_umb.o:  gwham.hpp densityofstate.hpp gwham_gromacs463_umb.cpp hamiltonian.hpp gmxxvgio.hpp fileio.hpp fileio_utils.hpp gnarray.hpp typedefs.hpp ensemble.o gmxmdpio.o 
	$(CC) $(CFLAGS) -c gwham_gromacs463_umb.cpp

gwham_gromacs463_umbfb.o:  gwham.hpp densityofstate.hpp gwham_gromacs463_umbfb.cpp hamiltonian.hpp gmxxvgio.hpp fileio.hpp fileio_utils.hpp gnarray.hpp typedefs.hpp ensemble.o gmxmdpio.o 
	$(CC) $(CFLAGS) -c gwham_gromacs463_umbfb.cpp

gmbar_gromacs463_umb.o:  gmbar.hpp densityofstate.hpp gmbar_gromacs463_umb.cpp hamiltonian.hpp gmxxvgio.hpp fileio.hpp fileio_utils.hpp gnarray.hpp typedefs.hpp ensemble.o gmxmdpio.o 
	$(CC) $(CFLAGS) -c gmbar_gromacs463_umb.cpp

gmbar.o:  gmbar.hpp gmbar.cpp 
	$(CC) $(CFLAGS) -c gmbar.cpp

mc.o: gwham.hpp densityofstate.hpp hamiltonian.hpp gnarray.hpp typedefs.hpp mc.hpp mc.cpp ensemble.o ran2.o fileio_utils.hpp
	$(CC) $(CFLAGS) -c mc.cpp

ensemble.o:  ensemble.cpp ensemble.hpp typedefs.hpp 
	$(CC) $(CFLAGS) -c ensemble.cpp

gmxpullpot.o:  gmxpullpot.cpp gmxpullpot.hpp typedefs.hpp 
	$(CC) $(CFLAGS) -c gmxpullpot.cpp

fileio_utils.o: fileio_utils.hpp fileio_utils.cpp
	$(CC) $(CFLAGS) -c fileio_utils.cpp

fileio.o: fileio_utils.hpp fileio_utils.cpp fileio.hpp fileio.cpp
	$(CC) $(CFLAGS) -c fileio.cpp 

gmxmdpio.o: fileio_utils.o fileio.o gmxmdpio.hpp gmxmdpio.cpp 
	$(CC) $(CFLAGS) -c gmxmdpio.cpp

gmxmdp.o: fileio_utils.o fileio.o mdp.hpp gmxmdp.hpp gmxmdp.cpp functor.hpp
	$(CC) $(CFLAGS) -c gmxmdp.cpp

gwham_gromacs463_umb : $(OBJSGMX463UMB) 
	$(CC) $(CFLAGS) -o gwham_gromacs463_umb $(OBJSGMX463UMB) $(LFLAGS) 

gwham_gromacs463_umbfb : $(OBJSGMX463UMBFB) 
	$(CC) $(CFLAGS) -o gwham_gromacs463_umbfb $(OBJSGMX463UMBFB) $(LFLAGS) 

gmbar_gromacs463_umb : $(OBJSMBARGMX463UMB) 
	$(CC) $(CFLAGS) -o gmbar_gromacs463_umb $(OBJSMBARGMX463UMB) $(LFLAGS) 

gwham_umb : $(OBJSUMB) 
	$(CC) $(CFLAGS) -o gwham_umb $(OBJSUMB) $(LFLAGS) 

mc : $(OBJSMC) 
	$(CC) $(CFLAGS) -o mc $(OBJSMC) $(LFLAGS)

#driver: bootstrap.o driver.o $(NR)/ran2.o $(NR)/locate.o
#	$(CC) $(CFLAGS) -o driver bootstrap.o driver.o ran2.o locate.o $(LFLAGS)

depend:
	makedepend -- $(CFLAGS) -- $(SRCS)

clean:
	rm -f gwham_umb gwham_gromacs463_umbfb gwham_gromacs463_umb gmbar_gromacs463_umbfb mc *.o .*.swp *~

remake:
	make clean
	make depend
	make all

tags: $(SRCS) $(HEADER_FILES) 
	ctags *
