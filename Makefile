# Linux flags
CC=g++
CFLAGS= -Wall -W -O3 
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

SRCS=ensemble.cpp  gmxmdpio.cpp fileio.cpp  gmxpullpot.cpp  gwham_gromacs463_umb.cpp mc.cpp ran2.c gmbar.cpp gmbar_gromacs463_umb.cpp
OBJSGMX463UMB=ensemble.o  gmxmdpio.o fileio.o  gmxpullpot.o  gwham_gromacs463_umb.o 
OBJSGMX463UMBFB=ensemble.o  gmxmdpio.o fileio.o  gmxpullpot.o  gwham_gromacs463_umbfb.o 
OBJSMBARGMX463UMB=ensemble.o  gmxmdpio.o fileio.o  gmxpullpot.o  gmbar_gromacs463_umb.o gmbar.o
OBJSMC=ensemble.o mc.o  
OBJSGMBAR=gmbar.o 
HEADER_FILES=densityofstate.hpp  ensemble.hpp  fileio.hpp  fileio_utils.hpp  gmxmdpio.hpp  gmxpullpot.hpp  gmxxvgio.hpp  gnarray.hpp  gwham.hpp  hamiltonian.hpp  typedefs.hpp mc.hpp gmbar.hpp trjsubtrj.hpp
all: gwham_gromacs463_umb gwham_gromacs463_umbfb mc gmbar_gromacs463_umb 

ran2.o: ran2.c
	$(CC) $(CFLAGS) -c ran2.c

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

fileio.o: fileio.hpp fileio.cpp
	$(CC) $(CFLAGS) -c fileio.cpp

gmxmdpio.o: fileio.o gmxmdpio.hpp gmxmdpio.cpp 
	$(CC) $(CFLAGS) -c gmxmdpio.cpp 

gwham_gromacs463_umb : $(OBJSGMX463UMB) 
	$(CC) $(CFLAGS) -o gwham_gromacs463_umb $(OBJSGMX463UMB) $(LFLAGS) 

gwham_gromacs463_umbfb : $(OBJSGMX463UMBFB) 
	$(CC) $(CFLAGS) -o gwham_gromacs463_umbfb $(OBJSGMX463UMBFB) $(LFLAGS) 

gmbar_gromacs463_umb : $(OBJSMBARGMX463UMB) 
	$(CC) $(CFLAGS) -o gmbar_gromacs463_umb $(OBJSMBARGMX463UMB) $(LFLAGS) 

mc : $(OBJSMC) 
	$(CC) $(CFLAGS) -o mc $(OBJSMC) $(LFLAGS)

#driver: bootstrap.o driver.o $(NR)/ran2.o $(NR)/locate.o
#	$(CC) $(CFLAGS) -o driver bootstrap.o driver.o ran2.o locate.o $(LFLAGS)

depend:
	makedepend -- $(CFLAGS) -- $(SRCS)

clean:
	rm -f gwham_gromacs463_umbfb gwham_gromacs463_umb gmbar_gromacs463_umbfb mc *.o .*.swp *~

remake:
	make clean
	make depend
	make all

tags: $(SRCS) $(HEADER_FILES) 
	ctags *

# DO NOT DELETE

ensemble.o: typedefs.hpp /usr/include/math.h /usr/include/features.h
ensemble.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
ensemble.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
ensemble.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
ensemble.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
ensemble.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
ensemble.o: /usr/include/bits/mathcalls.h ensemble.hpp /usr/include/stdlib.h
ensemble.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
ensemble.o: /usr/include/endian.h /usr/include/bits/endian.h
ensemble.o: /usr/include/bits/byteswap.h /usr/include/sys/types.h
ensemble.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
ensemble.o: /usr/include/time.h /usr/include/sys/select.h
ensemble.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
ensemble.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
ensemble.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
ensemble.o: /usr/include/stdio.h /usr/include/libio.h
ensemble.o: /usr/include/_G_config.h /usr/include/wchar.h
ensemble.o: /usr/include/bits/wchar.h /usr/include/xlocale.h
ensemble.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
gmxmdpio.o: typedefs.hpp /usr/include/math.h /usr/include/features.h
gmxmdpio.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
gmxmdpio.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
gmxmdpio.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
gmxmdpio.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
gmxmdpio.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
gmxmdpio.o: /usr/include/bits/mathcalls.h gmxmdpio.hpp gmxpullpot.hpp
gmxmdpio.o: hamiltonian.hpp ensemble.hpp fileio.hpp fileio_utils.hpp
gmxmdpio.o: /usr/include/string.h /usr/include/xlocale.h
gmxmdpio.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
gmxmdpio.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
gmxmdpio.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
gmxmdpio.o: /usr/include/sys/types.h /usr/include/bits/types.h
gmxmdpio.o: /usr/include/bits/typesizes.h /usr/include/time.h
gmxmdpio.o: /usr/include/sys/select.h /usr/include/bits/select.h
gmxmdpio.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
gmxmdpio.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
gmxmdpio.o: /usr/include/alloca.h /usr/include/stdio.h /usr/include/libio.h
gmxmdpio.o: /usr/include/_G_config.h /usr/include/wchar.h
gmxmdpio.o: /usr/include/bits/wchar.h /usr/include/bits/stdio_lim.h
gmxmdpio.o: /usr/include/bits/sys_errlist.h
fileio.o: typedefs.hpp /usr/include/math.h /usr/include/features.h
fileio.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
fileio.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
fileio.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
fileio.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
fileio.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
fileio.o: /usr/include/bits/mathcalls.h /usr/include/stdlib.h
fileio.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
fileio.o: /usr/include/endian.h /usr/include/bits/endian.h
fileio.o: /usr/include/bits/byteswap.h /usr/include/sys/types.h
fileio.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
fileio.o: /usr/include/time.h /usr/include/sys/select.h
fileio.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
fileio.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
fileio.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
fileio.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
fileio.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
fileio.o: /usr/include/bits/posix2_lim.h fileio.hpp fileio_utils.hpp
fileio.o: /usr/include/string.h /usr/include/xlocale.h
gmxpullpot.o: gmxpullpot.hpp typedefs.hpp /usr/include/math.h
gmxpullpot.o: /usr/include/features.h /usr/include/bits/predefs.h
gmxpullpot.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
gmxpullpot.o: /usr/include/gnu/stubs.h /usr/include/bits/huge_val.h
gmxpullpot.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
gmxpullpot.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
gmxpullpot.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
gwham_gromacs463_umb.o: typedefs.hpp /usr/include/math.h
gwham_gromacs463_umb.o: /usr/include/features.h /usr/include/bits/predefs.h
gwham_gromacs463_umb.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
gwham_gromacs463_umb.o: /usr/include/gnu/stubs.h /usr/include/bits/huge_val.h
gwham_gromacs463_umb.o: /usr/include/bits/huge_valf.h
gwham_gromacs463_umb.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
gwham_gromacs463_umb.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
gwham_gromacs463_umb.o: /usr/include/bits/mathcalls.h hamiltonian.hpp
gwham_gromacs463_umb.o: ensemble.hpp gwham.hpp densityofstate.hpp
gwham_gromacs463_umb.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
gwham_gromacs463_umb.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
gwham_gromacs463_umb.o: /usr/include/bits/endian.h
gwham_gromacs463_umb.o: /usr/include/bits/byteswap.h /usr/include/sys/types.h
gwham_gromacs463_umb.o: /usr/include/bits/types.h
gwham_gromacs463_umb.o: /usr/include/bits/typesizes.h /usr/include/time.h
gwham_gromacs463_umb.o: /usr/include/sys/select.h /usr/include/bits/select.h
gwham_gromacs463_umb.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
gwham_gromacs463_umb.o: /usr/include/sys/sysmacros.h
gwham_gromacs463_umb.o: /usr/include/bits/pthreadtypes.h
gwham_gromacs463_umb.o: /usr/include/alloca.h /usr/include/stdio.h
gwham_gromacs463_umb.o: /usr/include/libio.h /usr/include/_G_config.h
gwham_gromacs463_umb.o: /usr/include/wchar.h /usr/include/bits/wchar.h
gwham_gromacs463_umb.o: /usr/include/xlocale.h /usr/include/bits/stdio_lim.h
gwham_gromacs463_umb.o: /usr/include/bits/sys_errlist.h gnarray.hpp
gwham_gromacs463_umb.o: fileio.hpp fileio_utils.hpp /usr/include/string.h
gwham_gromacs463_umb.o: gmxxvgio.hpp gmxmdpio.hpp gmxpullpot.hpp
mc.o: ensemble.hpp typedefs.hpp /usr/include/math.h /usr/include/features.h
mc.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
mc.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
mc.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
mc.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
mc.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
mc.o: /usr/include/bits/mathcalls.h gnarray.hpp /usr/include/stdio.h
mc.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
mc.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
mc.o: /usr/include/bits/wchar.h /usr/include/xlocale.h
mc.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h gwham.hpp
mc.o: densityofstate.hpp hamiltonian.hpp /usr/include/stdlib.h
mc.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
mc.o: /usr/include/endian.h /usr/include/bits/endian.h
mc.o: /usr/include/bits/byteswap.h /usr/include/sys/types.h
mc.o: /usr/include/time.h /usr/include/sys/select.h
mc.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
mc.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
mc.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h mc.hpp ran2.c
mc.o: fileio_utils.hpp /usr/include/string.h
gmbar.o: typedefs.hpp /usr/include/math.h /usr/include/features.h
gmbar.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
gmbar.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
gmbar.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
gmbar.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
gmbar.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
gmbar.o: /usr/include/bits/mathcalls.h gmbar.hpp /usr/include/stdio.h
gmbar.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
gmbar.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
gmbar.o: /usr/include/bits/wchar.h /usr/include/xlocale.h
gmbar.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
gmbar.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
gmbar.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
gmbar.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
gmbar.o: /usr/include/sys/types.h /usr/include/time.h
gmbar.o: /usr/include/sys/select.h /usr/include/bits/select.h
gmbar.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
gmbar.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
gmbar.o: /usr/include/alloca.h
gmbar_gromacs463_umb.o: typedefs.hpp /usr/include/math.h
gmbar_gromacs463_umb.o: /usr/include/features.h /usr/include/bits/predefs.h
gmbar_gromacs463_umb.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
gmbar_gromacs463_umb.o: /usr/include/gnu/stubs.h /usr/include/bits/huge_val.h
gmbar_gromacs463_umb.o: /usr/include/bits/huge_valf.h
gmbar_gromacs463_umb.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
gmbar_gromacs463_umb.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
gmbar_gromacs463_umb.o: /usr/include/bits/mathcalls.h hamiltonian.hpp
gmbar_gromacs463_umb.o: ensemble.hpp gwham.hpp densityofstate.hpp
gmbar_gromacs463_umb.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
gmbar_gromacs463_umb.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
gmbar_gromacs463_umb.o: /usr/include/bits/endian.h
gmbar_gromacs463_umb.o: /usr/include/bits/byteswap.h /usr/include/sys/types.h
gmbar_gromacs463_umb.o: /usr/include/bits/types.h
gmbar_gromacs463_umb.o: /usr/include/bits/typesizes.h /usr/include/time.h
gmbar_gromacs463_umb.o: /usr/include/sys/select.h /usr/include/bits/select.h
gmbar_gromacs463_umb.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
gmbar_gromacs463_umb.o: /usr/include/sys/sysmacros.h
gmbar_gromacs463_umb.o: /usr/include/bits/pthreadtypes.h
gmbar_gromacs463_umb.o: /usr/include/alloca.h /usr/include/stdio.h
gmbar_gromacs463_umb.o: /usr/include/libio.h /usr/include/_G_config.h
gmbar_gromacs463_umb.o: /usr/include/wchar.h /usr/include/bits/wchar.h
gmbar_gromacs463_umb.o: /usr/include/xlocale.h /usr/include/bits/stdio_lim.h
gmbar_gromacs463_umb.o: /usr/include/bits/sys_errlist.h gmbar.hpp gnarray.hpp
gmbar_gromacs463_umb.o: fileio.hpp fileio_utils.hpp /usr/include/string.h
gmbar_gromacs463_umb.o: gmxxvgio.hpp gmxmdpio.hpp gmxpullpot.hpp
