# Copyright (C) 2009 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

##########################################################################
#    You can modify this example makefile to fit for your own program.   #
#    Usually, you only need to change the four CHANGEME entries below.   #
##########################################################################

# CHANGEME: This should be the name of your executable
EXE = TutorialCpp

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS = TutorialCpp_main.o \
	TutorialCpp_nlp.o

# CHANGEME: Additional libraries
ADDLIBS =

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS =

##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile Ipopt.   #
##########################################################################


# C++ Compiler command
CXX = g++

# C++ Compiler options
CXXFLAGS = -O2 -DNDEBUG 

# additional C++ Compiler options for linking
CXXLINKFLAGS =  -Wl,--rpath -Wl,/home/hdang/Schreibtisch/HiWiSport/coinbrew-master/dist/lib

prefix=/home/hdang/Schreibtisch/HiWiSport/coinbrew-master/dist
exec_prefix=${prefix}

# Include directories
INCL = `PKG_CONFIG_PATH=/home/hdang/Schreibtisch/HiWiSport/coinbrew-master/dist/lib/pkgconfig: pkg-config --cflags ipopt` $(ADDINCFLAGS)
#INCL = -I${prefix}/include/coin -I/home/hdang/Schreibtisch/HiWiSport/coinbrew-master/dist/include/coin-or/mumps -I/home/hdang/Schreibtisch/HiWiSport/coinbrew-master/dist/include/coin-or/metis -I/home/jhogg/src/coinhsl/bin-linux_x86_64/coinhsl-linux-x86_64-2014.01.10/include     $(ADDINCFLAGS)

# Linker flags
LIBS = `PKG_CONFIG_PATH=/home/hdang/Schreibtisch/HiWiSport/coinbrew-master/dist/lib/pkgconfig: pkg-config --libs ipopt`
#LIBS = -L${exec_prefix}/lib -lipopt -L/home/hdang/Schreibtisch/HiWiSport/coinbrew-master/dist/lib -L/home/jhogg/src/coinhsl/bin-linux_x86_64/coinhsl-linux-x86_64-2014.01.10/lib -lcoinhsl -lcoinmumps -llapack -lblas    -lm  -ldl

all: $(EXE)

.SUFFIXES: .cpp .o

$(EXE): $(OBJS)
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $(OBJS) $(ADDLIBS) $(LIBS)

clean:
	rm -rf $(EXE) $(OBJS) ipopt.out

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ $<
