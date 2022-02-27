#---------------------------------------------------------------------------
#Code below searches for cplex install dir. If you know it, you can erase all
#lines below (up to first comment) and just put a line as follows.
#CPLEX=your/cplex/folder/
#Otherwise let it check for several folders and it will set the first existing
#one. If all fail, the compilation will halt before trying to build main
ifneq ($(wildcard /data/cplex126/),)  #check if this folder exists
    CPLEX   = /data/cplex126/
    CPLEXVER="cplex126"
endif 
ifneq ($(wildcard /opt/ibm/ILOG/),)            
    #Cplex is often installed in a subdir of above folder.
    #take last subdir inside (ex, CPLEX_Studio...) as it can be the last version
    CPLEXVER= $(shell ls /opt/ibm/ILOG/|tail -1)
    CPLEX   = /opt/ibm/ILOG/${CPLEXVER}
endif 
#If possible I always take version 12.10
ifneq ($(wildcard /opt/ibm/ILOG/CPLEX_Studio1210/),) 
    CPLEX   = /opt/ibm/ILOG/CPLEX_Studio1210/
    CPLEXVER="CPLEX_Studio126"
endif 
CPLEXVERNO_ALLDIGITS=$(shell echo ${CPLEXVER}|sed "s/[^0-9]//g")  #remove non digits
CPLEXVERNO=$(shell echo ${CPLEXVERNO_ALLDIGITS}|cut -c 1-3)       #take first 3 letters (eg. 126)

#----------------------------------------------------------------------------
#The rest depends on the cplex install dire defined above

CPLEXDIR      = $(CPLEX)/cplex/
CONCERTDIR    = $(CPLEX)/concert/

# ---------------------------------------------------------------------
# Compiler selection, code optimization, debug and warning options
# ---------------------------------------------------------------------
#use DNDEBUG to remove asserts and turn off some clog debug messages
CCC = g++ 
CCFLAGS =  -O3 -flto -m64 -Wall -Wextra -fomit-frame-pointer -funroll-loops -Wno-ignored-attributes -g 
#maybe add  -fomit-frame-pointer -funroll-loops 
#maybe add -Wextra and -Wpedantic
#you could add -flto for link-time optimizaton (e.g., in-lining function 
#from a .o to another .o). However, the speed-up can be about 1 or 2% or less
#maybe -Ofast that enables -ffast-math. 
#For me, the following can be good enough instead of -Ofast. However, it
#assumes things like floating point multiplication is associative, untrue.
#-funsafe-math-optimizations -ffinite-math-only
#also -Wshadow


# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------
LIBFORMAT  = static_pic
#cplex related
#use SYSTEM=x86-64_sles10_4.1 for cplex versions <12.6
#SYSTEM     = x86-64_linux
SYSTEM=$(shell ls $(CPLEXDIR)/lib/|grep x86)

CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CCLNFLAGS =" " 
CCLNFLAGSCPLEX = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lrt -lpthread -ldl
               #-m32 -pthread -lm (add more when needed)

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CCFLAGSCPLEX = $(CCFLAGS) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)  -DIL_STD #DIL_STD: CPLEX specific macro

# ---------------------------------------------------------------------
# Main compilation and linking commands
# ---------------------------------------------------------------------
all:cplexdirexists test
cplexdirexists:
	if [ -d "$(CPLEX)" ]; then echo -e "\nCplex dir found. I'll compile main if not already compiled."; else echo -e "\n\nCplex dir not found!!!!! \nPlease set add a line like below on \nCPLEX=/path/to/your/cplex/install/folder/\nat the beginning of this Makefile"; exit 1; fi
test:CuttingPlanesEngine.o src/test.cpp 
	$(CCC) $(CCFLAGSCPLEX) CuttingPlanesEngine.o src/test.cpp -o test $(CCLNFLAGSCPLEX) 
CuttingPlanesEngine.o : src/CuttingPlanesEngine.cpp src/CuttingPlanesEngine.h
	$(CCC) -c $(CCFLAGSCPLEX) src/CuttingPlanesEngine.cpp -o CuttingPlanesEngine.o 

clean: cleanobj
	rm -f main *.lp *.log *.err *.cout *.o
