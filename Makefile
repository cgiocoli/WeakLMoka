#####################################################################
#                                                                   #
#               Make file for MOKA                                  #
#                                                                   #
#                         cgiocoli@gmail.com                        #
#####################################################################

# executable name
PROG = $(HOME)/bin/WL_MOKA_1.4

MAIN = wlmain.cpp power2D.cpp cMZhao.cpp

# .cpp lib internal in MOKA
SOURCES = ../Moka/utilities.cpp \
          ../Moka/cosmology.cpp \
          ../Moka/halo.cpp \
          ../Moka/nfwHalo.cpp \
          ../Moka/nfwLens.cpp \
          ../Moka/hernq_Halo.cpp \
          ../Moka/hernq_Lens.cpp \
          ../Moka/jaffe_Halo.cpp \
          ../Moka/jaffe_Lens.cpp \
          ../Moka/sisHalo.cpp \
          ../Moka/sisLens.cpp 

# gsl, cfitsio, CCfits, fftw
LIBS = -L/home/cgiocoli/lib/gsl-1.13/lib  -lgslcblas -lgsl \
       -L/home/cgiocoli/lib/cfitsio/lib \
       -L/home/cgiocoli/lib/CCfits/lib -lCCfits -lcfitsio \
       -L/home/cgiocoli/lib/fftw-3.2.2/lib -lfftw3 -lm 

# gsl, cfitsio, CCfits, fftw  
ALLFLAGS = -I/home/cgiocoli/lib/gsl-1.13/include/gsl \
	   -I/home/cgiocoli/lib/gsl-1.13/include \
           -I/home/cgiocoli/lib/cfitsio/include \
           -I/home/cgiocoli/lib/CCfits/include \
           -I/home/cgiocoli/lib/fftw-3.2.2/include

# 
DEBUG = -g 
# compiler -Wall # use this to optimize -O2 
#CC = gcc -stdlib=libstdc++ -std=c++11 
CC = /usr/bin/gcc -std=c++0x
#
RM = rm -f -r
#
OBJ = $(SOURCES:.cpp=.o)
#

CFLAGS=$(ALLFLAGS) $(DEBUG)
#
CLEAR = clear

default: moka
moka: 
	$(CC) -c $(SOURCES) $(CFLAGS)
	ar r libmoka.a *.o
	$(CC) $(MAIN) -L. -lmoka $(CFLAGS) -o $(PROG) $(LIBS) 

main: 
	$(CC) $(MAIN) -L. -lmoka $(CFLAGS) -o $(PROG) $(LIBS) 

clean:
	$(RM) $(PROG) $(OBJ) *.o *.a

lib:
	$(CC) -c $(SOURCES) $(CFLAGS)
	ar r libmoka.a $(OBJ)

clall:
	$(RM) $(PROG) $(OBJ) libmoka.a *~ *.o

clall0:
	$(RM) $(PROG) $(OBJ) libmoka.a *~ html *.o

new: clean default

help:
	$(CLEAR)	
	@echo 
	@echo make - compile extlib, lib and main linking lib
	@echo make main - compile main linking lib
	@echo make lib - compile only internal lib
	@echo make clean - rm .o and executable
	@echo make clall - rm .o .a *.~ and executable
	@echo make clall0 - rm .o .a *.~ executable and html doxygen
	@echo 
	@echo - - - Moka written by Carlo Giocoli	
	@echo "                   " Matthias Bartelmann
	@echo "                   " Massimo Meneghetti
	@echo "          " please contact them for more
	@echo "          " detalis: carlo.giocoli@unibo.it
	@echo "                   " bartelmann@uni-heidelberg.de
	@echo "                   " massimo.meneghetti@oabo.inaf.it
	@echo "                   " margarita.petkova@unibo.it
	@echo
	@echo  website: http://cgiocoli.wordpress.com/research-interests/moka
	@echo
