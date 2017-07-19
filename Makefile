FC     	= gfortran
LC     	= $(FC)
EXE    	= h0cf3
FITSDIR = /home/wcardona/lib
LIBFITS = cfitsio
INCDIR	= /home/wcardona/software/Healpix/Healpix_3.31/include
IDIR	= #/home/wilmar/additional-software/Healpix_3.00/include
INFGSL  = /home/wcardona/include/fgsl
LIBDIR	= /home/wcardona/lib	
LDIR	= /home/wcardona/software/Healpix/Healpix_3.31/lib
#F_FL   	= -O3 -I$(INCDIR) -I$(IDIR) -DGFORTRAN -fno-second-underscore -fopenmp -fPIC -g
F_FL   	= -O3 -Wall -I$(INCDIR) -I$(INFGSL) -DGFORTRAN -fno-second-underscore -fopenmp -fPIC -g 
LIB_FL 	= -L$(LDIR) -L$(FITSDIR) -lhealpix -lhpxgif -lhpxgif -l$(LIBFITS) #-L$(LIBDIR) #-lfgsl #-llapack -lblas  -lranlib -lrnglib    -Wl,
#####################
OBJ   =  arrays.o fiducial.o functions.o h0cf3.o

def:	$(OBJ) $(OBJNR) $(OBJODE)
	$(LC) $(F_FL) $(OBJ) $(OBJNR) $(OBJODE) -o $(EXE)  $(LIB_FL)

%.o:	%.f90
	$(FC) $(F_FL) -c $<

%.o:	%.F90
	$(FC) $(F_FL) -c $<

%.o:	%.f
	$(FC) $(F_FL) -c $<

clean :
	rm -f *.o *.mod *.ini *~  fort.* *.out $(EXE)

### put dependencies here ###

h0cf3.o :	arrays.o functions.o fiducial.o
