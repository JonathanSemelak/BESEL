# Default NETCDF directory. Override with `make NETCDF_DIR=/your/dir`
NETCDF_DIR := /usr

NETCDFALL := -I$(NETCDF_DIR)/include -L$(NETCDF_DIR)/lib -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -lnetcdf -lnetcdf

FC := gfortran

SRC := readandget.f90 writeall.f90 opt.f90 statistics.f90 neb.f90 profiles.f90 fire.f90
OBJ := $(SRC:.f90=.o)

.PHONY: all clean

all: feneb bandbuilder extractor getdihe getnrmsd integrator segments getmaxgrad

feneb: $(OBJ) feneb.f90
	$(FC) $(OBJ) feneb.f90 -o $@ $(NETCDFALL)

bandbuilder: readandget.o opt.o statistics.o neb.o writeall.o profiles.o bandbuilder.f90
	$(FC) readandget.o opt.o statistics.o neb.o writeall.o profiles.o bandbuilder.f90 -o $@ $(NETCDFALL)

extractor: readandget.o extractcoordtofile.f90
	$(FC) readandget.o extractcoordtofile.f90 -o $@ $(NETCDFALL)

getdihe: readandget.o getdihe.f90
	$(FC) readandget.o getdihe.f90 -o $@ $(NETCDFALL)

getnrmsd: readandget.o getnrmsd.f90
	$(FC) readandget.o getnrmsd.f90 -o $@ $(NETCDFALL)

integrator: readandget.o profiles.o integrator.f90
	$(FC) readandget.o profiles.o integrator.f90 -o $@ $(NETCDFALL)

segments: readandget.o freenergysegments.f90
	$(FC) readandget.o freenergysegments.f90 -o $@ $(NETCDFALL)

getmaxgrad: getmaxgrad.f90
	$(FC) -o getmaxgrad getmaxgrad.f90

%.o: %.f90
	$(FC) -c $< -o $@ $(NETCDFALL)

clean:
	rm -f *.o feneb bandbuilder extractor getdihe getnrmsd integrator segments getmaxgrad

