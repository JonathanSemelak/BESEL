# Default NETCDF directory. Override with `make NETCDF_DIR=/your/dir`
NETCDF_DIR := /usr

NETCDFALL := -I$(NETCDF_DIR)/include -L$(NETCDF_DIR)/lib -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -lnetcdf -lnetcdf

FC := gfortran

SRC := readandget.f90 writeall.f90 opt.f90 statistics.f90 neb.f90 profiles.f90 fire.f90
OBJ := $(SRC:.f90=.o)
BIN_DIR := bin

.PHONY: all clean

all: feneb bandbuilder extractor getdihe getnrmsd integrator segments getmaxgrad
	@echo "FENEB installation is complete! 🎉"
	@echo "Don't forget to export the bin path: export PATH=\$$PATH:$(PWD)/$(BIN_DIR)"

feneb: $(OBJ) feneb.f90
	@echo "Compiling feneb..."
	$(FC) $(OBJ) feneb.f90 -o $(BIN_DIR)/$@ $(NETCDFALL)

bandbuilder: readandget.o opt.o statistics.o neb.o writeall.o profiles.o bandbuilder.f90
	@echo "Compiling bandbuilder..."
	$(FC) readandget.o opt.o statistics.o neb.o writeall.o profiles.o bandbuilder.f90 -o $(BIN_DIR)/$@ $(NETCDFALL)

extractor: readandget.o extractcoordtofile.f90
	@echo "Compiling extractor..."
	$(FC) readandget.o extractcoordtofile.f90 -o $(BIN_DIR)/$@ $(NETCDFALL)

getdihe: readandget.o getdihe.f90
	@echo "Compiling getdihe..."
	$(FC) readandget.o getdihe.f90 -o $(BIN_DIR)/$@ $(NETCDFALL)

getnrmsd: readandget.o getnrmsd.f90
	@echo "Compiling getnrmsd..."
	$(FC) readandget.o getnrmsd.f90 -o $(BIN_DIR)/$@ $(NETCDFALL)

integrator: readandget.o profiles.o integrator.f90
	@echo "Compiling integrator..."
	$(FC) readandget.o profiles.o integrator.f90 -o $(BIN_DIR)/$@ $(NETCDFALL)

segments: readandget.o freenergysegments.f90
	@echo "Compiling segments..."
	$(FC) readandget.o freenergysegments.f90 -o $(BIN_DIR)/$@ $(NETCDFALL)

getmaxgrad: getmaxgrad.f90
	@echo "Compiling getmaxgrad..."
	$(FC) -o $(BIN_DIR)/getmaxgrad getmaxgrad.f90

%.o: %.f90
	@echo "Compiling $<..."
	$(FC) -c $< -o $@ $(NETCDFALL)

clean:
	@echo "Cleaning up..."
	rm -f *.o $(BIN_DIR)/*

