#NETCDF_DIR=/home/fpedron/prog/Amber18

NETCDFALL=$(echo /usr/bin/nf-config --fflags --flibs)

#echo $NETCDFALL

gfortran -c readandget.f90 `$NETCDFALL`
gfortran -c opt.f90
gfortran -c neb.f90
gfortran -c writeall.f90
gfortran -c profiles.f90
gfortran -c statistics.f90

gfortran readandget.o writeall.o opt.o neb.o profiles.o statistics.o feneb.f90 -o feneb `$NETCDFALL`

gfortran readandget.o opt.o neb.o writeall.o profiles.o bandbuilder.f90 -o bandbuilder `$NETCDFALL`

gfortran readandget.o extractcoordtofile.f90 -o extractor `$NETCDFALL`

gfortran getdihe.f90 -o getdihe
