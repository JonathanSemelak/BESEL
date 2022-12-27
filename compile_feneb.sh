NETCDFALL=$(echo /usr/bin/nf-config --fflags --flibs)

gfortran -c readandget.f90 `$NETCDFALL`
gfortran -c opt.f90
gfortran -c neb.f90
gfortran -c writeall.f90
gfortran -c profiles.f90
gfortran -c statistics.f90
gfortran -c fire.f90

gfortran readandget.o writeall.o opt.o statistics.o neb.o profiles.o fire.o feneb.f90 -o feneb `$NETCDFALL`

gfortran readandget.o opt.o statistics.o neb.o writeall.o profiles.o bandbuilder.f90 -o bandbuilder `$NETCDFALL`

gfortran readandget.o extractcoordtofile.f90 -o extractor `$NETCDFALL`

gfortran readandget.o getdihe.f90 -o getdihe `$NETCDFALL`

gfortran readandget.o getnrmsd.f90 -o getnrmsd `$NETCDFALL`

gfortran readandget.o profiles.o integrator.f90 -o integrator `$NETCDFALL`

gfortran readandget.o freenergysegments.f90 -o segments `$NETCDFALL`

