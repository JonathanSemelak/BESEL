#gfortran feneb.f90 get_all.f90 write_all.f90 read_input.f90 opt.f90 neb.f90 -o feneb `/usr/bin/nf-config --fflags --flibs`

NETCDF_DIR=/home/fpedron/prog/Amber18

NETCDFALL=$(echo $NETCDF_DIR/include -L$NETCDF_DIR/lib -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -lnetcdf -lnetcdf)

#echo $NETCDFALL

gfortran -c readandget.f90 -I $NETCDFALL
gfortran -c opt.f90
gfortran -c neb.f90
gfortran -c writeall.f90
#gfortran -c feneb.f90

gfortran readandget.o writeall.o opt.o neb.o feneb.f90 -o feneb -I $NETCDFALL


#gfortran readandget.f90 write_all.f90 opt.f90 neb.f90 feneb.f90 -o feneb -I$NETCDF_DIR/include -L$NETCDF_DIR/lib -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -lnetcdf -lnetcdf


#gfortran get_all.f90 write_all.f90 read_input.f90 opt.f90 neb.f90 feneb.f90 -o feneb `/usr/bin/nf-config --fflags --flibs`
#compile_feneb.sh  feneb.f90  get_all.f90  opt.f90         write_all.f90
#feneb             feneb.in   neb.f90      read_input.f90
#gfortran feneb.f90 -o feneb `/usr/bin/nf-config --fflags --flibs`
