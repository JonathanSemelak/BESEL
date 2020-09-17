gfortran feneb.f90 get_all.f90 write_all.f90 read_input.f90 opt.f90 neb.f90 -o feneb `/usr/bin/nf-config --fflags --flibs`

#gfortran get_all.f90 write_all.f90 read_input.f90 opt.f90 neb.f90 feneb.f90 -o feneb `/usr/bin/nf-config --fflags --flibs`
#compile_feneb.sh  feneb.f90  get_all.f90  opt.f90         write_all.f90
#feneb             feneb.in   neb.f90      read_input.f90
#gfortran feneb.f90 -o feneb `/usr/bin/nf-config --fflags --flibs`
