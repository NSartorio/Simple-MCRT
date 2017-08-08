#*************** Makefile created by Mike Wolff ****************************

#******************************** G77/Linux Fortran ************************
FC     =       mpif90
#EXTRA_OPT =     -mpentium -malign-double -fforce-mem -fforce-addr \
#                -ffast-math -funroll-all-loops
EXTRA_OPT =     -malign-double -fforce-addr \
                -ffast-math -funroll-all-loops
# May want to experiment by adding the extra optimization flags to get
## better runtime. But then again, maybe not.
#FFLAGS  =       -O2 $(EXTRA_OPT) -ffloat-store
FFLAGS  =       -O2
LDFLAGS = 
time_it         = get_cpu_sun

#******************************** PGI Fortran ************************
#FC      =       pgf77
#FFLAGS  =      -fast
#LDFLAGS =	-fast 
#time_it         = get_cpu_sun

#******************************** Sun Fortran ************************
#FC     =       f77
#FFLAGS  =      -fast -O
#LDFLAGS =	-fast -O
#time_it         = get_cpu_sun

#******************************** Lahey-Fujitsu lf95 ************************
#
#FC      =       lf95
#FFLAGS  =       --tpp --nsav -O --nwarn -c
#LDFLAGS =
#time_it         = get_cpu_sun

#****************************************************************************

#sources = emit.f90 gridset.f90 ionize.f90 ran2.f90 sources.f90 \
#          sources.f90 tauint.f90      

OBJSB = mcemit.f90 mcgridset.f90 mcionize.f90 ran2.f90 mcsources.f90 \
          mctauint.f90 mcion.f90
          
#files   = $(sources) makefile make.default CHANGES.txt LICENSE.txt README.txt  \
#          hosts/default problems


mcgrid:	$(OBJSB)
		$(FC) $(OBJSB) $(LDFLAGS) -o mcion

clean:;		/bin/rm -f90 *.o

