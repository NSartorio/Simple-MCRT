program mcion

implicit none
!(ntot,time,pres)

      include 'simplegrid.txt'
      include 'photon.txt'    

!***** Parameter declarations ****************************************
      integer nphotons, niter,nph,iseed,is,jcount,nscatt,totscatt,xl,yl
      integer xcell,ycell,zcell,tflag,aflag,i,j,k,mubin,iter
      integer pesc,nphot
      integer nsource   !&&&&
      integer nbub,d_absorb
      integer ierr, irank, isize,l
      real*8 kappa,albedo,hgg,pl,pc,sc,xmax,ymax,zmax,rimage,tema,temb
      real*8 pi,twopi,fourpi,g2,ximage,yimage,ph,wgt1,wgt,tau,tau1,imtot
      real*8 sigma,jfac,Q49,sigr0,sigd,alphaa,length,vcell
      real*8 zstars,sigdustR,pionize,pdust, Q490, time
      integer nphot_tot
      real*8 unitdens, unitmassdens,unitElow,UnitEhigh
      real*8 conversion, factor, rhosum
      real*8 xsource(1),ysource(1),zsource(1)   !&&&&
      real*8 lsource(1),rbub(1),lumtot   !&&&&

      real  ran2
!**************** BASIC PARAMETERS ***********************
 namelist /params/ nphotons,niter, nsource,tema,temb,iseed,pionize,&
     & xmax,ymax,zmax,Q49,length,unitdens,sigd,alphaa

 open(10,file='parameters.par',status='unknown')
   read(10,params)
 close(10)

! number of photons
!      nphotons=100000
!number of iterations
!      niter=10
!number of sources
!      nsource = 1   
!temperature of neutral gas
!      tema=500
!temperature of ionized gas
!      temb=8000.
!seed for ran2 random number generator
!      iseed=76865
      iseed=-abs(iseed)  ! must be negative for ran2
!probability of emission of an ionizing
!photon after absorption by an H atom
!      pionize=0.38
!size of the box in x, y and z axis
!      xmax=5.
!      ymax=5.
!      zmax=5.
!luminosity (in 10^49 photons p/ second)
!      Q49=4.26
!unit length = 1pc (in cm)
!      length=3.086e18
!unit density
!      unitdens=3.1d3 !in number of hydrogen atoms
      unitmassdens=unitdens*1.67d-24 !mass of H in g
!cross section for diffuse photons
!      sigd=6.3d-18 ! in cm^2 
      sigr0=sigd*length !cross section*length unit (cm^3)
!recombination coefficient
!      alphaa=5.25e-13 ! in /cm^3/s  for T=8000K
!2*pi  
      twopi=8.*atan(1.)
!cell volume
      vcell = (2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)   
!factor that multiplies the intensity for ionization fraction calc.
      factor = (1.d49*Q49*sigd)/(vcell*alphaa*length*length)

!******************	HEADER ******************************

print*, "Running MC Ionizing code with:"
print*, "number of iterations:", niter
print*, "number of photons:", nphotons 
print*, "luminosity:", Q49
print*, "simulation box of",xmax ,"pc"

!**************  Loop through x, y, and z to set up grid density
    do i=1,nxg
      do j=1,nyg
        do k=1,nzg
           open(unit=7, file="density1.dat")
           read(7,*), ntot(i,j,k)
           ntot(i,j,k) = ntot(i,j,k)*500
           rhokap(i,j,k)=0.
           tem(i,j,k) = 0.
        end do
      end do
    end do
    close(7)

!_--_--_Set up density grid_--_--_--_--_--_--_--_--_--_--_--_--_--
      call mcgridset(xface,yface,zface,rhokap,rhosum,sigr0,length,&
                                &xmax,ymax,zmax,nfrac,ntot)

!_--_--_Set up point source locations, luminosities, total luminosity
      call mcsources(xsource,ysource,zsource,lsource,lumtot)

      !print*,"lumtot", lumtot

!--------------------START ITERATIONS--------------------!


do iter = 1, niter     !**********LOOP over iterations

!  if (iter.ge.7) then
!     nphotons = 10*nphotons
!  end if
     
 
!factor that multiplies the intensity for ionization fraction calc.
  jfac = factor/nphotons

print*, jfac

    do is=1,nsource !**********LOOP over sources

!calculate number of photons per source
          nph=int(nphotons*lsource(is)/lumtot)
          jcount = 0
          do j=1,nph !**********LOOP over photons from each source

!print @ which point we are
           jcount=jcount+1
           if(mod(jcount,1000).eq.0) then
          print *, 'iteration ',iter,',',jcount,' photons completed'
           end if
           
!***** Release photon from point source *******************************
!*************** Cartesian Grid 
          xp = xsource(is)
          yp = ysource(is)
          zp = zsource(is)
!*************** Linear Grid
          xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
          ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
          zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
!emit photon
          call mcemit(nxp,nyp,nzp,twopi,iseed)

          sigma=0.44 ! cross section for source photons: flux averaged
          nscatt=0   ! it hasn't been scattered yet
          tflag=0    !it is within the grid and it is ionizing

!******* Find first scattering location for random tau.
          call mctauint(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,&
     &xface,yface,zface,rhokap,ntot,jmean,&
     &xcell,ycell,zcell,sigma,tflag,iseed)

!--------------------------start do while---------------------------------!
          do while(tflag.eq.0)

!_______________________________________________
!tflag signalizes when a photon no longer      |
!need to be tracked, that is, if a photon has  |
!exited the grid or it was re-emitted as       |
!a non-ionizing photon                         |
!----------------------------------------------!

!check if it is going to be ionizing or not
              if(ran2(iseed).le.pionize) then

!photon absorbed and re-emitted as ionizing photon: change
!cross section and re-emit isotropically
                 call mcemit(nxp,nyp,nzp,twopi,iseed)

                 sigma=1.  !cross section for diffuse photons
              else
!********* photon absorbed and re-emitted as a non-ionizing photon:terminate 
                 !print*, "terminate: photon re-emitted as non-ionizing"
                 tflag = 1
                 goto 100 ! terminate
              end if

!************ Find next scattering location
              call mctauint(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,&
     &xface,yface,zface,rhokap,ntot,jmean,&
     &xcell,ycell,zcell,sigma,tflag,iseed)

           end do	
!--------------------------end do while---------------------------------!
100       continue
          end do !**********END LOOP over photons
        end do   !**********END LOOP over sources


!**************update the opacity (rhokap)************************

!_________________________________________________________________
!use the sum of path lengths to compute ionizing fraction        |
!and adjust opacity accordingly (ie. the fraction of H ionized   |
!does not contribute to the opacity).                            |
!this will be the opacity grid the photons of the next           |
!iteration will see.                                                 |
!----------------------------------------------------------------!

       call mcionize(jmean,rhokap,nfrac,ntot,jfac,sigr0) 

      end do !**********END LOOP over iterations

!******calculate the temperature and pressure
      tem(:,:,:) = tema
      ntot(:,:,:) = ntot(:,:,:)/unitdens 
      pres(:,:,:) = ntot(:,:,:)*(tema/tema)         
       do i = 1,nxg	
         do j = 1,nyg
            do k = 1,nzg
!            if (nfrac(i,j,k) .lt. 0.5) then 
      tem(i,j,k) = tema+temb*(1.-nfrac(i,j,k))
      pres(i,j,k) = (2.-nfrac(i,j,k))*pres(i,j,k)*(tem(i,j,k)/tema)
!            endif
             end do
          end do
       end do

!output ntot, temperature and ionized fraction -------------------
      open(unit=10,file='ntot@50.dat')
              do i=1,nxg 
	       write(10,*)(ntot(i,k,50),k=1,nyg)
              end do
      close(10)
      open(unit=10,file='tem@50.dat')
           do i=1,nxg 
	       write(10,*)(tem(i,k,50),k=1,nyg)
           end do
      close(10)

      open(unit=10,file='nfrac@50.dat')
              do i=1,nxg 
	       write(10,*)(nfrac(i,k,50),k=1,nyg)
              end do
      close(10)
     

end program

