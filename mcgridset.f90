subroutine mcgridset(xface,yface,zface,rhokap,rhosum,sigr0,length,&
                                &xmax,ymax,zmax,nfrac,ntot)

implicit none

      include 'simplegrid.txt'

      real*8 xmax,ymax,zmax

      integer i,j,k
      real*8 x,y,z,rhosum,dV,taueq,taupole, sigr0, length

      print *, 'Setting up density grid....'

!**********  Linear Cartesian grid. Set up grid faces ****************

      do i=1,nxg+1
         xface(i)=(i-1)*2.*xmax/nxg
      end do
      do i=1,nyg+1
         yface(i)=(i-1)*2.*ymax/nyg
      end do
      do i=1,nzg+1
         zface(i)=(i-1)*2.*zmax/nzg
      end do

      rhosum=0.
!**************  Loop through x, y, and z to set up grid density.  
      do i=1,nxg
       do j=1,nyg
        do k=1,nzg
           x=xface(i)-xmax+xmax/nxg
           y=yface(j)-ymax+ymax/nyg
           z=zface(k)-zmax+zmax/nzg
           nfrac(i,j,k)=1.e-6
           rhokap(i,j,k)=ntot(i,j,k)*nfrac(i,j,k)*sigr0
           rhosum=rhosum+ntot(i,j,k)
        end do
       end do
      end do

      return
end
