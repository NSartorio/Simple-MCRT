subroutine mcsources(xsource,ysource,zsource,lsource,lumtot)

implicit none

      include 'simplegrid.txt'

 integer, parameter :: nsource=1
 real*8 xsource(nsource),ysource(nsource),zsource(nsource)
 real*8 lsource(nsource),lumtot
 integer i


!**** Set photon locations and luminosities

 xsource(1)=0.
 ysource(1)=0.
 zsource(1)=0.
 lsource(1)=1.


!**** Calculate total luminosity of all sources
 lumtot=0.
 do i=1,nsource
   lumtot=lumtot+lsource(i)
 end do

 return

end
