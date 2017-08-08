subroutine mcemit(nxp,nyp,nzp,twopi,iseed)

implicit none
!uses the monte carlo method to select a direction
!of travel for the photon and sets the photon
!at the position of the source.

 integer iseed
 real*8 twopi
 real*8 nxp,nyp,nzp,sint,cost,sinp,cosp,phi
 real  ran2

!***** emit photon isotropically from point source location

!angle theta (measured from z axis) between 0 and pi
 cost=2.*ran2(iseed)-1        !value btw -1 and 1
 sint=sqrt(1.-cost*cost)      !value btw 0 and 1

!angle phi between 0 and 2pi
 phi=twopi*ran2(iseed)
 cosp=cos(phi)
 sinp=sin(phi)

!***** Set photon direction cosines for direction of travel 
 nxp=sint*cosp  
 nyp=sint*sinp
 nzp=cost

 return
end
