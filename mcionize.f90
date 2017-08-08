subroutine mcionize(jmean,rhokap,nfrac,ntot,jfac,sigr0)

 implicit none



 include 'simplegrid.txt'

 integer i,j,k
 real*8 jfac,sigr0
 real*8 aa,bb,cc
 do i=1,nxg
  do j=1,nyg
   do k=1,nzg

 if((jmean(i,j,k).gt.0.).and.(ntot(i,j,k).gt.0.)) then
    aa=jmean(i,j,k)/ntot(i,j,k)/2.*jfac
    bb=2./aa
   if(bb.le.0.) then
    print *,'IONIZE: CRASH'
    stop
   endif

   cc=sqrt(bb+1.)

   nfrac(i,j,k)=1.+aa*(1.-cc)

 else

  nfrac(i,j,k)=1.  

 endif

 rhokap(i,j,k)=ntot(i,j,k)*nfrac(i,j,k)*sigr0 
 jmean(i,j,k) =0

   end do
  end do
 end do
 return
end
