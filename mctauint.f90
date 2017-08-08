subroutine mctauint(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,&
&xface,yface,zface,rhokap,ntot,jmean,&
&xcell,ycell,zcell,sigma,tflag,iseed)
 implicit none

 include 'simplegrid.txt'

 integer tflag,iseed,xcell,ycell,zcell,k,i
 real*8 xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,sigma
 real ran2

 integer celli,cellj,cellk
 real*8 tau,taurun,taucell,d,d1,dcell,xcur,ycur,zcur,dsx,dsy,dsz
 real*8 dx,dy,dz,smax,delta
! !print*, "calculate next interaction position"
 delta=1.e-3

!***** tflag=0 means photon is in envelope
 tflag=0

!**** generate random optical depth tau
 tau=-alog(ran2(iseed))
  
!set the current photon coordinates.
!note that the origin of the (xcur,ycur,zcur) system is at the 
!bottom corner of the grid.
 xcur=xp+xmax
 ycur=yp+ymax
 zcur=zp+zmax

!!print*, "photon position: xcur, ycur, zcur",xcur, ycur, zcur

 celli=xcell
 cellj=ycell
 cellk=zcell

!!print*, "cell location: celli, cellj, cellk",celli, cellj, cellk

!set the cumulative distance and optical depth (d and taurun) 
!along the photon path to zero. 
 taurun=0.

 10   continue

 d=0.

!***** calculate smax -- maximum distance photon can travel
!(that is distance to edges of the simulation box)

 if(nxp.gt.0.) then           !photon travels in +ve x dir.
    dsx=(2.*xmax-xcur)/nxp    !distance to +ve x-face
 elseif(nxp.lt.0.) then       !photon travels in -ve x dir
    dsx=-xcur/nxp             !distance to -ve x-face
 elseif(nxp.eq.0.) then       !photon doenst move in x dir
    dsx=1.e2*xmax             !set arbitrarily high dist.
 endif

 if(nyp.gt.0.) then
    dsy=(2.*ymax-ycur)/nyp
 elseif(nyp.lt.0.) then
    dsy=-ycur/nyp
 elseif(nyp.eq.0.) then
    dsy=1.e2*ymax
 endif

 if(nzp.gt.0.) then
    dsz=(2.*zmax-zcur)/nzp
 elseif(nzp.lt.0.) then
    dsz=-zcur/nzp
 elseif(nzp.eq.0.) then
    dsz=1.e2*zmax
 endif

!______________________________________________
!the minimum of dsx,dsy,dsz gives us the first |
!face of the simulation box the photon would   |
!escape from if it continued travelling in the |
!present direction                             |
!----------------------------------------------!

 smax=amin1(dsx,dsy,dsz)
!if distance to travel to leave the box is very
!small (smax<0.001) then terminate photon (tflag=1)

 if(smax.lt.1.e-3) then
    tflag=1
    return
 endif

! !print*, "maximum distance photon can travel is:", smax
! !print*, "photon is at position",xcur, ycur, zcur

  
!***** integrate through grid
!photon travels from cell to cell until it either 
!reaches the interaction position (tau = taurun)
!or it leaves the box (d > 0.999*smax)

 do while((taurun.lt.tau).and.(d.lt.(.999*smax)))  !------START DO WHILE---!
 if (celli.gt.nxg) then
    !print*, "x problem celli, cellj, cellk",celli, cellj, cellk
 end if
 if (cellj.gt.nyg) then
    !print*, "y problem celli, cellj, cellk",celli, cellj, cellk
 end if
 if (cellk.gt.nzg) then
    !print*, "z problem celli, cellj, cellk",celli, cellj, cellk
 end if
 !!print*, "random optical depth of interaction is:", tau



!***** find distance to next x, y, and z cell walls.  
!***** note that dx is not the x-distance, but the actual distance along 
!*****the direction of travel to the next x-face, and likewise for dy and dz.
  if(nxp.gt.0.) then              !photon travels in +ve x dir.
    dx=(xface(celli+1)-xcur)/nxp  !distance to +ve x-face of cell
    if(dx.lt.1.e-5) then          !if distance small jump to next cell
       xcur=xface(celli+1)
       celli=celli+1
       dx=(xface(celli+1)-xcur)/nxp !distance to new +ve x-face of cell
    endif
  elseif(nxp.lt.0.) then          !photon travels in -ve x dir.
    dx=(xface(celli)-xcur)/nxp    !distance to -ve x-face of cell
    if(dx.lt.1.e-5) then          !if distance small jump to previous cell
       xcur=xface(celli)
       celli=celli-1
       dx=(xface(celli)-xcur)/nxp !distance to new -ve x-face of cell
    endif
  elseif(nxp.eq.0.) then          !photon doesnt move in x dir.
    dx=1.e2*xmax                  !set arbitrarily large dist.
  endif

  if(nyp.gt.0.) then
    dy=(yface(cellj+1)-ycur)/nyp
    if(dy.lt.1.e-5) then
       ycur=yface(cellj+1)
       cellj=cellj+1
       dy=(yface(cellj+1)-ycur)/nyp
    endif
  elseif(nyp.lt.0.) then
    dy=(yface(cellj)-ycur)/nyp
    if(dy.lt.1.e-5) then
       ycur=yface(cellj)
       dy=(yface(cellj-1)-ycur)/nyp
       cellj=cellj-1
    endif
  elseif(nyp.eq.0.) then
    dy=1.e2*ymax
  endif

  if(nzp.gt.0.) then
    dz=(zface(cellk+1)-zcur)/nzp
    if(dz.lt.1.e-5) then
       zcur=zface(cellk+1)
       cellk=cellk+1
       dz=(zface(cellk+1)-zcur)/nzp
    endif
  elseif(nzp.lt.0.) then
    dz=(zface(cellk)-zcur)/nzp
    if(dz.lt.1.e-5) then
       zcur=zface(cellk)
       dz=(zface(cellk-1)-zcur)/nzp
       cellk=cellk-1
    endif
  elseif(nzp.eq.0.) then
    dz=1.e2*zmax
  endif

!***** distances are only zero if photon is on cell wall.  if it is 
!***** on cell wall then set to arbitrary large distance, since we will
!***** in fact hit another wall
         if( (dx.eq.0.) .or. ((abs(dx)).lt.(1.e-3)) ) dx=1.e2*xmax
         if( (dy.eq.0.) .or. ((abs(dy)).lt.(1.e-3)) ) dy=1.e2*ymax
         if( (dz.eq.0.) .or. ((abs(dz)).lt.(1.e-3)) ) dz=1.e2*zmax

!***** find distance to next cell wall -- minimum of dx, dy, and dz
!____________________________________________
!the minimum of dx,dy,dz gives us the first |
!face of the cell the photon would escape   |
!from if it continued travelling in the     |
!present direction                          |
!-------------------------------------------!

   dcell=amin1(dx,dy,dz)
   !!print*, "distance to travel within cell", dcell

   if(dcell.le.0.) then
     !!print *,'tauint2: dcell < 0'
     tflag=1
     return           
   endif
   if(dx.lt.0.) dcell=amin1(dy,dz)
   if(dy.lt.0.) dcell=amin1(dx,dz)
   if(dz.lt.0.) dcell=amin1(dx,dy)

!***** optical depth to next cell wall is 
!***** taucell= (distance to cell)*(opacity of current cell)
!         taucell=dcell*rhokap(celli,cellj,cellk)*sigma

   taucell=dcell*(rhokap(celli,cellj,cellk)*sigma)
   !!print*, "tau within cell", taucell
   !!print*, "rhokap, sigma", rhokap(celli,cellj,cellk), sigma

!***** if taurun+taucell>tau then photon will interact
!at a distance d+d1.  
 
    if((taurun+taucell).ge.tau) then
      
      d1=(tau-taurun)/(rhokap(celli,cellj,cellk)*sigma)
      d=d+d1
!update photon position 
      xcur=xcur+d1*nxp
      ycur=ycur+d1*nyp
      zcur=zcur+d1*nzp
!update pathlenght sum of the cell
      jmean(celli,cellj,cellk)=jmean(celli,cellj,cellk)&
            &+d1*sigma
!update cell
      celli=int(nxg*xcur/(2.*xmax))+1
      cellj=int(nyg*ycur/(2.*ymax))+1
      cellk=int(nzg*zcur/(2.*zmax))+1
!update taurun
!(this will make us exit do while loop)
      taurun=taurun+taucell

    else
!***** if taurun+taucell<tau then photon will not interact
!yet. Travel dcell and move to next cell.
      d=d+dcell
!update photon position 
      xcur=xcur+dcell*nxp
      ycur=ycur+dcell*nyp
      zcur=zcur+dcell*nzp
!update pathlenght sum of the cell
      jmean(celli,cellj,cellk)=jmean(celli,cellj,cellk)+dcell*sigma
!update cell     
      celli=int(nxg*xcur/(2.*xmax))+1
      cellj=int(nyg*ycur/(2.*ymax))+1
      cellk=int(nzg*zcur/(2.*zmax))+1
!update taurun
      taurun=taurun+taucell
    endif
    if (celli.gt.nxg) then
       !print*, "2x problem celli, cellj, cellk",celli, cellj, cellk
    end if
    if (cellj.gt.nyg) then
       !print*, "2y problem celli, cellj, cellk",celli, cellj, cellk
    end if
    if (cellk.gt.nzg) then
       !print*, "2z problem celli, cellj, cellk",celli, cellj, cellk
    end if

    !!print*, "optical depth travelled so far:", taurun

 end do                                       !------END DO WHILE---!

!***** calculate photon final position.  if it escapes envelope then
!***** set tflag=1.  if photon doesn't escape leave tflag=0 and update 
!***** photon position.
      if((d.ge.(.999*smax))) then
!          tflag=1
!           !print*, "terminate: photon out of the grid"
           if(zcur.lt.1.e-3.or.zcur.gt.0.995*2.*zmax) then
            tflag=1
           elseif(xcur.lt.1.e-3) then
            xcur=2.*xmax-delta
            celli=int(nxg*xcur/(2.*xmax))+1
            tflag=0
           elseif(xcur.gt.0.995*2.*xmax) then
            xcur=delta
            celli=int(nxg*xcur/(2.*xmax))+1
            tflag=0
           elseif(ycur.lt.1.e-3) then
            ycur=2.*ymax-delta
            cellj=int(nyg*ycur/(2.*ymax))+1
            tflag=0
           elseif(ycur.gt.0.995*2.*ymax) then
            ycur=delta
            cellj=int(nyg*ycur/(2.*ymax))+1
            tflag=0
           endif
!          !!print*, "tflag =1"
          tflag=1 ! uncomment this line for no repeating boundaries

         if(tflag.eq.0) then
             xp=xcur-xmax
             yp=ycur-ymax
             zp=zcur-zmax
!             !!print *,xcur,ycur,zcur,celli,cellj,cellk
!             !!print *,tau,taurun
!             !!print *, 'hello',xp,yp,zp
             goto 10
          endif
       if (celli.gt.nxg) then
         !print*, "3x problem celli, cellj, cellk",celli, cellj, cellk
       end if
       if (cellj.gt.nyg) then
         !print*, "3y problem celli, cellj, cellk",celli, cellj, cellk
      end if
      if (cellk.gt.nzg) then
         !print*, "3z problem celli, cellj, cellk",celli, cellj, cellk
      end if

      else
!update photon position
         xp=xp+d*nxp
         yp=yp+d*nyp
         zp=zp+d*nzp
         xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
         ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
         zcell=int(nzg*(zp+zmax)/(2.*zmax))+1

           if (xcell.gt.nxg) then
             !print*, "5x problem celli, cellj, cellk",xcell, ycell, zcell
           end if
           if (ycell.gt.nyg) then
             !print*, "5y problem celli, cellj, cellk",xcell, ycell, zcell
           end if
           if (zcell.gt.nzg) then
             !print*, "5z problem celli, cellj, cellk",xcell, ycell, zcell
           end if
      endif

 return
end
