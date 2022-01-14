!      Written by In-Ho Lee, KRISS, September 11, 2013.
       program genk
       implicit none
       real*8 dg,a1(3),a2(3),a3(3),scale0
       character*80 kname

       open(2,file='POSCAR',form='formatted')
       read(2,*)
       read(2,*) scale0
       read(2,*) a1(1),a1(2),a1(3)
       read(2,*) a2(1),a2(2),a2(3)
       read(2,*) a3(1),a3(2),a3(3)
       a1=a1*scale0
       a2=a2*scale0
       a3=a3*scale0
       close(2)

       kname='KPOINTS_012'
       dg=0.12d0
       call write_kpoints(dg,a1,a2,a3,kname)

       kname='KPOINTS_006'
       dg=0.060d0
       call write_kpoints(dg,a1,a2,a3,kname)

       kname='KPOINTS_003'
       dg=0.030d0
       call write_kpoints(dg,a1,a2,a3,kname)

       kname='KPOINTS_000'
       dg=0.000d0
       call write_kpoints(dg,a1,a2,a3,kname)
       stop
       end program genk
!
!      recycled from calypso code, September 11, 2013
       subroutine write_kpoints(dg,a1,a2,a3,kname)
       implicit none
       character*80 kname
       real*8 dg,a1(3),a2(3),a3(3)
       real*8 ra,rb,rc,ga,gb,gc,cosine_alpha,cosine_beta,cosine_gamma
       real*8 angle_alpha,angle_beta,angle_gamma
       real*8 r(3,3),g(3,3),c(3),volume,da,db,dc,pi
       integer ka,kb,kc
       integer iswitch

       iswitch=0
       if(dg <  0.d0) dg=0.06d0
       if(dg == 0.d0)then
       dg=0.02d0
       iswitch=1
                     endif
!
       pi=4.d0*atan(1.d0)
       r(1,:)=a1(:)
       r(2,:)=a2(:)
       r(3,:)=a3(:)
!-     Real Lattice Parameters ------------------------
       ra=sqrt(r(1,1)**2+r(1,2)**2+r(1,3)**2)
       rb=sqrt(r(2,1)**2+r(2,2)**2+r(2,3)**2)
       rc=sqrt(r(3,1)**2+r(3,2)**2+r(3,3)**2)
!-     Cell Angles ------------------------------------
       cosine_alpha=(r(2,1)*r(3,1)+r(2,2)*r(3,2)+r(2,3)*r(3,3))/rb/rc
       cosine_beta=(r(1,1)*r(3,1)+r(1,2)*r(3,2)+r(1,3)*r(3,3))/ra/rc
       cosine_gamma=(r(1,1)*r(2,1)+r(1,2)*r(2,2)+R(1,3)*r(2,3))/ra/rb
       angle_alpha=(acos(cosine_alpha)/pi)*180.d0
       angle_beta=(acos(cosine_beta)/pi)*180.d0
       angle_gamma=(acos(cosine_gamma)/pi)*180.d0
!-     Volume ------------------------------------------
       call cross(r(2,1),r(2,2),r(2,3),r(3,1),r(3,2),r(3,3),c(1),c(2),c(3))
       call dot(r(1,1),r(1,2),r(1,3),c(1),c(2),c(3),volume)
!-     Reciprocal Lattices -----------------------------
       g(1,:)=2.0d0*pi*c(:)/volume
       call cross(r(3,1),r(3,2),r(3,3),r(1,1),r(1,2),r(1,3),c(1),c(2),c(3))
       g(2,:)=2.0d0*pi*c(:)/volume
       call cross(r(1,1),r(1,2),r(1,3),r(2,1),r(2,2),r(2,3),c(1),c(2),c(3))
       g(3,:)=2.0d0*pi*c(:)/volume
!-     Reciprocal Lattice Parameters ------------------
       ga=sqrt(g(1,1)**2+g(1,2)**2+g(1,3)**2)
       gb=sqrt(g(2,1)**2+g(2,2)**2+g(2,3)**2)
       gc=sqrt(g(3,1)**2+g(3,2)**2+g(3,3)**2)
!-     Actual spacing in Monkhorst-Pack grid ----------
!-     Mesh parameters in Monkhorst-Pack grid ---------
       call kmfind(dg,ga,da,ka)
       call kmfind(dg,gb,db,kb)
       call kmfind(dg,gc,dc,kc)
       if(ka > 10) ka=10
       if(kb > 10) kb=10
       if(kc > 10) kc=10
       if(iswitch == 1)then
       ka=ka+1+ka/2
       kb=kb+1+kb/2
       kc=kc+1+kc/2
                       endif
       open(71,file=trim(kname),form='formatted')
       write(71,"(a1)") 'A'
       write(71,"(i1)")  0
       write(71,"(a1)") 'G'
       write(71,"(1x,i3,2x,i3,2x,i3)") ka,kb,kc
       write(71,"(1x,i2,2x,i2,2x,i2)") 0,0,0
!      write(71,"(1x,3f4.2)") 0.,0.,0.
       close(71)
       end subroutine write_kpoints

!--------------------------------------------------
!  subroutine to perform dot-product of 2 vectors (In1 and In2) and
!  write
!       the results into the scalar Out
!       written by Eric J. Simon, Cray Research, Inc. Fri Aug  2 1996
!
       subroutine dot(Inx1, Iny1, Inz1,Inx2, Iny2, Inz2, outt)
       implicit none
       real*8 Inx1, Iny1, Inz1, Inx2, Iny2, Inz2, outt
       outt = Inx1*Inx2 + Iny1*Iny2 + Inz1*Inz2
       end subroutine dot
!---------------------------------------------------
!  subroutine to perform cross-product of 2 vectors (In1 and In2) and
!  write
!       the results into a third vector (Out)
!       written by Eric J. Simon, Cray Research, Inc. Fri Aug  2 1996
!
!         i       j       k
!       In1x    In1y    In1z
!       In2x    In2y    In2z
!
!       In1XIn2 = In1Y*In2z-In1z*In2y, In1z*In2x-In1x*In2z,
!       In1x*In2y-In1y*In2x
!
       subroutine cross(In1x, In1y, In1z, In2x, In2y, In2z, Outx, Outy, Outz)
       implicit none
       real*8 In1x, In1y, In1z, In2x, In2y, In2z, Outx
       real*8 Outy, Outz

       Outx = In1y*In2z - In2y*In1z
       Outy = In2x*In1z - In1x*In2z
       Outz = In1x*In2y - In2x*In1y
       end subroutine cross
!----------------------------------------------------
! Subroutine to find out the MP Kmesh accordding to one quality of
! k-point separation.
!
!      recycled from calypso code, September 11, 2013
      subroutine kmfind(di,gi,dd,kd)
      implicit none
      integer kd
      real*8 gi,di,dd
      real*8 pi
      integer i

      pi=4.d0*atan(1.d0)
      kd=int(gi/di/2.d0/pi)
      if(kd == 0) kd=1
      dd=gi/kd/2.d0/pi
      if(dd .ge. di)then
      do i=1,10
      kd=kd+i
      dd=gi/kd/2.d0/pi
      if(dd .le. di) exit
      end do
      endif
      end subroutine kmfind
