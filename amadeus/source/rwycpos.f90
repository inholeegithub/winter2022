!234567890
!      List of Wyckoff positions. 
!      http://cci.lbl.gov/cctbx/cctbx_web.cgi
!      http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-wp-list
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine rwycpos(indexsg,pos,nwyc)
       implicit none
       integer indexsg,nwyc
       real*8 pos(3,nwyc)
       real*8, allocatable :: posdum(:,:),wrk11(:)
       integer, allocatable :: iwrk11(:)
       integer i,j
       real ranmar

       select case(indexsg)
       case(1)
       do i=1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(2)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.0d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.0d0
       pos(1,5)=0.5d0 ; pos(2,5)=0.5d0 ; pos(3,5)=0.0d0
       pos(1,6)=0.5d0 ; pos(2,6)=0.0d0 ; pos(3,6)=0.5d0
       pos(1,7)=0.0d0 ; pos(2,7)=0.5d0 ; pos(3,7)=0.5d0
       pos(1,8)=0.5d0 ; pos(2,8)=0.5d0 ; pos(3,8)=0.5d0
       do i=9,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(3)
       do i=1,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(4)
       do i=1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(5)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0    ; pos(2,i)=ranmar()  ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.0d0    ; pos(2,i)=ranmar()  ; pos(3,i)=0.5d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
 
       case(6)
       do i=1,int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(7)
       do i=1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(8)
       do i=1,int(nwyc/2)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(9)
       do i=1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(10)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.0d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.0d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.0d0
       pos(1,5)=0.5d0 ; pos(2,5)=0.5d0 ; pos(3,5)=0.0d0
       pos(1,6)=0.0d0 ; pos(2,6)=0.5d0 ; pos(3,6)=0.5d0
       pos(1,7)=0.5d0 ; pos(2,7)=0.0d0 ; pos(3,7)=0.5d0
       pos(1,8)=0.5d0 ; pos(2,8)=0.5d0 ; pos(3,8)=0.5d0
       do i=9,int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(11)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.0d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.0d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.5d0
       do i=5,int(nwyc/2)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(12)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.50d0 ; pos(3,2)=0.00d0
       pos(1,3)=0.00d0 ; pos(2,3)=0.00d0 ; pos(3,3)=0.50d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.50d0 ; pos(3,4)=0.50d0
       pos(1,5)=0.25d0 ; pos(2,5)=0.25d0 ; pos(3,5)=0.00d0
       pos(1,6)=0.25d0 ; pos(2,6)=0.25d0 ; pos(3,6)=0.50d0
       do i=7,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(13)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.0d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.0d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.0d0
       do i=5,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(14)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.0d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.0d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.5d0
       do i=5,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(15)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.50d0 ; pos(3,2)=0.00d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.25d0 ; pos(2,4)=0.25d0 ; pos(3,4)=0.50d0
       do i=5,int(nwyc/2)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(16)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.0d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.0d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.5d0
       pos(1,5)=0.5d0 ; pos(2,5)=0.5d0 ; pos(3,5)=0.0d0
       pos(1,6)=0.5d0 ; pos(2,6)=0.0d0 ; pos(3,6)=0.5d0
       pos(1,7)=0.0d0 ; pos(2,7)=0.5d0 ; pos(3,7)=0.5d0
       pos(1,8)=0.5d0 ; pos(2,8)=0.5d0 ; pos(3,8)=0.5d0
       do i=9,int(nwyc/13)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/13)+1,2*int(nwyc/13)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=2*int(nwyc/13)+1,3*int(nwyc/13)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/13)+1,4*int(nwyc/13)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/13)+1,5*int(nwyc/13)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=5*int(nwyc/13)+1,6*int(nwyc/13)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=6*int(nwyc/13)+1,7*int(nwyc/13)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=7*int(nwyc/13)+1,8*int(nwyc/13)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=8*int(nwyc/13)+1,9*int(nwyc/13)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=9*int(nwyc/13)+1,10*int(nwyc/13)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=10*int(nwyc/13)+1,11*int(nwyc/13)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=11*int(nwyc/13)+1,12*int(nwyc/13)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=12*int(nwyc/13)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(17)
       do i=1,int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(18)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(19)
       do i=1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(20)
       do i=1,int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(21)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.0d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.0d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.5d0
       do i=5,int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/8)+1,2*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=2*int(nwyc/8)+1,3*int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/8)+1,4*int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/8)+1,5*int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/8)+1,6*int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=6*int(nwyc/8)+1,7*int(nwyc/8)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=7*int(nwyc/8)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(22)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.00d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.25d0 ; pos(2,4)=0.25d0 ; pos(3,4)=0.75d0
       do i=5,int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.25d0
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(23)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.0d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.0d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.0d0
       do i=5,int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)+1
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(24)
       do i=1,int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(25)
       do i=1,int(nwyc/9)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/9)+1,2*int(nwyc/9)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/9)+1,3*int(nwyc/9)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/9)+1,4*int(nwyc/9)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/9)+1,5*int(nwyc/9)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/9)+1,6*int(nwyc/9)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=6*int(nwyc/9)+1,7*int(nwyc/9)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=7*int(nwyc/9)+1,8*int(nwyc/9)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=8*int(nwyc/9)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(26)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(27)
       do i=1,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(28)
       do i=1,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(29)
       do i=1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(30)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(31)
       do i=1,int(nwyc/2)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(32)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(33)
       do i=1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(34)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(35)
       do i=1,int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(36)
       do i=1,int(nwyc/2)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(37)
       do i=1,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(38)
       do i=1,int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(39)
       do i=1,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(40)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(41)
       do i=1,int(nwyc/2)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(42)
       do i=1,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(43)
       do i=1,int(nwyc/2)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(44)
       do i=1,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(45)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(46)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(47)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.0d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.0d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.5d0
       pos(1,5)=0.0d0 ; pos(2,5)=0.5d0 ; pos(3,5)=0.0d0
       pos(1,6)=0.5d0 ; pos(2,6)=0.5d0 ; pos(3,6)=0.0d0
       pos(1,7)=0.0d0 ; pos(2,7)=0.5d0 ; pos(3,7)=0.5d0
       pos(1,8)=0.5d0 ; pos(2,8)=0.5d0 ; pos(3,8)=0.5d0
       do i=9,int(nwyc/19)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/19)+1,2*int(nwyc/19)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=2*int(nwyc/19)+1,3*int(nwyc/19)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/19)+1,4*int(nwyc/19)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/19)+1,5*int(nwyc/19)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=5*int(nwyc/19)+1,6*int(nwyc/19)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=6*int(nwyc/19)+1,7*int(nwyc/19)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=7*int(nwyc/19)+1,8*int(nwyc/19)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=8*int(nwyc/19)+1,9*int(nwyc/19)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=9*int(nwyc/19)+1,10*int(nwyc/19)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=10*int(nwyc/19)+1,11*int(nwyc/19)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=11*int(nwyc/19)+1,12*int(nwyc/19)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=12*int(nwyc/19)+1,13*int(nwyc/19)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=13*int(nwyc/19)+1,14*int(nwyc/19)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=14*int(nwyc/19)+1,15*int(nwyc/19)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=15*int(nwyc/19)+1,16*int(nwyc/19)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=16*int(nwyc/19)+1,17*int(nwyc/19)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=17*int(nwyc/19)+1,18*int(nwyc/19)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=18*int(nwyc/19)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(48)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.50d0 ; pos(2,2)=0.50d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.75d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.25d0 ; pos(2,4)=0.25d0 ; pos(3,4)=0.75d0
       pos(1,5)=0.75d0 ; pos(2,5)=0.25d0 ; pos(3,5)=0.25d0
       pos(1,6)=0.25d0 ; pos(2,6)=0.25d0 ; pos(3,6)=0.25d0
       do i=7,int(nwyc/7)
       pos(1,i)=0.25d0 ; pos(2,i)=0.75d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.75d0
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.25d0
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(1,i)=0.75d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(49)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.00d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.00d0
       pos(1,5)=0.0d0 ; pos(2,5)=0.0d0 ; pos(3,5)=0.25d0
       pos(1,6)=0.5d0 ; pos(2,6)=0.0d0 ; pos(3,6)=0.25d0
       pos(1,7)=0.0d0 ; pos(2,7)=0.5d0 ; pos(3,7)=0.25d0
       pos(1,8)=0.5d0 ; pos(2,8)=0.5d0 ; pos(3,8)=0.25d0
       do i=9,int(nwyc/10)
       pos(1,i)=0.00d0 ; pos(2,i)=0.00d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/10)+1,2*int(nwyc/10)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/10)+1,3*int(nwyc/10)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/10)+1,4*int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/10)+1,5*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=5*int(nwyc/10)+1,6*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.25d0
       enddo
       do i=6*int(nwyc/10)+1,7*int(nwyc/10)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=7*int(nwyc/10)+1,8*int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=8*int(nwyc/10)+1,9*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=9*int(nwyc/10)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(50)
       pos(1,1)=0.25d0 ; pos(2,1)=0.25d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.75d0 ; pos(2,2)=0.25d0 ; pos(3,2)=0.0d0
       pos(1,3)=0.75d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.25d0 ; pos(2,4)=0.25d0 ; pos(3,4)=0.5d0
       pos(1,5)=0.00d0 ; pos(2,5)=0.00d0 ; pos(3,5)=0.0d0
       pos(1,6)=0.00d0 ; pos(2,6)=0.00d0 ; pos(3,6)=0.5d0
       do i=7,int(nwyc/7)
       pos(1,i)=0.25d0 ; pos(2,i)=0.75d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.5d0
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(51)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.0d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.0d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.5d0
       do i=5,int(nwyc/8)
       pos(1,i)=0.25d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/8)+1,2*int(nwyc/8)
       pos(1,i)=0.25d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/8)+1,3*int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/8)+1,4*int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/8)+1,5*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/8)+1,6*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=6*int(nwyc/8)+1,7*int(nwyc/8)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=7*int(nwyc/8)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(52)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       do i=3,int(nwyc/3)
       pos(1,i)=0.25d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(53)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.0d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.0d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.0d0
       do i=5,int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(54)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.0d0
       do i=3,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.25d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=0.25d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(55)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.0d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.5d0
       do i=5,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(56)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       do i=3,int(nwyc/3)
       pos(1,i)=0.25d0 ; pos(2,i)=0.75d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(57)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.0d0
       do i=3,int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(58)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.0d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.5d0
       do i=5,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(59)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       do i=3,int(nwyc/5)
       pos(1,i)=0.25d0 ; pos(2,i)=0.75d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(60)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.0d0
       do i=3,int(nwyc/2)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(61)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       do i=3,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(62)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       do i=3,int(nwyc/2)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(63)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.50d0 ; pos(3,2)=0.0d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.0d0
       do i=4,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(64)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.50d0 ; pos(2,2)=0.00d0 ; pos(3,2)=0.0d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.0d0
       do i=4,int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()  
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(65)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.50d0 ; pos(2,2)=0.00d0 ; pos(3,2)=0.0d0
       pos(1,3)=0.50d0 ; pos(2,3)=0.00d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.00d0 ; pos(3,4)=0.5d0
       pos(1,5)=0.25d0 ; pos(2,5)=0.25d0 ; pos(3,5)=0.0d0
       pos(1,6)=0.25d0 ; pos(2,6)=0.25d0 ; pos(3,6)=0.5d0
       do i=7,int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/12)+1,2*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=2*int(nwyc/12)+1,3*int(nwyc/12)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/12)+1,4*int(nwyc/12)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/12)+1,5*int(nwyc/12)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/12)+1,6*int(nwyc/12)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=6*int(nwyc/12)+1,7*int(nwyc/12)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=7*int(nwyc/12)+1,8*int(nwyc/12)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=8*int(nwyc/12)+1,9*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=9*int(nwyc/12)+1,10*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=10*int(nwyc/12)+1,11*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=11*int(nwyc/12)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(66)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.25d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.50d0 ; pos(3,2)=0.25d0
       pos(1,3)=0.00d0 ; pos(2,3)=0.00d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.50d0 ; pos(3,4)=0.00d0
       pos(1,5)=0.25d0 ; pos(2,5)=0.25d0 ; pos(3,5)=0.00d0
       pos(1,6)=0.25d0 ; pos(2,6)=0.75d0 ; pos(3,6)=0.00d0
       do i=7,int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(67)
       pos(1,1)=0.25d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.25d0 ; pos(2,2)=0.00d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.00d0 ; pos(2,3)=0.00d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.00d0 ; pos(3,4)=0.50d0
       pos(1,5)=0.25d0 ; pos(2,5)=0.25d0 ; pos(3,5)=0.00d0
       pos(1,6)=0.25d0 ; pos(2,6)=0.25d0 ; pos(3,6)=0.50d0
       do i=7,int(nwyc/9)
       pos(1,i)=0.0d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/9)+1,2*int(nwyc/9)
       pos(1,i)=0.25d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/9)+1,3*int(nwyc/9)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=3*int(nwyc/9)+1,4*int(nwyc/9)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/9)+1,5*int(nwyc/9)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=5*int(nwyc/9)+1,6*int(nwyc/9)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=6*int(nwyc/9)+1,7*int(nwyc/9)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=7*int(nwyc/9)+1,8*int(nwyc/9)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=8*int(nwyc/9)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(68)
       pos(1,1)=0.00d0 ; pos(2,1)=0.25d0 ; pos(3,1)=0.25d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.25d0 ; pos(3,2)=0.75d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.75d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.00d0 ; pos(3,4)=0.00d0
       do i=5,int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.25d0
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=0.25d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(69)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.00d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.00d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.25d0 ; pos(2,4)=0.00d0 ; pos(3,4)=0.25d0
       pos(1,5)=0.25d0 ; pos(2,5)=0.25d0 ; pos(3,5)=0.00d0
       pos(1,6)=0.25d0 ; pos(2,6)=0.25d0 ; pos(3,6)=0.25d0
       do i=7,int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/10)+1,2*int(nwyc/10)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/10)+1,3*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/10)+1,4*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.25d0
       enddo
       do i=4*int(nwyc/10)+1,5*int(nwyc/10)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=5*int(nwyc/10)+1,6*int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=6*int(nwyc/10)+1,7*int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=7*int(nwyc/10)+1,8*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=8*int(nwyc/10)+1,9*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=9*int(nwyc/10)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(70)
       pos(1,1)=0.125d0 ; pos(2,1)=0.125d0 ; pos(3,1)=0.125d0
       pos(1,2)=0.125d0 ; pos(2,2)=0.125d0 ; pos(3,2)=0.625d0
       pos(1,3)=0.000d0 ; pos(2,3)=0.000d0 ; pos(3,3)=0.000d0
       pos(1,4)=0.500d0 ; pos(2,4)=0.500d0 ; pos(3,4)=0.500d0
       do i=5,int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.125d0 ; pos(3,i)=0.125d0
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.125d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.125d0
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=0.125d0 ; pos(2,i)=0.125d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(71)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.50d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.50d0 ; pos(2,3)=0.50d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.50d0 ; pos(2,4)=0.00d0 ; pos(3,4)=0.50d0
       pos(1,5)=0.25d0 ; pos(2,5)=0.25d0 ; pos(3,5)=0.25d0
       do i=6,int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/10)+1,2*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/10)+1,3*int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/10)+1,4*int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/10)+1,5*int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/10)+1,6*int(nwyc/10)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=6*int(nwyc/10)+1,7*int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=7*int(nwyc/10)+1,8*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=8*int(nwyc/10)+1,9*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=9*int(nwyc/10)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(72)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.50d0 ; pos(2,2)=0.00d0 ; pos(3,2)=0.00d0
       pos(1,3)=0.50d0 ; pos(2,3)=0.00d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.00d0 ; pos(3,4)=0.25d0
       pos(1,5)=0.25d0 ; pos(2,5)=0.25d0 ; pos(3,5)=0.25d0
       do i=6,int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(73)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.25d0 ; pos(2,2)=0.25d0 ; pos(3,2)=0.25d0
       do i=3,int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(74)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.00d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.25d0 ; pos(2,4)=0.25d0 ; pos(3,4)=0.75d0
       do i=5,int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(75)
       do i=1,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(76)
       do i=1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(77)
       do i=1,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(78)
       do i=1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(79)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(80)
       do i=1,int(nwyc/2)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(81)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.0d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.5d0
       do i=5,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(82)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.75d0
       do i=5,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(83)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.0d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.5d0
       pos(1,5)=0.0d0 ; pos(2,5)=0.5d0 ; pos(3,5)=0.0d0
       pos(1,6)=0.0d0 ; pos(2,6)=0.5d0 ; pos(3,6)=0.5d0
       do i=7,int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(84)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.00d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.50d0
       pos(1,5)=0.0d0 ; pos(2,5)=0.0d0 ; pos(3,5)=0.25d0
       pos(1,6)=0.5d0 ; pos(2,6)=0.5d0 ; pos(3,6)=0.25d0
       do i=7,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(85)
       pos(1,1)=0.25d0 ; pos(2,1)=0.75d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.25d0 ; pos(2,2)=0.75d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.00d0 ; pos(2,3)=0.00d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.00d0 ; pos(3,4)=0.50d0
       do i=5,int(nwyc/3)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.25d0 ; pos(2,i)=0.75d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
 
       case(86)
       pos(1,1)=0.25d0 ; pos(2,1)=0.25d0 ; pos(3,1)=0.25d0
       pos(1,2)=0.25d0 ; pos(2,2)=0.25d0 ; pos(3,2)=0.75d0
       pos(1,3)=0.00d0 ; pos(2,3)=0.00d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.00d0 ; pos(3,4)=0.50d0
       do i=5,int(nwyc/3)
       pos(1,i)=0.75d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(87)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.00d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.00d0 ; pos(2,3)=0.50d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.50d0 ; pos(3,4)=0.25d0
       pos(1,5)=0.25d0 ; pos(2,5)=0.25d0 ; pos(3,5)=0.25d0
       do i=6,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(88)
       pos(1,1)=0.0d0 ; pos(2,1)=0.25d0 ; pos(3,1)=0.125d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.25d0 ; pos(3,2)=5.0d0/8.0d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.00d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.00d0 ; pos(3,4)=0.50d0
       do i=5,int(nwyc/2)
       pos(1,i)=0.0d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(89)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.0d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.5d0
       pos(1,5)=0.5d0 ; pos(2,5)=0.0d0 ; pos(3,5)=0.0d0
       pos(1,6)=0.5d0 ; pos(2,6)=0.0d0 ; pos(3,6)=0.5d0
       do i=7,int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/10)+1,2*int(nwyc/10)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/10)+1,3*int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/10)+1,4*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/10)+1,5*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/10)+1,6*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=6*int(nwyc/10)+1,7*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.5d0
       enddo
       do i=7*int(nwyc/10)+1,8*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=8*int(nwyc/10)+1,9*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=9*int(nwyc/10)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(90)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       do i=3,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.5d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(91)
       do i=1,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=3.0d0/8.0d0
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(92)
       do i=1,int(nwyc/2)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(93)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.00d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.50d0
       pos(1,5)=0.5d0 ; pos(2,5)=0.5d0 ; pos(3,5)=0.25d0
       pos(1,6)=0.0d0 ; pos(2,6)=0.0d0 ; pos(3,6)=0.25d0
       do i=7,int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/10)+1,2*int(nwyc/10)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/10)+1,3*int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/10)+1,4*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/10)+1,5*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/10)+1,6*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=6*int(nwyc/10)+1,7*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=7*int(nwyc/10)+1,8*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.75d0
       enddo
       do i=8*int(nwyc/10)+1,9*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.25d0
       enddo
       do i=9*int(nwyc/10)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(94)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       do i=3,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(95)
       do i=1,int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=5.0d0/8.0d0
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(96)
       do i=1,int(nwyc/2)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(97)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.25d0
       do i=5,int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i)+0.5d0 ; pos(3,i)=0.25d0
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(98)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       do i=3,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,2)=0.0d0
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,2)=0.0d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,2)=0.125d0
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(99)
       do i=1,int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ;  pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(100)
       do i=1,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i)+0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(101)
       do i=1,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(102)
       do i=1,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(103)
       do i=1,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(104)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(105)
       do i=1,int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(106)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(107)
       do i=1,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(108)
       do i=1,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i)+0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(109)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(110)
       do i=1,int(nwyc/2)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(111)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.0d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.0d0
       pos(1,5)=0.5d0 ; pos(2,5)=0.0d0 ; pos(3,5)=0.0d0
       pos(1,6)=0.5d0 ; pos(2,6)=0.0d0 ; pos(3,6)=0.5d0
       do i=7,int(nwyc/9)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/9)+1,2*int(nwyc/9)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/9)+1,3*int(nwyc/9)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/9)+1,4*int(nwyc/9)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/9)+1,5*int(nwyc/9)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/9)+1,6*int(nwyc/9)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=6*int(nwyc/9)+1,7*int(nwyc/9)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=7*int(nwyc/9)+1,8*int(nwyc/9)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=8*int(nwyc/9)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(112)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.25d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.25d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.25d0
       pos(1,5)=0.0d0 ; pos(2,5)=0.0d0 ; pos(3,5)=0.00d0
       pos(1,6)=0.5d0 ; pos(2,6)=0.5d0 ; pos(3,6)=0.00d0
       do i=7,int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/8)+1,2*int(nwyc/8)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/8)+1,3*int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/8)+1,4*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=4*int(nwyc/8)+1,5*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.25d0
       enddo
       do i=5*int(nwyc/8)+1,6*int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=6*int(nwyc/8)+1,7*int(nwyc/8)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=7*int(nwyc/8)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(113)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       do i=3,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i)+0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(114)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       do i=3,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(115)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.0d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.5d0
       do i=5,int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/8)+1,2*int(nwyc/8)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/8)+1,3*int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/8)+1,4*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/8)+1,5*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/8)+1,6*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=6*int(nwyc/8)+1,7*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=7*int(nwyc/8)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(116)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.25d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.25d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.0d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.00d0
       do i=5,int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.25d0
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.75d0
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
 
       case(117)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.50d0
       do i=5,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i)+0.5d0 ; pos(3,i)=0.5d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i)+0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(118)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.75d0
       do i=5,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0+pos(1,i) ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0-pos(1,i) ; pos(3,i)=0.25d0
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(119)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.75d0
       do i=5,int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0+pos(1,i) ; pos(3,i)=0.25d0
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(120)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.25d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.00d0
       do i=5,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i)+0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.25d0
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(121)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.25d0
       do i=5,int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(122)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       do i=3,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.125d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(123)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.0d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.5d0
       pos(1,5)=0.0d0 ; pos(2,5)=0.5d0 ; pos(3,5)=0.5d0
       pos(1,6)=0.0d0 ; pos(2,6)=0.5d0 ; pos(3,6)=0.0d0
       do i=7,int(nwyc/15)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/15)+1,2*int(nwyc/15)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/15)+1,3*int(nwyc/15)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/15)+1,4*int(nwyc/15)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/15)+1,5*int(nwyc/15)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/15)+1,6*int(nwyc/15)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=6*int(nwyc/15)+1,7*int(nwyc/15)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=7*int(nwyc/15)+1,8*int(nwyc/15)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=8*int(nwyc/15)+1,9*int(nwyc/15)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.5d0
       enddo
       do i=9*int(nwyc/15)+1,10*int(nwyc/15)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=10*int(nwyc/15)+1,11*int(nwyc/15)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=11*int(nwyc/15)+1,12*int(nwyc/15)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=12*int(nwyc/15)+1,13*int(nwyc/15)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=13*int(nwyc/15)+1,14*int(nwyc/15)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=14*int(nwyc/15)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(124)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.25d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.00d0
       pos(1,5)=0.0d0 ; pos(2,5)=0.5d0 ; pos(3,5)=0.00d0
       pos(1,6)=0.0d0 ; pos(2,6)=0.5d0 ; pos(3,6)=0.25d0
       do i=7,int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/8)+1,2*int(nwyc/8)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/8)+1,3*int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/8)+1,4*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=4*int(nwyc/8)+1,5*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.25d0
       enddo
       do i=5*int(nwyc/8)+1,6*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.25d0
       enddo
       do i=6*int(nwyc/8)+1,7*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=7*int(nwyc/8)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(125)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.00d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.75d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.75d0 ; pos(2,4)=0.25d0 ; pos(3,4)=0.0d0
       pos(1,5)=0.25d0 ; pos(2,5)=0.25d0 ; pos(3,5)=0.5d0
       pos(1,6)=0.25d0 ; pos(2,6)=0.25d0 ; pos(3,6)=0.0d0
       do i=7,int(nwyc/8)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/8)+1,2*int(nwyc/8)
       pos(1,i)=0.75d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/8)+1,3*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.5d0
       enddo
       do i=3*int(nwyc/8)+1,4*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/8)+1,5*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/8)+1,6*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=6*int(nwyc/8)+1,7*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=7*int(nwyc/8)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(126)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.25d0 ; pos(2,2)=0.75d0 ; pos(3,2)=0.00d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.75d0 ; pos(3,3)=0.75d0
       pos(1,4)=0.25d0 ; pos(2,4)=0.25d0 ; pos(3,4)=0.75d0
       pos(1,5)=0.25d0 ; pos(2,5)=0.25d0 ; pos(3,5)=0.25d0
       do i=6,int(nwyc/6)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=0.25d0 ; pos(2,i)=0.75d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.75d0 ; pos(3,i)=0.25d0
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.25d0
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(127)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.0d0
       do i=5,int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/8)+1,2*int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/8)+1,3*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i)+0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/8)+1,4*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i)+0.5d0 ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/8)+1,5*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/8)+1,6*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=6*int(nwyc/8)+1,7*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i)+0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=7*int(nwyc/8)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(128)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.00d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.00d0 ; pos(2,3)=0.50d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.50d0 ; pos(3,4)=0.25d0
       do i=5,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i)+0.5d0 ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(129)
       pos(1,1)=0.75d0 ; pos(2,1)=0.25d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.75d0 ; pos(2,2)=0.25d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.00d0 ; pos(2,3)=0.00d0 ; pos(3,3)=0.0d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.00d0 ; pos(2,4)=0.5d0
       do i=5,int(nwyc/7)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=0.75d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(130)
       pos(1,1)=0.75d0 ; pos(2,1)=0.25d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.75d0 ; pos(2,2)=0.25d0 ; pos(3,2)=0.25d0
       pos(1,3)=0.00d0 ; pos(2,3)=0.00d0 ; pos(3,3)=0.00d0
       do i=4,int(nwyc/4)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.75d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(131)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.00d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.50d0
       pos(1,5)=0.0d0 ; pos(2,5)=0.0d0 ; pos(3,5)=0.25d0
       pos(1,6)=0.5d0 ; pos(2,6)=0.5d0 ; pos(3,6)=0.25d0
       do i=7,int(nwyc/12)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/12)+1,2*int(nwyc/12)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/12)+1,3*int(nwyc/12)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/12)+1,4*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/12)+1,5*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/12)+1,6*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=6*int(nwyc/12)+1,7*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.5d0
       enddo
       do i=7*int(nwyc/12)+1,8*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.25d0
       enddo
       do i=8*int(nwyc/12)+1,9*int(nwyc/12)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=9*int(nwyc/12)+1,10*int(nwyc/12)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=10*int(nwyc/12)+1,11*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=11*int(nwyc/12)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(132)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.25d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.25d0
       pos(1,5)=0.0d0 ; pos(2,5)=0.5d0 ; pos(3,5)=0.25d0
       pos(1,6)=0.0d0 ; pos(2,6)=0.5d0 ; pos(3,6)=0.00d0
       do i=7,int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/10)+1,2*int(nwyc/10)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/10)+1,3*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/10)+1,4*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/10)+1,5*int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/10)+1,6*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=6*int(nwyc/10)+1,7*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.25d0
       enddo
       do i=7*int(nwyc/10)+1,8*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=8*int(nwyc/10)+1,9*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=9*int(nwyc/10)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(133)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0  ; pos(3,1)=0.00d0
       pos(1,2)=0.75d0 ; pos(2,2)=0.25d0  ; pos(3,2)=0.75d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.25d0  ; pos(3,3)=0.25d0
       pos(1,4)=0.75d0 ; pos(2,4)=0.25d0  ; pos(3,4)=0.00d0
       pos(1,5)=0.25d0 ; pos(2,5)=0.25d0  ; pos(3,5)=0.00d0
       do i=6,int(nwyc/6)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=0.75d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.25d0
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(134)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.00d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.25d0 ; pos(2,4)=0.25d0 ; pos(3,4)=0.25d0
       pos(1,5)=0.75d0 ; pos(2,5)=0.25d0 ; pos(3,5)=0.25d0
       pos(1,6)=0.25d0 ; pos(2,6)=0.75d0 ; pos(3,6)=0.25d0
       do i=7,int(nwyc/8)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/8)+1,2*int(nwyc/8)
       pos(1,i)=0.75d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/8)+1,3*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/8)+1,4*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.75d0
       enddo
       do i=4*int(nwyc/8)+1,5*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/8)+1,6*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=6*int(nwyc/8)+1,7*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=7*int(nwyc/8)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(135)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.25d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.25d0
       do i=5,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i)+0.5d0 ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(136)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.5d0 ; pos(3,4)=0.25d0
       do i=5,int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(137)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.75d0 ; pos(2,2)=0.25d0 ; pos(3,2)=0.25d0
       pos(1,3)=0.75d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.75d0
       do i=4,int(nwyc/5)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=0.75d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(138)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.75d0 ; pos(2,2)=0.25d0 ; pos(3,2)=0.75d0
       pos(1,3)=0.75d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.00d0 ; pos(3,4)=0.50d0
       do i=5,int(nwyc/6)
       pos(1,i)=0.25d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=0.75d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=0.5d0
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(139)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.00d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.00d0 ; pos(2,3)=0.50d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.50d0 ; pos(3,4)=0.25d0
       pos(1,5)=0.25d0 ; pos(2,5)=0.25d0 ; pos(3,5)=0.25d0
       do i=6,int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/10)+1,2*int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/10)+1,3*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/10)+1,4*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/10)+1,5*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=5*int(nwyc/10)+1,6*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i)+0.5 ; pos(3,i)=0.25d0
       enddo
       do i=6*int(nwyc/10)+1,7*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=7*int(nwyc/10)+1,8*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=8*int(nwyc/10)+1,9*int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=9*int(nwyc/10)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(140)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.25d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.00d0 ; pos(3,2)=0.00d0
       pos(1,3)=0.00d0 ; pos(2,3)=0.50d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.50d0 ; pos(3,4)=0.00d0
       pos(1,5)=0.25d0 ; pos(2,5)=0.25d0 ; pos(3,5)=0.25d0
       do i=6,int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/8)+1,2*int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/8)+1,3*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i)+0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/8)+1,4*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=0.25d0
       enddo
       do i=4*int(nwyc/8)+1,5*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=5*int(nwyc/8)+1,6*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=6*int(nwyc/8)+1,7*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i)+0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=7*int(nwyc/8)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(141)
       pos(1,1)=0.0d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.000d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.00d0 ; pos(3,2)=0.500d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.25d0 ; pos(3,3)=3.0d0/8.0d0
       pos(1,4)=0.0d0 ; pos(2,4)=0.75d0 ; pos(3,4)=0.125d0
       do i=5,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i)+0.25d0 ; pos(3,i)=7.0d0/8.0d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(142)
       pos(1,1)=0.0d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.000d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.25d0 ; pos(3,2)=0.125d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.25d0 ; pos(3,3)=3.0d0/8.0d0
       do i=4,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.25d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i)+0.25d0 ; pos(3,i)=1.0d0/8.0d0
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(143)
       do i=1,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=2.0d0/3.0d0 ; pos(2,i)=1.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(144)
       do i=1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(145)
       do i=1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(146)
       do i=1,int(nwyc/2)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(147)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.0d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.0d0
       do i=5,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(148)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.0d0 ; pos(3,3)=0.0d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.5d0
       do i=5,int(nwyc/2)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(149)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=1.0d0/3.0d0 ; pos(2,3)=2.0d0/3.0d0 ; pos(3,3)=0.0d0
       pos(1,4)=1.0d0/3.0d0 ; pos(2,4)=2.0d0/3.0d0 ; pos(3,4)=0.5d0
       pos(1,5)=2.0d0/3.0d0 ; pos(2,5)=1.0d0/3.0d0 ; pos(3,5)=0.0d0
       pos(1,6)=2.0d0/3.0d0 ; pos(2,6)=1.0d0/3.0d0 ; pos(3,6)=0.5d0
       do i=7,int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=2.0d0/3.0d0 ; pos(2,i)=1.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(150)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       do i=3,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(151)
       do i=1,int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=1.0d0/3.0d0
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=5.0d0/6.0d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(152)
       do i=1,int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=1.0d0/3.0d0
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=5.0d0/6.0d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(153)
       do i=1,int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=2.0d0/3.0d0
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=1.0d0/6.0d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(154)
       do i=1,int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=1.0d0/6.0d0
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=2.0d0/3.0d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(155)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       do i=3,int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(156)
       do i=1,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=2.0d0/3.0d0 ; pos(2,i)=1.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(157)
       do i=1,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(158)
       do i=1,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=2.0d0/3.0d0 ; pos(2,i)=1.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(159)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(160)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(161)
       do i=1,int(nwyc/2)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(162)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=1.0d0/3.0d0 ; pos(2,3)=2.0d0/3.0d0 ; pos(3,3)=0.0d0
       pos(1,4)=1.0d0/3.0d0 ; pos(2,4)=2.0d0/3.0d0 ; pos(3,4)=0.5d0
       pos(1,5)=0.5d0 ; pos(2,5)=0.0d0 ; pos(3,5)=0.0d0
       pos(1,6)=0.5d0 ; pos(2,6)=0.0d0 ; pos(3,6)=0.5d0
       do i=7,int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(163)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.25d0
       pos(1,3)=1.0d0/3.0d0 ; pos(2,3)=2.0d0/3.0d0 ; pos(3,3)=0.25d0
       pos(1,4)=2.0d0/3.0d0 ; pos(2,4)=1.0d0/3.0d0 ; pos(3,4)=0.25d0
       pos(1,5)=0.5d0 ; pos(2,5)=0.0d0 ; pos(3,5)=0.0d0
       do i=6,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(164)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.0d0 ; pos(3,3)=0.0d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.5d0
       do i=5,int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(165)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.25d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.0d0 ; pos(3,3)=0.00d0
       do i=4,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(166)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.0d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.0d0
       do i=5,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(167)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.25d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.0d0 ; pos(3,3)=0.00d0
       do i=4,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(168)
       do i=1,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(169)
       do i=1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(170)
       do i=1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(171)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(172)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=0.5d0 ; pos(2,i)=0.5d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(173)
       do i=1,int(nwyc/3)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(174)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=1.0d0/3.0d0 ; pos(2,3)=2.0d0/3.0d0 ; pos(3,3)=0.0d0
       pos(1,4)=1.0d0/3.0d0 ; pos(2,4)=2.0d0/3.0d0 ; pos(3,4)=0.5d0
       pos(1,5)=2.0d0/3.0d0 ; pos(2,5)=1.0d0/3.0d0 ; pos(3,5)=0.0d0
       pos(1,6)=2.0d0/3.0d0 ; pos(2,6)=1.0d0/3.0d0 ; pos(3,6)=0.5d0
       do i=7,int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=2.0d0/3.0d0 ; pos(2,i)=1.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(175)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=1.0d0/3.0d0 ; pos(2,3)=2.0d0/3.0d0 ; pos(3,3)=0.0d0
       pos(1,4)=1.0d0/3.0d0 ; pos(2,4)=2.0d0/3.0d0 ; pos(3,4)=0.5d0
       pos(1,5)=0.5d0 ; pos(2,5)=0.0d0 ; pos(3,5)=0.0d0
       pos(1,6)=0.5d0 ; pos(2,6)=0.0d0 ; pos(3,6)=0.5d0
       do i=7,int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar() 
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(176)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.25d0
       pos(1,3)=1.0d0/3.0d0 ; pos(2,3)=2.0d0/3.0d0 ; pos(3,3)=0.25d0
       pos(1,4)=2.0d0/3.0d0 ; pos(2,4)=1.0d0/3.0d0 ; pos(3,4)=0.25d0
       pos(1,5)=0.5d0 ; pos(2,5)=0.0d0 ; pos(3,5)=0.0d0
       do i=6,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(177)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=1.0d0/3.0d0 ; pos(2,3)=2.0d0/3.0d0 ; pos(3,3)=0.0d0
       pos(1,4)=1.0d0/3.0d0 ; pos(2,4)=2.0d0/3.0d0 ; pos(3,4)=0.5d0
       pos(1,5)=0.5d0 ; pos(2,5)=0.0d0 ; pos(3,5)=0.0d0
       pos(1,6)=0.5d0 ; pos(2,6)=0.0d0 ; pos(3,6)=0.5d0
       do i=7,int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/8)+1,2*int(nwyc/8)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/8)+1,3*int(nwyc/8)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/8)+1,4*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/8)+1,5*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/8)+1,6*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=6*int(nwyc/8)+1,7*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=0.5d0
       enddo
       do i=7*int(nwyc/8)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(178)
       do i=1,int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=2.0d0*pos(1,i) ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(179)
       do i=1,int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=2.0d0*pos(1,i) ; pos(3,i)=0.75d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(180)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.0d0 ; pos(3,3)=0.0d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.5d0
       do i=5,int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=2.0d0*pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=2.0d0*pos(1,i) ; pos(3,i)=0.5d0
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(181)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.5d0 ; pos(2,3)=0.0d0 ; pos(3,3)=0.0d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.5d0
       do i=5,int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=2.0d0*(pos(1,i)) ; pos(3,i)=0.0d0
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=2.0d0*(pos(1,i)) ; pos(3,i)=0.5d0
       enddo
       do i=6*int(nwyc/7),nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
 
       case(182)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.25d0
       pos(1,3)=1.0d0/3.0d0 ; pos(2,3)=2.0d0/3.0d0 ; pos(3,3)=0.25d0
       pos(1,4)=1.0d0/3.0d0 ; pos(2,4)=2.0d0/3.0d0 ; pos(3,4)=0.75d0
       do i=5,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=2.0d0*(pos(1,i)) ; pos(3,i)=0.25d0
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(183)
       do i=1,int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(184)
       do i=1,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(185)
       do i=1,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(186)
       do i=1,int(nwyc/4)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(187)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=1.0d0/3.0d0 ; pos(2,3)=2.0d0/3.0d0 ; pos(3,3)=0.0d0
       pos(1,4)=1.0d0/3.0d0 ; pos(2,4)=2.0d0/3.0d0 ; pos(3,4)=0.5d0
       pos(1,5)=2.0d0/3.0d0 ; pos(2,5)=1.0d0/3.0d0 ; pos(3,5)=0.0d0
       pos(1,6)=2.0d0/3.0d0 ; pos(2,6)=1.0d0/3.0d0 ; pos(3,6)=0.5d0
       do i=7,int(nwyc/9)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/9)+1,2*int(nwyc/9)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/9)+1,3*int(nwyc/9)
       pos(1,i)=2.0d0/3.0d0 ; pos(2,i)=1.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/9)+1,4*int(nwyc/9)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/9)+1,5*int(nwyc/9)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/9)+1,6*int(nwyc/9)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=6*int(nwyc/9)+1,7*int(nwyc/9)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=7*int(nwyc/9)+1,8*int(nwyc/9)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=8*int(nwyc/9),nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(188)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.25d0
       pos(1,3)=1.0d0/3.0d0 ; pos(2,3)=2.0d0/3.0d0 ; pos(3,3)=0.00d0
       pos(1,4)=1.0d0/3.0d0 ; pos(2,4)=2.0d0/3.0d0 ; pos(3,4)=0.25d0
       pos(1,5)=2.0d0/3.0d0 ; pos(2,5)=1.0d0/3.0d0 ; pos(3,5)=0.00d0
       pos(1,6)=2.0d0/3.0d0 ; pos(2,6)=1.0d0/3.0d0 ; pos(3,6)=0.25d0
       do i=7,int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=2.0d0/3.0d0 ; pos(2,i)=1.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=-pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(189)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=1.0d0/3.0d0 ; pos(2,3)=2.0d0/3.0d0 ; pos(3,3)=0.0d0
       pos(1,4)=1.0d0/3.0d0 ; pos(2,4)=2.0d0/3.0d0 ; pos(3,4)=0.5d0
       do i=5,int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/8)+1,2*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/8)+1,3*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=3*int(nwyc/8)+1,4*int(nwyc/8)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/8)+1,5*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/8)+1,6*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=6*int(nwyc/8)+1,7*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=7*int(nwyc/8)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(190)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.25d0
       pos(1,3)=1.0d0/3.0d0 ; pos(2,3)=2.0d0/3.0d0 ; pos(3,3)=0.25d0
       pos(1,4)=2.0d0/3.0d0 ; pos(2,4)=1.0d0/3.0d0 ; pos(3,4)=0.25d0
       do i=5,int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(191)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.5d0
       pos(1,3)=1.0d0/3.0d0 ; pos(2,3)=2.0d0/3.0d0 ; pos(3,3)=0.0d0
       pos(1,4)=1.0d0/3.0d0 ; pos(2,4)=2.0d0/3.0d0 ; pos(3,4)=0.5d0
       pos(1,5)=0.5d0 ; pos(2,5)=0.0d0 ; pos(3,5)=0.0d0
       pos(1,6)=0.5d0 ; pos(2,6)=0.0d0 ; pos(3,6)=0.5d0
       do i=7,int(nwyc/12)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/12)+1,2*int(nwyc/12)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/12)+1,3*int(nwyc/12)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/12)+1,4*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/12)+1,5*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/12)+1,6*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=2.0d0*pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=6*int(nwyc/12)+1,7*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=2.d0*pos(1,i) ; pos(3,i)=0.5d0
       enddo
       do i=7*int(nwyc/12)+1,8*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=8*int(nwyc/12)+1,9*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=2.d0*pos(1,i) ;  pos(3,i)=ranmar()
       enddo
       do i=9*int(nwyc/12)+1,10*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=10*int(nwyc/12)+1,11*int(nwyc/12)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0
       enddo
       do i=11*int(nwyc/12)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(192)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.25d0
       pos(1,3)=1.0d0/3.0d0 ; pos(2,3)=2.0d0/3.0d0 ; pos(3,3)=0.25d0
       pos(1,4)=1.0d0/3.0d0 ; pos(2,4)=2.0d0/3.0d0 ; pos(3,4)=0.00d0
       pos(1,5)=0.5d0 ; pos(2,5)=0.0d0 ; pos(3,5)=0.25d0
       pos(1,6)=0.5d0 ; pos(2,6)=0.0d0 ; pos(3,6)=0.00d0
       do i=7,int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=0.5d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=2.0d0*pos(1,i) ; pos(3,i)=0.25d0
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.0d0
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(193)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.25d0
       pos(1,3)=1.0d0/3.0d0 ; pos(2,3)=2.0d0/3.0d0 ; pos(3,3)=0.25d0
       pos(1,4)=1.0d0/3.0d0 ; pos(2,4)=2.0d0/3.0d0 ; pos(3,4)=0.00d0
       pos(1,5)=0.5d0 ; pos(2,5)=0.0d0 ; pos(3,5)=0.0d0
       do i=6,int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=2.0d0*pos(1,i) ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(194)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.25d0
       pos(1,3)=1.0d0/3.0d0 ; pos(2,3)=2.0d0/3.0d0 ; pos(3,3)=0.25d0
       pos(1,4)=1.0d0/3.0d0 ; pos(2,4)=2.0d0/3.0d0 ; pos(3,4)=0.75d0
       pos(1,5)=0.5d0 ; pos(2,5)=0.0d0 ; pos(3,5)=0.0d0
       do i=6,int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=0.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=1.0d0/3.0d0 ; pos(2,i)=2.0d0/3.0d0 ; pos(3,i)=ranmar()
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=2.d0*pos(1,i) ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=0.25d0
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=2.0d0*pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(195)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.0d0
       do i=5,int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(196)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.50d0 ; pos(2,2)=0.50d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.75d0 ; pos(2,4)=0.75d0 ; pos(3,4)=0.75d0
       do i=5,int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(197)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.0d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.5d0
       do i=3,1*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(198)
       do i=1,int(nwyc/2)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(199)
       do i=1,int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(200)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.0d0
       do i=5,int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/8)+1,2*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/8)+1,3*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/8)+1,4*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/8)+1,5*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.5d0
       enddo
       do i=5*int(nwyc/8)+1,6*int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=6*int(nwyc/8)+1,7*int(nwyc/8)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=7*int(nwyc/8)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
 
       case(201)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.25d0 ; pos(2,2)=0.25d0 ; pos(3,2)=0.25d0
       pos(1,3)=0.50d0 ; pos(2,3)=0.50d0 ; pos(3,3)=0.50d0
       pos(1,4)=0.25d0 ; pos(2,4)=0.75d0 ; pos(3,4)=0.75d0
       do i=5,int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.75d0 ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(202)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.50d0 ; pos(2,2)=0.50d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.25d0 ; pos(3,4)=0.25d0
       do i=5,int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.00d0 ; pos(3,i)=0.00d0
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(203)
       pos(1,1)=0.000d0 ; pos(2,1)=0.000d0 ; pos(3,1)=0.000d0
       pos(1,2)=0.500d0 ; pos(2,2)=0.500d0 ; pos(3,2)=0.500d0
       pos(1,3)=0.125d0 ; pos(2,3)=0.125d0 ; pos(3,3)=0.125d0
       pos(1,4)=0.625d0 ; pos(2,4)=0.625d0 ; pos(3,4)=0.625d0
       do i=5,int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=0.125d0 ; pos(3,i)=0.125d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(204)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.50d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.25d0
       do i=4,int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(205)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.5d0
       do i=3,int(nwyc/2)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=int(nwyc/2)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(206)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.25d0 ; pos(2,2)=0.25d0 ; pos(3,2)=0.25d0
       do i=3,int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(207)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.0d0
       do i=5,int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.5d0
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=pos(2,i)
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=pos(2,i)
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(208)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.25d0 ; pos(2,2)=0.25d0 ; pos(3,2)=0.25d0
       pos(1,3)=0.75d0 ; pos(2,3)=0.75d0 ; pos(3,3)=0.75d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.50d0 ; pos(3,4)=0.50d0
       pos(1,5)=0.25d0 ; pos(2,5)=0.00d0 ; pos(3,5)=0.50d0
       pos(1,6)=0.25d0 ; pos(2,6)=0.50d0 ; pos(3,6)=0.00d0
       do i=7,int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(2,i)=ranmar() ; pos(1,i)=0.25d0 ; pos(3,i)=-pos(2,i)+0.5d0
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(2,i)=ranmar() ; pos(1,i)=0.25d0 ; pos(3,i)=pos(2,i)+0.5d0
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
 
       case(209)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.50d0 ; pos(2,2)=0.50d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.25d0 ; pos(3,4)=0.25d0
       do i=5,int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=pos(2,i)
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=pos(2,i)
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.25d0
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(210)
       pos(1,1)=0.000d0 ; pos(2,1)=0.000d0 ; pos(3,1)=0.000d0
       pos(1,2)=0.500d0 ; pos(2,2)=0.500d0 ; pos(3,2)=0.500d0
       pos(1,3)=0.125d0 ; pos(2,3)=0.125d0 ; pos(3,3)=0.125d0
       pos(1,4)=0.625d0 ; pos(2,4)=0.625d0 ; pos(3,4)=0.625d0
       do i=5,int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(2,i)=ranmar() ; pos(1,i)=0.125d0 ; pos(3,i)=-pos(2,i)+0.25d0
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(211)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.50d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.25d0 ; pos(2,4)=0.50d0 ; pos(3,4)=0.00d0
       do i=5,int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=pos(2,i)
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(2,i)=ranmar() ; pos(1,i)=0.25d0 ; pos(3,i)=-pos(2,i)+0.5d0
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(212)
       pos(1,1)=0.125d0 ; pos(2,1)=0.125d0 ; pos(3,1)=0.125d0
       pos(1,2)=0.625d0 ; pos(2,2)=0.625d0 ; pos(3,2)=0.625d0
       do i=3,int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(2,i)=ranmar() ; pos(1,i)=0.125d0 ; pos(3,i)=-pos(2,i)+0.25d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(213)
       pos(1,1)=0.375d0 ; pos(2,1)=0.375d0 ; pos(3,1)=0.375d0
       pos(1,2)=0.875d0 ; pos(2,2)=0.875d0 ; pos(3,2)=0.875d0
       do i=3,int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(2,i)=ranmar() ; pos(1,i)=0.125d0 ; pos(3,i)=pos(2,i)+0.25d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(214)
       pos(1,1)=0.125d0 ; pos(2,1)=0.125d0 ; pos(3,1)=0.125d0
       pos(1,2)=0.875d0 ; pos(2,2)=0.875d0 ; pos(3,2)=0.875d0
       pos(1,3)=0.125d0 ; pos(2,3)=0.000d0 ; pos(3,3)=0.250d0
       pos(1,4)=0.625d0 ; pos(2,4)=0.000d0 ; pos(3,4)=0.250d0
       do i=5,int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(2,i)=ranmar() ; pos(1,i)=0.125d0 ; pos(3,i)=-pos(2,i)+0.25d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(2,i)=ranmar() ; pos(1,i)=0.125d0 ; pos(3,i)=pos(2,i)+0.25d0
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(215)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.0d0
       do i=5,int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.5d0
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(216)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.50d0 ; pos(2,2)=0.50d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.75d0 ; pos(2,4)=0.75d0 ; pos(3,4)=0.75d0
       do i=5,int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/5),nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(217)
       pos(1,1)=0.00d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.0d0
       do i=4,int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(218)
       pos(1,1)=0.00d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.0d0
       pos(1,4)=0.25d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.5d0
       do i=5,int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(219)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.25d0 ; pos(2,2)=0.25d0 ; pos(3,2)=0.25d0
       pos(1,3)=0.00d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.25d0 ; pos(2,4)=0.00d0 ; pos(3,4)=0.00d0
       do i=5,int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
 
       case(220)
       pos(1,1)=0.375d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.25d0
       pos(1,2)=0.875d0 ; pos(2,2)=0.0d0 ; pos(3,2)=0.25d0
       do i=3,int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/3)+1,2*int(nwyc/3)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/3)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(221)
       pos(1,1)=0.0d0 ; pos(2,1)=0.0d0 ; pos(3,1)=0.0d0
       pos(1,2)=0.5d0 ; pos(2,2)=0.5d0 ; pos(3,2)=0.5d0
       pos(1,3)=0.0d0 ; pos(2,3)=0.5d0 ; pos(3,3)=0.5d0
       pos(1,4)=0.5d0 ; pos(2,4)=0.0d0 ; pos(3,4)=0.0d0
       do i=5,int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/10)+1,2*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.5d0
       enddo
       do i=2*int(nwyc/10)+1,3*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=3*int(nwyc/10)+1,4*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=4*int(nwyc/10)+1,5*int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=pos(2,i)
       enddo
       do i=5*int(nwyc/10)+1,6*int(nwyc/10)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=pos(2,i)
       enddo
       do i=6*int(nwyc/10)+1,7*int(nwyc/10)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=7*int(nwyc/10)+1,8*int(nwyc/10)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=8*int(nwyc/10)+1,9*int(nwyc/10)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=9*int(nwyc/10)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(222)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.75d0 ; pos(2,2)=0.25d0 ; pos(3,2)=0.25d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.75d0 ; pos(3,4)=0.25d0
       do i=5,int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.25d0
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.75d0 ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=pos(2,i)
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(223)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.50d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.00d0 ; pos(3,3)=0.50d0
       pos(1,4)=0.25d0 ; pos(2,4)=0.50d0 ; pos(3,4)=0.00d0
       pos(1,5)=0.25d0 ; pos(2,5)=0.25d0 ; pos(3,5)=0.25d0
       do i=6,int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.5d0
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.5d0 ; pos(3,i)=0.0d0
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(2,i)=ranmar() ; pos(1,i)=0.25d0 ; pos(3,i)=pos(2,i)+0.5d0
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(224)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.25d0 ; pos(2,2)=0.25d0 ; pos(3,2)=0.25d0
       pos(1,3)=0.50d0 ; pos(2,3)=0.50d0 ; pos(3,3)=0.50d0
       pos(1,4)=0.25d0 ; pos(2,4)=0.75d0 ; pos(3,4)=0.75d0
       pos(1,5)=0.50d0 ; pos(2,5)=0.25d0 ; pos(3,5)=0.75d0
       do i=6,int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/7)+1,2*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/7)+1,3*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.75d0
       enddo
       do i=3*int(nwyc/7)+1,4*int(nwyc/7)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=-pos(2,i)
       enddo
       do i=4*int(nwyc/7)+1,5*int(nwyc/7)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=0.5d0+pos(2,i)
       enddo
       do i=5*int(nwyc/7)+1,6*int(nwyc/7)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=6*int(nwyc/7)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(225)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.50d0 ; pos(2,2)=0.50d0 ; pos(3,2)=0.50d0
       pos(1,3)=0.25d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.25d0 ; pos(3,4)=0.25d0
       do i=5,int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/8)+1,2*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=2*int(nwyc/8)+1,3*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.25d0
       enddo
       do i=3*int(nwyc/8)+1,4*int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=pos(2,i)
       enddo
       do i=4*int(nwyc/8)+1,5*int(nwyc/8)
       pos(1,i)=0.5d0 ; pos(2,i)=ranmar() ; pos(3,i)=pos(2,i)
       enddo
       do i=5*int(nwyc/8)+1,6*int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=6*int(nwyc/8)+1,7*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=7*int(nwyc/8),nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(226)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.25d0 ; pos(2,2)=0.25d0 ; pos(3,2)=0.25d0  
       pos(1,3)=0.25d0 ; pos(2,3)=0.00d0 ; pos(3,3)=0.00d0
       pos(1,4)=0.00d0 ; pos(2,4)=0.25d0 ; pos(3,4)=0.25d0
       do i=5,int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/6)+1,2*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=0.25d0 ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/6)+1,3*int(nwyc/6)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=3*int(nwyc/6)+1,4*int(nwyc/6)
       pos(1,i)=0.25d0 ; pos(2,i)=ranmar() ; pos(3,i)=pos(2,i)
       enddo
       do i=4*int(nwyc/6)+1,5*int(nwyc/6)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=5*int(nwyc/6)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(227)
       pos(1,1)=0.000d0 ; pos(2,1)=0.000d0 ; pos(3,1)=0.000d0
       pos(1,2)=0.500d0 ; pos(2,2)=0.500d0 ; pos(3,2)=0.500d0
       pos(1,3)=0.125d0 ; pos(2,3)=0.125d0 ; pos(3,3)=0.125d0
       pos(1,4)=0.375d0 ; pos(2,4)=0.375d0 ; pos(3,4)=0.375d0
       do i=5,int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/5)+1,2*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=0.125d0 ; pos(3,i)=0.125d0
       enddo
       do i=2*int(nwyc/5)+1,3*int(nwyc/5)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=3*int(nwyc/5)+1,4*int(nwyc/5)
       pos(1,i)=0.0d0  ; pos(2,i)=ranmar() ; pos(3,i)=-pos(2,i)
       enddo
       do i=4*int(nwyc/5)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(228)
       pos(1,1)=0.000d0 ; pos(2,1)=0.000d0 ; pos(3,1)=0.000d0
       pos(1,2)=0.125d0 ; pos(2,2)=0.125d0 ; pos(3,2)=0.125d0
       pos(1,3)=0.875d0 ; pos(2,3)=0.125d0 ; pos(3,3)=0.125d0
       pos(1,4)=0.250d0 ; pos(2,4)=0.250d0 ; pos(3,4)=0.250d0
       do i=5,int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.125d0 ; pos(3,i)=0.125d0
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(2,i)=ranmar() ; pos(1,i)=0.25d0 ; pos(3,i)=-pos(2,i)
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(229)
       pos(1,1)=0.00d0 ; pos(2,1)=0.00d0 ; pos(3,1)=0.00d0
       pos(1,2)=0.00d0 ; pos(2,2)=0.50d0 ; pos(3,2)=0.50d0 
       pos(1,3)=0.25d0 ; pos(2,3)=0.25d0 ; pos(3,3)=0.25d0
       pos(1,4)=0.25d0 ; pos(2,4)=0.00d0 ; pos(3,4)=0.50d0
       do i=5,int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.0d0
       enddo
       do i=1*int(nwyc/8)+1,2*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=2*int(nwyc/8)+1,3*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.50d0
       enddo
       do i=3*int(nwyc/8)+1,4*int(nwyc/8)
       pos(2,i)=ranmar() ; pos(1,i)=0.0d0 ; pos(3,i)=pos(2,i)
       enddo
       do i=4*int(nwyc/8)+1,5*int(nwyc/8)
       pos(2,i)=ranmar() ; pos(1,i)=0.25d0 ; pos(3,i)=-pos(2,i)+0.5d0
       enddo
       do i=5*int(nwyc/8)+1,6*int(nwyc/8)
       pos(1,i)=0.0d0 ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
       do i=6*int(nwyc/8)+1,7*int(nwyc/8)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=ranmar()
       enddo
       do i=7*int(nwyc/8)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo

       case(230)
       pos(1,1)=0.000d0 ; pos(2,1)=0.000d0 ; pos(3,1)=0.000d0
       pos(1,2)=0.125d0 ; pos(2,2)=0.125d0 ; pos(3,2)=0.125d0
       pos(1,3)=0.125d0 ; pos(2,3)=0.000d0 ; pos(3,3)=0.250d0
       pos(1,4)=0.375d0 ; pos(2,4)=0.000d0 ; pos(3,4)=0.250d0
       do i=5,int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=pos(1,i) ; pos(3,i)=pos(1,i)
       enddo
       do i=1*int(nwyc/4)+1,2*int(nwyc/4)
       pos(1,i)=ranmar() ; pos(2,i)=0.0d0 ; pos(3,i)=0.25d0
       enddo
       do i=2*int(nwyc/4)+1,3*int(nwyc/4)
       pos(2,i)=ranmar() ; pos(1,i)=0.125d0 ; pos(3,i)=-pos(2,i)+0.25d0
       enddo
       do i=3*int(nwyc/4)+1,nwyc
       pos(1,i)=ranmar() ; pos(2,i)=ranmar() ; pos(3,i)=ranmar()
       enddo
        end select
       do i=1,nwyc
       do j=1,3
       pos(j,i)=pos(j,i)-floor(pos(j,i))
       enddo
       enddo
       allocate(posdum(3,nwyc))
       posdum=pos
       allocate(wrk11(nwyc),iwrk11(nwyc))
       do i=1,nwyc
       wrk11(i)=ranmar()
       enddo
       call sortnr(nwyc,wrk11,iwrk11)
       do i=1,nwyc
       j=iwrk11(i) ; pos(:,i)=posdum(:,j)
       enddo
       deallocate(wrk11,iwrk11)
       deallocate(posdum)
       end
!234567890
