!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine gen_sg_lat(indexsg,volume,wmat)
       implicit none
       integer indexsg
       real*8 wmat(3,3),volume
       real*8 rlat(6),pi
       logical llattice,l2d

       l2d=.false. ; pi=4.0d0*atan(1.0d0)
       do while(.true.)
       call genrlat(rlat)
       select case(indexsg)
       case(1:2)
       rlat(4:6)=rlat(4:6)*pi/2.0d0
!  1:2, Triclinic, a /= b /= c and alpha /= beta /= gamma

       case(3:4)
       rlat(4)=pi/2.0d0 ; rlat(6)=pi/2.0d0 ; rlat(5)=rlat(5)*pi/2.0d0

       case(6:7)
       rlat(4)=pi/2.0d0 ; rlat(6)=pi/2.0d0 ; rlat(5)=rlat(5)*pi/2.0d0

       case(5)
       rlat(1)=rlat(2) ; rlat(4)=rlat(5) ; rlat(4:5)=rlat(4:5)*pi ; rlat(6)=rlat(6)*pi/2.0d0

       case(8:9)
       rlat(1)=rlat(2) ; rlat(4)=rlat(5) ; rlat(4:5)=rlat(4:5)*pi ; rlat(6)=rlat(6)*pi/2.0d0

       case(10:11)
       rlat(4)=pi/2.0d0 ; rlat(6)=pi/2.0d0 ; rlat(5)=rlat(5)*pi/2.0d0

       case(12)
       rlat(1)=rlat(2) ; rlat(4)=rlat(5) ; rlat(4:5)=rlat(4:5)*pi ; rlat(6)=rlat(6)*pi/2.0d0

       case(13:14)
       rlat(4)=pi/2.0d0 ; rlat(6)=pi/2.0d0 ; rlat(5)=rlat(5)*pi/2.0d0

       case(15)
       rlat(1)=rlat(2) ; rlat(4)=rlat(5) ; rlat(4:5)=rlat(4:5)*pi ; rlat(6)=rlat(6)*pi/2.0d0
!  3:15, Monoclinic, a /= b /= c and alpha = gamma = 90 deg,  beta /= 90 deg

       case(16:19)
       rlat(4)=0.5d0 ; rlat(5)=0.5d0 ; rlat(6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(20:21)
       rlat(1)=rlat(2) ; rlat(4)=0.5d0*pi ; rlat(5)=0.5d0*pi ; rlat(6)=rlat(6)*pi/2.0d0

       case(22)
       rlat(4)=acos((rlat(2)**2+rlat(3)**2-rlat(1)**2)/2.0d0/rlat(2)/rlat(3))
       rlat(5)=acos((rlat(1)**2+rlat(3)**2-rlat(2)**2)/2.0d0/rlat(1)/rlat(3))
       rlat(6)=acos((rlat(2)**2+rlat(1)**2-rlat(3)**2)/2.0d0/rlat(1)/rlat(2))

       case(23:24)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=rlat(4:6)*pi/2.0d0
       if(cos(rlat(4))+cos(rlat(5)) > 0.0d0)then
       rlat(6)=acos(-1.0d0+(cos(rlat(4))+cos(rlat(5))))
                                            else
       rlat(6)=acos(-1.0d0-cos(rlat(4))-cos(rlat(5)))
                                            endif

       case(25:34)
       rlat(4)=0.5d0 ; rlat(5)=0.5d0 ; rlat(6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(35:41)
       rlat(1)=rlat(2) ; rlat(4)=0.5d0*pi ; rlat(5)=0.5d0*pi ; rlat(6)=rlat(6)*pi/2.0d0

       case(42:43)
       rlat(4)=acos((rlat(2)**2+rlat(3)**2-rlat(1)**2)/2.0/rlat(2)/rlat(3))
       rlat(5)=acos((rlat(1)**2+rlat(3)**2-rlat(2)**2)/2.0/rlat(1)/rlat(3))
       rlat(6)=acos((rlat(2)**2+rlat(1)**2-rlat(3)**2)/2.0/rlat(1)/rlat(2))

       case(44:46)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=rlat(4:6)*pi/2.0d0
       if(cos(rlat(4))+cos(rlat(5)) > 0.0d0)then
       rlat(6)=acos(-1.0d0+(cos(rlat(4))+cos(rlat(5))))
                                            else
       rlat(6)=acos(-1.0d0-cos(rlat(4))-cos(rlat(5)))
                                            endif

       case(47:62)
       rlat(4)=0.5d0 ; rlat(5)=0.5d0 ; rlat(6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(63:68)
       rlat(1)=rlat(2) ; rlat(4)=0.5d0*pi ; rlat(5)=0.5d0*pi ; rlat(6)=rlat(6)*pi/2.0d0

       case(69:70)
       rlat(4)=acos((rlat(2)**2+rlat(3)**2-rlat(1)**2)/2.0/rlat(2)/rlat(3))
       rlat(5)=acos((rlat(1)**2+rlat(3)**2-rlat(2)**2)/2.0/rlat(1)/rlat(3))
       rlat(6)=acos((rlat(2)**2+rlat(1)**2-rlat(3)**2)/2.0/rlat(1)/rlat(2))

       case(71:74)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=rlat(4:6)*pi/2.0d0
       if(cos(rlat(4))+cos(rlat(5)) > 0.0d0)then
       rlat(6)=acos(-1.0d0+(cos(rlat(4))+cos(rlat(5))))
                                            else
       rlat(6)=acos(-1.0d0-cos(rlat(4))-cos(rlat(5)))
                                            endif
!   16:74, Orthorhombic, a /= b /= c and alpha = beta = gamma= 90 deg

       case(75:78)
       rlat(1)=rlat(2) ; rlat(4:6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(81)
       rlat(1)=rlat(2) ; rlat(4:6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(83:86)
       rlat(1)=rlat(2) ; rlat(4:6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(89:96)
       rlat(1)=rlat(2) ; rlat(4:6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(99:106)
       rlat(1)=rlat(2) ; rlat(4:6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(111:118)
       rlat(1)=rlat(2) ; rlat(4:6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(123:138)
       rlat(1)=rlat(2) ; rlat(4:6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(79:80)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4)=rlat(5) ; rlat(4:6)=rlat(4:6)*pi/2.0d0
       if(cos(rlat(4))+cos(rlat(5)) > 0.0d0)then
       rlat(6)=acos(-1.0d0+(cos(rlat(4))+cos(rlat(5))))
                                            else
       rlat(6)=acos(-1.0d0-cos(rlat(4))-cos(rlat(5)))
                                            endif

       case(82)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4)=rlat(5) ; rlat(4:6)=rlat(4:6)*pi/2.0d0
       if(cos(rlat(4))+cos(rlat(5)) > 0.0d0)then
       rlat(6)=acos(-1.0d0+(cos(rlat(4))+cos(rlat(5))))
                                            else
       rlat(6)=acos(-1.0d0-cos(rlat(4))-cos(rlat(5)))
                                            endif

       case(87:88)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4)=rlat(5) ; rlat(4:6)=rlat(4:6)*pi/2.0d0
       if(cos(rlat(4))+cos(rlat(5)) > 0.0d0)then
       rlat(6)=acos(-1.0d0+(cos(rlat(4))+cos(rlat(5))))
                                            else
       rlat(6)=acos(-1.0d0-cos(rlat(4))-cos(rlat(5)))
                                            endif

       case(97:98)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4)=rlat(5) ; rlat(4:6)=rlat(4:6)*pi/2.0d0
       if(cos(rlat(4))+cos(rlat(5)) > 0.0d0)then
       rlat(6)=acos(-1.0d0+(cos(rlat(4))+cos(rlat(5))))
                                            else
       rlat(6)=acos(-1.0d0-cos(rlat(4))-cos(rlat(5)))
                                            endif

       case(107:110)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4)=rlat(5) ; rlat(4:6)=rlat(4:6)*pi/2.0d0
       if(cos(rlat(4))+cos(rlat(5)) > 0.0d0)then
       rlat(6)=acos(-1.0d0+(cos(rlat(4))+cos(rlat(5))))
                                            else
       rlat(6)=acos(-1.0d0-cos(rlat(4))-cos(rlat(5)))
                                            endif

       case(119:122)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4)=rlat(5) ; rlat(4:6)=rlat(4:6)*pi/2.0d0
       if(cos(rlat(4))+cos(rlat(5)) > 0.0d0)then
       rlat(6)=acos(-1.0d0+(cos(rlat(4))+cos(rlat(5))))
                                            else
       rlat(6)=acos(-1.0d0-cos(rlat(4))-cos(rlat(5)))
                                            endif

       case(139:142)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4)=rlat(5) ; rlat(4:6)=rlat(4:6)*pi/2.0d0
       if(cos(rlat(4))+cos(rlat(5)) > 0.0d0)then
       rlat(6)=acos(-1.0d0+(cos(rlat(4))+cos(rlat(5))))
                                            else
       rlat(6)=acos(-1.0d0-cos(rlat(4))-cos(rlat(5)))
                                            endif
!  75:142, Tetragonal, a = b, a /= c and alpha = beta = gamma= 90 deg

       case(143:145)
       rlat(1)=rlat(2) ; rlat(5)=0.5d0 ; rlat(4)=rlat(5) ; rlat(6)=0.6666666666666d0 ; rlat(4:6)=rlat(4:6)*pi

       case(147)
       rlat(1)=rlat(2) ; rlat(5)=0.5d0 ; rlat(4)=rlat(5) ; rlat(6)=0.6666666666666d0 ; rlat(4:6)=rlat(4:6)*pi

       case(149:154)
       rlat(1)=rlat(2) ; rlat(5)=0.5d0 ; rlat(4)=rlat(5) ; rlat(6)=0.6666666666666d0 ; rlat(4:6)=rlat(4:6)*pi

       case(156:159)
       rlat(1)=rlat(2) ; rlat(5)=0.5d0 ; rlat(4)=rlat(5) ; rlat(6)=0.6666666666666d0 ; rlat(4:6)=rlat(4:6)*pi

       case(162:165)
       rlat(1)=rlat(2) ; rlat(5)=0.5d0 ; rlat(4)=rlat(5) ; rlat(6)=0.6666666666666d0 ; rlat(4:6)=rlat(4:6)*pi

       case(168:194)
       rlat(1)=rlat(2) ; rlat(5)=0.5d0 ; rlat(4)=rlat(5) ; rlat(6)=0.6666666666666d0 ; rlat(4:6)=rlat(4:6)*pi

       case(146)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4)=rlat(5) ; rlat(6)=rlat(5) ; rlat(4:6)=rlat(4:6)*pi/2.0d0

       case(148)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4)=rlat(5) ; rlat(6)=rlat(5) ; rlat(4:6)=rlat(4:6)*pi/2.0d0

       case(155)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4)=rlat(5) ; rlat(6)=rlat(5) ; rlat(4:6)=rlat(4:6)*pi/2.0d0

       case(160:161)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4)=rlat(5) ; rlat(6)=rlat(5) ; rlat(4:6)=rlat(4:6)*pi/2.0d0

       case(166:167)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4)=rlat(5) ; rlat(6)=rlat(5) ; rlat(4:6)=rlat(4:6)*pi/2.0d0
! 143:167, Trigonal, a = b= c  and alpha = beta = gamma /= 90 deg
! 168:194, Hexagonal, a = b /= c and alpha, beta, gamma = 120 deg 

       case(195)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(198)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(200:201)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(205)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(207:208)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(212:213)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(215)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(218)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(221:224)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=0.5d0 ; rlat(4:6)=rlat(4:6)*pi

       case(196)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=1.0d0/3.0d0 ; rlat(4:6)=rlat(4:6)*pi

       case(202:203)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=1.0d0/3.0d0 ; rlat(4:6)=rlat(4:6)*pi

       case(209:210)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=1.0d0/3.0d0 ; rlat(4:6)=rlat(4:6)*pi

       case(216)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=1.0d0/3.0d0 ; rlat(4:6)=rlat(4:6)*pi

       case(219)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=1.0d0/3.0d0 ; rlat(4:6)=rlat(4:6)*pi

       case(225:228)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=1.0d0/3.0d0 ; rlat(4:6)=rlat(4:6)*pi

       case(197)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=109.471220634491d0/180.d0 ; rlat(4:6)=rlat(4:6)*pi

       case(199)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=109.471220634491d0/180.d0 ; rlat(4:6)=rlat(4:6)*pi

       case(204)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=109.471220634491d0/180.d0 ; rlat(4:6)=rlat(4:6)*pi

       case(206)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=109.471220634491d0/180.d0 ; rlat(4:6)=rlat(4:6)*pi

       case(211)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=109.471220634491d0/180.d0 ; rlat(4:6)=rlat(4:6)*pi

       case(214)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=109.471220634491d0/180.d0 ; rlat(4:6)=rlat(4:6)*pi

       case(217)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=109.471220634491d0/180.d0 ; rlat(4:6)=rlat(4:6)*pi

       case(220)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=109.471220634491d0/180.d0 ; rlat(4:6)=rlat(4:6)*pi

       case(229)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=109.471220634491d0/180.d0 ; rlat(4:6)=rlat(4:6)*pi

       case(230)
       rlat(1)=rlat(2) ; rlat(3)=rlat(2) ; rlat(4:6)=109.471220634491d0/180.d0 ; rlat(4:6)=rlat(4:6)*pi
!  195:230, Cubic
       end select
       call latmatvol(rlat,wmat,volume)
       call llcheck(llattice,l2d,wmat)
       if(llattice) return
       enddo
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine latmatvol(rlat,wmat,volume)
       implicit none
       real*8 rlat(6),wmat(3,3),volume
       real*8 ylat(6),slat(6),zmat(3,3),tmq,tmp

       slat=rlat
       call latmat(slat,zmat,1)
       tmp=(zmat(1,2)*zmat(2,3)-zmat(1,3)*zmat(2,2))*zmat(3,1)  &
          +(zmat(1,3)*zmat(2,1)-zmat(1,1)*zmat(2,3))*zmat(3,2)  &
          +(zmat(1,1)*zmat(2,2)-zmat(1,2)*zmat(2,1))*zmat(3,3)
       tmq=volume/tmp ; tmq=tmq**(1.0d0/3.0d0)
       ylat(1)=rlat(1)*tmq ; ylat(2)=rlat(2)*tmq ; ylat(3)=rlat(3)*tmq
       ylat(4)=rlat(4) ; ylat(5)=rlat(5) ; ylat(6)=rlat(6)
       call latmat(ylat,wmat,1)
       return
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine latmat(rlat,wmat,ksign)
       implicit none
       integer ksign
       real*8 rlat(6),wmat(3,3)
       real*8 ra,rb,rc,cosinea,cosineb,cosinec
       real*8 epslat,tmr
       integer i,j

       epslat=1.0d-6
       if(ksign == 1)then
       wmat=0.0d0
       wmat(1,1)=rlat(1)
       wmat(2,1)=rlat(2)*cos(rlat(6))
       wmat(2,2)=rlat(2)*sin(rlat(6))
       wmat(3,1)=rlat(3)*cos(rlat(5))
       wmat(3,2)=rlat(3)*cos(rlat(4))*sin(rlat(6)) &
        -((rlat(3)*cos(rlat(5))-rlat(3)*cos(rlat(4))*cos(rlat(6)))/tan(rlat(6)))
       tmr=rlat(3)**2-wmat(3,1)**2-wmat(3,2)**2  ; if(tmr <= 1.0d-12) tmr=0.0d0
       wmat(3,3)=sqrt(tmr)
       do i=1,3
       do j=1,3
       if(abs(wmat(i,j)) < epslat) wmat(i,j)=0.0d0
       enddo
       enddo
                     else
       rlat=0.0d0
       ra=sqrt(wmat(1,1)**2+wmat(1,2)**2+wmat(1,3)**2)
       rb=sqrt(wmat(2,1)**2+wmat(2,2)**2+wmat(2,3)**2)
       rc=sqrt(wmat(3,1)**2+wmat(3,2)**2+wmat(3,3)**2)
       cosinea=(wmat(2,1)*wmat(3,1)+wmat(2,2)*wmat(3,2)+wmat(2,3)*wmat(3,3))/rb/rc
       cosineb=(wmat(1,1)*wmat(3,1)+wmat(1,2)*wmat(3,2)+wmat(1,3)*wmat(3,3))/rc/ra
       cosinec=(wmat(1,1)*wmat(2,1)+wmat(1,2)*wmat(2,2)+wmat(1,3)*wmat(2,3))/ra/rb
       rlat(1)=ra ; rlat(2)=rb ; rlat(3)=rc
       rlat(4)=acos(cosinea) ; rlat(5)=acos(cosineb) ; rlat(6)=acos(cosinec)
                     endif
       return
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine genrlat(rlat)
       implicit none
       real*8 rlat(6)
       real*8 tmp,tmq
       integer k
       real ranmar

       tmp=0.10d0 ; tmq=1.0d0/tmp
       do while(.true.)
       do k=1,6
       rlat(k)=ranmar()
       enddo
       if((abs(rlat(1)/rlat(2)) > tmp .or. abs(rlat(1)/rlat(2)) < tmq) .and. &
          (abs(rlat(1)/rlat(3)) > tmp .or. abs(rlat(1)/rlat(3)) < tmq) .and. &
          (abs(rlat(2)/rlat(3)) > tmp .or. abs(rlat(2)/rlat(3)) < tmq)) exit
       enddo
       return
       end
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine llcheck(llattice,l2d,zmat)
       implicit none
       logical llattice,l2d
       real*8 zmat(3,3)
       real*8 altm(3,3),ra,rb,rc,alpha,beta,gama,pi,tmp
       integer i,j
 
       altm=zmat
       pi=4.0d0*atan(1.0d0)
!---{
!
       do i=1,3
       do j=1,3
       if(isnan(altm(i,j)))then
       llattice=.false.
                           return
                           endif
       enddo
       enddo
!
!---}
       ra=sqrt(altm(1,1)**2+altm(1,2)**2+altm(1,3)**2)
       rb=sqrt(altm(2,1)**2+altm(2,2)**2+altm(2,3)**2)
       rc=sqrt(altm(3,1)**2+altm(3,2)**2+altm(3,3)**2)
       alpha=(altm(2,1)*altm(3,1)+altm(2,2)*altm(3,2)+altm(2,3)*altm(3,3))/rb/rc
       beta=(altm(1,1)*altm(3,1)+altm(1,2)*altm(3,2)+altm(1,3)*altm(3,3))/rc/ra
       gama=(altm(1,1)*altm(2,1)+altm(1,2)*altm(2,2)+altm(1,3)*altm(2,3))/ra/rb  
       tmp=180.0d0/pi
       alpha=tmp*acos(alpha) ; beta=tmp*acos(beta) ; gama=tmp*acos(gama)
       llattice=.true.
       if(.not. l2d)then
       if(ra < 1.20d0 .or. rb < 1.20d0 .or. rc < 1.20d0) llattice=.false.
       if(alpha < 20.0d0 .or. alpha  > 160.0d0) llattice=.false.
       if(beta  < 20.0d0 .or.  beta  > 160.0d0) llattice=.false.
       if(gama  < 20.0d0 .or.  gama  > 160.0d0) llattice=.false.
       if(ra/rb >  6.0d0 .or. ra/rb  <   0.3d0) llattice=.false.
       if(ra/rc >  6.0d0 .or. ra/rc  <   0.3d0) llattice=.false.
       if(rb/rc >  6.0d0 .or. rb/rc  <   0.3d0) llattice=.false.
                    else
       if(ra    < 1.20d0  .or.    rb  <   1.20d0) llattice=.false.
       if(alpha <  20.0d0 .or. alpha  >  160.0d0) llattice=.false.
       if(beta  <  20.0d0 .or.  beta  >  160.0d0) llattice=.false.
       if(gama  <  20.0d0 .or. gama   >  160.0d0) llattice=.false.
       if(ra/rb >   6.0d0 .or. ra/rb  <    0.3d0) llattice=.false.
                    endif
       return
       end subroutine llcheck
!234567890
