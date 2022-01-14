!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine xyz2rtp(x,y,z,r,theta,phi)
       implicit none
       real*8 x,y,z,r,theta,phi
       real*8 pi,tol

       pi=4.0d0*atan(1.0d0)
       tol=0.001d0
       r=sqrt(x*x+y*y+z*z)
       if(abs(x) < tol  .and. abs(y) < tol)then 
       if(abs(z)  < tol)then
       theta=0.0d0
       elseif(z < 0.0d0)then
       theta=pi
                        else
       theta=0.d0
                        endif
                                           else
       theta=acos(z/r)
                                           endif
       if(abs(x) < tol .and.  abs(y) < tol)then
       phi=0.0d0
       elseif(abs(x) < tol .and. y > 0.0d0)then
       phi=pi/2.0d0
       elseif(abs(x) < tol .and. y < 0.0d0)then
       phi=3.0d0*pi/2.0d0
       elseif(abs(y) < tol .and. x > 0.0d0)then
       phi=0.0d0
       elseif(abs(y) < tol .and. x < 0.0d0)then
       phi=pi
       elseif(x > 0.0d0 .and.    y > 0.0d0)then
       phi=atan(y/x)
       elseif(x < 0.0d0 .and.    y > 0.0d0)then
       phi=atan(y/x) + pi
       elseif(x < 0.0d0 .and.    y < 0.0d0)then
       phi=atan(y/x) + pi
       elseif(x > 0.0d0 .and.    y < 0.0d0)then
       phi=atan(y/x)+2.0d0*pi
                                           endif
       end subroutine xyz2rtp
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine rtp2xyz(r,theta,phi,x,y,z)
       implicit none
       real*8 x,y,z,r,theta,phi

       x=r*sin(theta)*cos(phi) ; y=r*sin(theta)*sin(phi) ; z=r*cos(theta)
       end subroutine rtp2xyz
!234567890
!      Written by In-Ho Lee, KRISS, September 11, 2013.
       subroutine sphhar(l,m,theta,phi,ylm)
       implicit none
       integer l,m
       real*8 theta,phi,sina,cosa,sinb,cosb,pi,c,d,twopi
       complex*16 i,ylm

       i=(0.0d0,1.0d0) ; pi=4.0d0*atan(1.0d0) ; twopi=2.0d0*pi
       sina=sin(theta)
       cosa=cos(theta)
     if(l==0)    then
  if(m==0) ylm=cmplx(0.5d0*sqrt(1.0d0/pi),0.0d0)
     elseif(l==1)then
  if(m==-1) ylm=0.5d0*sqrt(3.0d0/twopi)*exp(-i*phi)*sin(theta)
  if(m==1)  ylm=-0.5d0*sqrt(3.0d0/twopi)*exp(i*phi)*sin(theta)
  if(m==0)  ylm=0.5d0*sqrt(3.0d0/pi)*cos(theta)
     elseif(l==2)then
  if(m==-2) ylm=0.25d0*sqrt(15.0d0/twopi)*exp(-2.0d0*i*phi)*sin(theta)**2
  if(m==2)  ylm=0.25d0*sqrt(15.0d0/twopi)*exp(2.0d0*i*phi)*sina**2
  if(m==-1) ylm=0.5d0*sqrt(15.0d0/twopi)*exp(-i*phi)*sin(theta)*cos(theta) 
  if(m==1)  ylm=-0.5d0*sqrt(15.0d0/twopi)*exp(i*phi)*sin(theta)*cos(theta) 
  if(m==0)  ylm=0.25d0*sqrt(5.0d0/pi)*(3.0d0*(cos(theta))**2-1.0d0)
     elseif(l==3)then
  if(m==-3) ylm=1.0d0/8.0d0*sqrt(35.0d0/pi)*exp(-3.0d0*i*phi)*(sin(theta))**3
  if(m==3)  ylm=-1.0d0/8.0d0*sqrt(35.0d0/pi)*exp(3.0d0*i*phi)*(sin(theta))**3
  if(m==-2) ylm=1.0d0/4.0d0*sqrt(105.0d0/twopi)*exp(-2.0d0*i*phi)*(sin(theta))**2*cos(theta)
  if(m==2)  ylm=1.0d0/4.0d0*sqrt(105.0d0/twopi)*exp(2.0d0*i*phi)*(sin(theta))**2*cos(theta)
  if(m==-1) ylm=1.0d0/8.0d0*sqrt(21.0d0/pi)*exp(-i*phi)*sin(theta)*(5.0d0*(cos(theta))**2-1.0d0)
  if(m==1)  ylm=-1.0d0/8.0d0*sqrt(21.0d0/pi)*exp(i*phi)*sin(theta)*(5.0d0*(cos(theta))**2-1.0d0)
  if(m==0)  ylm=1.0d0/4.0d0*sqrt(7.0d0/pi)*(5.0d0*(cos(theta))**3-3.0d0*cos(theta))
     elseif(l==4)then
  if(m==-4) ylm=3.0d0/16.0d0*sqrt(35.0d0/twopi)*exp(-4.0*i*phi)*(sin(theta))**4
  if(m==4)  ylm=3.0d0/16.0d0*sqrt(35.0d0/twopi)*exp(4.0d0*i*phi)*sina**4
  if(m==-3) ylm=3.0d0/8.0d0*sqrt(35.0d0/pi)*exp(-3.0d0*i*phi)*(sin(theta))**3*cos(theta)
  if(m==3)  ylm=-3.d0/8.0d0*sqrt(35.0d0/pi)*exp(3.0d0*i*phi)*sina**3*cosa
  if(m==-2) ylm=3.0d0/8.0d0*sqrt(5.0d0/twopi)*exp(-2.0d0*i*phi)*(sin(theta))**2*(7.0d0*(cos(theta))**2-1.0d0)
  if(m==2)  ylm=3.0d0/8.0d0*sqrt(5.0d0/twopi)*exp(2.0d0*i*phi)*sina**2*(7.0d0*cosa**2-1.0d0)
  if(m==-1) ylm=3.0d0/8.0d0*sqrt(5.0d0/pi)*exp(-i*phi)*sin(theta)*(7.0d0*(cos(theta))**3-3.0d0*cos(theta))
  if(m==1)  ylm=-3.0d0/8.0d0*sqrt(5.0d0/pi)*exp(i*phi)*sina*(7.0d0*cosa**3-3.0d0*cosa)
  if(m==0)  ylm=3.0d0/16.0d0*sqrt(1.0d0/pi)*(35.0d0*cosa**4-30.0d0*cosa**2+3.0d0)
     elseif(l==5)then
  if(m==-5) ylm=3.0d0/32.0d0*sqrt(77.0d0/pi)*exp(-5.0d0*i*phi)*sina**5
  if(m==5)  ylm=-3.0d0/32.0d0*sqrt(77.0d0/pi)*exp(5.0d0*i*phi)*sina**5
  if(m==-4) ylm=3.0d0/16.0*sqrt(385.0d0/twopi)*exp(-4.0d0*i*phi)*sina**4*cosa
  if(m==4)  ylm=3.0d0/16.0d0*sqrt(385.0d0/twopi)*exp(4.0d0*i*phi)*sina**4*cosa
  if(m==-3) ylm=1.0d0/32.0d0*sqrt(385.0d0/pi)*exp(-3.0d0*i*phi)*sina**3*(9.0d0*cosa**2-1.0d0)
  if(m==3)  ylm=-1.0d0/32.0d0*sqrt(385.0d0/pi)*exp(3.0d0*i*phi)*sina**3*(9.0d0*cosa**2-1.0d0)
  if(m==-2) ylm=1.0d0/8.0d0*sqrt(1155.0d0/twopi)*exp(-2.0d0*i*phi)*sina**2*(3.0d0*cosa**3-cosa)
  if(m==2)  ylm=1.0d0/8.0d0*sqrt(1155.0d0/twopi)*exp(2.0d0*i*phi)*sina**2*(3.0d0*cosa**3-cosa)
  if(m==-1) ylm=1.d0/16.0d0*sqrt(165.0d0/twopi)*exp(-i*phi)*sina*(21.0d0*cosa**4-14.0d0*cosa**2+1.0d0)
  if(m==1)  ylm=-1.0d0/16.0*sqrt(165.0d0/twopi)*exp(i*phi)*sina*(21.0d0*cosa**4-14.0d0*cosa**2+1.0d0)
  if(m==0)  ylm=1.0d0/16.0d0*sqrt(11.0d0/pi)*(63.0d0*cosa**5-70.0d0*cosa**3+15.0d0*cosa)
     elseif(l==6)then
  if(m==-6) ylm=1.0d0/64.0d0*sqrt(3003.0d0/pi)*exp(-6.0d0*i*phi)*sina**6
  if(m==6)  ylm=1.0d0/64.0d0*sqrt(3003.0d0/pi)*exp(6.0d0*i*phi)*sina**6
  if(m==-5) ylm=3.0d0/32.0d0*sqrt(1001.0d0/pi)*exp(-5.0d0*i*phi)*sina**5*cosa
  if(m==5)  ylm=-3.0d0/32.0d0*sqrt(1001.0d0/pi)*exp(5.0d0*i*phi)*sina**5*cosa
  if(m==-4) ylm=3.0d0/32.0d0*sqrt(91.0d0/twopi)*exp(-4.0d0*i*phi)*sina**4*(11.0d0*cosa**2-1.0d0)
  if(m==4)  ylm=3.0d0/32.0d0*sqrt(91.0d0/twopi)*exp(4.0d0*i*phi)*sina**4*(11.0d0*cosa**2-1.0d0)
  if(m==-3) ylm=1.0d0/32.0d0*sqrt(1365.0d0/pi)*exp(-3.0d0*i*phi)*sina**3*(11.0d0*cosa**3-3.0d0*cosa)
  if(m==3)  ylm=-1.0d0/32.0d0*sqrt(1365.0d0/pi)*exp(3.0d0*i*phi)*sina**3*(11.0d0*cosa**3-3.0d0*cosa)
  if(m==-2) ylm=1.0d0/64.0d0*sqrt(1365.0d0/pi)*exp(-2.0d0*i*phi)*sina**2*(33.0d0*cosa**4-18.0d0*cosa**2+1.0d0)
  if(m==2)  ylm=1.0d0/64.0d0*sqrt(1365.0d0/pi)*exp(2.0d0*i*phi)*sina**2*(33.0d0*cosa**4-18.0d0*cosa**2+1.0d0)
  if(m==-1) ylm=1.0d0/16.0d0*sqrt(273.0d0/twopi)*exp(-i*phi)*sina*(33.0d0*cosa**5-30.0d0*cosa**3+5.0d0*cosa) 
  if(m==1)  ylm=-1.0d0/16.0*sqrt(273.0d0/twopi)*exp(i*phi)*sina*(33.0d0*cosa**5-30.0d0*cosa**3+5.0d0*cosa)
  if(m==0)  ylm=1.0d0/32.0d0*sqrt(13.0d0/pi)*(231.0d0*cosa**6-315.0d0*cosa**4+105.0d0*cosa**2-5.0d0)
     elseif(l==7)then
  if(m==-7) ylm=3.0d0/64.0d0*sqrt(715.0d0/twopi)*exp(-7.0d0*i*phi)*sina**7
  if(m==7)  ylm=-3.0d0/64.0d0*sqrt(715.0d0/twopi)*exp(7.0d0*i*phi)*sina**7
  if(m==-6) ylm=3.0d0/64.0d0*sqrt(5005.0d0/pi)*exp(-6.0d0*i*phi)*sina**6*cosa
  if(m==6)  ylm=3.0d0/64.0d0*sqrt(5005.0d0/pi)*exp(6.0d0*i*phi)*sina**6*cosa
  if(m==-5) ylm=3.0d0/64.0d0*sqrt(385.0d0/twopi)*exp(-5.0d0*i*phi)*sina**5*(13.0d0*cosa**2-1.0d0)
  if(m==5)  ylm=-3.0d0/64.0*sqrt(385.0d0/twopi)*exp(5.0d0*i*phi)*sina**5*(13.0d0*cosa**2-1.0d0)
  if(m==-4) ylm=3.0d0/32.0d0*sqrt(385.0d0/twopi)*exp(-4.0d0*i*phi)*sina**4*(13.0d0*cosa**3-3.0d0*cosa)
  if(m==4)  ylm=3.0d0/32.0d0*sqrt(385.0d0/twopi)*exp(4.0d0*i*phi)*sina**4*(13.0d0*cosa**3-3.0d0*cosa)
  if(m==-3) ylm=3.0d0/64.0d0*sqrt(35.0d0/twopi)*exp(-3.0d0*i*phi)*sina**3*(143.0d0*cosa**4-66.0d0*cosa**2+3.0d0)
  if(m==3)  ylm=-3.0d0/64.0d0*sqrt(35.0d0/twopi)*exp(3.0d0*i*phi)*sina**3*(143.0d0*cosa**4-66.0d0*cosa**2+3.0d0)
  if(m==-2) ylm=3.0d0/64.0d0*sqrt(35.0d0/pi)*exp(-2.0d0*i*phi)*sina**2*(143.0d0*cosa**5-110.0d0*cosa**3+15.0d0*cosa)
  if(m==2)  ylm=3.0d0/64.0d0*sqrt(35.0d0/pi)*exp(2.0d0*i*phi)*sina**2*(143.0d0*cosa**5-110.0d0*cosa**3+15.0d0*cosa)
  if(m==-1) ylm=1.0d0/64.0d0*sqrt(105.0d0/twopi)*exp(-i*phi)*sina*(429.0d0*cosa**6-495.0d0*cosa**4+135.0d0*cosa**2-5.0d0)
  if(m==1)  ylm=-1.0d0/64.0d0*sqrt(105.0d0/twopi)*exp(i*phi)*sina*(429.0d0*cosa**6-495.0d0*cosa**4+135.0d0*cosa**2-5.0d0)
  if(m==0)  ylm=1.0d0/32.0d0*sqrt(15.0d0/pi)*(429.0d0*cosa**7-693.0d0*cosa**5+315.0d0*cosa**3-35.0d0*cosa)
     elseif(l==8)then
  if(m==-8) ylm=3.0d0/256.0d0*sqrt(12155.0d0/twopi)*exp(-8.0d0*i*phi)*sina**8
  if(m==8)  ylm=3.0d0/256.0d0*sqrt(12155.0d0/twopi)*exp(8.0d0*i*phi)*sina**8
  if(m==-7) ylm=3.0d0/64.0d0*sqrt(12155.0d0/twopi)*exp(-7.0d0*i*phi)*sina**7*cosa
  if(m==7)  ylm=-3.0d0/64.0d0*sqrt(12155.0d0/twopi)*exp(7.0d0*i*phi)*sina**7*cosa
  if(m==-6) ylm=1.0d0/128.0d0*sqrt(7293.0d0/pi)*exp(-6.0d0*i*phi)*sina**6*(15.0d0*cosa**2-1.0d0)
  if(m==6)  ylm=1.0d0/128.0d0*sqrt(7293.0d0/pi)*exp(6.0d0*i*phi)*sina**6*(15.0d0*cosa**2-1.0d0)
  if(m==-5) ylm=3.0d0/64.0d0*sqrt(17017.0d0/twopi)*exp(-5.0d0*i*phi)*sina**5*(5.0d0*cosa**3-cosa)
  if(m==5)  ylm=-3.0d0/64.0d0*sqrt(17017.0d0/twopi)*exp(5.0d0*i*phi)*sina**5*(5.0d0*cosa**3-cosa)
  if(m==-4) ylm=3.0d0/128.0d0*sqrt(1309.0d0/twopi)*exp(-4.0d0*i*phi)*sina**4*(65.0d0*cosa**4-26.0d0*cosa**2+1.0d0)
  if(m==4)  ylm=3.0d0/128.0d0*sqrt(1309.0d0/twopi)*exp(4.0d0*i*phi)*sina**4*(65.0d0*cosa**4-26.0d0*cosa**2+1.0d0)
  if(m==-3) ylm=1.0d0/64.0d0*sqrt(19635.0d0/twopi)*exp(-3.0d0*i*phi)*sina**3*(39.0d0*cosa**5-26.0d0*cosa**3+3.0d0*cosa)
  if(m==3)  ylm=-1.0d0/64.0d0*sqrt(19635.0d0/twopi)*exp(3.0d0*i*phi)*sina**3*(39.0d0*cosa**5-26.0d0*cosa**3+3.0d0*cosa)
  if(m==-2) ylm=3.0d0/128.0d0*sqrt(595.0d0/pi)*exp(-2.0d0*i*phi)*sina**2*(143.0d0*cosa**6-143.0d0*cosa**4+33.0d0*cosa**2-1.0d0)
  if(m==2)  ylm=3.0d0/128.0d0*sqrt(595.0d0/pi)*exp(2.0d0*i*phi)*sina**2*(143.0d0*cosa**6-143.0d0*cosa**4+33.0d0*cosa**2-1.0d0)
  if(m==-1) ylm=3.0d0/64.0d0*sqrt(17.0d0/twopi)*exp(-i*phi)*sina*(715.0d0*cosa**7-1001.0d0*cosa**5+385.0d0*cosa**3-35.0d0*cosa)
  if(m==1)  ylm=-3.0d0/64.0d0*sqrt(17.0d0/twopi)*exp(i*phi)*sina*(715.0d0*cosa**7-1001.0d0*cosa**5+385.0d0*cosa**3-35.0d0*cosa)
  if(m==0)  ylm=1.0d0/256.0d0*sqrt(17.0d0/pi)*(6435.0d0*cosa**8-12012.0d0*cosa**6+6930.0d0*cosa**4-1260.0d0*cosa**2+35.0d0)
     elseif(l==9)then
  if(m==-9) ylm=1.0d0/512.0d0*sqrt(230945.0d0/pi)*exp(-9.0d0*i*phi)*sina**9
  if(m==9)  ylm=-1.0d0/512.0d0*sqrt(230945.0d0/pi)*exp(9.0d0*i*phi)*sina**9  
  if(m==-8) ylm=3.0d0/256.0d0*sqrt(230945.0d0/twopi)*exp(-8.0d0*i*phi)*sina**8*cosa
  if(m==8)  ylm=3.0d0/256.0d0*sqrt(230945.0d0/twopi)*exp(8.0d0*i*phi)*sina**8.0*cosa
  if(m==-7) ylm=3.0d0/512.0d0*sqrt(13585.0d0/pi)*exp(-7.0d0*i*phi)*sina**7*(17.0d0*cosa**2-1.0d0)
  if(m==7)  ylm=-3.0d0/512.0d0*sqrt(13585.0d0/pi)*exp(7.0d0*i*phi)*sina**7*(17.0d0*cosa**2-1.0d0)
  if(m==-6) ylm=1.0d0/128.0d0*sqrt(40755.0d0/pi)*exp(-6.0d0*i*phi)*sina**6*(17.0d0*cosa**3-3.0d0*cosa)
  if(m==6)  ylm=1.0d0/128.0d0*sqrt(40755.0d0/pi)*exp(6.0d0*i*phi)*sina**6*(17.0d0*cosa**3-3.0d0*cosa)
  if(m==-5) ylm=3.0d0/256.0d0*sqrt(2717.0d0/pi)*exp(-5.0d0*i*phi)*sina**5*(85.0d0*cosa**4-30.0d0*cosa**2+1.0d0)
  if(m==5)  ylm=-3.0d0/256.0d0*sqrt(2717.0d0/pi)*exp(5.0d0*i*phi)*sina**5*(85.0d0*cosa**4-30.0d0*cosa**2+1.0d0)
  if(m==-4) ylm=3.0d0/128.0d0*sqrt(95095.0d0/twopi)*exp(-4.0d0*i*phi)*sina**4*(17.0d0*cosa**5-10.0d0*cosa**3+cosa)
  if(m==4)  ylm=3.0d0/128.0d0*sqrt(95095.0d0/twopi)*exp(4.0d0*i*phi)*sina**4*(17.0d0*cosa**5-10.0d0*cosa**3+cosa)
  if(m==-3) ylm=1.0d0/256.0d0*sqrt(21945.0d0/pi)*exp(-3.0d0*i*phi)*sina**3*(221.0d0*cosa**6-195.0d0*cosa**4+39.0d0*cosa**2-1.0d0)
  if(m==3)  ylm=-1.0d0/256.0d0*sqrt(21945.0d0/pi)*exp(3.0d0*i*phi)*sina**3*(221.0d0*cosa**6-195.0d0*cosa**4+39.0d0*cosa**2-1.0d0)
  if(m==-2) ylm=3.0d0/128.0d0*sqrt(1045.0d0/pi)*exp(-2.0d0*i*phi)*sina**2*(221.0d0*cosa**7-273.0d0*cosa**5+91.0d0*cosa**3-7.0d0*cosa)
  if(m==2)  ylm=3.d0/128.0d0*sqrt(1045.0d0/pi)*exp(2.0d0*i*phi)*sina**2*(221.0d0*cosa**7-273.0d0*cosa**5+91.0d0*cosa**3-7.0d0*cosa)
  if(m==-1) ylm=3.0d0/256.0d0*sqrt(95.0d0/twopi)*exp(-i*phi)*sina*(2431.0d0*cosa**8-4004.0d0*cosa**6+2002.0d0*cosa**4-308.0d0*cosa**2+7.0d0)
  if(m==1)  ylm=-3.0d0/256.0d0*sqrt(95.0d0/twopi)*exp(i*phi)*sina*(2431.0d0*cosa**8-4004.0d0*cosa**6+2002.0d0*cosa**4-308.0d0*cosa**2+7.0d0)
  if(m==0)  ylm=1.0d0/256.0d0*sqrt(19.0d0/pi)*(12155.0d0*cosa**9-25740.0d0*cosa**7+18018.0d0*cosa**5-4620.0d0*cosa**3+315.0d0*cosa)  
    elseif(l==10)then
  if(m==-10) ylm=1.0d0/1024.0d0*sqrt(969969.0d0/pi)*exp(-10.0d0*i*phi)*sina**10
  if(m==10)  ylm=1.0d0/1024.0d0*sqrt(969969.0d0/pi)*exp(10.0d0*i*phi)*sina**10 
  if(m==-9)  ylm=1.0d0/512.0d0*sqrt(4849845.0d0/pi)*exp(-9.0d0*i*phi)*sina**9*cosa
  if(m==9)   ylm=-1.0d0/512.0d0*sqrt(4849845.0d0/pi)*exp(9.0d0*i*phi)*sina**9*cosa
  if(m==-8)  ylm=1.0d0/512.0d0*sqrt(255255.0d0/twopi)*exp(-8.0d0*i*phi)*sina**8*(19.0d0*cosa**2-1.0d0)
  if(m==8)   ylm=1.0d0/512.0d0*sqrt(255255.0d0/twopi)*exp(8.0d0*i*phi)*sina**8*(19.0d0*cosa**2-1.0d0)
  if(m==-7)  ylm=3.0d0/512.0d0*sqrt(85085.0d0/pi)*exp(-7.0d0*i*phi)*sina**7*(19.0d0*cosa**3-3.0d0*cosa)
  if(m==7)   ylm=-3.0d0/512.0d0*sqrt(85085.0d0/pi)*exp(7.0d0*i*phi)*sina**7*(19.0d0*cosa**3-3.0d0*cosa)
  if(m==-6)  ylm=3.0d0/1024.0d0*sqrt(5005.0d0/pi)*exp(-6.0d0*i*phi)*sina**6*(323.0d0*cosa**4-102.0d0*cosa**2+3.0d0)
  if(m==6)   ylm=3.0d0/1024.0d0*sqrt(5005.0d0/pi)*exp(6.0d0*i*phi)*sina**6*(323.0d0*cosa**4-102.0d0*cosa**2+3.0d0)
  if(m==-5)  ylm=3.0d0/256.0d0*sqrt(1001.0d0/pi)*exp(-5.0d0*i*phi)*sina**5*(323.0d0*cosa**5-170.0*cosa**3+15.0d0*cosa)
  if(m==5)   ylm=-3.0d0/256.0d0*sqrt(1001.0d0/pi)*exp(5.0d0*i*phi)*sina**5*(323.0d0*cosa**5-170.0d0*cosa**3+15.0d0*cosa)
  if(m==-4)  ylm=3.0d0/256.0d0*sqrt(5005.0d0/twopi)*exp(-4.0d0*i*phi)*sina**4*(323.0d0*cosa**6-255.0d0*cosa**4+45.0d0*cosa**2-1.0d0)
  if(m==4)   ylm=3.0d0/256.0d0*sqrt(5005.0d0/twopi)*exp(4.0d0*i*phi)*sina**4*(323.0d0*cosa**6-255.0d0*cosa**4+45.0*cosa**2-1.0d0)
  if(m==-3)  ylm=3.0d0/256.0d0*sqrt(5005.0d0/pi)*exp(-3.0d0*i*phi)*sina**3*(323.0d0*cosa**7-357.0*cosa**5+105.0d0*cosa**3-7.0*cosa)
  if(m==3)   ylm=-3.0d0/256.0d0*sqrt(5005.0d0/pi)*exp(3.0d0*i*phi)*sina**3*(323.0d0*cosa**7-357.0d0*cosa**5+105.0d0*cosa**3-7.0d0*cosa)
  if(m==-2)  ylm=3.0d0/512.0d0*sqrt(385.0d0/twopi)*exp(-2.0d0*i*phi)*sina**2*(4199.0d0*cosa**8-6188.0d0*cosa**6+2730.0d0*cosa**4-364.0d0*cosa**2+7.0d0)
  if(m==2)   ylm=3.0d0/512.0d0*sqrt(385.0d0/twopi)*exp(2.0d0*i*phi)*sina**2*(4199.0d0*cosa**8-6188.0d0*cosa**6+2730.0d0*cosa**4-364.0d0*cosa**2+7.0d0)
  if(m==-1)  ylm=1.0d0/256.0d0*sqrt(1155.0d0/twopi)*exp(-i*phi)*sina*(4199.0d0*cosa**9-7956.0*cosa**7+4914.0d0*cosa**5-1092.0*cosa**3+63.0d0*cosa)
  if(m==1)   ylm=-1.0d0/256.0d0*sqrt(1155.0d0/twopi)*exp(i*phi)*sina*(4199.0d0*cosa**9-7956.0*cosa**7+4914.0d0*cosa**5-1092.0*cosa**3+63.0d0*cosa) 
  if(m==0)   ylm=1.0d0/512.0d0*sqrt(21.0d0/pi)*(46189.0d0*cosa**10-109395.0d0*cosa**8+90090.0d0*cosa**6-30030.0d0*cosa**4+3465.0d0*cosa**2-63.0d0)
                 endif
       end subroutine sphhar
!234567890
