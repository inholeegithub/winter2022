      SUBROUTINE pearsn(x,y,n,r,prob,z)
      INTEGER n
      REAL*8 prob,r,z,x(n),y(n),TINY
!     PARAMETER (TINY=1.d-20)
      PARAMETER (TINY=1.d-15)
!U    USES betai
      INTEGER j
      REAL*8 ax,ay,df,sxx,sxy,syy,t,xt,yt,betai
      ax=0.d0
      ay=0.d0
      do 11 j=1,n
        ax=ax+x(j)
        ay=ay+y(j)
11    continue
      ax=ax/n
      ay=ay/n
      sxx=0.d0
      syy=0.d0
      sxy=0.d0
      do 12 j=1,n
        xt=x(j)-ax
        yt=y(j)-ay
        sxx=sxx+xt**2
        syy=syy+yt**2
        sxy=sxy+xt*yt
12    continue
      r=sxy/sqrt(sxx*syy)
      z=0.5d0*log(((1.d0+r)+TINY)/((1.d0-r)+TINY))
      df=n-2
      t=r*sqrt(df/(((1.-r)+TINY)*((1.+r)+TINY)))
      prob=betai(0.5*df,0.5d0,df/(df+t**2))
!     prob=erfcc(abs(z*sqrt(n-1.d0))/1.414213562d0)
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software .37. p 632
      FUNCTION betai(a,b,x)
      REAL*8 betai,a,b,x
!U    USES betacf,gammln
      REAL*8 bt,betacf,gammln
      if(x.lt.0..or.x.gt.1.)then
      write(6,*) 'bad argument x in betai'
                            stop
                            endif
      if(x.eq.0..or.x.eq.1.)then
        bt=0.
      else
        bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.d0-x))
      endif
      if(x.lt.(a+1.)/(a+b+2.))then
        betai=bt*betacf(a,b,x)/a
        return
      else
        betai=1.-bt*betacf(b,a,1.d0-x)/b
        return
      endif
      END
!  (C) Copr. 1986-92 Numerical Recipes Software .37.
      FUNCTION gammln(xx)
      REAL*8 gammln,xx
      INTEGER j
      real*8 ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, &
      24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
      -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software .37.
      FUNCTION betacf(a,b,x)
      INTEGER MAXIT
      REAL*8 betacf,a,b,x,EPS,FPMIN
      PARAMETER (MAXIT=100,EPS=3.d-7,FPMIN=1.d-30)
      INTEGER m,m2
      REAL*8 aa,c,d,del,h,qab,qam,qap
      qab=a+b
      qap=a+1.
      qam=a-1.
      c=1.
      d=1.-qab*x/qap
      if(abs(d).lt.FPMIN)d=FPMIN
      d=1./d
      h=d
      do 11 m=1,MAXIT
        m2=2*m
        aa=m*(b-m)*x/((qam+m2)*(a+m2))
        d=1.+aa*d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        h=h*d*c
        aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d=1.+aa*d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
  11  continue
!     pause 'a or b too big, or MAXIT too small in betacf'
      write(6,*) 'a or b too big, or MAXIT too small in betacf'
   1  betacf=h
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software .37.

      program abc
      implicit none
      integer n
      real*8 x(0:180),y(0:180)
      real*8 r,prob,z
      integer i
      do i=0,180
      x(i)=exp(-(i-60.)**2)
      y(i)=2.*exp(-(i-60.)**2)
      enddo
      n=180+1
      call pearsn(x,y,n,r,prob,z)
      write(6,*) r,prob,z
      end program abc

