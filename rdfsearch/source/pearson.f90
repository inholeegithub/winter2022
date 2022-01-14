!234567890
      subroutine pearsn(x,y,n,r,prob,z)
      integer n
      real*8 prob,r,z,x(n),y(n),tiny
      parameter (tiny=1.d-20)
!u    uses betai
      integer j
      real*8 ax,ay,df,sxx,sxy,syy,t,xt,yt,betai

      ax=0.d0
      ay=0.d0
      do 11 j=1,n
        ax=ax+x(j)
        ay=ay+y(j)
 11   continue
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
 12   continue
      r=sxy/sqrt(sxx*syy)
      z=0.5d0*log(((1.d0+r)+TINY)/((1.d0-r)+TINY))
      df=n-2
      t=r*sqrt(df/(((1.-r)+TINY)*((1.+r)+TINY)))
      if(.false.)then
      prob=betai(0.5*df,0.5d0,df/(df+t**2))
!     prob=erfcc(abs(z*sqrt(n-1.d0))/1.414213562d0)
                 endif
      end
      SUBROUTINE KENDL1(DATA1,DATA2,N,TAU,Z,PROB)
      implicit none
      integer n
      real*8 DATA1(N),DATA2(N),tau,z,prob
      real*8 var,a1,a2,aa
      integer is,n1,n2,j,k
      real*8, external :: erfcc

      N1=0
      N2=0
      IS=0
      DO 12 J=1,N-1
        DO 11 K=J+1,N
          A1=DATA1(J)-DATA1(K)
          A2=DATA2(J)-DATA2(K)
          AA=A1*A2
          IF(AA.NE.0.)THEN
            N1=N1+1
            N2=N2+1
            IF(AA.GT.0.)THEN
              IS=IS+1
            ELSE
              IS=IS-1
            ENDIF
          ELSE
            IF(A1.NE.0.)N1=N1+1
            IF(A2.NE.0.)N2=N2+1
          ENDIF
  11    CONTINUE
  12  CONTINUE
      TAU=FLOAT(IS)/SQRT(FLOAT(N1)*FLOAT(N2))
      VAR=(4.*N+10.)/(9.*N*(N-1.))
      Z=TAU/SQRT(VAR)
      PROB=ERFCC(ABS(Z)/1.4142136)
      RETURN
      END
      SUBROUTINE KENDL2(TAB,I,J,IP,JP,TAU,Z,PROB)
      implicit none
      integer ip,jp,i,j
      real*8 TAB(IP,JP),tau,z,prob
      integer nn,k,ki,kj,m1,m2,mm,li,lj,l
      real*8 en1,en2,s,points,pairs,var
      real*8, external :: erfcc

      EN1=0.
      EN2=0.
      S=0.
      NN=I*J
      POINTS=TAB(I,J)
      DO 12 K=0,NN-2
        KI=K/J
        KJ=K-J*KI
        POINTS=POINTS+TAB(KI+1,KJ+1)
        DO 11 L=K+1,NN-1
          LI=L/J
          LJ=L-J*LI
          M1=LI-KI
          M2=LJ-KJ
          MM=M1*M2
          PAIRS=TAB(KI+1,KJ+1)*TAB(LI+1,LJ+1)
          IF(MM.NE.0)THEN
            EN1=EN1+PAIRS
            EN2=EN2+PAIRS
            IF(MM.GT.0)THEN
              S=S+PAIRS
            ELSE
              S=S-PAIRS
            ENDIF
          ELSE
            IF(M1.NE.0)EN1=EN1+PAIRS
            IF(M2.NE.0)EN2=EN2+PAIRS
          ENDIF
  11    CONTINUE
  12  CONTINUE
      TAU=S/SQRT(EN1*EN2)
      VAR=(4.*POINTS+10.)/(9.*POINTS*(POINTS-1.))
      Z=TAU/SQRT(VAR)
      PROB=ERFCC(ABS(Z)/1.4142136)
      RETURN
      END
      SUBROUTINE SPEAR(DATA1,DATA2,N,WKSP1,WKSP2,D,ZD,PROBD,RS,PROBRS)
      implicit none
      integer n
      real*8 DATA1(N),DATA2(N),WKSP1(N),WKSP2(N)
      real*8 d,zd,probd,rs,probrs
      integer j
      real*8 t,en,en3n,aved,fac,vard,df,sf,sg
      real*8, external :: erfcc,betai

      DO 11 J=1,N
        WKSP1(J)=DATA1(J)
        WKSP2(J)=DATA2(J)
  11  CONTINUE
      CALL SORT2(N,WKSP1,WKSP2)
      CALL CRANK(N,WKSP1,SF)
      CALL SORT2(N,WKSP2,WKSP1)
      CALL CRANK(N,WKSP2,SG)
      D=0.
      DO 12 J=1,N
        D=D+(WKSP1(J)-WKSP2(J))**2
  12  CONTINUE
      EN=N
      EN3N=EN**3-EN
      AVED=EN3N/6.-(SF+SG)/12.
      FAC=(1.-SF/EN3N)*(1.-SG/EN3N)
      VARD=((EN-1.)*EN**2*(EN+1.)**2/36.)*FAC
      ZD=(D-AVED)/SQRT(VARD)
      PROBD=ERFCC(ABS(ZD)/1.4142136)
      RS=(1.-(6./EN3N)*(D+0.5*(SF+SG)))/FAC
      T=RS*SQRT((EN-2.)/((1.+RS)*(1.-RS)))
      DF=EN-2.
      PROBRS=BETAI(0.5*DF,0.5,DF/(DF+T**2))
      RETURN
      END
      SUBROUTINE SORT(N,RA)
      implicit none
      integer n
      real*8 RA(N)
      integer l,ir,j,i
      real*8 rra

      L=N/2+1
      IR=N
  10  CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
  20    IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END
      SUBROUTINE SORT2(N,RA,RB)
      implicit none
      integer n
      real*8 RA(N),RB(N)
      integer l,ir,j,i
      real*8 rra,rrb

      L=N/2+1
      IR=N
  10  CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
        ELSE
          RRA=RA(IR)
          RRB=RB(IR)
          RA(IR)=RA(1)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
  20    IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
      GO TO 10
      END
      SUBROUTINE CRANK(N,W,S)
      implicit none
      integer n
      real*8 W(N)
      real*8 s
      integer j,jt,ji
      real*8 t,rank

      S=0.
      J=1
  1   IF(J.LT.N)THEN
        IF(W(J+1).NE.W(J))THEN
          W(J)=J
          J=J+1
        ELSE
          DO 11 JT=J+1,N
            IF(W(JT).NE.W(J))GO TO 2
  11      CONTINUE
          JT=N+1
  2       RANK=0.5*(J+JT-1)
          DO 12 JI=J,JT-1
            W(JI)=RANK
  12      CONTINUE
          T=JT-J
          S=S+T**3-T
          J=JT
        ENDIF
      GO TO 1
      ENDIF
      IF(J.EQ.N)W(N)=N
      RETURN
      END
      SUBROUTINE SORT3(N,RA,RB,RC,WKSP,IWKSP)
      implicit none
      integer n
      integer IWKSP(N)
      real*8 RA(N),RB(N),RC(N),WKSP(N)
      integer j

      CALL INDEXX(N,RA,IWKSP)
      DO 11 J=1,N
        WKSP(J)=RA(J)
  11  CONTINUE
      DO 12 J=1,N
        RA(J)=WKSP(IWKSP(J))
  12  CONTINUE
      DO 13 J=1,N
        WKSP(J)=RB(J)
  13  CONTINUE
      DO 14 J=1,N
        RB(J)=WKSP(IWKSP(J))
  14  CONTINUE
      DO 15 J=1,N
        WKSP(J)=RC(J)
  15  CONTINUE
      DO 16 J=1,N
        RC(J)=WKSP(IWKSP(J))
  16  CONTINUE
      RETURN
      END
      SUBROUTINE INDEXX(N,ARRIN,INDX)
      implicit none
      integer n
      integer INDX(N)
      real*8 ARRIN(N)
      integer j,ir,l,i,indxt
      real*8 q

      DO 11 J=1,N
        INDX(J)=J
  11  CONTINUE
      L=N/2+1
      IR=N
  10  CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
  20    IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END
      FUNCTION ERFCC(X)
      implicit none
      real*8 x
      real*8 erfcc
      real*8 z,t

      Z=ABS(X)      
      T=1./(1.+0.5*Z)
      ERFCC=T*EXP(-Z*Z-1.26551223+T*(1.00002368+T*(.37409196+ &
          T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+ &
          T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
      IF (X.LT.0.) ERFCC=2.-ERFCC
      RETURN
      END
!  (C) Copr. 1986-92 Numerical Recipes Software .37. p 632
      function betai(a,b,x)
      real*8 betai,a,b,x
!U    USES betacf,gammln
      real*8 bt,betacf,gammln

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
      end
!  (C) Copr. 1986-92 Numerical Recipes Software .37.
      function gammln(xx)
      real*8 gammln,xx
      integer j
      real*8 ser,stp,tmp,x,y,cof(6)
      save cof,stp
      data cof,stp/76.18009172947146d0,-86.50532032941677d0, &
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
 11   continue
      gammln=tmp+log(stp*ser/x)
      return
      end
!  (C) Copr. 1986-92 Numerical Recipes Software .37.
      function betacf(a,b,x)
      integer maxit
      real*8 betacf,a,b,x,EPS,FPMIN
      parameter (maxit=100,eps=3.d-7,fpmin=1.d-30)
      integer m,m2
      real*8 aa,c,d,del,h,qab,qam,qap

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
 11   continue
!     pause 'a or b too big, or MAXIT too small in betacf'
      write(6,*) 'a or b too big, or MAXIT too small in betacf'
  1   betacf=h
      return
      end
!  (C) Copr. 1986-92 Numerical Recipes Software .37. 
      subroutine pearson(xdata,ydata,ndata,rr)
      implicit none
      integer ndata
      real*8 xdata(ndata),ydata(ndata),rr
      real*8 prob,zz 

      call pearsn(xdata,ydata,ndata,rr,prob,zz)
!     write(6,'(f8.4)') rr
!     write(6,*) rr,rr**2, ' r, r^2'
!     write(6,*) prob, ' prob'
!     write(6,*) zz, ' z'
      end
!234567890
