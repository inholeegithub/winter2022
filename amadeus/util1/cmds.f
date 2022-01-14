C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                         C
C  Carry out a PRINCIPAL COORDINATES ANALYSIS
C                                                                         C
C              (CLASSICAL MULTIDIMENSIONAL SCALING).
C                                                                         C
C                                                                         C
C  To call:   CALL PCA(N,A,IPRINT,W1,W2,A2,IERR)    where
C                                                                         C
C                                                                         C
C                                                                         C
C  N     : integer number of objects.
C                                                                         C
C  A     : input distances array, real, of dimensions N by N.
C                                                                         C
C          On output, A contains in the first 7 columns the projections
C                                                                         C
C          of the objects on the first 7 principal components.
C                                                                         C
C  IPRINT: print options.
C                                                                         C
C          = 3: full printing of items calculated.
C                                                                         C
C          = 2: printing of everything except the input distance matrix.
C                                                                         C
C          Otherwise: no printing.
C                                                                         C
C  W1,W2 : real vectors of dimension M (see called routines for use).
C                                                                         C
C          On output, W1 contains the eigenvalues (in increasing order
C          of                                                             C
C          magnitude).
C                                                                         C
C  A2    : real array of dimensions M * M (see called routines for use).
C                                                                         C
C  IERR  : error indicator (normally zero).
C                                                                         C
C                                                                         C
C                                                                         C
C  Inputs here are N, A, IPRINT (and IERR).
C                                                                         C
C  Output information is contained in A, and W1.
C                                                                         C
C  All printed outputs are carried out in easily recognizable
C  subroutines                                                            C
C  called from the first subroutine following.
C                                                                         C
C                                                                         C
C                                                                         C
C  F. Murtagh, ST-ECF/ESA/ESO, Garching-bei-Muenchen, January 1986.
C                                                                         C
C                                                                         C
C-------------------------------------------------------------------------C
        SUBROUTINE CMDS(N,A,IPRINT,W,FV1,Z,IERR)
        REAL    A(N,N), W(N), FV1(N), Z(N,N)
C
        IF (IPRINT.EQ.3) CALL OUTMAT(N,A)
C
        TOT = 0.0
        DO 200 I1 = 1, N
           W(I1) = 0.0
           DO 100 I2 = 1, N
              W(I1) = W(I1) + A(I1,I2)
              TOT = TOT + A(I1,I2)
  100      CONTINUE
           W(I1) = W(I1)/FLOAT(N)
  200   CONTINUE
        TOT = TOT/(FLOAT(N)*FLOAT(N))
C
        DO 300 I1 = 1, N
           DO 300 I2 = 1, N
              A(I1,I2) = -0.5 * (A(I1,I2) - W(I1) - W(I2) - TOT)
  300   CONTINUE
C
C          Carry out the eigenreduction.
C
        N2 = N
        CALL TRED2(N,N2,A,W,FV1,Z)
        CALL TQL2(N,N2,W,FV1,Z,IERR)
        IF (IERR.NE.0) GOTO 9000
C
C          Output eigenvalues and eigenvectors.
C
        IF (IPRINT.GE.2) CALL OUTEVL(N,W)
        IF (IPRINT.GE.2) CALL OUTEVC(N,Z)
C
C          Determine projections and output them.
C
        CALL PROJN(N,W,A,Z,FV1)
        IF (IPRINT.GE.2) CALL OUTPR(N,A)
C
 9000   RETURN  
        END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                         C
C Reduce a real, symmetric matrix to a symmetric, tridiagonal matrix.
C                                                                         C
C                                                                         C
C To call:    CALL TRED2(NM,N,A,D,E,Z)    where
C                                                                         C
C                                                                         C 
C NM = row dimension of A and Z;
C                                                                         C
C N = order of matrix A (will always be <= NM);
C                                                                         C
C A = symmetric matrix of order N to be reduced to tridiagonal form;
C                                                                         C
C D = vector of dim. N containing, on output, diagonal elts. of trid.
C                                                                         C
C     matrix.
C                                                                         C
C E = working vector of dim. at least N-1 to contain subdiagonal elts.;
C                                                                         C
C Z = matrix of dims. NM by N containing, on output, orthogonal
C                                                                         C
C     transformation matrix producing the reduction.
C                                                                         C
C                                                                         C
C Normally a call to TQL2 will follow the call to TRED2 in order to
C                                                                         C
C produce all eigenvectors and eigenvalues of matrix A.
C                                                                         C
C                                                                         C
C Algorithm used: Martin et al., Num. Math. 11, 181-195, 1968.
C                                                                         C
C                                                                         C
C Reference: Smith et al., Matrix Eigensystem Routines - EISPACK
C                                                                         C
C Guide, Lecture Notes in Computer Science 6, Springer-Verlag, 1976,
C                                                                         C
C pp. 489-494.
C                                                                         C
C                                                                         C
C F. Murtagh, ST-ECF/ESA/ESO, Garching-bei-Muenchen, January 1986.
C                                                                         C
C                                                                         C
C-------------------------------------------------------------------------C
        SUBROUTINE TRED2(NM,N,A,D,E,Z)
        REAL A(NM,N),D(N),E(N),Z(NM,N)
C
        DO 100 I = 1, N
           DO 100 J = 1, I
              Z(I,J) = A(I,J)
  100   CONTINUE
        IF (N.EQ.1) GOTO 320
        DO 300 II = 2, N
           I = N + 2 - II
           L = I - 1
           H = 0.0
           SCALE = 0.0
           IF (L.LT.2) GOTO 130
           DO 120 K = 1, L
              SCALE = SCALE + ABS(Z(I,K))
  120      CONTINUE
           IF (SCALE.NE.0.0) GOTO 140
  130      E(I) = Z(I,L)
           GOTO 290
  140      DO 150 K = 1, L
              Z(I,K) = Z(I,K)/SCALE
              H = H + Z(I,K)*Z(I,K)
  150      CONTINUE
C
           F = Z(I,L)
           G = -SIGN(SQRT(H),F)
           E(I) = SCALE * G
           H = H - F * G
           Z(I,L) = F - G
           F = 0.0
C
           DO 240 J = 1, L
              Z(J,I) = Z(I,J)/H
              G = 0.0
C             Form element of A*U.
              DO 180 K = 1, J
                 G = G + Z(J,K)*Z(I,K)
  180         CONTINUE
              JP1 = J + 1
              IF (L.LT.JP1) GOTO 220
              DO 200 K = JP1, L
                 G = G + Z(K,J)*Z(I,K)
  200         CONTINUE
C             Form element of P where P = I - U U' / H .
  220         E(J) = G/H
              F = F + E(J) * Z(I,J)
  240      CONTINUE
           HH = F/(H + H)
C          Form reduced A.
           DO 260 J = 1, L
              F = Z(I,J)
              G = E(J) - HH * F
              E(J) = G
              DO 250 K = 1, J
                 Z(J,K) = Z(J,K) - F*E(K) - G*Z(I,K)
  250         CONTINUE
  260      CONTINUE
  290      D(I) = H
  300   CONTINUE
  320   D(1) = 0.0
        E(1) = 0.0
C       Accumulation of transformation matrices.
        DO 500 I = 1, N
           L = I - 1
           IF (D(I).EQ.0.0) GOTO 380
           DO 360 J = 1, L
              G = 0.0
              DO 340 K = 1, L
                 G = G + Z(I,K) * Z(K,J)
  340         CONTINUE
              DO 350 K = 1, L
                 Z(K,J) = Z(K,J) - G * Z(K,I)
  350         CONTINUE
  360      CONTINUE
  380      D(I) = Z(I,I)
           Z(I,I) = 1.0
           IF (L.LT.1) GOTO 500
           DO 400 J = 1, L
              Z(I,J) = 0.0
              Z(J,I) = 0.0
  400      CONTINUE
  500   CONTINUE
C
        RETURN
        END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                         C
C Determine eigenvalues and eigenvectors of a symmetric,
C                                                                         C 
C tridiagonal matrix.
C                                                                         C
C                                                                         C
C To call:    CALL TQL2(NM,N,D,E,Z,IERR)    where
C                                                                         C
C                                                                         C
C NM = row dimension of Z;
C                                                                         C
C N = order of matrix Z;
C                                                                         C
C D = vector of dim. N containing, on output, eigenvalues;
C                                                                         C
C E = working vector of dim. at least N-1;
C                                                                         C
C Z = matrix of dims. NM by N containing, on output, eigenvectors;
C                                                                         C
C IERR = error, normally 0, but 1 if no convergence.
C                                                                         C
C                                                                         C
C Normally the call to TQL2 will be preceded by a call to TRED2 in
C                                                                         C
C order to set up the tridiagonal matrix.
C                                                                         C
C                                                                         C
C Algorithm used: QL method of Bowdler et al., Num. Math. 11,
C                                                                         C
C 293-306, 1968.
C                                                                         C
C                                                                         C
C Reference: Smith et al., Matrix Eigensystem Routines - EISPACK
C                                                                         C
C Guide, Lecture Notes in Computer Science 6, Springer-Verlag, 1976,
C                                                                         C
C pp. 468-474.
C                                                                         C
C                                                                         C
C F. Murtagh, ST-ECF/ESA/ESO, Garching-bei-Muenchen, January 1986.
C                                                                         C
C                                                                         C
C-------------------------------------------------------------------------C
        SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
        REAL    D(N), E(N), Z(NM,N)
        DATA    EPS/1.E-12/
C
        IERR = 0
        IF (N.EQ.1) GOTO 1001
        DO 100 I = 2, N
           E(I-1) = E(I)
  100   CONTINUE
        F = 0.0
        B = 0.0
        E(N) = 0.0
C
        DO 240 L = 1, N
           J = 0
           H = EPS * (ABS(D(L)) + ABS(E(L)))
           IF (B.LT.H) B = H
C          Look for small sub-diagonal element.
           DO 110 M = L, N
              IF (ABS(E(M)).LE.B) GOTO 120
C             E(N) is always 0, so there is no exit through the 
C             bottom of the loop.
  110      CONTINUE
  120      IF (M.EQ.L) GOTO 220
  130      IF (J.EQ.30) GOTO 1000
           J = J + 1
C          Form shift.
           L1 = L + 1
           G = D(L)
           P = (D(L1)-G)/(2.0*E(L))
           R = SQRT(P*P+1.0)
           D(L) = E(L)/(P+SIGN(R,P))
           H = G-D(L)
C
           DO 140 I = L1, N
              D(I) = D(I) - H
  140      CONTINUE
C
           F = F + H
C          QL transformation.
           P = D(M)
           C = 1.0
           S = 0.0
           MML = M - L
C
           DO 200 II = 1, MML
              I = M - II
              G = C * E(I)
              H = C * P
              IF (ABS(P).LT.ABS(E(I))) GOTO 150
              C = E(I)/P
              R = SQRT(C*C+1.0)
              E(I+1) = S * P * R
              S = C/R
              C = 1.0/R
              GOTO 160
  150         C = P/E(I)
              R = SQRT(C*C+1.0)
              E(I+1) = S * E(I) * R
              S = 1.0/R
              C = C * S
  160         P = C * D(I) - S * G
              D(I+1) = H + S * (C * G + S * D(I))
C             Form vector.
              DO 180 K = 1, N
                 H = Z(K,I+1)
                 Z(K,I+1) = S * Z(K,I) + C * H
                 Z(K,I) = C * Z(K,I) - S * H
  180         CONTINUE
  200      CONTINUE
           E(L) = S * P
           D(L) = C * P
           IF (ABS(E(L)).GT.B) GOTO 130
  220      D(L) = D(L) + F
  240   CONTINUE
C
C       Order eigenvectors and eigenvalues.
        DO 300 II = 2, N
           I = II - 1
           K = I
           P = D(I)
           DO 260 J = II, N
              IF (D(J).GE.P) GOTO 260
              K = J
              P = D(J)
  260      CONTINUE
           IF (K.EQ.I) GOTO 300
           D(K) = D(I)
           D(I) = P
           DO 280 J = 1, N
              P = Z(J,I)
              Z(J,I) = Z(J,K)
              Z(J,K) = P
  280      CONTINUE
  300   CONTINUE
C
        GOTO 1001
C       Set error - no convergence to an eigenvalue after 30 iterations.
 1000   IERR = L
 1001   RETURN
        END  

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                         C
C  Output array.
C                                                                         C
C                                                                         C
C-------------------------------------------------------------------------C
        SUBROUTINE OUTMAT(N,ARRAY)
        DIMENSION ARRAY(N,N)
C
        DO 100 K1 = 1, N
           WRITE (6,1000) (ARRAY(K1,K2),K2=1,N)
  100   CONTINUE
C
 1000   FORMAT(10(2X,F8.4))
        RETURN
        END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                         C
C  Output eigenvalues in order of decreasing value.
C                                                                         C
C  Ignore first (trivial) eigenvalue.
C                                                                         C
C                                                                         C
C-------------------------------------------------------------------------C
        SUBROUTINE OUTEVL(NVALS,VALS)
        DIMENSION       VALS(NVALS)
C
        TOT = 0.0
        DO 100 K = 1, NVALS-1
           TOT = TOT + VALS(K)
  100   CONTINUE
C
        WRITE (6,1000)
        CUM = 0.0
        K = NVALS 
        WRITE (6,1010)
        WRITE (6,1020)
  200   CONTINUE
        K = K - 1
        CUM = CUM + VALS(K)
        VPC = VALS(K) * 100.0 / TOT
        VCPC = CUM * 100.0 / TOT
        WRITE (6,1030) VALS(K),VPC,VCPC
        IF (K.GT.1) GOTO 200
C
        RETURN
 1000   FORMAT(1H0,'EIGENVALUES FOLLOW.',/)
 1010   FORMAT
     X  (' Eigenvalues       As Percentages    Cumul. Percentages')
 1020   FORMAT
     X  (' -----------       --------------    ------------------')
 1030   FORMAT(F10.4,9X,F10.4,10X,F10.4)
        END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                         C
C         Output FIRST SEVEN eigenvectors associated with eigenvalues
C                                                                         C
C         in descending order.
C                                                                         C
C         Ignore first (trivial) eigenvector.
C                                                                         C
C                                                                         C
C-------------------------------------------------------------------------C
        SUBROUTINE OUTEVC(NDIM,VECS)
        DIMENSION       VECS(NDIM,NDIM)
C
        NUM = MIN0(NDIM,7)
        WRITE (6,1000)
        WRITE (6,1010)
        WRITE (6,1020)
        DO 100 K1 = 1, NDIM
        WRITE (6,1030) K1,(VECS(K1,NDIM-K2),K2=1,NUM)
  100   CONTINUE
C
        RETURN
 1000   FORMAT(1H0,'EIGENVECTORS FOLLOW.',/)
 1010   FORMAT('  VBLE.   EV-1    EV-2    EV-3    EV-4    EV-5    EV-6 
     X   EV-7')
 1020   FORMAT(' ------  ------  ------  ------  ------  ------  ------  
     X------')
 1030   FORMAT(I5,2X,7F8.4)
        END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                         C
C  Output projections of objects on first 7 principal components.
C                                                                         C
C                                                                         C
C-------------------------------------------------------------------------C
        SUBROUTINE OUTPR(N,PRJNS)
        REAL    PRJNS(N,N)
C
        NUM = MIN0(N,7)
        WRITE (6,1000)
        WRITE (6,1010)
        WRITE (6,1020)
        DO 100 K = 1, N
           WRITE (6,1030) K,(PRJNS(K,J),J=1,NUM)
  100   CONTINUE
C
 1000   FORMAT(1H0,'PROJECTIONS OF OBJECTS FOLLOW.',/)
 1010   FORMAT('  VBLE.  PROJ-1  PROJ-2  PROJ-3  PROJ-4  PROJ-5  PROJ-6
     X  PROJ-7')
 1020   FORMAT(' ------  ------  ------  ------  ------  ------  ------
     X  ------')
 1030   FORMAT(I5,2X,7F8.4)
        RETURN
        END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                         C
C  Determine projections of objects on 7 principal components.
C                                                                         C
C  Ignore first (trivial) axis.
C                                                                         C
C                                                                         C
C-------------------------------------------------------------------------C
        SUBROUTINE PROJN(N,EVALS,A,Z,VEC)
        REAL    EVALS(N), A(N,N), Z(N,N), VEC(N)
C
        NUM = MIN0(N,7)
        DO 300 J1 = 1, N
           DO 50 L = 1, N
              VEC(L) = A(J1,L)
   50      CONTINUE
           DO 200 J2 = 1, NUM
              A(J1,J2) = 0.0
              DO 100 J3 = 1, N
                 A(J1,J2) = A(J1,J2) + VEC(J3)*Z(J3,N-J2)
  100         CONTINUE
              IF (EVALS(N-J2).GT.0.0) A(J1,J2) = 
     X                                  A(J1,J2)/SQRT(EVALS(N-J2))
              IF (EVALS(N-J2).LE.0.0) A(J1,J2) = 0.0 
  200      CONTINUE
  300   CONTINUE
C
        RETURN
        END
