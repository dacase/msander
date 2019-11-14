c-----------------------------------------------------------------------
      SUBROUTINE matinv(A,N,D,L,M)
      implicit double precision (a-h,o-z)
C
C     ----- STANDARD IBM MATRIX INVERSION ROUTINE -----
c
c     arguments:
c       a:  square matrix of dimension nxn
c       d:  resultant determinant
c       l:  work vector of length n
c       m:  work vector of length n
C
      DIMENSION A(*),L(*),M(*)
C
C     ----- SEARCH FOR LARGEST ELEMENT -----
C
      D = 1.0d0
      NK = -N
      DO 80 K = 1,N
      NK = NK+N
      L(K) = K
      M(K) = K
      KK = NK+K
      BIGA = A(KK)
      DO 25 J = K,N
         IZ = N*(J-1)
         DO 20 I = K,N
            IJ = IZ+I
            IF( ABS(BIGA)- ABS(A(IJ)) .ge. 0.d0) cycle
   15       BIGA = A(IJ)
            L(K) = I
            M(K) = J
   20    CONTINUE
   25 CONTINUE
C
C     ----- INTERCHANGE ROWS -----
C
      J = L(K)
      IF(J-K .le. 0) go to 35
      KI = K-N
      DO I = 1,N
         KI = KI+N
         HOLD = -A(KI)
         JI = KI-K+J
         A(KI) = A(JI)
         A(JI) = HOLD
      END DO
C
C     ----- INTERCHANGE COLUMNS -----
C
   35 I = M(K)
      IF(I-K .le. 0) go to 45
      JP = N*(I-1)
      DO J = 1,N
         JK = NK+J
         JI = JP+J
         HOLD = -A(JK)
         A(JK) = A(JI)
         A(JI) = HOLD
      END DO
C
C     ----- DIVIDE COLUMN BY MINUS PIVOT -----
C
   45 IF(BIGA .ne. 0.d0) go to 48
      D = 0.0d0
      GO TO 150
   48 DO I = 1,N
         IF(I-K .eq. 0) cycle
         IK = NK+I
         A(IK) = A(IK)/(-BIGA)
      END DO
C
C     ----- REDUCE MATRIX -----
C
      DO I = 1,N
         IK = NK+I
         HOLD = A(IK)
         IJ = I-N
         DO J = 1,N
            IJ = IJ+N
            IF(I-K .eq. 0) cycle
            IF(J-K .eq. 0) cycle
            KJ = IJ-I+K
            A(IJ) = HOLD*A(KJ)+A(IJ)
         END DO
      END DO
C
C     ----- DIVIDE ROW BY PIVOT -----
C
      KJ = K-N
      DO J = 1,N
         KJ = KJ+N
         IF(J-K .eq. 0) cycle
         A(KJ) = A(KJ)/BIGA
      END DO
C
C     ----- PRODUCT OF PIVOTS -----
C
      D = D*BIGA
C
C     ----- REPLACE PIVOT BY RECIPROCAL -----
C
      A(KK) = 1.0d0/BIGA
   80 CONTINUE
C
C     ----- FINAL ROW AND COLUMN INTERCHANGE -----
C
      K = N
  100 K = (K-1)
      IF(K .le. 0) go to 150
      I = L(K)
      IF(I-K .le. 0) go to 120
      JQ = N*(K-1)
      JR = N*(I-1)
      DO J = 1,N
         JK = JQ+J
         HOLD = A(JK)
         JI = JR+J
         A(JK) = -A(JI)
         A(JI) = HOLD
      END DO
  120 J = M(K)
      IF(J-K .le. 0) go to 100
      KI = K-N
      DO I = 1,N
         KI = KI+N
         HOLD = A(KI)
         JI = KI-K+J
         A(KI) = -A(JI)
         A(JI) = HOLD
      END DO
      GOTO 100
  150 RETURN
      END
