C   23/11/92 211231325  MEMBER NAME  LPC8     (EISPACK.)    FVS
      SUBROUTINE CGECO(A,LDA,N,IPVT,RCOND,Z)
      INTEGER LDA,N,IPVT(1)
      COMPLEX A(LDA,1),Z(1)
	REAL RCOND
C
C     CGECO FACTORS A COMPLEX MATRIX BY GAUSSIAN ELIMINATION
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, CGEFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW CGECO BY CGESL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW CGECO BY CGESL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW CGECO BY CGEDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW CGECO BY CGEDI.
C
C     ON ENTRY
C
C        A       COMPLEX(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       COMPLEX(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK CGEFA
C     BLAS CAXPY,CDOTC,CSSCAL,SCASUM
C     FORTRAN ABS,AIMAG,AMAX1,CMPLX,CONJG,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX CDOTC,EK,T,WK,WKM
      REAL ANORM,S,SCASUM,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
C
      COMPLEX ZDUM,ZDUM1,ZDUM2,CSIGN1
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
      CSIGN1(ZDUM1,ZDUM2) = CABS1(ZDUM1)*(ZDUM2/CABS1(ZDUM2))
C
C     COMPUTE 1-NORM OF A
C
      ANORM = 0.0E0
      DO 10 J = 1, N
         ANORM = AMAX1(ANORM,SCASUM(N,A(1,J),1))
   10 CONTINUE
C
C     FACTOR
C
      CALL CGEFA(A,LDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  CTRANS(A)*Y = E .
C     CTRANS(A)  IS THE CONJUGATE TRANSPOSE OF A .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(U)*W = E .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C     SOLVE CTRANS(U)*W = E
C
      EK = (1.0E0,0.0E0)
      DO 20 J = 1, N
         Z(J) = (0.0E0,0.0E0)
   20 CONTINUE
      DO 100 K = 1, N
         IF (CABS1(Z(K)) .NE. 0.0E0) EK = CSIGN1(EK,-Z(K))
         IF (CABS1(EK-Z(K)) .LE. CABS1(A(K,K))) GO TO 30
            S = CABS1(A(K,K))/CABS1(EK-Z(K))
            CALL CSSCAL(N,S,Z,1)
            EK = CMPLX(S,0.0E0)*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = CABS1(WK)
         SM = CABS1(WKM)
         IF (CABS1(A(K,K)) .EQ. 0.0E0) GO TO 40
            WK = WK/CONJG(A(K,K))
            WKM = WKM/CONJG(A(K,K))
         GO TO 50
   40    CONTINUE
            WK = (1.0E0,0.0E0)
            WKM = (1.0E0,0.0E0)
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + CABS1(Z(J)+WKM*CONJG(A(K,J)))
               Z(J) = Z(J) + WK*CONJG(A(K,J))
               S = S + CABS1(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*CONJG(A(K,J))
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
C     SOLVE CTRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + CDOTC(N-K,A(K+1,K),1,Z(K+1),1)
         IF (CABS1(Z(K)) .LE. 1.0E0) GO TO 110
            S = 1.0E0/CABS1(Z(K))
            CALL CSSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
      YNORM = 1.0E0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL CAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (CABS1(Z(K)) .LE. 1.0E0) GO TO 130
            S = 1.0E0/CABS1(Z(K))
            CALL CSSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (CABS1(Z(K)) .LE. CABS1(A(K,K))) GO TO 150
            S = CABS1(A(K,K))/CABS1(Z(K))
            CALL CSSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (CABS1(A(K,K)) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)
         IF (CABS1(A(K,K)) .EQ. 0.0E0) Z(K) = (1.0E0,0.0E0)
         T = -Z(K)
         CALL CAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
      SUBROUTINE CGEFA(A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(1),INFO
      COMPLEX A(LDA,1)
C
C     CGEFA FACTORS A COMPLEX MATRIX BY GAUSSIAN ELIMINATION.
C
C     CGEFA IS USUALLY CALLED BY CGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR CGECO) = (1 + 9/N)*(TIME FOR CGEFA) .
C
C     ON ENTRY
C
C        A       COMPLEX(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT CGESL OR CGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN CGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CSCAL,ICAMAX
C     FORTRAN ABS,AIMAG,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX T
      INTEGER ICAMAX,J,K,KP1,L,NM1
C
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = ICAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (CABS1(A(L,K)) .EQ. 0.0E0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -(1.0E0,0.0E0)/A(K,K)
            CALL CSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL CAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (CABS1(A(N,N)) .EQ. 0.0E0) INFO = N
      RETURN
      END
      SUBROUTINE CGESL(A,LDA,N,IPVT,B,JOB)
      INTEGER LDA,N,IPVT(1),JOB
      COMPLEX A(LDA,1),B(1)
C
C     CGESL SOLVES THE COMPLEX SYSTEM
C     A * X = B  OR  CTRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY CGECO OR CGEFA.
C
C     ON ENTRY
C
C        A       COMPLEX(LDA, N)
C                THE OUTPUT FROM CGECO OR CGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM CGECO OR CGEFA.
C
C        B       COMPLEX(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  CTRANS(A)*X = B  WHERE
C                            CTRANS(A)  IS THE CONJUGATE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF CGECO HAS SET RCOND .GT. 0.0
C        OR CGEFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL CGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL CGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CDOTC
C     FORTRAN CONJG
C
C     INTERNAL VARIABLES
C
      COMPLEX CDOTC,T
      INTEGER K,KB,L,NM1
C
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL CAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL CAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  CTRANS(A) * X = B
C        FIRST SOLVE  CTRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = CDOTC(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/CONJG(A(K,K))
   60    CONTINUE
C
C        NOW SOLVE CTRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + CDOTC(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE CGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
      INTEGER LDA,N,IPVT(1),JOB
      COMPLEX A(LDA,1),DET(2),WORK(1)
C
C     CGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
C     USING THE FACTORS COMPUTED BY CGECO OR CGEFA.
C
C     ON ENTRY
C
C        A       COMPLEX(LDA, N)
C                THE OUTPUT FROM CGECO OR CGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM CGECO OR CGEFA.
C
C        WORK    COMPLEX(N)
C                WORK VECTOR.  CONTENTS DESTROYED.
C
C        JOB     INTEGER
C                = 11   BOTH DETERMINANT AND INVERSE.
C                = 01   INVERSE ONLY.
C                = 10   DETERMINANT ONLY.
C
C     ON RETURN
C
C        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE UNCHANGED.
C
C        DET     COMPLEX(2)
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE NOT REFERENCED.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. CABS1(DET(1)) .LT. 10.0
C                OR  DET(1) .EQ. 0.0 .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
C        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
C        AND IF CGECO HAS SET RCOND .GT. 0.0 OR CGEFA HAS SET
C        INFO .EQ. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CSCAL,CSWAP
C     FORTRAN ABS,AIMAG,CMPLX,MOD,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX T
      REAL TEN
      INTEGER I,J,K,KB,KP1,L,NM1
C
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = (1.0E0,0.0E0)
         DET(2) = (0.0E0,0.0E0)
         TEN = 10.0E0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
C        ...EXIT
            IF (CABS1(DET(1)) .EQ. 0.0E0) GO TO 60
   10       IF (CABS1(DET(1)) .GE. 1.0E0) GO TO 20
               DET(1) = CMPLX(TEN,0.0E0)*DET(1)
               DET(2) = DET(2) - (1.0E0,0.0E0)
            GO TO 10
   20       CONTINUE
   30       IF (CABS1(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/CMPLX(TEN,0.0E0)
               DET(2) = DET(2) + (1.0E0,0.0E0)
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = (1.0E0,0.0E0)/A(K,K)
            T = -A(K,K)
            CALL CSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = (0.0E0,0.0E0)
               CALL CAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = (0.0E0,0.0E0)
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL CAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL CSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
      SUBROUTINE CGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
      INTEGER LDA,N,ML,MU,IPVT(1)
      COMPLEX ABD(LDA,1),Z(1)
      REAL RCOND
C
C     CGBCO FACTORS A COMPLEX BAND MATRIX BY GAUSSIAN
C     ELIMINATION AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, CGBFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW CGBCO BY CGBSL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW CGBCO BY CGBSL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW CGBCO BY CGBDI.
C
C     ON ENTRY
C
C        ABD     COMPLEX(LDA, N)
C                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
C                ML+1 THROUGH 2*ML+MU+1 OF  ABD .
C                SEE THE COMMENTS BELOW FOR DETAILS.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C                LDA MUST BE .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.
C
C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C                0 .LE. ML .LT. N .
C
C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C                0 .LE. MU .LT. N .
C                MORE EFFICIENT IF  ML .LE. MU .
C
C     ON RETURN
C
C        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       COMPLEX(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     BAND STORAGE
C
C           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT
C           WILL SET UP THE INPUT.
C
C                   ML = (BAND WIDTH BELOW THE DIAGONAL)
C                   MU = (BAND WIDTH ABOVE THE DIAGONAL)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-MU)
C                      I2 = MIN0(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
C           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
C           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
C           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
C           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
C           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.
C
C     EXAMPLE..  IF THE ORIGINAL MATRIX IS
C
C           11 12 13  0  0  0
C           21 22 23 24  0  0
C            0 32 33 34 35  0
C            0  0 43 44 45 46
C            0  0  0 54 55 56
C            0  0  0  0 65 66
C
C      THEN  N = 6, ML = 1, MU = 2, LDA .GE. 5  AND ABD SHOULD CONTAIN
C
C            *  *  *  +  +  +  , * = NOT USED
C            *  * 13 24 35 46  , + = USED FOR PIVOTING
C            * 12 23 34 45 56
C           11 22 33 44 55 66
C           21 32 43 54 65  *
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK CGBFA
C     BLAS CAXPY,CDOTC,CSSCAL,SCASUM
C     FORTRAN ABS,AIMAG,AMAX1,CMPLX,CONJG,MAX0,MIN0,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX CDOTC,EK,T,WK,WKM
      REAL ANORM,S,SCASUM,SM,YNORM
      INTEGER IS,INFO,J,JU,K,KB,KP1,L,LA,LM,LZ,M,MM
C
      COMPLEX ZDUM,ZDUM1,ZDUM2,CSIGN1
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
      CSIGN1(ZDUM1,ZDUM2) = CABS1(ZDUM1)*(ZDUM2/CABS1(ZDUM2))
C
C     COMPUTE 1-NORM OF A
C
      ANORM = 0.0E0
      L = ML + 1
      IS = L + MU
      DO 10 J = 1, N
         ANORM = AMAX1(ANORM,SCASUM(L,ABD(IS,J),1))
         IF (IS .GT. ML + 1) IS = IS - 1
         IF (J .LE. MU) L = L + 1
         IF (J .GE. N - ML) L = L - 1
   10 CONTINUE
C
C     FACTOR
C
      CALL CGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  CTRANS(A)*Y = E .
C     CTRANS(A)  IS THE CONJUGATE TRANSPOSE OF A .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(U)*W = E .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C     SOLVE CTRANS(U)*W = E
C
      EK = (1.0E0,0.0E0)
      DO 20 J = 1, N
         Z(J) = (0.0E0,0.0E0)
   20 CONTINUE
      M = ML + MU + 1
      JU = 0
      DO 100 K = 1, N
         IF (CABS1(Z(K)) .NE. 0.0E0) EK = CSIGN1(EK,-Z(K))
         IF (CABS1(EK-Z(K)) .LE. CABS1(ABD(M,K))) GO TO 30
            S = CABS1(ABD(M,K))/CABS1(EK-Z(K))
            CALL CSSCAL(N,S,Z,1)
            EK = CMPLX(S,0.0E0)*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = CABS1(WK)
         SM = CABS1(WKM)
         IF (CABS1(ABD(M,K)) .EQ. 0.0E0) GO TO 40
            WK = WK/CONJG(ABD(M,K))
            WKM = WKM/CONJG(ABD(M,K))
         GO TO 50
   40    CONTINUE
            WK = (1.0E0,0.0E0)
            WKM = (1.0E0,0.0E0)
   50    CONTINUE
         KP1 = K + 1
         JU = MIN0(MAX0(JU,MU+IPVT(K)),N)
         MM = M
         IF (KP1 .GT. JU) GO TO 90
            DO 60 J = KP1, JU
               MM = MM - 1
               SM = SM + CABS1(Z(J)+WKM*CONJG(ABD(MM,J)))
               Z(J) = Z(J) + WK*CONJG(ABD(MM,J))
               S = S + CABS1(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               MM = M
               DO 70 J = KP1, JU
                  MM = MM - 1
                  Z(J) = Z(J) + T*CONJG(ABD(MM,J))
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
C     SOLVE CTRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         LM = MIN0(ML,N-K)
         IF (K .LT. N) Z(K) = Z(K) + CDOTC(LM,ABD(M+1,K),1,Z(K+1),1)
         IF (CABS1(Z(K)) .LE. 1.0E0) GO TO 110
            S = 1.0E0/CABS1(Z(K))
            CALL CSSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
      YNORM = 1.0E0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         LM = MIN0(ML,N-K)
         IF (K .LT. N) CALL CAXPY(LM,T,ABD(M+1,K),1,Z(K+1),1)
         IF (CABS1(Z(K)) .LE. 1.0E0) GO TO 130
            S = 1.0E0/CABS1(Z(K))
            CALL CSSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = W
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (CABS1(Z(K)) .LE. CABS1(ABD(M,K))) GO TO 150
            S = CABS1(ABD(M,K))/CABS1(Z(K))
            CALL CSSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (CABS1(ABD(M,K)) .NE. 0.0E0) Z(K) = Z(K)/ABD(M,K)
         IF (CABS1(ABD(M,K)) .EQ. 0.0E0) Z(K) = (1.0E0,0.0E0)
         LM = MIN0(K,M) - 1
         LA = M - LM
         LZ = K - LM
         T = -Z(K)
         CALL CAXPY(LM,T,ABD(LA,K),1,Z(LZ),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
      SUBROUTINE CGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)
      INTEGER LDA,N,ML,MU,IPVT(1),INFO
      COMPLEX ABD(LDA,1)
C
C     CGBFA FACTORS A COMPLEX BAND MATRIX BY ELIMINATION.
C
C     CGBFA IS USUALLY CALLED BY CGBCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C
C     ON ENTRY
C
C        ABD     COMPLEX(LDA, N)
C                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
C                ML+1 THROUGH 2*ML+MU+1 OF  ABD .
C                SEE THE COMMENTS BELOW FOR DETAILS.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C                LDA MUST BE .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.
C
C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C                0 .LE. ML .LT. N .
C
C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C                0 .LE. MU .LT. N .
C                MORE EFFICIENT IF  ML .LE. MU .
C     ON RETURN
C
C        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT CGBSL WILL DIVIDE BY ZERO IF
C                     CALLED.  USE  RCOND  IN CGBCO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     BAND STORAGE
C
C           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT
C           WILL SET UP THE INPUT.
C
C                   ML = (BAND WIDTH BELOW THE DIAGONAL)
C                   MU = (BAND WIDTH ABOVE THE DIAGONAL)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-MU)
C                      I2 = MIN0(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
C           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
C           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
C           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
C           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
C           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CSCAL,ICAMAX
C     FORTRAN ABS,AIMAG,MAX0,MIN0,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX T
      INTEGER I,ICAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
C
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
      M = ML + MU + 1
      INFO = 0
C
C     ZERO INITIAL FILL-IN COLUMNS
C
      J0 = MU + 2
      J1 = MIN0(N,M) - 1
      IF (J1 .LT. J0) GO TO 30
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = (0.0E0,0.0E0)
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
      JZ = J1
      JU = 0
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 130
      DO 120 K = 1, NM1
         KP1 = K + 1
C
C        ZERO NEXT FILL-IN COLUMN
C
         JZ = JZ + 1
         IF (JZ .GT. N) GO TO 50
         IF (ML .LT. 1) GO TO 50
            DO 40 I = 1, ML
               ABD(I,JZ) = (0.0E0,0.0E0)
   40       CONTINUE
   50    CONTINUE
C
C        FIND L = PIVOT INDEX
C
         LM = MIN0(ML,N-K)
         L = ICAMAX(LM+1,ABD(M,K),1) + M - 1
         IPVT(K) = L + K - M
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (CABS1(ABD(L,K)) .EQ. 0.0E0) GO TO 100
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. M) GO TO 60
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
   60       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -(1.0E0,0.0E0)/ABD(M,K)
            CALL CSCAL(LM,T,ABD(M+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            JU = MIN0(MAX0(JU,MU+IPVT(K)),N)
            MM = M
            IF (JU .LT. KP1) GO TO 90
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .EQ. MM) GO TO 70
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
   70          CONTINUE
               CALL CAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
   80       CONTINUE
   90       CONTINUE
         GO TO 110
  100    CONTINUE
            INFO = K
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE
      IPVT(N) = N
      IF (CABS1(ABD(M,N)) .EQ. 0.0E0) INFO = N
      RETURN
      END
      SUBROUTINE CGBSL(ABD,LDA,N,ML,MU,IPVT,B,JOB)
      INTEGER LDA,N,ML,MU,IPVT(1),JOB
      COMPLEX ABD(LDA,1),B(1)
C
C     CGBSL SOLVES THE COMPLEX BAND SYSTEM
C     A * X = B  OR  CTRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY CGBCO OR CGBFA.
C
C     ON ENTRY
C
C        ABD     COMPLEX(LDA, N)
C                THE OUTPUT FROM CGBCO OR CGBFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C
C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.
C
C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C
C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM CGBCO OR CGBFA.
C
C        B       COMPLEX(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  CTRANS(A)*X = B , WHERE
C                            CTRANS(A)  IS THE CONJUGATE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF CGBCO HAS SET RCOND .GT. 0.0
C        OR CGBFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL CGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL CGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CDOTC
C     FORTRAN CONJG,MIN0
C
C     INTERNAL VARIABLES
C
      COMPLEX CDOTC,T
      INTEGER K,KB,L,LA,LB,LM,M,NM1
C
      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE L*Y = B
C
         IF (ML .EQ. 0) GO TO 30
         IF (NM1 .LT. 1) GO TO 30
            DO 20 K = 1, NM1
               LM = MIN0(ML,N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .EQ. K) GO TO 10
                  B(L) = B(K)
                  B(K) = T
   10          CONTINUE
               CALL CAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
   20       CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/ABD(M,K)
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL CAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  CTRANS(A) * X = B
C        FIRST SOLVE  CTRANS(U)*Y = B
C
         DO 60 K = 1, N
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = CDOTC(LM,ABD(LA,K),1,B(LB),1)
            B(K) = (B(K) - T)/CONJG(ABD(M,K))
   60    CONTINUE
C
C        NOW SOLVE CTRANS(L)*X = Y
C
         IF (ML .EQ. 0) GO TO 90
         IF (NM1 .LT. 1) GO TO 90
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN0(ML,N-K)
               B(K) = B(K) + CDOTC(LM,ABD(M+1,K),1,B(K+1),1)
               L = IPVT(K)
               IF (L .EQ. K) GO TO 70
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE CGBDI(ABD,LDA,N,ML,MU,IPVT,DET)
      INTEGER LDA,N,ML,MU,IPVT(1)
      COMPLEX ABD(LDA,1),DET(2)
C
C     CGBDI COMPUTES THE DETERMINANT OF A BAND MATRIX
C     USING THE FACTORS COMPUTED BY CGBCO OR CGBFA.
C     IF THE INVERSE IS NEEDED, USE CGBSL  N  TIMES.
C
C     ON ENTRY
C
C        ABD     COMPLEX(LDA, N)
C                THE OUTPUT FROM CGBCO OR CGBFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C
C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.
C
C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C
C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM CGBCO OR CGBFA.
C
C     ON RETURN
C
C        DET     COMPLEX(2)
C                DETERMINANT OF ORIGINAL MATRIX.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. CABS1(DET(1)) .LT. 10.0
C                OR  DET(1) = 0.0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     FORTRAN ABS,AIMAG,CMPLX,REAL
C
C     INTERNAL VARIABLES
C
      REAL TEN
      INTEGER I,M
C
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
      M = ML + MU + 1
      DET(1) = (1.0E0,0.0E0)
      DET(2) = (0.0E0,0.0E0)
      TEN = 10.0E0
      DO 50 I = 1, N
         IF (IPVT(I) .NE. I) DET(1) = -DET(1)
         DET(1) = ABD(M,I)*DET(1)
C     ...EXIT
         IF (CABS1(DET(1)) .EQ. 0.0E0) GO TO 60
   10    IF (CABS1(DET(1)) .GE. 1.0E0) GO TO 20
            DET(1) = CMPLX(TEN,0.0E0)*DET(1)
            DET(2) = DET(2) - (1.0E0,0.0E0)
         GO TO 10
   20    CONTINUE
   30    IF (CABS1(DET(1)) .LT. TEN) GO TO 40
            DET(1) = DET(1)/CMPLX(TEN,0.0E0)
            DET(2) = DET(2) + (1.0E0,0.0E0)
         GO TO 30
   40    CONTINUE
   50 CONTINUE
   60 CONTINUE
      RETURN
      END
      SUBROUTINE CPOCO(A,LDA,N,RCOND,Z,INFO)
      INTEGER LDA,N,INFO
      COMPLEX A(LDA,1),Z(1)
      REAL RCOND
C
C     CPOCO FACTORS A COMPLEX HERMITIAN POSITIVE DEFINITE MATRIX
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, CPOFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW CPOCO BY CPOSL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW CPOCO BY CPOSL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW CPOCO BY CPODI.
C     TO COMPUTE  INVERSE(A) , FOLLOW CPOCO BY CPODI.
C
C     ON ENTRY
C
C        A       COMPLEX(LDA, N)
C                THE HERMITIAN MATRIX TO BE FACTORED.  ONLY THE
C                DIAGONAL AND UPPER TRIANGLE ARE USED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A =
C                CTRANS(R)*R WHERE  CTRANS(R)  IS THE CONJUGATE
C                TRANSPOSE.  THE STRICT LOWER TRIANGLE IS UNALTERED.
C                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED.
C
C        Z       COMPLEX(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C                IF  INFO .NE. 0 , Z  IS UNCHANGED.
C
C        INFO    INTEGER
C                = 0  FOR NORMAL RETURN.
C                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
C                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK CPOFA
C     BLAS CAXPY,CDOTC,CSSCAL,SCASUM
C     FORTRAN ABS,AIMAG,AMAX1,CMPLX,CONJG,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX CDOTC,EK,T,WK,WKM
      REAL ANORM,S,SCASUM,SM,YNORM
      INTEGER I,J,JM1,K,KB,KP1
C
      COMPLEX ZDUM,ZDUM2,CSIGN1
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
      CSIGN1(ZDUM,ZDUM2) = CABS1(ZDUM)*(ZDUM2/CABS1(ZDUM2))
C
C     FIND NORM OF A USING ONLY UPPER HALF
C
      DO 30 J = 1, N
         Z(J) = CMPLX(SCASUM(J,A(1,J),1),0.0E0)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = CMPLX(REAL(Z(I))+CABS1(A(I,J)),0.0E0)
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0E0
      DO 40 J = 1, N
         ANORM = AMAX1(ANORM,REAL(Z(J)))
   40 CONTINUE
C
C     FACTOR
C
      CALL CPOFA(A,LDA,N,INFO)
      IF (INFO .NE. 0) GO TO 180
C
C        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C        GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(R)*W = E .
C        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C        SOLVE CTRANS(R)*W = E
C
         EK = (1.0E0,0.0E0)
         DO 50 J = 1, N
            Z(J) = (0.0E0,0.0E0)
   50    CONTINUE
         DO 110 K = 1, N
            IF (CABS1(Z(K)) .NE. 0.0E0) EK = CSIGN1(EK,-Z(K))
            IF (CABS1(EK-Z(K)) .LE. REAL(A(K,K))) GO TO 60
               S = REAL(A(K,K))/CABS1(EK-Z(K))
               CALL CSSCAL(N,S,Z,1)
               EK = CMPLX(S,0.0E0)*EK
   60       CONTINUE
            WK = EK - Z(K)
            WKM = -EK - Z(K)
            S = CABS1(WK)
            SM = CABS1(WKM)
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
            KP1 = K + 1
            IF (KP1 .GT. N) GO TO 100
               DO 70 J = KP1, N
                  SM = SM + CABS1(Z(J)+WKM*CONJG(A(K,J)))
                  Z(J) = Z(J) + WK*CONJG(A(K,J))
                  S = S + CABS1(Z(J))
   70          CONTINUE
               IF (S .GE. SM) GO TO 90
                  T = WKM - WK
                  WK = WKM
                  DO 80 J = KP1, N
                     Z(J) = Z(J) + T*CONJG(A(K,J))
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
            Z(K) = WK
  110    CONTINUE
         S = 1.0E0/SCASUM(N,Z,1)
         CALL CSSCAL(N,S,Z,1)
C
C        SOLVE R*Y = W
C
         DO 130 KB = 1, N
            K = N + 1 - KB
            IF (CABS1(Z(K)) .LE. REAL(A(K,K))) GO TO 120
               S = REAL(A(K,K))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
  120       CONTINUE
            Z(K) = Z(K)/A(K,K)
            T = -Z(K)
            CALL CAXPY(K-1,T,A(1,K),1,Z(1),1)
  130    CONTINUE
         S = 1.0E0/SCASUM(N,Z,1)
         CALL CSSCAL(N,S,Z,1)
C
         YNORM = 1.0E0
C
C        SOLVE CTRANS(R)*V = Y
C
         DO 150 K = 1, N
            Z(K) = Z(K) - CDOTC(K-1,A(1,K),1,Z(1),1)
            IF (CABS1(Z(K)) .LE. REAL(A(K,K))) GO TO 140
               S = REAL(A(K,K))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
               YNORM = S*YNORM
  140       CONTINUE
            Z(K) = Z(K)/A(K,K)
  150    CONTINUE
         S = 1.0E0/SCASUM(N,Z,1)
         CALL CSSCAL(N,S,Z,1)
         YNORM = S*YNORM
C
C        SOLVE R*Z = V
C
         DO 170 KB = 1, N
            K = N + 1 - KB
            IF (CABS1(Z(K)) .LE. REAL(A(K,K))) GO TO 160
               S = REAL(A(K,K))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
               YNORM = S*YNORM
  160       CONTINUE
            Z(K) = Z(K)/A(K,K)
            T = -Z(K)
            CALL CAXPY(K-1,T,A(1,K),1,Z(1),1)
  170    CONTINUE
C        MAKE ZNORM = 1.0
         S = 1.0E0/SCASUM(N,Z,1)
         CALL CSSCAL(N,S,Z,1)
         YNORM = S*YNORM
C
         IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
         IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
  180 CONTINUE
      RETURN
      END
      SUBROUTINE CPOFA(A,LDA,N,INFO)
      INTEGER LDA,N,INFO
      COMPLEX A(LDA,1)
C
C     CPOFA FACTORS A COMPLEX HERMITIAN POSITIVE DEFINITE MATRIX.
C
C     CPOFA IS USUALLY CALLED BY CPOCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR CPOCO) = (1 + 18/N)*(TIME FOR CPOFA) .
C
C     ON ENTRY
C
C        A       COMPLEX(LDA, N)
C                THE HERMITIAN MATRIX TO BE FACTORED.  ONLY THE
C                DIAGONAL AND UPPER TRIANGLE ARE USED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A =
C                CTRANS(R)*R WHERE  CTRANS(R)  IS THE CONJUGATE
C                TRANSPOSE.  THE STRICT LOWER TRIANGLE IS UNALTERED.
C                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
C
C        INFO    INTEGER
C                = 0  FOR NORMAL RETURN.
C                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
C                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CDOTC
C     FORTRAN AIMAG,CMPLX,CONJG,REAL,SQRT
C
C     INTERNAL VARIABLES
C
      COMPLEX CDOTC,T
      REAL S
      INTEGER J,JM1,K
C     BEGIN BLOCK WITH ...EXITS TO 40
C
C
         DO 30 J = 1, N
            INFO = J
            S = 0.0E0
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               T = A(K,J) - CDOTC(K-1,A(1,K),1,A(1,J),1)
               T = T/A(K,K)
               A(K,J) = T
               S = S + REAL(T*CONJG(T))
   10       CONTINUE
   20       CONTINUE
            S = REAL(A(J,J)) - S
C     ......EXIT
            IF (S .LE. 0.0E0 .OR. AIMAG(A(J,J)) .NE. 0.0E0) GO TO 40
            A(J,J) = CMPLX(SQRT(S),0.0E0)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
      SUBROUTINE CPOSL(A,LDA,N,B)
      INTEGER LDA,N
      COMPLEX A(LDA,1),B(1)
C
C     CPOSL SOLVES THE COMPLEX HERMITIAN POSITIVE DEFINITE SYSTEM
C     A * X = B
C     USING THE FACTORS COMPUTED BY CPOCO OR CPOFA.
C
C     ON ENTRY
C
C        A       COMPLEX(LDA, N)
C                THE OUTPUT FROM CPOCO OR CPOFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        B       COMPLEX(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES
C        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE
C        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED
C        CORRECTLY AND  INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL CPOCO(A,LDA,N,RCOND,Z,INFO)
C           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL CPOSL(A,LDA,N,C(1,J))
C        10 CONTINUE
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CDOTC
C
C     INTERNAL VARIABLES
C
      COMPLEX CDOTC,T
      INTEGER K,KB
C
C     SOLVE CTRANS(R)*Y = B
C
      DO 10 K = 1, N
         T = CDOTC(K-1,A(1,K),1,B(1),1)
         B(K) = (B(K) - T)/A(K,K)
   10 CONTINUE
C
C     SOLVE R*X = Y
C
      DO 20 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/A(K,K)
         T = -B(K)
         CALL CAXPY(K-1,T,A(1,K),1,B(1),1)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE CPODI(A,LDA,N,DET,JOB)
      INTEGER LDA,N,JOB
      COMPLEX A(LDA,1)
      REAL DET(2)
C
C     CPODI COMPUTES THE DETERMINANT AND INVERSE OF A CERTAIN
C     COMPLEX HERMITIAN POSITIVE DEFINITE MATRIX (SEE BELOW)
C     USING THE FACTORS COMPUTED BY CPOCO, CPOFA OR CQRDC.
C
C     ON ENTRY
C
C        A       COMPLEX(LDA, N)
C                THE OUTPUT  A  FROM CPOCO OR CPOFA
C                OR THE OUTPUT  X  FROM CQRDC.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        JOB     INTEGER
C                = 11   BOTH DETERMINANT AND INVERSE.
C                = 01   INVERSE ONLY.
C                = 10   DETERMINANT ONLY.
C
C     ON RETURN
C
C        A       IF CPOCO OR CPOFA WAS USED TO FACTOR  A  THEN
C                CPODI PRODUCES THE UPPER HALF OF INVERSE(A) .
C                IF CQRDC WAS USED TO DECOMPOSE  X  THEN
C                CPODI PRODUCES THE UPPER HALF OF INVERSE(CTRANS(X)*X)
C                WHERE CTRANS(X) IS THE CONJUGATE TRANSPOSE.
C                ELEMENTS OF  A  BELOW THE DIAGONAL ARE UNCHANGED.
C                IF THE UNITS DIGIT OF JOB IS ZERO,  A  IS UNCHANGED.
C
C        DET     REAL(2)
C                DETERMINANT OF  A  OR OF  CTRANS(X)*X  IF REQUESTED.
C                OTHERWISE NOT REFERENCED.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. DET(1) .LT. 10.0
C                OR  DET(1) .EQ. 0.0 .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
C        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
C        AND IF CPOCO OR CPOFA HAS SET INFO .EQ. 0 .
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CSCAL
C     FORTRAN CONJG,MOD,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX T
      REAL S
      INTEGER I,J,JM1,K,KP1
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0E0
         DET(2) = 0.0E0
         S = 10.0E0
         DO 50 I = 1, N
            DET(1) = REAL(A(I,I))**2*DET(1)
C        ...EXIT
            IF (DET(1) .EQ. 0.0E0) GO TO 60
   10       IF (DET(1) .GE. 1.0E0) GO TO 20
               DET(1) = S*DET(1)
               DET(2) = DET(2) - 1.0E0
            GO TO 10
   20       CONTINUE
   30       IF (DET(1) .LT. S) GO TO 40
               DET(1) = DET(1)/S
               DET(2) = DET(2) + 1.0E0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(R)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 140
         DO 100 K = 1, N
            A(K,K) = (1.0E0,0.0E0)/A(K,K)
            T = -A(K,K)
            CALL CSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = (0.0E0,0.0E0)
               CALL CAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM  INVERSE(R) * CTRANS(INVERSE(R))
C
         DO 130 J = 1, N
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 120
            DO 110 K = 1, JM1
               T = CONJG(A(K,J))
               CALL CAXPY(K,T,A(1,J),1,A(1,K),1)
  110       CONTINUE
  120       CONTINUE
            T = CONJG(A(J,J))
            CALL CSCAL(J,T,A(1,J),1)
  130    CONTINUE
  140 CONTINUE
      RETURN
      END
      SUBROUTINE CPPCO(AP,N,RCOND,Z,INFO)
      INTEGER N,INFO
      COMPLEX AP(1),Z(1)
      REAL RCOND
C
C     CPPCO FACTORS A COMPLEX HERMITIAN POSITIVE DEFINITE MATRIX
C     STORED IN PACKED FORM
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, CPPFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW CPPCO BY CPPSL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW CPPCO BY CPPSL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW CPPCO BY CPPDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW CPPCO BY CPPDI.
C
C     ON ENTRY
C
C        AP      COMPLEX (N*(N+1)/2)
C                THE PACKED FORM OF A HERMITIAN MATRIX  A .  THE
C                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY
C                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
C                SEE COMMENTS BELOW FOR DETAILS.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        AP      AN UPPER TRIANGULAR MATRIX  R , STORED IN PACKED
C                FORM, SO THAT  A = CTRANS(R)*R .
C                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED.
C
C        Z       COMPLEX(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS SINGULAR TO WORKING PRECISION, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C                IF  INFO .NE. 0 , Z  IS UNCHANGED.
C
C        INFO    INTEGER
C                = 0  FOR NORMAL RETURN.
C                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
C                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
C
C     PACKED STORAGE
C
C          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
C          TRIANGLE OF A HERMITIAN MATRIX.
C
C                K = 0
C                DO 20 J = 1, N
C                   DO 10 I = 1, J
C                      K = K + 1
C                      AP(K) = A(I,J)
C             10    CONTINUE
C             20 CONTINUE
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK CPPFA
C     BLAS CAXPY,CDOTC,CSSCAL,SCASUM
C     FORTRAN ABS,AIMAG,AMAX1,CMPLX,CONJG,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX CDOTC,EK,T,WK,WKM
      REAL ANORM,S,SCASUM,SM,YNORM
      INTEGER I,IJ,J,JM1,J1,K,KB,KJ,KK,KP1
C
      COMPLEX ZDUM,ZDUM2,CSIGN1
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
      CSIGN1(ZDUM,ZDUM2) = CABS1(ZDUM)*(ZDUM2/CABS1(ZDUM2))
C
C     FIND NORM OF A
C
      J1 = 1
      DO 30 J = 1, N
         Z(J) = CMPLX(SCASUM(J,AP(J1),1),0.0E0)
         IJ = J1
         J1 = J1 + J
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = CMPLX(REAL(Z(I))+CABS1(AP(IJ)),0.0E0)
            IJ = IJ + 1
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0E0
      DO 40 J = 1, N
         ANORM = AMAX1(ANORM,REAL(Z(J)))
   40 CONTINUE
C
C     FACTOR
C
      CALL CPPFA(AP,N,INFO)
      IF (INFO .NE. 0) GO TO 180
C
C        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C        GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(R)*W = E .
C        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C        SOLVE CTRANS(R)*W = E
C
         EK = (1.0E0,0.0E0)
         DO 50 J = 1, N
            Z(J) = (0.0E0,0.0E0)
   50    CONTINUE
         KK = 0
         DO 110 K = 1, N
            KK = KK + K
            IF (CABS1(Z(K)) .NE. 0.0E0) EK = CSIGN1(EK,-Z(K))
            IF (CABS1(EK-Z(K)) .LE. REAL(AP(KK))) GO TO 60
               S = REAL(AP(KK))/CABS1(EK-Z(K))
               CALL CSSCAL(N,S,Z,1)
               EK = CMPLX(S,0.0E0)*EK
   60       CONTINUE
            WK = EK - Z(K)
            WKM = -EK - Z(K)
            S = CABS1(WK)
            SM = CABS1(WKM)
            WK = WK/AP(KK)
            WKM = WKM/AP(KK)
            KP1 = K + 1
            KJ = KK + K
            IF (KP1 .GT. N) GO TO 100
               DO 70 J = KP1, N
                  SM = SM + CABS1(Z(J)+WKM*CONJG(AP(KJ)))
                  Z(J) = Z(J) + WK*CONJG(AP(KJ))
                  S = S + CABS1(Z(J))
                  KJ = KJ + J
   70          CONTINUE
               IF (S .GE. SM) GO TO 90
                  T = WKM - WK
                  WK = WKM
                  KJ = KK + K
                  DO 80 J = KP1, N
                     Z(J) = Z(J) + T*CONJG(AP(KJ))
                     KJ = KJ + J
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
            Z(K) = WK
  110    CONTINUE
         S = 1.0E0/SCASUM(N,Z,1)
         CALL CSSCAL(N,S,Z,1)
C
C        SOLVE R*Y = W
C
         DO 130 KB = 1, N
            K = N + 1 - KB
            IF (CABS1(Z(K)) .LE. REAL(AP(KK))) GO TO 120
               S = REAL(AP(KK))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
  120       CONTINUE
            Z(K) = Z(K)/AP(KK)
            KK = KK - K
            T = -Z(K)
            CALL CAXPY(K-1,T,AP(KK+1),1,Z(1),1)
  130    CONTINUE
         S = 1.0E0/SCASUM(N,Z,1)
         CALL CSSCAL(N,S,Z,1)
C
         YNORM = 1.0E0
C
C        SOLVE CTRANS(R)*V = Y
C
         DO 150 K = 1, N
            Z(K) = Z(K) - CDOTC(K-1,AP(KK+1),1,Z(1),1)
            KK = KK + K
            IF (CABS1(Z(K)) .LE. REAL(AP(KK))) GO TO 140
               S = REAL(AP(KK))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
               YNORM = S*YNORM
  140       CONTINUE
            Z(K) = Z(K)/AP(KK)
  150    CONTINUE
         S = 1.0E0/SCASUM(N,Z,1)
         CALL CSSCAL(N,S,Z,1)
         YNORM = S*YNORM
C
C        SOLVE R*Z = V
C
         DO 170 KB = 1, N
            K = N + 1 - KB
            IF (CABS1(Z(K)) .LE. REAL(AP(KK))) GO TO 160
               S = REAL(AP(KK))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
               YNORM = S*YNORM
  160       CONTINUE
            Z(K) = Z(K)/AP(KK)
            KK = KK - K
            T = -Z(K)
            CALL CAXPY(K-1,T,AP(KK+1),1,Z(1),1)
  170    CONTINUE
C        MAKE ZNORM = 1.0
         S = 1.0E0/SCASUM(N,Z,1)
         CALL CSSCAL(N,S,Z,1)
         YNORM = S*YNORM
C
         IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
         IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
  180 CONTINUE
      RETURN
      END
      SUBROUTINE CPPFA(AP,N,INFO)
      INTEGER N,INFO
      COMPLEX AP(1)
C
C     CPPFA FACTORS A COMPLEX HERMITIAN POSITIVE DEFINITE MATRIX
C     STORED IN PACKED FORM.
C
C     CPPFA IS USUALLY CALLED BY CPPCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR CPPCO) = (1 + 18/N)*(TIME FOR CPPFA) .
C
C     ON ENTRY
C
C        AP      COMPLEX (N*(N+1)/2)
C                THE PACKED FORM OF A HERMITIAN MATRIX  A .  THE
C                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY
C                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
C                SEE COMMENTS BELOW FOR DETAILS.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        AP      AN UPPER TRIANGULAR MATRIX  R , STORED IN PACKED
C                FORM, SO THAT  A = CTRANS(R)*R .
C
C        INFO    INTEGER
C                = 0  FOR NORMAL RETURN.
C                = K  IF THE LEADING MINOR OF ORDER  K  IS NOT
C                     POSITIVE DEFINITE.
C
C
C     PACKED STORAGE
C
C          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
C          TRIANGLE OF A HERMITIAN MATRIX.
C
C                K = 0
C                DO 20 J = 1, N
C                   DO 10 I = 1, J
C                      K = K + 1
C                      AP(K) = A(I,J)
C             10    CONTINUE
C             20 CONTINUE
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CDOTC
C     FORTRAN AIMAG,CMPLX,CONJG,REAL,SQRT
C
C     INTERNAL VARIABLES
C
      COMPLEX CDOTC,T
      REAL S
      INTEGER J,JJ,JM1,K,KJ,KK
C     BEGIN BLOCK WITH ...EXITS TO 40
C
C
         JJ = 0
         DO 30 J = 1, N
            INFO = J
            S = 0.0E0
            JM1 = J - 1
            KJ = JJ
            KK = 0
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               KJ = KJ + 1
               T = AP(KJ) - CDOTC(K-1,AP(KK+1),1,AP(JJ+1),1)
               KK = KK + K
               T = T/AP(KK)
               AP(KJ) = T
               S = S + REAL(T*CONJG(T))
   10       CONTINUE
   20       CONTINUE
            JJ = JJ + J
            S = REAL(AP(JJ)) - S
C     ......EXIT
            IF (S .LE. 0.0E0 .OR. AIMAG(AP(JJ)) .NE. 0.0E0) GO TO 40
            AP(JJ) = CMPLX(SQRT(S),0.0E0)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
      SUBROUTINE CPPSL(AP,N,B)
      INTEGER N
      COMPLEX AP(1),B(1)
C
C     CPPSL SOLVES THE COMPLEX HERMITIAN POSITIVE DEFINITE SYSTEM
C     A * X = B
C     USING THE FACTORS COMPUTED BY CPPCO OR CPPFA.
C
C     ON ENTRY
C
C        AP      COMPLEX (N*(N+1)/2)
C                THE OUTPUT FROM CPPCO OR CPPFA.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        B       COMPLEX(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES
C        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE
C        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED
C        CORRECTLY AND  INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL CPPCO(AP,N,RCOND,Z,INFO)
C           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL CPPSL(AP,N,C(1,J))
C        10 CONTINUE
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CDOTC
C
C     INTERNAL VARIABLES
C
      COMPLEX CDOTC,T
      INTEGER K,KB,KK
C
      KK = 0
      DO 10 K = 1, N
         T = CDOTC(K-1,AP(KK+1),1,B(1),1)
         KK = KK + K
         B(K) = (B(K) - T)/AP(KK)
   10 CONTINUE
      DO 20 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/AP(KK)
         KK = KK - K
         T = -B(K)
         CALL CAXPY(K-1,T,AP(KK+1),1,B(1),1)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE CPPDI(AP,N,DET,JOB)
      INTEGER N,JOB
      COMPLEX AP(1)
      REAL DET(2)
C
C     CPPDI COMPUTES THE DETERMINANT AND INVERSE
C     OF A COMPLEX HERMITIAN POSITIVE DEFINITE MATRIX
C     USING THE FACTORS COMPUTED BY CPPCO OR CPPFA .
C
C     ON ENTRY
C
C        AP      COMPLEX (N*(N+1)/2)
C                THE OUTPUT FROM CPPCO OR CPPFA.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        JOB     INTEGER
C                = 11   BOTH DETERMINANT AND INVERSE.
C                = 01   INVERSE ONLY.
C                = 10   DETERMINANT ONLY.
C
C     ON RETURN
C
C        AP      THE UPPER TRIANGULAR HALF OF THE INVERSE .
C                THE STRICT LOWER TRIANGLE IS UNALTERED.
C
C        DET     REAL(2)
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE NOT REFERENCED.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. DET(1) .LT. 10.0
C                OR  DET(1) .EQ. 0.0 .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
C        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
C        AND IF CPOCO OR CPOFA HAS SET INFO .EQ. 0 .
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CSCAL
C     FORTRAN CONJG,MOD,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX T
      REAL S
      INTEGER I,II,J,JJ,JM1,J1,K,KJ,KK,KP1,K1
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0E0
         DET(2) = 0.0E0
         S = 10.0E0
         II = 0
         DO 50 I = 1, N
            II = II + I
            DET(1) = REAL(AP(II))**2*DET(1)
C        ...EXIT
            IF (DET(1) .EQ. 0.0E0) GO TO 60
   10       IF (DET(1) .GE. 1.0E0) GO TO 20
               DET(1) = S*DET(1)
               DET(2) = DET(2) - 1.0E0
            GO TO 10
   20       CONTINUE
   30       IF (DET(1) .LT. S) GO TO 40
               DET(1) = DET(1)/S
               DET(2) = DET(2) + 1.0E0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(R)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 140
         KK = 0
         DO 100 K = 1, N
            K1 = KK + 1
            KK = KK + K
            AP(KK) = (1.0E0,0.0E0)/AP(KK)
            T = -AP(KK)
            CALL CSCAL(K-1,T,AP(K1),1)
            KP1 = K + 1
            J1 = KK + 1
            KJ = KK + K
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = AP(KJ)
               AP(KJ) = (0.0E0,0.0E0)
               CALL CAXPY(K,T,AP(K1),1,AP(J1),1)
               J1 = J1 + J
               KJ = KJ + J
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM  INVERSE(R) * CTRANS(INVERSE(R))
C
         JJ = 0
         DO 130 J = 1, N
            J1 = JJ + 1
            JJ = JJ + J
            JM1 = J - 1
            K1 = 1
            KJ = J1
            IF (JM1 .LT. 1) GO TO 120
            DO 110 K = 1, JM1
               T = CONJG(AP(KJ))
               CALL CAXPY(K,T,AP(J1),1,AP(K1),1)
               K1 = K1 + K
               KJ = KJ + 1
  110       CONTINUE
  120       CONTINUE
            T = CONJG(AP(JJ))
            CALL CSCAL(J,T,AP(J1),1)
  130    CONTINUE
  140 CONTINUE
      RETURN
      END
      SUBROUTINE CPBCO(ABD,LDA,N,M,RCOND,Z,INFO)
      INTEGER LDA,N,M,INFO
      COMPLEX ABD(LDA,1),Z(1)
      REAL RCOND
C
C     CPBCO FACTORS A COMPLEX HERMITIAN POSITIVE DEFINITE MATRIX
C     STORED IN BAND FORM AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, CPBFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW CPBCO BY CPBSL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW CPBCO BY CPBSL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW CPBCO BY CPBDI.
C
C     ON ENTRY
C
C        ABD     COMPLEX(LDA, N)
C                THE MATRIX TO BE FACTORED.  THE COLUMNS OF THE UPPER
C                TRIANGLE ARE STORED IN THE COLUMNS OF ABD AND THE
C                DIAGONALS OF THE UPPER TRIANGLE ARE STORED IN THE
C                ROWS OF ABD .  SEE THE COMMENTS BELOW FOR DETAILS.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C                LDA MUST BE .GE. M + 1 .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        M       INTEGER
C                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C                0 .LE. M .LT. N .
C
C     ON RETURN
C
C        ABD     AN UPPER TRIANGULAR MATRIX  R , STORED IN BAND
C                FORM, SO THAT  A = CTRANS(R)*R .
C                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED.
C
C        Z       COMPLEX(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS SINGULAR TO WORKING PRECISION, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C                IF  INFO .NE. 0 , Z  IS UNCHANGED.
C
C        INFO    INTEGER
C                = 0  FOR NORMAL RETURN.
C                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
C                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
C
C     BAND STORAGE
C
C           IF  A  IS A HERMITIAN POSITIVE DEFINITE BAND MATRIX,
C           THE FOLLOWING PROGRAM SEGMENT WILL SET UP THE INPUT.
C
C                   M = (BAND WIDTH ABOVE DIAGONAL)
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-M)
C                      DO 10 I = I1, J
C                         K = I-J+M+1
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           THIS USES  M + 1  ROWS OF  A , EXCEPT FOR THE  M BY M
C           UPPER LEFT TRIANGLE, WHICH IS IGNORED.
C
C     EXAMPLE..  IF THE ORIGINAL MATRIX IS
C
C           11 12 13  0  0  0
C           12 22 23 24  0  0
C           13 23 33 34 35  0
C            0 24 34 44 45 46
C            0  0 35 45 55 56
C            0  0  0 46 56 66
C
C     THEN  N = 6 , M = 2  AND  ABD  SHOULD CONTAIN
C
C            *  * 13 24 35 46
C            * 12 23 34 45 56
C           11 22 33 44 55 66
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK CPBFA
C     BLAS CAXPY,CDOTC,CSSCAL,SCASUM
C     FORTRAN ABS,AIMAG,AMAX1,CMPLX,CONJG,MAX0,MIN0,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX CDOTC,EK,T,WK,WKM
      REAL ANORM,S,SCASUM,SM,YNORM
      INTEGER I,J,J2,K,KB,KP1,L,LA,LB,LM,MU
C
      COMPLEX ZDUM,ZDUM2,CSIGN1
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
      CSIGN1(ZDUM,ZDUM2) = CABS1(ZDUM)*(ZDUM2/CABS1(ZDUM2))
C
C     FIND NORM OF A
C
      DO 30 J = 1, N
         L = MIN0(J,M+1)
         MU = MAX0(M+2-J,1)
         Z(J) = CMPLX(SCASUM(L,ABD(MU,J),1),0.0E0)
         K = J - L
         IF (M .LT. MU) GO TO 20
         DO 10 I = MU, M
            K = K + 1
            Z(K) = CMPLX(REAL(Z(K))+CABS1(ABD(I,J)),0.0E0)
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0E0
      DO 40 J = 1, N
         ANORM = AMAX1(ANORM,REAL(Z(J)))
   40 CONTINUE
C
C     FACTOR
C
      CALL CPBFA(ABD,LDA,N,M,INFO)
      IF (INFO .NE. 0) GO TO 180
C
C        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C        GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(R)*W = E .
C        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C        SOLVE CTRANS(R)*W = E
C
         EK = (1.0E0,0.0E0)
         DO 50 J = 1, N
            Z(J) = (0.0E0,0.0E0)
   50    CONTINUE
         DO 110 K = 1, N
            IF (CABS1(Z(K)) .NE. 0.0E0) EK = CSIGN1(EK,-Z(K))
            IF (CABS1(EK-Z(K)) .LE. REAL(ABD(M+1,K))) GO TO 60
               S = REAL(ABD(M+1,K))/CABS1(EK-Z(K))
               CALL CSSCAL(N,S,Z,1)
               EK = CMPLX(S,0.0E0)*EK
   60       CONTINUE
            WK = EK - Z(K)
            WKM = -EK - Z(K)
            S = CABS1(WK)
            SM = CABS1(WKM)
            WK = WK/ABD(M+1,K)
            WKM = WKM/ABD(M+1,K)
            KP1 = K + 1
            J2 = MIN0(K+M,N)
            I = M + 1
            IF (KP1 .GT. J2) GO TO 100
               DO 70 J = KP1, J2
                  I = I - 1
                  SM = SM + CABS1(Z(J)+WKM*CONJG(ABD(I,J)))
                  Z(J) = Z(J) + WK*CONJG(ABD(I,J))
                  S = S + CABS1(Z(J))
   70          CONTINUE
               IF (S .GE. SM) GO TO 90
                  T = WKM - WK
                  WK = WKM
                  I = M + 1
                  DO 80 J = KP1, J2
                     I = I - 1
                     Z(J) = Z(J) + T*CONJG(ABD(I,J))
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
            Z(K) = WK
  110    CONTINUE
         S = 1.0E0/SCASUM(N,Z,1)
         CALL CSSCAL(N,S,Z,1)
C
C        SOLVE  R*Y = W
C
         DO 130 KB = 1, N
            K = N + 1 - KB
            IF (CABS1(Z(K)) .LE. REAL(ABD(M+1,K))) GO TO 120
               S = REAL(ABD(M+1,K))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
  120       CONTINUE
            Z(K) = Z(K)/ABD(M+1,K)
            LM = MIN0(K-1,M)
            LA = M + 1 - LM
            LB = K - LM
            T = -Z(K)
            CALL CAXPY(LM,T,ABD(LA,K),1,Z(LB),1)
  130    CONTINUE
         S = 1.0E0/SCASUM(N,Z,1)
         CALL CSSCAL(N,S,Z,1)
C
         YNORM = 1.0E0
C
C        SOLVE CTRANS(R)*V = Y
C
         DO 150 K = 1, N
            LM = MIN0(K-1,M)
            LA = M + 1 - LM
            LB = K - LM
            Z(K) = Z(K) - CDOTC(LM,ABD(LA,K),1,Z(LB),1)
            IF (CABS1(Z(K)) .LE. REAL(ABD(M+1,K))) GO TO 140
               S = REAL(ABD(M+1,K))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
               YNORM = S*YNORM
  140       CONTINUE
            Z(K) = Z(K)/ABD(M+1,K)
  150    CONTINUE
         S = 1.0E0/SCASUM(N,Z,1)
         CALL CSSCAL(N,S,Z,1)
         YNORM = S*YNORM
C
C        SOLVE  R*Z = W
C
         DO 170 KB = 1, N
            K = N + 1 - KB
            IF (CABS1(Z(K)) .LE. REAL(ABD(M+1,K))) GO TO 160
               S = REAL(ABD(M+1,K))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
               YNORM = S*YNORM
  160       CONTINUE
            Z(K) = Z(K)/ABD(M+1,K)
            LM = MIN0(K-1,M)
            LA = M + 1 - LM
            LB = K - LM
            T = -Z(K)
            CALL CAXPY(LM,T,ABD(LA,K),1,Z(LB),1)
  170    CONTINUE
C        MAKE ZNORM = 1.0
         S = 1.0E0/SCASUM(N,Z,1)
         CALL CSSCAL(N,S,Z,1)
         YNORM = S*YNORM
C
         IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
         IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
  180 CONTINUE
      RETURN
      END
      SUBROUTINE CPBFA(ABD,LDA,N,M,INFO)
      INTEGER LDA,N,M,INFO
      COMPLEX ABD(LDA,1)
C
C     CPBFA FACTORS A COMPLEX HERMITIAN POSITIVE DEFINITE MATRIX
C     STORED IN BAND FORM.
C
C     CPBFA IS USUALLY CALLED BY CPBCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C
C     ON ENTRY
C
C        ABD     COMPLEX(LDA, N)
C                THE MATRIX TO BE FACTORED.  THE COLUMNS OF THE UPPER
C                TRIANGLE ARE STORED IN THE COLUMNS OF ABD AND THE
C                DIAGONALS OF THE UPPER TRIANGLE ARE STORED IN THE
C                ROWS OF ABD .  SEE THE COMMENTS BELOW FOR DETAILS.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C                LDA MUST BE .GE. M + 1 .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        M       INTEGER
C                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C                0 .LE. M .LT. N .
C
C     ON RETURN
C
C        ABD     AN UPPER TRIANGULAR MATRIX  R , STORED IN BAND
C                FORM, SO THAT  A = CTRANS(R)*R .
C
C        INFO    INTEGER
C                = 0  FOR NORMAL RETURN.
C                = K  IF THE LEADING MINOR OF ORDER  K  IS NOT
C                     POSITIVE DEFINITE.
C
C     BAND STORAGE
C
C           IF  A  IS A HERMITIAN POSITIVE DEFINITE BAND MATRIX,
C           THE FOLLOWING PROGRAM SEGMENT WILL SET UP THE INPUT.
C
C                   M = (BAND WIDTH ABOVE DIAGONAL)
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-M)
C                      DO 10 I = I1, J
C                         K = I-J+M+1
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CDOTC
C     FORTRAN AIMAG,CMPLX,CONJG,MAX0,REAL,SQRT
C
C     INTERNAL VARIABLES
C
      COMPLEX CDOTC,T
      REAL S
      INTEGER IK,J,JK,K,MU
C     BEGIN BLOCK WITH ...EXITS TO 40
C
C
         DO 30 J = 1, N
            INFO = J
            S = 0.0E0
            IK = M + 1
            JK = MAX0(J-M,1)
            MU = MAX0(M+2-J,1)
            IF (M .LT. MU) GO TO 20
            DO 10 K = MU, M
               T = ABD(K,J) - CDOTC(K-MU,ABD(IK,JK),1,ABD(MU,J),1)
               T = T/ABD(M+1,JK)
               ABD(K,J) = T
               S = S + REAL(T*CONJG(T))
               IK = IK - 1
               JK = JK + 1
   10       CONTINUE
   20       CONTINUE
            S = REAL(ABD(M+1,J)) - S
C     ......EXIT
            IF (S .LE. 0.0E0 .OR. AIMAG(ABD(M+1,J)) .NE. 0.0E0)
     *         GO TO 40
            ABD(M+1,J) = CMPLX(SQRT(S),0.0E0)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
      SUBROUTINE CPBSL(ABD,LDA,N,M,B)
      INTEGER LDA,N,M
      COMPLEX ABD(LDA,1),B(1)
C
C     CPBSL SOLVES THE COMPLEX HERMITIAN POSITIVE DEFINITE BAND
C     SYSTEM  A*X = B
C     USING THE FACTORS COMPUTED BY CPBCO OR CPBFA.
C
C     ON ENTRY
C
C        ABD     COMPLEX(LDA, N)
C                THE OUTPUT FROM CPBCO OR CPBFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        M       INTEGER
C                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C
C        B       COMPLEX(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES
C        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE
C        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED
C        CORRECTLY AND  INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL CPBCO(ABD,LDA,N,RCOND,Z,INFO)
C           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL CPBSL(ABD,LDA,N,C(1,J))
C        10 CONTINUE
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CDOTC
C     FORTRAN MIN0
C
C     INTERNAL VARIABLES
C
      COMPLEX CDOTC,T
      INTEGER K,KB,LA,LB,LM
C
C     SOLVE CTRANS(R)*Y = B
C
      DO 10 K = 1, N
         LM = MIN0(K-1,M)
         LA = M + 1 - LM
         LB = K - LM
         T = CDOTC(LM,ABD(LA,K),1,B(LB),1)
         B(K) = (B(K) - T)/ABD(M+1,K)
   10 CONTINUE
C
C     SOLVE R*X = Y
C
      DO 20 KB = 1, N
         K = N + 1 - KB
         LM = MIN0(K-1,M)
         LA = M + 1 - LM
         LB = K - LM
         B(K) = B(K)/ABD(M+1,K)
         T = -B(K)
         CALL CAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE CPBDI(ABD,LDA,N,M,DET)
      INTEGER LDA,N,M
      COMPLEX ABD(LDA,1)
      REAL DET(2)
C
C     CPBDI COMPUTES THE DETERMINANT
C     OF A COMPLEX HERMITIAN POSITIVE DEFINITE BAND MATRIX
C     USING THE FACTORS COMPUTED BY CPBCO OR CPBFA.
C     IF THE INVERSE IS NEEDED, USE CPBSL  N  TIMES.
C
C     ON ENTRY
C
C        ABD     COMPLEX(LDA, N)
C                THE OUTPUT FROM CPBCO OR CPBFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        M       INTEGER
C                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C
C     ON RETURN
C
C        DET     REAL(2)
C                DETERMINANT OF ORIGINAL MATRIX IN THE FORM
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. DET(1) .LT. 10.0
C                OR  DET(1) .EQ. 0.0 .
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     FORTRAN REAL
C
C     INTERNAL VARIABLES
C
      REAL S
      INTEGER I
C
C     COMPUTE DETERMINANT
C
      DET(1) = 1.0E0
      DET(2) = 0.0E0
      S = 10.0E0
      DO 50 I = 1, N
         DET(1) = REAL(ABD(M+1,I))**2*DET(1)
C     ...EXIT
         IF (DET(1) .EQ. 0.0E0) GO TO 60
   10    IF (DET(1) .GE. 1.0E0) GO TO 20
            DET(1) = S*DET(1)
            DET(2) = DET(2) - 1.0E0
         GO TO 10
   20    CONTINUE
   30    IF (DET(1) .LT. S) GO TO 40
            DET(1) = DET(1)/S
            DET(2) = DET(2) + 1.0E0
         GO TO 30
   40    CONTINUE
   50 CONTINUE
   60 CONTINUE
      RETURN
      END
      SUBROUTINE CSICO(A,LDA,N,KPVT,RCOND,Z)
      INTEGER LDA,N,KPVT(1)
      COMPLEX A(LDA,1),Z(1)
      REAL RCOND
C
C     CSICO FACTORS A COMPLEX SYMMETRIC MATRIX BY ELIMINATION WITH
C     SYMMETRIC PIVOTING AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, CSIFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW CSICO BY CSISL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW CSICO BY CSISL.
C     TO COMPUTE  INVERSE(A) , FOLLOW CSICO BY CSIDI.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW CSICO BY CSIDI.
C
C     ON ENTRY
C
C        A       COMPLEX(LDA, N)
C                THE SYMMETRIC MATRIX TO BE FACTORED.
C                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
C                WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
C                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
C                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
C
C        KPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       COMPLEX(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK CSIFA
C     BLAS CAXPY,CDOTU,CSSCAL,SCASUM
C     FORTRAN ABS,AIMAG,AMAX1,CMPLX,IABS,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX AK,AKM1,BK,BKM1,CDOTU,DENOM,EK,T
      REAL ANORM,S,SCASUM,YNORM
      INTEGER I,INFO,J,JM1,K,KP,KPS,KS
C
      COMPLEX ZDUM,ZDUM2,CSIGN1
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
      CSIGN1(ZDUM,ZDUM2) = CABS1(ZDUM)*(ZDUM2/CABS1(ZDUM2))
C
C     FIND NORM OF A USING ONLY UPPER HALF
C
      DO 30 J = 1, N
         Z(J) = CMPLX(SCASUM(J,A(1,J),1),0.0E0)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = CMPLX(REAL(Z(I))+CABS1(A(I,J)),0.0E0)
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0E0
      DO 40 J = 1, N
         ANORM = AMAX1(ANORM,REAL(Z(J)))
   40 CONTINUE
C
C     FACTOR
C
      CALL CSIFA(A,LDA,N,KPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C     SOLVE U*D*W = E
C
      EK = (1.0E0,0.0E0)
      DO 50 J = 1, N
         Z(J) = (0.0E0,0.0E0)
   50 CONTINUE
      K = N
   60 IF (K .EQ. 0) GO TO 120
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         KP = IABS(KPVT(K))
         KPS = K + 1 - KS
         IF (KP .EQ. KPS) GO TO 70
            T = Z(KPS)
            Z(KPS) = Z(KP)
            Z(KP) = T
   70    CONTINUE
         IF (CABS1(Z(K)) .NE. 0.0E0) EK = CSIGN1(EK,Z(K))
         Z(K) = Z(K) + EK
         CALL CAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)
         IF (KS .EQ. 1) GO TO 80
            IF (CABS1(Z(K-1)) .NE. 0.0E0) EK = CSIGN1(EK,Z(K-1))
            Z(K-1) = Z(K-1) + EK
            CALL CAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)
   80    CONTINUE
         IF (KS .EQ. 2) GO TO 100
            IF (CABS1(Z(K)) .LE. CABS1(A(K,K))) GO TO 90
               S = CABS1(A(K,K))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
               EK = CMPLX(S,0.0E0)*EK
   90       CONTINUE
            IF (CABS1(A(K,K)) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)
            IF (CABS1(A(K,K)) .EQ. 0.0E0) Z(K) = (1.0E0,0.0E0)
         GO TO 110
  100    CONTINUE
            AK = A(K,K)/A(K-1,K)
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = Z(K)/A(K-1,K)
            BKM1 = Z(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0E0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  110    CONTINUE
         K = K - KS
      GO TO 60
  120 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
C     SOLVE TRANS(U)*Y = W
C
      K = 1
  130 IF (K .GT. N) GO TO 160
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 150
            Z(K) = Z(K) + CDOTU(K-1,A(1,K),1,Z(1),1)
            IF (KS .EQ. 2)
     *         Z(K+1) = Z(K+1) + CDOTU(K-1,A(1,K+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 140
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  140       CONTINUE
  150    CONTINUE
         K = K + KS
      GO TO 130
  160 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
      YNORM = 1.0E0
C
C     SOLVE U*D*V = Y
C
      K = N
  170 IF (K .EQ. 0) GO TO 230
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. KS) GO TO 190
            KP = IABS(KPVT(K))
            KPS = K + 1 - KS
            IF (KP .EQ. KPS) GO TO 180
               T = Z(KPS)
               Z(KPS) = Z(KP)
               Z(KP) = T
  180       CONTINUE
            CALL CAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)
            IF (KS .EQ. 2) CALL CAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)
  190    CONTINUE
         IF (KS .EQ. 2) GO TO 210
            IF (CABS1(Z(K)) .LE. CABS1(A(K,K))) GO TO 200
               S = CABS1(A(K,K))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
               YNORM = S*YNORM
  200       CONTINUE
            IF (CABS1(A(K,K)) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)
            IF (CABS1(A(K,K)) .EQ. 0.0E0) Z(K) = (1.0E0,0.0E0)
         GO TO 220
  210    CONTINUE
            AK = A(K,K)/A(K-1,K)
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = Z(K)/A(K-1,K)
            BKM1 = Z(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0E0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  220    CONTINUE
         K = K - KS
      GO TO 170
  230 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE TRANS(U)*Z = V
C
      K = 1
  240 IF (K .GT. N) GO TO 270
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 260
            Z(K) = Z(K) + CDOTU(K-1,A(1,K),1,Z(1),1)
            IF (KS .EQ. 2)
     *         Z(K+1) = Z(K+1) + CDOTU(K-1,A(1,K+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 250
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  250       CONTINUE
  260    CONTINUE
         K = K + KS
      GO TO 240
  270 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
      SUBROUTINE CSIFA(A,LDA,N,KPVT,INFO)
      INTEGER LDA,N,KPVT(1),INFO
      COMPLEX A(LDA,1)
C
C     CSIFA FACTORS A COMPLEX SYMMETRIC MATRIX BY ELIMINATION
C     WITH SYMMETRIC PIVOTING.
C
C     TO SOLVE  A*X = B , FOLLOW CSIFA BY CSISL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW CSIFA BY CSISL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW CSIFA BY CSIDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW CSIFA BY CSIDI.
C
C     ON ENTRY
C
C        A       COMPLEX(LDA,N)
C                THE SYMMETRIC MATRIX TO BE FACTORED.
C                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
C                WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
C                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
C                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
C
C        KPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS
C                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,
C                     BUT IT DOES INDICATE THAT CSISL OR CSIDI MAY
C                     DIVIDE BY ZERO IF CALLED.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CSWAP,ICAMAX
C     FORTRAN ABS,AIMAG,AMAX1,REAL,SQRT
C
C     INTERNAL VARIABLES
C
      COMPLEX AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
      REAL ABSAKK,ALPHA,COLMAX,ROWMAX
      INTEGER IMAX,IMAXP1,J,JJ,JMAX,K,KM1,KM2,KSTEP,ICAMAX
      LOGICAL SWAP
C
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
C     INITIALIZE
C
C     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
      ALPHA = (1.0E0 + SQRT(17.0E0))/8.0E0
C
      INFO = 0
C
C     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
C
      K = N
   10 CONTINUE
C
C        LEAVE THE LOOP IF K=0 OR K=1.
C
C     ...EXIT
         IF (K .EQ. 0) GO TO 200
         IF (K .GT. 1) GO TO 20
            KPVT(1) = 1
            IF (CABS1(A(1,1)) .EQ. 0.0E0) INFO = 1
C     ......EXIT
            GO TO 200
   20    CONTINUE
C
C        THIS SECTION OF CODE DETERMINES THE KIND OF
C        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
C        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
C        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
C        REQUIRED.
C
         KM1 = K - 1
         ABSAKK = CABS1(A(K,K))
C
C        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C        COLUMN K.
C
         IMAX = ICAMAX(K-1,A(1,K),1)
         COLMAX = CABS1(A(IMAX,K))
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30
            KSTEP = 1
            SWAP = .FALSE.
         GO TO 90
   30    CONTINUE
C
C           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C           ROW IMAX.
C
            ROWMAX = 0.0E0
            IMAXP1 = IMAX + 1
            DO 40 J = IMAXP1, K
               ROWMAX = AMAX1(ROWMAX,CABS1(A(IMAX,J)))
   40       CONTINUE
            IF (IMAX .EQ. 1) GO TO 50
               JMAX = ICAMAX(IMAX-1,A(1,IMAX),1)
               ROWMAX = AMAX1(ROWMAX,CABS1(A(JMAX,IMAX)))
   50       CONTINUE
            IF (CABS1(A(IMAX,IMAX)) .LT. ALPHA*ROWMAX) GO TO 60
               KSTEP = 1
               SWAP = .TRUE.
            GO TO 80
   60       CONTINUE
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70
               KSTEP = 1
               SWAP = .FALSE.
            GO TO 80
   70       CONTINUE
               KSTEP = 2
               SWAP = IMAX .NE. KM1
   80       CONTINUE
   90    CONTINUE
         IF (AMAX1(ABSAKK,COLMAX) .NE. 0.0E0) GO TO 100
C
C           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
C
            KPVT(K) = K
            INFO = K
         GO TO 190
  100    CONTINUE
         IF (KSTEP .EQ. 2) GO TO 140
C
C           1 X 1 PIVOT BLOCK.
C
            IF (.NOT.SWAP) GO TO 120
C
C              PERFORM AN INTERCHANGE.
C
               CALL CSWAP(IMAX,A(1,IMAX),1,A(1,K),1)
               DO 110 JJ = IMAX, K
                  J = K + IMAX - JJ
                  T = A(J,K)
                  A(J,K) = A(IMAX,J)
                  A(IMAX,J) = T
  110          CONTINUE
  120       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            DO 130 JJ = 1, KM1
               J = K - JJ
               MULK = -A(J,K)/A(K,K)
               T = MULK
               CALL CAXPY(J,T,A(1,K),1,A(1,J),1)
               A(J,K) = MULK
  130       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = K
            IF (SWAP) KPVT(K) = IMAX
         GO TO 190
  140    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (.NOT.SWAP) GO TO 160
C
C              PERFORM AN INTERCHANGE.
C
               CALL CSWAP(IMAX,A(1,IMAX),1,A(1,K-1),1)
               DO 150 JJ = IMAX, KM1
                  J = KM1 + IMAX - JJ
                  T = A(J,K-1)
                  A(J,K-1) = A(IMAX,J)
                  A(IMAX,J) = T
  150          CONTINUE
               T = A(K-1,K)
               A(K-1,K) = A(IMAX,K)
               A(IMAX,K) = T
  160       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            KM2 = K - 2
            IF (KM2 .EQ. 0) GO TO 180
               AK = A(K,K)/A(K-1,K)
               AKM1 = A(K-1,K-1)/A(K-1,K)
               DENOM = 1.0E0 - AK*AKM1
               DO 170 JJ = 1, KM2
                  J = KM1 - JJ
                  BK = A(J,K)/A(K-1,K)
                  BKM1 = A(J,K-1)/A(K-1,K)
                  MULK = (AKM1*BK - BKM1)/DENOM
                  MULKM1 = (AK*BKM1 - BK)/DENOM
                  T = MULK
                  CALL CAXPY(J,T,A(1,K),1,A(1,J),1)
                  T = MULKM1
                  CALL CAXPY(J,T,A(1,K-1),1,A(1,J),1)
                  A(J,K) = MULK
                  A(J,K-1) = MULKM1
  170          CONTINUE
  180       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = 1 - K
            IF (SWAP) KPVT(K) = -IMAX
            KPVT(K-1) = KPVT(K)
  190    CONTINUE
         K = K - KSTEP
      GO TO 10
  200 CONTINUE
      RETURN
      END
      SUBROUTINE CSISL(A,LDA,N,KPVT,B)
      INTEGER LDA,N,KPVT(1)
      COMPLEX A(LDA,1),B(1)
C
C     CSISL SOLVES THE COMPLEX SYMMETRIC SYSTEM
C     A * X = B
C     USING THE FACTORS COMPUTED BY CSIFA.
C
C     ON ENTRY
C
C        A       COMPLEX(LDA,N)
C                THE OUTPUT FROM CSIFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        KPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM CSIFA.
C
C        B       COMPLEX(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO MAY OCCUR IF  CSICO  HAS SET RCOND .EQ. 0.0
C        OR  CSIFA  HAS SET INFO .NE. 0  .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL CSIFA(A,LDA,N,KPVT,INFO)
C           IF (INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL CSISL(A,LDA,N,KPVT,C(1,J))
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CDOTU
C     FORTRAN IABS
C
C     INTERNAL VARIABLES.
C
      COMPLEX AK,AKM1,BK,BKM1,CDOTU,DENOM,TEMP
      INTEGER K,KP
C
C     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
C     D INVERSE TO B.
C
      K = N
   10 IF (K .EQ. 0) GO TO 80
         IF (KPVT(K) .LT. 0) GO TO 40
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 30
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 20
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
   20          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL CAXPY(K-1,B(K),A(1,K),1,B(1),1)
   30       CONTINUE
C
C           APPLY D INVERSE.
C
            B(K) = B(K)/A(K,K)
            K = K - 1
         GO TO 70
   40    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (K .EQ. 2) GO TO 60
               KP = IABS(KPVT(K))
               IF (KP .EQ. K - 1) GO TO 50
C
C                 INTERCHANGE.
C
                  TEMP = B(K-1)
                  B(K-1) = B(KP)
                  B(KP) = TEMP
   50          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL CAXPY(K-2,B(K),A(1,K),1,B(1),1)
               CALL CAXPY(K-2,B(K-1),A(1,K-1),1,B(1),1)
   60       CONTINUE
C
C           APPLY D INVERSE.
C
            AK = A(K,K)/A(K-1,K)
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = B(K)/A(K-1,K)
            BKM1 = B(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0E0
            B(K) = (AKM1*BK - BKM1)/DENOM
            B(K-1) = (AK*BKM1 - BK)/DENOM
            K = K - 2
   70    CONTINUE
      GO TO 10
   80 CONTINUE
C
C     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
C
      K = 1
   90 IF (K .GT. N) GO TO 160
         IF (KPVT(K) .LT. 0) GO TO 120
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 110
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + CDOTU(K-1,A(1,K),1,B(1),1)
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 100
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  100          CONTINUE
  110       CONTINUE
            K = K + 1
         GO TO 150
  120    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 140
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + CDOTU(K-1,A(1,K),1,B(1),1)
               B(K+1) = B(K+1) + CDOTU(K-1,A(1,K+1),1,B(1),1)
               KP = IABS(KPVT(K))
               IF (KP .EQ. K) GO TO 130
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  130          CONTINUE
  140       CONTINUE
            K = K + 2
  150    CONTINUE
      GO TO 90
  160 CONTINUE
      RETURN
      END
      SUBROUTINE CSIDI(A,LDA,N,KPVT,DET,WORK,JOB)
      INTEGER LDA,N,JOB
      COMPLEX A(LDA,1),DET(2),WORK(1)
      INTEGER KPVT(1)
C
C     CSIDI COMPUTES THE DETERMINANT AND INVERSE
C     OF A COMPLEX SYMMETRIC MATRIX USING THE FACTORS FROM CSIFA.
C
C     ON ENTRY
C
C        A       COMPLEX(LDA,N)
C                THE OUTPUT FROM CSIFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY A.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX A.
C
C        KPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM CSIFA.
C
C        WORK    COMPLEX(N)
C                WORK VECTOR.  CONTENTS DESTROYED.
C
C        JOB     INTEGER
C                JOB HAS THE DECIMAL EXPANSION  AB  WHERE
C                   IF  B .NE. 0, THE INVERSE IS COMPUTED,
C                   IF  A .NE. 0, THE DETERMINANT IS COMPUTED,
C
C                FOR EXAMPLE, JOB = 11  GIVES BOTH.
C
C     ON RETURN
C
C        VARIABLES NOT REQUESTED BY JOB ARE NOT USED.
C
C        A      CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF
C               THE ORIGINAL MATRIX.  THE STRICT LOWER TRIANGLE
C               IS NEVER REFERENCED.
C
C        DET    COMPLEX(2)
C               DETERMINANT OF ORIGINAL MATRIX.
C               DETERMINANT = DET(1) * 10.0**DET(2)
C               WITH 1.0 .LE. ABS(DET(1)) .LT. 10.0
C               OR DET(1) = 0.0.
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO MAY OCCUR IF THE INVERSE IS REQUESTED
C        AND  CSICO  HAS SET RCOND .EQ. 0.0
C        OR  CSIFA  HAS SET  INFO .NE. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CCOPY,CDOTU,CSWAP
C     FORTRAN ABS,CMPLX,IABS,MOD,REAL
C
C     INTERNAL VARIABLES.
C
      COMPLEX AK,AKP1,AKKP1,CDOTU,D,T,TEMP
      REAL TEN
      INTEGER J,JB,K,KM1,KS,KSTEP
      LOGICAL NOINV,NODET
C
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
      NOINV = MOD(JOB,10) .EQ. 0
      NODET = MOD(JOB,100)/10 .EQ. 0
C
      IF (NODET) GO TO 100
         DET(1) = (1.0E0,0.0E0)
         DET(2) = (0.0E0,0.0E0)
         TEN = 10.0E0
         T = (0.0E0,0.0E0)
         DO 90 K = 1, N
            D = A(K,K)
C
C           CHECK IF 1 BY 1
C
            IF (KPVT(K) .GT. 0) GO TO 30
C
C              2 BY 2 BLOCK
C              USE DET (D  T)  =  (D/T * C - T) * T
C                      (T  C)
C              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
C              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
C
               IF (CABS1(T) .NE. 0.0E0) GO TO 10
                  T = A(K,K+1)
                  D = (D/T)*A(K+1,K+1) - T
               GO TO 20
   10          CONTINUE
                  D = T
                  T = (0.0E0,0.0E0)
   20          CONTINUE
   30       CONTINUE
C
            DET(1) = D*DET(1)
            IF (CABS1(DET(1)) .EQ. 0.0E0) GO TO 80
   40          IF (CABS1(DET(1)) .GE. 1.0E0) GO TO 50
                  DET(1) = CMPLX(TEN,0.0E0)*DET(1)
                  DET(2) = DET(2) - (1.0E0,0.0E0)
               GO TO 40
   50          CONTINUE
   60          IF (CABS1(DET(1)) .LT. TEN) GO TO 70
                  DET(1) = DET(1)/CMPLX(TEN,0.0E0)
                  DET(2) = DET(2) + (1.0E0,0.0E0)
               GO TO 60
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
C
C     COMPUTE INVERSE(A)
C
      IF (NOINV) GO TO 230
         K = 1
  110    IF (K .GT. N) GO TO 220
            KM1 = K - 1
            IF (KPVT(K) .LT. 0) GO TO 140
C
C              1 BY 1
C
               A(K,K) = (1.0E0,0.0E0)/A(K,K)
               IF (KM1 .LT. 1) GO TO 130
                  CALL CCOPY(KM1,A(1,K),1,WORK,1)
                  DO 120 J = 1, KM1
                     A(J,K) = CDOTU(J,A(1,J),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  120             CONTINUE
                  A(K,K) = A(K,K) + CDOTU(KM1,WORK,1,A(1,K),1)
  130          CONTINUE
               KSTEP = 1
            GO TO 180
  140       CONTINUE
C
C              2 BY 2
C
               T = A(K,K+1)
               AK = A(K,K)/T
               AKP1 = A(K+1,K+1)/T
               AKKP1 = A(K,K+1)/T
               D = T*(AK*AKP1 - (1.0E0,0.0E0))
               A(K,K) = AKP1/D
               A(K+1,K+1) = AK/D
               A(K,K+1) = -AKKP1/D
               IF (KM1 .LT. 1) GO TO 170
                  CALL CCOPY(KM1,A(1,K+1),1,WORK,1)
                  DO 150 J = 1, KM1
                     A(J,K+1) = CDOTU(J,A(1,J),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),A(1,J),1,A(1,K+1),1)
  150             CONTINUE
                  A(K+1,K+1) = A(K+1,K+1)
     *                         + CDOTU(KM1,WORK,1,A(1,K+1),1)
                  A(K,K+1) = A(K,K+1) + CDOTU(KM1,A(1,K),1,A(1,K+1),1)
                  CALL CCOPY(KM1,A(1,K),1,WORK,1)
                  DO 160 J = 1, KM1
                     A(J,K) = CDOTU(J,A(1,J),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  160             CONTINUE
                  A(K,K) = A(K,K) + CDOTU(KM1,WORK,1,A(1,K),1)
  170          CONTINUE
               KSTEP = 2
  180       CONTINUE
C
C           SWAP
C
            KS = IABS(KPVT(K))
            IF (KS .EQ. K) GO TO 210
               CALL CSWAP(KS,A(1,KS),1,A(1,K),1)
               DO 190 JB = KS, K
                  J = K + KS - JB
                  TEMP = A(J,K)
                  A(J,K) = A(KS,J)
                  A(KS,J) = TEMP
  190          CONTINUE
               IF (KSTEP .EQ. 1) GO TO 200
                  TEMP = A(KS,K+1)
                  A(KS,K+1) = A(K,K+1)
                  A(K,K+1) = TEMP
  200          CONTINUE
  210       CONTINUE
            K = K + KSTEP
         GO TO 110
  220    CONTINUE
  230 CONTINUE
      RETURN
      END
      SUBROUTINE CSPCO(AP,N,KPVT,RCOND,Z)
      INTEGER N,KPVT(1)
      COMPLEX AP(1),Z(1)
      REAL RCOND
C
C     CSPCO FACTORS A COMPLEX SYMMETRIC MATRIX STORED IN PACKED
C     FORM BY ELIMINATION WITH SYMMETRIC PIVOTING AND ESTIMATES
C     THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, CSPFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW CSPCO BY CSPSL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW CSPCO BY CSPSL.
C     TO COMPUTE  INVERSE(A) , FOLLOW CSPCO BY CSPDI.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW CSPCO BY CSPDI.
C
C     ON ENTRY
C
C        AP      COMPLEX (N*(N+1)/2)
C                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE
C                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY
C                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
C                SEE COMMENTS BELOW FOR DETAILS.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
C                WERE USED TO OBTAIN IT STORED IN PACKED FORM.
C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
C                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
C                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
C
C        KPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       COMPLEX(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     PACKED STORAGE
C
C          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
C          TRIANGLE OF A SYMMETRIC MATRIX.
C
C                K = 0
C                DO 20 J = 1, N
C                   DO 10 I = 1, J
C                      K = K + 1
C                      AP(K) = A(I,J)
C             10    CONTINUE
C             20 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK CSPFA
C     BLAS CAXPY,CDOTU,CSSCAL,SCASUM
C     FORTRAN ABS,AIMAG,AMAX1,CMPLX,IABS,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX AK,AKM1,BK,BKM1,CDOTU,DENOM,EK,T
      REAL ANORM,S,SCASUM,YNORM
      INTEGER I,IJ,IK,IKM1,IKP1,INFO,J,JM1,J1
      INTEGER K,KK,KM1K,KM1KM1,KP,KPS,KS
C
      COMPLEX ZDUM,ZDUM2,CSIGN1
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
      CSIGN1(ZDUM,ZDUM2) = CABS1(ZDUM)*(ZDUM2/CABS1(ZDUM2))
C
C     FIND NORM OF A USING ONLY UPPER HALF
C
      J1 = 1
      DO 30 J = 1, N
         Z(J) = CMPLX(SCASUM(J,AP(J1),1),0.0E0)
         IJ = J1
         J1 = J1 + J
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = CMPLX(REAL(Z(I))+CABS1(AP(IJ)),0.0E0)
            IJ = IJ + 1
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0E0
      DO 40 J = 1, N
         ANORM = AMAX1(ANORM,REAL(Z(J)))
   40 CONTINUE
C
C     FACTOR
C
      CALL CSPFA(AP,N,KPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C     SOLVE U*D*W = E
C
      EK = (1.0E0,0.0E0)
      DO 50 J = 1, N
         Z(J) = (0.0E0,0.0E0)
   50 CONTINUE
      K = N
      IK = (N*(N - 1))/2
   60 IF (K .EQ. 0) GO TO 120
         KK = IK + K
         IKM1 = IK - (K - 1)
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         KP = IABS(KPVT(K))
         KPS = K + 1 - KS
         IF (KP .EQ. KPS) GO TO 70
            T = Z(KPS)
            Z(KPS) = Z(KP)
            Z(KP) = T
   70    CONTINUE
         IF (CABS1(Z(K)) .NE. 0.0E0) EK = CSIGN1(EK,Z(K))
         Z(K) = Z(K) + EK
         CALL CAXPY(K-KS,Z(K),AP(IK+1),1,Z(1),1)
         IF (KS .EQ. 1) GO TO 80
            IF (CABS1(Z(K-1)) .NE. 0.0E0) EK = CSIGN1(EK,Z(K-1))
            Z(K-1) = Z(K-1) + EK
            CALL CAXPY(K-KS,Z(K-1),AP(IKM1+1),1,Z(1),1)
   80    CONTINUE
         IF (KS .EQ. 2) GO TO 100
            IF (CABS1(Z(K)) .LE. CABS1(AP(KK))) GO TO 90
               S = CABS1(AP(KK))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
               EK = CMPLX(S,0.0E0)*EK
   90       CONTINUE
            IF (CABS1(AP(KK)) .NE. 0.0E0) Z(K) = Z(K)/AP(KK)
            IF (CABS1(AP(KK)) .EQ. 0.0E0) Z(K) = (1.0E0,0.0E0)
         GO TO 110
  100    CONTINUE
            KM1K = IK + K - 1
            KM1KM1 = IKM1 + K - 1
            AK = AP(KK)/AP(KM1K)
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = Z(K)/AP(KM1K)
            BKM1 = Z(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0E0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  110    CONTINUE
         K = K - KS
         IK = IK - K
         IF (KS .EQ. 2) IK = IK - (K + 1)
      GO TO 60
  120 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
C     SOLVE TRANS(U)*Y = W
C
      K = 1
      IK = 0
  130 IF (K .GT. N) GO TO 160
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 150
            Z(K) = Z(K) + CDOTU(K-1,AP(IK+1),1,Z(1),1)
            IKP1 = IK + K
            IF (KS .EQ. 2)
     *         Z(K+1) = Z(K+1) + CDOTU(K-1,AP(IKP1+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 140
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  140       CONTINUE
  150    CONTINUE
         IK = IK + K
         IF (KS .EQ. 2) IK = IK + (K + 1)
         K = K + KS
      GO TO 130
  160 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
      YNORM = 1.0E0
C
C     SOLVE U*D*V = Y
C
      K = N
      IK = N*(N - 1)/2
  170 IF (K .EQ. 0) GO TO 230
         KK = IK + K
         IKM1 = IK - (K - 1)
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. KS) GO TO 190
            KP = IABS(KPVT(K))
            KPS = K + 1 - KS
            IF (KP .EQ. KPS) GO TO 180
               T = Z(KPS)
               Z(KPS) = Z(KP)
               Z(KP) = T
  180       CONTINUE
            CALL CAXPY(K-KS,Z(K),AP(IK+1),1,Z(1),1)
            IF (KS .EQ. 2) CALL CAXPY(K-KS,Z(K-1),AP(IKM1+1),1,Z(1),1)
  190    CONTINUE
         IF (KS .EQ. 2) GO TO 210
            IF (CABS1(Z(K)) .LE. CABS1(AP(KK))) GO TO 200
               S = CABS1(AP(KK))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
               YNORM = S*YNORM
  200       CONTINUE
            IF (CABS1(AP(KK)) .NE. 0.0E0) Z(K) = Z(K)/AP(KK)
            IF (CABS1(AP(KK)) .EQ. 0.0E0) Z(K) = (1.0E0,0.0E0)
         GO TO 220
  210    CONTINUE
            KM1K = IK + K - 1
            KM1KM1 = IKM1 + K - 1
            AK = AP(KK)/AP(KM1K)
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = Z(K)/AP(KM1K)
            BKM1 = Z(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0E0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  220    CONTINUE
         K = K - KS
         IK = IK - K
         IF (KS .EQ. 2) IK = IK - (K + 1)
      GO TO 170
  230 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE TRANS(U)*Z = V
C
      K = 1
      IK = 0
  240 IF (K .GT. N) GO TO 270
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 260
            Z(K) = Z(K) + CDOTU(K-1,AP(IK+1),1,Z(1),1)
            IKP1 = IK + K
            IF (KS .EQ. 2)
     *         Z(K+1) = Z(K+1) + CDOTU(K-1,AP(IKP1+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 250
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  250       CONTINUE
  260    CONTINUE
         IK = IK + K
         IF (KS .EQ. 2) IK = IK + (K + 1)
         K = K + KS
      GO TO 240
  270 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
      SUBROUTINE CSPFA(AP,N,KPVT,INFO)
      INTEGER N,KPVT(1),INFO
      COMPLEX AP(1)
C
C     CSPFA FACTORS A COMPLEX SYMMETRIC MATRIX STORED IN
C     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING.
C
C     TO SOLVE  A*X = B , FOLLOW CSPFA BY CSPSL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW CSPFA BY CSPSL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW CSPFA BY CSPDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW CSPFA BY CSPDI.
C
C     ON ENTRY
C
C        AP      COMPLEX (N*(N+1)/2)
C                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE
C                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY
C                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
C                SEE COMMENTS BELOW FOR DETAILS.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
C                WERE USED TO OBTAIN IT STORED IN PACKED FORM.
C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
C                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
C                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
C
C        KPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS
C                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,
C                     BUT IT DOES INDICATE THAT CSPSL OR CSPDI MAY
C                     DIVIDE BY ZERO IF CALLED.
C
C     PACKED STORAGE
C
C          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
C          TRIANGLE OF A SYMMETRIC MATRIX.
C
C                K = 0
C                DO 20 J = 1, N
C                   DO 10 I = 1, J
C                      K = K + 1
C                      AP(K)  = A(I,J)
C             10    CONTINUE
C             20 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CSWAP,ICAMAX
C     FORTRAN ABS,AIMAG,AMAX1,REAL,SQRT
C
C     INTERNAL VARIABLES
C
      COMPLEX AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
      REAL ABSAKK,ALPHA,COLMAX,ROWMAX
      INTEGER ICAMAX,IJ,IJJ,IK,IKM1,IM,IMAX,IMAXP1,IMIM,IMJ,IMK
      INTEGER J,JJ,JK,JKM1,JMAX,JMIM,K,KK,KM1,KM1K,KM1KM1,KM2,KSTEP
      LOGICAL SWAP
C
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
C     INITIALIZE
C
C     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
      ALPHA = (1.0E0 + SQRT(17.0E0))/8.0E0
C
      INFO = 0
C
C     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
C
      K = N
      IK = (N*(N - 1))/2
   10 CONTINUE
C
C        LEAVE THE LOOP IF K=0 OR K=1.
C
C     ...EXIT
         IF (K .EQ. 0) GO TO 200
         IF (K .GT. 1) GO TO 20
            KPVT(1) = 1
            IF (CABS1(AP(1)) .EQ. 0.0E0) INFO = 1
C     ......EXIT
            GO TO 200
   20    CONTINUE
C
C        THIS SECTION OF CODE DETERMINES THE KIND OF
C        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
C        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
C        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
C        REQUIRED.
C
         KM1 = K - 1
         KK = IK + K
         ABSAKK = CABS1(AP(KK))
C
C        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C        COLUMN K.
C
         IMAX = ICAMAX(K-1,AP(IK+1),1)
         IMK = IK + IMAX
         COLMAX = CABS1(AP(IMK))
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30
            KSTEP = 1
            SWAP = .FALSE.
         GO TO 90
   30    CONTINUE
C
C           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C           ROW IMAX.
C
            ROWMAX = 0.0E0
            IMAXP1 = IMAX + 1
            IM = IMAX*(IMAX - 1)/2
            IMJ = IM + 2*IMAX
            DO 40 J = IMAXP1, K
               ROWMAX = AMAX1(ROWMAX,CABS1(AP(IMJ)))
               IMJ = IMJ + J
   40       CONTINUE
            IF (IMAX .EQ. 1) GO TO 50
               JMAX = ICAMAX(IMAX-1,AP(IM+1),1)
               JMIM = JMAX + IM
               ROWMAX = AMAX1(ROWMAX,CABS1(AP(JMIM)))
   50       CONTINUE
            IMIM = IMAX + IM
            IF (CABS1(AP(IMIM)) .LT. ALPHA*ROWMAX) GO TO 60
               KSTEP = 1
               SWAP = .TRUE.
            GO TO 80
   60       CONTINUE
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70
               KSTEP = 1
               SWAP = .FALSE.
            GO TO 80
   70       CONTINUE
               KSTEP = 2
               SWAP = IMAX .NE. KM1
   80       CONTINUE
   90    CONTINUE
         IF (AMAX1(ABSAKK,COLMAX) .NE. 0.0E0) GO TO 100
C
C           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
C
            KPVT(K) = K
            INFO = K
         GO TO 190
  100    CONTINUE
         IF (KSTEP .EQ. 2) GO TO 140
C
C           1 X 1 PIVOT BLOCK.
C
            IF (.NOT.SWAP) GO TO 120
C
C              PERFORM AN INTERCHANGE.
C
               CALL CSWAP(IMAX,AP(IM+1),1,AP(IK+1),1)
               IMJ = IK + IMAX
               DO 110 JJ = IMAX, K
                  J = K + IMAX - JJ
                  JK = IK + J
                  T = AP(JK)
                  AP(JK) = AP(IMJ)
                  AP(IMJ) = T
                  IMJ = IMJ - (J - 1)
  110          CONTINUE
  120       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            IJ = IK - (K - 1)
            DO 130 JJ = 1, KM1
               J = K - JJ
               JK = IK + J
               MULK = -AP(JK)/AP(KK)
               T = MULK
               CALL CAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
               IJJ = IJ + J
               AP(JK) = MULK
               IJ = IJ - (J - 1)
  130       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = K
            IF (SWAP) KPVT(K) = IMAX
         GO TO 190
  140    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            KM1K = IK + K - 1
            IKM1 = IK - (K - 1)
            IF (.NOT.SWAP) GO TO 160
C
C              PERFORM AN INTERCHANGE.
C
               CALL CSWAP(IMAX,AP(IM+1),1,AP(IKM1+1),1)
               IMJ = IKM1 + IMAX
               DO 150 JJ = IMAX, KM1
                  J = KM1 + IMAX - JJ
                  JKM1 = IKM1 + J
                  T = AP(JKM1)
                  AP(JKM1) = AP(IMJ)
                  AP(IMJ) = T
                  IMJ = IMJ - (J - 1)
  150          CONTINUE
               T = AP(KM1K)
               AP(KM1K) = AP(IMK)
               AP(IMK) = T
  160       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            KM2 = K - 2
            IF (KM2 .EQ. 0) GO TO 180
               AK = AP(KK)/AP(KM1K)
               KM1KM1 = IKM1 + K - 1
               AKM1 = AP(KM1KM1)/AP(KM1K)
               DENOM = 1.0E0 - AK*AKM1
               IJ = IK - (K - 1) - (K - 2)
               DO 170 JJ = 1, KM2
                  J = KM1 - JJ
                  JK = IK + J
                  BK = AP(JK)/AP(KM1K)
                  JKM1 = IKM1 + J
                  BKM1 = AP(JKM1)/AP(KM1K)
                  MULK = (AKM1*BK - BKM1)/DENOM
                  MULKM1 = (AK*BKM1 - BK)/DENOM
                  T = MULK
                  CALL CAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
                  T = MULKM1
                  CALL CAXPY(J,T,AP(IKM1+1),1,AP(IJ+1),1)
                  AP(JK) = MULK
                  AP(JKM1) = MULKM1
                  IJJ = IJ + J
                  IJ = IJ - (J - 1)
  170          CONTINUE
  180       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = 1 - K
            IF (SWAP) KPVT(K) = -IMAX
            KPVT(K-1) = KPVT(K)
  190    CONTINUE
         IK = IK - (K - 1)
         IF (KSTEP .EQ. 2) IK = IK - (K - 2)
         K = K - KSTEP
      GO TO 10
  200 CONTINUE
      RETURN
      END
      SUBROUTINE CSPSL(AP,N,KPVT,B)
      INTEGER N,KPVT(1)
      COMPLEX AP(1),B(1)
C
C     CSISL SOLVES THE COMPLEX SYMMETRIC SYSTEM
C     A * X = B
C     USING THE FACTORS COMPUTED BY CSPFA.
C
C     ON ENTRY
C
C        AP      COMPLEX(N*(N+1)/2)
C                THE OUTPUT FROM CSPFA.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        KPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM CSPFA.
C
C        B       COMPLEX(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO MAY OCCUR IF  CSPCO  HAS SET RCOND .EQ. 0.0
C        OR  CSPFA  HAS SET INFO .NE. 0  .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL CSPFA(AP,N,KPVT,INFO)
C           IF (INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL CSPSL(AP,N,KPVT,C(1,J))
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CDOTU
C     FORTRAN IABS
C
C     INTERNAL VARIABLES.
C
      COMPLEX AK,AKM1,BK,BKM1,CDOTU,DENOM,TEMP
      INTEGER IK,IKM1,IKP1,K,KK,KM1K,KM1KM1,KP
C
C     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
C     D INVERSE TO B.
C
      K = N
      IK = (N*(N - 1))/2
   10 IF (K .EQ. 0) GO TO 80
         KK = IK + K
         IF (KPVT(K) .LT. 0) GO TO 40
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 30
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 20
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
   20          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL CAXPY(K-1,B(K),AP(IK+1),1,B(1),1)
   30       CONTINUE
C
C           APPLY D INVERSE.
C
            B(K) = B(K)/AP(KK)
            K = K - 1
            IK = IK - K
         GO TO 70
   40    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IKM1 = IK - (K - 1)
            IF (K .EQ. 2) GO TO 60
               KP = IABS(KPVT(K))
               IF (KP .EQ. K - 1) GO TO 50
C
C                 INTERCHANGE.
C
                  TEMP = B(K-1)
                  B(K-1) = B(KP)
                  B(KP) = TEMP
   50          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL CAXPY(K-2,B(K),AP(IK+1),1,B(1),1)
               CALL CAXPY(K-2,B(K-1),AP(IKM1+1),1,B(1),1)
   60       CONTINUE
C
C           APPLY D INVERSE.
C
            KM1K = IK + K - 1
            KK = IK + K
            AK = AP(KK)/AP(KM1K)
            KM1KM1 = IKM1 + K - 1
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = B(K)/AP(KM1K)
            BKM1 = B(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0E0
            B(K) = (AKM1*BK - BKM1)/DENOM
            B(K-1) = (AK*BKM1 - BK)/DENOM
            K = K - 2
            IK = IK - (K + 1) - K
   70    CONTINUE
      GO TO 10
   80 CONTINUE
C
C     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
C
      K = 1
      IK = 0
   90 IF (K .GT. N) GO TO 160
         IF (KPVT(K) .LT. 0) GO TO 120
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 110
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + CDOTU(K-1,AP(IK+1),1,B(1),1)
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 100
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  100          CONTINUE
  110       CONTINUE
            IK = IK + K
            K = K + 1
         GO TO 150
  120    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 140
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + CDOTU(K-1,AP(IK+1),1,B(1),1)
               IKP1 = IK + K
               B(K+1) = B(K+1) + CDOTU(K-1,AP(IKP1+1),1,B(1),1)
               KP = IABS(KPVT(K))
               IF (KP .EQ. K) GO TO 130
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  130          CONTINUE
  140       CONTINUE
            IK = IK + K + K + 1
            K = K + 2
  150    CONTINUE
      GO TO 90
  160 CONTINUE
      RETURN
      END
      SUBROUTINE CSPDI(AP,N,KPVT,DET,WORK,JOB)
      INTEGER N,JOB
      COMPLEX AP(1),WORK(1),DET(2)
      INTEGER KPVT(1)
C
C     CSPDI COMPUTES THE DETERMINANT AND INVERSE
C     OF A COMPLEX SYMMETRIC MATRIX USING THE FACTORS FROM CSPFA,
C     WHERE THE MATRIX IS STORED IN PACKED FORM.
C
C     ON ENTRY
C
C        AP      COMPLEX (N*(N+1)/2)
C                THE OUTPUT FROM CSPFA.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX A.
C
C        KPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM CSPFA.
C
C        WORK    COMPLEX(N)
C                WORK VECTOR.  CONTENTS IGNORED.
C
C        JOB     INTEGER
C                JOB HAS THE DECIMAL EXPANSION  AB  WHERE
C                   IF  B .NE. 0, THE INVERSE IS COMPUTED,
C                   IF  A .NE. 0, THE DETERMINANT IS COMPUTED,
C
C                FOR EXAMPLE, JOB = 11  GIVES BOTH.
C
C     ON RETURN
C
C        VARIABLES NOT REQUESTED BY JOB ARE NOT USED.
C
C        AP     CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF
C               THE ORIGINAL MATRIX, STORED IN PACKED FORM.
C               THE COLUMNS OF THE UPPER TRIANGLE ARE STORED
C               SEQUENTIALLY IN A ONE-DIMENSIONAL ARRAY.
C
C        DET    COMPLEX(2)
C               DETERMINANT OF ORIGINAL MATRIX.
C               DETERMINANT = DET(1) * 10.0**DET(2)
C               WITH 1.0 .LE. ABS(DET(1)) .LT. 10.0
C               OR DET(1) = 0.0.
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INVERSE IS REQUESTED
C        AND  CSPCO  HAS SET RCOND .EQ. 0.0
C        OR  CSPFA  HAS SET  INFO .NE. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CCOPY,CDOTU,CSWAP
C     FORTRAN ABS,CMPLX,IABS,MOD,REAL
C
C     INTERNAL VARIABLES.
C
      COMPLEX AK,AKKP1,AKP1,CDOTU,D,T,TEMP
      REAL TEN
      INTEGER IJ,IK,IKP1,IKS,J,JB,JK,JKP1
      INTEGER K,KK,KKP1,KM1,KS,KSJ,KSKP1,KSTEP
      LOGICAL NOINV,NODET
C
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
      NOINV = MOD(JOB,10) .EQ. 0
      NODET = MOD(JOB,100)/10 .EQ. 0
C
      IF (NODET) GO TO 110
         DET(1) = (1.0E0,0.0E0)
         DET(2) = (0.0E0,0.0E0)
         TEN = 10.0E0
         T = (0.0E0,0.0E0)
         IK = 0
         DO 100 K = 1, N
            KK = IK + K
            D = AP(KK)
C
C           CHECK IF 1 BY 1
C
            IF (KPVT(K) .GT. 0) GO TO 30
C
C              2 BY 2 BLOCK
C              USE DET (D  T)  =  (D/T * C - T) * T
C                      (T  C)
C              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
C              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
C
               IF (CABS1(T) .NE. 0.0E0) GO TO 10
                  IKP1 = IK + K
                  KKP1 = IKP1 + K
                  T = AP(KKP1)
                  D = (D/T)*AP(KKP1+1) - T
               GO TO 20
   10          CONTINUE
                  D = T
                  T = (0.0E0,0.0E0)
   20          CONTINUE
   30       CONTINUE
C
            IF (NODET) GO TO 90
               DET(1) = D*DET(1)
               IF (CABS1(DET(1)) .EQ. 0.0E0) GO TO 80
   40             IF (CABS1(DET(1)) .GE. 1.0E0) GO TO 50
                     DET(1) = CMPLX(TEN,0.0E0)*DET(1)
                     DET(2) = DET(2) - (1.0E0,0.0E0)
                  GO TO 40
   50             CONTINUE
   60             IF (CABS1(DET(1)) .LT. TEN) GO TO 70
                     DET(1) = DET(1)/CMPLX(TEN,0.0E0)
                     DET(2) = DET(2) + (1.0E0,0.0E0)
                  GO TO 60
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
            IK = IK + K
  100    CONTINUE
  110 CONTINUE
C
C     COMPUTE INVERSE(A)
C
      IF (NOINV) GO TO 240
         K = 1
         IK = 0
  120    IF (K .GT. N) GO TO 230
            KM1 = K - 1
            KK = IK + K
            IKP1 = IK + K
            IF (KPVT(K) .LT. 0) GO TO 150
C
C              1 BY 1
C
               AP(KK) = (1.0E0,0.0E0)/AP(KK)
               IF (KM1 .LT. 1) GO TO 140
                  CALL CCOPY(KM1,AP(IK+1),1,WORK,1)
                  IJ = 0
                  DO 130 J = 1, KM1
                     JK = IK + J
                     AP(JK) = CDOTU(J,AP(IJ+1),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
                     IJ = IJ + J
  130             CONTINUE
                  AP(KK) = AP(KK) + CDOTU(KM1,WORK,1,AP(IK+1),1)
  140          CONTINUE
               KSTEP = 1
            GO TO 190
  150       CONTINUE
C
C              2 BY 2
C
               KKP1 = IKP1 + K
               T = AP(KKP1)
               AK = AP(KK)/T
               AKP1 = AP(KKP1+1)/T
               AKKP1 = AP(KKP1)/T
               D = T*(AK*AKP1 - (1.0E0,0.0E0))
               AP(KK) = AKP1/D
               AP(KKP1+1) = AK/D
               AP(KKP1) = -AKKP1/D
               IF (KM1 .LT. 1) GO TO 180
                  CALL CCOPY(KM1,AP(IKP1+1),1,WORK,1)
                  IJ = 0
                  DO 160 J = 1, KM1
                     JKP1 = IKP1 + J
                     AP(JKP1) = CDOTU(J,AP(IJ+1),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IKP1+1),1)
                     IJ = IJ + J
  160             CONTINUE
                  AP(KKP1+1) = AP(KKP1+1)
     *                         + CDOTU(KM1,WORK,1,AP(IKP1+1),1)
                  AP(KKP1) = AP(KKP1)
     *                       + CDOTU(KM1,AP(IK+1),1,AP(IKP1+1),1)
                  CALL CCOPY(KM1,AP(IK+1),1,WORK,1)
                  IJ = 0
                  DO 170 J = 1, KM1
                     JK = IK + J
                     AP(JK) = CDOTU(J,AP(IJ+1),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
                     IJ = IJ + J
  170             CONTINUE
                  AP(KK) = AP(KK) + CDOTU(KM1,WORK,1,AP(IK+1),1)
  180          CONTINUE
               KSTEP = 2
  190       CONTINUE
C
C           SWAP
C
            KS = IABS(KPVT(K))
            IF (KS .EQ. K) GO TO 220
               IKS = (KS*(KS - 1))/2
               CALL CSWAP(KS,AP(IKS+1),1,AP(IK+1),1)
               KSJ = IK + KS
               DO 200 JB = KS, K
                  J = K + KS - JB
                  JK = IK + J
                  TEMP = AP(JK)
                  AP(JK) = AP(KSJ)
                  AP(KSJ) = TEMP
                  KSJ = KSJ - (J - 1)
  200          CONTINUE
               IF (KSTEP .EQ. 1) GO TO 210
                  KSKP1 = IKP1 + KS
                  TEMP = AP(KSKP1)
                  AP(KSKP1) = AP(KKP1)
                  AP(KKP1) = TEMP
  210          CONTINUE
  220       CONTINUE
            IK = IK + K
            IF (KSTEP .EQ. 2) IK = IK + K + 1
            K = K + KSTEP
         GO TO 120
  230    CONTINUE
  240 CONTINUE
      RETURN
      END
      SUBROUTINE CHICO(A,LDA,N,KPVT,RCOND,Z)
      INTEGER LDA,N,KPVT(1)
      COMPLEX A(LDA,1),Z(1)
      REAL RCOND
C
C     CHICO FACTORS A COMPLEX HERMITIAN MATRIX BY ELIMINATION WITH
C     SYMMETRIC PIVOTING AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, CHIFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW CHICO BY CHISL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW CHICO BY CHISL.
C     TO COMPUTE  INVERSE(A) , FOLLOW CHICO BY CHIDI.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW CHICO BY CHIDI.
C     TO COMPUTE  INERTIA(A), FOLLOW CHICO BY CHIDI.
C
C     ON ENTRY
C
C        A       COMPLEX(LDA, N)
C                THE HERMITIAN MATRIX TO BE FACTORED.
C                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     OUTPUT
C
C        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
C                WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*CTRANS(U)
C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
C                UPPER TRIANGULAR MATRICES , CTRANS(U) IS THE
C                CONJUGATE TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
C
C        KPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       COMPLEX(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK CHIFA
C     BLAS CAXPY,CDOTC,CSSCAL,SCASUM
C     FORTRAN ABS,AIMAG,AMAX1,CMPLX,CONJG,IABS,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX AK,AKM1,BK,BKM1,CDOTC,DENOM,EK,T
      REAL ANORM,S,SCASUM,YNORM
      INTEGER I,INFO,J,JM1,K,KP,KPS,KS
C
      COMPLEX ZDUM,ZDUM2,CSIGN1
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
      CSIGN1(ZDUM,ZDUM2) = CABS1(ZDUM)*(ZDUM2/CABS1(ZDUM2))
C
C     FIND NORM OF A USING ONLY UPPER HALF
C
      DO 30 J = 1, N
         Z(J) = CMPLX(SCASUM(J,A(1,J),1),0.0E0)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = CMPLX(REAL(Z(I))+CABS1(A(I,J)),0.0E0)
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0E0
      DO 40 J = 1, N
         ANORM = AMAX1(ANORM,REAL(Z(J)))
   40 CONTINUE
C
C     FACTOR
C
      CALL CHIFA(A,LDA,N,KPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C     SOLVE U*D*W = E
C
      EK = (1.0E0,0.0E0)
      DO 50 J = 1, N
         Z(J) = (0.0E0,0.0E0)
   50 CONTINUE
      K = N
   60 IF (K .EQ. 0) GO TO 120
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         KP = IABS(KPVT(K))
         KPS = K + 1 - KS
         IF (KP .EQ. KPS) GO TO 70
            T = Z(KPS)
            Z(KPS) = Z(KP)
            Z(KP) = T
   70    CONTINUE
         IF (CABS1(Z(K)) .NE. 0.0E0) EK = CSIGN1(EK,Z(K))
         Z(K) = Z(K) + EK
         CALL CAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)
         IF (KS .EQ. 1) GO TO 80
            IF (CABS1(Z(K-1)) .NE. 0.0E0) EK = CSIGN1(EK,Z(K-1))
            Z(K-1) = Z(K-1) + EK
            CALL CAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)
   80    CONTINUE
         IF (KS .EQ. 2) GO TO 100
            IF (CABS1(Z(K)) .LE. CABS1(A(K,K))) GO TO 90
               S = CABS1(A(K,K))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
               EK = CMPLX(S,0.0E0)*EK
   90       CONTINUE
            IF (CABS1(A(K,K)) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)
            IF (CABS1(A(K,K)) .EQ. 0.0E0) Z(K) = (1.0E0,0.0E0)
         GO TO 110
  100    CONTINUE
            AK = A(K,K)/CONJG(A(K-1,K))
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = Z(K)/CONJG(A(K-1,K))
            BKM1 = Z(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0E0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  110    CONTINUE
         K = K - KS
      GO TO 60
  120 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
C     SOLVE CTRANS(U)*Y = W
C
      K = 1
  130 IF (K .GT. N) GO TO 160
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 150
            Z(K) = Z(K) + CDOTC(K-1,A(1,K),1,Z(1),1)
            IF (KS .EQ. 2)
     *         Z(K+1) = Z(K+1) + CDOTC(K-1,A(1,K+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 140
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  140       CONTINUE
  150    CONTINUE
         K = K + KS
      GO TO 130
  160 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
      YNORM = 1.0E0
C
C     SOLVE U*D*V = Y
C
      K = N
  170 IF (K .EQ. 0) GO TO 230
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. KS) GO TO 190
            KP = IABS(KPVT(K))
            KPS = K + 1 - KS
            IF (KP .EQ. KPS) GO TO 180
               T = Z(KPS)
               Z(KPS) = Z(KP)
               Z(KP) = T
  180       CONTINUE
            CALL CAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)
            IF (KS .EQ. 2) CALL CAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)
  190    CONTINUE
         IF (KS .EQ. 2) GO TO 210
            IF (CABS1(Z(K)) .LE. CABS1(A(K,K))) GO TO 200
               S = CABS1(A(K,K))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
               YNORM = S*YNORM
  200       CONTINUE
            IF (CABS1(A(K,K)) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)
            IF (CABS1(A(K,K)) .EQ. 0.0E0) Z(K) = (1.0E0,0.0E0)
         GO TO 220
  210    CONTINUE
            AK = A(K,K)/CONJG(A(K-1,K))
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = Z(K)/CONJG(A(K-1,K))
            BKM1 = Z(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0E0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  220    CONTINUE
         K = K - KS
      GO TO 170
  230 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE CTRANS(U)*Z = V
C
      K = 1
  240 IF (K .GT. N) GO TO 270
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 260
            Z(K) = Z(K) + CDOTC(K-1,A(1,K),1,Z(1),1)
            IF (KS .EQ. 2)
     *         Z(K+1) = Z(K+1) + CDOTC(K-1,A(1,K+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 250
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  250       CONTINUE
  260    CONTINUE
         K = K + KS
      GO TO 240
  270 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
      SUBROUTINE CHIFA(A,LDA,N,KPVT,INFO)
      INTEGER LDA,N,KPVT(1),INFO
      COMPLEX A(LDA,1)
C
C     CHIFA FACTORS A COMPLEX HERMITIAN MATRIX BY ELIMINATION
C     WITH SYMMETRIC PIVOTING.
C
C     TO SOLVE  A*X = B , FOLLOW CHIFA BY CHISL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW CHIFA BY CHISL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW CHIFA BY CHIDI.
C     TO COMPUTE  INERTIA(A) , FOLLOW CHIFA BY CHIDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW CHIFA BY CHIDI.
C
C     ON ENTRY
C
C        A       COMPLEX(LDA,N)
C                THE HERMITIAN MATRIX TO BE FACTORED.
C                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
C                WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*CTRANS(U)
C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
C                UPPER TRIANGULAR MATRICES , CTRANS(U) IS THE
C                CONJUGATE TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
C
C        KPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS
C                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,
C                     BUT IT DOES INDICATE THAT CHISL OR CHIDI MAY
C                     DIVIDE BY ZERO IF CALLED.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CSWAP,ICAMAX
C     FORTRAN ABS,AIMAG,AMAX1,CMPLX,CONJG,REAL,SQRT
C
C     INTERNAL VARIABLES
C
      COMPLEX AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
      REAL ABSAKK,ALPHA,COLMAX,ROWMAX
      INTEGER IMAX,IMAXP1,J,JJ,JMAX,K,KM1,KM2,KSTEP,ICAMAX
      LOGICAL SWAP
C
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
C     INITIALIZE
C
C     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
      ALPHA = (1.0E0 + SQRT(17.0E0))/8.0E0
C
      INFO = 0
C
C     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
C
      K = N
   10 CONTINUE
C
C        LEAVE THE LOOP IF K=0 OR K=1.
C
C     ...EXIT
         IF (K .EQ. 0) GO TO 200
         IF (K .GT. 1) GO TO 20
            KPVT(1) = 1
            IF (CABS1(A(1,1)) .EQ. 0.0E0) INFO = 1
C     ......EXIT
            GO TO 200
   20    CONTINUE
C
C        THIS SECTION OF CODE DETERMINES THE KIND OF
C        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
C        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
C        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
C        REQUIRED.
C
         KM1 = K - 1
         ABSAKK = CABS1(A(K,K))
C
C        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C        COLUMN K.
C
         IMAX = ICAMAX(K-1,A(1,K),1)
         COLMAX = CABS1(A(IMAX,K))
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30
            KSTEP = 1
            SWAP = .FALSE.
         GO TO 90
   30    CONTINUE
C
C           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C           ROW IMAX.
C
            ROWMAX = 0.0E0
            IMAXP1 = IMAX + 1
            DO 40 J = IMAXP1, K
               ROWMAX = AMAX1(ROWMAX,CABS1(A(IMAX,J)))
   40       CONTINUE
            IF (IMAX .EQ. 1) GO TO 50
               JMAX = ICAMAX(IMAX-1,A(1,IMAX),1)
               ROWMAX = AMAX1(ROWMAX,CABS1(A(JMAX,IMAX)))
   50       CONTINUE
            IF (CABS1(A(IMAX,IMAX)) .LT. ALPHA*ROWMAX) GO TO 60
               KSTEP = 1
               SWAP = .TRUE.
            GO TO 80
   60       CONTINUE
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70
               KSTEP = 1
               SWAP = .FALSE.
            GO TO 80
   70       CONTINUE
               KSTEP = 2
               SWAP = IMAX .NE. KM1
   80       CONTINUE
   90    CONTINUE
         IF (AMAX1(ABSAKK,COLMAX) .NE. 0.0E0) GO TO 100
C
C           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
C
            KPVT(K) = K
            INFO = K
         GO TO 190
  100    CONTINUE
         IF (KSTEP .EQ. 2) GO TO 140
C
C           1 X 1 PIVOT BLOCK.
C
            IF (.NOT.SWAP) GO TO 120
C
C              PERFORM AN INTERCHANGE.
C
               CALL CSWAP(IMAX,A(1,IMAX),1,A(1,K),1)
               DO 110 JJ = IMAX, K
                  J = K + IMAX - JJ
                  T = CONJG(A(J,K))
                  A(J,K) = CONJG(A(IMAX,J))
                  A(IMAX,J) = T
  110          CONTINUE
  120       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            DO 130 JJ = 1, KM1
               J = K - JJ
               MULK = -A(J,K)/A(K,K)
               T = CONJG(MULK)
               CALL CAXPY(J,T,A(1,K),1,A(1,J),1)
               A(J,J) = CMPLX(REAL(A(J,J)),0.0E0)
               A(J,K) = MULK
  130       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = K
            IF (SWAP) KPVT(K) = IMAX
         GO TO 190
  140    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (.NOT.SWAP) GO TO 160
C
C              PERFORM AN INTERCHANGE.
C
               CALL CSWAP(IMAX,A(1,IMAX),1,A(1,K-1),1)
               DO 150 JJ = IMAX, KM1
                  J = KM1 + IMAX - JJ
                  T = CONJG(A(J,K-1))
                  A(J,K-1) = CONJG(A(IMAX,J))
                  A(IMAX,J) = T
  150          CONTINUE
               T = A(K-1,K)
               A(K-1,K) = A(IMAX,K)
               A(IMAX,K) = T
  160       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            KM2 = K - 2
            IF (KM2 .EQ. 0) GO TO 180
               AK = A(K,K)/A(K-1,K)
               AKM1 = A(K-1,K-1)/CONJG(A(K-1,K))
               DENOM = 1.0E0 - AK*AKM1
               DO 170 JJ = 1, KM2
                  J = KM1 - JJ
                  BK = A(J,K)/A(K-1,K)
                  BKM1 = A(J,K-1)/CONJG(A(K-1,K))
                  MULK = (AKM1*BK - BKM1)/DENOM
                  MULKM1 = (AK*BKM1 - BK)/DENOM
                  T = CONJG(MULK)
                  CALL CAXPY(J,T,A(1,K),1,A(1,J),1)
                  T = CONJG(MULKM1)
                  CALL CAXPY(J,T,A(1,K-1),1,A(1,J),1)
                  A(J,K) = MULK
                  A(J,K-1) = MULKM1
                  A(J,J) = CMPLX(REAL(A(J,J)),0.0E0)
  170          CONTINUE
  180       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = 1 - K
            IF (SWAP) KPVT(K) = -IMAX
            KPVT(K-1) = KPVT(K)
  190    CONTINUE
         K = K - KSTEP
      GO TO 10
  200 CONTINUE
      RETURN
      END
      SUBROUTINE CHISL(A,LDA,N,KPVT,B)
      INTEGER LDA,N,KPVT(1)
      COMPLEX A(LDA,1),B(1)
C
C     CHISL SOLVES THE COMPLEX HERMITIAN SYSTEM
C     A * X = B
C     USING THE FACTORS COMPUTED BY CHIFA.
C
C     ON ENTRY
C
C        A       COMPLEX(LDA,N)
C                THE OUTPUT FROM CHIFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        KPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM CHIFA.
C
C        B       COMPLEX(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO MAY OCCUR IF  CHICO  HAS SET RCOND .EQ. 0.0
C        OR  CHIFA  HAS SET INFO .NE. 0  .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL CHIFA(A,LDA,N,KPVT,INFO)
C           IF (INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL CHISL(A,LDA,N,KPVT,C(1,J))
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CDOTC
C     FORTRAN CONJG,IABS
C
C     INTERNAL VARIABLES.
C
      COMPLEX AK,AKM1,BK,BKM1,CDOTC,DENOM,TEMP
      INTEGER K,KP
C
C     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
C     D INVERSE TO B.
C
      K = N
   10 IF (K .EQ. 0) GO TO 80
         IF (KPVT(K) .LT. 0) GO TO 40
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 30
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 20
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
   20          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL CAXPY(K-1,B(K),A(1,K),1,B(1),1)
   30       CONTINUE
C
C           APPLY D INVERSE.
C
            B(K) = B(K)/A(K,K)
            K = K - 1
         GO TO 70
   40    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (K .EQ. 2) GO TO 60
               KP = IABS(KPVT(K))
               IF (KP .EQ. K - 1) GO TO 50
C
C                 INTERCHANGE.
C
                  TEMP = B(K-1)
                  B(K-1) = B(KP)
                  B(KP) = TEMP
   50          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL CAXPY(K-2,B(K),A(1,K),1,B(1),1)
               CALL CAXPY(K-2,B(K-1),A(1,K-1),1,B(1),1)
   60       CONTINUE
C
C           APPLY D INVERSE.
C
            AK = A(K,K)/CONJG(A(K-1,K))
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = B(K)/CONJG(A(K-1,K))
            BKM1 = B(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0E0
            B(K) = (AKM1*BK - BKM1)/DENOM
            B(K-1) = (AK*BKM1 - BK)/DENOM
            K = K - 2
   70    CONTINUE
      GO TO 10
   80 CONTINUE
C
C     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
C
      K = 1
   90 IF (K .GT. N) GO TO 160
         IF (KPVT(K) .LT. 0) GO TO 120
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 110
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + CDOTC(K-1,A(1,K),1,B(1),1)
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 100
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  100          CONTINUE
  110       CONTINUE
            K = K + 1
         GO TO 150
  120    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 140
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + CDOTC(K-1,A(1,K),1,B(1),1)
               B(K+1) = B(K+1) + CDOTC(K-1,A(1,K+1),1,B(1),1)
               KP = IABS(KPVT(K))
               IF (KP .EQ. K) GO TO 130
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  130          CONTINUE
  140       CONTINUE
            K = K + 2
  150    CONTINUE
      GO TO 90
  160 CONTINUE
      RETURN
      END
      SUBROUTINE CHIDI(A,LDA,N,KPVT,DET,INERT,WORK,JOB)
      INTEGER LDA,N,JOB
      COMPLEX A(LDA,1),WORK(1)
      REAL DET(2)
      INTEGER KPVT(1),INERT(3)
C
C     CHIDI COMPUTES THE DETERMINANT, INERTIA AND INVERSE
C     OF A COMPLEX HERMITIAN MATRIX USING THE FACTORS FROM CHIFA.
C
C     ON ENTRY
C
C        A       COMPLEX(LDA,N)
C                THE OUTPUT FROM CHIFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY A.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX A.
C
C        KPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM CHIFA.
C
C        WORK    COMPLEX(N)
C                WORK VECTOR.  CONTENTS DESTROYED.
C
C        JOB     INTEGER
C                JOB HAS THE DECIMAL EXPANSION  ABC  WHERE
C                   IF  C .NE. 0, THE INVERSE IS COMPUTED,
C                   IF  B .NE. 0, THE DETERMINANT IS COMPUTED,
C                   IF  A .NE. 0, THE INERTIA IS COMPUTED.
C
C                FOR EXAMPLE, JOB = 111  GIVES ALL THREE.
C
C     ON RETURN
C
C        VARIABLES NOT REQUESTED BY JOB ARE NOT USED.
C
C        A      CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF
C               THE ORIGINAL MATRIX.  THE STRICT LOWER TRIANGLE
C               IS NEVER REFERENCED.
C
C        DET    REAL(2)
C               DETERMINANT OF ORIGINAL MATRIX.
C               DETERMINANT = DET(1) * 10.0**DET(2)
C               WITH 1.0 .LE. ABS(DET(1)) .LT. 10.0
C               OR DET(1) = 0.0.
C
C        INERT  INTEGER(3)
C               THE INERTIA OF THE ORIGINAL MATRIX.
C               INERT(1)  =  NUMBER OF POSITIVE EIGENVALUES.
C               INERT(2)  =  NUMBER OF NEGATIVE EIGENVALUES.
C               INERT(3)  =  NUMBER OF ZERO EIGENVALUES.
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO MAY OCCUR IF THE INVERSE IS REQUESTED
C        AND  CHICO  HAS SET RCOND .EQ. 0.0
C        OR  CHIFA  HAS SET  INFO .NE. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CCOPY,CDOTC,CSWAP
C     FORTRAN ABS,CABS,CMPLX,CONJG,IABS,MOD,REAL
C
C     INTERNAL VARIABLES.
C
      COMPLEX AKKP1,CDOTC,TEMP
      REAL TEN,D,T,AK,AKP1
      INTEGER J,JB,K,KM1,KS,KSTEP
      LOGICAL NOINV,NODET,NOERT
C
      NOINV = MOD(JOB,10) .EQ. 0
      NODET = MOD(JOB,100)/10 .EQ. 0
      NOERT = MOD(JOB,1000)/100 .EQ. 0
C
      IF (NODET .AND. NOERT) GO TO 140
         IF (NOERT) GO TO 10
            INERT(1) = 0
            INERT(2) = 0
            INERT(3) = 0
   10    CONTINUE
         IF (NODET) GO TO 20
            DET(1) = 1.0E0
            DET(2) = 0.0E0
            TEN = 10.0E0
   20    CONTINUE
         T = 0.0E0
         DO 130 K = 1, N
            D = REAL(A(K,K))
C
C           CHECK IF 1 BY 1
C
            IF (KPVT(K) .GT. 0) GO TO 50
C
C              2 BY 2 BLOCK
C              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = CABS(S)
C                      (S  C)
C              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
C              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
C
               IF (T .NE. 0.0E0) GO TO 30
                  T = CABS(A(K,K+1))
                  D = (D/T)*REAL(A(K+1,K+1)) - T
               GO TO 40
   30          CONTINUE
                  D = T
                  T = 0.0E0
   40          CONTINUE
   50       CONTINUE
C
            IF (NOERT) GO TO 60
               IF (D .GT. 0.0E0) INERT(1) = INERT(1) + 1
               IF (D .LT. 0.0E0) INERT(2) = INERT(2) + 1
               IF (D .EQ. 0.0E0) INERT(3) = INERT(3) + 1
   60       CONTINUE
C
            IF (NODET) GO TO 120
               DET(1) = D*DET(1)
               IF (DET(1) .EQ. 0.0E0) GO TO 110
   70             IF (ABS(DET(1)) .GE. 1.0E0) GO TO 80
                     DET(1) = TEN*DET(1)
                     DET(2) = DET(2) - 1.0E0
                  GO TO 70
   80             CONTINUE
   90             IF (ABS(DET(1)) .LT. TEN) GO TO 100
                     DET(1) = DET(1)/TEN
                     DET(2) = DET(2) + 1.0E0
                  GO TO 90
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
C
C     COMPUTE INVERSE(A)
C
      IF (NOINV) GO TO 270
         K = 1
  150    IF (K .GT. N) GO TO 260
            KM1 = K - 1
            IF (KPVT(K) .LT. 0) GO TO 180
C
C              1 BY 1
C
               A(K,K) = CMPLX(1.0E0/REAL(A(K,K)),0.0E0)
               IF (KM1 .LT. 1) GO TO 170
                  CALL CCOPY(KM1,A(1,K),1,WORK,1)
                  DO 160 J = 1, KM1
                     A(J,K) = CDOTC(J,A(1,J),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  160             CONTINUE
                  A(K,K) = A(K,K)
     *                     + CMPLX(REAL(CDOTC(KM1,WORK,1,A(1,K),1)),
     *                             0.0E0)
  170          CONTINUE
               KSTEP = 1
            GO TO 220
  180       CONTINUE
C
C              2 BY 2
C
               T = CABS(A(K,K+1))
               AK = REAL(A(K,K))/T
               AKP1 = REAL(A(K+1,K+1))/T
               AKKP1 = A(K,K+1)/T
               D = T*(AK*AKP1 - 1.0E0)
               A(K,K) = CMPLX(AKP1/D,0.0E0)
               A(K+1,K+1) = CMPLX(AK/D,0.0E0)
               A(K,K+1) = -AKKP1/D
               IF (KM1 .LT. 1) GO TO 210
                  CALL CCOPY(KM1,A(1,K+1),1,WORK,1)
                  DO 190 J = 1, KM1
                     A(J,K+1) = CDOTC(J,A(1,J),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),A(1,J),1,A(1,K+1),1)
  190             CONTINUE
                  A(K+1,K+1) = A(K+1,K+1)
     *                         + CMPLX(REAL(CDOTC(KM1,WORK,1,A(1,K+1),
     *                                            1)),0.0E0)
                  A(K,K+1) = A(K,K+1) + CDOTC(KM1,A(1,K),1,A(1,K+1),1)
                  CALL CCOPY(KM1,A(1,K),1,WORK,1)
                  DO 200 J = 1, KM1
                     A(J,K) = CDOTC(J,A(1,J),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  200             CONTINUE
                  A(K,K) = A(K,K)
     *                     + CMPLX(REAL(CDOTC(KM1,WORK,1,A(1,K),1)),
     *                             0.0E0)
  210          CONTINUE
               KSTEP = 2
  220       CONTINUE
C
C           SWAP
C
            KS = IABS(KPVT(K))
            IF (KS .EQ. K) GO TO 250
               CALL CSWAP(KS,A(1,KS),1,A(1,K),1)
               DO 230 JB = KS, K
                  J = K + KS - JB
                  TEMP = CONJG(A(J,K))
                  A(J,K) = CONJG(A(KS,J))
                  A(KS,J) = TEMP
  230          CONTINUE
               IF (KSTEP .EQ. 1) GO TO 240
                  TEMP = A(KS,K+1)
                  A(KS,K+1) = A(K,K+1)
                  A(K,K+1) = TEMP
  240          CONTINUE
  250       CONTINUE
            K = K + KSTEP
         GO TO 150
  260    CONTINUE
  270 CONTINUE
      RETURN
      END
      SUBROUTINE CHPCO(AP,N,KPVT,RCOND,Z)
      INTEGER N,KPVT(1)
      COMPLEX AP(1),Z(1)
      REAL RCOND
C
C     CHPCO FACTORS A COMPLEX HERMITIAN MATRIX STORED IN PACKED
C     FORM BY ELIMINATION WITH SYMMETRIC PIVOTING AND ESTIMATES
C     THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, CHPFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW CHPCO BY CHPSL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW CHPCO BY CHPSL.
C     TO COMPUTE  INVERSE(A) , FOLLOW CHPCO BY CHPDI.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW CHPCO BY CHPDI.
C     TO COMPUTE  INERTIA(A), FOLLOW CHPCO BY CHPDI.
C
C     ON ENTRY
C
C        AP      COMPLEX (N*(N+1)/2)
C                THE PACKED FORM OF A HERMITIAN MATRIX  A .  THE
C                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY
C                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
C                SEE COMMENTS BELOW FOR DETAILS.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     OUTPUT
C
C        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
C                WERE USED TO OBTAIN IT STORED IN PACKED FORM.
C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*CTRANS(U)
C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
C                UPPER TRIANGULAR MATRICES , CTRANS(U) IS THE
C                CONJUGATE TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
C
C        KPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       COMPLEX(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     PACKED STORAGE
C
C          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
C          TRIANGLE OF A HERMITIAN MATRIX.
C
C                K = 0
C                DO 20 J = 1, N
C                   DO 10 I = 1, J
C                      K = K + 1
C                      AP(K) = A(I,J)
C             10    CONTINUE
C             20 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK CHPFA
C     BLAS CAXPY,CDOTC,CSSCAL,SCASUM
C     FORTRAN ABS,AIMAG,AMAX1,CMPLX,CONJG,IABS,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX AK,AKM1,BK,BKM1,CDOTC,DENOM,EK,T
      REAL ANORM,S,SCASUM,YNORM
      INTEGER I,IJ,IK,IKM1,IKP1,INFO,J,JM1,J1
      INTEGER K,KK,KM1K,KM1KM1,KP,KPS,KS
C
      COMPLEX ZDUM,ZDUM2,CSIGN1
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
      CSIGN1(ZDUM,ZDUM2) = CABS1(ZDUM)*(ZDUM2/CABS1(ZDUM2))
C
C     FIND NORM OF A USING ONLY UPPER HALF
C
      J1 = 1
      DO 30 J = 1, N
         Z(J) = CMPLX(SCASUM(J,AP(J1),1),0.0E0)
         IJ = J1
         J1 = J1 + J
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = CMPLX(REAL(Z(I))+CABS1(AP(IJ)),0.0E0)
            IJ = IJ + 1
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0E0
      DO 40 J = 1, N
         ANORM = AMAX1(ANORM,REAL(Z(J)))
   40 CONTINUE
C
C     FACTOR
C
      CALL CHPFA(AP,N,KPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C     SOLVE U*D*W = E
C
      EK = (1.0E0,0.0E0)
      DO 50 J = 1, N
         Z(J) = (0.0E0,0.0E0)
   50 CONTINUE
      K = N
      IK = (N*(N - 1))/2
   60 IF (K .EQ. 0) GO TO 120
         KK = IK + K
         IKM1 = IK - (K - 1)
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         KP = IABS(KPVT(K))
         KPS = K + 1 - KS
         IF (KP .EQ. KPS) GO TO 70
            T = Z(KPS)
            Z(KPS) = Z(KP)
            Z(KP) = T
   70    CONTINUE
         IF (CABS1(Z(K)) .NE. 0.0E0) EK = CSIGN1(EK,Z(K))
         Z(K) = Z(K) + EK
         CALL CAXPY(K-KS,Z(K),AP(IK+1),1,Z(1),1)
         IF (KS .EQ. 1) GO TO 80
            IF (CABS1(Z(K-1)) .NE. 0.0E0) EK = CSIGN1(EK,Z(K-1))
            Z(K-1) = Z(K-1) + EK
            CALL CAXPY(K-KS,Z(K-1),AP(IKM1+1),1,Z(1),1)
   80    CONTINUE
         IF (KS .EQ. 2) GO TO 100
            IF (CABS1(Z(K)) .LE. CABS1(AP(KK))) GO TO 90
               S = CABS1(AP(KK))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
               EK = CMPLX(S,0.0E0)*EK
   90       CONTINUE
            IF (CABS1(AP(KK)) .NE. 0.0E0) Z(K) = Z(K)/AP(KK)
            IF (CABS1(AP(KK)) .EQ. 0.0E0) Z(K) = (1.0E0,0.0E0)
         GO TO 110
  100    CONTINUE
            KM1K = IK + K - 1
            KM1KM1 = IKM1 + K - 1
            AK = AP(KK)/CONJG(AP(KM1K))
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = Z(K)/CONJG(AP(KM1K))
            BKM1 = Z(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0E0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  110    CONTINUE
         K = K - KS
         IK = IK - K
         IF (KS .EQ. 2) IK = IK - (K + 1)
      GO TO 60
  120 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
C     SOLVE CTRANS(U)*Y = W
C
      K = 1
      IK = 0
  130 IF (K .GT. N) GO TO 160
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 150
            Z(K) = Z(K) + CDOTC(K-1,AP(IK+1),1,Z(1),1)
            IKP1 = IK + K
            IF (KS .EQ. 2)
     *         Z(K+1) = Z(K+1) + CDOTC(K-1,AP(IKP1+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 140
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  140       CONTINUE
  150    CONTINUE
         IK = IK + K
         IF (KS .EQ. 2) IK = IK + (K + 1)
         K = K + KS
      GO TO 130
  160 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
      YNORM = 1.0E0
C
C     SOLVE U*D*V = Y
C
      K = N
      IK = N*(N - 1)/2
  170 IF (K .EQ. 0) GO TO 230
         KK = IK + K
         IKM1 = IK - (K - 1)
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. KS) GO TO 190
            KP = IABS(KPVT(K))
            KPS = K + 1 - KS
            IF (KP .EQ. KPS) GO TO 180
               T = Z(KPS)
               Z(KPS) = Z(KP)
               Z(KP) = T
  180       CONTINUE
            CALL CAXPY(K-KS,Z(K),AP(IK+1),1,Z(1),1)
            IF (KS .EQ. 2) CALL CAXPY(K-KS,Z(K-1),AP(IKM1+1),1,Z(1),1)
  190    CONTINUE
         IF (KS .EQ. 2) GO TO 210
            IF (CABS1(Z(K)) .LE. CABS1(AP(KK))) GO TO 200
               S = CABS1(AP(KK))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
               YNORM = S*YNORM
  200       CONTINUE
            IF (CABS1(AP(KK)) .NE. 0.0E0) Z(K) = Z(K)/AP(KK)
            IF (CABS1(AP(KK)) .EQ. 0.0E0) Z(K) = (1.0E0,0.0E0)
         GO TO 220
  210    CONTINUE
            KM1K = IK + K - 1
            KM1KM1 = IKM1 + K - 1
            AK = AP(KK)/CONJG(AP(KM1K))
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = Z(K)/CONJG(AP(KM1K))
            BKM1 = Z(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0E0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  220    CONTINUE
         K = K - KS
         IK = IK - K
         IF (KS .EQ. 2) IK = IK - (K + 1)
      GO TO 170
  230 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE CTRANS(U)*Z = V
C
      K = 1
      IK = 0
  240 IF (K .GT. N) GO TO 270
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 260
            Z(K) = Z(K) + CDOTC(K-1,AP(IK+1),1,Z(1),1)
            IKP1 = IK + K
            IF (KS .EQ. 2)
     *         Z(K+1) = Z(K+1) + CDOTC(K-1,AP(IKP1+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 250
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  250       CONTINUE
  260    CONTINUE
         IK = IK + K
         IF (KS .EQ. 2) IK = IK + (K + 1)
         K = K + KS
      GO TO 240
  270 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
      SUBROUTINE CHPFA(AP,N,KPVT,INFO)
      INTEGER N,KPVT(1),INFO
      COMPLEX AP(1)
C
C     CHPFA FACTORS A COMPLEX HERMITIAN MATRIX STORED IN
C     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING.
C
C     TO SOLVE  A*X = B , FOLLOW CHPFA BY CHPSL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW CHPFA BY CHPSL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW CHPFA BY CHPDI.
C     TO COMPUTE  INERTIA(A) , FOLLOW CHPFA BY CHPDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW CHPFA BY CHPDI.
C
C     ON ENTRY
C
C        AP      COMPLEX (N*(N+1)/2)
C                THE PACKED FORM OF A HERMITIAN MATRIX  A .  THE
C                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY
C                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
C                SEE COMMENTS BELOW FOR DETAILS.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     OUTPUT
C
C        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
C                WERE USED TO OBTAIN IT STORED IN PACKED FORM.
C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*CTRANS(U)
C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
C                UPPER TRIANGULAR MATRICES , CTRANS(U) IS THE
C                CONJUGATE TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
C
C        KPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS
C                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,
C                     BUT IT DOES INDICATE THAT CHPSL OR CHPDI MAY
C                     DIVIDE BY ZERO IF CALLED.
C
C     PACKED STORAGE
C
C          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
C          TRIANGLE OF A HERMITIAN MATRIX.
C
C                K = 0
C                DO 20 J = 1, N
C                   DO 10 I = 1, J
C                      K = K + 1
C                      AP(K)  = A(I,J)
C             10    CONTINUE
C             20 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CSWAP,ICAMAX
C     FORTRAN ABS,AIMAG,AMAX1,CMPLX,CONJG,REAL,SQRT
C
C     INTERNAL VARIABLES
C
      COMPLEX AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
      REAL ABSAKK,ALPHA,COLMAX,ROWMAX
      INTEGER ICAMAX,IJ,IJJ,IK,IKM1,IM,IMAX,IMAXP1,IMIM,IMJ,IMK
      INTEGER J,JJ,JK,JKM1,JMAX,JMIM,K,KK,KM1,KM1K,KM1KM1,KM2,KSTEP
      LOGICAL SWAP
C
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
C     INITIALIZE
C
C     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
      ALPHA = (1.0E0 + SQRT(17.0E0))/8.0E0
C
      INFO = 0
C
C     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
C
      K = N
      IK = (N*(N - 1))/2
   10 CONTINUE
C
C        LEAVE THE LOOP IF K=0 OR K=1.
C
C     ...EXIT
         IF (K .EQ. 0) GO TO 200
         IF (K .GT. 1) GO TO 20
            KPVT(1) = 1
            IF (CABS1(AP(1)) .EQ. 0.0E0) INFO = 1
C     ......EXIT
            GO TO 200
   20    CONTINUE
C
C        THIS SECTION OF CODE DETERMINES THE KIND OF
C        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
C        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
C        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
C        REQUIRED.
C
         KM1 = K - 1
         KK = IK + K
         ABSAKK = CABS1(AP(KK))
C
C        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C        COLUMN K.
C
         IMAX = ICAMAX(K-1,AP(IK+1),1)
         IMK = IK + IMAX
         COLMAX = CABS1(AP(IMK))
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30
            KSTEP = 1
            SWAP = .FALSE.
         GO TO 90
   30    CONTINUE
C
C           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C           ROW IMAX.
C
            ROWMAX = 0.0E0
            IMAXP1 = IMAX + 1
            IM = IMAX*(IMAX - 1)/2
            IMJ = IM + 2*IMAX
            DO 40 J = IMAXP1, K
               ROWMAX = AMAX1(ROWMAX,CABS1(AP(IMJ)))
               IMJ = IMJ + J
   40       CONTINUE
            IF (IMAX .EQ. 1) GO TO 50
               JMAX = ICAMAX(IMAX-1,AP(IM+1),1)
               JMIM = JMAX + IM
               ROWMAX = AMAX1(ROWMAX,CABS1(AP(JMIM)))
   50       CONTINUE
            IMIM = IMAX + IM
            IF (CABS1(AP(IMIM)) .LT. ALPHA*ROWMAX) GO TO 60
               KSTEP = 1
               SWAP = .TRUE.
            GO TO 80
   60       CONTINUE
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70
               KSTEP = 1
               SWAP = .FALSE.
            GO TO 80
   70       CONTINUE
               KSTEP = 2
               SWAP = IMAX .NE. KM1
   80       CONTINUE
   90    CONTINUE
         IF (AMAX1(ABSAKK,COLMAX) .NE. 0.0E0) GO TO 100
C
C           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
C
            KPVT(K) = K
            INFO = K
         GO TO 190
  100    CONTINUE
         IF (KSTEP .EQ. 2) GO TO 140
C
C           1 X 1 PIVOT BLOCK.
C
            IF (.NOT.SWAP) GO TO 120
C
C              PERFORM AN INTERCHANGE.
C
               CALL CSWAP(IMAX,AP(IM+1),1,AP(IK+1),1)
               IMJ = IK + IMAX
               DO 110 JJ = IMAX, K
                  J = K + IMAX - JJ
                  JK = IK + J
                  T = CONJG(AP(JK))
                  AP(JK) = CONJG(AP(IMJ))
                  AP(IMJ) = T
                  IMJ = IMJ - (J - 1)
  110          CONTINUE
  120       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            IJ = IK - (K - 1)
            DO 130 JJ = 1, KM1
               J = K - JJ
               JK = IK + J
               MULK = -AP(JK)/AP(KK)
               T = CONJG(MULK)
               CALL CAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
               IJJ = IJ + J
               AP(IJJ) = CMPLX(REAL(AP(IJJ)),0.0E0)
               AP(JK) = MULK
               IJ = IJ - (J - 1)
  130       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = K
            IF (SWAP) KPVT(K) = IMAX
         GO TO 190
  140    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            KM1K = IK + K - 1
            IKM1 = IK - (K - 1)
            IF (.NOT.SWAP) GO TO 160
C
C              PERFORM AN INTERCHANGE.
C
               CALL CSWAP(IMAX,AP(IM+1),1,AP(IKM1+1),1)
               IMJ = IKM1 + IMAX
               DO 150 JJ = IMAX, KM1
                  J = KM1 + IMAX - JJ
                  JKM1 = IKM1 + J
                  T = CONJG(AP(JKM1))
                  AP(JKM1) = CONJG(AP(IMJ))
                  AP(IMJ) = T
                  IMJ = IMJ - (J - 1)
  150          CONTINUE
               T = AP(KM1K)
               AP(KM1K) = AP(IMK)
               AP(IMK) = T
  160       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            KM2 = K - 2
            IF (KM2 .EQ. 0) GO TO 180
               AK = AP(KK)/AP(KM1K)
               KM1KM1 = IKM1 + K - 1
               AKM1 = AP(KM1KM1)/CONJG(AP(KM1K))
               DENOM = 1.0E0 - AK*AKM1
               IJ = IK - (K - 1) - (K - 2)
               DO 170 JJ = 1, KM2
                  J = KM1 - JJ
                  JK = IK + J
                  BK = AP(JK)/AP(KM1K)
                  JKM1 = IKM1 + J
                  BKM1 = AP(JKM1)/CONJG(AP(KM1K))
                  MULK = (AKM1*BK - BKM1)/DENOM
                  MULKM1 = (AK*BKM1 - BK)/DENOM
                  T = CONJG(MULK)
                  CALL CAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
                  T = CONJG(MULKM1)
                  CALL CAXPY(J,T,AP(IKM1+1),1,AP(IJ+1),1)
                  AP(JK) = MULK
                  AP(JKM1) = MULKM1
                  IJJ = IJ + J
                  AP(IJJ) = CMPLX(REAL(AP(IJJ)),0.0E0)
                  IJ = IJ - (J - 1)
  170          CONTINUE
  180       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = 1 - K
            IF (SWAP) KPVT(K) = -IMAX
            KPVT(K-1) = KPVT(K)
  190    CONTINUE
         IK = IK - (K - 1)
         IF (KSTEP .EQ. 2) IK = IK - (K - 2)
         K = K - KSTEP
      GO TO 10
  200 CONTINUE
      RETURN
      END
      SUBROUTINE CHPSL(AP,N,KPVT,B)
      INTEGER N,KPVT(1)
      COMPLEX AP(1),B(1)
C
C     CHISL SOLVES THE COMPLEX HERMITIAN SYSTEM
C     A * X = B
C     USING THE FACTORS COMPUTED BY CHPFA.
C
C     ON ENTRY
C
C        AP      COMPLEX(N*(N+1)/2)
C                THE OUTPUT FROM CHPFA.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        KPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM CHPFA.
C
C        B       COMPLEX(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO MAY OCCUR IF  CHPCO  HAS SET RCOND .EQ. 0.0
C        OR  CHPFA  HAS SET INFO .NE. 0  .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL CHPFA(AP,N,KPVT,INFO)
C           IF (INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL CHPSL(AP,N,KPVT,C(1,J))
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CDOTC
C     FORTRAN CONJG,IABS
C
C     INTERNAL VARIABLES.
C
      COMPLEX AK,AKM1,BK,BKM1,CDOTC,DENOM,TEMP
      INTEGER IK,IKM1,IKP1,K,KK,KM1K,KM1KM1,KP
C
C     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
C     D INVERSE TO B.
C
      K = N
      IK = (N*(N - 1))/2
   10 IF (K .EQ. 0) GO TO 80
         KK = IK + K
         IF (KPVT(K) .LT. 0) GO TO 40
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 30
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 20
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
   20          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL CAXPY(K-1,B(K),AP(IK+1),1,B(1),1)
   30       CONTINUE
C
C           APPLY D INVERSE.
C
            B(K) = B(K)/AP(KK)
            K = K - 1
            IK = IK - K
         GO TO 70
   40    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IKM1 = IK - (K - 1)
            IF (K .EQ. 2) GO TO 60
               KP = IABS(KPVT(K))
               IF (KP .EQ. K - 1) GO TO 50
C
C                 INTERCHANGE.
C
                  TEMP = B(K-1)
                  B(K-1) = B(KP)
                  B(KP) = TEMP
   50          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL CAXPY(K-2,B(K),AP(IK+1),1,B(1),1)
               CALL CAXPY(K-2,B(K-1),AP(IKM1+1),1,B(1),1)
   60       CONTINUE
C
C           APPLY D INVERSE.
C
            KM1K = IK + K - 1
            KK = IK + K
            AK = AP(KK)/CONJG(AP(KM1K))
            KM1KM1 = IKM1 + K - 1
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = B(K)/CONJG(AP(KM1K))
            BKM1 = B(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0E0
            B(K) = (AKM1*BK - BKM1)/DENOM
            B(K-1) = (AK*BKM1 - BK)/DENOM
            K = K - 2
            IK = IK - (K + 1) - K
   70    CONTINUE
      GO TO 10
   80 CONTINUE
C
C     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
C
      K = 1
      IK = 0
   90 IF (K .GT. N) GO TO 160
         IF (KPVT(K) .LT. 0) GO TO 120
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 110
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + CDOTC(K-1,AP(IK+1),1,B(1),1)
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 100
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  100          CONTINUE
  110       CONTINUE
            IK = IK + K
            K = K + 1
         GO TO 150
  120    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 140
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + CDOTC(K-1,AP(IK+1),1,B(1),1)
               IKP1 = IK + K
               B(K+1) = B(K+1) + CDOTC(K-1,AP(IKP1+1),1,B(1),1)
               KP = IABS(KPVT(K))
               IF (KP .EQ. K) GO TO 130
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  130          CONTINUE
  140       CONTINUE
            IK = IK + K + K + 1
            K = K + 2
  150    CONTINUE
      GO TO 90
  160 CONTINUE
      RETURN
      END
      SUBROUTINE CHPDI(AP,N,KPVT,DET,INERT,WORK,JOB)
      INTEGER N,JOB
      COMPLEX AP(1),WORK(1)
      REAL DET(2)
      INTEGER KPVT(1),INERT(3)
C
C     CHPDI COMPUTES THE DETERMINANT, INERTIA AND INVERSE
C     OF A COMPLEX HERMITIAN MATRIX USING THE FACTORS FROM CHPFA,
C     WHERE THE MATRIX IS STORED IN PACKED FORM.
C
C     ON ENTRY
C
C        AP      COMPLEX (N*(N+1)/2)
C                THE OUTPUT FROM CHPFA.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX A.
C
C        KPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM CHPFA.
C
C        WORK    COMPLEX(N)
C                WORK VECTOR.  CONTENTS IGNORED.
C
C        JOB     INTEGER
C                JOB HAS THE DECIMAL EXPANSION  ABC  WHERE
C                   IF  C .NE. 0, THE INVERSE IS COMPUTED,
C                   IF  B .NE. 0, THE DETERMINANT IS COMPUTED,
C                   IF  A .NE. 0, THE INERTIA IS COMPUTED.
C
C                FOR EXAMPLE, JOB = 111  GIVES ALL THREE.
C
C     ON RETURN
C
C        VARIABLES NOT REQUESTED BY JOB ARE NOT USED.
C
C        AP     CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF
C               THE ORIGINAL MATRIX, STORED IN PACKED FORM.
C               THE COLUMNS OF THE UPPER TRIANGLE ARE STORED
C               SEQUENTIALLY IN A ONE-DIMENSIONAL ARRAY.
C
C        DET    REAL(2)
C               DETERMINANT OF ORIGINAL MATRIX.
C               DETERMINANT = DET(1) * 10.0**DET(2)
C               WITH 1.0 .LE. ABS(DET(1)) .LT. 10.0
C               OR DET(1) = 0.0.
C
C        INERT  INTEGER(3)
C               THE INERTIA OF THE ORIGINAL MATRIX.
C               INERT(1)  =  NUMBER OF POSITIVE EIGENVALUES.
C               INERT(2)  =  NUMBER OF NEGATIVE EIGENVALUES.
C               INERT(3)  =  NUMBER OF ZERO EIGENVALUES.
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INVERSE IS REQUESTED
C        AND  CHPCO  HAS SET RCOND .EQ. 0.0
C        OR  CHPFA  HAS SET  INFO .NE. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CCOPY,CDOTC,CSWAP
C     FORTRAN ABS,CABS,CMPLX,CONJG,IABS,MOD,REAL
C
C     INTERNAL VARIABLES.
C
      COMPLEX AKKP1,CDOTC,TEMP
      REAL TEN,D,T,AK,AKP1
      INTEGER IJ,IK,IKP1,IKS,J,JB,JK,JKP1
      INTEGER K,KK,KKP1,KM1,KS,KSJ,KSKP1,KSTEP
      LOGICAL NOINV,NODET,NOERT
C
      NOINV = MOD(JOB,10) .EQ. 0
      NODET = MOD(JOB,100)/10 .EQ. 0
      NOERT = MOD(JOB,1000)/100 .EQ. 0
C
      IF (NODET .AND. NOERT) GO TO 140
         IF (NOERT) GO TO 10
            INERT(1) = 0
            INERT(2) = 0
            INERT(3) = 0
   10    CONTINUE
         IF (NODET) GO TO 20
            DET(1) = 1.0E0
            DET(2) = 0.0E0
            TEN = 10.0E0
   20    CONTINUE
         T = 0.0E0
         IK = 0
         DO 130 K = 1, N
            KK = IK + K
            D = REAL(AP(KK))
C
C           CHECK IF 1 BY 1
C
            IF (KPVT(K) .GT. 0) GO TO 50
C
C              2 BY 2 BLOCK
C              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = CABS(S)
C                      (S  C)
C              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
C              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
C
               IF (T .NE. 0.0E0) GO TO 30
                  IKP1 = IK + K
                  KKP1 = IKP1 + K
                  T = CABS(AP(KKP1))
                  D = (D/T)*REAL(AP(KKP1+1)) - T
               GO TO 40
   30          CONTINUE
                  D = T
                  T = 0.0E0
   40          CONTINUE
   50       CONTINUE
C
            IF (NOERT) GO TO 60
               IF (D .GT. 0.0E0) INERT(1) = INERT(1) + 1
               IF (D .LT. 0.0E0) INERT(2) = INERT(2) + 1
               IF (D .EQ. 0.0E0) INERT(3) = INERT(3) + 1
   60       CONTINUE
C
            IF (NODET) GO TO 120
               DET(1) = D*DET(1)
               IF (DET(1) .EQ. 0.0E0) GO TO 110
   70             IF (ABS(DET(1)) .GE. 1.0E0) GO TO 80
                     DET(1) = TEN*DET(1)
                     DET(2) = DET(2) - 1.0E0
                  GO TO 70
   80             CONTINUE
   90             IF (ABS(DET(1)) .LT. TEN) GO TO 100
                     DET(1) = DET(1)/TEN
                     DET(2) = DET(2) + 1.0E0
                  GO TO 90
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
            IK = IK + K
  130    CONTINUE
  140 CONTINUE
C
C     COMPUTE INVERSE(A)
C
      IF (NOINV) GO TO 270
         K = 1
         IK = 0
  150    IF (K .GT. N) GO TO 260
            KM1 = K - 1
            KK = IK + K
            IKP1 = IK + K
            KKP1 = IKP1 + K
            IF (KPVT(K) .LT. 0) GO TO 180
C
C              1 BY 1
C
               AP(KK) = CMPLX(1.0E0/REAL(AP(KK)),0.0E0)
               IF (KM1 .LT. 1) GO TO 170
                  CALL CCOPY(KM1,AP(IK+1),1,WORK,1)
                  IJ = 0
                  DO 160 J = 1, KM1
                     JK = IK + J
                     AP(JK) = CDOTC(J,AP(IJ+1),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
                     IJ = IJ + J
  160             CONTINUE
                  AP(KK) = AP(KK)
     *                     + CMPLX(REAL(CDOTC(KM1,WORK,1,AP(IK+1),1)),
     *                             0.0E0)
  170          CONTINUE
               KSTEP = 1
            GO TO 220
  180       CONTINUE
C
C              2 BY 2
C
               T = CABS(AP(KKP1))
               AK = REAL(AP(KK))/T
               AKP1 = REAL(AP(KKP1+1))/T
               AKKP1 = AP(KKP1)/T
               D = T*(AK*AKP1 - 1.0E0)
               AP(KK) = CMPLX(AKP1/D,0.0E0)
               AP(KKP1+1) = CMPLX(AK/D,0.0E0)
               AP(KKP1) = -AKKP1/D
               IF (KM1 .LT. 1) GO TO 210
                  CALL CCOPY(KM1,AP(IKP1+1),1,WORK,1)
                  IJ = 0
                  DO 190 J = 1, KM1
                     JKP1 = IKP1 + J
                     AP(JKP1) = CDOTC(J,AP(IJ+1),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IKP1+1),1)
                     IJ = IJ + J
  190             CONTINUE
                  AP(KKP1+1) = AP(KKP1+1)
     *                         + CMPLX(REAL(CDOTC(KM1,WORK,1,
     *                                            AP(IKP1+1),1)),0.0E0)
                  AP(KKP1) = AP(KKP1)
     *                       + CDOTC(KM1,AP(IK+1),1,AP(IKP1+1),1)
                  CALL CCOPY(KM1,AP(IK+1),1,WORK,1)
                  IJ = 0
                  DO 200 J = 1, KM1
                     JK = IK + J
                     AP(JK) = CDOTC(J,AP(IJ+1),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
                     IJ = IJ + J
  200             CONTINUE
                  AP(KK) = AP(KK)
     *                     + CMPLX(REAL(CDOTC(KM1,WORK,1,AP(IK+1),1)),
     *                             0.0E0)
  210          CONTINUE
               KSTEP = 2
  220       CONTINUE
C
C           SWAP
C
            KS = IABS(KPVT(K))
            IF (KS .EQ. K) GO TO 250
               IKS = (KS*(KS - 1))/2
               CALL CSWAP(KS,AP(IKS+1),1,AP(IK+1),1)
               KSJ = IK + KS
               DO 230 JB = KS, K
                  J = K + KS - JB
                  JK = IK + J
                  TEMP = CONJG(AP(JK))
                  AP(JK) = CONJG(AP(KSJ))
                  AP(KSJ) = TEMP
                  KSJ = KSJ - (J - 1)
  230          CONTINUE
               IF (KSTEP .EQ. 1) GO TO 240
                  KSKP1 = IKP1 + KS
                  TEMP = AP(KSKP1)
                  AP(KSKP1) = AP(KKP1)
                  AP(KKP1) = TEMP
  240          CONTINUE
  250       CONTINUE
            IK = IK + K
            IF (KSTEP .EQ. 2) IK = IK + K + 1
            K = K + KSTEP
         GO TO 150
  260    CONTINUE
  270 CONTINUE
      RETURN
      END
      SUBROUTINE CTRCO(T,LDT,N,RCOND,Z,JOB)
      INTEGER LDT,N,JOB
      COMPLEX T(LDT,1),Z(1)
      REAL RCOND
C
C     CTRCO ESTIMATES THE CONDITION OF A COMPLEX TRIANGULAR MATRIX.
C
C     ON ENTRY
C
C        T       COMPLEX(LDT,N)
C                T CONTAINS THE TRIANGULAR MATRIX. THE ZERO
C                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND
C                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE
C                USED TO STORE OTHER INFORMATION.
C
C        LDT     INTEGER
C                LDT IS THE LEADING DIMENSION OF THE ARRAY T.
C
C        N       INTEGER
C                N IS THE ORDER OF THE SYSTEM.
C
C        JOB     INTEGER
C                = 0         T  IS LOWER TRIANGULAR.
C                = NONZERO   T  IS UPPER TRIANGULAR.
C
C     ON RETURN
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  T .
C                FOR THE SYSTEM  T*X = B , RELATIVE PERTURBATIONS
C                IN  T  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  T  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       COMPLEX(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  T  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CSSCAL,SCASUM
C     FORTRAN ABS,AIMAG,AMAX1,CMPLX,CONJG,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX W,WK,WKM,EK
      REAL TNORM,YNORM,S,SM,SCASUM
      INTEGER I1,J,J1,J2,K,KK,L
      LOGICAL LOWER
      COMPLEX ZDUM,ZDUM1,ZDUM2,CSIGN1
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
      CSIGN1(ZDUM1,ZDUM2) = CABS1(ZDUM1)*(ZDUM2/CABS1(ZDUM2))
C
      LOWER = JOB .EQ. 0
C
C     COMPUTE 1-NORM OF T
C
      TNORM = 0.0E0
      DO 10 J = 1, N
         L = J
         IF (LOWER) L = N + 1 - J
         I1 = 1
         IF (LOWER) I1 = J
         TNORM = AMAX1(TNORM,SCASUM(L,T(I1,J),1))
   10 CONTINUE
C
C     RCOND = 1/(NORM(T)*(ESTIMATE OF NORM(INVERSE(T)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  T*Z = Y  AND  CTRANS(T)*Y = E .
C     CTRANS(T)  IS THE CONJUGATE TRANSPOSE OF T .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF Y .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C     SOLVE CTRANS(T)*Y = E
C
      EK = (1.0E0,0.0E0)
      DO 20 J = 1, N
         Z(J) = (0.0E0,0.0E0)
   20 CONTINUE
      DO 100 KK = 1, N
         K = KK
         IF (LOWER) K = N + 1 - KK
         IF (CABS1(Z(K)) .NE. 0.0E0) EK = CSIGN1(EK,-Z(K))
         IF (CABS1(EK-Z(K)) .LE. CABS1(T(K,K))) GO TO 30
            S = CABS1(T(K,K))/CABS1(EK-Z(K))
            CALL CSSCAL(N,S,Z,1)
            EK = CMPLX(S,0.0E0)*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = CABS1(WK)
         SM = CABS1(WKM)
         IF (CABS1(T(K,K)) .EQ. 0.0E0) GO TO 40
            WK = WK/CONJG(T(K,K))
            WKM = WKM/CONJG(T(K,K))
         GO TO 50
   40    CONTINUE
            WK = (1.0E0,0.0E0)
            WKM = (1.0E0,0.0E0)
   50    CONTINUE
         IF (KK .EQ. N) GO TO 90
            J1 = K + 1
            IF (LOWER) J1 = 1
            J2 = N
            IF (LOWER) J2 = K - 1
            DO 60 J = J1, J2
               SM = SM + CABS1(Z(J)+WKM*CONJG(T(K,J)))
               Z(J) = Z(J) + WK*CONJG(T(K,J))
               S = S + CABS1(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               W = WKM - WK
               WK = WKM
               DO 70 J = J1, J2
                  Z(J) = Z(J) + W*CONJG(T(K,J))
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
      YNORM = 1.0E0
C
C     SOLVE T*Z = Y
C
      DO 130 KK = 1, N
         K = N + 1 - KK
         IF (LOWER) K = KK
         IF (CABS1(Z(K)) .LE. CABS1(T(K,K))) GO TO 110
            S = CABS1(T(K,K))/CABS1(Z(K))
            CALL CSSCAL(N,S,Z,1)
            YNORM = S*YNORM
  110    CONTINUE
         IF (CABS1(T(K,K)) .NE. 0.0E0) Z(K) = Z(K)/T(K,K)
         IF (CABS1(T(K,K)) .EQ. 0.0E0) Z(K) = (1.0E0,0.0E0)
         I1 = 1
         IF (LOWER) I1 = K + 1
         IF (KK .GE. N) GO TO 120
            W = -Z(K)
            CALL CAXPY(N-KK,W,T(I1,K),1,Z(I1),1)
  120    CONTINUE
  130 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (TNORM .NE. 0.0E0) RCOND = YNORM/TNORM
      IF (TNORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
      SUBROUTINE CTRSL(T,LDT,N,B,JOB,INFO)
      INTEGER LDT,N,JOB,INFO
      COMPLEX T(LDT,1),B(1)
C
C
C     CTRSL SOLVES SYSTEMS OF THE FORM
C
C                   T * X = B
C     OR
C                   CTRANS(T) * X = B
C
C     WHERE T IS A TRIANGULAR MATRIX OF ORDER N. HERE CTRANS(T)
C     DENOTES THE CONJUGATE TRANSPOSE OF THE MATRIX T.
C
C     ON ENTRY
C
C         T         COMPLEX(LDT,N)
C                   T CONTAINS THE MATRIX OF THE SYSTEM. THE ZERO
C                   ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND
C                   THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE
C                   USED TO STORE OTHER INFORMATION.
C
C         LDT       INTEGER
C                   LDT IS THE LEADING DIMENSION OF THE ARRAY T.
C
C         N         INTEGER
C                   N IS THE ORDER OF THE SYSTEM.
C
C         B         COMPLEX(N).
C                   B CONTAINS THE RIGHT HAND SIDE OF THE SYSTEM.
C
C         JOB       INTEGER
C                   JOB SPECIFIES WHAT KIND OF SYSTEM IS TO BE SOLVED.
C                   IF JOB IS
C
C                        00   SOLVE T*X=B, T LOWER TRIANGULAR,
C                        01   SOLVE T*X=B, T UPPER TRIANGULAR,
C                        10   SOLVE CTRANS(T)*X=B, T LOWER TRIANGULAR,
C                        11   SOLVE CTRANS(T)*X=B, T UPPER TRIANGULAR.
C
C     ON RETURN
C
C         B         B CONTAINS THE SOLUTION, IF INFO .EQ. 0.
C                   OTHERWISE B IS UNALTERED.
C
C         INFO      INTEGER
C                   INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR.
C                   OTHERWISE INFO CONTAINS THE INDEX OF
C                   THE FIRST ZERO DIAGONAL ELEMENT OF T.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     G. W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CDOTC
C     FORTRAN ABS,AIMAG,CONJG,MOD,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX CDOTC,TEMP
      INTEGER CASE,J,JJ
C
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C     BEGIN BLOCK PERMITTING ...EXITS TO 150
C
C        CHECK FOR ZERO DIAGONAL ELEMENTS.
C
         DO 10 INFO = 1, N
C     ......EXIT
            IF (CABS1(T(INFO,INFO)) .EQ. 0.0E0) GO TO 150
   10    CONTINUE
         INFO = 0
C
C        DETERMINE THE TASK AND GO TO IT.
C
         CASE = 1
         IF (MOD(JOB,10) .NE. 0) CASE = 2
         IF (MOD(JOB,100)/10 .NE. 0) CASE = CASE + 2
         GO TO (20,50,80,110), CASE
C
C        SOLVE T*X=B FOR T LOWER TRIANGULAR
C
   20    CONTINUE
            B(1) = B(1)/T(1,1)
            IF (N .LT. 2) GO TO 40
            DO 30 J = 2, N
               TEMP = -B(J-1)
               CALL CAXPY(N-J+1,TEMP,T(J,J-1),1,B(J),1)
               B(J) = B(J)/T(J,J)
   30       CONTINUE
   40       CONTINUE
         GO TO 140
C
C        SOLVE T*X=B FOR T UPPER TRIANGULAR.
C
   50    CONTINUE
            B(N) = B(N)/T(N,N)
            IF (N .LT. 2) GO TO 70
            DO 60 JJ = 2, N
               J = N - JJ + 1
               TEMP = -B(J+1)
               CALL CAXPY(J,TEMP,T(1,J+1),1,B(1),1)
               B(J) = B(J)/T(J,J)
   60       CONTINUE
   70       CONTINUE
         GO TO 140
C
C        SOLVE CTRANS(T)*X=B FOR T LOWER TRIANGULAR.
C
   80    CONTINUE
            B(N) = B(N)/CONJG(T(N,N))
            IF (N .LT. 2) GO TO 100
            DO 90 JJ = 2, N
               J = N - JJ + 1
               B(J) = B(J) - CDOTC(JJ-1,T(J+1,J),1,B(J+1),1)
               B(J) = B(J)/CONJG(T(J,J))
   90       CONTINUE
  100       CONTINUE
         GO TO 140
C
C        SOLVE CTRANS(T)*X=B FOR T UPPER TRIANGULAR.
C
  110    CONTINUE
            B(1) = B(1)/CONJG(T(1,1))
            IF (N .LT. 2) GO TO 130
            DO 120 J = 2, N
               B(J) = B(J) - CDOTC(J-1,T(1,J),1,B(1),1)
               B(J) = B(J)/CONJG(T(J,J))
  120       CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
      SUBROUTINE CTRDI(T,LDT,N,DET,JOB,INFO)
      INTEGER LDT,N,JOB,INFO
      COMPLEX T(LDT,1),DET(2)
C
C     CTRDI COMPUTES THE DETERMINANT AND INVERSE OF A COMPLEX
C     TRIANGULAR MATRIX.
C
C     ON ENTRY
C
C        T       COMPLEX(LDT,N)
C                T CONTAINS THE TRIANGULAR MATRIX. THE ZERO
C                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND
C                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE
C                USED TO STORE OTHER INFORMATION.
C
C        LDT     INTEGER
C                LDT IS THE LEADING DIMENSION OF THE ARRAY T.
C
C        N       INTEGER
C                N IS THE ORDER OF THE SYSTEM.
C
C        JOB     INTEGER
C                = 010       NO DET, INVERSE OF LOWER TRIANGULAR.
C                = 011       NO DET, INVERSE OF UPPER TRIANGULAR.
C                = 100       DET, NO INVERSE.
C                = 110       DET, INVERSE OF LOWER TRIANGULAR.
C                = 111       DET, INVERSE OF UPPER TRIANGULAR.
C
C     ON RETURN
C
C        T       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE UNCHANGED.
C
C        DET     COMPLEX(2)
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE NOT REFERENCED.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. CABS1(DET(1)) .LT. 10.0
C                OR  DET(1) .EQ. 0.0 .
C
C        INFO    INTEGER
C                INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR
C                AND THE INVERSE IS REQUESTED.
C                OTHERWISE INFO CONTAINS THE INDEX OF
C                A ZERO DIAGONAL ELEMENT OF T.
C
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CSCAL
C     FORTRAN ABS,AIMAG,CMPLX,MOD,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX TEMP
      REAL TEN
      INTEGER I,J,K,KB,KM1,KP1
C
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C     BEGIN BLOCK PERMITTING ...EXITS TO 180
C
C        COMPUTE DETERMINANT
C
         IF (JOB/100 .EQ. 0) GO TO 70
            DET(1) = (1.0E0,0.0E0)
            DET(2) = (0.0E0,0.0E0)
            TEN = 10.0E0
            DO 50 I = 1, N
               DET(1) = T(I,I)*DET(1)
C           ...EXIT
               IF (CABS1(DET(1)) .EQ. 0.0E0) GO TO 60
   10          IF (CABS1(DET(1)) .GE. 1.0E0) GO TO 20
                  DET(1) = CMPLX(TEN,0.0E0)*DET(1)
                  DET(2) = DET(2) - (1.0E0,0.0E0)
               GO TO 10
   20          CONTINUE
   30          IF (CABS1(DET(1)) .LT. TEN) GO TO 40
                  DET(1) = DET(1)/CMPLX(TEN,0.0E0)
                  DET(2) = DET(2) + (1.0E0,0.0E0)
               GO TO 30
   40          CONTINUE
   50       CONTINUE
   60       CONTINUE
   70    CONTINUE
C
C        COMPUTE INVERSE OF UPPER TRIANGULAR
C
         IF (MOD(JOB/10,10) .EQ. 0) GO TO 170
            IF (MOD(JOB,10) .EQ. 0) GO TO 120
C              BEGIN BLOCK PERMITTING ...EXITS TO 110
                  DO 100 K = 1, N
                     INFO = K
C              ......EXIT
                     IF (CABS1(T(K,K)) .EQ. 0.0E0) GO TO 110
                     T(K,K) = (1.0E0,0.0E0)/T(K,K)
                     TEMP = -T(K,K)
                     CALL CSCAL(K-1,TEMP,T(1,K),1)
                     KP1 = K + 1
                     IF (N .LT. KP1) GO TO 90
                     DO 80 J = KP1, N
                        TEMP = T(K,J)
                        T(K,J) = (0.0E0,0.0E0)
                        CALL CAXPY(K,TEMP,T(1,K),1,T(1,J),1)
   80                CONTINUE
   90                CONTINUE
  100             CONTINUE
                  INFO = 0
  110          CONTINUE
            GO TO 160
  120       CONTINUE
C
C              COMPUTE INVERSE OF LOWER TRIANGULAR
C
               DO 150 KB = 1, N
                  K = N + 1 - KB
                  INFO = K
C     ............EXIT
                  IF (CABS1(T(K,K)) .EQ. 0.0E0) GO TO 180
                  T(K,K) = (1.0E0,0.0E0)/T(K,K)
                  TEMP = -T(K,K)
                  IF (K .NE. N) CALL CSCAL(N-K,TEMP,T(K+1,K),1)
                  KM1 = K - 1
                  IF (KM1 .LT. 1) GO TO 140
                  DO 130 J = 1, KM1
                     TEMP = T(K,J)
                     T(K,J) = (0.0E0,0.0E0)
                     CALL CAXPY(N-K+1,TEMP,T(K,K),1,T(K,J),1)
  130             CONTINUE
  140             CONTINUE
  150          CONTINUE
               INFO = 0
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE
      RETURN
      END
      SUBROUTINE CGTSL(N,C,D,E,B,INFO)
      INTEGER N,INFO
      COMPLEX C(1),D(1),E(1),B(1)
C
C     CGTSL GIVEN A GENERAL TRIDIAGONAL MATRIX AND A RIGHT HAND
C     SIDE WILL FIND THE SOLUTION.
C
C     ON ENTRY
C
C        N       INTEGER
C                IS THE ORDER OF THE TRIDIAGONAL MATRIX.
C
C        C       COMPLEX(N)
C                IS THE SUBDIAGONAL OF THE TRIDIAGONAL MATRIX.
C                C(2) THROUGH C(N) SHOULD CONTAIN THE SUBDIAGONAL.
C                ON OUTPUT C IS DESTROYED.
C
C        D       COMPLEX(N)
C                IS THE DIAGONAL OF THE TRIDIAGONAL MATRIX.
C                ON OUTPUT D IS DESTROYED.
C
C        E       COMPLEX(N)
C                IS THE SUPERDIAGONAL OF THE TRIDIAGONAL MATRIX.
C                E(1) THROUGH E(N-1) SHOULD CONTAIN THE SUPERDIAGONAL.
C                ON OUTPUT E IS DESTROYED.
C
C        B       COMPLEX(N)
C                IS THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       IS THE SOLUTION VECTOR.
C
C        INFO    INTEGER
C                = 0 NORMAL VALUE.
C                = K IF THE K-TH ELEMENT OF THE DIAGONAL BECOMES
C                    EXACTLY ZERO.  THE SUBROUTINE RETURNS WHEN
C                    THIS IS DETECTED.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.
C
C     NO EXTERNALS
C     FORTRAN ABS,AIMAG,REAL
C
C     INTERNAL VARIABLES
C
      INTEGER K,KB,KP1,NM1,NM2
      COMPLEX T
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C     BEGIN BLOCK PERMITTING ...EXITS TO 100
C
         INFO = 0
         C(1) = D(1)
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 40
            D(1) = E(1)
            E(1) = (0.0E0,0.0E0)
            E(N) = (0.0E0,0.0E0)
C
            DO 30 K = 1, NM1
               KP1 = K + 1
C
C              FIND THE LARGEST OF THE TWO ROWS
C
               IF (CABS1(C(KP1)) .LT. CABS1(C(K))) GO TO 10
C
C                 INTERCHANGE ROW
C
                  T = C(KP1)
                  C(KP1) = C(K)
                  C(K) = T
                  T = D(KP1)
                  D(KP1) = D(K)
                  D(K) = T
                  T = E(KP1)
                  E(KP1) = E(K)
                  E(K) = T
                  T = B(KP1)
                  B(KP1) = B(K)
                  B(K) = T
   10          CONTINUE
C
C              ZERO ELEMENTS
C
               IF (CABS1(C(K)) .NE. 0.0E0) GO TO 20
                  INFO = K
C     ............EXIT
                  GO TO 100
   20          CONTINUE
               T = -C(KP1)/C(K)
               C(KP1) = D(KP1) + T*D(K)
               D(KP1) = E(KP1) + T*E(K)
               E(KP1) = (0.0E0,0.0E0)
               B(KP1) = B(KP1) + T*B(K)
   30       CONTINUE
   40    CONTINUE
         IF (CABS1(C(N)) .NE. 0.0E0) GO TO 50
            INFO = N
         GO TO 90
   50    CONTINUE
C
C           BACK SOLVE
C
            NM2 = N - 2
            B(N) = B(N)/C(N)
            IF (N .EQ. 1) GO TO 80
               B(NM1) = (B(NM1) - D(NM1)*B(N))/C(NM1)
               IF (NM2 .LT. 1) GO TO 70
               DO 60 KB = 1, NM2
                  K = NM2 - KB + 1
                  B(K) = (B(K) - D(K)*B(K+1) - E(K)*B(K+2))/C(K)
   60          CONTINUE
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE CPTSL(N,D,E,B)
      INTEGER N
      COMPLEX D(1),E(1),B(1)
C
C     CPTSL GIVEN A POSITIVE DEFINITE TRIDIAGONAL MATRIX AND A RIGHT
C     HAND SIDE WILL FIND THE SOLUTION.
C
C     ON ENTRY
C
C        N        INTEGER
C                 IS THE ORDER OF THE TRIDIAGONAL MATRIX.
C
C        D        COMPLEX(N)
C                 IS THE DIAGONAL OF THE TRIDIAGONAL MATRIX.
C                 ON OUTPUT D IS DESTROYED.
C
C        E        COMPLEX(N)
C                 IS THE OFFDIAGONAL OF THE TRIDIAGONAL MATRIX.
C                 E(1) THROUGH E(N-1) SHOULD CONTAIN THE
C                 OFFDIAGONAL.
C
C        B        COMPLEX(N)
C                 IS THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B        CONTAINS THE SOULTION.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.
C
C     NO EXTERNALS
C     FORTRAN CONJG,MOD
C
C     INTERNAL VARIABLES
C
      INTEGER K,KBM1,KE,KF,KP1,NM1,NM1D2
      COMPLEX T1,T2
C
C     CHECK FOR 1 X 1 CASE
C
      IF (N .NE. 1) GO TO 10
         B(1) = B(1)/D(1)
      GO TO 70
   10 CONTINUE
         NM1 = N - 1
         NM1D2 = NM1/2
         IF (N .EQ. 2) GO TO 30
            KBM1 = N - 1
C
C           ZERO TOP HALF OF SUBDIAGONAL AND BOTTOM HALF OF
C           SUPERDIAGONAL
C
            DO 20 K = 1, NM1D2
               T1 = CONJG(E(K))/D(K)
               D(K+1) = D(K+1) - T1*E(K)
               B(K+1) = B(K+1) - T1*B(K)
               T2 = E(KBM1)/D(KBM1+1)
               D(KBM1) = D(KBM1) - T2*CONJG(E(KBM1))
               B(KBM1) = B(KBM1) - T2*B(KBM1+1)
               KBM1 = KBM1 - 1
   20       CONTINUE
   30    CONTINUE
         KP1 = NM1D2 + 1
C
C        CLEAN UP FOR POSSIBLE 2 X 2 BLOCK AT CENTER
C
         IF (MOD(N,2) .NE. 0) GO TO 40
            T1 = CONJG(E(KP1))/D(KP1)
            D(KP1+1) = D(KP1+1) - T1*E(KP1)
            B(KP1+1) = B(KP1+1) - T1*B(KP1)
            KP1 = KP1 + 1
   40    CONTINUE
C
C        BACK SOLVE STARTING AT THE CENTER, GOING TOWARDS THE TOP
C        AND BOTTOM
C
         B(KP1) = B(KP1)/D(KP1)
         IF (N .EQ. 2) GO TO 60
            K = KP1 - 1
            KE = KP1 + NM1D2 - 1
            DO 50 KF = KP1, KE
               B(K) = (B(K) - E(K)*B(K+1))/D(K)
               B(KF+1) = (B(KF+1) - CONJG(E(KF))*B(KF))/D(KF+1)
               K = K - 1
   50       CONTINUE
   60    CONTINUE
         IF (MOD(N,2) .EQ. 0) B(1) = (B(1) - E(1)*B(2))/D(1)
   70 CONTINUE
      RETURN
      END
      SUBROUTINE CCHDC(A,LDA,P,WORK,JPVT,JOB,INFO)
      INTEGER LDA,P,JPVT(1),JOB,INFO
      COMPLEX A(LDA,1),WORK(1)
C
C     CCHDC COMPUTES THE CHOLESKY DECOMPOSITION OF A POSITIVE DEFINITE
C     MATRIX.  A PIVOTING OPTION ALLOWS THE USER TO ESTIMATE THE
C     CONDITION OF A POSITIVE DEFINITE MATRIX OR DETERMINE THE RANK
C     OF A POSITIVE SEMIDEFINITE MATRIX.
C
C     ON ENTRY
C
C         A      COMPLEX(LDA,P).
C                A CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO
C                BE COMPUTED.  ONLT THE UPPER HALF OF A NEED BE STORED.
C                THE LOWER PART OF THE ARRAY A IS NOT REFERENCED.
C
C         LDA    INTEGER.
C                LDA IS THE LEADING DIMENSION OF THE ARRAY A.
C
C         P      INTEGER.
C                P IS THE ORDER OF THE MATRIX.
C
C         WORK   COMPLEX.
C                WORK IS A WORK ARRAY.
C
C         JPVT   INTEGER(P).
C                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION
C                OF THE PIVOT ELEMENTS, IF PIVOTING HAS BEEN REQUESTED.
C                EACH DIAGONAL ELEMENT A(K,K)
C                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE
C                VALUE OF JPVT(K).
C
C                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL
C                                      ELEMENT.
C
C                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE ELEMENT.
C
C                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL ELEMENT.
C
C                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL ELEMENTS
C                ARE MOVED BY SYMMETRIC ROW AND COLUMN INTERCHANGES TO
C                THE BEGINNING OF THE ARRAY A AND FINAL
C                ELEMENTS TO THE END.  BOTH INITIAL AND FINAL ELEMENTS
C                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY
C                FREE ELEMENTS ARE MOVED.  AT THE K-TH STAGE OF THE
C                REDUCTION, IF A(K,K) IS OCCUPIED BY A FREE ELEMENT
C                IT IS INTERCHANGED WITH THE LARGEST FREE ELEMENT
C                A(L,L) WITH L .GE. K.  JPVT IS NOT REFERENCED IF
C                JOB .EQ. 0.
C
C        JOB     INTEGER.
C                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING.
C                IF JOB .EQ. 0, NO PIVOTING IS DONE.
C                IF JOB .NE. 0, PIVOTING IS DONE.
C
C     ON RETURN
C
C         A      A CONTAINS IN ITS UPPER HALF THE CHOLESKY FACTOR
C                OF THE MATRIX A AS IT HAS BEEN PERMUTED BY PIVOTING.
C
C         JPVT   JPVT(J) CONTAINS THE INDEX OF THE DIAGONAL ELEMENT
C                OF A THAT WAS MOVED INTO THE J-TH POSITION,
C                PROVIDED PIVOTING WAS REQUESTED.
C
C         INFO   CONTAINS THE INDEX OF THE LAST POSITIVE DIAGONAL
C                ELEMENT OF THE CHOLESKY FACTOR.
C
C     FOR POSITIVE DEFINITE MATRICES INFO = P IS THE NORMAL RETURN.
C     FOR PIVOTING WITH POSITIVE SEMIDEFINITE MATRICES INFO WILL
C     IN GENERAL BE LESS THAN P.  HOWEVER, INFO MAY BE GREATER THAN
C     THE RANK OF A, SINCE ROUNDING ERROR CAN CAUSE AN OTHERWISE ZERO
C     ELEMENT TO BE POSITIVE. INDEFINITE SYSTEMS WILL ALWAYS CAUSE
C     INFO TO BE LESS THAN P.
C
C     LINPACK. THIS VERSION DATED 03/19/79 .
C     J.J. DONGARRA AND G.W. STEWART, ARGONNE NATIONAL LABORATORY AND
C     UNIVERSITY OF MARYLAND.
C
C
C     BLAS CAXPY,CSWAP
C     FORTRAN SQRT,REAL,CONJG
C
C     INTERNAL VARIABLES
C
      INTEGER PU,PL,PLP1,I,J,JP,JT,K,KB,KM1,KP1,L,MAXL
      COMPLEX TEMP
      REAL MAXDIA
      LOGICAL SWAPK,NEGK
C
      PL = 1
      PU = 0
      INFO = P
      IF (JOB .EQ. 0) GO TO 160
C
C        PIVOTING HAS BEEN REQUESTED. REARRANGE THE
C        THE ELEMENTS ACCORDING TO JPVT.
C
         DO 70 K = 1, P
            SWAPK = JPVT(K) .GT. 0
            NEGK = JPVT(K) .LT. 0
            JPVT(K) = K
            IF (NEGK) JPVT(K) = -JPVT(K)
            IF (.NOT.SWAPK) GO TO 60
               IF (K .EQ. PL) GO TO 50
                  CALL CSWAP(PL-1,A(1,K),1,A(1,PL),1)
                  TEMP = A(K,K)
                  A(K,K) = A(PL,PL)
                  A(PL,PL) = TEMP
                  A(PL,K) = CONJG(A(PL,K))
                  PLP1 = PL + 1
                  IF (P .LT. PLP1) GO TO 40
                  DO 30 J = PLP1, P
                     IF (J .GE. K) GO TO 10
                        TEMP = CONJG(A(PL,J))
                        A(PL,J) = CONJG(A(J,K))
                        A(J,K) = TEMP
                     GO TO 20
   10                CONTINUE
                     IF (J .EQ. K) GO TO 20
                        TEMP = A(K,J)
                        A(K,J) = A(PL,J)
                        A(PL,J) = TEMP
   20                CONTINUE
   30             CONTINUE
   40             CONTINUE
                  JPVT(K) = JPVT(PL)
                  JPVT(PL) = K
   50          CONTINUE
               PL = PL + 1
   60       CONTINUE
   70    CONTINUE
         PU = P
         IF (P .LT. PL) GO TO 150
         DO 140 KB = PL, P
            K = P - KB + PL
            IF (JPVT(K) .GE. 0) GO TO 130
               JPVT(K) = -JPVT(K)
               IF (PU .EQ. K) GO TO 120
                  CALL CSWAP(K-1,A(1,K),1,A(1,PU),1)
                  TEMP = A(K,K)
                  A(K,K) = A(PU,PU)
                  A(PU,PU) = TEMP
                  A(K,PU) = CONJG(A(K,PU))
                  KP1 = K + 1
                  IF (P .LT. KP1) GO TO 110
                  DO 100 J = KP1, P
                     IF (J .GE. PU) GO TO 80
                        TEMP = CONJG(A(K,J))
                        A(K,J) = CONJG(A(J,PU))
                        A(J,PU) = TEMP
                     GO TO 90
   80                CONTINUE
                     IF (J .EQ. PU) GO TO 90
                        TEMP = A(K,J)
                        A(K,J) = A(PU,J)
                        A(PU,J) = TEMP
   90                CONTINUE
  100             CONTINUE
  110             CONTINUE
                  JT = JPVT(K)
                  JPVT(K) = JPVT(PU)
                  JPVT(PU) = JT
  120          CONTINUE
               PU = PU - 1
  130       CONTINUE
  140    CONTINUE
  150    CONTINUE
  160 CONTINUE
      DO 270 K = 1, P
C
C        REDUCTION LOOP.
C
         MAXDIA = REAL(A(K,K))
         KP1 = K + 1
         MAXL = K
C
C        DETERMINE THE PIVOT ELEMENT.
C
         IF (K .LT. PL .OR. K .GE. PU) GO TO 190
            DO 180 L = KP1, PU
               IF (REAL(A(L,L)) .LE. MAXDIA) GO TO 170
                  MAXDIA = REAL(A(L,L))
                  MAXL = L
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
C
C        QUIT IF THE PIVOT ELEMENT IS NOT POSITIVE.
C
         IF (MAXDIA .GT. 0.0E0) GO TO 200
            INFO = K - 1
C     ......EXIT
            GO TO 280
  200    CONTINUE
         IF (K .EQ. MAXL) GO TO 210
C
C           START THE PIVOTING AND UPDATE JPVT.
C
            KM1 = K - 1
            CALL CSWAP(KM1,A(1,K),1,A(1,MAXL),1)
            A(MAXL,MAXL) = A(K,K)
            A(K,K) = CMPLX(MAXDIA,0.0E0)
            JP = JPVT(MAXL)
            JPVT(MAXL) = JPVT(K)
            JPVT(K) = JP
            A(K,MAXL) = CONJG(A(K,MAXL))
  210    CONTINUE
C
C        REDUCTION STEP. PIVOTING IS CONTAINED ACROSS THE ROWS.
C
         WORK(K) = CMPLX(SQRT(REAL(A(K,K))),0.0E0)
         A(K,K) = WORK(K)
         IF (P .LT. KP1) GO TO 260
         DO 250 J = KP1, P
            IF (K .EQ. MAXL) GO TO 240
               IF (J .GE. MAXL) GO TO 220
                  TEMP = CONJG(A(K,J))
                  A(K,J) = CONJG(A(J,MAXL))
                  A(J,MAXL) = TEMP
               GO TO 230
  220          CONTINUE
               IF (J .EQ. MAXL) GO TO 230
                  TEMP = A(K,J)
                  A(K,J) = A(MAXL,J)
                  A(MAXL,J) = TEMP
  230          CONTINUE
  240       CONTINUE
            A(K,J) = A(K,J)/WORK(K)
            WORK(J) = CONJG(A(K,J))
            TEMP = -A(K,J)
            CALL CAXPY(J-K,TEMP,WORK(KP1),1,A(KP1,J),1)
  250    CONTINUE
  260    CONTINUE
  270 CONTINUE
  280 CONTINUE
      RETURN
      END
      SUBROUTINE CCHUD(R,LDR,P,X,Z,LDZ,NZ,Y,RHO,C,S)
      INTEGER LDR,P,LDZ,NZ
      REAL RHO(1),C(1)
      COMPLEX R(LDR,1),X(1),Z(LDZ,1),Y(1),S(1)
C
C     CCHUD UPDATES AN AUGMENTED CHOLESKY DECOMPOSITION OF THE
C     TRIANGULAR PART OF AN AUGMENTED QR DECOMPOSITION.  SPECIFICALLY,
C     GIVEN AN UPPER TRIANGULAR MATRIX R OF ORDER P, A ROW VECTOR
C     X, A COLUMN VECTOR Z, AND A SCALAR Y, CCHUD DETERMINES A
C     UNTIARY MATRIX U AND A SCALAR ZETA SUCH THAT
C
C
C                              (R  Z)     (RR   ZZ )
C                         U  * (    )  =  (        ) ,
C                              (X  Y)     ( 0  ZETA)
C
C     WHERE RR IS UPPER TRIANGULAR.  IF R AND Z HAVE BEEN
C     OBTAINED FROM THE FACTORIZATION OF A LEAST SQUARES
C     PROBLEM, THEN RR AND ZZ ARE THE FACTORS CORRESPONDING TO
C     THE PROBLEM WITH THE OBSERVATION (X,Y) APPENDED.  IN THIS
C     CASE, IF RHO IS THE NORM OF THE RESIDUAL VECTOR, THEN THE
C     NORM OF THE RESIDUAL VECTOR OF THE UPDATED PROBLEM IS
C     SQRT(RHO**2 + ZETA**2).  CCHUD WILL SIMULTANEOUSLY UPDATE
C     SEVERAL TRIPLETS (Z,Y,RHO).
C     FOR A LESS TERSE DESCRIPTION OF WHAT CCHUD DOES AND HOW
C     IT MAY BE APPLIED SEE THE LINPACK GUIDE.
C
C     THE MATRIX U IS DETERMINED AS THE PRODUCT U(P)*...*U(1),
C     WHERE U(I) IS A ROTATION IN THE (I,P+1) PLANE OF THE
C     FORM
C
C                       (     C(I)      S(I) )
C                       (                    ) .
C                       ( -CONJG(S(I))  C(I) )
C
C     THE ROTATIONS ARE CHOSEN SO THAT C(I) IS REAL.
C
C     ON ENTRY
C
C         R      COMPLEX(LDR,P), WHERE LDR .GE. P.
C                R CONTAINS THE UPPER TRIANGULAR MATRIX
C                THAT IS TO BE UPDATED.  THE PART OF R
C                BELOW THE DIAGONAL IS NOT REFERENCED.
C
C         LDR    INTEGER.
C                LDR IS THE LEADING DIMENSION OF THE ARRAY R.
C
C         P      INTEGER.
C                P IS THE ORDER OF THE MATRIX R.
C
C         X      COMPLEX(P).
C                X CONTAINS THE ROW TO BE ADDED TO R.  X IS
C                NOT ALTERED BY CCHUD.
C
C         Z      COMPLEX(LDZ,NZ), WHERE LDZ .GE. P.
C                Z IS AN ARRAY CONTAINING NZ P-VECTORS TO
C                BE UPDATED WITH R.
C
C         LDZ    INTEGER.
C                LDZ IS THE LEADING DIMENSION OF THE ARRAY Z.
C
C         NZ     INTEGER.
C                NZ IS THE NUMBER OF VECTORS TO BE UPDATED
C                NZ MAY BE ZERO, IN WHICH CASE Z, Y, AND RHO
C                ARE NOT REFERENCED.
C
C         Y      COMPLEX(NZ).
C                Y CONTAINS THE SCALARS FOR UPDATING THE VECTORS
C                Z.  Y IS NOT ALTERED BY CCHUD.
C
C         RHO    REAL(NZ).
C                RHO CONTAINS THE NORMS OF THE RESIDUAL
C                VECTORS THAT ARE TO BE UPDATED.  IF RHO(J)
C                IS NEGATIVE, IT IS LEFT UNALTERED.
C
C     ON RETURN
C
C         RC
C         RHO    CONTAIN THE UPDATED QUANTITIES.
C         Z
C
C         C      REAL(P).
C                C CONTAINS THE COSINES OF THE TRANSFORMING
C                ROTATIONS.
C
C         S      COMPLEX(P).
C                S CONTAINS THE SINES OF THE TRANSFORMING
C                ROTATIONS.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     CCHUD USES THE FOLLOWING FUNCTIONS AND SUBROUTINES.
C
C     EXTENDED BLAS CROTG
C     FORTRAN CONJG,SQRT
C
      INTEGER I,J,JM1
      REAL AZETA,SCALE
      COMPLEX T,XJ,ZETA
C
C     UPDATE R.
C
      DO 30 J = 1, P
         XJ = X(J)
C
C        APPLY THE PREVIOUS ROTATIONS.
C
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            T = C(I)*R(I,J) + S(I)*XJ
            XJ = C(I)*XJ - CONJG(S(I))*R(I,J)
            R(I,J) = T
   10    CONTINUE
   20    CONTINUE
C
C        COMPUTE THE NEXT ROTATION.
C
         CALL CROTG(R(J,J),XJ,C(J),S(J))
   30 CONTINUE
C
C     IF REQUIRED, UPDATE Z AND RHO.
C
      IF (NZ .LT. 1) GO TO 70
      DO 60 J = 1, NZ
         ZETA = Y(J)
         DO 40 I = 1, P
            T = C(I)*Z(I,J) + S(I)*ZETA
            ZETA = C(I)*ZETA - CONJG(S(I))*Z(I,J)
            Z(I,J) = T
   40    CONTINUE
         AZETA = CABS(ZETA)
         IF (AZETA .EQ. 0.0E0 .OR. RHO(J) .LT. 0.0E0) GO TO 50
            SCALE = AZETA + RHO(J)
            RHO(J) = SCALE*SQRT((AZETA/SCALE)**2+(RHO(J)/SCALE)**2)
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      RETURN
      END
      SUBROUTINE CCHDD(R,LDR,P,X,Z,LDZ,NZ,Y,RHO,C,S,INFO)
      INTEGER LDR,P,LDZ,NZ,INFO
      COMPLEX R(LDR,1),X(1),Z(LDZ,1),Y(1),S(1)
      REAL RHO(1),C(1)
C
C     CCHDD DOWNDATES AN AUGMENTED CHOLESKY DECOMPOSITION OR THE
C     TRIANGULAR FACTOR OF AN AUGMENTED QR DECOMPOSITION.
C     SPECIFICALLY, GIVEN AN UPPER TRIANGULAR MATRIX R OF ORDER P,  A
C     ROW VECTOR X, A COLUMN VECTOR Z, AND A SCALAR Y, CCHDD
C     DETERMINEDS A UNITARY MATRIX U AND A SCALAR ZETA SUCH THAT
C
C                        (R   Z )     (RR  ZZ)
C                    U * (      )  =  (      ) ,
C                        (0 ZETA)     ( X   Y)
C
C     WHERE RR IS UPPER TRIANGULAR.  IF R AND Z HAVE BEEN OBTAINED
C     FROM THE FACTORIZATION OF A LEAST SQUARES PROBLEM, THEN
C     RR AND ZZ ARE THE FACTORS CORRESPONDING TO THE PROBLEM
C     WITH THE OBSERVATION (X,Y) REMOVED.  IN THIS CASE, IF RHO
C     IS THE NORM OF THE RESIDUAL VECTOR, THEN THE NORM OF
C     THE RESIDUAL VECTOR OF THE DOWNDATED PROBLEM IS
C     SQRT(RHO**2 - ZETA**2). CCHDD WILL SIMULTANEOUSLY DOWNDATE
C     SEVERAL TRIPLETS (Z,Y,RHO) ALONG WITH R.
C     FOR A LESS TERSE DESCRIPTION OF WHAT CCHDD DOES AND HOW
C     IT MAY BE APPLIED, SEE THE LINPACK GUIDE.
C
C     THE MATRIX U IS DETERMINED AS THE PRODUCT U(1)*...*U(P)
C     WHERE U(I) IS A ROTATION IN THE (P+1,I)-PLANE OF THE
C     FORM
C
C                       ( C(I)  -CONJG(S(I)) )
C                       (                    ) .
C                       ( S(I)       C(I)    )
C
C     THE ROTATIONS ARE CHOSEN SO THAT C(I) IS REAL.
C
C     THE USER IS WARNED THAT A GIVEN DOWNDATING PROBLEM MAY
C     BE IMPOSSIBLE TO ACCOMPLISH OR MAY PRODUCE
C     INACCURATE RESULTS.  FOR EXAMPLE, THIS CAN HAPPEN
C     IF X IS NEAR A VECTOR WHOSE REMOVAL WILL REDUCE THE
C     RANK OF R.  BEWARE.
C
C     ON ENTRY
C
C         R      COMPLEX(LDR,P), WHERE LDR .GE. P.
C                R CONTAINS THE UPPER TRIANGULAR MATRIX
C                THAT IS TO BE DOWNDATED.  THE PART OF  R
C                BELOW THE DIAGONAL IS NOT REFERENCED.
C
C         LDR    INTEGER.
C                LDR IS THE LEADING DIMENSION FO THE ARRAY R.
C
C         P      INTEGER.
C                P IS THE ORDER OF THE MATRIX R.
C
C         X      COMPLEX(P).
C                X CONTAINS THE ROW VECTOR THAT IS TO
C                BE REMOVED FROM R.  X IS NOT ALTERED BY CCHDD.
C
C         Z      COMPLEX(LDZ,NZ), WHERE LDZ .GE. P.
C                Z IS AN ARRAY OF NZ P-VECTORS WHICH
C                ARE TO BE DOWNDATED ALONG WITH R.
C
C         LDZ    INTEGER.
C                LDZ IS THE LEADING DIMENSION OF THE ARRAY Z.
C
C         NZ     INTEGER.
C                NZ IS THE NUMBER OF VECTORS TO BE DOWNDATED
C                NZ MAY BE ZERO, IN WHICH CASE Z, Y, AND RHO
C                ARE NOT REFERENCED.
C
C         Y      COMPLEX(NZ).
C                Y CONTAINS THE SCALARS FOR THE DOWNDATING
C                OF THE VECTORS Z.  Y IS NOT ALTERED BY CCHDD.
C
C         RHO    REAL(NZ).
C                RHO CONTAINS THE NORMS OF THE RESIDUAL
C                VECTORS THAT ARE TO BE DOWNDATED.
C
C     ON RETURN
C
C         R
C         Z      CONTAIN THE DOWNDATED QUANTITIES.
C         RHO
C
C         C      REAL(P).
C                C CONTAINS THE COSINES OF THE TRANSFORMING
C                ROTATIONS.
C
C         S      COMPLEX(P).
C                S CONTAINS THE SINES OF THE TRANSFORMING
C                ROTATIONS.
C
C         INFO   INTEGER.
C                INFO IS SET AS FOLLOWS.
C
C                   INFO = 0  IF THE ENTIRE DOWNDATING
C                             WAS SUCCESSFUL.
C
C                   INFO =-1  IF R COULD NOT BE DOWNDATED.
C                             IN THIS CASE, ALL QUANTITIES
C                             ARE LEFT UNALTERED.
C
C                   INFO = 1  IF SOME RHO COULD NOT BE
C                             DOWNDATED.  THE OFFENDING RHOS ARE
C                             SET TO -1.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     CCHDD USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
C
C     FORTRAN AIMAG,CABS,CONJG
C     BLAS CDOTC, SCNRM2
C
      INTEGER I,II,J
      REAL A,ALPHA,AZETA,NORM,SCNRM2
      COMPLEX CDOTC,T,ZETA,B,XX
C
C     SOLVE THE SYSTEM CTRANS(R)*A = X, PLACING THE RESULT
C     IN THE ARRAY S.
C
      INFO = 0
      S(1) = CONJG(X(1))/CONJG(R(1,1))
      IF (P .LT. 2) GO TO 20
      DO 10 J = 2, P
         S(J) = CONJG(X(J)) - CDOTC(J-1,R(1,J),1,S,1)
         S(J) = S(J)/CONJG(R(J,J))
   10 CONTINUE
   20 CONTINUE
      NORM = SCNRM2(P,S,1)
      IF (NORM .LT. 1.0E0) GO TO 30
         INFO = -1
      GO TO 120
   30 CONTINUE
         ALPHA = SQRT(1.0E0-NORM**2)
C
C        DETERMINE THE TRANSFORMATIONS.
C
         DO 40 II = 1, P
            I = P - II + 1
            SCALE = ALPHA + CABS(S(I))
            A = ALPHA/SCALE
            B = S(I)/SCALE
            NORM = SQRT(A**2+REAL(B)**2+AIMAG(B)**2)
            C(I) = A/NORM
            S(I) = CONJG(B)/NORM
            ALPHA = SCALE*NORM
   40    CONTINUE
C
C        APPLY THE TRANSFORMATIONS TO R.
C
         DO 60 J = 1, P
            XX = (0.0E0,0.0E0)
            DO 50 II = 1, J
               I = J - II + 1
               T = C(I)*XX + S(I)*R(I,J)
               R(I,J) = C(I)*R(I,J) - CONJG(S(I))*XX
               XX = T
   50       CONTINUE
   60    CONTINUE
C
C        IF REQUIRED, DOWNDATE Z AND RHO.
C
         IF (NZ .LT. 1) GO TO 110
         DO 100 J = 1, NZ
            ZETA = Y(J)
            DO 70 I = 1, P
               Z(I,J) = (Z(I,J) - CONJG(S(I))*ZETA)/C(I)
               ZETA = C(I)*ZETA - S(I)*Z(I,J)
   70       CONTINUE
            AZETA = CABS(ZETA)
            IF (AZETA .LE. RHO(J)) GO TO 80
               INFO = 1
               RHO(J) = -1.0E0
            GO TO 90
   80       CONTINUE
               RHO(J) = RHO(J)*SQRT(1.0E0-(AZETA/RHO(J))**2)
   90       CONTINUE
  100    CONTINUE
  110    CONTINUE
  120 CONTINUE
      RETURN
      END
      SUBROUTINE CCHEX(R,LDR,P,K,L,Z,LDZ,NZ,C,S,JOB)
      INTEGER LDR,P,K,L,LDZ,NZ,JOB
      COMPLEX R(LDR,1),Z(LDZ,1),S(1)
      REAL C(1)
C
C     CCHEX UPDATES THE CHOLESKY FACTORIZATION
C
C                   A = CTRANS(R)*R
C
C     OF A POSITIVE DEFINITE MATRIX A OF ORDER P UNDER DIAGONAL
C     PERMUTATIONS OF THE FORM
C
C                   TRANS(E)*A*E
C
C     WHERE E IS A PERMUTATION MATRIX.  SPECIFICALLY, GIVEN
C     AN UPPER TRIANGULAR MATRIX R AND A PERMUTATION MATRIX
C     E (WHICH IS SPECIFIED BY K, L, AND JOB), CCHEX DETERMINES
C     A UNITARY MATRIX U SUCH THAT
C
C                           U*R*E = RR,
C
C     WHERE RR IS UPPER TRIANGULAR.  AT THE USERS OPTION, THE
C     TRANSFORMATION U WILL BE MULTIPLIED INTO THE ARRAY Z.
C     IF A = CTRANS(X)*X, SO THAT R IS THE TRIANGULAR PART OF THE
C     QR FACTORIZATION OF X, THEN RR IS THE TRIANGULAR PART OF THE
C     QR FACTORIZATION OF X*E, I.E. X WITH ITS COLUMNS PERMUTED.
C     FOR A LESS TERSE DESCRIPTION OF WHAT CCHEX DOES AND HOW
C     IT MAY BE APPLIED, SEE THE LINPACK GUIDE.
C
C     THE MATRIX Q IS DETERMINED AS THE PRODUCT U(L-K)*...*U(1)
C     OF PLANE ROTATIONS OF THE FORM
C
C                           (    C(I)       S(I) )
C                           (                    ) ,
C                           ( -CONJG(S(I))  C(I) )
C
C     WHERE C(I) IS REAL, THE ROWS THESE ROTATIONS OPERATE ON
C     ARE DESCRIBED BELOW.
C
C     THERE ARE TWO TYPES OF PERMUTATIONS, WHICH ARE DETERMINED
C     BY THE VALUE OF JOB.
C
C     1. RIGHT CIRCULAR SHIFT (JOB = 1).
C
C         THE COLUMNS ARE REARRANGED IN THE FOLLOWING ORDER.
C
C                1,...,K-1,L,K,K+1,...,L-1,L+1,...,P.
C
C         U IS THE PRODUCT OF L-K ROTATIONS U(I), WHERE U(I)
C         ACTS IN THE (L-I,L-I+1)-PLANE.
C
C     2. LEFT CIRCULAR SHIFT (JOB = 2).
C         THE COLUMNS ARE REARRANGED IN THE FOLLOWING ORDER
C
C                1,...,K-1,K+1,K+2,...,L,K,L+1,...,P.
C
C         U IS THE PRODUCT OF L-K ROTATIONS U(I), WHERE U(I)
C         ACTS IN THE (K+I-1,K+I)-PLANE.
C
C     ON ENTRY
C
C         R      COMPLEX(LDR,P), WHERE LDR.GE.P.
C                R CONTAINS THE UPPER TRIANGULAR FACTOR
C                THAT IS TO BE UPDATED.  ELEMENTS OF R
C                BELOW THE DIAGONAL ARE NOT REFERENCED.
C
C         LDR    INTEGER.
C                LDR IS THE LEADING DIMENSION OF THE ARRAY R.
C
C         P      INTEGER.
C                P IS THE ORDER OF THE MATRIX R.
C
C         K      INTEGER.
C                K IS THE FIRST COLUMN TO BE PERMUTED.
C
C         L      INTEGER.
C                L IS THE LAST COLUMN TO BE PERMUTED.
C                L MUST BE STRICTLY GREATER THAN K.
C
C         Z      COMPLEX(LDZ,NZ), WHERE LDZ.GE.P.
C                Z IS AN ARRAY OF NZ P-VECTORS INTO WHICH THE
C                TRANSFORMATION U IS MULTIPLIED.  Z IS
C                NOT REFERENCED IF NZ = 0.
C
C         LDZ    INTEGER.
C                LDZ IS THE LEADING DIMENSION OF THE ARRAY Z.
C
C         NZ     INTEGER.
C                NZ IS THE NUMBER OF COLUMNS OF THE MATRIX Z.
C
C         JOB    INTEGER.
C                JOB DETERMINES THE TYPE OF PERMUTATION.
C                       JOB = 1  RIGHT CIRCULAR SHIFT.
C                       JOB = 2  LEFT CIRCULAR SHIFT.
C
C     ON RETURN
C
C         R      CONTAINS THE UPDATED FACTOR.
C
C         Z      CONTAINS THE UPDATED MATRIX Z.
C
C         C      REAL(P).
C                C CONTAINS THE COSINES OF THE TRANSFORMING ROTATIONS.
C
C         S      COMPLEX(P).
C                S CONTAINS THE SINES OF THE TRANSFORMING ROTATIONS.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     CCHEX USES THE FOLLOWING FUNCTIONS AND SUBROUTINES.
C
C     EXTENDED BLAS CROTG
C     FORTRAN MIN0
C
      INTEGER I,II,IL,IU,J,JJ,KM1,KP1,LMK,LM1
      COMPLEX RJP1J,T
C
C     INITIALIZE
C
      KM1 = K - 1
      KP1 = K + 1
      LMK = L - K
      LM1 = L - 1
C
C     PERFORM THE APPROPRIATE TASK.
C
      GO TO (10,130), JOB
C
C     RIGHT CIRCULAR SHIFT.
C
   10 CONTINUE
C
C        REORDER THE COLUMNS.
C
         DO 20 I = 1, L
            II = L - I + 1
            S(I) = R(II,L)
   20    CONTINUE
         DO 40 JJ = K, LM1
            J = LM1 - JJ + K
            DO 30 I = 1, J
               R(I,J+1) = R(I,J)
   30       CONTINUE
            R(J+1,J+1) = (0.0E0,0.0E0)
   40    CONTINUE
         IF (K .EQ. 1) GO TO 60
            DO 50 I = 1, KM1
               II = L - I + 1
               R(I,K) = S(II)
   50       CONTINUE
   60    CONTINUE
C
C        CALCULATE THE ROTATIONS.
C
         T = S(1)
         DO 70 I = 1, LMK
            CALL CROTG(S(I+1),T,C(I),S(I))
            T = S(I+1)
   70    CONTINUE
         R(K,K) = T
         DO 90 J = KP1, P
            IL = MAX0(1,L-J+1)
            DO 80 II = IL, LMK
               I = L - II
               T = C(II)*R(I,J) + S(II)*R(I+1,J)
               R(I+1,J) = C(II)*R(I+1,J) - CONJG(S(II))*R(I,J)
               R(I,J) = T
   80       CONTINUE
   90    CONTINUE
C
C        IF REQUIRED, APPLY THE TRANSFORMATIONS TO Z.
C
         IF (NZ .LT. 1) GO TO 120
         DO 110 J = 1, NZ
            DO 100 II = 1, LMK
               I = L - II
               T = C(II)*Z(I,J) + S(II)*Z(I+1,J)
               Z(I+1,J) = C(II)*Z(I+1,J) - CONJG(S(II))*Z(I,J)
               Z(I,J) = T
  100       CONTINUE
  110    CONTINUE
  120    CONTINUE
      GO TO 260
C
C     LEFT CIRCULAR SHIFT
C
  130 CONTINUE
C
C        REORDER THE COLUMNS
C
         DO 140 I = 1, K
            II = LMK + I
            S(II) = R(I,K)
  140    CONTINUE
         DO 160 J = K, LM1
            DO 150 I = 1, J
               R(I,J) = R(I,J+1)
  150       CONTINUE
            JJ = J - KM1
            S(JJ) = R(J+1,J+1)
  160    CONTINUE
         DO 170 I = 1, K
            II = LMK + I
            R(I,L) = S(II)
  170    CONTINUE
         DO 180 I = KP1, L
            R(I,L) = (0.0E0,0.0E0)
  180    CONTINUE
C
C        REDUCTION LOOP.
C
         DO 220 J = K, P
            IF (J .EQ. K) GO TO 200
C
C              APPLY THE ROTATIONS.
C
               IU = MIN0(J-1,L-1)
               DO 190 I = K, IU
                  II = I - K + 1
                  T = C(II)*R(I,J) + S(II)*R(I+1,J)
                  R(I+1,J) = C(II)*R(I+1,J) - CONJG(S(II))*R(I,J)
                  R(I,J) = T
  190          CONTINUE
  200       CONTINUE
            IF (J .GE. L) GO TO 210
               JJ = J - K + 1
               T = S(JJ)
               CALL CROTG(R(J,J),T,C(JJ),S(JJ))
  210       CONTINUE
  220    CONTINUE
C
C        APPLY THE ROTATIONS TO Z.
C
         IF (NZ .LT. 1) GO TO 250
         DO 240 J = 1, NZ
            DO 230 I = K, LM1
               II = I - KM1
               T = C(II)*Z(I,J) + S(II)*Z(I+1,J)
               Z(I+1,J) = C(II)*Z(I+1,J) - CONJG(S(II))*Z(I,J)
               Z(I,J) = T
  230       CONTINUE
  240    CONTINUE
  250    CONTINUE
  260 CONTINUE
      RETURN
      END
      SUBROUTINE CQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,JOB)
      INTEGER LDX,N,P,JOB
      INTEGER JPVT(1)
      COMPLEX X(LDX,1),QRAUX(1),WORK(1)
C
C     CQRDC USES HOUSEHOLDER TRANSFORMATIONS TO COMPUTE THE QR
C     FACTORIZATION OF AN N BY P MATRIX X.  COLUMN PIVOTING
C     BASED ON THE 2-NORMS OF THE REDUCED COLUMNS MAY BE
C     PERFORMED AT THE USERS OPTION.
C
C     ON ENTRY
C
C        X       COMPLEX(LDX,P), WHERE LDX .GE. N.
C                X CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO BE
C                COMPUTED.
C
C        LDX     INTEGER.
C                LDX IS THE LEADING DIMENSION OF THE ARRAY X.
C
C        N       INTEGER.
C                N IS THE NUMBER OF ROWS OF THE MATRIX X.
C
C        P       INTEGER.
C                P IS THE NUMBER OF COLUMNS OF THE MATRIX X.
C
C        JPVT    INTEGER(P).
C                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION
C                OF THE PIVOT COLUMNS.  THE K-TH COLUMN X(K) OF X
C                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE
C                VALUE OF JPVT(K).
C
C                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL
C                                      COLUMN.
C
C                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE COLUMN.
C
C                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL COLUMN.
C
C                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL COLUMNS
C                ARE MOVED TO THE BEGINNING OF THE ARRAY X AND FINAL
C                COLUMNS TO THE END.  BOTH INITIAL AND FINAL COLUMNS
C                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY
C                FREE COLUMNS ARE MOVED.  AT THE K-TH STAGE OF THE
C                REDUCTION, IF X(K) IS OCCUPIED BY A FREE COLUMN
C                IT IS INTERCHANGED WITH THE FREE COLUMN OF LARGEST
C                REDUCED NORM.  JPVT IS NOT REFERENCED IF
C                JOB .EQ. 0.
C
C        WORK    COMPLEX(P).
C                WORK IS A WORK ARRAY.  WORK IS NOT REFERENCED IF
C                JOB .EQ. 0.
C
C        JOB     INTEGER.
C                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING.
C                IF JOB .EQ. 0, NO PIVOTING IS DONE.
C                IF JOB .NE. 0, PIVOTING IS DONE.
C
C     ON RETURN
C
C        X       X CONTAINS IN ITS UPPER TRIANGLE THE UPPER
C                TRIANGULAR MATRIX R OF THE QR FACTORIZATION.
C                BELOW ITS DIAGONAL X CONTAINS INFORMATION FROM
C                WHICH THE UNITARY PART OF THE DECOMPOSITION
C                CAN BE RECOVERED.  NOTE THAT IF PIVOTING HAS
C                BEEN REQUESTED, THE DECOMPOSITION IS NOT THAT
C                OF THE ORIGINAL MATRIX X BUT THAT OF X
C                WITH ITS COLUMNS PERMUTED AS DESCRIBED BY JPVT.
C
C        QRAUX   COMPLEX(P).
C                QRAUX CONTAINS FURTHER INFORMATION REQUIRED TO RECOVER
C                THE UNITARY PART OF THE DECOMPOSITION.
C
C        JPVT    JPVT(K) CONTAINS THE INDEX OF THE COLUMN OF THE
C                ORIGINAL MATRIX THAT HAS BEEN INTERCHANGED INTO
C                THE K-TH COLUMN, IF PIVOTING WAS REQUESTED.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     CQRDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
C
C     BLAS CAXPY,CDOTC,CSCAL,CSWAP,SCNRM2
C     FORTRAN ABS,AIMAG,AMAX1,CABS,CMPLX,CSQRT,MIN0,REAL
C
C     INTERNAL VARIABLES
C
      INTEGER J,JP,L,LP1,LUP,MAXJ,PL,PU
      REAL MAXNRM,SCNRM2,TT
      COMPLEX CDOTC,NRMXL,T
      LOGICAL NEGJ,SWAPJ
C
      COMPLEX CSIGN,ZDUM,ZDUM1,ZDUM2
      REAL CABS1
      CSIGN(ZDUM1,ZDUM2) = CABS(ZDUM1)*(ZDUM2/CABS(ZDUM2))
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
      PL = 1
      PU = 0
      IF (JOB .EQ. 0) GO TO 60
C
C        PIVOTING HAS BEEN REQUESTED.  REARRANGE THE COLUMNS
C        ACCORDING TO JPVT.
C
         DO 20 J = 1, P
            SWAPJ = JPVT(J) .GT. 0
            NEGJ = JPVT(J) .LT. 0
            JPVT(J) = J
            IF (NEGJ) JPVT(J) = -J
            IF (.NOT.SWAPJ) GO TO 10
               IF (J .NE. PL) CALL CSWAP(N,X(1,PL),1,X(1,J),1)
               JPVT(J) = JPVT(PL)
               JPVT(PL) = J
               PL = PL + 1
   10       CONTINUE
   20    CONTINUE
         PU = P
         DO 50 JJ = 1, P
            J = P - JJ + 1
            IF (JPVT(J) .GE. 0) GO TO 40
               JPVT(J) = -JPVT(J)
               IF (J .EQ. PU) GO TO 30
                  CALL CSWAP(N,X(1,PU),1,X(1,J),1)
                  JP = JPVT(PU)
                  JPVT(PU) = JPVT(J)
                  JPVT(J) = JP
   30          CONTINUE
               PU = PU - 1
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
C
C     COMPUTE THE NORMS OF THE FREE COLUMNS.
C
      IF (PU .LT. PL) GO TO 80
      DO 70 J = PL, PU
         QRAUX(J) = CMPLX(SCNRM2(N,X(1,J),1),0.0E0)
         WORK(J) = QRAUX(J)
   70 CONTINUE
   80 CONTINUE
C
C     PERFORM THE HOUSEHOLDER REDUCTION OF X.
C
      LUP = MIN0(N,P)
      DO 200 L = 1, LUP
         IF (L .LT. PL .OR. L .GE. PU) GO TO 120
C
C           LOCATE THE COLUMN OF LARGEST NORM AND BRING IT
C           INTO THE PIVOT POSITION.
C
            MAXNRM = 0.0E0
            MAXJ = L
            DO 100 J = L, PU
               IF (REAL(QRAUX(J)) .LE. MAXNRM) GO TO 90
                  MAXNRM = REAL(QRAUX(J))
                  MAXJ = J
   90          CONTINUE
  100       CONTINUE
            IF (MAXJ .EQ. L) GO TO 110
               CALL CSWAP(N,X(1,L),1,X(1,MAXJ),1)
               QRAUX(MAXJ) = QRAUX(L)
               WORK(MAXJ) = WORK(L)
               JP = JPVT(MAXJ)
               JPVT(MAXJ) = JPVT(L)
               JPVT(L) = JP
  110       CONTINUE
  120    CONTINUE
         QRAUX(L) = (0.0E0,0.0E0)
         IF (L .EQ. N) GO TO 190
C
C           COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L.
C
            NRMXL = CMPLX(SCNRM2(N-L+1,X(L,L),1),0.0E0)
            IF (CABS1(NRMXL) .EQ. 0.0E0) GO TO 180
               IF (CABS1(X(L,L)) .NE. 0.0E0)
     *            NRMXL = CSIGN(NRMXL,X(L,L))
               CALL CSCAL(N-L+1,(1.0E0,0.0E0)/NRMXL,X(L,L),1)
               X(L,L) = (1.0E0,0.0E0) + X(L,L)
C
C              APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS,
C              UPDATING THE NORMS.
C
               LP1 = L + 1
               IF (P .LT. LP1) GO TO 170
               DO 160 J = LP1, P
                  T = -CDOTC(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
                  CALL CAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
                  IF (J .LT. PL .OR. J .GT. PU) GO TO 150
                  IF (CABS1(QRAUX(J)) .EQ. 0.0E0) GO TO 150
                     TT = 1.0E0 - (CABS(X(L,J))/REAL(QRAUX(J)))**2
                     TT = AMAX1(TT,0.0E0)
                     T = CMPLX(TT,0.0E0)
                     TT = 1.0E0
     *                    + 0.05E0*TT*(REAL(QRAUX(J))/REAL(WORK(J)))**2
                     IF (TT .EQ. 1.0E0) GO TO 130
                        QRAUX(J) = QRAUX(J)*CSQRT(T)
                     GO TO 140
  130                CONTINUE
                        QRAUX(J) = CMPLX(SCNRM2(N-L,X(L+1,J),1),0.0E0)
                        WORK(J) = QRAUX(J)
  140                CONTINUE
  150             CONTINUE
  160          CONTINUE
  170          CONTINUE
C
C              SAVE THE TRANSFORMATION.
C
               QRAUX(L) = X(L,L)
               X(L,L) = -NRMXL
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
      RETURN
      END
      SUBROUTINE CQRSL(X,LDX,N,K,QRAUX,Y,QY,QTY,B,RSD,XB,JOB,INFO)
      INTEGER LDX,N,K,JOB,INFO
      COMPLEX X(LDX,1),QRAUX(1),Y(1),QY(1),QTY(1),B(1),RSD(1),XB(1)
C
C     CQRSL APPLIES THE OUTPUT OF CQRDC TO COMPUTE COORDINATE
C     TRANSFORMATIONS, PROJECTIONS, AND LEAST SQUARES SOLUTIONS.
C     FOR K .LE. MIN(N,P), LET XK BE THE MATRIX
C
C            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K)))
C
C     FORMED FROM COLUMNNS JPVT(1), ... ,JPVT(K) OF THE ORIGINAL
C     N X P MATRIX X THAT WAS INPUT TO CQRDC (IF NO PIVOTING WAS
C     DONE, XK CONSISTS OF THE FIRST K COLUMNS OF X IN THEIR
C     ORIGINAL ORDER).  CQRDC PRODUCES A FACTORED UNITARY MATRIX Q
C     AND AN UPPER TRIANGULAR MATRIX R SUCH THAT
C
C              XK = Q * (R)
C                       (0)
C
C     THIS INFORMATION IS CONTAINED IN CODED FORM IN THE ARRAYS
C     X AND QRAUX.
C
C     ON ENTRY
C
C        X      COMPLEX(LDX,P).
C               X CONTAINS THE OUTPUT OF CQRDC.
C
C        LDX    INTEGER.
C               LDX IS THE LEADING DIMENSION OF THE ARRAY X.
C
C        N      INTEGER.
C               N IS THE NUMBER OF ROWS OF THE MATRIX XK.  IT MUST
C               HAVE THE SAME VALUE AS N IN CQRDC.
C
C        K      INTEGER.
C               K IS THE NUMBER OF COLUMNS OF THE MATRIX XK.  K
C               MUST NNOT BE GREATER THAN MIN(N,P), WHERE P IS THE
C               SAME AS IN THE CALLING SEQUENCE TO CQRDC.
C
C        QRAUX  COMPLEX(P).
C               QRAUX CONTAINS THE AUXILIARY OUTPUT FROM CQRDC.
C
C        Y      COMPLEX(N)
C               Y CONTAINS AN N-VECTOR THAT IS TO BE MANIPULATED
C               BY CQRSL.
C
C        JOB    INTEGER.
C               JOB SPECIFIES WHAT IS TO BE COMPUTED.  JOB HAS
C               THE DECIMAL EXPANSION ABCDE, WITH THE FOLLOWING
C               MEANING.
C
C                    IF A.NE.0, COMPUTE QY.
C                    IF B,C,D, OR E .NE. 0, COMPUTE QTY.
C                    IF C.NE.0, COMPUTE B.
C                    IF D.NE.0, COMPUTE RSD.
C                    IF E.NE.0, COMPUTE XB.
C
C               NOTE THAT A REQUEST TO COMPUTE B, RSD, OR XB
C               AUTOMATICALLY TRIGGERS THE COMPUTATION OF QTY, FOR
C               WHICH AN ARRAY MUST BE PROVIDED IN THE CALLING
C               SEQUENCE.
C
C     ON RETURN
C
C        QY     COMPLEX(N).
C               QY CONNTAINS Q*Y, IF ITS COMPUTATION HAS BEEN
C               REQUESTED.
C
C        QTY    COMPLEX(N).
C               QTY CONTAINS CTRANS(Q)*Y, IF ITS COMPUTATION HAS
C               BEEN REQUESTED.  HERE CTRANS(Q) IS THE CONJUGATE
C               TRANSPOSE OF THE MATRIX Q.
C
C        B      COMPLEX(K)
C               B CONTAINS THE SOLUTION OF THE LEAST SQUARES PROBLEM
C
C                    MINIMIZE NORM2(Y - XK*B),
C
C               IF ITS COMPUTATION HAS BEEN REQUESTED.  (NOTE THAT
C               IF PIVOTING WAS REQUESTED IN CQRDC, THE J-TH
C               COMPONENT OF B WILL BE ASSOCIATED WITH COLUMN JPVT(J)
C               OF THE ORIGINAL MATRIX X THAT WAS INPUT INTO CQRDC.)
C
C        RSD    COMPLEX(N).
C               RSD CONTAINS THE LEAST SQUARES RESIDUAL Y - XK*B,
C               IF ITS COMPUTATION HAS BEEN REQUESTED.  RSD IS
C               ALSO THE ORTHOGONAL PROJECTION OF Y ONTO THE
C               ORTHOGONAL COMPLEMENT OF THE COLUMN SPACE OF XK.
C
C        XB     COMPLEX(N).
C               XB CONTAINS THE LEAST SQUARES APPROXIMATION XK*B,
C               IF ITS COMPUTATION HAS BEEN REQUESTED.  XB IS ALSO
C               THE ORTHOGONAL PROJECTION OF Y ONTO THE COLUMN SPACE
C               OF X.
C
C        INFO   INTEGER.
C               INFO IS ZERO UNLESS THE COMPUTATION OF B HAS
C               BEEN REQUESTED AND R IS EXACTLY SINGULAR.  IN
C               THIS CASE, INFO IS THE INDEX OF THE FIRST ZERO
C               DIAGONAL ELEMENT OF R AND B IS LEFT UNALTERED.
C
C     THE PARAMETERS QY, QTY, B, RSD, AND XB ARE NOT REFERENCED
C     IF THEIR COMPUTATION IS NOT REQUESTED AND IN THIS CASE
C     CAN BE REPLACED BY DUMMY VARIABLES IN THE CALLING PROGRAM.
C     TO SAVE STORAGE, THE USER MAY IN SOME CASES USE THE SAME
C     ARRAY FOR DIFFERENT PARAMETERS IN THE CALLING SEQUENCE.  A
C     FREQUENTLY OCCURING EXAMPLE IS WHEN ONE WISHES TO COMPUTE
C     ANY OF B, RSD, OR XB AND DOES NOT NEED Y OR QTY.  IN THIS
C     CASE ONE MAY IDENTIFY Y, QTY, AND ONE OF B, RSD, OR XB, WHILE
C     PROVIDING SEPARATE ARRAYS FOR ANYTHING ELSE THAT IS TO BE
C     COMPUTED.  THUS THE CALLING SEQUENCE
C
C          CALL CQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO)
C
C     WILL RESULT IN THE COMPUTATION OF B AND RSD, WITH RSD
C     OVERWRITING Y.  MORE GENERALLY, EACH ITEM IN THE FOLLOWING
C     LIST CONTAINS GROUPS OF PERMISSIBLE IDENTIFICATIONS FOR
C     A SINGLE CALLINNG SEQUENCE.
C
C          1. (Y,QTY,B) (RSD) (XB) (QY)
C
C          2. (Y,QTY,RSD) (B) (XB) (QY)
C
C          3. (Y,QTY,XB) (B) (RSD) (QY)
C
C          4. (Y,QY) (QTY,B) (RSD) (XB)
C
C          5. (Y,QY) (QTY,RSD) (B) (XB)
C
C          6. (Y,QY) (QTY,XB) (B) (RSD)
C
C     IN ANY GROUP THE VALUE RETURNED IN THE ARRAY ALLOCATED TO
C     THE GROUP CORRESPONDS TO THE LAST MEMBER OF THE GROUP.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     CQRSL USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
C
C     BLAS CAXPY,CCOPY,CDOTC
C     FORTRAN ABS,AIMAG,MIN0,MOD,REAL
C
C     INTERNAL VARIABLES
C
      INTEGER I,J,JJ,JU,KP1
      COMPLEX CDOTC,T,TEMP
      LOGICAL CB,CQY,CQTY,CR,CXB
C
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
C     SET INFO FLAG.
C
      INFO = 0
C
C     DETERMINE WHAT IS TO BE COMPUTED.
C
      CQY = JOB/10000 .NE. 0
      CQTY = MOD(JOB,10000) .NE. 0
      CB = MOD(JOB,1000)/100 .NE. 0
      CR = MOD(JOB,100)/10 .NE. 0
      CXB = MOD(JOB,10) .NE. 0
      JU = MIN0(K,N-1)
C
C     SPECIAL ACTION WHEN N=1.
C
      IF (JU .NE. 0) GO TO 40
         IF (CQY) QY(1) = Y(1)
         IF (CQTY) QTY(1) = Y(1)
         IF (CXB) XB(1) = Y(1)
         IF (.NOT.CB) GO TO 30
            IF (CABS1(X(1,1)) .NE. 0.0E0) GO TO 10
               INFO = 1
            GO TO 20
   10       CONTINUE
               B(1) = Y(1)/X(1,1)
   20       CONTINUE
   30    CONTINUE
         IF (CR) RSD(1) = (0.0E0,0.0E0)
      GO TO 250
   40 CONTINUE
C
C        SET UP TO COMPUTE QY OR QTY.
C
         IF (CQY) CALL CCOPY(N,Y,1,QY,1)
         IF (CQTY) CALL CCOPY(N,Y,1,QTY,1)
         IF (.NOT.CQY) GO TO 70
C
C           COMPUTE QY.
C
            DO 60 JJ = 1, JU
               J = JU - JJ + 1
               IF (CABS1(QRAUX(J)) .EQ. 0.0E0) GO TO 50
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -CDOTC(N-J+1,X(J,J),1,QY(J),1)/X(J,J)
                  CALL CAXPY(N-J+1,T,X(J,J),1,QY(J),1)
                  X(J,J) = TEMP
   50          CONTINUE
   60       CONTINUE
   70    CONTINUE
         IF (.NOT.CQTY) GO TO 100
C
C           COMPUTE CTRANS(Q)*Y.
C
            DO 90 J = 1, JU
               IF (CABS1(QRAUX(J)) .EQ. 0.0E0) GO TO 80
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -CDOTC(N-J+1,X(J,J),1,QTY(J),1)/X(J,J)
                  CALL CAXPY(N-J+1,T,X(J,J),1,QTY(J),1)
                  X(J,J) = TEMP
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        SET UP TO COMPUTE B, RSD, OR XB.
C
         IF (CB) CALL CCOPY(K,QTY,1,B,1)
         KP1 = K + 1
         IF (CXB) CALL CCOPY(K,QTY,1,XB,1)
         IF (CR .AND. K .LT. N) CALL CCOPY(N-K,QTY(KP1),1,RSD(KP1),1)
         IF (.NOT.CXB .OR. KP1 .GT. N) GO TO 120
            DO 110 I = KP1, N
               XB(I) = (0.0E0,0.0E0)
  110       CONTINUE
  120    CONTINUE
         IF (.NOT.CR) GO TO 140
            DO 130 I = 1, K
               RSD(I) = (0.0E0,0.0E0)
  130       CONTINUE
  140    CONTINUE
         IF (.NOT.CB) GO TO 190
C
C           COMPUTE B.
C
            DO 170 JJ = 1, K
               J = K - JJ + 1
               IF (CABS1(X(J,J)) .NE. 0.0E0) GO TO 150
                  INFO = J
C           ......EXIT
                  GO TO 180
  150          CONTINUE
               B(J) = B(J)/X(J,J)
               IF (J .EQ. 1) GO TO 160
                  T = -B(J)
                  CALL CAXPY(J-1,T,X(1,J),1,B,1)
  160          CONTINUE
  170       CONTINUE
  180       CONTINUE
  190    CONTINUE
         IF (.NOT.CR .AND. .NOT.CXB) GO TO 240
C
C           COMPUTE RSD OR XB AS REQUIRED.
C
            DO 230 JJ = 1, JU
               J = JU - JJ + 1
               IF (CABS1(QRAUX(J)) .EQ. 0.0E0) GO TO 220
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  IF (.NOT.CR) GO TO 200
                     T = -CDOTC(N-J+1,X(J,J),1,RSD(J),1)/X(J,J)
                     CALL CAXPY(N-J+1,T,X(J,J),1,RSD(J),1)
  200             CONTINUE
                  IF (.NOT.CXB) GO TO 210
                     T = -CDOTC(N-J+1,X(J,J),1,XB(J),1)/X(J,J)
                     CALL CAXPY(N-J+1,T,X(J,J),1,XB(J),1)
  210             CONTINUE
                  X(J,J) = TEMP
  220          CONTINUE
  230       CONTINUE
  240    CONTINUE
  250 CONTINUE
      RETURN
      END
      SUBROUTINE CSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,JOB,INFO)
      INTEGER LDX,N,P,LDU,LDV,JOB,INFO
      COMPLEX X(LDX,1),S(1),E(1),U(LDU,1),V(LDV,1),WORK(1)
C
C
C     CSVDC IS A SUBROUTINE TO REDUCE A COMPLEX NXP MATRIX X BY
C     UNITARY TRANSFORMATIONS U AND V TO DIAGONAL FORM.  THE
C     DIAGONAL ELEMENTS S(I) ARE THE SINGULAR VALUES OF X.  THE
C     COLUMNS OF U ARE THE CORRESPONDING LEFT SINGULAR VECTORS,
C     AND THE COLUMNS OF V THE RIGHT SINGULAR VECTORS.
C
C     ON ENTRY
C
C         X         COMPLEX(LDX,P), WHERE LDX.GE.N.
C                   X CONTAINS THE MATRIX WHOSE SINGULAR VALUE
C                   DECOMPOSITION IS TO BE COMPUTED.  X IS
C                   DESTROYED BY CSVDC.
C
C         LDX       INTEGER.
C                   LDX IS THE LEADING DIMENSION OF THE ARRAY X.
C
C         N         INTEGER.
C                   N IS THE NUMBER OF COLUMNS OF THE MATRIX X.
C
C         P         INTEGER.
C                   P IS THE NUMBER OF ROWS OF THE MATRIX X.
C
C         LDU       INTEGER.
C                   LDU IS THE LEADING DIMENSION OF THE ARRAY U
C                   (SEE BELOW).
C
C         LDV       INTEGER.
C                   LDV IS THE LEADING DIMENSION OF THE ARRAY V
C                   (SEE BELOW).
C
C         WORK      COMPLEX(N).
C                   WORK IS A SCRATCH ARRAY.
C
C         JOB       INTEGER.
C                   JOB CONTROLS THE COMPUTATION OF THE SINGULAR
C                   VECTORS.  IT HAS THE DECIMAL EXPANSION AB
C                   WITH THE FOLLOWING MEANING
C
C                        A.EQ.0    DO NOT COMPUTE THE LEFT SINGULAR
C                                  VECTORS.
C                        A.EQ.1    RETURN THE N LEFT SINGULAR VECTORS
C                                  IN U.
C                        A.GE.2    RETURNS THE FIRST MIN(N,P)
C                                  LEFT SINGULAR VECTORS IN U.
C                        B.EQ.0    DO NOT COMPUTE THE RIGHT SINGULAR
C                                  VECTORS.
C                        B.EQ.1    RETURN THE RIGHT SINGULAR VECTORS
C                                  IN V.
C
C     ON RETURN
C
C         S         COMPLEX(MM), WHERE MM=MIN(N+1,P).
C                   THE FIRST MIN(N,P) ENTRIES OF S CONTAIN THE
C                   SINGULAR VALUES OF X ARRANGED IN DESCENDING
C                   ORDER OF MAGNITUDE.
C
C         E         COMPLEX(P).
C                   E ORDINARILY CONTAINS ZEROS.  HOWEVER SEE THE
C                   DISCUSSION OF INFO FOR EXCEPTIONS.
C
C         U         COMPLEX(LDU,K), WHERE LDU.GE.N.  IF JOBA.EQ.1 THEN
C                                   K.EQ.N, IF JOBA.GE.2 THEN
C                                   K.EQ.MIN(N,P).
C                   U CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.
C                   U IS NOT REFERENCED IF JOBA.EQ.0.  IF N.LE.P
C                   OR IF JOBA.GT.2, THEN U MAY BE IDENTIFIED WITH X
C                   IN THE SUBROUTINE CALL.
C
C         V         COMPLEX(LDV,P), WHERE LDV.GE.P.
C                   V CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.
C                   V IS NOT REFERENCED IF JOBB.EQ.0.  IF P.LE.N,
C                   THEN V MAY BE IDENTIFIED WHTH X IN THE
C                   SUBROUTINE CALL.
C
C         INFO      INTEGER.
C                   THE SINGULAR VALUES (AND THEIR CORRESPONDING
C                   SINGULAR VECTORS) S(INFO+1),S(INFO+2),...,S(M)
C                   ARE CORRECT (HERE M=MIN(N,P)).  THUS IF
C                   INFO.EQ.0, ALL THE SINGULAR VALUES AND THEIR
C                   VECTORS ARE CORRECT.  IN ANY EVENT, THE MATRIX
C                   B = CTRANS(U)*X*V IS THE BIDIAGONAL MATRIX
C                   WITH THE ELEMENTS OF S ON ITS DIAGONAL AND THE
C                   ELEMENTS OF E ON ITS SUPER-DIAGONAL (CTRANS(U)
C                   IS THE CONJUGATE-TRANSPOSE OF U).  THUS THE
C                   SINGULAR VALUES OF X AND B ARE THE SAME.
C
C     LINPACK. THIS VERSION DATED 03/19/79 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     CSVDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
C
C     EXTERNAL CSROT
C     BLAS CAXPY,CDOTC,CSCAL,CSWAP,SCNRM2,SROTG
C     FORTRAN ABS,AIMAG,AMAX1,CABS,CMPLX
C     FORTRAN CONJG,MAX0,MIN0,MOD,REAL,SQRT
C
C     INTERNAL VARIABLES
C
      INTEGER I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,MAXIT,
     *        MM,MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1
      COMPLEX CDOTC,T,R
      REAL B,C,CS,EL,EMM1,F,G,SCNRM2,SCALE,SHIFT,SL,SM,SN,SMM1,T1,TEST,
     *     ZTEST
      LOGICAL WANTU,WANTV
C
      COMPLEX CSIGN,ZDUM,ZDUM1,ZDUM2
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
      CSIGN(ZDUM1,ZDUM2) = CABS(ZDUM1)*(ZDUM2/CABS(ZDUM2))
C
C     SET THE MAXIMUM NUMBER OF ITERATIONS.
C
      MAXIT = 30
C
C     DETERMINE WHAT IS TO BE COMPUTED.
C
      WANTU = .FALSE.
      WANTV = .FALSE.
      JOBU = MOD(JOB,100)/10
      NCU = N
      IF (JOBU .GT. 1) NCU = MIN0(N,P)
      IF (JOBU .NE. 0) WANTU = .TRUE.
      IF (MOD(JOB,10) .NE. 0) WANTV = .TRUE.
C
C     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS
C     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.
C
      INFO = 0
      NCT = MIN0(N-1,P)
      NRT = MAX0(0,MIN0(P-2,N))
      LU = MAX0(NCT,NRT)
      IF (LU .LT. 1) GO TO 170
      DO 160 L = 1, LU
         LP1 = L + 1
         IF (L .GT. NCT) GO TO 20
C
C           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND
C           PLACE THE L-TH DIAGONAL IN S(L).
C
            S(L) = CMPLX(SCNRM2(N-L+1,X(L,L),1),0.0E0)
            IF (CABS1(S(L)) .EQ. 0.0E0) GO TO 10
               IF (CABS1(X(L,L)) .NE. 0.0E0) S(L) = CSIGN(S(L),X(L,L))
               CALL CSCAL(N-L+1,1.0E0/S(L),X(L,L),1)
               X(L,L) = (1.0E0,0.0E0) + X(L,L)
   10       CONTINUE
            S(L) = -S(L)
   20    CONTINUE
         IF (P .LT. LP1) GO TO 50
         DO 40 J = LP1, P
            IF (L .GT. NCT) GO TO 30
            IF (CABS1(S(L)) .EQ. 0.0E0) GO TO 30
C
C              APPLY THE TRANSFORMATION.
C
               T = -CDOTC(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
               CALL CAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
   30       CONTINUE
C
C           PLACE THE L-TH ROW OF X INTO  E FOR THE
C           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
C
            E(J) = CONJG(X(L,J))
   40    CONTINUE
   50    CONTINUE
         IF (.NOT.WANTU .OR. L .GT. NCT) GO TO 70
C
C           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK
C           MULTIPLICATION.
C
            DO 60 I = L, N
               U(I,L) = X(I,L)
   60       CONTINUE
   70    CONTINUE
         IF (L .GT. NRT) GO TO 150
C
C           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE
C           L-TH SUPER-DIAGONAL IN E(L).
C
            E(L) = CMPLX(SCNRM2(P-L,E(LP1),1),0.0E0)
            IF (CABS1(E(L)) .EQ. 0.0E0) GO TO 80
               IF (CABS1(E(LP1)) .NE. 0.0E0) E(L) = CSIGN(E(L),E(LP1))
               CALL CSCAL(P-L,1.0E0/E(L),E(LP1),1)
               E(LP1) = (1.0E0,0.0E0) + E(LP1)
   80       CONTINUE
            E(L) = -CONJG(E(L))
            IF (LP1 .GT. N .OR. CABS1(E(L)) .EQ. 0.0E0) GO TO 120
C
C              APPLY THE TRANSFORMATION.
C
               DO 90 I = LP1, N
                  WORK(I) = (0.0E0,0.0E0)
   90          CONTINUE
               DO 100 J = LP1, P
                  CALL CAXPY(N-L,E(J),X(LP1,J),1,WORK(LP1),1)
  100          CONTINUE
               DO 110 J = LP1, P
                  CALL CAXPY(N-L,CONJG(-E(J)/E(LP1)),WORK(LP1),1,
     *                       X(LP1,J),1)
  110          CONTINUE
  120       CONTINUE
            IF (.NOT.WANTV) GO TO 140
C
C              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT
C              BACK MULTIPLICATION.
C
               DO 130 I = LP1, P
                  V(I,L) = E(I)
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE
  170 CONTINUE
C
C     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.
C
      M = MIN0(P,N+1)
      NCTP1 = NCT + 1
      NRTP1 = NRT + 1
      IF (NCT .LT. P) S(NCTP1) = X(NCTP1,NCTP1)
      IF (N .LT. M) S(M) = (0.0E0,0.0E0)
      IF (NRTP1 .LT. M) E(NRTP1) = X(NRTP1,M)
      E(M) = (0.0E0,0.0E0)
C
C     IF REQUIRED, GENERATE U.
C
      IF (.NOT.WANTU) GO TO 300
         IF (NCU .LT. NCTP1) GO TO 200
         DO 190 J = NCTP1, NCU
            DO 180 I = 1, N
               U(I,J) = (0.0E0,0.0E0)
  180       CONTINUE
            U(J,J) = (1.0E0,0.0E0)
  190    CONTINUE
  200    CONTINUE
         IF (NCT .LT. 1) GO TO 290
         DO 280 LL = 1, NCT
            L = NCT - LL + 1
            IF (CABS1(S(L)) .EQ. 0.0E0) GO TO 250
               LP1 = L + 1
               IF (NCU .LT. LP1) GO TO 220
               DO 210 J = LP1, NCU
                  T = -CDOTC(N-L+1,U(L,L),1,U(L,J),1)/U(L,L)
                  CALL CAXPY(N-L+1,T,U(L,L),1,U(L,J),1)
  210          CONTINUE
  220          CONTINUE
               CALL CSCAL(N-L+1,(-1.0E0,0.0E0),U(L,L),1)
               U(L,L) = (1.0E0,0.0E0) + U(L,L)
               LM1 = L - 1
               IF (LM1 .LT. 1) GO TO 240
               DO 230 I = 1, LM1
                  U(I,L) = (0.0E0,0.0E0)
  230          CONTINUE
  240          CONTINUE
            GO TO 270
  250       CONTINUE
               DO 260 I = 1, N
                  U(I,L) = (0.0E0,0.0E0)
  260          CONTINUE
               U(L,L) = (1.0E0,0.0E0)
  270       CONTINUE
  280    CONTINUE
  290    CONTINUE
  300 CONTINUE
C
C     IF IT IS REQUIRED, GENERATE V.
C
      IF (.NOT.WANTV) GO TO 350
         DO 340 LL = 1, P
            L = P - LL + 1
            LP1 = L + 1
            IF (L .GT. NRT) GO TO 320
            IF (CABS1(E(L)) .EQ. 0.0E0) GO TO 320
               DO 310 J = LP1, P
                  T = -CDOTC(P-L,V(LP1,L),1,V(LP1,J),1)/V(LP1,L)
                  CALL CAXPY(P-L,T,V(LP1,L),1,V(LP1,J),1)
  310          CONTINUE
  320       CONTINUE
            DO 330 I = 1, P
               V(I,L) = (0.0E0,0.0E0)
  330       CONTINUE
            V(L,L) = (1.0E0,0.0E0)
  340    CONTINUE
  350 CONTINUE
C
C     TRANSFORM S AND E SO THAT THEY ARE REAL.
C
      DO 380 I = 1, M
         IF (CABS1(S(I)) .EQ. 0.0E0) GO TO 360
            T = CMPLX(CABS(S(I)),0.0E0)
            R = S(I)/T
            S(I) = T
            IF (I .LT. M) E(I) = E(I)/R
            IF (WANTU) CALL CSCAL(N,R,U(1,I),1)
  360    CONTINUE
C     ...EXIT
         IF (I .EQ. M) GO TO 390
         IF (CABS1(E(I)) .EQ. 0.0E0) GO TO 370
            T = CMPLX(CABS(E(I)),0.0E0)
            R = T/E(I)
            E(I) = T
            S(I+1) = S(I+1)*R
            IF (WANTV) CALL CSCAL(P,R,V(1,I+1),1)
  370    CONTINUE
  380 CONTINUE
  390 CONTINUE
C
C     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
C
      MM = M
      ITER = 0
  400 CONTINUE
C
C        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.
C
C     ...EXIT
         IF (M .EQ. 0) GO TO 660
C
C        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET
C        FLAG AND RETURN.
C
         IF (ITER .LT. MAXIT) GO TO 410
            INFO = M
C     ......EXIT
            GO TO 660
  410    CONTINUE
C
C        THIS SECTION OF THE PROGRAM INSPECTS FOR
C        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON
C        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS.
C
C           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M
C           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M
C           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND
C                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP).
C           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).
C
         DO 430 LL = 1, M
            L = M - LL
C        ...EXIT
            IF (L .EQ. 0) GO TO 440
            TEST = CABS(S(L)) + CABS(S(L+1))
            ZTEST = TEST + CABS(E(L))
            IF (ZTEST .NE. TEST) GO TO 420
               E(L) = (0.0E0,0.0E0)
C        ......EXIT
               GO TO 440
  420       CONTINUE
  430    CONTINUE
  440    CONTINUE
         IF (L .NE. M - 1) GO TO 450
            KASE = 4
         GO TO 520
  450    CONTINUE
            LP1 = L + 1
            MP1 = M + 1
            DO 470 LLS = LP1, MP1
               LS = M - LLS + LP1
C           ...EXIT
               IF (LS .EQ. L) GO TO 480
               TEST = 0.0E0
               IF (LS .NE. M) TEST = TEST + CABS(E(LS))
               IF (LS .NE. L + 1) TEST = TEST + CABS(E(LS-1))
               ZTEST = TEST + CABS(S(LS))
               IF (ZTEST .NE. TEST) GO TO 460
                  S(LS) = (0.0E0,0.0E0)
C           ......EXIT
                  GO TO 480
  460          CONTINUE
  470       CONTINUE
  480       CONTINUE
            IF (LS .NE. L) GO TO 490
               KASE = 3
            GO TO 510
  490       CONTINUE
            IF (LS .NE. M) GO TO 500
               KASE = 1
            GO TO 510
  500       CONTINUE
               KASE = 2
               L = LS
  510       CONTINUE
  520    CONTINUE
         L = L + 1
C
C        PERFORM THE TASK INDICATED BY KASE.
C
         GO TO (530, 560, 580, 610), KASE
C
C        DEFLATE NEGLIGIBLE S(M).
C
  530    CONTINUE
            MM1 = M - 1
            F = REAL(E(M-1))
            E(M-1) = (0.0E0,0.0E0)
            DO 550 KK = L, MM1
               K = MM1 - KK + L
               T1 = REAL(S(K))
               CALL SROTG(T1,F,CS,SN)
               S(K) = CMPLX(T1,0.0E0)
               IF (K .EQ. L) GO TO 540
                  F = -SN*REAL(E(K-1))
                  E(K-1) = CS*E(K-1)
  540          CONTINUE
               IF (WANTV) CALL CSROT(P,V(1,K),1,V(1,M),1,CS,SN)
  550       CONTINUE
         GO TO 650
C
C        SPLIT AT NEGLIGIBLE S(L).
C
  560    CONTINUE
            F = REAL(E(L-1))
            E(L-1) = (0.0E0,0.0E0)
            DO 570 K = L, M
               T1 = REAL(S(K))
               CALL SROTG(T1,F,CS,SN)
               S(K) = CMPLX(T1,0.0E0)
               F = -SN*REAL(E(K))
               E(K) = CS*E(K)
               IF (WANTU) CALL CSROT(N,U(1,K),1,U(1,L-1),1,CS,SN)
  570       CONTINUE
         GO TO 650
C
C        PERFORM ONE QR STEP.
C
  580    CONTINUE
C
C           CALCULATE THE SHIFT.
C
            SCALE = AMAX1(CABS(S(M)),CABS(S(M-1)),CABS(E(M-1)),
     *                    CABS(S(L)),CABS(E(L)))
            SM = REAL(S(M))/SCALE
            SMM1 = REAL(S(M-1))/SCALE
            EMM1 = REAL(E(M-1))/SCALE
            SL = REAL(S(L))/SCALE
            EL = REAL(E(L))/SCALE
            B = ((SMM1 + SM)*(SMM1 - SM) + EMM1**2)/2.0E0
            C = (SM*EMM1)**2
            SHIFT = 0.0E0
            IF (B .EQ. 0.0E0 .AND. C .EQ. 0.0E0) GO TO 590
               SHIFT = SQRT(B**2+C)
               IF (B .LT. 0.0E0) SHIFT = -SHIFT
               SHIFT = C/(B + SHIFT)
  590       CONTINUE
            F = (SL + SM)*(SL - SM) - SHIFT
            G = SL*EL
C
C           CHASE ZEROS.
C
            MM1 = M - 1
            DO 600 K = L, MM1
               CALL SROTG(F,G,CS,SN)
               IF (K .NE. L) E(K-1) = CMPLX(F,0.0E0)
               F = CS*REAL(S(K)) + SN*REAL(E(K))
               E(K) = CS*E(K) - SN*S(K)
               G = SN*REAL(S(K+1))
               S(K+1) = CS*S(K+1)
               IF (WANTV) CALL CSROT(P,V(1,K),1,V(1,K+1),1,CS,SN)
               CALL SROTG(F,G,CS,SN)
               S(K) = CMPLX(F,0.0E0)
               F = CS*REAL(E(K)) + SN*REAL(S(K+1))
               S(K+1) = -SN*E(K) + CS*S(K+1)
               G = SN*REAL(E(K+1))
               E(K+1) = CS*E(K+1)
               IF (WANTU .AND. K .LT. N)
     *            CALL CSROT(N,U(1,K),1,U(1,K+1),1,CS,SN)
  600       CONTINUE
            E(M-1) = CMPLX(F,0.0E0)
            ITER = ITER + 1
         GO TO 650
C
C        CONVERGENCE.
C
  610    CONTINUE
C
C           MAKE THE SINGULAR VALUE  POSITIVE
C
            IF (REAL(S(L)) .GE. 0.0E0) GO TO 620
               S(L) = -S(L)
               IF (WANTV) CALL CSCAL(P,(-1.0E0,0.0E0),V(1,L),1)
  620       CONTINUE
C
C           ORDER THE SINGULAR VALUE.
C
  630       IF (L .EQ. MM) GO TO 640
C           ...EXIT
               IF (REAL(S(L)) .GE. REAL(S(L+1))) GO TO 640
               T = S(L)
               S(L) = S(L+1)
               S(L+1) = T
               IF (WANTV .AND. L .LT. P)
     *            CALL CSWAP(P,V(1,L),1,V(1,L+1),1)
               IF (WANTU .AND. L .LT. N)
     *            CALL CSWAP(N,U(1,L),1,U(1,L+1),1)
               L = L + 1
            GO TO 630
  640       CONTINUE
            ITER = 0
            M = M - 1
  650    CONTINUE
      GO TO 400
  660 CONTINUE
      RETURN
      END
