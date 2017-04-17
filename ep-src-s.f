      SUBROUTINE CDIV(AR,AI,BR,BI,CR,CI)
      REAL AR,AI,BR,BI,CR,CI
C
C     COMPLEX DIVISION, (CR,CI) = (AR,AI)/(BR,BI)
C
      REAL S,ARS,AIS,BRS,BIS
      S = ABS(BR) + ABS(BI)
      ARS = AR/S
      AIS = AI/S
      BRS = BR/S
      BIS = BI/S
      S = BRS**2 + BIS**2
      CR = (ARS*BRS + AIS*BIS)/S
      CI = (AIS*BRS - ARS*BIS)/S
      RETURN
      END
      SUBROUTINE CSROOT(XR,XI,YR,YI)
      REAL XR,XI,YR,YI
C
C     (YR,YI) = COMPLEX SQRT(XR,XI)
C     BRANCH CHOSEN SO THAT YR .GE. 0.0 AND SIGN(YI) .EQ. SIGN(XI)
C
      REAL S,TR,TI,PYTHAG
      TR = XR
      TI = XI
      S = SQRT(0.5E0*(PYTHAG(TR,TI) + ABS(TR)))
      IF (TR .GE. 0.0E0) YR = S
      IF (TI .LT. 0.0E0) S = -S
      IF (TR .LE. 0.0E0) YI = S
      IF (TR .LT. 0.0E0) YR = 0.5E0*(TI/YI)
      IF (TR .GT. 0.0E0) YI = 0.5E0*(TI/YR)
      RETURN
      END
      REAL FUNCTION EPSLON (X)
      REAL X
C
C     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.
C
      REAL A,B,C,EPS
C
C     THIS PROGRAM SHOULD FUNCTION PROPERLY ON ALL SYSTEMS
C     SATISFYING THE FOLLOWING TWO ASSUMPTIONS,
C        1.  THE BASE USED IN REPRESENTING FLOATING POINT
C            NUMBERS IS NOT A POWER OF THREE.
C        2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO
C            THE ACCURACY USED IN FLOATING POINT VARIABLES
C            THAT ARE STORED IN MEMORY.
C     THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO
C     FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING
C     ASSUMPTION 2.
C     UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT,
C            A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS,
C            B  HAS A ZERO FOR ITS LAST BIT OR DIGIT,
C            C  IS NOT EXACTLY EQUAL TO ONE,
C            EPS  MEASURES THE SEPARATION OF 1.0 FROM
C                 THE NEXT LARGER FLOATING POINT NUMBER.
C     THE DEVELOPERS OF EISPACK WOULD APPRECIATE BEING INFORMED
C     ABOUT ANY SYSTEMS WHERE THESE ASSUMPTIONS DO NOT HOLD.
C
C     THIS VERSION DATED 4/6/83.
C
      A = 4.0E0/3.0E0
   10 B = A - 1.0E0
      C = B + B + B
      EPS = ABS(C-1.0E0)
      IF (EPS .EQ. 0.0E0) GO TO 10
      EPSLON = EPS*ABS(X)
      RETURN
      END
      REAL FUNCTION PYTHAG(A,B)
      REAL A,B
C
C     FINDS SQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
C
      REAL P,R,S,T,U
      P = AMAX1(ABS(A),ABS(B))
      IF (P .EQ. 0.0E0) GO TO 20
      R = (AMIN1(ABS(A),ABS(B))/P)**2
   10 CONTINUE
         T = 4.0E0 + R
         IF (T .EQ. 4.0E0) GO TO 20
         S = R/T
         U = 1.0E0 + 2.0E0*S
         P = U*P
         R = (S/U)**2 * R
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
      SUBROUTINE BAKVEC(NM,N,T,E,M,Z,IERR)
C
      INTEGER I,J,M,N,NM,IERR
      REAL T(NM,3),E(N),Z(NM,M)
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A NONSYMMETRIC
C     TRIDIAGONAL MATRIX BY BACK TRANSFORMING THOSE OF THE
C     CORRESPONDING SYMMETRIC MATRIX DETERMINED BY  FIGI.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        T CONTAINS THE NONSYMMETRIC MATRIX.  ITS SUBDIAGONAL IS
C          STORED IN THE LAST N-1 POSITIONS OF THE FIRST COLUMN,
C          ITS DIAGONAL IN THE N POSITIONS OF THE SECOND COLUMN,
C          AND ITS SUPERDIAGONAL IN THE FIRST N-1 POSITIONS OF
C          THE THIRD COLUMN.  T(1,1) AND T(N,3) ARE ARBITRARY.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE SYMMETRIC
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
C
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        T IS UNALTERED.
C
C        E IS DESTROYED.
C
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS
C          IN ITS FIRST M COLUMNS.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          2*N+I      IF E(I) IS ZERO WITH T(I,1) OR T(I-1,3) NON-ZERO.
C                     IN THIS CASE, THE SYMMETRIC MATRIX IS NOT SIMILAR
C                     TO THE ORIGINAL MATRIX, AND THE EIGENVECTORS
C                     CANNOT BE FOUND BY THIS PROGRAM.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (M .EQ. 0) GO TO 1001
      E(1) = 1.0E0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
         IF (E(I) .NE. 0.0E0) GO TO 80
         IF (T(I,1) .NE. 0.0E0 .OR. T(I-1,3) .NE. 0.0E0) GO TO 1000
         E(I) = 1.0E0
         GO TO 100
   80    E(I) = E(I-1) * E(I) / T(I-1,3)
  100 CONTINUE
C
      DO 120 J = 1, M
C
         DO 120 I = 2, N
         Z(I,J) = Z(I,J) * E(I)
  120 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- EIGENVECTORS CANNOT BE
C                FOUND BY THIS PROGRAM ..........
 1000 IERR = 2 * N + I
 1001 RETURN
      END
      SUBROUTINE BALANC(NM,N,A,LOW,IGH,SCALE)
C
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      REAL A(NM,N),SCALE(N)
      REAL C,F,G,R,S,B2,RADIX
      LOGICAL NOCONV
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BALANCE,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     THIS SUBROUTINE BALANCES A REAL MATRIX AND ISOLATES
C     EIGENVALUES WHENEVER POSSIBLE.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS THE INPUT MATRIX TO BE BALANCED.
C
C     ON OUTPUT
C
C        A CONTAINS THE BALANCED MATRIX.
C
C        LOW AND IGH ARE TWO INTEGERS SUCH THAT A(I,J)
C          IS EQUAL TO ZERO IF
C           (1) I IS GREATER THAN J AND
C           (2) J=1,...,LOW-1 OR I=IGH+1,...,N.
C
C        SCALE CONTAINS INFORMATION DETERMINING THE
C           PERMUTATIONS AND SCALING FACTORS USED.
C
C     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH
C     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED
C     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS
C     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN
C        SCALE(J) = P(J),    FOR J = 1,...,LOW-1
C                 = D(J,J),      J = LOW,...,IGH
C                 = P(J)         J = IGH+1,...,N.
C     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,
C     THEN 1 TO LOW-1.
C
C     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.
C
C     THE ALGOL PROCEDURE EXC CONTAINED IN BALANCE APPEARS IN
C     BALANC  IN LINE.  (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS
C     K,L HAVE BEEN REVERSED.)
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      RADIX = 16.0E0
C
      B2 = RADIX * RADIX
      K = 1
      L = N
      GO TO 100
C     .......... IN-LINE PROCEDURE FOR ROW AND
C                COLUMN EXCHANGE ..........
   20 SCALE(M) = J
      IF (J .EQ. M) GO TO 50
C
      DO 30 I = 1, L
         F = A(I,J)
         A(I,J) = A(I,M)
         A(I,M) = F
   30 CONTINUE
C
      DO 40 I = K, N
         F = A(J,I)
         A(J,I) = A(M,I)
         A(M,I) = F
   40 CONTINUE
C
   50 GO TO (80,130), IEXC
C     .......... SEARCH FOR ROWS ISOLATING AN EIGENVALUE
C                AND PUSH THEM DOWN ..........
   80 IF (L .EQ. 1) GO TO 280
      L = L - 1
C     .......... FOR J=L STEP -1 UNTIL 1 DO -- ..........
  100 DO 120 JJ = 1, L
         J = L + 1 - JJ
C
         DO 110 I = 1, L
            IF (I .EQ. J) GO TO 110
            IF (A(J,I) .NE. 0.0E0) GO TO 120
  110    CONTINUE
C
         M = L
         IEXC = 1
         GO TO 20
  120 CONTINUE
C
      GO TO 140
C     .......... SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
C                AND PUSH THEM LEFT ..........
  130 K = K + 1
C
  140 DO 170 J = K, L
C
         DO 150 I = K, L
            IF (I .EQ. J) GO TO 150
            IF (A(I,J) .NE. 0.0E0) GO TO 170
  150    CONTINUE
C
         M = K
         IEXC = 2
         GO TO 20
  170 CONTINUE
C     .......... NOW BALANCE THE SUBMATRIX IN ROWS K TO L ..........
      DO 180 I = K, L
  180 SCALE(I) = 1.0E0
C     .......... ITERATIVE LOOP FOR NORM REDUCTION ..........
  190 NOCONV = .FALSE.
C
      DO 270 I = K, L
         C = 0.0E0
         R = 0.0E0
C
         DO 200 J = K, L
            IF (J .EQ. I) GO TO 200
            C = C + ABS(A(J,I))
            R = R + ABS(A(I,J))
  200    CONTINUE
C     .......... GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW ..........
         IF (C .EQ. 0.0E0 .OR. R .EQ. 0.0E0) GO TO 270
         G = R / RADIX
         F = 1.0E0
         S = C + R
  210    IF (C .GE. G) GO TO 220
         F = F * RADIX
         C = C * B2
         GO TO 210
  220    G = R * RADIX
  230    IF (C .LT. G) GO TO 240
         F = F / RADIX
         C = C / B2
         GO TO 230
C     .......... NOW BALANCE ..........
  240    IF ((C + R) / F .GE. 0.95E0 * S) GO TO 270
         G = 1.0E0 / F
         SCALE(I) = SCALE(I) * F
         NOCONV = .TRUE.
C
         DO 250 J = K, N
  250    A(I,J) = A(I,J) * G
C
         DO 260 J = 1, L
  260    A(J,I) = A(J,I) * F
C
  270 CONTINUE
C
      IF (NOCONV) GO TO 190
C
  280 LOW = K
      IGH = L
      RETURN
      END
      SUBROUTINE BALBAK(NM,N,LOW,IGH,SCALE,M,Z)
C
      INTEGER I,J,K,M,N,II,NM,IGH,LOW
      REAL SCALE(N),Z(NM,M)
      REAL S
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BALBAK,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL GENERAL
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     BALANCED MATRIX DETERMINED BY  BALANC.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY  BALANC.
C
C        SCALE CONTAINS INFORMATION DETERMINING THE PERMUTATIONS
C          AND SCALING FACTORS USED BY  BALANC.
C
C        M IS THE NUMBER OF COLUMNS OF Z TO BE BACK TRANSFORMED.
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGEN-
C          VECTORS TO BE BACK TRANSFORMED IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE
C          TRANSFORMED EIGENVECTORS IN ITS FIRST M COLUMNS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
      IF (IGH .EQ. LOW) GO TO 120
C
      DO 110 I = LOW, IGH
         S = SCALE(I)
C     .......... LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
C                IF THE FOREGOING STATEMENT IS REPLACED BY
C                S=1.0E0/SCALE(I). ..........
         DO 100 J = 1, M
  100    Z(I,J) = Z(I,J) * S
C
  110 CONTINUE
C     ......... FOR I=LOW-1 STEP -1 UNTIL 1,
C               IGH+1 STEP 1 UNTIL N DO -- ..........
  120 DO 140 II = 1, N
         I = II
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 140
         IF (I .LT. LOW) I = LOW - II
         K = SCALE(I)
         IF (K .EQ. I) GO TO 140
C
         DO 130 J = 1, M
            S = Z(I,J)
            Z(I,J) = Z(K,J)
            Z(K,J) = S
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE BANDR(NM,N,MB,A,D,E,E2,MATZ,Z)
C
      INTEGER J,K,L,N,R,I1,I2,J1,J2,KR,MB,MR,M1,NM,N2,R1,UGL,MAXL,MAXR
      REAL A(NM,MB),D(N),E(N),E2(N),Z(NM,N)
      REAL G,U,B1,B2,C2,F1,F2,S2,DMIN,DMINRT
      LOGICAL MATZ
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BANDRD,
C     NUM. MATH. 12, 231-241(1968) BY SCHWARZ.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 273-283(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC BAND MATRIX
C     TO A SYMMETRIC TRIDIAGONAL MATRIX USING AND OPTIONALLY
C     ACCUMULATING ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        MB IS THE (HALF) BAND WIDTH OF THE MATRIX, DEFINED AS THE
C          NUMBER OF ADJACENT DIAGONALS, INCLUDING THE PRINCIPAL
C          DIAGONAL, REQUIRED TO SPECIFY THE NON-ZERO PORTION OF THE
C          LOWER TRIANGLE OF THE MATRIX.
C
C        A CONTAINS THE LOWER TRIANGLE OF THE SYMMETRIC BAND INPUT
C          MATRIX STORED AS AN N BY MB ARRAY.  ITS LOWEST SUBDIAGONAL
C          IS STORED IN THE LAST N+1-MB POSITIONS OF THE FIRST COLUMN,
C          ITS NEXT SUBDIAGONAL IN THE LAST N+2-MB POSITIONS OF THE
C          SECOND COLUMN, FURTHER SUBDIAGONALS SIMILARLY, AND FINALLY
C          ITS PRINCIPAL DIAGONAL IN THE N POSITIONS OF THE LAST COLUMN.
C          CONTENTS OF STORAGES NOT PART OF THE MATRIX ARE ARBITRARY.
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE TRANSFORMATION MATRIX IS
C          TO BE ACCUMULATED, AND TO .FALSE. OTHERWISE.
C
C     ON OUTPUT
C
C        A HAS BEEN DESTROYED, EXCEPT FOR ITS LAST TWO COLUMNS WHICH
C          CONTAIN A COPY OF THE TRIDIAGONAL MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX PRODUCED IN
C          THE REDUCTION IF MATZ HAS BEEN SET TO .TRUE.  OTHERWISE, Z
C          IS NOT REFERENCED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      DMIN = 2.0E0**(-64)
      DMINRT = 2.0E0**(-32)
C     .......... INITIALIZE DIAGONAL SCALING MATRIX ..........
      DO 30 J = 1, N
   30 D(J) = 1.0E0
C
      IF (.NOT. MATZ) GO TO 60
C
      DO 50 J = 1, N
C
         DO 40 K = 1, N
   40    Z(J,K) = 0.0E0
C
         Z(J,J) = 1.0E0
   50 CONTINUE
C
   60 M1 = MB - 1
      IF (M1 - 1) 900, 800, 70
   70 N2 = N - 2
C
      DO 700 K = 1, N2
         MAXR = MIN0(M1,N-K)
C     .......... FOR R=MAXR STEP -1 UNTIL 2 DO -- ..........
         DO 600 R1 = 2, MAXR
            R = MAXR + 2 - R1
            KR = K + R
            MR = MB - R
            G = A(KR,MR)
            A(KR-1,1) = A(KR-1,MR+1)
            UGL = K
C
            DO 500 J = KR, N, M1
               J1 = J - 1
               J2 = J1 - 1
               IF (G .EQ. 0.0E0) GO TO 600
               B1 = A(J1,1) / G
               B2 = B1 * D(J1) / D(J)
               S2 = 1.0E0 / (1.0E0 + B1 * B2)
               IF (S2 .GE. 0.5E0 ) GO TO 450
               B1 = G / A(J1,1)
               B2 = B1 * D(J) / D(J1)
               C2 = 1.0E0 - S2
               D(J1) = C2 * D(J1)
               D(J) = C2 * D(J)
               F1 = 2.0E0 * A(J,M1)
               F2 = B1 * A(J1,MB)
               A(J,M1) = -B2 * (B1 * A(J,M1) - A(J,MB)) - F2 + A(J,M1)
               A(J1,MB) = B2 * (B2 * A(J,MB) + F1) + A(J1,MB)
               A(J,MB) = B1 * (F2 - F1) + A(J,MB)
C
               DO 200 L = UGL, J2
                  I2 = MB - J + L
                  U = A(J1,I2+1) + B2 * A(J,I2)
                  A(J,I2) = -B1 * A(J1,I2+1) + A(J,I2)
                  A(J1,I2+1) = U
  200          CONTINUE
C
               UGL = J
               A(J1,1) = A(J1,1) + B2 * G
               IF (J .EQ. N) GO TO 350
               MAXL = MIN0(M1,N-J1)
C
               DO 300 L = 2, MAXL
                  I1 = J1 + L
                  I2 = MB - L
                  U = A(I1,I2) + B2 * A(I1,I2+1)
                  A(I1,I2+1) = -B1 * A(I1,I2) + A(I1,I2+1)
                  A(I1,I2) = U
  300          CONTINUE
C
               I1 = J + M1
               IF (I1 .GT. N) GO TO 350
               G = B2 * A(I1,1)
  350          IF (.NOT. MATZ) GO TO 500
C
               DO 400 L = 1, N
                  U = Z(L,J1) + B2 * Z(L,J)
                  Z(L,J) = -B1 * Z(L,J1) + Z(L,J)
                  Z(L,J1) = U
  400          CONTINUE
C
               GO TO 500
C
  450          U = D(J1)
               D(J1) = S2 * D(J)
               D(J) = S2 * U
               F1 = 2.0E0 * A(J,M1)
               F2 = B1 * A(J,MB)
               U = B1 * (F2 - F1) + A(J1,MB)
               A(J,M1) = B2 * (B1 * A(J,M1) - A(J1,MB)) + F2 - A(J,M1)
               A(J1,MB) = B2 * (B2 * A(J1,MB) + F1) + A(J,MB)
               A(J,MB) = U
C
               DO 460 L = UGL, J2
                  I2 = MB - J + L
                  U = B2 * A(J1,I2+1) + A(J,I2)
                  A(J,I2) = -A(J1,I2+1) + B1 * A(J,I2)
                  A(J1,I2+1) = U
  460          CONTINUE
C
               UGL = J
               A(J1,1) = B2 * A(J1,1) + G
               IF (J .EQ. N) GO TO 480
               MAXL = MIN0(M1,N-J1)
C
               DO 470 L = 2, MAXL
                  I1 = J1 + L
                  I2 = MB - L
                  U = B2 * A(I1,I2) + A(I1,I2+1)
                  A(I1,I2+1) = -A(I1,I2) + B1 * A(I1,I2+1)
                  A(I1,I2) = U
  470          CONTINUE
C
               I1 = J + M1
               IF (I1 .GT. N) GO TO 480
               G = A(I1,1)
               A(I1,1) = B1 * A(I1,1)
  480          IF (.NOT. MATZ) GO TO 500
C
               DO 490 L = 1, N
                  U = B2 * Z(L,J1) + Z(L,J)
                  Z(L,J) = -Z(L,J1) + B1 * Z(L,J)
                  Z(L,J1) = U
  490          CONTINUE
C
  500       CONTINUE
C
  600    CONTINUE
C
         IF (MOD(K,64) .NE. 0) GO TO 700
C     .......... RESCALE TO AVOID UNDERFLOW OR OVERFLOW ..........
         DO 650 J = K, N
            IF (D(J) .GE. DMIN) GO TO 650
            MAXL = MAX0(1,MB+1-J)
C
            DO 610 L = MAXL, M1
  610       A(J,L) = DMINRT * A(J,L)
C
            IF (J .EQ. N) GO TO 630
            MAXL = MIN0(M1,N-J)
C
            DO 620 L = 1, MAXL
               I1 = J + L
               I2 = MB - L
               A(I1,I2) = DMINRT * A(I1,I2)
  620       CONTINUE
C
  630       IF (.NOT. MATZ) GO TO 645
C
            DO 640 L = 1, N
  640       Z(L,J) = DMINRT * Z(L,J)
C
  645       A(J,MB) = DMIN * A(J,MB)
            D(J) = D(J) / DMIN
  650    CONTINUE
C
  700 CONTINUE
C     .......... FORM SQUARE ROOT OF SCALING MATRIX ..........
  800 DO 810 J = 2, N
  810 E(J) = SQRT(D(J))
C
      IF (.NOT. MATZ) GO TO 840
C
      DO 830 J = 1, N
C
         DO 820 K = 2, N
  820    Z(J,K) = E(K) * Z(J,K)
C
  830 CONTINUE
C
  840 U = 1.0E0
C
      DO 850 J = 2, N
         A(J,M1) = U * E(J) * A(J,M1)
         U = E(J)
         E2(J) = A(J,M1) ** 2
         A(J,MB) = D(J) * A(J,MB)
         D(J) = A(J,MB)
         E(J) = A(J,M1)
  850 CONTINUE
C
      D(1) = A(1,MB)
      E(1) = 0.0E0
      E2(1) = 0.0E0
      GO TO 1001
C
  900 DO 950 J = 1, N
         D(J) = A(J,MB)
         E(J) = 0.0E0
         E2(J) = 0.0E0
  950 CONTINUE
C
 1001 RETURN
      END
      SUBROUTINE BANDV(NM,N,MBW,A,E21,M,W,Z,IERR,NV,RV,RV6)
C
      INTEGER I,J,K,M,N,R,II,IJ,JJ,KJ,MB,M1,NM,NV,IJ1,ITS,KJ1,MBW,M21,
     X        IERR,MAXJ,MAXK,GROUP
      REAL A(NM,MBW),W(M),Z(NM,M),RV(NV),RV6(N)
      REAL U,V,UK,XU,X0,X1,E21,EPS2,EPS3,EPS4,NORM,ORDER,
     X       EPSLON,PYTHAG
C
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A REAL SYMMETRIC
C     BAND MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES, USING INVERSE
C     ITERATION.  THE SUBROUTINE MAY ALSO BE USED TO SOLVE SYSTEMS
C     OF LINEAR EQUATIONS WITH A SYMMETRIC OR NON-SYMMETRIC BAND
C     COEFFICIENT MATRIX.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        MBW IS THE NUMBER OF COLUMNS OF THE ARRAY A USED TO STORE THE
C          BAND MATRIX.  IF THE MATRIX IS SYMMETRIC, MBW IS ITS (HALF)
C          BAND WIDTH, DENOTED MB AND DEFINED AS THE NUMBER OF ADJACENT
C          DIAGONALS, INCLUDING THE PRINCIPAL DIAGONAL, REQUIRED TO
C          SPECIFY THE NON-ZERO PORTION OF THE LOWER TRIANGLE OF THE
C          MATRIX.  IF THE SUBROUTINE IS BEING USED TO SOLVE SYSTEMS
C          OF LINEAR EQUATIONS AND THE COEFFICIENT MATRIX IS NOT
C          SYMMETRIC, IT MUST HOWEVER HAVE THE SAME NUMBER OF ADJACENT
C          DIAGONALS ABOVE THE MAIN DIAGONAL AS BELOW, AND IN THIS
C          CASE, MBW=2*MB-1.
C
C        A CONTAINS THE LOWER TRIANGLE OF THE SYMMETRIC BAND INPUT
C          MATRIX STORED AS AN N BY MB ARRAY.  ITS LOWEST SUBDIAGONAL
C          IS STORED IN THE LAST N+1-MB POSITIONS OF THE FIRST COLUMN,
C          ITS NEXT SUBDIAGONAL IN THE LAST N+2-MB POSITIONS OF THE
C          SECOND COLUMN, FURTHER SUBDIAGONALS SIMILARLY, AND FINALLY
C          ITS PRINCIPAL DIAGONAL IN THE N POSITIONS OF COLUMN MB.
C          IF THE SUBROUTINE IS BEING USED TO SOLVE SYSTEMS OF LINEAR
C          EQUATIONS AND THE COEFFICIENT MATRIX IS NOT SYMMETRIC, A IS
C          N BY 2*MB-1 INSTEAD WITH LOWER TRIANGLE AS ABOVE AND WITH
C          ITS FIRST SUPERDIAGONAL STORED IN THE FIRST N-1 POSITIONS OF
C          COLUMN MB+1, ITS SECOND SUPERDIAGONAL IN THE FIRST N-2
C          POSITIONS OF COLUMN MB+2, FURTHER SUPERDIAGONALS SIMILARLY,
C          AND FINALLY ITS HIGHEST SUPERDIAGONAL IN THE FIRST N+1-MB
C          POSITIONS OF THE LAST COLUMN.
C          CONTENTS OF STORAGES NOT PART OF THE MATRIX ARE ARBITRARY.
C
C        E21 SPECIFIES THE ORDERING OF THE EIGENVALUES AND CONTAINS
C            0.0E0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR
C            2.0E0 IF THE EIGENVALUES ARE IN DESCENDING ORDER.
C          IF THE SUBROUTINE IS BEING USED TO SOLVE SYSTEMS OF LINEAR
C          EQUATIONS, E21 SHOULD BE SET TO 1.0E0 IF THE COEFFICIENT
C          MATRIX IS SYMMETRIC AND TO -1.0E0 IF NOT.
C
C        M IS THE NUMBER OF SPECIFIED EIGENVALUES OR THE NUMBER OF
C          SYSTEMS OF LINEAR EQUATIONS.
C
C        W CONTAINS THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER.
C          IF THE SUBROUTINE IS BEING USED TO SOLVE SYSTEMS OF LINEAR
C          EQUATIONS (A-W(R)*I)*X(R)=B(R), WHERE I IS THE IDENTITY
C          MATRIX, W(R) SHOULD BE SET ACCORDINGLY, FOR R=1,2,...,M.
C
C        Z CONTAINS THE CONSTANT MATRIX COLUMNS (B(R),R=1,2,...,M), IF
C          THE SUBROUTINE IS USED TO SOLVE SYSTEMS OF LINEAR EQUATIONS.
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER RV
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT.
C
C     ON OUTPUT
C
C        A AND W ARE UNALTERED.
C
C        Z CONTAINS THE ASSOCIATED SET OF ORTHOGONAL EIGENVECTORS.
C          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO.  IF THE
C          SUBROUTINE IS USED TO SOLVE SYSTEMS OF LINEAR EQUATIONS,
C          Z CONTAINS THE SOLUTION MATRIX COLUMNS (X(R),R=1,2,...,M).
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          -R         IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
C                     EIGENVALUE FAILS TO CONVERGE, OR IF THE R-TH
C                     SYSTEM OF LINEAR EQUATIONS IS NEARLY SINGULAR.
C
C        RV AND RV6 ARE TEMPORARY STORAGE ARRAYS.  NOTE THAT RV IS
C          OF DIMENSION AT LEAST N*(2*MB-1).  IF THE SUBROUTINE
C          IS BEING USED TO SOLVE SYSTEMS OF LINEAR EQUATIONS, THE
C          DETERMINANT (UP TO SIGN) OF A-W(M)*I IS AVAILABLE, UPON
C          RETURN, AS THE PRODUCT OF THE FIRST N ELEMENTS OF RV.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (M .EQ. 0) GO TO 1001
      MB = MBW
      IF (E21 .LT. 0.0E0) MB = (MBW + 1) / 2
      M1 = MB - 1
      M21 = M1 + MB
      ORDER = 1.0E0 - ABS(E21)
C     .......... FIND VECTORS BY INVERSE ITERATION ..........
      DO 920 R = 1, M
         ITS = 1
         X1 = W(R)
         IF (R .NE. 1) GO TO 100
C     .......... COMPUTE NORM OF MATRIX ..........
         NORM = 0.0E0
C
         DO 60 J = 1, MB
            JJ = MB + 1 - J
            KJ = JJ + M1
            IJ = 1
            V = 0.0E0
C
            DO 40 I = JJ, N
               V = V + ABS(A(I,J))
               IF (E21 .GE. 0.0E0) GO TO 40
               V = V + ABS(A(IJ,KJ))
               IJ = IJ + 1
   40       CONTINUE
C
            NORM = AMAX1(NORM,V)
   60    CONTINUE
C
         IF (E21 .LT. 0.0E0) NORM = 0.5E0 * NORM
C     .......... EPS2 IS THE CRITERION FOR GROUPING,
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                ROOTS ARE MODIFIED BY EPS3,
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ..........
         IF (NORM .EQ. 0.0E0) NORM = 1.0E0
         EPS2 = 1.0E-3 * NORM * ABS(ORDER)
         EPS3 = EPSLON(NORM)
         UK = N
         UK = SQRT(UK)
         EPS4 = UK * EPS3
   80    GROUP = 0
         GO TO 120
C     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
  100    IF (ABS(X1-X0) .GE. EPS2) GO TO 80
         GROUP = GROUP + 1
         IF (ORDER * (X1 - X0) .LE. 0.0E0) X1 = X0 + ORDER * EPS3
C     .......... EXPAND MATRIX, SUBTRACT EIGENVALUE,
C                AND INITIALIZE VECTOR ..........
  120    DO 200 I = 1, N
            IJ = I + MIN0(0,I-M1) * N
            KJ = IJ + MB * N
            IJ1 = KJ + M1 * N
            IF (M1 .EQ. 0) GO TO 180
C
            DO 150 J = 1, M1
               IF (IJ .GT. M1) GO TO 125
               IF (IJ .GT. 0) GO TO 130
               RV(IJ1) = 0.0E0
               IJ1 = IJ1 + N
               GO TO 130
  125          RV(IJ) = A(I,J)
  130          IJ = IJ + N
               II = I + J
               IF (II .GT. N) GO TO 150
               JJ = MB - J
               IF (E21 .GE. 0.0E0) GO TO 140
               II = I
               JJ = MB + J
  140          RV(KJ) = A(II,JJ)
               KJ = KJ + N
  150       CONTINUE
C
  180       RV(IJ) = A(I,MB) - X1
            RV6(I) = EPS4
            IF (ORDER .EQ. 0.0E0) RV6(I) = Z(I,R)
  200    CONTINUE
C
         IF (M1 .EQ. 0) GO TO 600
C     .......... ELIMINATION WITH INTERCHANGES ..........
         DO 580 I = 1, N
            II = I + 1
            MAXK = MIN0(I+M1-1,N)
            MAXJ = MIN0(N-I,M21-2) * N
C
            DO 360 K = I, MAXK
               KJ1 = K
               J = KJ1 + N
               JJ = J + MAXJ
C
               DO 340 KJ = J, JJ, N
                  RV(KJ1) = RV(KJ)
                  KJ1 = KJ
  340          CONTINUE
C
               RV(KJ1) = 0.0E0
  360       CONTINUE
C
            IF (I .EQ. N) GO TO 580
            U = 0.0E0
            MAXK = MIN0(I+M1,N)
            MAXJ = MIN0(N-II,M21-2) * N
C
            DO 450 J = I, MAXK
               IF (ABS(RV(J)) .LT. ABS(U)) GO TO 450
               U = RV(J)
               K = J
  450       CONTINUE
C
            J = I + N
            JJ = J + MAXJ
            IF (K .EQ. I) GO TO 520
            KJ = K
C
            DO 500 IJ = I, JJ, N
               V = RV(IJ)
               RV(IJ) = RV(KJ)
               RV(KJ) = V
               KJ = KJ + N
  500       CONTINUE
C
            IF (ORDER .NE. 0.0E0) GO TO 520
            V = RV6(I)
            RV6(I) = RV6(K)
            RV6(K) = V
  520       IF (U .EQ. 0.0E0) GO TO 580
C
            DO 560 K = II, MAXK
               V = RV(K) / U
               KJ = K
C
               DO 540 IJ = J, JJ, N
                  KJ = KJ + N
                  RV(KJ) = RV(KJ) - V * RV(IJ)
  540          CONTINUE
C
               IF (ORDER .EQ. 0.0E0) RV6(K) = RV6(K) - V * RV6(I)
  560       CONTINUE
C
  580    CONTINUE
C     .......... BACK SUBSTITUTION
C                FOR I=N STEP -1 UNTIL 1 DO -- ..........
  600    DO 630 II = 1, N
            I = N + 1 - II
            MAXJ = MIN0(II,M21)
            IF (MAXJ .EQ. 1) GO TO 620
            IJ1 = I
            J = IJ1 + N
            JJ = J + (MAXJ - 2) * N
C
            DO 610 IJ = J, JJ, N
               IJ1 = IJ1 + 1
               RV6(I) = RV6(I) - RV(IJ) * RV6(IJ1)
  610       CONTINUE
C
  620       V = RV(I)
            IF (ABS(V) .GE. EPS3) GO TO 625
C     .......... SET ERROR -- NEARLY SINGULAR LINEAR SYSTEM ..........
            IF (ORDER .EQ. 0.0E0) IERR = -R
            V = SIGN(EPS3,V)
  625       RV6(I) = RV6(I) / V
  630    CONTINUE
C
         XU = 1.0E0
         IF (ORDER .EQ. 0.0E0) GO TO 870
C     .......... ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                MEMBERS OF GROUP ..........
         IF (GROUP .EQ. 0) GO TO 700
C
         DO 680 JJ = 1, GROUP
            J = R - GROUP - 1 + JJ
            XU = 0.0E0
C
            DO 640 I = 1, N
  640       XU = XU + RV6(I) * Z(I,J)
C
            DO 660 I = 1, N
  660       RV6(I) = RV6(I) - XU * Z(I,J)
C
  680    CONTINUE
C
  700    NORM = 0.0E0
C
         DO 720 I = 1, N
  720    NORM = NORM + ABS(RV6(I))
C
         IF (NORM .GE. 0.1E0) GO TO 840
C     .......... IN-LINE PROCEDURE FOR CHOOSING
C                A NEW STARTING VECTOR ..........
         IF (ITS .GE. N) GO TO 830
         ITS = ITS + 1
         XU = EPS4 / (UK + 1.0E0)
         RV6(1) = EPS4
C
         DO 760 I = 2, N
  760    RV6(I) = XU
C
         RV6(ITS) = RV6(ITS) - EPS4 * UK
         GO TO 600
C     .......... SET ERROR -- NON-CONVERGED EIGENVECTOR ..........
  830    IERR = -R
         XU = 0.0E0
         GO TO 870
C     .......... NORMALIZE SO THAT SUM OF SQUARES IS
C                1 AND EXPAND TO FULL ORDER ..........
  840    U = 0.0E0
C
         DO 860 I = 1, N
  860    U = PYTHAG(U,RV6(I))
C
         XU = 1.0E0 / U
C
  870    DO 900 I = 1, N
  900    Z(I,R) = RV6(I) * XU
C
         X0 = X1
  920 CONTINUE
C
 1001 RETURN
      END
      SUBROUTINE BISECT(N,EPS1,D,E,E2,LB,UB,MM,M,W,IND,IERR,RV4,RV5)
C
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,MM,M1,M2,TAG,IERR,ISTURM
      REAL D(N),E(N),E2(N),W(MM),RV4(N),RV5(N)
      REAL U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,TST1,TST2,EPSLON
      INTEGER IND(MM)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE BISECTION TECHNIQUE
C     IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVALUES OF A TRIDIAGONAL
C     SYMMETRIC MATRIX WHICH LIE IN A SPECIFIED INTERVAL,
C     USING BISECTION.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        EPS1 IS AN ABSOLUTE ERROR TOLERANCE FOR THE COMPUTED
C          EIGENVALUES.  IF THE INPUT EPS1 IS NON-POSITIVE,
C          IT IS RESET FOR EACH SUBMATRIX TO A DEFAULT VALUE,
C          NAMELY, MINUS THE PRODUCT OF THE RELATIVE MACHINE
C          PRECISION AND THE 1-NORM OF THE SUBMATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY.
C
C        LB AND UB DEFINE THE INTERVAL TO BE SEARCHED FOR EIGENVALUES.
C          IF LB IS NOT LESS THAN UB, NO EIGENVALUES WILL BE FOUND.
C
C        MM SHOULD BE SET TO AN UPPER BOUND FOR THE NUMBER OF
C          EIGENVALUES IN THE INTERVAL.  WARNING. IF MORE THAN
C          MM EIGENVALUES ARE DETERMINED TO LIE IN THE INTERVAL,
C          AN ERROR RETURN IS MADE WITH NO EIGENVALUES FOUND.
C
C     ON OUTPUT
C
C        EPS1 IS UNALTERED UNLESS IT HAS BEEN RESET TO ITS
C          (LAST) DEFAULT VALUE.
C
C        D AND E ARE UNALTERED.
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS ALSO SET TO ZERO.
C
C        M IS THE NUMBER OF EIGENVALUES DETERMINED TO LIE IN (LB,UB).
C
C        W CONTAINS THE M EIGENVALUES IN ASCENDING ORDER.
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC..
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          3*N+1      IF M EXCEEDS MM.
C
C        RV4 AND RV5 ARE TEMPORARY STORAGE ARRAYS.
C
C     THE ALGOL PROCEDURE STURMCNT CONTAINED IN TRISTURM
C     APPEARS IN BISECT IN-LINE.
C
C     NOTE THAT SUBROUTINE TQL1 OR IMTQL1 IS GENERALLY FASTER THAN
C     BISECT, IF MORE THAN N/4 EIGENVALUES ARE TO BE FOUND.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      TAG = 0
      T1 = LB
      T2 = UB
C     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES ..........
      DO 40 I = 1, N
         IF (I .EQ. 1) GO TO 20
         TST1 = ABS(D(I)) + ABS(D(I-1))
         TST2 = TST1 + ABS(E(I))
         IF (TST2 .GT. TST1) GO TO 40
   20    E2(I) = 0.0E0
   40 CONTINUE
C     .......... DETERMINE THE NUMBER OF EIGENVALUES
C                IN THE INTERVAL ..........
      P = 1
      Q = N
      X1 = UB
      ISTURM = 1
      GO TO 320
   60 M = S
      X1 = LB
      ISTURM = 2
      GO TO 320
   80 M = M - S
      IF (M .GT. MM) GO TO 980
      Q = 0
      R = 0
C     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
C                INTERVAL BY THE GERSCHGORIN BOUNDS ..........
  100 IF (R .EQ. M) GO TO 1001
      TAG = TAG + 1
      P = Q + 1
      XU = D(P)
      X0 = D(P)
      U = 0.0E0
C
      DO 120 Q = P, N
         X1 = U
         U = 0.0E0
         V = 0.0E0
         IF (Q .EQ. N) GO TO 110
         U = ABS(E(Q+1))
         V = E2(Q+1)
  110    XU = AMIN1(D(Q)-(X1+U),XU)
         X0 = AMAX1(D(Q)+(X1+U),X0)
         IF (V .EQ. 0.0E0) GO TO 140
  120 CONTINUE
C
  140 X1 = EPSLON(AMAX1(ABS(XU),ABS(X0)))
      IF (EPS1 .LE. 0.0E0) EPS1 = -X1
      IF (P .NE. Q) GO TO 180
C     .......... CHECK FOR ISOLATED ROOT WITHIN INTERVAL ..........
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 900
  180 X1 = X1 * (Q - P + 1)
      LB = AMAX1(T1,XU-X1)
      UB = AMIN1(T2,X0+X1)
      X1 = LB
      ISTURM = 3
      GO TO 320
  200 M1 = S + 1
      X1 = UB
      ISTURM = 4
      GO TO 320
  220 M2 = S
      IF (M1 .GT. M2) GO TO 940
C     .......... FIND ROOTS BY BISECTION ..........
      X0 = UB
      ISTURM = 5
C
      DO 240 I = M1, M2
         RV5(I) = UB
         RV4(I) = LB
  240 CONTINUE
C     .......... LOOP FOR K-TH EIGENVALUE
C                FOR K=M2 STEP -1 UNTIL M1 DO --
C                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) ..........
      K = M2
  250    XU = LB
C     .......... FOR I=K STEP -1 UNTIL M1 DO -- ..........
         DO 260 II = M1, K
            I = M1 + K - II
            IF (XU .GE. RV4(I)) GO TO 260
            XU = RV4(I)
            GO TO 280
  260    CONTINUE
C
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K)
C     .......... NEXT BISECTION STEP ..........
  300    X1 = (XU + X0) * 0.5E0
         IF ((X0 - XU) .LE. ABS(EPS1)) GO TO 420
         TST1 = 2.0E0 * (ABS(XU) + ABS(X0))
         TST2 = TST1 + (X0 - XU)
         IF (TST2 .EQ. TST1) GO TO 420
C     .......... IN-LINE PROCEDURE FOR STURM SEQUENCE ..........
  320    S = P - 1
         U = 1.0E0
C
         DO 340 I = P, Q
            IF (U .NE. 0.0E0) GO TO 325
            V = ABS(E(I)) / EPSLON(1.0E0)
            IF (E2(I) .EQ. 0.0E0) V = 0.0E0
            GO TO 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
            IF (U .LT. 0.0E0) S = S + 1
  340    CONTINUE
C
         GO TO (60,80,200,220,360), ISTURM
C     .......... REFINE INTERVALS ..........
  360    IF (S .GE. K) GO TO 400
         XU = X1
         IF (S .GE. M1) GO TO 380
         RV4(M1) = X1
         GO TO 300
  380    RV4(S+1) = X1
         IF (RV5(S) .GT. X1) RV5(S) = X1
         GO TO 300
  400    X0 = X1
         GO TO 300
C     .......... K-TH EIGENVALUE FOUND ..........
  420    RV5(K) = X1
      K = K - 1
      IF (K .GE. M1) GO TO 250
C     .......... ORDER EIGENVALUES TAGGED WITH THEIR
C                SUBMATRIX ASSOCIATIONS ..........
  900 S = R
      R = R + M2 - M1 + 1
      J = 1
      K = M1
C
      DO 920 L = 1, R
         IF (J .GT. S) GO TO 910
         IF (K .GT. M2) GO TO 940
         IF (RV5(K) .GE. W(L)) GO TO 915
C
         DO 905 II = J, S
            I = L + S - II
            W(I+1) = W(I)
            IND(I+1) = IND(I)
  905    CONTINUE
C
  910    W(L) = RV5(K)
         IND(L) = TAG
         K = K + 1
         GO TO 920
  915    J = J + 1
  920 CONTINUE
C
  940 IF (Q .LT. N) GO TO 100
      GO TO 1001
C     .......... SET ERROR -- UNDERESTIMATE OF NUMBER OF
C                EIGENVALUES IN INTERVAL ..........
  980 IERR = 3 * N + 1
 1001 LB = T1
      UB = T2
      RETURN
      END
      SUBROUTINE BQR(NM,N,MB,A,T,R,IERR,NV,RV)
C
      INTEGER I,J,K,L,M,N,II,IK,JK,JM,KJ,KK,KM,LL,MB,MK,MN,MZ,
     X        M1,M2,M3,M4,NI,NM,NV,ITS,KJ1,M21,M31,IERR,IMULT
      REAL A(NM,MB),RV(NV)
      REAL F,G,Q,R,S,T,TST1,TST2,SCALE,PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BQR,
C     NUM. MATH. 16, 85-92(1970) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 266-272(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUE OF SMALLEST (USUALLY)
C     MAGNITUDE OF A REAL SYMMETRIC BAND MATRIX USING THE
C     QR ALGORITHM WITH SHIFTS OF ORIGIN.  CONSECUTIVE CALLS
C     CAN BE MADE TO FIND FURTHER EIGENVALUES.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        MB IS THE (HALF) BAND WIDTH OF THE MATRIX, DEFINED AS THE
C          NUMBER OF ADJACENT DIAGONALS, INCLUDING THE PRINCIPAL
C          DIAGONAL, REQUIRED TO SPECIFY THE NON-ZERO PORTION OF THE
C          LOWER TRIANGLE OF THE MATRIX.
C
C        A CONTAINS THE LOWER TRIANGLE OF THE SYMMETRIC BAND INPUT
C          MATRIX STORED AS AN N BY MB ARRAY.  ITS LOWEST SUBDIAGONAL
C          IS STORED IN THE LAST N+1-MB POSITIONS OF THE FIRST COLUMN,
C          ITS NEXT SUBDIAGONAL IN THE LAST N+2-MB POSITIONS OF THE
C          SECOND COLUMN, FURTHER SUBDIAGONALS SIMILARLY, AND FINALLY
C          ITS PRINCIPAL DIAGONAL IN THE N POSITIONS OF THE LAST COLUMN.
C          CONTENTS OF STORAGES NOT PART OF THE MATRIX ARE ARBITRARY.
C          ON A SUBSEQUENT CALL, ITS OUTPUT CONTENTS FROM THE PREVIOUS
C          CALL SHOULD BE PASSED.
C
C        T SPECIFIES THE SHIFT (OF EIGENVALUES) APPLIED TO THE DIAGONAL
C          OF A IN FORMING THE INPUT MATRIX. WHAT IS ACTUALLY DETERMINED
C          IS THE EIGENVALUE OF A+TI (I IS THE IDENTITY MATRIX) NEAREST
C          TO T.  ON A SUBSEQUENT CALL, THE OUTPUT VALUE OF T FROM THE
C          PREVIOUS CALL SHOULD BE PASSED IF THE NEXT NEAREST EIGENVALUE
C          IS SOUGHT.
C
C        R SHOULD BE SPECIFIED AS ZERO ON THE FIRST CALL, AND AS ITS
C          OUTPUT VALUE FROM THE PREVIOUS CALL ON A SUBSEQUENT CALL.
C          IT IS USED TO DETERMINE WHEN THE LAST ROW AND COLUMN OF
C          THE TRANSFORMED BAND MATRIX CAN BE REGARDED AS NEGLIGIBLE.
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER RV
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT.
C
C     ON OUTPUT
C
C        A CONTAINS THE TRANSFORMED BAND MATRIX.  THE MATRIX A+TI
C          DERIVED FROM THE OUTPUT PARAMETERS IS SIMILAR TO THE
C          INPUT A+TI TO WITHIN ROUNDING ERRORS.  ITS LAST ROW AND
C          COLUMN ARE NULL (IF IERR IS ZERO).
C
C        T CONTAINS THE COMPUTED EIGENVALUE OF A+TI (IF IERR IS ZERO).
C
C        R CONTAINS THE MAXIMUM OF ITS INPUT VALUE AND THE NORM OF THE
C          LAST COLUMN OF THE INPUT MATRIX A.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          N          IF THE EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C        RV IS A TEMPORARY STORAGE ARRAY OF DIMENSION AT LEAST
C          (2*MB**2+4*MB-3).  THE FIRST (3*MB-2) LOCATIONS CORRESPOND
C          TO THE ALGOL ARRAY B, THE NEXT (2*MB-1) LOCATIONS CORRESPOND
C          TO THE ALGOL ARRAY H, AND THE FINAL (2*MB**2-MB) LOCATIONS
C          CORRESPOND TO THE MB BY (2*MB-1) ALGOL ARRAY U.
C
C     NOTE. FOR A SUBSEQUENT CALL, N SHOULD BE REPLACED BY N-1, BUT
C     MB SHOULD NOT BE ALTERED EVEN WHEN IT EXCEEDS THE CURRENT N.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      M1 = MIN0(MB,N)
      M = M1 - 1
      M2 = M + M
      M21 = M2 + 1
      M3 = M21 + M
      M31 = M3 + 1
      M4 = M31 + M2
      MN = M + N
      MZ = MB - M1
      ITS = 0
C     .......... TEST FOR CONVERGENCE ..........
   40 G = A(N,MB)
      IF (M .EQ. 0) GO TO 360
      F = 0.0E0
C
      DO 50 K = 1, M
         MK = K + MZ
         F = F + ABS(A(N,MK))
   50 CONTINUE
C
      IF (ITS .EQ. 0 .AND. F .GT. R) R = F
      TST1 = R
      TST2 = TST1 + F
      IF (TST2 .LE. TST1) GO TO 360
      IF (ITS .EQ. 30) GO TO 1000
      ITS = ITS + 1
C     .......... FORM SHIFT FROM BOTTOM 2 BY 2 MINOR ..........
      IF (F .GT. 0.25E0 * R .AND. ITS .LT. 5) GO TO 90
      F = A(N,MB-1)
      IF (F .EQ. 0.0E0) GO TO 70
      Q = (A(N-1,MB) - G) / (2.0E0 * F)
      S = PYTHAG(Q,1.0E0)
      G = G - F / (Q + SIGN(S,Q))
   70 T = T + G
C
      DO 80 I = 1, N
   80 A(I,MB) = A(I,MB) - G
C
   90 DO 100 K = M31, M4
  100 RV(K) = 0.0E0
C
      DO 350 II = 1, MN
         I = II - M
         NI = N - II
         IF (NI .LT. 0) GO TO 230
C     .......... FORM COLUMN OF SHIFTED MATRIX A-G*I ..........
         L = MAX0(1,2-I)
C
         DO 110 K = 1, M3
  110    RV(K) = 0.0E0
C
         DO 120 K = L, M1
            KM = K + M
            MK = K + MZ
            RV(KM) = A(II,MK)
  120    CONTINUE
C
         LL = MIN0(M,NI)
         IF (LL .EQ. 0) GO TO 135
C
         DO 130 K = 1, LL
            KM = K + M21
            IK = II + K
            MK = MB - K
            RV(KM) = A(IK,MK)
  130    CONTINUE
C     .......... PRE-MULTIPLY WITH HOUSEHOLDER REFLECTIONS ..........
  135    LL = M2
         IMULT = 0
C     .......... MULTIPLICATION PROCEDURE ..........
  140    KJ = M4 - M1
C
         DO 170 J = 1, LL
            KJ = KJ + M1
            JM = J + M3
            IF (RV(JM) .EQ. 0.0E0) GO TO 170
            F = 0.0E0
C
            DO 150 K = 1, M1
               KJ = KJ + 1
               JK = J + K - 1
               F = F + RV(KJ) * RV(JK)
  150       CONTINUE
C
            F = F / RV(JM)
            KJ = KJ - M1
C
            DO 160 K = 1, M1
               KJ = KJ + 1
               JK = J + K - 1
               RV(JK) = RV(JK) - RV(KJ) * F
  160       CONTINUE
C
            KJ = KJ - M1
  170    CONTINUE
C
         IF (IMULT .NE. 0) GO TO 280
C     .......... HOUSEHOLDER REFLECTION ..........
         F = RV(M21)
         S = 0.0E0
         RV(M4) = 0.0E0
         SCALE = 0.0E0
C
         DO 180 K = M21, M3
  180    SCALE = SCALE + ABS(RV(K))
C
         IF (SCALE .EQ. 0.0E0) GO TO 210
C
         DO 190 K = M21, M3
  190    S = S + (RV(K)/SCALE)**2
C
         S = SCALE * SCALE * S
         G = -SIGN(SQRT(S),F)
         RV(M21) = G
         RV(M4) = S - F * G
         KJ = M4 + M2 * M1 + 1
         RV(KJ) = F - G
C
         DO 200 K = 2, M1
            KJ = KJ + 1
            KM = K + M2
            RV(KJ) = RV(KM)
  200    CONTINUE
C     .......... SAVE COLUMN OF TRIANGULAR FACTOR R ..........
  210    DO 220 K = L, M1
            KM = K + M
            MK = K + MZ
            A(II,MK) = RV(KM)
  220    CONTINUE
C
  230    L = MAX0(1,M1+1-I)
         IF (I .LE. 0) GO TO 300
C     .......... PERFORM ADDITIONAL STEPS ..........
         DO 240 K = 1, M21
  240    RV(K) = 0.0E0
C
         LL = MIN0(M1,NI+M1)
C     .......... GET ROW OF TRIANGULAR FACTOR R ..........
         DO 250 KK = 1, LL
            K = KK - 1
            KM = K + M1
            IK = I + K
            MK = MB - K
            RV(KM) = A(IK,MK)
  250    CONTINUE
C     .......... POST-MULTIPLY WITH HOUSEHOLDER REFLECTIONS ..........
         LL = M1
         IMULT = 1
         GO TO 140
C     .......... STORE COLUMN OF NEW A MATRIX ..........
  280    DO 290 K = L, M1
            MK = K + MZ
            A(I,MK) = RV(K)
  290    CONTINUE
C     .......... UPDATE HOUSEHOLDER REFLECTIONS ..........
  300    IF (L .GT. 1) L = L - 1
         KJ1 = M4 + L * M1
C
         DO 320 J = L, M2
            JM = J + M3
            RV(JM) = RV(JM+1)
C
            DO 320 K = 1, M1
               KJ1 = KJ1 + 1
               KJ = KJ1 - M1
               RV(KJ) = RV(KJ1)
  320    CONTINUE
C
  350 CONTINUE
C
      GO TO 40
C     .......... CONVERGENCE ..........
  360 T = T + G
C
      DO 380 I = 1, N
  380 A(I,MB) = A(I,MB) - G
C
      DO 400 K = 1, M1
         MK = K + MZ
         A(N,MK) = 0.0E0
  400 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = N
 1001 RETURN
      END
      SUBROUTINE CBABK2(NM,N,LOW,IGH,SCALE,M,ZR,ZI)
C
      INTEGER I,J,K,M,N,II,NM,IGH,LOW
      REAL SCALE(N),ZR(NM,M),ZI(NM,M)
      REAL S
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE
C     CBABK2, WHICH IS A COMPLEX VERSION OF BALBAK,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX GENERAL
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     BALANCED MATRIX DETERMINED BY  CBAL.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY  CBAL.
C
C        SCALE CONTAINS INFORMATION DETERMINING THE PERMUTATIONS
C          AND SCALING FACTORS USED BY  CBAL.
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVECTORS TO BE
C          BACK TRANSFORMED IN THEIR FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE TRANSFORMED EIGENVECTORS
C          IN THEIR FIRST M COLUMNS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
      IF (IGH .EQ. LOW) GO TO 120
C
      DO 110 I = LOW, IGH
         S = SCALE(I)
C     .......... LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
C                IF THE FOREGOING STATEMENT IS REPLACED BY
C                S=1.0E0/SCALE(I). ..........
         DO 100 J = 1, M
            ZR(I,J) = ZR(I,J) * S
            ZI(I,J) = ZI(I,J) * S
  100    CONTINUE
C
  110 CONTINUE
C     .......... FOR I=LOW-1 STEP -1 UNTIL 1,
C                IGH+1 STEP 1 UNTIL N DO -- ..........
  120 DO 140 II = 1, N
         I = II
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 140
         IF (I .LT. LOW) I = LOW - II
         K = SCALE(I)
         IF (K .EQ. I) GO TO 140
C
         DO 130 J = 1, M
            S = ZR(I,J)
            ZR(I,J) = ZR(K,J)
            ZR(K,J) = S
            S = ZI(I,J)
            ZI(I,J) = ZI(K,J)
            ZI(K,J) = S
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE CBAL(NM,N,AR,AI,LOW,IGH,SCALE)
C
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      REAL AR(NM,N),AI(NM,N),SCALE(N)
      REAL C,F,G,R,S,B2,RADIX
      LOGICAL NOCONV
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE
C     CBALANCE, WHICH IS A COMPLEX VERSION OF BALANCE,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     THIS SUBROUTINE BALANCES A COMPLEX MATRIX AND ISOLATES
C     EIGENVALUES WHENEVER POSSIBLE.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX MATRIX TO BE BALANCED.
C
C     ON OUTPUT
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE BALANCED MATRIX.
C
C        LOW AND IGH ARE TWO INTEGERS SUCH THAT AR(I,J) AND AI(I,J)
C          ARE EQUAL TO ZERO IF
C           (1) I IS GREATER THAN J AND
C           (2) J=1,...,LOW-1 OR I=IGH+1,...,N.
C
C        SCALE CONTAINS INFORMATION DETERMINING THE
C           PERMUTATIONS AND SCALING FACTORS USED.
C
C     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH
C     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED
C     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS
C     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN
C        SCALE(J) = P(J),    FOR J = 1,...,LOW-1
C                 = D(J,J)       J = LOW,...,IGH
C                 = P(J)         J = IGH+1,...,N.
C     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,
C     THEN 1 TO LOW-1.
C
C     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.
C
C     THE ALGOL PROCEDURE EXC CONTAINED IN CBALANCE APPEARS IN
C     CBAL  IN LINE.  (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS
C     K,L HAVE BEEN REVERSED.)
C
C     ARITHMETIC IS REAL THROUGHOUT.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      RADIX = 16.0E0
C
      B2 = RADIX * RADIX
      K = 1
      L = N
      GO TO 100
C     .......... IN-LINE PROCEDURE FOR ROW AND
C                COLUMN EXCHANGE ..........
   20 SCALE(M) = J
      IF (J .EQ. M) GO TO 50
C
      DO 30 I = 1, L
         F = AR(I,J)
         AR(I,J) = AR(I,M)
         AR(I,M) = F
         F = AI(I,J)
         AI(I,J) = AI(I,M)
         AI(I,M) = F
   30 CONTINUE
C
      DO 40 I = K, N
         F = AR(J,I)
         AR(J,I) = AR(M,I)
         AR(M,I) = F
         F = AI(J,I)
         AI(J,I) = AI(M,I)
         AI(M,I) = F
   40 CONTINUE
C
   50 GO TO (80,130), IEXC
C     .......... SEARCH FOR ROWS ISOLATING AN EIGENVALUE
C                AND PUSH THEM DOWN ..........
   80 IF (L .EQ. 1) GO TO 280
      L = L - 1
C     .......... FOR J=L STEP -1 UNTIL 1 DO -- ..........
  100 DO 120 JJ = 1, L
         J = L + 1 - JJ
C
         DO 110 I = 1, L
            IF (I .EQ. J) GO TO 110
            IF (AR(J,I) .NE. 0.0E0 .OR. AI(J,I) .NE. 0.0E0) GO TO 120
  110    CONTINUE
C
         M = L
         IEXC = 1
         GO TO 20
  120 CONTINUE
C
      GO TO 140
C     .......... SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
C                AND PUSH THEM LEFT ..........
  130 K = K + 1
C
  140 DO 170 J = K, L
C
         DO 150 I = K, L
            IF (I .EQ. J) GO TO 150
            IF (AR(I,J) .NE. 0.0E0 .OR. AI(I,J) .NE. 0.0E0) GO TO 170
  150    CONTINUE
C
         M = K
         IEXC = 2
         GO TO 20
  170 CONTINUE
C     .......... NOW BALANCE THE SUBMATRIX IN ROWS K TO L ..........
      DO 180 I = K, L
  180 SCALE(I) = 1.0E0
C     .......... ITERATIVE LOOP FOR NORM REDUCTION ..........
  190 NOCONV = .FALSE.
C
      DO 270 I = K, L
         C = 0.0E0
         R = 0.0E0
C
         DO 200 J = K, L
            IF (J .EQ. I) GO TO 200
            C = C + ABS(AR(J,I)) + ABS(AI(J,I))
            R = R + ABS(AR(I,J)) + ABS(AI(I,J))
  200    CONTINUE
C     .......... GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW ..........
         IF (C .EQ. 0.0E0 .OR. R .EQ. 0.0E0) GO TO 270
         G = R / RADIX
         F = 1.0E0
         S = C + R
  210    IF (C .GE. G) GO TO 220
         F = F * RADIX
         C = C * B2
         GO TO 210
  220    G = R * RADIX
  230    IF (C .LT. G) GO TO 240
         F = F / RADIX
         C = C / B2
         GO TO 230
C     .......... NOW BALANCE ..........
  240    IF ((C + R) / F .GE. 0.95E0 * S) GO TO 270
         G = 1.0E0 / F
         SCALE(I) = SCALE(I) * F
         NOCONV = .TRUE.
C
         DO 250 J = K, N
            AR(I,J) = AR(I,J) * G
            AI(I,J) = AI(I,J) * G
  250    CONTINUE
C
         DO 260 J = 1, L
            AR(J,I) = AR(J,I) * F
            AI(J,I) = AI(J,I) * F
  260    CONTINUE
C
  270 CONTINUE
C
      IF (NOCONV) GO TO 190
C
  280 LOW = K
      IGH = L
      RETURN
      END
      SUBROUTINE CG(NM,N,AR,AI,WR,WI,MATZ,ZR,ZI,FV1,FV2,FV3,IERR)
C
      INTEGER N,NM,IS1,IS2,IERR,MATZ
      REAL AR(NM,N),AI(NM,N),WR(N),WI(N),ZR(NM,N),ZI(NM,N),
     X       FV1(N),FV2(N),FV3(N)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     OF A COMPLEX GENERAL MATRIX.
C
C     ON INPUT
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT.
C
C        N  IS THE ORDER OF THE MATRIX  A=(AR,AI).
C
C        AR  AND  AI  CONTAIN THE REAL AND IMAGINARY PARTS,
C        RESPECTIVELY, OF THE COMPLEX GENERAL MATRIX.
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT
C
C        WR  AND  WI  CONTAIN THE REAL AND IMAGINARY PARTS,
C        RESPECTIVELY, OF THE EIGENVALUES.
C
C        ZR  AND  ZI  CONTAIN THE REAL AND IMAGINARY PARTS,
C        RESPECTIVELY, OF THE EIGENVECTORS IF MATZ IS NOT ZERO.
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR COMQR
C           AND COMQR2.  THE NORMAL COMPLETION CODE IS ZERO.
C
C        FV1, FV2, AND  FV3  ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 CALL  CBAL(NM,N,AR,AI,IS1,IS2,FV1)
      CALL  CORTH(NM,N,IS1,IS2,AR,AI,FV2,FV3)
      IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  COMQR(NM,N,IS1,IS2,AR,AI,WR,WI,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  COMQR2(NM,N,IS1,IS2,FV2,FV3,AR,AI,WR,WI,ZR,ZI,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  CBABK2(NM,N,IS1,IS2,FV1,N,ZR,ZI)
   50 RETURN
      END
      SUBROUTINE CH(NM,N,AR,AI,W,MATZ,ZR,ZI,FV1,FV2,FM1,IERR)
C
      INTEGER I,J,N,NM,IERR,MATZ
      REAL AR(NM,N),AI(NM,N),W(N),ZR(NM,N),ZI(NM,N),
     X       FV1(N),FV2(N),FM1(2,N)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     OF A COMPLEX HERMITIAN MATRIX.
C
C     ON INPUT
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT.
C
C        N  IS THE ORDER OF THE MATRIX  A=(AR,AI).
C
C        AR  AND  AI  CONTAIN THE REAL AND IMAGINARY PARTS,
C        RESPECTIVELY, OF THE COMPLEX HERMITIAN MATRIX.
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT
C
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.
C
C        ZR  AND  ZI  CONTAIN THE REAL AND IMAGINARY PARTS,
C        RESPECTIVELY, OF THE EIGENVECTORS IF MATZ IS NOT ZERO.
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT
C           AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO.
C
C        FV1, FV2, AND  FM1  ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 CALL  HTRIDI(NM,N,AR,AI,W,FV1,FV2,FM1)
      IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 DO 40 I = 1, N
C
         DO 30 J = 1, N
            ZR(J,I) = 0.0E0
   30    CONTINUE
C
         ZR(I,I) = 1.0E0
   40 CONTINUE
C
      CALL  TQL2(NM,N,W,FV1,ZR,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  HTRIBK(NM,N,AR,AI,FM1,N,ZR,ZI)
   50 RETURN
      END
      SUBROUTINE CINVIT(NM,N,AR,AI,WR,WI,SELECT,MM,M,ZR,ZI,
     X                  IERR,RM1,RM2,RV1,RV2)
C
      INTEGER I,J,K,M,N,S,II,MM,MP,NM,UK,IP1,ITS,KM1,IERR
      REAL AR(NM,N),AI(NM,N),WR(N),WI(N),ZR(NM,MM),
     X       ZI(NM,MM),RM1(N,N),RM2(N,N),RV1(N),RV2(N)
      REAL X,Y,EPS3,NORM,NORMV,EPSLON,GROWTO,ILAMBD,PYTHAG,
     X       RLAMBD,UKROOT
      LOGICAL SELECT(N)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE CX INVIT
C     BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP. VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A COMPLEX UPPER
C     HESSENBERG MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
C     USING INVERSE ITERATION.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE HESSENBERG MATRIX.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY,
C          OF THE EIGENVALUES OF THE MATRIX.  THE EIGENVALUES MUST BE
C          STORED IN A MANNER IDENTICAL TO THAT OF SUBROUTINE  COMLR,
C          WHICH RECOGNIZES POSSIBLE SPLITTING OF THE MATRIX.
C
C        SELECT SPECIFIES THE EIGENVECTORS TO BE FOUND.  THE
C          EIGENVECTOR CORRESPONDING TO THE J-TH EIGENVALUE IS
C          SPECIFIED BY SETTING SELECT(J) TO .TRUE..
C
C        MM SHOULD BE SET TO AN UPPER BOUND FOR THE NUMBER OF
C          EIGENVECTORS TO BE FOUND.
C
C     ON OUTPUT
C
C        AR, AI, WI, AND SELECT ARE UNALTERED.
C
C        WR MAY HAVE BEEN ALTERED SINCE CLOSE EIGENVALUES ARE PERTURBED
C          SLIGHTLY IN SEARCHING FOR INDEPENDENT EIGENVECTORS.
C
C        M IS THE NUMBER OF EIGENVECTORS ACTUALLY FOUND.
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY,
C          OF THE EIGENVECTORS.  THE EIGENVECTORS ARE NORMALIZED
C          SO THAT THE COMPONENT OF LARGEST MAGNITUDE IS 1.
C          ANY VECTOR WHICH FAILS THE ACCEPTANCE TEST IS SET TO ZERO.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          -(2*N+1)   IF MORE THAN MM EIGENVECTORS HAVE BEEN SPECIFIED,
C          -K         IF THE ITERATION CORRESPONDING TO THE K-TH
C                     VALUE FAILS,
C          -(N+K)     IF BOTH ERROR SITUATIONS OCCUR.
C
C        RM1, RM2, RV1, AND RV2 ARE TEMPORARY STORAGE ARRAYS.
C
C     THE ALGOL PROCEDURE GUESSVEC APPEARS IN CINVIT IN LINE.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      UK = 0
      S = 1
C
      DO 980 K = 1, N
         IF (.NOT. SELECT(K)) GO TO 980
         IF (S .GT. MM) GO TO 1000
         IF (UK .GE. K) GO TO 200
C     .......... CHECK FOR POSSIBLE SPLITTING ..........
         DO 120 UK = K, N
            IF (UK .EQ. N) GO TO 140
            IF (AR(UK+1,UK) .EQ. 0.0E0 .AND. AI(UK+1,UK) .EQ. 0.0E0)
     X         GO TO 140
  120    CONTINUE
C     .......... COMPUTE INFINITY NORM OF LEADING UK BY UK
C                (HESSENBERG) MATRIX ..........
  140    NORM = 0.0E0
         MP = 1
C
         DO 180 I = 1, UK
            X = 0.0E0
C
            DO 160 J = MP, UK
  160       X = X + PYTHAG(AR(I,J),AI(I,J))
C
            IF (X .GT. NORM) NORM = X
            MP = I
  180    CONTINUE
C     .......... EPS3 REPLACES ZERO PIVOT IN DECOMPOSITION
C                AND CLOSE ROOTS ARE MODIFIED BY EPS3 ..........
         IF (NORM .EQ. 0.0E0) NORM = 1.0E0
         EPS3 = EPSLON(NORM)
C     .......... GROWTO IS THE CRITERION FOR GROWTH ..........
         UKROOT = UK
         UKROOT = SQRT(UKROOT)
         GROWTO = 0.1E0 / UKROOT
  200    RLAMBD = WR(K)
         ILAMBD = WI(K)
         IF (K .EQ. 1) GO TO 280
         KM1 = K - 1
         GO TO 240
C     .......... PERTURB EIGENVALUE IF IT IS CLOSE
C                TO ANY PREVIOUS EIGENVALUE ..........
  220    RLAMBD = RLAMBD + EPS3
C     .......... FOR I=K-1 STEP -1 UNTIL 1 DO -- ..........
  240    DO 260 II = 1, KM1
            I = K - II
            IF (SELECT(I) .AND. ABS(WR(I)-RLAMBD) .LT. EPS3 .AND.
     X         ABS(WI(I)-ILAMBD) .LT. EPS3) GO TO 220
  260    CONTINUE
C
         WR(K) = RLAMBD
C     .......... FORM UPPER HESSENBERG (AR,AI)-(RLAMBD,ILAMBD)*I
C                AND INITIAL COMPLEX VECTOR ..........
  280    MP = 1
C
         DO 320 I = 1, UK
C
            DO 300 J = MP, UK
               RM1(I,J) = AR(I,J)
               RM2(I,J) = AI(I,J)
  300       CONTINUE
C
            RM1(I,I) = RM1(I,I) - RLAMBD
            RM2(I,I) = RM2(I,I) - ILAMBD
            MP = I
            RV1(I) = EPS3
  320    CONTINUE
C     .......... TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
C                REPLACING ZERO PIVOTS BY EPS3 ..........
         IF (UK .EQ. 1) GO TO 420
C
         DO 400 I = 2, UK
            MP = I - 1
            IF (PYTHAG(RM1(I,MP),RM2(I,MP)) .LE.
     X          PYTHAG(RM1(MP,MP),RM2(MP,MP))) GO TO 360
C
            DO 340 J = MP, UK
               Y = RM1(I,J)
               RM1(I,J) = RM1(MP,J)
               RM1(MP,J) = Y
               Y = RM2(I,J)
               RM2(I,J) = RM2(MP,J)
               RM2(MP,J) = Y
  340       CONTINUE
C
  360       IF (RM1(MP,MP) .EQ. 0.0E0 .AND. RM2(MP,MP) .EQ. 0.0E0)
     X         RM1(MP,MP) = EPS3
            CALL CDIV(RM1(I,MP),RM2(I,MP),RM1(MP,MP),RM2(MP,MP),X,Y)
            IF (X .EQ. 0.0E0 .AND. Y .EQ. 0.0E0) GO TO 400
C
            DO 380 J = I, UK
               RM1(I,J) = RM1(I,J) - X * RM1(MP,J) + Y * RM2(MP,J)
               RM2(I,J) = RM2(I,J) - X * RM2(MP,J) - Y * RM1(MP,J)
  380       CONTINUE
C
  400    CONTINUE
C
  420    IF (RM1(UK,UK) .EQ. 0.0E0 .AND. RM2(UK,UK) .EQ. 0.0E0)
     X      RM1(UK,UK) = EPS3
         ITS = 0
C     .......... BACK SUBSTITUTION
C                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
  660    DO 720 II = 1, UK
            I = UK + 1 - II
            X = RV1(I)
            Y = 0.0E0
            IF (I .EQ. UK) GO TO 700
            IP1 = I + 1
C
            DO 680 J = IP1, UK
               X = X - RM1(I,J) * RV1(J) + RM2(I,J) * RV2(J)
               Y = Y - RM1(I,J) * RV2(J) - RM2(I,J) * RV1(J)
  680       CONTINUE
C
  700       CALL CDIV(X,Y,RM1(I,I),RM2(I,I),RV1(I),RV2(I))
  720    CONTINUE
C     .......... ACCEPTANCE TEST FOR EIGENVECTOR
C                AND NORMALIZATION ..........
         ITS = ITS + 1
         NORM = 0.0E0
         NORMV = 0.0E0
C
         DO 780 I = 1, UK
            X = PYTHAG(RV1(I),RV2(I))
            IF (NORMV .GE. X) GO TO 760
            NORMV = X
            J = I
  760       NORM = NORM + X
  780    CONTINUE
C
         IF (NORM .LT. GROWTO) GO TO 840
C     .......... ACCEPT VECTOR ..........
         X = RV1(J)
         Y = RV2(J)
C
         DO 820 I = 1, UK
            CALL CDIV(RV1(I),RV2(I),X,Y,ZR(I,S),ZI(I,S))
  820    CONTINUE
C
         IF (UK .EQ. N) GO TO 940
         J = UK + 1
         GO TO 900
C     .......... IN-LINE PROCEDURE FOR CHOOSING
C                A NEW STARTING VECTOR ..........
  840    IF (ITS .GE. UK) GO TO 880
         X = UKROOT
         Y = EPS3 / (X + 1.0E0)
         RV1(1) = EPS3
C
         DO 860 I = 2, UK
  860    RV1(I) = Y
C
         J = UK - ITS + 1
         RV1(J) = RV1(J) - EPS3 * X
         GO TO 660
C     .......... SET ERROR -- UNACCEPTED EIGENVECTOR ..........
  880    J = 1
         IERR = -K
C     .......... SET REMAINING VECTOR COMPONENTS TO ZERO ..........
  900    DO 920 I = J, N
            ZR(I,S) = 0.0E0
            ZI(I,S) = 0.0E0
  920    CONTINUE
C
  940    S = S + 1
  980 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- UNDERESTIMATE OF EIGENVECTOR
C                SPACE REQUIRED ..........
 1000 IF (IERR .NE. 0) IERR = IERR - N
      IF (IERR .EQ. 0) IERR = -(2 * N + 1)
 1001 M = S - 1
      RETURN
      END
      SUBROUTINE COMBAK(NM,LOW,IGH,AR,AI,INT,M,ZR,ZI)
C
      INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
      REAL AR(NM,IGH),AI(NM,IGH),ZR(NM,M),ZI(NM,M)
      REAL XR,XI
      INTEGER INT(IGH)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE COMBAK,
C     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX GENERAL
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     UPPER HESSENBERG MATRIX DETERMINED BY  COMHES.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1 AND IGH EQUAL TO THE ORDER OF THE MATRIX.
C
C        AR AND AI CONTAIN THE MULTIPLIERS WHICH WERE USED IN THE
C          REDUCTION BY  COMHES  IN THEIR LOWER TRIANGLES
C          BELOW THE SUBDIAGONAL.
C
C        INT CONTAINS INFORMATION ON THE ROWS AND COLUMNS
C          INTERCHANGED IN THE REDUCTION BY  COMHES.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVECTORS TO BE
C          BACK TRANSFORMED IN THEIR FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE TRANSFORMED EIGENVECTORS
C          IN THEIR FIRST M COLUMNS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 140 MM = KP1, LA
         MP = LOW + IGH - MM
         MP1 = MP + 1
C
         DO 110 I = MP1, IGH
            XR = AR(I,MP-1)
            XI = AI(I,MP-1)
            IF (XR .EQ. 0.0E0 .AND. XI .EQ. 0.0E0) GO TO 110
C
            DO 100 J = 1, M
               ZR(I,J) = ZR(I,J) + XR * ZR(MP,J) - XI * ZI(MP,J)
               ZI(I,J) = ZI(I,J) + XR * ZI(MP,J) + XI * ZR(MP,J)
  100       CONTINUE
C
  110    CONTINUE
C
         I = INT(MP)
         IF (I .EQ. MP) GO TO 140
C
         DO 130 J = 1, M
            XR = ZR(I,J)
            ZR(I,J) = ZR(MP,J)
            ZR(MP,J) = XR
            XI = ZI(I,J)
            ZI(I,J) = ZI(MP,J)
            ZI(MP,J) = XI
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE COMHES(NM,N,LOW,IGH,AR,AI,INT)
C
      INTEGER I,J,M,N,LA,NM,IGH,KP1,LOW,MM1,MP1
      REAL AR(NM,N),AI(NM,N)
      REAL XR,XI,YR,YI
      INTEGER INT(IGH)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE COMHES,
C     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     GIVEN A COMPLEX GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     STABILIZED ELEMENTARY SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX INPUT MATRIX.
C
C     ON OUTPUT
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE HESSENBERG MATRIX.  THE
C          MULTIPLIERS WHICH WERE USED IN THE REDUCTION
C          ARE STORED IN THE REMAINING TRIANGLES UNDER THE
C          HESSENBERG MATRIX.
C
C        INT CONTAINS INFORMATION ON THE ROWS AND COLUMNS
C          INTERCHANGED IN THE REDUCTION.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C
      DO 180 M = KP1, LA
         MM1 = M - 1
         XR = 0.0E0
         XI = 0.0E0
         I = M
C
         DO 100 J = M, IGH
            IF (ABS(AR(J,MM1)) + ABS(AI(J,MM1))
     X         .LE. ABS(XR) + ABS(XI)) GO TO 100
            XR = AR(J,MM1)
            XI = AI(J,MM1)
            I = J
  100    CONTINUE
C
         INT(M) = I
         IF (I .EQ. M) GO TO 130
C     .......... INTERCHANGE ROWS AND COLUMNS OF AR AND AI ..........
         DO 110 J = MM1, N
            YR = AR(I,J)
            AR(I,J) = AR(M,J)
            AR(M,J) = YR
            YI = AI(I,J)
            AI(I,J) = AI(M,J)
            AI(M,J) = YI
  110    CONTINUE
C
         DO 120 J = 1, IGH
            YR = AR(J,I)
            AR(J,I) = AR(J,M)
            AR(J,M) = YR
            YI = AI(J,I)
            AI(J,I) = AI(J,M)
            AI(J,M) = YI
  120    CONTINUE
C     .......... END INTERCHANGE ..........
  130    IF (XR .EQ. 0.0E0 .AND. XI .EQ. 0.0E0) GO TO 180
         MP1 = M + 1
C
         DO 160 I = MP1, IGH
            YR = AR(I,MM1)
            YI = AI(I,MM1)
            IF (YR .EQ. 0.0E0 .AND. YI .EQ. 0.0E0) GO TO 160
            CALL CDIV(YR,YI,XR,XI,YR,YI)
            AR(I,MM1) = YR
            AI(I,MM1) = YI
C
            DO 140 J = M, N
               AR(I,J) = AR(I,J) - YR * AR(M,J) + YI * AI(M,J)
               AI(I,J) = AI(I,J) - YR * AI(M,J) - YI * AR(M,J)
  140       CONTINUE
C
            DO 150 J = 1, IGH
               AR(J,M) = AR(J,M) + YR * AR(J,I) - YI * AI(J,I)
               AI(J,M) = AI(J,M) + YR * AI(J,I) + YI * AR(J,I)
  150       CONTINUE
C
  160    CONTINUE
C
  180 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE COMLR(NM,N,LOW,IGH,HR,HI,WR,WI,IERR)
C
      INTEGER I,J,L,M,N,EN,LL,MM,NM,IGH,IM1,ITN,ITS,LOW,MP1,ENM1,IERR
      REAL HR(NM,N),HI(NM,N),WR(N),WI(N)
      REAL SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,TST1,TST2
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE COMLR,
C     NUM. MATH. 12, 369-376(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 396-403(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A COMPLEX
C     UPPER HESSENBERG MATRIX BY THE MODIFIED LR METHOD.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        HR AND HI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX.
C          THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAIN THE
C          MULTIPLIERS WHICH WERE USED IN THE REDUCTION BY  COMHES,
C          IF PERFORMED.
C
C     ON OUTPUT
C
C        THE UPPER HESSENBERG PORTIONS OF HR AND HI HAVE BEEN
C          DESTROYED.  THEREFORE, THEY MUST BE SAVED BEFORE
C          CALLING  COMLR  IF SUBSEQUENT CALCULATION OF
C          EIGENVECTORS IS TO BE PERFORMED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  IF AN ERROR
C          EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C     CALLS CSROOT FOR COMPLEX SQUARE ROOT.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
C     .......... STORE ROOTS ISOLATED BY CBAL ..........
      DO 200 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 200
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  200 CONTINUE
C
      EN = IGH
      TR = 0.0E0
      TI = 0.0E0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUE ..........
  220 IF (EN .LT. LOW) GO TO 1001
      ITS = 0
      ENM1 = EN - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW E0 -- ..........
  240 DO 260 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 300
         TST1 = ABS(HR(L-1,L-1)) + ABS(HI(L-1,L-1))
     X            + ABS(HR(L,L)) + ABS(HI(L,L))
         TST2 = TST1 + ABS(HR(L,L-1)) + ABS(HI(L,L-1))
         IF (TST2 .EQ. TST1) GO TO 300
  260 CONTINUE
C     .......... FORM SHIFT ..........
  300 IF (L .EQ. EN) GO TO 660
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 320
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1) - HI(ENM1,EN) * HI(EN,ENM1)
      XI = HR(ENM1,EN) * HI(EN,ENM1) + HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR .EQ. 0.0E0 .AND. XI .EQ. 0.0E0) GO TO 340
      YR = (HR(ENM1,ENM1) - SR) / 2.0E0
      YI = (HI(ENM1,ENM1) - SI) / 2.0E0
      CALL CSROOT(YR**2-YI**2+XR,2.0E0*YR*YI+XI,ZZR,ZZI)
      IF (YR * ZZR + YI * ZZI .GE. 0.0E0) GO TO 310
      ZZR = -ZZR
      ZZI = -ZZI
  310 CALL CDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI)
      SR = SR - XR
      SI = SI - XI
      GO TO 340
C     .......... FORM EXCEPTIONAL SHIFT ..........
  320 SR = ABS(HR(EN,ENM1)) + ABS(HR(ENM1,EN-2))
      SI = ABS(HI(EN,ENM1)) + ABS(HI(ENM1,EN-2))
C
  340 DO 360 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  360 CONTINUE
C
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      ITN = ITN - 1
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS ..........
      XR = ABS(HR(ENM1,ENM1)) + ABS(HI(ENM1,ENM1))
      YR = ABS(HR(EN,ENM1)) + ABS(HI(EN,ENM1))
      ZZR = ABS(HR(EN,EN)) + ABS(HI(EN,EN))
C     .......... FOR M=EN-1 STEP -1 UNTIL L DO -- ..........
      DO 380 MM = L, ENM1
         M = ENM1 + L - MM
         IF (M .EQ. L) GO TO 420
         YI = YR
         YR = ABS(HR(M,M-1)) + ABS(HI(M,M-1))
         XI = ZZR
         ZZR = XR
         XR = ABS(HR(M-1,M-1)) + ABS(HI(M-1,M-1))
         TST1 = ZZR / YI * (ZZR + XR + XI)
         TST2 = TST1 + YR
         IF (TST2 .EQ. TST1) GO TO 420
  380 CONTINUE
C     .......... TRIANGULAR DECOMPOSITION H=L*R ..........
  420 MP1 = M + 1
C
      DO 520 I = MP1, EN
         IM1 = I - 1
         XR = HR(IM1,IM1)
         XI = HI(IM1,IM1)
         YR = HR(I,IM1)
         YI = HI(I,IM1)
         IF (ABS(XR) + ABS(XI) .GE. ABS(YR) + ABS(YI)) GO TO 460
C     .......... INTERCHANGE ROWS OF HR AND HI ..........
         DO 440 J = IM1, EN
            ZZR = HR(IM1,J)
            HR(IM1,J) = HR(I,J)
            HR(I,J) = ZZR
            ZZI = HI(IM1,J)
            HI(IM1,J) = HI(I,J)
            HI(I,J) = ZZI
  440    CONTINUE
C
         CALL CDIV(XR,XI,YR,YI,ZZR,ZZI)
         WR(I) = 1.0E0
         GO TO 480
  460    CALL CDIV(YR,YI,XR,XI,ZZR,ZZI)
         WR(I) = -1.0E0
  480    HR(I,IM1) = ZZR
         HI(I,IM1) = ZZI
C
         DO 500 J = I, EN
            HR(I,J) = HR(I,J) - ZZR * HR(IM1,J) + ZZI * HI(IM1,J)
            HI(I,J) = HI(I,J) - ZZR * HI(IM1,J) - ZZI * HR(IM1,J)
  500    CONTINUE
C
  520 CONTINUE
C     .......... COMPOSITION R*L=H ..........
      DO 640 J = MP1, EN
         XR = HR(J,J-1)
         XI = HI(J,J-1)
         HR(J,J-1) = 0.0E0
         HI(J,J-1) = 0.0E0
C     .......... INTERCHANGE COLUMNS OF HR AND HI,
C                IF NECESSARY ..........
         IF (WR(J) .LE. 0.0E0) GO TO 580
C
         DO 540 I = L, J
            ZZR = HR(I,J-1)
            HR(I,J-1) = HR(I,J)
            HR(I,J) = ZZR
            ZZI = HI(I,J-1)
            HI(I,J-1) = HI(I,J)
            HI(I,J) = ZZI
  540    CONTINUE
C
  580    DO 600 I = L, J
            HR(I,J-1) = HR(I,J-1) + XR * HR(I,J) - XI * HI(I,J)
            HI(I,J-1) = HI(I,J-1) + XR * HI(I,J) + XI * HR(I,J)
  600    CONTINUE
C
  640 CONTINUE
C
      GO TO 240
C     .......... A ROOT FOUND ..........
  660 WR(EN) = HR(EN,EN) + TR
      WI(EN) = HI(EN,EN) + TI
      EN = ENM1
      GO TO 220
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END
      SUBROUTINE COMLR2(NM,N,LOW,IGH,INT,HR,HI,WR,WI,ZR,ZI,IERR)
C
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NM,NN,IGH,IM1,IP1,
     X        ITN,ITS,LOW,MP1,ENM1,IEND,IERR
      REAL HR(NM,N),HI(NM,N),WR(N),WI(N),ZR(NM,N),ZI(NM,N)
      REAL SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,TST1,TST2
      INTEGER INT(IGH)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE COMLR2,
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A COMPLEX UPPER HESSENBERG MATRIX BY THE MODIFIED LR
C     METHOD.  THE EIGENVECTORS OF A COMPLEX GENERAL MATRIX
C     CAN ALSO BE FOUND IF  COMHES  HAS BEEN USED TO REDUCE
C     THIS GENERAL MATRIX TO HESSENBERG FORM.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        INT CONTAINS INFORMATION ON THE ROWS AND COLUMNS INTERCHANGED
C          IN THE REDUCTION BY  COMHES, IF PERFORMED.  ONLY ELEMENTS
C          LOW THROUGH IGH ARE USED.  IF THE EIGENVECTORS OF THE HESSEN-
C          BERG MATRIX ARE DESIRED, SET INT(J)=J FOR THESE ELEMENTS.
C
C        HR AND HI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX.
C          THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAIN THE
C          MULTIPLIERS WHICH WERE USED IN THE REDUCTION BY  COMHES,
C          IF PERFORMED.  IF THE EIGENVECTORS OF THE HESSENBERG
C          MATRIX ARE DESIRED, THESE ELEMENTS MUST BE SET TO ZERO.
C
C     ON OUTPUT
C
C        THE UPPER HESSENBERG PORTIONS OF HR AND HI HAVE BEEN
C          DESTROYED, BUT THE LOCATION HR(1,1) CONTAINS THE NORM
C          OF THE TRIANGULARIZED MATRIX.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  IF AN ERROR
C          EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVECTORS.  THE EIGENVECTORS
C          ARE UNNORMALIZED.  IF AN ERROR EXIT IS MADE, NONE OF
C          THE EIGENVECTORS HAS BEEN FOUND.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C     CALLS CSROOT FOR COMPLEX SQUARE ROOT.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
C     .......... INITIALIZE EIGENVECTOR MATRIX ..........
      DO 100 I = 1, N
C
         DO 100 J = 1, N
            ZR(I,J) = 0.0E0
            ZI(I,J) = 0.0E0
            IF (I .EQ. J) ZR(I,J) = 1.0E0
  100 CONTINUE
C     .......... FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS
C                FROM THE INFORMATION LEFT BY COMHES ..........
      IEND = IGH - LOW - 1
      IF (IEND .LE. 0) GO TO 180
C     .......... FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 160 II = 1, IEND
         I = IGH - II
         IP1 = I + 1
C
         DO 120 K = IP1, IGH
            ZR(K,I) = HR(K,I-1)
            ZI(K,I) = HI(K,I-1)
  120    CONTINUE
C
         J = INT(I)
         IF (I .EQ. J) GO TO 160
C
         DO 140 K = I, IGH
            ZR(I,K) = ZR(J,K)
            ZI(I,K) = ZI(J,K)
            ZR(J,K) = 0.0E0
            ZI(J,K) = 0.0E0
  140    CONTINUE
C
         ZR(J,I) = 1.0E0
  160 CONTINUE
C     .......... STORE ROOTS ISOLATED BY CBAL ..........
  180 DO 200 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 200
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  200 CONTINUE
C
      EN = IGH
      TR = 0.0E0
      TI = 0.0E0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUE ..........
  220 IF (EN .LT. LOW) GO TO 680
      ITS = 0
      ENM1 = EN - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
  240 DO 260 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 300
         TST1 = ABS(HR(L-1,L-1)) + ABS(HI(L-1,L-1))
     X            + ABS(HR(L,L)) + ABS(HI(L,L))
         TST2 = TST1 + ABS(HR(L,L-1)) + ABS(HI(L,L-1))
         IF (TST2 .EQ. TST1) GO TO 300
  260 CONTINUE
C     .......... FORM SHIFT ..........
  300 IF (L .EQ. EN) GO TO 660
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 320
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1) - HI(ENM1,EN) * HI(EN,ENM1)
      XI = HR(ENM1,EN) * HI(EN,ENM1) + HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR .EQ. 0.0E0 .AND. XI .EQ. 0.0E0) GO TO 340
      YR = (HR(ENM1,ENM1) - SR) / 2.0E0
      YI = (HI(ENM1,ENM1) - SI) / 2.0E0
      CALL CSROOT(YR**2-YI**2+XR,2.0E0*YR*YI+XI,ZZR,ZZI)
      IF (YR * ZZR + YI * ZZI .GE. 0.0E0) GO TO 310
      ZZR = -ZZR
      ZZI = -ZZI
  310 CALL CDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI)
      SR = SR - XR
      SI = SI - XI
      GO TO 340
C     .......... FORM EXCEPTIONAL SHIFT ..........
  320 SR = ABS(HR(EN,ENM1)) + ABS(HR(ENM1,EN-2))
      SI = ABS(HI(EN,ENM1)) + ABS(HI(ENM1,EN-2))
C
  340 DO 360 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  360 CONTINUE
C
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      ITN = ITN - 1
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS ..........
      XR = ABS(HR(ENM1,ENM1)) + ABS(HI(ENM1,ENM1))
      YR = ABS(HR(EN,ENM1)) + ABS(HI(EN,ENM1))
      ZZR = ABS(HR(EN,EN)) + ABS(HI(EN,EN))
C     .......... FOR M=EN-1 STEP -1 UNTIL L DO -- ..........
      DO 380 MM = L, ENM1
         M = ENM1 + L - MM
         IF (M .EQ. L) GO TO 420
         YI = YR
         YR = ABS(HR(M,M-1)) + ABS(HI(M,M-1))
         XI = ZZR
         ZZR = XR
         XR = ABS(HR(M-1,M-1)) + ABS(HI(M-1,M-1))
         TST1 = ZZR / YI * (ZZR + XR + XI)
         TST2 = TST1 + YR
         IF (TST2 .EQ. TST1) GO TO 420
  380 CONTINUE
C     .......... TRIANGULAR DECOMPOSITION H=L*R ..........
  420 MP1 = M + 1
C
      DO 520 I = MP1, EN
         IM1 = I - 1
         XR = HR(IM1,IM1)
         XI = HI(IM1,IM1)
         YR = HR(I,IM1)
         YI = HI(I,IM1)
         IF (ABS(XR) + ABS(XI) .GE. ABS(YR) + ABS(YI)) GO TO 460
C     .......... INTERCHANGE ROWS OF HR AND HI ..........
         DO 440 J = IM1, N
            ZZR = HR(IM1,J)
            HR(IM1,J) = HR(I,J)
            HR(I,J) = ZZR
            ZZI = HI(IM1,J)
            HI(IM1,J) = HI(I,J)
            HI(I,J) = ZZI
  440    CONTINUE
C
         CALL CDIV(XR,XI,YR,YI,ZZR,ZZI)
         WR(I) = 1.0E0
         GO TO 480
  460    CALL CDIV(YR,YI,XR,XI,ZZR,ZZI)
         WR(I) = -1.0E0
  480    HR(I,IM1) = ZZR
         HI(I,IM1) = ZZI
C
         DO 500 J = I, N
            HR(I,J) = HR(I,J) - ZZR * HR(IM1,J) + ZZI * HI(IM1,J)
            HI(I,J) = HI(I,J) - ZZR * HI(IM1,J) - ZZI * HR(IM1,J)
  500    CONTINUE
C
  520 CONTINUE
C     .......... COMPOSITION R*L=H ..........
      DO 640 J = MP1, EN
         XR = HR(J,J-1)
         XI = HI(J,J-1)
         HR(J,J-1) = 0.0E0
         HI(J,J-1) = 0.0E0
C     .......... INTERCHANGE COLUMNS OF HR, HI, ZR, AND ZI,
C                IF NECESSARY ..........
         IF (WR(J) .LE. 0.0E0) GO TO 580
C
         DO 540 I = 1, J
            ZZR = HR(I,J-1)
            HR(I,J-1) = HR(I,J)
            HR(I,J) = ZZR
            ZZI = HI(I,J-1)
            HI(I,J-1) = HI(I,J)
            HI(I,J) = ZZI
  540    CONTINUE
C
         DO 560 I = LOW, IGH
            ZZR = ZR(I,J-1)
            ZR(I,J-1) = ZR(I,J)
            ZR(I,J) = ZZR
            ZZI = ZI(I,J-1)
            ZI(I,J-1) = ZI(I,J)
            ZI(I,J) = ZZI
  560    CONTINUE
C
  580    DO 600 I = 1, J
            HR(I,J-1) = HR(I,J-1) + XR * HR(I,J) - XI * HI(I,J)
            HI(I,J-1) = HI(I,J-1) + XR * HI(I,J) + XI * HR(I,J)
  600    CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
         DO 620 I = LOW, IGH
            ZR(I,J-1) = ZR(I,J-1) + XR * ZR(I,J) - XI * ZI(I,J)
            ZI(I,J-1) = ZI(I,J-1) + XR * ZI(I,J) + XI * ZR(I,J)
  620    CONTINUE
C
  640 CONTINUE
C
      GO TO 240
C     .......... A ROOT FOUND ..........
  660 HR(EN,EN) = HR(EN,EN) + TR
      WR(EN) = HR(EN,EN)
      HI(EN,EN) = HI(EN,EN) + TI
      WI(EN) = HI(EN,EN)
      EN = ENM1
      GO TO 220
C     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM ..........
  680 NORM = 0.0E0
C
      DO 720 I = 1, N
C
         DO 720 J = I, N
            TR = ABS(HR(I,J)) + ABS(HI(I,J))
            IF (TR .GT. NORM) NORM = TR
  720 CONTINUE
C
      HR(1,1) = NORM
      IF (N .EQ. 1 .OR. NORM .EQ. 0.0E0) GO TO 1001
C     .......... FOR EN=N STEP -1 UNTIL 2 DO -- ..........
      DO 800 NN = 2, N
         EN = N + 2 - NN
         XR = WR(EN)
         XI = WI(EN)
         HR(EN,EN) = 1.0E0
         HI(EN,EN) = 0.0E0
         ENM1 = EN - 1
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
         DO 780 II = 1, ENM1
            I = EN - II
            ZZR = 0.0E0
            ZZI = 0.0E0
            IP1 = I + 1
C
            DO 740 J = IP1, EN
               ZZR = ZZR + HR(I,J) * HR(J,EN) - HI(I,J) * HI(J,EN)
               ZZI = ZZI + HR(I,J) * HI(J,EN) + HI(I,J) * HR(J,EN)
  740       CONTINUE
C
            YR = XR - WR(I)
            YI = XI - WI(I)
            IF (YR .NE. 0.0E0 .OR. YI .NE. 0.0E0) GO TO 765
               TST1 = NORM
               YR = TST1
  760          YR = 0.01E0 * YR
               TST2 = NORM + YR
               IF (TST2 .GT. TST1) GO TO 760
  765       CONTINUE
            CALL CDIV(ZZR,ZZI,YR,YI,HR(I,EN),HI(I,EN))
C     .......... OVERFLOW CONTROL ..........
            TR = ABS(HR(I,EN)) + ABS(HI(I,EN))
            IF (TR .EQ. 0.0E0) GO TO 780
            TST1 = TR
            TST2 = TST1 + 1.0E0/TST1
            IF (TST2 .GT. TST1) GO TO 780
            DO 770 J = I, EN
               HR(J,EN) = HR(J,EN)/TR
               HI(J,EN) = HI(J,EN)/TR
  770       CONTINUE
C
  780    CONTINUE
C
  800 CONTINUE
C     .......... END BACKSUBSTITUTION ..........
      ENM1 = N - 1
C     .......... VECTORS OF ISOLATED ROOTS ..........
      DO  840 I = 1, ENM1
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840
         IP1 = I + 1
C
         DO 820 J = IP1, N
            ZR(I,J) = HR(I,J)
            ZI(I,J) = HI(I,J)
  820    CONTINUE
C
  840 CONTINUE
C     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW+1 DO -- ..........
      DO 880 JJ = LOW, ENM1
         J = N + LOW - JJ
         M = MIN0(J,IGH)
C
         DO 880 I = LOW, IGH
            ZZR = 0.0E0
            ZZI = 0.0E0
C
            DO 860 K = LOW, M
               ZZR = ZZR + ZR(I,K) * HR(K,J) - ZI(I,K) * HI(K,J)
               ZZI = ZZI + ZR(I,K) * HI(K,J) + ZI(I,K) * HR(K,J)
  860       CONTINUE
C
            ZR(I,J) = ZZR
            ZI(I,J) = ZZI
  880 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END
      SUBROUTINE COMQR(NM,N,LOW,IGH,HR,HI,WR,WI,IERR)
C
      INTEGER I,J,L,N,EN,LL,NM,IGH,ITN,ITS,LOW,LP1,ENM1,IERR
      REAL HR(NM,N),HI(NM,N),WR(N),WI(N)
      REAL SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,TST1,TST2,
     X       PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF A UNITARY ANALOGUE OF THE
C     ALGOL PROCEDURE  COMLR, NUM. MATH. 12, 369-376(1968) BY MARTIN
C     AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 396-403(1971).
C     THE UNITARY ANALOGUE SUBSTITUTES THE QR ALGORITHM OF FRANCIS
C     (COMP. JOUR. 4, 332-345(1962)) FOR THE LR ALGORITHM.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A COMPLEX
C     UPPER HESSENBERG MATRIX BY THE QR METHOD.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        HR AND HI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX.
C          THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAIN
C          INFORMATION ABOUT THE UNITARY TRANSFORMATIONS USED IN
C          THE REDUCTION BY  CORTH, IF PERFORMED.
C
C     ON OUTPUT
C
C        THE UPPER HESSENBERG PORTIONS OF HR AND HI HAVE BEEN
C          DESTROYED.  THEREFORE, THEY MUST BE SAVED BEFORE
C          CALLING  COMQR  IF SUBSEQUENT CALCULATION OF
C          EIGENVECTORS IS TO BE PERFORMED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  IF AN ERROR
C          EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C     CALLS CSROOT FOR COMPLEX SQUARE ROOT.
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (LOW .EQ. IGH) GO TO 180
C     .......... CREATE REAL SUBDIAGONAL ELEMENTS ..........
      L = LOW + 1
C
      DO 170 I = L, IGH
         LL = MIN0(I+1,IGH)
         IF (HI(I,I-1) .EQ. 0.0E0) GO TO 170
         NORM = PYTHAG(HR(I,I-1),HI(I,I-1))
         YR = HR(I,I-1) / NORM
         YI = HI(I,I-1) / NORM
         HR(I,I-1) = NORM
         HI(I,I-1) = 0.0E0
C
         DO 155 J = I, IGH
            SI = YR * HI(I,J) - YI * HR(I,J)
            HR(I,J) = YR * HR(I,J) + YI * HI(I,J)
            HI(I,J) = SI
  155    CONTINUE
C
         DO 160 J = LOW, LL
            SI = YR * HI(J,I) + YI * HR(J,I)
            HR(J,I) = YR * HR(J,I) - YI * HI(J,I)
            HI(J,I) = SI
  160    CONTINUE
C
  170 CONTINUE
C     .......... STORE ROOTS ISOLATED BY CBAL ..........
  180 DO 200 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 200
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  200 CONTINUE
C
      EN = IGH
      TR = 0.0E0
      TI = 0.0E0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUE ..........
  220 IF (EN .LT. LOW) GO TO 1001
      ITS = 0
      ENM1 = EN - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW E0 -- ..........
  240 DO 260 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 300
         TST1 = ABS(HR(L-1,L-1)) + ABS(HI(L-1,L-1))
     X            + ABS(HR(L,L)) + ABS(HI(L,L))
         TST2 = TST1 + ABS(HR(L,L-1))
         IF (TST2 .EQ. TST1) GO TO 300
  260 CONTINUE
C     .......... FORM SHIFT ..........
  300 IF (L .EQ. EN) GO TO 660
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 320
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1)
      XI = HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR .EQ. 0.0E0 .AND. XI .EQ. 0.0E0) GO TO 340
      YR = (HR(ENM1,ENM1) - SR) / 2.0E0
      YI = (HI(ENM1,ENM1) - SI) / 2.0E0
      CALL CSROOT(YR**2-YI**2+XR,2.0E0*YR*YI+XI,ZZR,ZZI)
      IF (YR * ZZR + YI * ZZI .GE. 0.0E0) GO TO 310
      ZZR = -ZZR
      ZZI = -ZZI
  310 CALL CDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI)
      SR = SR - XR
      SI = SI - XI
      GO TO 340
C     .......... FORM EXCEPTIONAL SHIFT ..........
  320 SR = ABS(HR(EN,ENM1)) + ABS(HR(ENM1,EN-2))
      SI = 0.0E0
C
  340 DO 360 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  360 CONTINUE
C
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      ITN = ITN - 1
C     .......... REDUCE TO TRIANGLE (ROWS) ..........
      LP1 = L + 1
C
      DO 500 I = LP1, EN
         SR = HR(I,I-1)
         HR(I,I-1) = 0.0E0
         NORM = PYTHAG(PYTHAG(HR(I-1,I-1),HI(I-1,I-1)),SR)
         XR = HR(I-1,I-1) / NORM
         WR(I-1) = XR
         XI = HI(I-1,I-1) / NORM
         WI(I-1) = XI
         HR(I-1,I-1) = NORM
         HI(I-1,I-1) = 0.0E0
         HI(I,I-1) = SR / NORM
C
         DO 490 J = I, EN
            YR = HR(I-1,J)
            YI = HI(I-1,J)
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            HR(I-1,J) = XR * YR + XI * YI + HI(I,I-1) * ZZR
            HI(I-1,J) = XR * YI - XI * YR + HI(I,I-1) * ZZI
            HR(I,J) = XR * ZZR - XI * ZZI - HI(I,I-1) * YR
            HI(I,J) = XR * ZZI + XI * ZZR - HI(I,I-1) * YI
  490    CONTINUE
C
  500 CONTINUE
C
      SI = HI(EN,EN)
      IF (SI .EQ. 0.0E0) GO TO 540
      NORM = PYTHAG(HR(EN,EN),SI)
      SR = HR(EN,EN) / NORM
      SI = SI / NORM
      HR(EN,EN) = NORM
      HI(EN,EN) = 0.0E0
C     .......... INVERSE OPERATION (COLUMNS) ..........
  540 DO 600 J = LP1, EN
         XR = WR(J-1)
         XI = WI(J-1)
C
         DO 580 I = L, J
            YR = HR(I,J-1)
            YI = 0.0E0
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            IF (I .EQ. J) GO TO 560
            YI = HI(I,J-1)
            HI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
  560       HR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            HR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            HI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  580    CONTINUE
C
  600 CONTINUE
C
      IF (SI .EQ. 0.0E0) GO TO 240
C
      DO 630 I = L, EN
         YR = HR(I,EN)
         YI = HI(I,EN)
         HR(I,EN) = SR * YR - SI * YI
         HI(I,EN) = SR * YI + SI * YR
  630 CONTINUE
C
      GO TO 240
C     .......... A ROOT FOUND ..........
  660 WR(EN) = HR(EN,EN) + TR
      WI(EN) = HI(EN,EN) + TI
      EN = ENM1
      GO TO 220
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END
      SUBROUTINE COMQR2(NM,N,LOW,IGH,ORTR,ORTI,HR,HI,WR,WI,ZR,ZI,IERR)
C
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,NM,NN,IGH,IP1,
     X        ITN,ITS,LOW,LP1,ENM1,IEND,IERR
      REAL HR(NM,N),HI(NM,N),WR(N),WI(N),ZR(NM,N),ZI(NM,N),
     X       ORTR(IGH),ORTI(IGH)
      REAL SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,TST1,TST2,
     X       PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF A UNITARY ANALOGUE OF THE
C     ALGOL PROCEDURE  COMLR2, NUM. MATH. 16, 181-204(1970) BY PETERS
C     AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C     THE UNITARY ANALOGUE SUBSTITUTES THE QR ALGORITHM OF FRANCIS
C     (COMP. JOUR. 4, 332-345(1962)) FOR THE LR ALGORITHM.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A COMPLEX UPPER HESSENBERG MATRIX BY THE QR
C     METHOD.  THE EIGENVECTORS OF A COMPLEX GENERAL MATRIX
C     CAN ALSO BE FOUND IF  CORTH  HAS BEEN USED TO REDUCE
C     THIS GENERAL MATRIX TO HESSENBERG FORM.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        ORTR AND ORTI CONTAIN INFORMATION ABOUT THE UNITARY TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  CORTH, IF PERFORMED.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.  IF THE EIGENVECTORS
C          OF THE HESSENBERG MATRIX ARE DESIRED, SET ORTR(J) AND
C          ORTI(J) TO 0.0E0 FOR THESE ELEMENTS.
C
C        HR AND HI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX.
C          THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAIN FURTHER
C          INFORMATION ABOUT THE TRANSFORMATIONS WHICH WERE USED IN THE
C          REDUCTION BY  CORTH, IF PERFORMED.  IF THE EIGENVECTORS OF
C          THE HESSENBERG MATRIX ARE DESIRED, THESE ELEMENTS MAY BE
C          ARBITRARY.
C
C     ON OUTPUT
C
C        ORTR, ORTI, AND THE UPPER HESSENBERG PORTIONS OF HR AND HI
C          HAVE BEEN DESTROYED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  IF AN ERROR
C          EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVECTORS.  THE EIGENVECTORS
C          ARE UNNORMALIZED.  IF AN ERROR EXIT IS MADE, NONE OF
C          THE EIGENVECTORS HAS BEEN FOUND.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C     CALLS CSROOT FOR COMPLEX SQUARE ROOT.
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
C     .......... INITIALIZE EIGENVECTOR MATRIX ..........
      DO 101 J = 1, N
C
         DO 100 I = 1, N
            ZR(I,J) = 0.0E0
            ZI(I,J) = 0.0E0
  100    CONTINUE
         ZR(J,J) = 1.0E0
  101 CONTINUE
C     .......... FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS
C                FROM THE INFORMATION LEFT BY CORTH ..........
      IEND = IGH - LOW - 1
      IF (IEND) 180, 150, 105
C     .......... FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
  105 DO 140 II = 1, IEND
         I = IGH - II
         IF (ORTR(I) .EQ. 0.0E0 .AND. ORTI(I) .EQ. 0.0E0) GO TO 140
         IF (HR(I,I-1) .EQ. 0.0E0 .AND. HI(I,I-1) .EQ. 0.0E0) GO TO 140
C     .......... NORM BELOW IS NEGATIVE OF H FORMED IN CORTH ..........
         NORM = HR(I,I-1) * ORTR(I) + HI(I,I-1) * ORTI(I)
         IP1 = I + 1
C
         DO 110 K = IP1, IGH
            ORTR(K) = HR(K,I-1)
            ORTI(K) = HI(K,I-1)
  110    CONTINUE
C
         DO 130 J = I, IGH
            SR = 0.0E0
            SI = 0.0E0
C
            DO 115 K = I, IGH
               SR = SR + ORTR(K) * ZR(K,J) + ORTI(K) * ZI(K,J)
               SI = SI + ORTR(K) * ZI(K,J) - ORTI(K) * ZR(K,J)
  115       CONTINUE
C
            SR = SR / NORM
            SI = SI / NORM
C
            DO 120 K = I, IGH
               ZR(K,J) = ZR(K,J) + SR * ORTR(K) - SI * ORTI(K)
               ZI(K,J) = ZI(K,J) + SR * ORTI(K) + SI * ORTR(K)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C     .......... CREATE REAL SUBDIAGONAL ELEMENTS ..........
  150 L = LOW + 1
C
      DO 170 I = L, IGH
         LL = MIN0(I+1,IGH)
         IF (HI(I,I-1) .EQ. 0.0E0) GO TO 170
         NORM = PYTHAG(HR(I,I-1),HI(I,I-1))
         YR = HR(I,I-1) / NORM
         YI = HI(I,I-1) / NORM
         HR(I,I-1) = NORM
         HI(I,I-1) = 0.0E0
C
         DO 155 J = I, N
            SI = YR * HI(I,J) - YI * HR(I,J)
            HR(I,J) = YR * HR(I,J) + YI * HI(I,J)
            HI(I,J) = SI
  155    CONTINUE
C
         DO 160 J = 1, LL
            SI = YR * HI(J,I) + YI * HR(J,I)
            HR(J,I) = YR * HR(J,I) - YI * HI(J,I)
            HI(J,I) = SI
  160    CONTINUE
C
         DO 165 J = LOW, IGH
            SI = YR * ZI(J,I) + YI * ZR(J,I)
            ZR(J,I) = YR * ZR(J,I) - YI * ZI(J,I)
            ZI(J,I) = SI
  165    CONTINUE
C
  170 CONTINUE
C     .......... STORE ROOTS ISOLATED BY CBAL ..........
  180 DO 200 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 200
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  200 CONTINUE
C
      EN = IGH
      TR = 0.0E0
      TI = 0.0E0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUE ..........
  220 IF (EN .LT. LOW) GO TO 680
      ITS = 0
      ENM1 = EN - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
  240 DO 260 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 300
         TST1 = ABS(HR(L-1,L-1)) + ABS(HI(L-1,L-1))
     X            + ABS(HR(L,L)) + ABS(HI(L,L))
         TST2 = TST1 + ABS(HR(L,L-1))
         IF (TST2 .EQ. TST1) GO TO 300
  260 CONTINUE
C     .......... FORM SHIFT ..........
  300 IF (L .EQ. EN) GO TO 660
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 320
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1)
      XI = HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR .EQ. 0.0E0 .AND. XI .EQ. 0.0E0) GO TO 340
      YR = (HR(ENM1,ENM1) - SR) / 2.0E0
      YI = (HI(ENM1,ENM1) - SI) / 2.0E0
      CALL CSROOT(YR**2-YI**2+XR,2.0E0*YR*YI+XI,ZZR,ZZI)
      IF (YR * ZZR + YI * ZZI .GE. 0.0E0) GO TO 310
      ZZR = -ZZR
      ZZI = -ZZI
  310 CALL CDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI)
      SR = SR - XR
      SI = SI - XI
      GO TO 340
C     .......... FORM EXCEPTIONAL SHIFT ..........
  320 SR = ABS(HR(EN,ENM1)) + ABS(HR(ENM1,EN-2))
      SI = 0.0E0
C
  340 DO 360 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  360 CONTINUE
C
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      ITN = ITN - 1
C     .......... REDUCE TO TRIANGLE (ROWS) ..........
      LP1 = L + 1
C
      DO 500 I = LP1, EN
         SR = HR(I,I-1)
         HR(I,I-1) = 0.0E0
         NORM = PYTHAG(PYTHAG(HR(I-1,I-1),HI(I-1,I-1)),SR)
         XR = HR(I-1,I-1) / NORM
         WR(I-1) = XR
         XI = HI(I-1,I-1) / NORM
         WI(I-1) = XI
         HR(I-1,I-1) = NORM
         HI(I-1,I-1) = 0.0E0
         HI(I,I-1) = SR / NORM
C
         DO 490 J = I, N
            YR = HR(I-1,J)
            YI = HI(I-1,J)
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            HR(I-1,J) = XR * YR + XI * YI + HI(I,I-1) * ZZR
            HI(I-1,J) = XR * YI - XI * YR + HI(I,I-1) * ZZI
            HR(I,J) = XR * ZZR - XI * ZZI - HI(I,I-1) * YR
            HI(I,J) = XR * ZZI + XI * ZZR - HI(I,I-1) * YI
  490    CONTINUE
C
  500 CONTINUE
C
      SI = HI(EN,EN)
      IF (SI .EQ. 0.0E0) GO TO 540
      NORM = PYTHAG(HR(EN,EN),SI)
      SR = HR(EN,EN) / NORM
      SI = SI / NORM
      HR(EN,EN) = NORM
      HI(EN,EN) = 0.0E0
      IF (EN .EQ. N) GO TO 540
      IP1 = EN + 1
C
      DO 520 J = IP1, N
         YR = HR(EN,J)
         YI = HI(EN,J)
         HR(EN,J) = SR * YR + SI * YI
         HI(EN,J) = SR * YI - SI * YR
  520 CONTINUE
C     .......... INVERSE OPERATION (COLUMNS) ..........
  540 DO 600 J = LP1, EN
         XR = WR(J-1)
         XI = WI(J-1)
C
         DO 580 I = 1, J
            YR = HR(I,J-1)
            YI = 0.0E0
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            IF (I .EQ. J) GO TO 560
            YI = HI(I,J-1)
            HI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
  560       HR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            HR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            HI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  580    CONTINUE
C
         DO 590 I = LOW, IGH
            YR = ZR(I,J-1)
            YI = ZI(I,J-1)
            ZZR = ZR(I,J)
            ZZI = ZI(I,J)
            ZR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            ZI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
            ZR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            ZI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  590    CONTINUE
C
  600 CONTINUE
C
      IF (SI .EQ. 0.0E0) GO TO 240
C
      DO 630 I = 1, EN
         YR = HR(I,EN)
         YI = HI(I,EN)
         HR(I,EN) = SR * YR - SI * YI
         HI(I,EN) = SR * YI + SI * YR
  630 CONTINUE
C
      DO 640 I = LOW, IGH
         YR = ZR(I,EN)
         YI = ZI(I,EN)
         ZR(I,EN) = SR * YR - SI * YI
         ZI(I,EN) = SR * YI + SI * YR
  640 CONTINUE
C
      GO TO 240
C     .......... A ROOT FOUND ..........
  660 HR(EN,EN) = HR(EN,EN) + TR
      WR(EN) = HR(EN,EN)
      HI(EN,EN) = HI(EN,EN) + TI
      WI(EN) = HI(EN,EN)
      EN = ENM1
      GO TO 220
C     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM ..........
  680 NORM = 0.0E0
C
      DO 720 I = 1, N
C
         DO 720 J = I, N
            TR = ABS(HR(I,J)) + ABS(HI(I,J))
            IF (TR .GT. NORM) NORM = TR
  720 CONTINUE
C
      IF (N .EQ. 1 .OR. NORM .EQ. 0.0E0) GO TO 1001
C     .......... FOR EN=N STEP -1 UNTIL 2 DO -- ..........
      DO 800 NN = 2, N
         EN = N + 2 - NN
         XR = WR(EN)
         XI = WI(EN)
         HR(EN,EN) = 1.0E0
         HI(EN,EN) = 0.0E0
         ENM1 = EN - 1
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
         DO 780 II = 1, ENM1
            I = EN - II
            ZZR = 0.0E0
            ZZI = 0.0E0
            IP1 = I + 1
C
            DO 740 J = IP1, EN
               ZZR = ZZR + HR(I,J) * HR(J,EN) - HI(I,J) * HI(J,EN)
               ZZI = ZZI + HR(I,J) * HI(J,EN) + HI(I,J) * HR(J,EN)
  740       CONTINUE
C
            YR = XR - WR(I)
            YI = XI - WI(I)
            IF (YR .NE. 0.0E0 .OR. YI .NE. 0.0E0) GO TO 765
               TST1 = NORM
               YR = TST1
  760          YR = 0.01E0 * YR
               TST2 = NORM + YR
               IF (TST2 .GT. TST1) GO TO 760
  765       CONTINUE
            CALL CDIV(ZZR,ZZI,YR,YI,HR(I,EN),HI(I,EN))
C     .......... OVERFLOW CONTROL ..........
            TR = ABS(HR(I,EN)) + ABS(HI(I,EN))
            IF (TR .EQ. 0.0E0) GO TO 780
            TST1 = TR
            TST2 = TST1 + 1.0E0/TST1
            IF (TST2 .GT. TST1) GO TO 780
            DO 770 J = I, EN
               HR(J,EN) = HR(J,EN)/TR
               HI(J,EN) = HI(J,EN)/TR
  770       CONTINUE
C
  780    CONTINUE
C
  800 CONTINUE
C     .......... END BACKSUBSTITUTION ..........
      ENM1 = N - 1
C     .......... VECTORS OF ISOLATED ROOTS ..........
      DO  840 I = 1, ENM1
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840
         IP1 = I + 1
C
         DO 820 J = IP1, N
            ZR(I,J) = HR(I,J)
            ZI(I,J) = HI(I,J)
  820    CONTINUE
C
  840 CONTINUE
C     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW+1 DO -- ..........
      DO 880 JJ = LOW, ENM1
         J = N + LOW - JJ
         M = MIN0(J,IGH)
C
         DO 880 I = LOW, IGH
            ZZR = 0.0E0
            ZZI = 0.0E0
C
            DO 860 K = LOW, M
               ZZR = ZZR + ZR(I,K) * HR(K,J) - ZI(I,K) * HI(K,J)
               ZZI = ZZI + ZR(I,K) * HI(K,J) + ZI(I,K) * HR(K,J)
  860       CONTINUE
C
            ZR(I,J) = ZZR
            ZI(I,J) = ZZI
  880 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END
      SUBROUTINE CORTB(NM,LOW,IGH,AR,AI,ORTR,ORTI,M,ZR,ZI)
C
      INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
      REAL AR(NM,IGH),AI(NM,IGH),ORTR(IGH),ORTI(IGH),
     X       ZR(NM,M),ZI(NM,M)
      REAL H,GI,GR
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE ORTBAK, NUM. MATH. 12, 349-368(1968)
C     BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX GENERAL
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     UPPER HESSENBERG MATRIX DETERMINED BY  CORTH.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1 AND IGH EQUAL TO THE ORDER OF THE MATRIX.
C
C        AR AND AI CONTAIN INFORMATION ABOUT THE UNITARY
C          TRANSFORMATIONS USED IN THE REDUCTION BY  CORTH
C          IN THEIR STRICT LOWER TRIANGLES.
C
C        ORTR AND ORTI CONTAIN FURTHER INFORMATION ABOUT THE
C          TRANSFORMATIONS USED IN THE REDUCTION BY  CORTH.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C        M IS THE NUMBER OF COLUMNS OF ZR AND ZI TO BE BACK TRANSFORMED.
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVECTORS TO BE
C          BACK TRANSFORMED IN THEIR FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE TRANSFORMED EIGENVECTORS
C          IN THEIR FIRST M COLUMNS.
C
C        ORTR AND ORTI HAVE BEEN ALTERED.
C
C     NOTE THAT CORTB PRESERVES VECTOR EUCLIDEAN NORMS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 140 MM = KP1, LA
         MP = LOW + IGH - MM
         IF (AR(MP,MP-1) .EQ. 0.0E0 .AND. AI(MP,MP-1) .EQ. 0.0E0)
     X      GO TO 140
C     .......... H BELOW IS NEGATIVE OF H FORMED IN CORTH ..........
         H = AR(MP,MP-1) * ORTR(MP) + AI(MP,MP-1) * ORTI(MP)
         MP1 = MP + 1
C
         DO 100 I = MP1, IGH
            ORTR(I) = AR(I,MP-1)
            ORTI(I) = AI(I,MP-1)
  100    CONTINUE
C
         DO 130 J = 1, M
            GR = 0.0E0
            GI = 0.0E0
C
            DO 110 I = MP, IGH
               GR = GR + ORTR(I) * ZR(I,J) + ORTI(I) * ZI(I,J)
               GI = GI + ORTR(I) * ZI(I,J) - ORTI(I) * ZR(I,J)
  110       CONTINUE
C
            GR = GR / H
            GI = GI / H
C
            DO 120 I = MP, IGH
               ZR(I,J) = ZR(I,J) + GR * ORTR(I) - GI * ORTI(I)
               ZI(I,J) = ZI(I,J) + GR * ORTI(I) + GI * ORTR(I)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE CORTH(NM,N,LOW,IGH,AR,AI,ORTR,ORTI)
C
      INTEGER I,J,M,N,II,JJ,LA,MP,NM,IGH,KP1,LOW
      REAL AR(NM,N),AI(NM,N),ORTR(IGH),ORTI(IGH)
      REAL F,G,H,FI,FR,SCALE,PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE ORTHES, NUM. MATH. 12, 349-368(1968)
C     BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     GIVEN A COMPLEX GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     UNITARY SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX INPUT MATRIX.
C
C     ON OUTPUT
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE HESSENBERG MATRIX.  INFORMATION
C          ABOUT THE UNITARY TRANSFORMATIONS USED IN THE REDUCTION
C          IS STORED IN THE REMAINING TRIANGLES UNDER THE
C          HESSENBERG MATRIX.
C
C        ORTR AND ORTI CONTAIN FURTHER INFORMATION ABOUT THE
C          TRANSFORMATIONS.  ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C
      DO 180 M = KP1, LA
         H = 0.0E0
         ORTR(M) = 0.0E0
         ORTI(M) = 0.0E0
         SCALE = 0.0E0
C     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ..........
         DO 90 I = M, IGH
   90    SCALE = SCALE + ABS(AR(I,M-1)) + ABS(AI(I,M-1))
C
         IF (SCALE .EQ. 0.0E0) GO TO 180
         MP = M + IGH
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
         DO 100 II = M, IGH
            I = MP - II
            ORTR(I) = AR(I,M-1) / SCALE
            ORTI(I) = AI(I,M-1) / SCALE
            H = H + ORTR(I) * ORTR(I) + ORTI(I) * ORTI(I)
  100    CONTINUE
C
         G = SQRT(H)
         F = PYTHAG(ORTR(M),ORTI(M))
         IF (F .EQ. 0.0E0) GO TO 103
         H = H + F * G
         G = G / F
         ORTR(M) = (1.0E0 + G) * ORTR(M)
         ORTI(M) = (1.0E0 + G) * ORTI(M)
         GO TO 105
C
  103    ORTR(M) = G
         AR(M,M-1) = SCALE
C     .......... FORM (I-(U*UT)/H) * A ..........
  105    DO 130 J = M, N
            FR = 0.0E0
            FI = 0.0E0
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
            DO 110 II = M, IGH
               I = MP - II
               FR = FR + ORTR(I) * AR(I,J) + ORTI(I) * AI(I,J)
               FI = FI + ORTR(I) * AI(I,J) - ORTI(I) * AR(I,J)
  110       CONTINUE
C
            FR = FR / H
            FI = FI / H
C
            DO 120 I = M, IGH
               AR(I,J) = AR(I,J) - FR * ORTR(I) + FI * ORTI(I)
               AI(I,J) = AI(I,J) - FR * ORTI(I) - FI * ORTR(I)
  120       CONTINUE
C
  130    CONTINUE
C     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ..........
         DO 160 I = 1, IGH
            FR = 0.0E0
            FI = 0.0E0
C     .......... FOR J=IGH STEP -1 UNTIL M DO -- ..........
            DO 140 JJ = M, IGH
               J = MP - JJ
               FR = FR + ORTR(J) * AR(I,J) - ORTI(J) * AI(I,J)
               FI = FI + ORTR(J) * AI(I,J) + ORTI(J) * AR(I,J)
  140       CONTINUE
C
            FR = FR / H
            FI = FI / H
C
            DO 150 J = M, IGH
               AR(I,J) = AR(I,J) - FR * ORTR(J) - FI * ORTI(J)
               AI(I,J) = AI(I,J) + FR * ORTI(J) - FI * ORTR(J)
  150       CONTINUE
C
  160    CONTINUE
C
         ORTR(M) = SCALE * ORTR(M)
         ORTI(M) = SCALE * ORTI(M)
         AR(M,M-1) = -G * AR(M,M-1)
         AI(M,M-1) = -G * AI(M,M-1)
  180 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE ELMBAK(NM,LOW,IGH,A,INT,M,Z)
C
      INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
      REAL A(NM,IGH),Z(NM,M)
      REAL X
      INTEGER INT(IGH)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ELMBAK,
C     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL GENERAL
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     UPPER HESSENBERG MATRIX DETERMINED BY  ELMHES.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1 AND IGH EQUAL TO THE ORDER OF THE MATRIX.
C
C        A CONTAINS THE MULTIPLIERS WHICH WERE USED IN THE
C          REDUCTION BY  ELMHES  IN ITS LOWER TRIANGLE
C          BELOW THE SUBDIAGONAL.
C
C        INT CONTAINS INFORMATION ON THE ROWS AND COLUMNS
C          INTERCHANGED IN THE REDUCTION BY  ELMHES.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C        M IS THE NUMBER OF COLUMNS OF Z TO BE BACK TRANSFORMED.
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGEN-
C          VECTORS TO BE BACK TRANSFORMED IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE
C          TRANSFORMED EIGENVECTORS IN ITS FIRST M COLUMNS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 140 MM = KP1, LA
         MP = LOW + IGH - MM
         MP1 = MP + 1
C
         DO 110 I = MP1, IGH
            X = A(I,MP-1)
            IF (X .EQ. 0.0E0) GO TO 110
C
            DO 100 J = 1, M
  100       Z(I,J) = Z(I,J) + X * Z(MP,J)
C
  110    CONTINUE
C
         I = INT(MP)
         IF (I .EQ. MP) GO TO 140
C
         DO 130 J = 1, M
            X = Z(I,J)
            Z(I,J) = Z(MP,J)
            Z(MP,J) = X
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE ELMHES(NM,N,LOW,IGH,A,INT)
C
      INTEGER I,J,M,N,LA,NM,IGH,KP1,LOW,MM1,MP1
      REAL A(NM,N)
      REAL X,Y
      INTEGER INT(IGH)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ELMHES,
C     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     GIVEN A REAL GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     STABILIZED ELEMENTARY SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        A CONTAINS THE INPUT MATRIX.
C
C     ON OUTPUT
C
C        A CONTAINS THE HESSENBERG MATRIX.  THE MULTIPLIERS
C          WHICH WERE USED IN THE REDUCTION ARE STORED IN THE
C          REMAINING TRIANGLE UNDER THE HESSENBERG MATRIX.
C
C        INT CONTAINS INFORMATION ON THE ROWS AND COLUMNS
C          INTERCHANGED IN THE REDUCTION.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C
      DO 180 M = KP1, LA
         MM1 = M - 1
         X = 0.0E0
         I = M
C
         DO 100 J = M, IGH
            IF (ABS(A(J,MM1)) .LE. ABS(X)) GO TO 100
            X = A(J,MM1)
            I = J
  100    CONTINUE
C
         INT(M) = I
         IF (I .EQ. M) GO TO 130
C     .......... INTERCHANGE ROWS AND COLUMNS OF A ..........
         DO 110 J = MM1, N
            Y = A(I,J)
            A(I,J) = A(M,J)
            A(M,J) = Y
  110    CONTINUE
C
         DO 120 J = 1, IGH
            Y = A(J,I)
            A(J,I) = A(J,M)
            A(J,M) = Y
  120    CONTINUE
C     .......... END INTERCHANGE ..........
  130    IF (X .EQ. 0.0E0) GO TO 180
         MP1 = M + 1
C
         DO 160 I = MP1, IGH
            Y = A(I,MM1)
            IF (Y .EQ. 0.0E0) GO TO 160
            Y = Y / X
            A(I,MM1) = Y
C
            DO 140 J = M, N
  140       A(I,J) = A(I,J) - Y * A(M,J)
C
            DO 150 J = 1, IGH
  150       A(J,M) = A(J,M) + Y * A(J,I)
C
  160    CONTINUE
C
  180 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE ELTRAN(NM,N,LOW,IGH,A,INT,Z)
C
      INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
      REAL A(NM,IGH),Z(NM,N)
      INTEGER INT(IGH)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ELMTRANS,
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     THIS SUBROUTINE ACCUMULATES THE STABILIZED ELEMENTARY
C     SIMILARITY TRANSFORMATIONS USED IN THE REDUCTION OF A
C     REAL GENERAL MATRIX TO UPPER HESSENBERG FORM BY  ELMHES.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        A CONTAINS THE MULTIPLIERS WHICH WERE USED IN THE
C          REDUCTION BY  ELMHES  IN ITS LOWER TRIANGLE
C          BELOW THE SUBDIAGONAL.
C
C        INT CONTAINS INFORMATION ON THE ROWS AND COLUMNS
C          INTERCHANGED IN THE REDUCTION BY  ELMHES.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     ON OUTPUT
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  ELMHES.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
C     .......... INITIALIZE Z TO IDENTITY MATRIX ..........
      DO 80 J = 1, N
C
         DO 60 I = 1, N
   60    Z(I,J) = 0.0E0
C
         Z(J,J) = 1.0E0
   80 CONTINUE
C
      KL = IGH - LOW - 1
      IF (KL .LT. 1) GO TO 200
C     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 140 MM = 1, KL
         MP = IGH - MM
         MP1 = MP + 1
C
         DO 100 I = MP1, IGH
  100    Z(I,MP) = A(I,MP-1)
C
         I = INT(MP)
         IF (I .EQ. MP) GO TO 140
C
         DO 130 J = MP, IGH
            Z(MP,J) = Z(I,J)
            Z(I,J) = 0.0E0
  130    CONTINUE
C
         Z(I,MP) = 1.0E0
  140 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE FIGI(NM,N,T,D,E,E2,IERR)
C
      INTEGER I,N,NM,IERR
      REAL T(NM,3),D(N),E(N),E2(N)
C
C     GIVEN A NONSYMMETRIC TRIDIAGONAL MATRIX SUCH THAT THE PRODUCTS
C     OF CORRESPONDING PAIRS OF OFF-DIAGONAL ELEMENTS ARE ALL
C     NON-NEGATIVE, THIS SUBROUTINE REDUCES IT TO A SYMMETRIC
C     TRIDIAGONAL MATRIX WITH THE SAME EIGENVALUES.  IF, FURTHER,
C     A ZERO PRODUCT ONLY OCCURS WHEN BOTH FACTORS ARE ZERO,
C     THE REDUCED MATRIX IS SIMILAR TO THE ORIGINAL MATRIX.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        T CONTAINS THE INPUT MATRIX.  ITS SUBDIAGONAL IS
C          STORED IN THE LAST N-1 POSITIONS OF THE FIRST COLUMN,
C          ITS DIAGONAL IN THE N POSITIONS OF THE SECOND COLUMN,
C          AND ITS SUPERDIAGONAL IN THE FIRST N-1 POSITIONS OF
C          THE THIRD COLUMN.  T(1,1) AND T(N,3) ARE ARBITRARY.
C
C     ON OUTPUT
C
C        T IS UNALTERED.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE SYMMETRIC MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE SYMMETRIC
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS NOT SET.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          N+I        IF T(I,1)*T(I-1,3) IS NEGATIVE,
C          -(3*N+I)   IF T(I,1)*T(I-1,3) IS ZERO WITH ONE FACTOR
C                     NON-ZERO.  IN THIS CASE, THE EIGENVECTORS OF
C                     THE SYMMETRIC MATRIX ARE NOT SIMPLY RELATED
C                     TO THOSE OF  T  AND SHOULD NOT BE SOUGHT.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
C
      DO 100 I = 1, N
         IF (I .EQ. 1) GO TO 90
         E2(I) = T(I,1) * T(I-1,3)
         IF (E2(I)) 1000, 60, 80
   60    IF (T(I,1) .EQ. 0.0E0 .AND. T(I-1,3) .EQ. 0.0E0) GO TO 80
C     .......... SET ERROR -- PRODUCT OF SOME PAIR OF OFF-DIAGONAL
C                ELEMENTS IS ZERO WITH ONE MEMBER NON-ZERO ..........
         IERR = -(3 * N + I)
   80    E(I) = SQRT(E2(I))
   90    D(I) = T(I,2)
  100 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- PRODUCT OF SOME PAIR OF OFF-DIAGONAL
C                ELEMENTS IS NEGATIVE ..........
 1000 IERR = N + I
 1001 RETURN
      END
      SUBROUTINE FIGI2(NM,N,T,D,E,Z,IERR)
C
      INTEGER I,J,N,NM,IERR
      REAL T(NM,3),D(N),E(N),Z(NM,N)
      REAL H
C
C     GIVEN A NONSYMMETRIC TRIDIAGONAL MATRIX SUCH THAT THE PRODUCTS
C     OF CORRESPONDING PAIRS OF OFF-DIAGONAL ELEMENTS ARE ALL
C     NON-NEGATIVE, AND ZERO ONLY WHEN BOTH FACTORS ARE ZERO, THIS
C     SUBROUTINE REDUCES IT TO A SYMMETRIC TRIDIAGONAL MATRIX
C     USING AND ACCUMULATING DIAGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        T CONTAINS THE INPUT MATRIX.  ITS SUBDIAGONAL IS
C          STORED IN THE LAST N-1 POSITIONS OF THE FIRST COLUMN,
C          ITS DIAGONAL IN THE N POSITIONS OF THE SECOND COLUMN,
C          AND ITS SUPERDIAGONAL IN THE FIRST N-1 POSITIONS OF
C          THE THIRD COLUMN.  T(1,1) AND T(N,3) ARE ARBITRARY.
C
C     ON OUTPUT
C
C        T IS UNALTERED.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE SYMMETRIC MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE SYMMETRIC
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS NOT SET.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN
C          THE REDUCTION.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          N+I        IF T(I,1)*T(I-1,3) IS NEGATIVE,
C          2*N+I      IF T(I,1)*T(I-1,3) IS ZERO WITH
C                     ONE FACTOR NON-ZERO.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
C
      DO 100 I = 1, N
C
         DO 50 J = 1, N
   50    Z(I,J) = 0.0E0
C
         IF (I .EQ. 1) GO TO 70
         H = T(I,1) * T(I-1,3)
         IF (H) 900, 60, 80
   60    IF (T(I,1) .NE. 0.0E0 .OR. T(I-1,3) .NE. 0.0E0) GO TO 1000
         E(I) = 0.0E0
   70    Z(I,I) = 1.0E0
         GO TO 90
   80    E(I) = SQRT(H)
         Z(I,I) = Z(I-1,I-1) * E(I) / T(I-1,3)
   90    D(I) = T(I,2)
  100 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- PRODUCT OF SOME PAIR OF OFF-DIAGONAL
C                ELEMENTS IS NEGATIVE ..........
  900 IERR = N + I
      GO TO 1001
C     .......... SET ERROR -- PRODUCT OF SOME PAIR OF OFF-DIAGONAL
C                ELEMENTS IS ZERO WITH ONE MEMBER NON-ZERO ..........
 1000 IERR = 2 * N + I
 1001 RETURN
      END
      SUBROUTINE HQR(NM,N,LOW,IGH,H,WR,WI,IERR)
C
      INTEGER I,J,K,L,M,N,EN,LL,MM,NA,NM,IGH,ITN,ITS,LOW,MP2,ENM2,IERR
      REAL H(NM,N),WR(N),WI(N)
      REAL P,Q,R,S,T,W,X,Y,ZZ,NORM,TST1,TST2
      LOGICAL NOTLAS
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR,
C     NUM. MATH. 14, 219-231(1970) BY MARTIN, PETERS, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 359-371(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A REAL
C     UPPER HESSENBERG MATRIX BY THE QR METHOD.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        H CONTAINS THE UPPER HESSENBERG MATRIX.  INFORMATION ABOUT
C          THE TRANSFORMATIONS USED IN THE REDUCTION TO HESSENBERG
C          FORM BY  ELMHES  OR  ORTHES, IF PERFORMED, IS STORED
C          IN THE REMAINING TRIANGLE UNDER THE HESSENBERG MATRIX.
C
C     ON OUTPUT
C
C        H HAS BEEN DESTROYED.  THEREFORE, IT MUST BE SAVED
C          BEFORE CALLING  HQR  IF SUBSEQUENT CALCULATION AND
C          BACK TRANSFORMATION OF EIGENVECTORS IS TO BE PERFORMED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES
C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS
C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE
C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      NORM = 0.0E0
      K = 1
C     .......... STORE ROOTS ISOLATED BY BALANC
C                AND COMPUTE MATRIX NORM ..........
      DO 50 I = 1, N
C
         DO 40 J = K, N
   40    NORM = NORM + ABS(H(I,J))
C
         K = I
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50
         WR(I) = H(I,I)
         WI(I) = 0.0E0
   50 CONTINUE
C
      EN = IGH
      T = 0.0E0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUES ..........
   60 IF (EN .LT. LOW) GO TO 1001
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
   70 DO 80 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 100
         S = ABS(H(L-1,L-1)) + ABS(H(L,L))
         IF (S .EQ. 0.0E0) S = NORM
         TST1 = S
         TST2 = TST1 + ABS(H(L,L-1))
         IF (TST2 .EQ. TST1) GO TO 100
   80 CONTINUE
C     .......... FORM SHIFT ..........
  100 X = H(EN,EN)
      IF (L .EQ. EN) GO TO 270
      Y = H(NA,NA)
      W = H(EN,NA) * H(NA,EN)
      IF (L .EQ. NA) GO TO 280
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130
C     .......... FORM EXCEPTIONAL SHIFT ..........
      T = T + X
C
      DO 120 I = LOW, EN
  120 H(I,I) = H(I,I) - X
C
      S = ABS(H(EN,NA)) + ABS(H(NA,ENM2))
      X = 0.75E0 * S
      Y = X
      W = -0.4375E0 * S * S
  130 ITS = ITS + 1
      ITN = ITN - 1
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS.
C                FOR M=EN-2 STEP -1 UNTIL L DO -- ..........
      DO 140 MM = L, ENM2
         M = ENM2 + L - MM
         ZZ = H(M,M)
         R = X - ZZ
         S = Y - ZZ
         P = (R * S - W) / H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - ZZ - R - S
         R = H(M+2,M+1)
         S = ABS(P) + ABS(Q) + ABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF (M .EQ. L) GO TO 150
         TST1 = ABS(P)*(ABS(H(M-1,M-1)) + ABS(ZZ) + ABS(H(M+1,M+1)))
         TST2 = TST1 + ABS(H(M,M-1))*(ABS(Q) + ABS(R))
         IF (TST2 .EQ. TST1) GO TO 150
  140 CONTINUE
C
  150 MP2 = M + 2
C
      DO 160 I = MP2, EN
         H(I,I-2) = 0.0E0
         IF (I .EQ. MP2) GO TO 160
         H(I,I-3) = 0.0E0
  160 CONTINUE
C     .......... DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C                COLUMNS M TO EN ..........
      DO 260 K = M, NA
         NOTLAS = K .NE. NA
         IF (K .EQ. M) GO TO 170
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.0E0
         IF (NOTLAS) R = H(K+2,K-1)
         X = ABS(P) + ABS(Q) + ABS(R)
         IF (X .EQ. 0.0E0) GO TO 260
         P = P / X
         Q = Q / X
         R = R / X
  170    S = SIGN(SQRT(P*P+Q*Q+R*R),P)
         IF (K .EQ. M) GO TO 180
         H(K,K-1) = -S * X
         GO TO 190
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)
  190    P = P + S
         X = P / S
         Y = Q / S
         ZZ = R / S
         Q = Q / P
         R = R / P
         IF (NOTLAS) GO TO 225
C     .......... ROW MODIFICATION ..........
         DO 200 J = K, N
            P = H(K,J) + Q * H(K+1,J)
            H(K,J) = H(K,J) - P * X
            H(K+1,J) = H(K+1,J) - P * Y
  200    CONTINUE
C
         J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
         DO 210 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1)
            H(I,K) = H(I,K) - P
            H(I,K+1) = H(I,K+1) - P * Q
  210    CONTINUE
         GO TO 255
  225    CONTINUE
C     .......... ROW MODIFICATION ..........
         DO 230 J = K, N
            P = H(K,J) + Q * H(K+1,J) + R * H(K+2,J)
            H(K,J) = H(K,J) - P * X
            H(K+1,J) = H(K+1,J) - P * Y
            H(K+2,J) = H(K+2,J) - P * ZZ
  230    CONTINUE
C
         J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
         DO 240 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1) + ZZ * H(I,K+2)
            H(I,K) = H(I,K) - P
            H(I,K+1) = H(I,K+1) - P * Q
            H(I,K+2) = H(I,K+2) - P * R
  240    CONTINUE
  255    CONTINUE
C
  260 CONTINUE
C
      GO TO 70
C     .......... ONE ROOT FOUND ..........
  270 WR(EN) = X + T
      WI(EN) = 0.0E0
      EN = NA
      GO TO 60
C     .......... TWO ROOTS FOUND ..........
  280 P = (Y - X) / 2.0E0
      Q = P * P + W
      ZZ = SQRT(ABS(Q))
      X = X + T
      IF (Q .LT. 0.0E0) GO TO 320
C     .......... REAL PAIR ..........
      ZZ = P + SIGN(ZZ,P)
      WR(NA) = X + ZZ
      WR(EN) = WR(NA)
      IF (ZZ .NE. 0.0E0) WR(EN) = X - W / ZZ
      WI(NA) = 0.0E0
      WI(EN) = 0.0E0
      GO TO 330
C     .......... COMPLEX PAIR ..........
  320 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = ZZ
      WI(EN) = -ZZ
  330 EN = ENM2
      GO TO 60
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END
      SUBROUTINE HQR2(NM,N,LOW,IGH,H,WR,WI,Z,IERR)
C
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN,
     X        IGH,ITN,ITS,LOW,MP2,ENM2,IERR
      REAL H(NM,N),WR(N),WI(N),Z(NM,N)
      REAL P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,TST1,TST2
      LOGICAL NOTLAS
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR2,
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A REAL UPPER HESSENBERG MATRIX BY THE QR METHOD.  THE
C     EIGENVECTORS OF A REAL GENERAL MATRIX CAN ALSO BE FOUND
C     IF  ELMHES  AND  ELTRAN  OR  ORTHES  AND  ORTRAN  HAVE
C     BEEN USED TO REDUCE THIS GENERAL MATRIX TO HESSENBERG FORM
C     AND TO ACCUMULATE THE SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        H CONTAINS THE UPPER HESSENBERG MATRIX.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED BY  ELTRAN
C          AFTER THE REDUCTION BY  ELMHES, OR BY  ORTRAN  AFTER THE
C          REDUCTION BY  ORTHES, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE HESSENBERG MATRIX ARE DESIRED, Z MUST CONTAIN THE
C          IDENTITY MATRIX.
C
C     ON OUTPUT
C
C        H HAS BEEN DESTROYED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES
C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS
C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE
C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
C          IF THE I-TH EIGENVALUE IS REAL, THE I-TH COLUMN OF Z
C          CONTAINS ITS EIGENVECTOR.  IF THE I-TH EIGENVALUE IS COMPLEX
C          WITH POSITIVE IMAGINARY PART, THE I-TH AND (I+1)-TH
C          COLUMNS OF Z CONTAIN THE REAL AND IMAGINARY PARTS OF ITS
C          EIGENVECTOR.  THE EIGENVECTORS ARE UNNORMALIZED.  IF AN
C          ERROR EXIT IS MADE, NONE OF THE EIGENVECTORS HAS BEEN FOUND.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      NORM = 0.0E0
      K = 1
C     .......... STORE ROOTS ISOLATED BY BALANC
C                AND COMPUTE MATRIX NORM ..........
      DO 50 I = 1, N
C
         DO 40 J = K, N
   40    NORM = NORM + ABS(H(I,J))
C
         K = I
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50
         WR(I) = H(I,I)
         WI(I) = 0.0E0
   50 CONTINUE
C
      EN = IGH
      T = 0.0E0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUES ..........
   60 IF (EN .LT. LOW) GO TO 340
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
   70 DO 80 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 100
         S = ABS(H(L-1,L-1)) + ABS(H(L,L))
         IF (S .EQ. 0.0E0) S = NORM
         TST1 = S
         TST2 = TST1 + ABS(H(L,L-1))
         IF (TST2 .EQ. TST1) GO TO 100
   80 CONTINUE
C     .......... FORM SHIFT ..........
  100 X = H(EN,EN)
      IF (L .EQ. EN) GO TO 270
      Y = H(NA,NA)
      W = H(EN,NA) * H(NA,EN)
      IF (L .EQ. NA) GO TO 280
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130
C     .......... FORM EXCEPTIONAL SHIFT ..........
      T = T + X
C
      DO 120 I = LOW, EN
  120 H(I,I) = H(I,I) - X
C
      S = ABS(H(EN,NA)) + ABS(H(NA,ENM2))
      X = 0.75E0 * S
      Y = X
      W = -0.4375E0 * S * S
  130 ITS = ITS + 1
      ITN = ITN - 1
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS.
C                FOR M=EN-2 STEP -1 UNTIL L DO -- ..........
      DO 140 MM = L, ENM2
         M = ENM2 + L - MM
         ZZ = H(M,M)
         R = X - ZZ
         S = Y - ZZ
         P = (R * S - W) / H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - ZZ - R - S
         R = H(M+2,M+1)
         S = ABS(P) + ABS(Q) + ABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF (M .EQ. L) GO TO 150
         TST1 = ABS(P)*(ABS(H(M-1,M-1)) + ABS(ZZ) + ABS(H(M+1,M+1)))
         TST2 = TST1 + ABS(H(M,M-1))*(ABS(Q) + ABS(R))
         IF (TST2 .EQ. TST1) GO TO 150
  140 CONTINUE
C
  150 MP2 = M + 2
C
      DO 160 I = MP2, EN
         H(I,I-2) = 0.0E0
         IF (I .EQ. MP2) GO TO 160
         H(I,I-3) = 0.0E0
  160 CONTINUE
C     .......... DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C                COLUMNS M TO EN ..........
      DO 260 K = M, NA
         NOTLAS = K .NE. NA
         IF (K .EQ. M) GO TO 170
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.0E0
         IF (NOTLAS) R = H(K+2,K-1)
         X = ABS(P) + ABS(Q) + ABS(R)
         IF (X .EQ. 0.0E0) GO TO 260
         P = P / X
         Q = Q / X
         R = R / X
  170    S = SIGN(SQRT(P*P+Q*Q+R*R),P)
         IF (K .EQ. M) GO TO 180
         H(K,K-1) = -S * X
         GO TO 190
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)
  190    P = P + S
         X = P / S
         Y = Q / S
         ZZ = R / S
         Q = Q / P
         R = R / P
         IF (NOTLAS) GO TO 225
C     .......... ROW MODIFICATION ..........
         DO 200 J = K, N
            P = H(K,J) + Q * H(K+1,J)
            H(K,J) = H(K,J) - P * X
            H(K+1,J) = H(K+1,J) - P * Y
  200    CONTINUE
C
         J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
         DO 210 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1)
            H(I,K) = H(I,K) - P
            H(I,K+1) = H(I,K+1) - P * Q
  210    CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
         DO 220 I = LOW, IGH
            P = X * Z(I,K) + Y * Z(I,K+1)
            Z(I,K) = Z(I,K) - P
            Z(I,K+1) = Z(I,K+1) - P * Q
  220    CONTINUE
         GO TO 255
  225    CONTINUE
C     .......... ROW MODIFICATION ..........
         DO 230 J = K, N
            P = H(K,J) + Q * H(K+1,J) + R * H(K+2,J)
            H(K,J) = H(K,J) - P * X
            H(K+1,J) = H(K+1,J) - P * Y
            H(K+2,J) = H(K+2,J) - P * ZZ
  230    CONTINUE
C
         J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
         DO 240 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1) + ZZ * H(I,K+2)
            H(I,K) = H(I,K) - P
            H(I,K+1) = H(I,K+1) - P * Q
            H(I,K+2) = H(I,K+2) - P * R
  240    CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
         DO 250 I = LOW, IGH
            P = X * Z(I,K) + Y * Z(I,K+1) + ZZ * Z(I,K+2)
            Z(I,K) = Z(I,K) - P
            Z(I,K+1) = Z(I,K+1) - P * Q
            Z(I,K+2) = Z(I,K+2) - P * R
  250    CONTINUE
  255    CONTINUE
C
  260 CONTINUE
C
      GO TO 70
C     .......... ONE ROOT FOUND ..........
  270 H(EN,EN) = X + T
      WR(EN) = H(EN,EN)
      WI(EN) = 0.0E0
      EN = NA
      GO TO 60
C     .......... TWO ROOTS FOUND ..........
  280 P = (Y - X) / 2.0E0
      Q = P * P + W
      ZZ = SQRT(ABS(Q))
      H(EN,EN) = X + T
      X = H(EN,EN)
      H(NA,NA) = Y + T
      IF (Q .LT. 0.0E0) GO TO 320
C     .......... REAL PAIR ..........
      ZZ = P + SIGN(ZZ,P)
      WR(NA) = X + ZZ
      WR(EN) = WR(NA)
      IF (ZZ .NE. 0.0E0) WR(EN) = X - W / ZZ
      WI(NA) = 0.0E0
      WI(EN) = 0.0E0
      X = H(EN,NA)
      S = ABS(X) + ABS(ZZ)
      P = X / S
      Q = ZZ / S
      R = SQRT(P*P+Q*Q)
      P = P / R
      Q = Q / R
C     .......... ROW MODIFICATION ..........
      DO 290 J = NA, N
         ZZ = H(NA,J)
         H(NA,J) = Q * ZZ + P * H(EN,J)
         H(EN,J) = Q * H(EN,J) - P * ZZ
  290 CONTINUE
C     .......... COLUMN MODIFICATION ..........
      DO 300 I = 1, EN
         ZZ = H(I,NA)
         H(I,NA) = Q * ZZ + P * H(I,EN)
         H(I,EN) = Q * H(I,EN) - P * ZZ
  300 CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
      DO 310 I = LOW, IGH
         ZZ = Z(I,NA)
         Z(I,NA) = Q * ZZ + P * Z(I,EN)
         Z(I,EN) = Q * Z(I,EN) - P * ZZ
  310 CONTINUE
C
      GO TO 330
C     .......... COMPLEX PAIR ..........
  320 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = ZZ
      WI(EN) = -ZZ
  330 EN = ENM2
      GO TO 60
C     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM ..........
  340 IF (NORM .EQ. 0.0E0) GO TO 1001
C     .......... FOR EN=N STEP -1 UNTIL 1 DO -- ..........
      DO 800 NN = 1, N
         EN = N + 1 - NN
         P = WR(EN)
         Q = WI(EN)
         NA = EN - 1
         IF (Q) 710, 600, 800
C     .......... REAL VECTOR ..........
  600    M = EN
         H(EN,EN) = 1.0E0
         IF (NA .EQ. 0) GO TO 800
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
         DO 700 II = 1, NA
            I = EN - II
            W = H(I,I) - P
            R = 0.0E0
C
            DO 610 J = M, EN
  610       R = R + H(I,J) * H(J,EN)
C
            IF (WI(I) .GE. 0.0E0) GO TO 630
            ZZ = W
            S = R
            GO TO 700
  630       M = I
            IF (WI(I) .NE. 0.0E0) GO TO 640
            T = W
            IF (T .NE. 0.0E0) GO TO 635
               TST1 = NORM
               T = TST1
  632          T = 0.01E0 * T
               TST2 = NORM + T
               IF (TST2 .GT. TST1) GO TO 632
  635       H(I,EN) = -R / T
            GO TO 680
C     .......... SOLVE REAL EQUATIONS ..........
  640       X = H(I,I+1)
            Y = H(I+1,I)
            Q = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I)
            T = (X * S - ZZ * R) / Q
            H(I,EN) = T
            IF (ABS(X) .LE. ABS(ZZ)) GO TO 650
            H(I+1,EN) = (-R - W * T) / X
            GO TO 680
  650       H(I+1,EN) = (-S - Y * T) / ZZ
C
C     .......... OVERFLOW CONTROL ..........
  680       T = ABS(H(I,EN))
            IF (T .EQ. 0.0E0) GO TO 700
            TST1 = T
            TST2 = TST1 + 1.0E0/TST1
            IF (TST2 .GT. TST1) GO TO 700
            DO 690 J = I, EN
               H(J,EN) = H(J,EN)/T
  690       CONTINUE
C
  700    CONTINUE
C     .......... END REAL VECTOR ..........
         GO TO 800
C     .......... COMPLEX VECTOR ..........
  710    M = NA
C     .......... LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
C                EIGENVECTOR MATRIX IS TRIANGULAR ..........
         IF (ABS(H(EN,NA)) .LE. ABS(H(NA,EN))) GO TO 720
         H(NA,NA) = Q / H(EN,NA)
         H(NA,EN) = -(H(EN,EN) - P) / H(EN,NA)
         GO TO 730
  720    CALL CDIV(0.0E0,-H(NA,EN),H(NA,NA)-P,Q,H(NA,NA),H(NA,EN))
  730    H(EN,NA) = 0.0E0
         H(EN,EN) = 1.0E0
         ENM2 = NA - 1
         IF (ENM2 .EQ. 0) GO TO 800
C     .......... FOR I=EN-2 STEP -1 UNTIL 1 DO -- ..........
         DO 795 II = 1, ENM2
            I = NA - II
            W = H(I,I) - P
            RA = 0.0E0
            SA = 0.0E0
C
            DO 760 J = M, EN
               RA = RA + H(I,J) * H(J,NA)
               SA = SA + H(I,J) * H(J,EN)
  760       CONTINUE
C
            IF (WI(I) .GE. 0.0E0) GO TO 770
            ZZ = W
            R = RA
            S = SA
            GO TO 795
  770       M = I
            IF (WI(I) .NE. 0.0E0) GO TO 780
            CALL CDIV(-RA,-SA,W,Q,H(I,NA),H(I,EN))
            GO TO 790
C     .......... SOLVE COMPLEX EQUATIONS ..........
  780       X = H(I,I+1)
            Y = H(I+1,I)
            VR = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I) - Q * Q
            VI = (WR(I) - P) * 2.0E0 * Q
            IF (VR .NE. 0.0E0 .OR. VI .NE. 0.0E0) GO TO 784
               TST1 = NORM * (ABS(W) + ABS(Q) + ABS(X)
     X                      + ABS(Y) + ABS(ZZ))
               VR = TST1
  783          VR = 0.01E0 * VR
               TST2 = TST1 + VR
               IF (TST2 .GT. TST1) GO TO 783
  784       CALL CDIV(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA,VR,VI,
     X                H(I,NA),H(I,EN))
            IF (ABS(X) .LE. ABS(ZZ) + ABS(Q)) GO TO 785
            H(I+1,NA) = (-RA - W * H(I,NA) + Q * H(I,EN)) / X
            H(I+1,EN) = (-SA - W * H(I,EN) - Q * H(I,NA)) / X
            GO TO 790
  785       CALL CDIV(-R-Y*H(I,NA),-S-Y*H(I,EN),ZZ,Q,
     X                H(I+1,NA),H(I+1,EN))
C
C     .......... OVERFLOW CONTROL ..........
  790       T = AMAX1(ABS(H(I,NA)), ABS(H(I,EN)))
            IF (T .EQ. 0.0E0) GO TO 795
            TST1 = T
            TST2 = TST1 + 1.0E0/TST1
            IF (TST2 .GT. TST1) GO TO 795
            DO 792 J = I, EN
               H(J,NA) = H(J,NA)/T
               H(J,EN) = H(J,EN)/T
  792       CONTINUE
C
  795    CONTINUE
C     .......... END COMPLEX VECTOR ..........
  800 CONTINUE
C     .......... END BACK SUBSTITUTION.
C                VECTORS OF ISOLATED ROOTS ..........
      DO 840 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840
C
         DO 820 J = I, N
  820    Z(I,J) = H(I,J)
C
  840 CONTINUE
C     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW DO -- ..........
      DO 880 JJ = LOW, N
         J = N + LOW - JJ
         M = MIN0(J,IGH)
C
         DO 880 I = LOW, IGH
            ZZ = 0.0E0
C
            DO 860 K = LOW, M
  860       ZZ = ZZ + Z(I,K) * H(K,J)
C
            Z(I,J) = ZZ
  880 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END
      SUBROUTINE HTRIB3(NM,N,A,TAU,M,ZR,ZI)
C
      INTEGER I,J,K,L,M,N,NM
      REAL A(NM,N),TAU(2,N),ZR(NM,M),ZI(NM,M)
      REAL H,S,SI
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE TRBAK3, NUM. MATH. 11, 181-195(1968)
C     BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX HERMITIAN
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     REAL SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  HTRID3.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS INFORMATION ABOUT THE UNITARY TRANSFORMATIONS
C          USED IN THE REDUCTION BY  HTRID3.
C
C        TAU CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
C
C        ZR CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE TRANSFORMED EIGENVECTORS
C          IN THEIR FIRST M COLUMNS.
C
C     NOTE THAT THE LAST COMPONENT OF EACH RETURNED VECTOR
C     IS REAL AND THAT VECTOR EUCLIDEAN NORMS ARE PRESERVED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
C     .......... TRANSFORM THE EIGENVECTORS OF THE REAL SYMMETRIC
C                TRIDIAGONAL MATRIX TO THOSE OF THE HERMITIAN
C                TRIDIAGONAL MATRIX. ..........
      DO 50 K = 1, N
C
         DO 50 J = 1, M
            ZI(K,J) = -ZR(K,J) * TAU(2,K)
            ZR(K,J) = ZR(K,J) * TAU(1,K)
   50 CONTINUE
C
      IF (N .EQ. 1) GO TO 200
C     .......... RECOVER AND APPLY THE HOUSEHOLDER MATRICES ..........
      DO 140 I = 2, N
         L = I - 1
         H = A(I,I)
         IF (H .EQ. 0.0E0) GO TO 140
C
         DO 130 J = 1, M
            S = 0.0E0
            SI = 0.0E0
C
            DO 110 K = 1, L
               S = S + A(I,K) * ZR(K,J) - A(K,I) * ZI(K,J)
               SI = SI + A(I,K) * ZI(K,J) + A(K,I) * ZR(K,J)
  110       CONTINUE
C     .......... DOUBLE DIVISIONS AVOID POSSIBLE UNDERFLOW ..........
            S = (S / H) / H
            SI = (SI / H) / H
C
            DO 120 K = 1, L
               ZR(K,J) = ZR(K,J) - S * A(I,K) - SI * A(K,I)
               ZI(K,J) = ZI(K,J) - SI * A(I,K) + S * A(K,I)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE HTRIBK(NM,N,AR,AI,TAU,M,ZR,ZI)
C
      INTEGER I,J,K,L,M,N,NM
      REAL AR(NM,N),AI(NM,N),TAU(2,N),ZR(NM,M),ZI(NM,M)
      REAL H,S,SI
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE TRBAK1, NUM. MATH. 11, 181-195(1968)
C     BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX HERMITIAN
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     REAL SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  HTRIDI.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        AR AND AI CONTAIN INFORMATION ABOUT THE UNITARY TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  HTRIDI  IN THEIR
C          FULL LOWER TRIANGLES EXCEPT FOR THE DIAGONAL OF AR.
C
C        TAU CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
C
C        ZR CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE TRANSFORMED EIGENVECTORS
C          IN THEIR FIRST M COLUMNS.
C
C     NOTE THAT THE LAST COMPONENT OF EACH RETURNED VECTOR
C     IS REAL AND THAT VECTOR EUCLIDEAN NORMS ARE PRESERVED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
C     .......... TRANSFORM THE EIGENVECTORS OF THE REAL SYMMETRIC
C                TRIDIAGONAL MATRIX TO THOSE OF THE HERMITIAN
C                TRIDIAGONAL MATRIX. ..........
      DO 50 K = 1, N
C
         DO 50 J = 1, M
            ZI(K,J) = -ZR(K,J) * TAU(2,K)
            ZR(K,J) = ZR(K,J) * TAU(1,K)
   50 CONTINUE
C
      IF (N .EQ. 1) GO TO 200
C     .......... RECOVER AND APPLY THE HOUSEHOLDER MATRICES ..........
      DO 140 I = 2, N
         L = I - 1
         H = AI(I,I)
         IF (H .EQ. 0.0E0) GO TO 140
C
         DO 130 J = 1, M
            S = 0.0E0
            SI = 0.0E0
C
            DO 110 K = 1, L
               S = S + AR(I,K) * ZR(K,J) - AI(I,K) * ZI(K,J)
               SI = SI + AR(I,K) * ZI(K,J) + AI(I,K) * ZR(K,J)
  110       CONTINUE
C     .......... DOUBLE DIVISIONS AVOID POSSIBLE UNDERFLOW ..........
            S = (S / H) / H
            SI = (SI / H) / H
C
            DO 120 K = 1, L
               ZR(K,J) = ZR(K,J) - S * AR(I,K) - SI * AI(I,K)
               ZI(K,J) = ZI(K,J) - SI * AR(I,K) + S * AI(I,K)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE HTRID3(NM,N,A,D,E,E2,TAU)
C
      INTEGER I,J,K,L,N,II,NM,JM1,JP1
      REAL A(NM,N),D(N),E(N),E2(N),TAU(2,N)
      REAL F,G,H,FI,GI,HH,SI,SCALE,PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE TRED3, NUM. MATH. 11, 181-195(1968)
C     BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A COMPLEX HERMITIAN MATRIX, STORED AS
C     A SINGLE SQUARE ARRAY, TO A REAL SYMMETRIC TRIDIAGONAL MATRIX
C     USING UNITARY SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS THE LOWER TRIANGLE OF THE COMPLEX HERMITIAN INPUT
C          MATRIX.  THE REAL PARTS OF THE MATRIX ELEMENTS ARE STORED
C          IN THE FULL LOWER TRIANGLE OF A, AND THE IMAGINARY PARTS
C          ARE STORED IN THE TRANSPOSED POSITIONS OF THE STRICT UPPER
C          TRIANGLE OF A.  NO STORAGE IS REQUIRED FOR THE ZERO
C          IMAGINARY PARTS OF THE DIAGONAL ELEMENTS.
C
C     ON OUTPUT
C
C        A CONTAINS INFORMATION ABOUT THE UNITARY TRANSFORMATIONS
C          USED IN THE REDUCTION.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE THE TRIDIAGONAL MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C        TAU CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      TAU(1,N) = 1.0E0
      TAU(2,N) = 0.0E0
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0E0
         SCALE = 0.0E0
         IF (L .LT. 1) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(A(I,K)) + ABS(A(K,I))
C
         IF (SCALE .NE. 0.0E0) GO TO 140
         TAU(1,L) = 1.0E0
         TAU(2,L) = 0.0E0
  130    E(I) = 0.0E0
         E2(I) = 0.0E0
         GO TO 290
C
  140    DO 150 K = 1, L
            A(I,K) = A(I,K) / SCALE
            A(K,I) = A(K,I) / SCALE
            H = H + A(I,K) * A(I,K) + A(K,I) * A(K,I)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         G = SQRT(H)
         E(I) = SCALE * G
         F = PYTHAG(A(I,L),A(L,I))
C     .......... FORM NEXT DIAGONAL ELEMENT OF MATRIX T ..........
         IF (F .EQ. 0.0E0) GO TO 160
         TAU(1,L) = (A(L,I) * TAU(2,I) - A(I,L) * TAU(1,I)) / F
         SI = (A(I,L) * TAU(2,I) + A(L,I) * TAU(1,I)) / F
         H = H + F * G
         G = 1.0E0 + G / F
         A(I,L) = G * A(I,L)
         A(L,I) = G * A(L,I)
         IF (L .EQ. 1) GO TO 270
         GO TO 170
  160    TAU(1,L) = -TAU(1,I)
         SI = TAU(2,I)
         A(I,L) = G
  170    F = 0.0E0
C
         DO 240 J = 1, L
            G = 0.0E0
            GI = 0.0E0
            IF (J .EQ. 1) GO TO 190
            JM1 = J - 1
C     .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, JM1
               G = G + A(J,K) * A(I,K) + A(K,J) * A(K,I)
               GI = GI - A(J,K) * A(K,I) + A(K,J) * A(I,K)
  180       CONTINUE
C
  190       G = G + A(J,J) * A(I,J)
            GI = GI - A(J,J) * A(J,I)
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
               G = G + A(K,J) * A(I,K) - A(J,K) * A(K,I)
               GI = GI - A(K,J) * A(K,I) - A(J,K) * A(I,K)
  200       CONTINUE
C     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
            TAU(2,J) = GI / H
            F = F + E(J) * A(I,J) - TAU(2,J) * A(J,I)
  240    CONTINUE
C
         HH = F / (H + H)
C     .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = A(I,J)
            G = E(J) - HH * F
            E(J) = G
            FI = -A(J,I)
            GI = TAU(2,J) - HH * FI
            TAU(2,J) = -GI
            A(J,J) = A(J,J) - 2.0E0 * (F * G + FI * GI)
            IF (J .EQ. 1) GO TO 260
            JM1 = J - 1
C
            DO 250 K = 1, JM1
               A(J,K) = A(J,K) - F * E(K) - G * A(I,K)
     X                         + FI * TAU(2,K) + GI * A(K,I)
               A(K,J) = A(K,J) - F * TAU(2,K) - G * A(K,I)
     X                         - FI * E(K) - GI * A(I,K)
  250       CONTINUE
C
  260    CONTINUE
C
  270    DO 280 K = 1, L
            A(I,K) = SCALE * A(I,K)
            A(K,I) = SCALE * A(K,I)
  280    CONTINUE
C
         TAU(2,L) = -SI
  290    D(I) = A(I,I)
         A(I,I) = SCALE * SQRT(H)
  300 CONTINUE
C
      RETURN
      END
      SUBROUTINE HTRIDI(NM,N,AR,AI,D,E,E2,TAU)
C
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL AR(NM,N),AI(NM,N),D(N),E(N),E2(N),TAU(2,N)
      REAL F,G,H,FI,GI,HH,SI,SCALE,PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE TRED1, NUM. MATH. 11, 181-195(1968)
C     BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A COMPLEX HERMITIAN MATRIX
C     TO A REAL SYMMETRIC TRIDIAGONAL MATRIX USING
C     UNITARY SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX HERMITIAN INPUT MATRIX.
C          ONLY THE LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT
C
C        AR AND AI CONTAIN INFORMATION ABOUT THE UNITARY TRANS-
C          FORMATIONS USED IN THE REDUCTION IN THEIR FULL LOWER
C          TRIANGLES.  THEIR STRICT UPPER TRIANGLES AND THE
C          DIAGONAL OF AR ARE UNALTERED.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE THE TRIDIAGONAL MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C        TAU CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      TAU(1,N) = 1.0E0
      TAU(2,N) = 0.0E0
C
      DO 100 I = 1, N
  100 D(I) = AR(I,I)
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0E0
         SCALE = 0.0E0
         IF (L .LT. 1) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(AR(I,K)) + ABS(AI(I,K))
C
         IF (SCALE .NE. 0.0E0) GO TO 140
         TAU(1,L) = 1.0E0
         TAU(2,L) = 0.0E0
  130    E(I) = 0.0E0
         E2(I) = 0.0E0
         GO TO 290
C
  140    DO 150 K = 1, L
            AR(I,K) = AR(I,K) / SCALE
            AI(I,K) = AI(I,K) / SCALE
            H = H + AR(I,K) * AR(I,K) + AI(I,K) * AI(I,K)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         G = SQRT(H)
         E(I) = SCALE * G
         F = PYTHAG(AR(I,L),AI(I,L))
C     .......... FORM NEXT DIAGONAL ELEMENT OF MATRIX T ..........
         IF (F .EQ. 0.0E0) GO TO 160
         TAU(1,L) = (AI(I,L) * TAU(2,I) - AR(I,L) * TAU(1,I)) / F
         SI = (AR(I,L) * TAU(2,I) + AI(I,L) * TAU(1,I)) / F
         H = H + F * G
         G = 1.0E0 + G / F
         AR(I,L) = G * AR(I,L)
         AI(I,L) = G * AI(I,L)
         IF (L .EQ. 1) GO TO 270
         GO TO 170
  160    TAU(1,L) = -TAU(1,I)
         SI = TAU(2,I)
         AR(I,L) = G
  170    F = 0.0E0
C
         DO 240 J = 1, L
            G = 0.0E0
            GI = 0.0E0
C     .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, J
               G = G + AR(J,K) * AR(I,K) + AI(J,K) * AI(I,K)
               GI = GI - AR(J,K) * AI(I,K) + AI(J,K) * AR(I,K)
  180       CONTINUE
C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
               G = G + AR(K,J) * AR(I,K) - AI(K,J) * AI(I,K)
               GI = GI - AR(K,J) * AI(I,K) - AI(K,J) * AR(I,K)
  200       CONTINUE
C     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
            TAU(2,J) = GI / H
            F = F + E(J) * AR(I,J) - TAU(2,J) * AI(I,J)
  240    CONTINUE
C
         HH = F / (H + H)
C     .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = AR(I,J)
            G = E(J) - HH * F
            E(J) = G
            FI = -AI(I,J)
            GI = TAU(2,J) - HH * FI
            TAU(2,J) = -GI
C
            DO 260 K = 1, J
               AR(J,K) = AR(J,K) - F * E(K) - G * AR(I,K)
     X                           + FI * TAU(2,K) + GI * AI(I,K)
               AI(J,K) = AI(J,K) - F * TAU(2,K) - G * AI(I,K)
     X                           - FI * E(K) - GI * AR(I,K)
  260    CONTINUE
C
  270    DO 280 K = 1, L
            AR(I,K) = SCALE * AR(I,K)
            AI(I,K) = SCALE * AI(I,K)
  280    CONTINUE
C
         TAU(2,L) = -SI
  290    HH = D(I)
         D(I) = AR(I,I)
         AR(I,I) = HH
         AI(I,I) = SCALE * SQRT(H)
  300 CONTINUE
C
      RETURN
      END
      SUBROUTINE IMTQL1(N,D,E,IERR)
C
      INTEGER I,J,L,M,N,II,MML,IERR
      REAL D(N),E(N)
      REAL B,C,F,G,P,R,S,TST1,TST2,PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL1,
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE IMPLICIT QL METHOD.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES.
C
C        E HAS BEEN DESTROYED.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      E(N) = 0.0E0
C
      DO 290 L = 1, N
         J = 0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
            TST1 = ABS(D(M)) + ABS(D(M+1))
            TST2 = TST1 + ABS(E(M))
            IF (TST2 .EQ. TST1) GO TO 120
  110    CONTINUE
C
  120    P = D(L)
         IF (M .EQ. L) GO TO 215
         IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         G = (D(L+1) - P) / (2.0E0 * E(L))
         R = PYTHAG(G,1.0E0)
         G = D(M) - P + E(L) / (G + SIGN(R,G))
         S = 1.0E0
         C = 1.0E0
         P = 0.0E0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            R = PYTHAG(F,G)
            E(I+1) = R
            IF (R .EQ. 0.0E0) GO TO 210
            S = F / R
            C = G / R
            G = D(I+1) - P
            R = (D(I) - G) * S + 2.0E0 * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
  200    CONTINUE
C
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.0E0
         GO TO 105
C     .......... RECOVER FROM UNDERFLOW ..........
  210    D(I+1) = D(I+1) - P
         E(M) = 0.0E0
         GO TO 105
C     .......... ORDER EIGENVALUES ..........
  215    IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      SUBROUTINE IMTQL2(NM,N,D,E,Z,IERR)
C
      INTEGER I,J,K,L,M,N,II,NM,MML,IERR
      REAL D(N),E(N),Z(NM,N)
      REAL B,C,F,G,P,R,S,TST1,TST2,PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL2,
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1.
C
C        E HAS BEEN DESTROYED.
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      E(N) = 0.0E0
C
      DO 240 L = 1, N
         J = 0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
            TST1 = ABS(D(M)) + ABS(D(M+1))
            TST2 = TST1 + ABS(E(M))
            IF (TST2 .EQ. TST1) GO TO 120
  110    CONTINUE
C
  120    P = D(L)
         IF (M .EQ. L) GO TO 240
         IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         G = (D(L+1) - P) / (2.0E0 * E(L))
         R = PYTHAG(G,1.0E0)
         G = D(M) - P + E(L) / (G + SIGN(R,G))
         S = 1.0E0
         C = 1.0E0
         P = 0.0E0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            R = PYTHAG(F,G)
            E(I+1) = R
            IF (R .EQ. 0.0E0) GO TO 210
            S = F / R
            C = G / R
            G = D(I+1) - P
            R = (D(I) - G) * S + 2.0E0 * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
C     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               F = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * F
               Z(K,I) = C * Z(K,I) - S * F
  180       CONTINUE
C
  200    CONTINUE
C
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.0E0
         GO TO 105
C     .......... RECOVER FROM UNDERFLOW ..........
  210    D(I+1) = D(I+1) - P
         E(M) = 0.0E0
         GO TO 105
  240 CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      SUBROUTINE IMTQLV(N,D,E,E2,W,IND,IERR,RV1)
C
      INTEGER I,J,K,L,M,N,II,MML,TAG,IERR
      REAL D(N),E(N),E2(N),W(N),RV1(N)
      REAL B,C,F,G,P,R,S,TST1,TST2,PYTHAG
      INTEGER IND(N)
C
C     THIS SUBROUTINE IS A VARIANT OF  IMTQL1  WHICH IS A TRANSLATION OF
C     ALGOL PROCEDURE IMTQL1, NUM. MATH. 12, 377-383(1968) BY MARTIN AND
C     WILKINSON, AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC TRIDIAGONAL
C     MATRIX BY THE IMPLICIT QL METHOD AND ASSOCIATES WITH THEM
C     THEIR CORRESPONDING SUBMATRIX INDICES.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY.
C
C     ON OUTPUT
C
C        D AND E ARE UNALTERED.
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS ALSO SET TO ZERO.
C
C        W CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES.
C
C        IND CONTAINS THE SUBMATRIX INDICES ASSOCIATED WITH THE
C          CORRESPONDING EIGENVALUES IN W -- 1 FOR EIGENVALUES
C          BELONGING TO THE FIRST SUBMATRIX FROM THE TOP,
C          2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC..
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C        RV1 IS A TEMPORARY STORAGE ARRAY.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      K = 0
      TAG = 0
C
      DO 100 I = 1, N
         W(I) = D(I)
         IF (I .NE. 1) RV1(I-1) = E(I)
  100 CONTINUE
C
      E2(1) = 0.0E0
      RV1(N) = 0.0E0
C
      DO 290 L = 1, N
         J = 0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
            TST1 = ABS(W(M)) + ABS(W(M+1))
            TST2 = TST1 + ABS(RV1(M))
            IF (TST2 .EQ. TST1) GO TO 120
C     .......... GUARD AGAINST UNDERFLOWED ELEMENT OF E2 ..........
            IF (E2(M+1) .EQ. 0.0E0) GO TO 125
  110    CONTINUE
C
  120    IF (M .LE. K) GO TO 130
         IF (M .NE. N) E2(M+1) = 0.0E0
  125    K = M
         TAG = TAG + 1
  130    P = W(L)
         IF (M .EQ. L) GO TO 215
         IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         G = (W(L+1) - P) / (2.0E0 * RV1(L))
         R = PYTHAG(G,1.0E0)
         G = W(M) - P + RV1(L) / (G + SIGN(R,G))
         S = 1.0E0
         C = 1.0E0
         P = 0.0E0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            F = S * RV1(I)
            B = C * RV1(I)
            R = PYTHAG(F,G)
            RV1(I+1) = R
            IF (R .EQ. 0.0E0) GO TO 210
            S = F / R
            C = G / R
            G = W(I+1) - P
            R = (W(I) - G) * S + 2.0E0 * C * B
            P = S * R
            W(I+1) = G + P
            G = C * R - B
  200    CONTINUE
C
         W(L) = W(L) - P
         RV1(L) = G
         RV1(M) = 0.0E0
         GO TO 105
C     .......... RECOVER FROM UNDERFLOW ..........
  210    W(I+1) = W(I+1) - P
         RV1(M) = 0.0E0
         GO TO 105
C     .......... ORDER EIGENVALUES ..........
  215    IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. W(I-1)) GO TO 270
            W(I) = W(I-1)
            IND(I) = IND(I-1)
  230    CONTINUE
C
  250    I = 1
  270    W(I) = P
         IND(I) = TAG
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      SUBROUTINE INVIT(NM,N,A,WR,WI,SELECT,MM,M,Z,IERR,RM1,RV1,RV2)
C
      INTEGER I,J,K,L,M,N,S,II,IP,MM,MP,NM,NS,N1,UK,IP1,ITS,KM1,IERR
      REAL A(NM,N),WR(N),WI(N),Z(NM,MM),RM1(N,N),
     X       RV1(N),RV2(N)
      REAL T,W,X,Y,EPS3,NORM,NORMV,EPSLON,GROWTO,ILAMBD,
     X       PYTHAG,RLAMBD,UKROOT
      LOGICAL SELECT(N)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE INVIT
C     BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A REAL UPPER
C     HESSENBERG MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
C     USING INVERSE ITERATION.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS THE HESSENBERG MATRIX.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY,
C          OF THE EIGENVALUES OF THE MATRIX.  THE EIGENVALUES MUST BE
C          STORED IN A MANNER IDENTICAL TO THAT OF SUBROUTINE  HQR,
C          WHICH RECOGNIZES POSSIBLE SPLITTING OF THE MATRIX.
C
C        SELECT SPECIFIES THE EIGENVECTORS TO BE FOUND. THE
C          EIGENVECTOR CORRESPONDING TO THE J-TH EIGENVALUE IS
C          SPECIFIED BY SETTING SELECT(J) TO .TRUE..
C
C        MM SHOULD BE SET TO AN UPPER BOUND FOR THE NUMBER OF
C          COLUMNS REQUIRED TO STORE THE EIGENVECTORS TO BE FOUND.
C          NOTE THAT TWO COLUMNS ARE REQUIRED TO STORE THE
C          EIGENVECTOR CORRESPONDING TO A COMPLEX EIGENVALUE.
C
C     ON OUTPUT
C
C        A AND WI ARE UNALTERED.
C
C        WR MAY HAVE BEEN ALTERED SINCE CLOSE EIGENVALUES ARE PERTURBED
C          SLIGHTLY IN SEARCHING FOR INDEPENDENT EIGENVECTORS.
C
C        SELECT MAY HAVE BEEN ALTERED.  IF THE ELEMENTS CORRESPONDING
C          TO A PAIR OF CONJUGATE COMPLEX EIGENVALUES WERE EACH
C          INITIALLY SET TO .TRUE., THE PROGRAM RESETS THE SECOND OF
C          THE TWO ELEMENTS TO .FALSE..
C
C        M IS THE NUMBER OF COLUMNS ACTUALLY USED TO STORE
C          THE EIGENVECTORS.
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
C          IF THE NEXT SELECTED EIGENVALUE IS REAL, THE NEXT COLUMN
C          OF Z CONTAINS ITS EIGENVECTOR.  IF THE EIGENVALUE IS
C          COMPLEX, THE NEXT TWO COLUMNS OF Z CONTAIN THE REAL AND
C          IMAGINARY PARTS OF ITS EIGENVECTOR.  THE EIGENVECTORS ARE
C          NORMALIZED SO THAT THE COMPONENT OF LARGEST MAGNITUDE IS 1.
C          ANY VECTOR WHICH FAILS THE ACCEPTANCE TEST IS SET TO ZERO.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          -(2*N+1)   IF MORE THAN MM COLUMNS OF Z ARE NECESSARY
C                     TO STORE THE EIGENVECTORS CORRESPONDING TO
C                     THE SPECIFIED EIGENVALUES.
C          -K         IF THE ITERATION CORRESPONDING TO THE K-TH
C                     VALUE FAILS,
C          -(N+K)     IF BOTH ERROR SITUATIONS OCCUR.
C
C        RM1, RV1, AND RV2 ARE TEMPORARY STORAGE ARRAYS.  NOTE THAT RM1
C          IS SQUARE OF DIMENSION N BY N AND, AUGMENTED BY TWO COLUMNS
C          OF Z, IS THE TRANSPOSE OF THE CORRESPONDING ALGOL B ARRAY.
C
C     THE ALGOL PROCEDURE GUESSVEC APPEARS IN INVIT IN LINE.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      UK = 0
      S = 1
C     .......... IP = 0, REAL EIGENVALUE
C                     1, FIRST OF CONJUGATE COMPLEX PAIR
C                    -1, SECOND OF CONJUGATE COMPLEX PAIR ..........
      IP = 0
      N1 = N - 1
C
      DO 980 K = 1, N
         IF (WI(K) .EQ. 0.0E0 .OR. IP .LT. 0) GO TO 100
         IP = 1
         IF (SELECT(K) .AND. SELECT(K+1)) SELECT(K+1) = .FALSE.
  100    IF (.NOT. SELECT(K)) GO TO 960
         IF (WI(K) .NE. 0.0E0) S = S + 1
         IF (S .GT. MM) GO TO 1000
         IF (UK .GE. K) GO TO 200
C     .......... CHECK FOR POSSIBLE SPLITTING ..........
         DO 120 UK = K, N
            IF (UK .EQ. N) GO TO 140
            IF (A(UK+1,UK) .EQ. 0.0E0) GO TO 140
  120    CONTINUE
C     .......... COMPUTE INFINITY NORM OF LEADING UK BY UK
C                (HESSENBERG) MATRIX ..........
  140    NORM = 0.0E0
         MP = 1
C
         DO 180 I = 1, UK
            X = 0.0E0
C
            DO 160 J = MP, UK
  160       X = X + ABS(A(I,J))
C
            IF (X .GT. NORM) NORM = X
            MP = I
  180    CONTINUE
C     .......... EPS3 REPLACES ZERO PIVOT IN DECOMPOSITION
C                AND CLOSE ROOTS ARE MODIFIED BY EPS3 ..........
         IF (NORM .EQ. 0.0E0) NORM = 1.0E0
         EPS3 = EPSLON(NORM)
C     .......... GROWTO IS THE CRITERION FOR THE GROWTH ..........
         UKROOT = UK
         UKROOT = SQRT(UKROOT)
         GROWTO = 0.1E0 / UKROOT
  200    RLAMBD = WR(K)
         ILAMBD = WI(K)
         IF (K .EQ. 1) GO TO 280
         KM1 = K - 1
         GO TO 240
C     .......... PERTURB EIGENVALUE IF IT IS CLOSE
C                TO ANY PREVIOUS EIGENVALUE ..........
  220    RLAMBD = RLAMBD + EPS3
C     .......... FOR I=K-1 STEP -1 UNTIL 1 DO -- ..........
  240    DO 260 II = 1, KM1
            I = K - II
            IF (SELECT(I) .AND. ABS(WR(I)-RLAMBD) .LT. EPS3 .AND.
     X         ABS(WI(I)-ILAMBD) .LT. EPS3) GO TO 220
  260    CONTINUE
C
         WR(K) = RLAMBD
C     .......... PERTURB CONJUGATE EIGENVALUE TO MATCH ..........
         IP1 = K + IP
         WR(IP1) = RLAMBD
C     .......... FORM UPPER HESSENBERG A-RLAMBD*I (TRANSPOSED)
C                AND INITIAL REAL VECTOR ..........
  280    MP = 1
C
         DO 320 I = 1, UK
C
            DO 300 J = MP, UK
  300       RM1(J,I) = A(I,J)
C
            RM1(I,I) = RM1(I,I) - RLAMBD
            MP = I
            RV1(I) = EPS3
  320    CONTINUE
C
         ITS = 0
         IF (ILAMBD .NE. 0.0E0) GO TO 520
C     .......... REAL EIGENVALUE.
C                TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
C                REPLACING ZERO PIVOTS BY EPS3 ..........
         IF (UK .EQ. 1) GO TO 420
C
         DO 400 I = 2, UK
            MP = I - 1
            IF (ABS(RM1(MP,I)) .LE. ABS(RM1(MP,MP))) GO TO 360
C
            DO 340 J = MP, UK
               Y = RM1(J,I)
               RM1(J,I) = RM1(J,MP)
               RM1(J,MP) = Y
  340       CONTINUE
C
  360       IF (RM1(MP,MP) .EQ. 0.0E0) RM1(MP,MP) = EPS3
            X = RM1(MP,I) / RM1(MP,MP)
            IF (X .EQ. 0.0E0) GO TO 400
C
            DO 380 J = I, UK
  380       RM1(J,I) = RM1(J,I) - X * RM1(J,MP)
C
  400    CONTINUE
C
  420    IF (RM1(UK,UK) .EQ. 0.0E0) RM1(UK,UK) = EPS3
C     .......... BACK SUBSTITUTION FOR REAL VECTOR
C                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
  440    DO 500 II = 1, UK
            I = UK + 1 - II
            Y = RV1(I)
            IF (I .EQ. UK) GO TO 480
            IP1 = I + 1
C
            DO 460 J = IP1, UK
  460       Y = Y - RM1(J,I) * RV1(J)
C
  480       RV1(I) = Y / RM1(I,I)
  500    CONTINUE
C
         GO TO 740
C     .......... COMPLEX EIGENVALUE.
C                TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
C                REPLACING ZERO PIVOTS BY EPS3.  STORE IMAGINARY
C                PARTS IN UPPER TRIANGLE STARTING AT (1,3) ..........
  520    NS = N - S
         Z(1,S-1) = -ILAMBD
         Z(1,S) = 0.0E0
         IF (N .EQ. 2) GO TO 550
         RM1(1,3) = -ILAMBD
         Z(1,S-1) = 0.0E0
         IF (N .EQ. 3) GO TO 550
C
         DO 540 I = 4, N
  540    RM1(1,I) = 0.0E0
C
  550    DO 640 I = 2, UK
            MP = I - 1
            W = RM1(MP,I)
            IF (I .LT. N) T = RM1(MP,I+1)
            IF (I .EQ. N) T = Z(MP,S-1)
            X = RM1(MP,MP) * RM1(MP,MP) + T * T
            IF (W * W .LE. X) GO TO 580
            X = RM1(MP,MP) / W
            Y = T / W
            RM1(MP,MP) = W
            IF (I .LT. N) RM1(MP,I+1) = 0.0E0
            IF (I .EQ. N) Z(MP,S-1) = 0.0E0
C
            DO 560 J = I, UK
               W = RM1(J,I)
               RM1(J,I) = RM1(J,MP) - X * W
               RM1(J,MP) = W
               IF (J .LT. N1) GO TO 555
               L = J - NS
               Z(I,L) = Z(MP,L) - Y * W
               Z(MP,L) = 0.0E0
               GO TO 560
  555          RM1(I,J+2) = RM1(MP,J+2) - Y * W
               RM1(MP,J+2) = 0.0E0
  560       CONTINUE
C
            RM1(I,I) = RM1(I,I) - Y * ILAMBD
            IF (I .LT. N1) GO TO 570
            L = I - NS
            Z(MP,L) = -ILAMBD
            Z(I,L) = Z(I,L) + X * ILAMBD
            GO TO 640
  570       RM1(MP,I+2) = -ILAMBD
            RM1(I,I+2) = RM1(I,I+2) + X * ILAMBD
            GO TO 640
  580       IF (X .NE. 0.0E0) GO TO 600
            RM1(MP,MP) = EPS3
            IF (I .LT. N) RM1(MP,I+1) = 0.0E0
            IF (I .EQ. N) Z(MP,S-1) = 0.0E0
            T = 0.0E0
            X = EPS3 * EPS3
  600       W = W / X
            X = RM1(MP,MP) * W
            Y = -T * W
C
            DO 620 J = I, UK
               IF (J .LT. N1) GO TO 610
               L = J - NS
               T = Z(MP,L)
               Z(I,L) = -X * T - Y * RM1(J,MP)
               GO TO 615
  610          T = RM1(MP,J+2)
               RM1(I,J+2) = -X * T - Y * RM1(J,MP)
  615          RM1(J,I) = RM1(J,I) - X * RM1(J,MP) + Y * T
  620       CONTINUE
C
            IF (I .LT. N1) GO TO 630
            L = I - NS
            Z(I,L) = Z(I,L) - ILAMBD
            GO TO 640
  630       RM1(I,I+2) = RM1(I,I+2) - ILAMBD
  640    CONTINUE
C
         IF (UK .LT. N1) GO TO 650
         L = UK - NS
         T = Z(UK,L)
         GO TO 655
  650    T = RM1(UK,UK+2)
  655    IF (RM1(UK,UK) .EQ. 0.0E0 .AND. T .EQ. 0.0E0) RM1(UK,UK) = EPS3
C     .......... BACK SUBSTITUTION FOR COMPLEX VECTOR
C                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
  660    DO 720 II = 1, UK
            I = UK + 1 - II
            X = RV1(I)
            Y = 0.0E0
            IF (I .EQ. UK) GO TO 700
            IP1 = I + 1
C
            DO 680 J = IP1, UK
               IF (J .LT. N1) GO TO 670
               L = J - NS
               T = Z(I,L)
               GO TO 675
  670          T = RM1(I,J+2)
  675          X = X - RM1(J,I) * RV1(J) + T * RV2(J)
               Y = Y - RM1(J,I) * RV2(J) - T * RV1(J)
  680       CONTINUE
C
  700       IF (I .LT. N1) GO TO 710
            L = I - NS
            T = Z(I,L)
            GO TO 715
  710       T = RM1(I,I+2)
  715       CALL CDIV(X,Y,RM1(I,I),T,RV1(I),RV2(I))
  720    CONTINUE
C     .......... ACCEPTANCE TEST FOR REAL OR COMPLEX
C                EIGENVECTOR AND NORMALIZATION ..........
  740    ITS = ITS + 1
         NORM = 0.0E0
         NORMV = 0.0E0
C
         DO 780 I = 1, UK
            IF (ILAMBD .EQ. 0.0E0) X = ABS(RV1(I))
            IF (ILAMBD .NE. 0.0E0) X = PYTHAG(RV1(I),RV2(I))
            IF (NORMV .GE. X) GO TO 760
            NORMV = X
            J = I
  760       NORM = NORM + X
  780    CONTINUE
C
         IF (NORM .LT. GROWTO) GO TO 840
C     .......... ACCEPT VECTOR ..........
         X = RV1(J)
         IF (ILAMBD .EQ. 0.0E0) X = 1.0E0 / X
         IF (ILAMBD .NE. 0.0E0) Y = RV2(J)
C
         DO 820 I = 1, UK
            IF (ILAMBD .NE. 0.0E0) GO TO 800
            Z(I,S) = RV1(I) * X
            GO TO 820
  800       CALL CDIV(RV1(I),RV2(I),X,Y,Z(I,S-1),Z(I,S))
  820    CONTINUE
C
         IF (UK .EQ. N) GO TO 940
         J = UK + 1
         GO TO 900
C     .......... IN-LINE PROCEDURE FOR CHOOSING
C                A NEW STARTING VECTOR ..........
  840    IF (ITS .GE. UK) GO TO 880
         X = UKROOT
         Y = EPS3 / (X + 1.0E0)
         RV1(1) = EPS3
C
         DO 860 I = 2, UK
  860    RV1(I) = Y
C
         J = UK - ITS + 1
         RV1(J) = RV1(J) - EPS3 * X
         IF (ILAMBD .EQ. 0.0E0) GO TO 440
         GO TO 660
C     .......... SET ERROR -- UNACCEPTED EIGENVECTOR ..........
  880    J = 1
         IERR = -K
C     .......... SET REMAINING VECTOR COMPONENTS TO ZERO ..........
  900    DO 920 I = J, N
            Z(I,S) = 0.0E0
            IF (ILAMBD .NE. 0.0E0) Z(I,S-1) = 0.0E0
  920    CONTINUE
C
  940    S = S + 1
  960    IF (IP .EQ. (-1)) IP = 0
         IF (IP .EQ. 1) IP = -1
  980 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- UNDERESTIMATE OF EIGENVECTOR
C                SPACE REQUIRED ..........
 1000 IF (IERR .NE. 0) IERR = IERR - N
      IF (IERR .EQ. 0) IERR = -(2 * N + 1)
 1001 M = S - 1 - IABS(IP)
      RETURN
      END
      SUBROUTINE MINFIT(NM,M,N,A,W,IP,B,IERR,RV1)
C
      INTEGER I,J,K,L,M,N,II,IP,I1,KK,K1,LL,L1,M1,NM,ITS,IERR
      REAL A(NM,N),W(N),B(NM,IP),RV1(N)
      REAL C,F,G,H,S,X,Y,Z,TST1,TST2,SCALE,PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE MINFIT,
C     NUM. MATH. 14, 403-420(1970) BY GOLUB AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
C
C     THIS SUBROUTINE DETERMINES, TOWARDS THE SOLUTION OF THE LINEAR
C                                                        T
C     SYSTEM AX=B, THE SINGULAR VALUE DECOMPOSITION A=USV  OF A REAL
C                                         T
C     M BY N RECTANGULAR MATRIX, FORMING U B RATHER THAN U.  HOUSEHOLDER
C     BIDIAGONALIZATION AND A VARIANT OF THE QR ALGORITHM ARE USED.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.  NOTE THAT NM MUST BE AT LEAST
C          AS LARGE AS THE MAXIMUM OF M AND N.
C
C        M IS THE NUMBER OF ROWS OF A AND B.
C
C        N IS THE NUMBER OF COLUMNS OF A AND THE ORDER OF V.
C
C        A CONTAINS THE RECTANGULAR COEFFICIENT MATRIX OF THE SYSTEM.
C
C        IP IS THE NUMBER OF COLUMNS OF B.  IP CAN BE ZERO.
C
C        B CONTAINS THE CONSTANT COLUMN MATRIX OF THE SYSTEM
C          IF IP IS NOT ZERO.  OTHERWISE B IS NOT REFERENCED.
C
C     ON OUTPUT
C
C        A HAS BEEN OVERWRITTEN BY THE MATRIX V (ORTHOGONAL) OF THE
C          DECOMPOSITION IN ITS FIRST N ROWS AND COLUMNS.  IF AN
C          ERROR EXIT IS MADE, THE COLUMNS OF V CORRESPONDING TO
C          INDICES OF CORRECT SINGULAR VALUES SHOULD BE CORRECT.
C
C        W CONTAINS THE N (NON-NEGATIVE) SINGULAR VALUES OF A (THE
C          DIAGONAL ELEMENTS OF S).  THEY ARE UNORDERED.  IF AN
C          ERROR EXIT IS MADE, THE SINGULAR VALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,IERR+2,...,N.
C
C                                   T
C        B HAS BEEN OVERWRITTEN BY U B.  IF AN ERROR EXIT IS MADE,
C                       T
C          THE ROWS OF U B CORRESPONDING TO INDICES OF CORRECT
C          SINGULAR VALUES SHOULD BE CORRECT.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          K          IF THE K-TH SINGULAR VALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C        RV1 IS A TEMPORARY STORAGE ARRAY.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
C     .......... HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM ..........
      G = 0.0E0
      SCALE = 0.0E0
      X = 0.0E0
C
      DO 300 I = 1, N
         L = I + 1
         RV1(I) = SCALE * G
         G = 0.0E0
         S = 0.0E0
         SCALE = 0.0E0
         IF (I .GT. M) GO TO 210
C
         DO 120 K = I, M
  120    SCALE = SCALE + ABS(A(K,I))
C
         IF (SCALE .EQ. 0.0E0) GO TO 210
C
         DO 130 K = I, M
            A(K,I) = A(K,I) / SCALE
            S = S + A(K,I)**2
  130    CONTINUE
C
         F = A(I,I)
         G = -SIGN(SQRT(S),F)
         H = F * G - S
         A(I,I) = F - G
         IF (I .EQ. N) GO TO 160
C
         DO 150 J = L, N
            S = 0.0E0
C
            DO 140 K = I, M
  140       S = S + A(K,I) * A(K,J)
C
            F = S / H
C
            DO 150 K = I, M
               A(K,J) = A(K,J) + F * A(K,I)
  150    CONTINUE
C
  160    IF (IP .EQ. 0) GO TO 190
C
         DO 180 J = 1, IP
            S = 0.0E0
C
            DO 170 K = I, M
  170       S = S + A(K,I) * B(K,J)
C
            F = S / H
C
            DO 180 K = I, M
               B(K,J) = B(K,J) + F * A(K,I)
  180    CONTINUE
C
  190    DO 200 K = I, M
  200    A(K,I) = SCALE * A(K,I)
C
  210    W(I) = SCALE * G
         G = 0.0E0
         S = 0.0E0
         SCALE = 0.0E0
         IF (I .GT. M .OR. I .EQ. N) GO TO 290
C
         DO 220 K = L, N
  220    SCALE = SCALE + ABS(A(I,K))
C
         IF (SCALE .EQ. 0.0E0) GO TO 290
C
         DO 230 K = L, N
            A(I,K) = A(I,K) / SCALE
            S = S + A(I,K)**2
  230    CONTINUE
C
         F = A(I,L)
         G = -SIGN(SQRT(S),F)
         H = F * G - S
         A(I,L) = F - G
C
         DO 240 K = L, N
  240    RV1(K) = A(I,K) / H
C
         IF (I .EQ. M) GO TO 270
C
         DO 260 J = L, M
            S = 0.0E0
C
            DO 250 K = L, N
  250       S = S + A(J,K) * A(I,K)
C
            DO 260 K = L, N
               A(J,K) = A(J,K) + S * RV1(K)
  260    CONTINUE
C
  270    DO 280 K = L, N
  280    A(I,K) = SCALE * A(I,K)
C
  290    X = AMAX1(X,ABS(W(I))+ABS(RV1(I)))
  300 CONTINUE
C     .......... ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS.
C                FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 400 II = 1, N
         I = N + 1 - II
         IF (I .EQ. N) GO TO 390
         IF (G .EQ. 0.0E0) GO TO 360
C
         DO 320 J = L, N
C     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
  320    A(J,I) = (A(I,J) / A(I,L)) / G
C
         DO 350 J = L, N
            S = 0.0E0
C
            DO 340 K = L, N
  340       S = S + A(I,K) * A(K,J)
C
            DO 350 K = L, N
               A(K,J) = A(K,J) + S * A(K,I)
  350    CONTINUE
C
  360    DO 380 J = L, N
            A(I,J) = 0.0E0
            A(J,I) = 0.0E0
  380    CONTINUE
C
  390    A(I,I) = 1.0E0
         G = RV1(I)
         L = I
  400 CONTINUE
C
      IF (M .GE. N .OR. IP .EQ. 0) GO TO 510
      M1 = M + 1
C
      DO 500 I = M1, N
C
         DO 500 J = 1, IP
            B(I,J) = 0.0E0
  500 CONTINUE
C     .......... DIAGONALIZATION OF THE BIDIAGONAL FORM ..........
  510 TST1 = X
C     .......... FOR K=N STEP -1 UNTIL 1 DO -- ..........
      DO 700 KK = 1, N
         K1 = N - KK
         K = K1 + 1
         ITS = 0
C     .......... TEST FOR SPLITTING.
C                FOR L=K STEP -1 UNTIL 1 DO -- ..........
  520    DO 530 LL = 1, K
            L1 = K - LL
            L = L1 + 1
            TST2 = TST1 + ABS(RV1(L))
            IF (TST2 .EQ. TST1) GO TO 565
C     .......... RV1(1) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
            TST2 = TST1 + ABS(W(L1))
            IF (TST2 .EQ. TST1) GO TO 540
  530    CONTINUE
C     .......... CANCELLATION OF RV1(L) IF L GREATER THAN 1 ..........
  540    C = 0.0E0
         S = 1.0E0
C
         DO 560 I = L, K
            F = S * RV1(I)
            RV1(I) = C * RV1(I)
            TST2 = TST1 + ABS(F)
            IF (TST2 .EQ. TST1) GO TO 565
            G = W(I)
            H = PYTHAG(F,G)
            W(I) = H
            C = G / H
            S = -F / H
            IF (IP .EQ. 0) GO TO 560
C
            DO 550 J = 1, IP
               Y = B(L1,J)
               Z = B(I,J)
               B(L1,J) = Y * C + Z * S
               B(I,J) = -Y * S + Z * C
  550       CONTINUE
C
  560    CONTINUE
C     .......... TEST FOR CONVERGENCE ..........
  565    Z = W(K)
         IF (L .EQ. K) GO TO 650
C     .......... SHIFT FROM BOTTOM 2 BY 2 MINOR ..........
         IF (ITS .EQ. 30) GO TO 1000
         ITS = ITS + 1
         X = W(L)
         Y = W(K1)
         G = RV1(K1)
         H = RV1(K)
         F = 0.5E0 * (((G + Z) / H) * ((G - Z) / Y) + Y / H - H / Y)
         G = PYTHAG(F,1.0E0)
         F = X - (Z / X) * Z + (H / X) * (Y / (F + SIGN(G,F)) - H)
C     .......... NEXT QR TRANSFORMATION ..........
         C = 1.0E0
         S = 1.0E0
C
         DO 600 I1 = L, K1
            I = I1 + 1
            G = RV1(I)
            Y = W(I)
            H = S * G
            G = C * G
            Z = PYTHAG(F,H)
            RV1(I1) = Z
            C = F / Z
            S = H / Z
            F = X * C + G * S
            G = -X * S + G * C
            H = Y * S
            Y = Y * C
C
            DO 570 J = 1, N
               X = A(J,I1)
               Z = A(J,I)
               A(J,I1) = X * C + Z * S
               A(J,I) = -X * S + Z * C
  570       CONTINUE
C
            Z = PYTHAG(F,H)
            W(I1) = Z
C     .......... ROTATION CAN BE ARBITRARY IF Z IS ZERO ..........
            IF (Z .EQ. 0.0E0) GO TO 580
            C = F / Z
            S = H / Z
  580       F = C * G + S * Y
            X = -S * G + C * Y
            IF (IP .EQ. 0) GO TO 600
C
            DO 590 J = 1, IP
               Y = B(I1,J)
               Z = B(I,J)
               B(I1,J) = Y * C + Z * S
               B(I,J) = -Y * S + Z * C
  590       CONTINUE
C
  600    CONTINUE
C
         RV1(L) = 0.0E0
         RV1(K) = F
         W(K) = X
         GO TO 520
C     .......... CONVERGENCE ..........
  650    IF (Z .GE. 0.0E0) GO TO 700
C     .......... W(K) IS MADE NON-NEGATIVE ..........
         W(K) = -Z
C
         DO 690 J = 1, N
  690    A(J,K) = -A(J,K)
C
  700 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO A
C                SINGULAR VALUE AFTER 30 ITERATIONS ..........
 1000 IERR = K
 1001 RETURN
      END
      SUBROUTINE ORTBAK(NM,LOW,IGH,A,ORT,M,Z)
C
      INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
      REAL A(NM,IGH),ORT(IGH),Z(NM,M)
      REAL G
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ORTBAK,
C     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL GENERAL
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     UPPER HESSENBERG MATRIX DETERMINED BY  ORTHES.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1 AND IGH EQUAL TO THE ORDER OF THE MATRIX.
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  ORTHES
C          IN ITS STRICT LOWER TRIANGLE.
C
C        ORT CONTAINS FURTHER INFORMATION ABOUT THE TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  ORTHES.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C        M IS THE NUMBER OF COLUMNS OF Z TO BE BACK TRANSFORMED.
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGEN-
C          VECTORS TO BE BACK TRANSFORMED IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE
C          TRANSFORMED EIGENVECTORS IN ITS FIRST M COLUMNS.
C
C        ORT HAS BEEN ALTERED.
C
C     NOTE THAT ORTBAK PRESERVES VECTOR EUCLIDEAN NORMS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 140 MM = KP1, LA
         MP = LOW + IGH - MM
         IF (A(MP,MP-1) .EQ. 0.0E0) GO TO 140
         MP1 = MP + 1
C
         DO 100 I = MP1, IGH
  100    ORT(I) = A(I,MP-1)
C
         DO 130 J = 1, M
            G = 0.0E0
C
            DO 110 I = MP, IGH
  110       G = G + ORT(I) * Z(I,J)
C     .......... DIVISOR BELOW IS NEGATIVE OF H FORMED IN ORTHES.
C                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            G = (G / ORT(MP)) / A(MP,MP-1)
C
            DO 120 I = MP, IGH
  120       Z(I,J) = Z(I,J) + G * ORT(I)
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE ORTHES(NM,N,LOW,IGH,A,ORT)
C
      INTEGER I,J,M,N,II,JJ,LA,MP,NM,IGH,KP1,LOW
      REAL A(NM,N),ORT(IGH)
      REAL F,G,H,SCALE
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ORTHES,
C     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     GIVEN A REAL GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        A CONTAINS THE INPUT MATRIX.
C
C     ON OUTPUT
C
C        A CONTAINS THE HESSENBERG MATRIX.  INFORMATION ABOUT
C          THE ORTHOGONAL TRANSFORMATIONS USED IN THE REDUCTION
C          IS STORED IN THE REMAINING TRIANGLE UNDER THE
C          HESSENBERG MATRIX.
C
C        ORT CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C
      DO 180 M = KP1, LA
         H = 0.0E0
         ORT(M) = 0.0E0
         SCALE = 0.0E0
C     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ..........
         DO 90 I = M, IGH
   90    SCALE = SCALE + ABS(A(I,M-1))
C
         IF (SCALE .EQ. 0.0E0) GO TO 180
         MP = M + IGH
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
         DO 100 II = M, IGH
            I = MP - II
            ORT(I) = A(I,M-1) / SCALE
            H = H + ORT(I) * ORT(I)
  100    CONTINUE
C
         G = -SIGN(SQRT(H),ORT(M))
         H = H - ORT(M) * G
         ORT(M) = ORT(M) - G
C     .......... FORM (I-(U*UT)/H) * A ..........
         DO 130 J = M, N
            F = 0.0E0
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
            DO 110 II = M, IGH
               I = MP - II
               F = F + ORT(I) * A(I,J)
  110       CONTINUE
C
            F = F / H
C
            DO 120 I = M, IGH
  120       A(I,J) = A(I,J) - F * ORT(I)
C
  130    CONTINUE
C     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ..........
         DO 160 I = 1, IGH
            F = 0.0E0
C     .......... FOR J=IGH STEP -1 UNTIL M DO -- ..........
            DO 140 JJ = M, IGH
               J = MP - JJ
               F = F + ORT(J) * A(I,J)
  140       CONTINUE
C
            F = F / H
C
            DO 150 J = M, IGH
  150       A(I,J) = A(I,J) - F * ORT(J)
C
  160    CONTINUE
C
         ORT(M) = SCALE * ORT(M)
         A(M,M-1) = SCALE * G
  180 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE ORTRAN(NM,N,LOW,IGH,A,ORT,Z)
C
      INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
      REAL A(NM,IGH),ORT(IGH),Z(NM,N)
      REAL G
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ORTRANS,
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     THIS SUBROUTINE ACCUMULATES THE ORTHOGONAL SIMILARITY
C     TRANSFORMATIONS USED IN THE REDUCTION OF A REAL GENERAL
C     MATRIX TO UPPER HESSENBERG FORM BY  ORTHES.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  ORTHES
C          IN ITS STRICT LOWER TRIANGLE.
C
C        ORT CONTAINS FURTHER INFORMATION ABOUT THE TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  ORTHES.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     ON OUTPUT
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  ORTHES.
C
C        ORT HAS BEEN ALTERED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
C     .......... INITIALIZE Z TO IDENTITY MATRIX ..........
      DO 80 J = 1, N
C
         DO 60 I = 1, N
   60    Z(I,J) = 0.0E0
C
         Z(J,J) = 1.0E0
   80 CONTINUE
C
      KL = IGH - LOW - 1
      IF (KL .LT. 1) GO TO 200
C     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 140 MM = 1, KL
         MP = IGH - MM
         IF (A(MP,MP-1) .EQ. 0.0E0) GO TO 140
         MP1 = MP + 1
C
         DO 100 I = MP1, IGH
  100    ORT(I) = A(I,MP-1)
C
         DO 130 J = MP, IGH
            G = 0.0E0
C
            DO 110 I = MP, IGH
  110       G = G + ORT(I) * Z(I,J)
C     .......... DIVISOR BELOW IS NEGATIVE OF H FORMED IN ORTHES.
C                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            G = (G / ORT(MP)) / A(MP,MP-1)
C
            DO 120 I = MP, IGH
  120       Z(I,J) = Z(I,J) + G * ORT(I)
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE QZHES(NM,N,A,B,MATZ,Z)
C
      INTEGER I,J,K,L,N,LB,L1,NM,NK1,NM1,NM2
      REAL A(NM,N),B(NM,N),Z(NM,N)
      REAL R,S,T,U1,U2,V1,V2,RHO
      LOGICAL MATZ
C
C     THIS SUBROUTINE IS THE FIRST STEP OF THE QZ ALGORITHM
C     FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF REAL GENERAL MATRICES AND
C     REDUCES ONE OF THEM TO UPPER HESSENBERG FORM AND THE OTHER
C     TO UPPER TRIANGULAR FORM USING ORTHOGONAL TRANSFORMATIONS.
C     IT IS USUALLY FOLLOWED BY  QZIT,  QZVAL  AND, POSSIBLY,  QZVEC.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRICES.
C
C        A CONTAINS A REAL GENERAL MATRIX.
C
C        B CONTAINS A REAL GENERAL MATRIX.
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C     ON OUTPUT
C
C        A HAS BEEN REDUCED TO UPPER HESSENBERG FORM.  THE ELEMENTS
C          BELOW THE FIRST SUBDIAGONAL HAVE BEEN SET TO ZERO.
C
C        B HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO.
C
C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS IF
C          MATZ HAS BEEN SET TO .TRUE.  OTHERWISE, Z IS NOT REFERENCED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
C     .......... INITIALIZE Z ..........
      IF (.NOT. MATZ) GO TO 10
C
      DO 3 J = 1, N
C
         DO 2 I = 1, N
            Z(I,J) = 0.0E0
    2    CONTINUE
C
         Z(J,J) = 1.0E0
    3 CONTINUE
C     .......... REDUCE B TO UPPER TRIANGULAR FORM ..........
   10 IF (N .LE. 1) GO TO 170
      NM1 = N - 1
C
      DO 100 L = 1, NM1
         L1 = L + 1
         S = 0.0E0
C
         DO 20 I = L1, N
            S = S + ABS(B(I,L))
   20    CONTINUE
C
         IF (S .EQ. 0.0E0) GO TO 100
         S = S + ABS(B(L,L))
         R = 0.0E0
C
         DO 25 I = L, N
            B(I,L) = B(I,L) / S
            R = R + B(I,L)**2
   25    CONTINUE
C
         R = SIGN(SQRT(R),B(L,L))
         B(L,L) = B(L,L) + R
         RHO = R * B(L,L)
C
         DO 50 J = L1, N
            T = 0.0E0
C
            DO 30 I = L, N
               T = T + B(I,L) * B(I,J)
   30       CONTINUE
C
            T = -T / RHO
C
            DO 40 I = L, N
               B(I,J) = B(I,J) + T * B(I,L)
   40       CONTINUE
C
   50    CONTINUE
C
         DO 80 J = 1, N
            T = 0.0E0
C
            DO 60 I = L, N
               T = T + B(I,L) * A(I,J)
   60       CONTINUE
C
            T = -T / RHO
C
            DO 70 I = L, N
               A(I,J) = A(I,J) + T * B(I,L)
   70       CONTINUE
C
   80    CONTINUE
C
         B(L,L) = -S * R
C
         DO 90 I = L1, N
            B(I,L) = 0.0E0
   90    CONTINUE
C
  100 CONTINUE
C     .......... REDUCE A TO UPPER HESSENBERG FORM, WHILE
C                KEEPING B TRIANGULAR ..........
      IF (N .EQ. 2) GO TO 170
      NM2 = N - 2
C
      DO 160 K = 1, NM2
         NK1 = NM1 - K
C     .......... FOR L=N-1 STEP -1 UNTIL K+1 DO -- ..........
         DO 150 LB = 1, NK1
            L = N - LB
            L1 = L + 1
C     .......... ZERO A(L+1,K) ..........
            S = ABS(A(L,K)) + ABS(A(L1,K))
            IF (S .EQ. 0.0E0) GO TO 150
            U1 = A(L,K) / S
            U2 = A(L1,K) / S
            R = SIGN(SQRT(U1*U1+U2*U2),U1)
            V1 =  -(U1 + R) / R
            V2 = -U2 / R
            U2 = V2 / V1
C
            DO 110 J = K, N
               T = A(L,J) + U2 * A(L1,J)
               A(L,J) = A(L,J) + T * V1
               A(L1,J) = A(L1,J) + T * V2
  110       CONTINUE
C
            A(L1,K) = 0.0E0
C
            DO 120 J = L, N
               T = B(L,J) + U2 * B(L1,J)
               B(L,J) = B(L,J) + T * V1
               B(L1,J) = B(L1,J) + T * V2
  120       CONTINUE
C     .......... ZERO B(L+1,L) ..........
            S = ABS(B(L1,L1)) + ABS(B(L1,L))
            IF (S .EQ. 0.0E0) GO TO 150
            U1 = B(L1,L1) / S
            U2 = B(L1,L) / S
            R = SIGN(SQRT(U1*U1+U2*U2),U1)
            V1 =  -(U1 + R) / R
            V2 = -U2 / R
            U2 = V2 / V1
C
            DO 130 I = 1, L1
               T = B(I,L1) + U2 * B(I,L)
               B(I,L1) = B(I,L1) + T * V1
               B(I,L) = B(I,L) + T * V2
  130       CONTINUE
C
            B(L1,L) = 0.0E0
C
            DO 140 I = 1, N
               T = A(I,L1) + U2 * A(I,L)
               A(I,L1) = A(I,L1) + T * V1
               A(I,L) = A(I,L) + T * V2
  140       CONTINUE
C
            IF (.NOT. MATZ) GO TO 150
C
            DO 145 I = 1, N
               T = Z(I,L1) + U2 * Z(I,L)
               Z(I,L1) = Z(I,L1) + T * V1
               Z(I,L) = Z(I,L) + T * V2
  145       CONTINUE
C
  150    CONTINUE
C
  160 CONTINUE
C
  170 RETURN
      END
      SUBROUTINE QZIT(NM,N,A,B,EPS1,MATZ,Z,IERR)
C
      INTEGER I,J,K,L,N,EN,K1,K2,LD,LL,L1,NA,NM,ISH,ITN,ITS,KM1,LM1,
     X        ENM2,IERR,LOR1,ENORN
      REAL A(NM,N),B(NM,N),Z(NM,N)
      REAL R,S,T,A1,A2,A3,EP,SH,U1,U2,U3,V1,V2,V3,ANI,A11,
     X       A12,A21,A22,A33,A34,A43,A44,BNI,B11,B12,B22,B33,B34,
     X       B44,EPSA,EPSB,EPS1,ANORM,BNORM,EPSLON
      LOGICAL MATZ,NOTLAS
C
C     THIS SUBROUTINE IS THE SECOND STEP OF THE QZ ALGORITHM
C     FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART,
C     AS MODIFIED IN TECHNICAL NOTE NASA TN D-7305(1973) BY WARD.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF REAL MATRICES, ONE OF THEM
C     IN UPPER HESSENBERG FORM AND THE OTHER IN UPPER TRIANGULAR FORM.
C     IT REDUCES THE HESSENBERG MATRIX TO QUASI-TRIANGULAR FORM USING
C     ORTHOGONAL TRANSFORMATIONS WHILE MAINTAINING THE TRIANGULAR FORM
C     OF THE OTHER MATRIX.  IT IS USUALLY PRECEDED BY  QZHES  AND
C     FOLLOWED BY  QZVAL  AND, POSSIBLY,  QZVEC.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRICES.
C
C        A CONTAINS A REAL UPPER HESSENBERG MATRIX.
C
C        B CONTAINS A REAL UPPER TRIANGULAR MATRIX.
C
C        EPS1 IS A TOLERANCE USED TO DETERMINE NEGLIGIBLE ELEMENTS.
C          EPS1 = 0.0 (OR NEGATIVE) MAY BE INPUT, IN WHICH CASE AN
C          ELEMENT WILL BE NEGLECTED ONLY IF IT IS LESS THAN ROUNDOFF
C          ERROR TIMES THE NORM OF ITS MATRIX.  IF THE INPUT EPS1 IS
C          POSITIVE, THEN AN ELEMENT WILL BE CONSIDERED NEGLIGIBLE
C          IF IT IS LESS THAN EPS1 TIMES THE NORM OF ITS MATRIX.  A
C          POSITIVE VALUE OF EPS1 MAY RESULT IN FASTER EXECUTION,
C          BUT LESS ACCURATE RESULTS.
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C        Z CONTAINS, IF MATZ HAS BEEN SET TO .TRUE., THE
C          TRANSFORMATION MATRIX PRODUCED IN THE REDUCTION
C          BY  QZHES, IF PERFORMED, OR ELSE THE IDENTITY MATRIX.
C          IF MATZ HAS BEEN SET TO .FALSE., Z IS NOT REFERENCED.
C
C     ON OUTPUT
C
C        A HAS BEEN REDUCED TO QUASI-TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE FIRST SUBDIAGONAL ARE STILL ZERO AND NO TWO
C          CONSECUTIVE SUBDIAGONAL ELEMENTS ARE NONZERO.
C
C        B IS STILL IN UPPER TRIANGULAR FORM, ALTHOUGH ITS ELEMENTS
C          HAVE BEEN ALTERED.  THE LOCATION B(N,1) IS USED TO STORE
C          EPS1 TIMES THE NORM OF B FOR LATER USE BY  QZVAL  AND  QZVEC.
C
C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS
C          (FOR BOTH STEPS) IF MATZ HAS BEEN SET TO .TRUE..
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
C     .......... COMPUTE EPSA,EPSB ..........
      ANORM = 0.0E0
      BNORM = 0.0E0
C
      DO 30 I = 1, N
         ANI = 0.0E0
         IF (I .NE. 1) ANI = ABS(A(I,I-1))
         BNI = 0.0E0
C
         DO 20 J = I, N
            ANI = ANI + ABS(A(I,J))
            BNI = BNI + ABS(B(I,J))
   20    CONTINUE
C
         IF (ANI .GT. ANORM) ANORM = ANI
         IF (BNI .GT. BNORM) BNORM = BNI
   30 CONTINUE
C
      IF (ANORM .EQ. 0.0E0) ANORM = 1.0E0
      IF (BNORM .EQ. 0.0E0) BNORM = 1.0E0
      EP = EPS1
      IF (EP .GT. 0.0E0) GO TO 50
C     .......... USE ROUNDOFF LEVEL IF EPS1 IS ZERO ..........
      EP = EPSLON(1.0E0)
   50 EPSA = EP * ANORM
      EPSB = EP * BNORM
C     .......... REDUCE A TO QUASI-TRIANGULAR FORM, WHILE
C                KEEPING B TRIANGULAR ..........
      LOR1 = 1
      ENORN = N
      EN = N
      ITN = 30*N
C     .......... BEGIN QZ STEP ..........
   60 IF (EN .LE. 2) GO TO 1001
      IF (.NOT. MATZ) ENORN = EN
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
   70 ISH = 2
C     .......... CHECK FOR CONVERGENCE OR REDUCIBILITY.
C                FOR L=EN STEP -1 UNTIL 1 DO -- ..........
      DO 80 LL = 1, EN
         LM1 = EN - LL
         L = LM1 + 1
         IF (L .EQ. 1) GO TO 95
         IF (ABS(A(L,LM1)) .LE. EPSA) GO TO 90
   80 CONTINUE
C
   90 A(L,LM1) = 0.0E0
      IF (L .LT. NA) GO TO 95
C     .......... 1-BY-1 OR 2-BY-2 BLOCK ISOLATED ..........
      EN = LM1
      GO TO 60
C     .......... CHECK FOR SMALL TOP OF B ..........
   95 LD = L
  100 L1 = L + 1
      B11 = B(L,L)
      IF (ABS(B11) .GT. EPSB) GO TO 120
      B(L,L) = 0.0E0
      S = ABS(A(L,L)) + ABS(A(L1,L))
      U1 = A(L,L) / S
      U2 = A(L1,L) / S
      R = SIGN(SQRT(U1*U1+U2*U2),U1)
      V1 = -(U1 + R) / R
      V2 = -U2 / R
      U2 = V2 / V1
C
      DO 110 J = L, ENORN
         T = A(L,J) + U2 * A(L1,J)
         A(L,J) = A(L,J) + T * V1
         A(L1,J) = A(L1,J) + T * V2
         T = B(L,J) + U2 * B(L1,J)
         B(L,J) = B(L,J) + T * V1
         B(L1,J) = B(L1,J) + T * V2
  110 CONTINUE
C
      IF (L .NE. 1) A(L,LM1) = -A(L,LM1)
      LM1 = L
      L = L1
      GO TO 90
  120 A11 = A(L,L) / B11
      A21 = A(L1,L) / B11
      IF (ISH .EQ. 1) GO TO 140
C     .......... ITERATION STRATEGY ..........
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .EQ. 10) GO TO 155
C     .......... DETERMINE TYPE OF SHIFT ..........
      B22 = B(L1,L1)
      IF (ABS(B22) .LT. EPSB) B22 = EPSB
      B33 = B(NA,NA)
      IF (ABS(B33) .LT. EPSB) B33 = EPSB
      B44 = B(EN,EN)
      IF (ABS(B44) .LT. EPSB) B44 = EPSB
      A33 = A(NA,NA) / B33
      A34 = A(NA,EN) / B44
      A43 = A(EN,NA) / B33
      A44 = A(EN,EN) / B44
      B34 = B(NA,EN) / B44
      T = 0.5E0 * (A43 * B34 - A33 - A44)
      R = T * T + A34 * A43 - A33 * A44
      IF (R .LT. 0.0E0) GO TO 150
C     .......... DETERMINE SINGLE SHIFT ZEROTH COLUMN OF A ..........
      ISH = 1
      R = SQRT(R)
      SH = -T + R
      S = -T - R
      IF (ABS(S-A44) .LT. ABS(SH-A44)) SH = S
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS OF A.
C                FOR L=EN-2 STEP -1 UNTIL LD DO -- ..........
      DO 130 LL = LD, ENM2
         L = ENM2 + LD - LL
         IF (L .EQ. LD) GO TO 140
         LM1 = L - 1
         L1 = L + 1
         T = A(L,L)
         IF (ABS(B(L,L)) .GT. EPSB) T = T - SH * B(L,L)
         IF (ABS(A(L,LM1)) .LE. ABS(T/A(L1,L)) * EPSA) GO TO 100
  130 CONTINUE
C
  140 A1 = A11 - SH
      A2 = A21
      IF (L .NE. LD) A(L,LM1) = -A(L,LM1)
      GO TO 160
C     .......... DETERMINE DOUBLE SHIFT ZEROTH COLUMN OF A ..........
  150 A12 = A(L,L1) / B22
      A22 = A(L1,L1) / B22
      B12 = B(L,L1) / B22
      A1 = ((A33 - A11) * (A44 - A11) - A34 * A43 + A43 * B34 * A11)
     X     / A21 + A12 - A11 * B12
      A2 = (A22 - A11) - A21 * B12 - (A33 - A11) - (A44 - A11)
     X     + A43 * B34
      A3 = A(L1+1,L1) / B22
      GO TO 160
C     .......... AD HOC SHIFT ..........
  155 A1 = 0.0E0
      A2 = 1.0E0
      A3 = 1.1605E0
  160 ITS = ITS + 1
      ITN = ITN - 1
      IF (.NOT. MATZ) LOR1 = LD
C     .......... MAIN LOOP ..........
      DO 260 K = L, NA
         NOTLAS = K .NE. NA .AND. ISH .EQ. 2
         K1 = K + 1
         K2 = K + 2
         KM1 = MAX0(K-1,L)
         LL = MIN0(EN,K1+ISH)
         IF (NOTLAS) GO TO 190
C     .......... ZERO A(K+1,K-1) ..........
         IF (K .EQ. L) GO TO 170
         A1 = A(K,KM1)
         A2 = A(K1,KM1)
  170    S = ABS(A1) + ABS(A2)
         IF (S .EQ. 0.0E0) GO TO 70
         U1 = A1 / S
         U2 = A2 / S
         R = SIGN(SQRT(U1*U1+U2*U2),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         U2 = V2 / V1
C
         DO 180 J = KM1, ENORN
            T = A(K,J) + U2 * A(K1,J)
            A(K,J) = A(K,J) + T * V1
            A(K1,J) = A(K1,J) + T * V2
            T = B(K,J) + U2 * B(K1,J)
            B(K,J) = B(K,J) + T * V1
            B(K1,J) = B(K1,J) + T * V2
  180    CONTINUE
C
         IF (K .NE. L) A(K1,KM1) = 0.0E0
         GO TO 240
C     .......... ZERO A(K+1,K-1) AND A(K+2,K-1) ..........
  190    IF (K .EQ. L) GO TO 200
         A1 = A(K,KM1)
         A2 = A(K1,KM1)
         A3 = A(K2,KM1)
  200    S = ABS(A1) + ABS(A2) + ABS(A3)
         IF (S .EQ. 0.0E0) GO TO 260
         U1 = A1 / S
         U2 = A2 / S
         U3 = A3 / S
         R = SIGN(SQRT(U1*U1+U2*U2+U3*U3),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         V3 = -U3 / R
         U2 = V2 / V1
         U3 = V3 / V1
C
         DO 210 J = KM1, ENORN
            T = A(K,J) + U2 * A(K1,J) + U3 * A(K2,J)
            A(K,J) = A(K,J) + T * V1
            A(K1,J) = A(K1,J) + T * V2
            A(K2,J) = A(K2,J) + T * V3
            T = B(K,J) + U2 * B(K1,J) + U3 * B(K2,J)
            B(K,J) = B(K,J) + T * V1
            B(K1,J) = B(K1,J) + T * V2
            B(K2,J) = B(K2,J) + T * V3
  210    CONTINUE
C
         IF (K .EQ. L) GO TO 220
         A(K1,KM1) = 0.0E0
         A(K2,KM1) = 0.0E0
C     .......... ZERO B(K+2,K+1) AND B(K+2,K) ..........
  220    S = ABS(B(K2,K2)) + ABS(B(K2,K1)) + ABS(B(K2,K))
         IF (S .EQ. 0.0E0) GO TO 240
         U1 = B(K2,K2) / S
         U2 = B(K2,K1) / S
         U3 = B(K2,K) / S
         R = SIGN(SQRT(U1*U1+U2*U2+U3*U3),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         V3 = -U3 / R
         U2 = V2 / V1
         U3 = V3 / V1
C
         DO 230 I = LOR1, LL
            T = A(I,K2) + U2 * A(I,K1) + U3 * A(I,K)
            A(I,K2) = A(I,K2) + T * V1
            A(I,K1) = A(I,K1) + T * V2
            A(I,K) = A(I,K) + T * V3
            T = B(I,K2) + U2 * B(I,K1) + U3 * B(I,K)
            B(I,K2) = B(I,K2) + T * V1
            B(I,K1) = B(I,K1) + T * V2
            B(I,K) = B(I,K) + T * V3
  230    CONTINUE
C
         B(K2,K) = 0.0E0
         B(K2,K1) = 0.0E0
         IF (.NOT. MATZ) GO TO 240
C
         DO 235 I = 1, N
            T = Z(I,K2) + U2 * Z(I,K1) + U3 * Z(I,K)
            Z(I,K2) = Z(I,K2) + T * V1
            Z(I,K1) = Z(I,K1) + T * V2
            Z(I,K) = Z(I,K) + T * V3
  235    CONTINUE
C     .......... ZERO B(K+1,K) ..........
  240    S = ABS(B(K1,K1)) + ABS(B(K1,K))
         IF (S .EQ. 0.0E0) GO TO 260
         U1 = B(K1,K1) / S
         U2 = B(K1,K) / S
         R = SIGN(SQRT(U1*U1+U2*U2),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         U2 = V2 / V1
C
         DO 250 I = LOR1, LL
            T = A(I,K1) + U2 * A(I,K)
            A(I,K1) = A(I,K1) + T * V1
            A(I,K) = A(I,K) + T * V2
            T = B(I,K1) + U2 * B(I,K)
            B(I,K1) = B(I,K1) + T * V1
            B(I,K) = B(I,K) + T * V2
  250    CONTINUE
C
         B(K1,K) = 0.0E0
         IF (.NOT. MATZ) GO TO 260
C
         DO 255 I = 1, N
            T = Z(I,K1) + U2 * Z(I,K)
            Z(I,K1) = Z(I,K1) + T * V1
            Z(I,K) = Z(I,K) + T * V2
  255    CONTINUE
C
  260 CONTINUE
C     .......... END QZ STEP ..........
      GO TO 70
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
C     .......... SAVE EPSB FOR USE BY QZVAL AND QZVEC ..........
 1001 IF (N .GT. 1) B(N,1) = EPSB
      RETURN
      END
      SUBROUTINE QZVAL(NM,N,A,B,ALFR,ALFI,BETA,MATZ,Z)
C
      INTEGER I,J,N,EN,NA,NM,NN,ISW
      REAL A(NM,N),B(NM,N),ALFR(N),ALFI(N),BETA(N),Z(NM,N)
      REAL C,D,E,R,S,T,AN,A1,A2,BN,CQ,CZ,DI,DR,EI,TI,TR,U1,
     X       U2,V1,V2,A1I,A11,A12,A2I,A21,A22,B11,B12,B22,SQI,SQR,
     X       SSI,SSR,SZI,SZR,A11I,A11R,A12I,A12R,A22I,A22R,EPSB
      LOGICAL MATZ
C
C     THIS SUBROUTINE IS THE THIRD STEP OF THE QZ ALGORITHM
C     FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF REAL MATRICES, ONE OF THEM
C     IN QUASI-TRIANGULAR FORM AND THE OTHER IN UPPER TRIANGULAR FORM.
C     IT REDUCES THE QUASI-TRIANGULAR MATRIX FURTHER, SO THAT ANY
C     REMAINING 2-BY-2 BLOCKS CORRESPOND TO PAIRS OF COMPLEX
C     EIGENVALUES, AND RETURNS QUANTITIES WHOSE RATIOS GIVE THE
C     GENERALIZED EIGENVALUES.  IT IS USUALLY PRECEDED BY  QZHES
C     AND  QZIT  AND MAY BE FOLLOWED BY  QZVEC.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRICES.
C
C        A CONTAINS A REAL UPPER QUASI-TRIANGULAR MATRIX.
C
C        B CONTAINS A REAL UPPER TRIANGULAR MATRIX.  IN ADDITION,
C          LOCATION B(N,1) CONTAINS THE TOLERANCE QUANTITY (EPSB)
C          COMPUTED AND SAVED IN  QZIT.
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C        Z CONTAINS, IF MATZ HAS BEEN SET TO .TRUE., THE
C          TRANSFORMATION MATRIX PRODUCED IN THE REDUCTIONS BY QZHES
C          AND QZIT, IF PERFORMED, OR ELSE THE IDENTITY MATRIX.
C          IF MATZ HAS BEEN SET TO .FALSE., Z IS NOT REFERENCED.
C
C     ON OUTPUT
C
C        A HAS BEEN REDUCED FURTHER TO A QUASI-TRIANGULAR MATRIX
C          IN WHICH ALL NONZERO SUBDIAGONAL ELEMENTS CORRESPOND TO
C          PAIRS OF COMPLEX EIGENVALUES.
C
C        B IS STILL IN UPPER TRIANGULAR FORM, ALTHOUGH ITS ELEMENTS
C          HAVE BEEN ALTERED.  B(N,1) IS UNALTERED.
C
C        ALFR AND ALFI CONTAIN THE REAL AND IMAGINARY PARTS OF THE
C          DIAGONAL ELEMENTS OF THE TRIANGULAR MATRIX THAT WOULD BE
C          OBTAINED IF A WERE REDUCED COMPLETELY TO TRIANGULAR FORM
C          BY UNITARY TRANSFORMATIONS.  NON-ZERO VALUES OF ALFI OCCUR
C          IN PAIRS, THE FIRST MEMBER POSITIVE AND THE SECOND NEGATIVE.
C
C        BETA CONTAINS THE DIAGONAL ELEMENTS OF THE CORRESPONDING B,
C          NORMALIZED TO BE REAL AND NON-NEGATIVE.  THE GENERALIZED
C          EIGENVALUES ARE THEN THE RATIOS ((ALFR+I*ALFI)/BETA).
C
C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS
C          (FOR ALL THREE STEPS) IF MATZ HAS BEEN SET TO .TRUE.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      EPSB = B(N,1)
      ISW = 1
C     .......... FIND EIGENVALUES OF QUASI-TRIANGULAR MATRICES.
C                FOR EN=N STEP -1 UNTIL 1 DO -- ..........
      DO 510 NN = 1, N
         EN = N + 1 - NN
         NA = EN - 1
         IF (ISW .EQ. 2) GO TO 505
         IF (EN .EQ. 1) GO TO 410
         IF (A(EN,NA) .NE. 0.0E0) GO TO 420
C     .......... 1-BY-1 BLOCK, ONE REAL ROOT ..........
  410    ALFR(EN) = A(EN,EN)
         IF (B(EN,EN) .LT. 0.0E0) ALFR(EN) = -ALFR(EN)
         BETA(EN) = ABS(B(EN,EN))
         ALFI(EN) = 0.0E0
         GO TO 510
C     .......... 2-BY-2 BLOCK ..........
  420    IF (ABS(B(NA,NA)) .LE. EPSB) GO TO 455
         IF (ABS(B(EN,EN)) .GT. EPSB) GO TO 430
         A1 = A(EN,EN)
         A2 = A(EN,NA)
         BN = 0.0E0
         GO TO 435
  430    AN = ABS(A(NA,NA)) + ABS(A(NA,EN)) + ABS(A(EN,NA))
     X      + ABS(A(EN,EN))
         BN = ABS(B(NA,NA)) + ABS(B(NA,EN)) + ABS(B(EN,EN))
         A11 = A(NA,NA) / AN
         A12 = A(NA,EN) / AN
         A21 = A(EN,NA) / AN
         A22 = A(EN,EN) / AN
         B11 = B(NA,NA) / BN
         B12 = B(NA,EN) / BN
         B22 = B(EN,EN) / BN
         E = A11 / B11
         EI = A22 / B22
         S = A21 / (B11 * B22)
         T = (A22 - E * B22) / B22
         IF (ABS(E) .LE. ABS(EI)) GO TO 431
         E = EI
         T = (A11 - E * B11) / B11
  431    C = 0.5E0 * (T - S * B12)
         D = C * C + S * (A12 - E * B12)
         IF (D .LT. 0.0E0) GO TO 480
C     .......... TWO REAL ROOTS.
C                ZERO BOTH A(EN,NA) AND B(EN,NA) ..........
         E = E + (C + SIGN(SQRT(D),C))
         A11 = A11 - E * B11
         A12 = A12 - E * B12
         A22 = A22 - E * B22
         IF (ABS(A11) + ABS(A12) .LT.
     X       ABS(A21) + ABS(A22)) GO TO 432
         A1 = A12
         A2 = A11
         GO TO 435
  432    A1 = A22
         A2 = A21
C     .......... CHOOSE AND APPLY REAL Z ..........
  435    S = ABS(A1) + ABS(A2)
         U1 = A1 / S
         U2 = A2 / S
         R = SIGN(SQRT(U1*U1+U2*U2),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         U2 = V2 / V1
C
         DO 440 I = 1, EN
            T = A(I,EN) + U2 * A(I,NA)
            A(I,EN) = A(I,EN) + T * V1
            A(I,NA) = A(I,NA) + T * V2
            T = B(I,EN) + U2 * B(I,NA)
            B(I,EN) = B(I,EN) + T * V1
            B(I,NA) = B(I,NA) + T * V2
  440    CONTINUE
C
         IF (.NOT. MATZ) GO TO 450
C
         DO 445 I = 1, N
            T = Z(I,EN) + U2 * Z(I,NA)
            Z(I,EN) = Z(I,EN) + T * V1
            Z(I,NA) = Z(I,NA) + T * V2
  445    CONTINUE
C
  450    IF (BN .EQ. 0.0E0) GO TO 475
         IF (AN .LT. ABS(E) * BN) GO TO 455
         A1 = B(NA,NA)
         A2 = B(EN,NA)
         GO TO 460
  455    A1 = A(NA,NA)
         A2 = A(EN,NA)
C     .......... CHOOSE AND APPLY REAL Q ..........
  460    S = ABS(A1) + ABS(A2)
         IF (S .EQ. 0.0E0) GO TO 475
         U1 = A1 / S
         U2 = A2 / S
         R = SIGN(SQRT(U1*U1+U2*U2),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         U2 = V2 / V1
C
         DO 470 J = NA, N
            T = A(NA,J) + U2 * A(EN,J)
            A(NA,J) = A(NA,J) + T * V1
            A(EN,J) = A(EN,J) + T * V2
            T = B(NA,J) + U2 * B(EN,J)
            B(NA,J) = B(NA,J) + T * V1
            B(EN,J) = B(EN,J) + T * V2
  470    CONTINUE
C
  475    A(EN,NA) = 0.0E0
         B(EN,NA) = 0.0E0
         ALFR(NA) = A(NA,NA)
         ALFR(EN) = A(EN,EN)
         IF (B(NA,NA) .LT. 0.0E0) ALFR(NA) = -ALFR(NA)
         IF (B(EN,EN) .LT. 0.0E0) ALFR(EN) = -ALFR(EN)
         BETA(NA) = ABS(B(NA,NA))
         BETA(EN) = ABS(B(EN,EN))
         ALFI(EN) = 0.0E0
         ALFI(NA) = 0.0E0
         GO TO 505
C     .......... TWO COMPLEX ROOTS ..........
  480    E = E + C
         EI = SQRT(-D)
         A11R = A11 - E * B11
         A11I = EI * B11
         A12R = A12 - E * B12
         A12I = EI * B12
         A22R = A22 - E * B22
         A22I = EI * B22
         IF (ABS(A11R) + ABS(A11I) + ABS(A12R) + ABS(A12I) .LT.
     X       ABS(A21) + ABS(A22R) + ABS(A22I)) GO TO 482
         A1 = A12R
         A1I = A12I
         A2 = -A11R
         A2I = -A11I
         GO TO 485
  482    A1 = A22R
         A1I = A22I
         A2 = -A21
         A2I = 0.0E0
C     .......... CHOOSE COMPLEX Z ..........
  485    CZ = SQRT(A1*A1+A1I*A1I)
         IF (CZ .EQ. 0.0E0) GO TO 487
         SZR = (A1 * A2 + A1I * A2I) / CZ
         SZI = (A1 * A2I - A1I * A2) / CZ
         R = SQRT(CZ*CZ+SZR*SZR+SZI*SZI)
         CZ = CZ / R
         SZR = SZR / R
         SZI = SZI / R
         GO TO 490
  487    SZR = 1.0E0
         SZI = 0.0E0
  490    IF (AN .LT. (ABS(E) + EI) * BN) GO TO 492
         A1 = CZ * B11 + SZR * B12
         A1I = SZI * B12
         A2 = SZR * B22
         A2I = SZI * B22
         GO TO 495
  492    A1 = CZ * A11 + SZR * A12
         A1I = SZI * A12
         A2 = CZ * A21 + SZR * A22
         A2I = SZI * A22
C     .......... CHOOSE COMPLEX Q ..........
  495    CQ = SQRT(A1*A1+A1I*A1I)
         IF (CQ .EQ. 0.0E0) GO TO 497
         SQR = (A1 * A2 + A1I * A2I) / CQ
         SQI = (A1 * A2I - A1I * A2) / CQ
         R = SQRT(CQ*CQ+SQR*SQR+SQI*SQI)
         CQ = CQ / R
         SQR = SQR / R
         SQI = SQI / R
         GO TO 500
  497    SQR = 1.0E0
         SQI = 0.0E0
C     .......... COMPUTE DIAGONAL ELEMENTS THAT WOULD RESULT
C                IF TRANSFORMATIONS WERE APPLIED ..........
  500    SSR = SQR * SZR + SQI * SZI
         SSI = SQR * SZI - SQI * SZR
         I = 1
         TR = CQ * CZ * A11 + CQ * SZR * A12 + SQR * CZ * A21
     X      + SSR * A22
         TI = CQ * SZI * A12 - SQI * CZ * A21 + SSI * A22
         DR = CQ * CZ * B11 + CQ * SZR * B12 + SSR * B22
         DI = CQ * SZI * B12 + SSI * B22
         GO TO 503
  502    I = 2
         TR = SSR * A11 - SQR * CZ * A12 - CQ * SZR * A21
     X      + CQ * CZ * A22
         TI = -SSI * A11 - SQI * CZ * A12 + CQ * SZI * A21
         DR = SSR * B11 - SQR * CZ * B12 + CQ * CZ * B22
         DI = -SSI * B11 - SQI * CZ * B12
  503    T = TI * DR - TR * DI
         J = NA
         IF (T .LT. 0.0E0) J = EN
         R = SQRT(DR*DR+DI*DI)
         BETA(J) = BN * R
         ALFR(J) = AN * (TR * DR + TI * DI) / R
         ALFI(J) = AN * T / R
         IF (I .EQ. 1) GO TO 502
  505    ISW = 3 - ISW
  510 CONTINUE
      B(N,1) = EPSB
C
      RETURN
      END
      SUBROUTINE QZVEC(NM,N,A,B,ALFR,ALFI,BETA,Z)
C
      INTEGER I,J,K,M,N,EN,II,JJ,NA,NM,NN,ISW,ENM2
      REAL A(NM,N),B(NM,N),ALFR(N),ALFI(N),BETA(N),Z(NM,N)
      REAL D,Q,R,S,T,W,X,Y,DI,DR,RA,RR,SA,TI,TR,T1,T2,W1,X1,
     X       ZZ,Z1,ALFM,ALMI,ALMR,BETM,EPSB
C
C     THIS SUBROUTINE IS THE OPTIONAL FOURTH STEP OF THE QZ ALGORITHM
C     FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF REAL MATRICES, ONE OF THEM IN
C     QUASI-TRIANGULAR FORM (IN WHICH EACH 2-BY-2 BLOCK CORRESPONDS TO
C     A PAIR OF COMPLEX EIGENVALUES) AND THE OTHER IN UPPER TRIANGULAR
C     FORM.  IT COMPUTES THE EIGENVECTORS OF THE TRIANGULAR PROBLEM AND
C     TRANSFORMS THE RESULTS BACK TO THE ORIGINAL COORDINATE SYSTEM.
C     IT IS USUALLY PRECEDED BY  QZHES,  QZIT, AND  QZVAL.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRICES.
C
C        A CONTAINS A REAL UPPER QUASI-TRIANGULAR MATRIX.
C
C        B CONTAINS A REAL UPPER TRIANGULAR MATRIX.  IN ADDITION,
C          LOCATION B(N,1) CONTAINS THE TOLERANCE QUANTITY (EPSB)
C          COMPUTED AND SAVED IN  QZIT.
C
C        ALFR, ALFI, AND BETA  ARE VECTORS WITH COMPONENTS WHOSE
C          RATIOS ((ALFR+I*ALFI)/BETA) ARE THE GENERALIZED
C          EIGENVALUES.  THEY ARE USUALLY OBTAINED FROM  QZVAL.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTIONS BY  QZHES,  QZIT, AND  QZVAL, IF PERFORMED.
C          IF THE EIGENVECTORS OF THE TRIANGULAR PROBLEM ARE
C          DESIRED, Z MUST CONTAIN THE IDENTITY MATRIX.
C
C     ON OUTPUT
C
C        A IS UNALTERED.  ITS SUBDIAGONAL ELEMENTS PROVIDE INFORMATION
C           ABOUT THE STORAGE OF THE COMPLEX EIGENVECTORS.
C
C        B HAS BEEN DESTROYED.
C
C        ALFR, ALFI, AND BETA ARE UNALTERED.
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
C          IF ALFI(I) .EQ. 0.0, THE I-TH EIGENVALUE IS REAL AND
C            THE I-TH COLUMN OF Z CONTAINS ITS EIGENVECTOR.
C          IF ALFI(I) .NE. 0.0, THE I-TH EIGENVALUE IS COMPLEX.
C            IF ALFI(I) .GT. 0.0, THE EIGENVALUE IS THE FIRST OF
C              A COMPLEX PAIR AND THE I-TH AND (I+1)-TH COLUMNS
C              OF Z CONTAIN ITS EIGENVECTOR.
C            IF ALFI(I) .LT. 0.0, THE EIGENVALUE IS THE SECOND OF
C              A COMPLEX PAIR AND THE (I-1)-TH AND I-TH COLUMNS
C              OF Z CONTAIN THE CONJUGATE OF ITS EIGENVECTOR.
C          EACH EIGENVECTOR IS NORMALIZED SO THAT THE MODULUS
C          OF ITS LARGEST COMPONENT IS 1.0 .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      EPSB = B(N,1)
      ISW = 1
C     .......... FOR EN=N STEP -1 UNTIL 1 DO -- ..........
      DO 800 NN = 1, N
         EN = N + 1 - NN
         NA = EN - 1
         IF (ISW .EQ. 2) GO TO 795
         IF (ALFI(EN) .NE. 0.0E0) GO TO 710
C     .......... REAL VECTOR ..........
         M = EN
         B(EN,EN) = 1.0E0
         IF (NA .EQ. 0) GO TO 800
         ALFM = ALFR(M)
         BETM = BETA(M)
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
         DO 700 II = 1, NA
            I = EN - II
            W = BETM * A(I,I) - ALFM * B(I,I)
            R = 0.0E0
C
            DO 610 J = M, EN
  610       R = R + (BETM * A(I,J) - ALFM * B(I,J)) * B(J,EN)
C
            IF (I .EQ. 1 .OR. ISW .EQ. 2) GO TO 630
            IF (BETM * A(I,I-1) .EQ. 0.0E0) GO TO 630
            ZZ = W
            S = R
            GO TO 690
  630       M = I
            IF (ISW .EQ. 2) GO TO 640
C     .......... REAL 1-BY-1 BLOCK ..........
            T = W
            IF (W .EQ. 0.0E0) T = EPSB
            B(I,EN) = -R / T
            GO TO 700
C     .......... REAL 2-BY-2 BLOCK ..........
  640       X = BETM * A(I,I+1) - ALFM * B(I,I+1)
            Y = BETM * A(I+1,I)
            Q = W * ZZ - X * Y
            T = (X * S - ZZ * R) / Q
            B(I,EN) = T
            IF (ABS(X) .LE. ABS(ZZ)) GO TO 650
            B(I+1,EN) = (-R - W * T) / X
            GO TO 690
  650       B(I+1,EN) = (-S - Y * T) / ZZ
  690       ISW = 3 - ISW
  700    CONTINUE
C     .......... END REAL VECTOR ..........
         GO TO 800
C     .......... COMPLEX VECTOR ..........
  710    M = NA
         ALMR = ALFR(M)
         ALMI = ALFI(M)
         BETM = BETA(M)
C     .......... LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
C                EIGENVECTOR MATRIX IS TRIANGULAR ..........
         Y = BETM * A(EN,NA)
         B(NA,NA) = -ALMI * B(EN,EN) / Y
         B(NA,EN) = (ALMR * B(EN,EN) - BETM * A(EN,EN)) / Y
         B(EN,NA) = 0.0E0
         B(EN,EN) = 1.0E0
         ENM2 = NA - 1
         IF (ENM2 .EQ. 0) GO TO 795
C     .......... FOR I=EN-2 STEP -1 UNTIL 1 DO -- ..........
         DO 790 II = 1, ENM2
            I = NA - II
            W = BETM * A(I,I) - ALMR * B(I,I)
            W1 = -ALMI * B(I,I)
            RA = 0.0E0
            SA = 0.0E0
C
            DO 760 J = M, EN
               X = BETM * A(I,J) - ALMR * B(I,J)
               X1 = -ALMI * B(I,J)
               RA = RA + X * B(J,NA) - X1 * B(J,EN)
               SA = SA + X * B(J,EN) + X1 * B(J,NA)
  760       CONTINUE
C
            IF (I .EQ. 1 .OR. ISW .EQ. 2) GO TO 770
            IF (BETM * A(I,I-1) .EQ. 0.0E0) GO TO 770
            ZZ = W
            Z1 = W1
            R = RA
            S = SA
            ISW = 2
            GO TO 790
  770       M = I
            IF (ISW .EQ. 2) GO TO 780
C     .......... COMPLEX 1-BY-1 BLOCK ..........
            TR = -RA
            TI = -SA
  773       DR = W
            DI = W1
C     .......... COMPLEX DIVIDE (T1,T2) = (TR,TI) / (DR,DI) ..........
  775       IF (ABS(DI) .GT. ABS(DR)) GO TO 777
            RR = DI / DR
            D = DR + DI * RR
            T1 = (TR + TI * RR) / D
            T2 = (TI - TR * RR) / D
            GO TO (787,782), ISW
  777       RR = DR / DI
            D = DR * RR + DI
            T1 = (TR * RR + TI) / D
            T2 = (TI * RR - TR) / D
            GO TO (787,782), ISW
C     .......... COMPLEX 2-BY-2 BLOCK ..........
  780       X = BETM * A(I,I+1) - ALMR * B(I,I+1)
            X1 = -ALMI * B(I,I+1)
            Y = BETM * A(I+1,I)
            TR = Y * RA - W * R + W1 * S
            TI = Y * SA - W * S - W1 * R
            DR = W * ZZ - W1 * Z1 - X * Y
            DI = W * Z1 + W1 * ZZ - X1 * Y
            IF (DR .EQ. 0.0E0 .AND. DI .EQ. 0.0E0) DR = EPSB
            GO TO 775
  782       B(I+1,NA) = T1
            B(I+1,EN) = T2
            ISW = 1
            IF (ABS(Y) .GT. ABS(W) + ABS(W1)) GO TO 785
            TR = -RA - X * B(I+1,NA) + X1 * B(I+1,EN)
            TI = -SA - X * B(I+1,EN) - X1 * B(I+1,NA)
            GO TO 773
  785       T1 = (-R - ZZ * B(I+1,NA) + Z1 * B(I+1,EN)) / Y
            T2 = (-S - ZZ * B(I+1,EN) - Z1 * B(I+1,NA)) / Y
  787       B(I,NA) = T1
            B(I,EN) = T2
  790    CONTINUE
C     .......... END COMPLEX VECTOR ..........
  795    ISW = 3 - ISW
  800 CONTINUE
C     .......... END BACK SUBSTITUTION.
C                TRANSFORM TO ORIGINAL COORDINATE SYSTEM.
C                FOR J=N STEP -1 UNTIL 1 DO -- ..........
      DO 880 JJ = 1, N
         J = N + 1 - JJ
C
         DO 880 I = 1, N
            ZZ = 0.0E0
C
            DO 860 K = 1, J
  860       ZZ = ZZ + Z(I,K) * B(K,J)
C
            Z(I,J) = ZZ
  880 CONTINUE
C     .......... NORMALIZE SO THAT MODULUS OF LARGEST
C                COMPONENT OF EACH VECTOR IS 1.
C                (ISW IS 1 INITIALLY FROM BEFORE) ..........
      DO 950 J = 1, N
         D = 0.0E0
         IF (ISW .EQ. 2) GO TO 920
         IF (ALFI(J) .NE. 0.0E0) GO TO 945
C
         DO 890 I = 1, N
            IF (ABS(Z(I,J)) .GT. D) D = ABS(Z(I,J))
  890    CONTINUE
C
         DO 900 I = 1, N
  900    Z(I,J) = Z(I,J) / D
C
         GO TO 950
C
  920    DO 930 I = 1, N
            R = ABS(Z(I,J-1)) + ABS(Z(I,J))
            IF (R .NE. 0.0E0) R = R * SQRT((Z(I,J-1)/R)**2
     X                                     +(Z(I,J)/R)**2)
            IF (R .GT. D) D = R
  930    CONTINUE
C
         DO 940 I = 1, N
            Z(I,J-1) = Z(I,J-1) / D
            Z(I,J) = Z(I,J) / D
  940    CONTINUE
C
  945    ISW = 3 - ISW
  950 CONTINUE
C
      RETURN
      END
      SUBROUTINE RATQR(N,EPS1,D,E,E2,M,W,IND,BD,TYPE,IDEF,IERR)
C
      INTEGER I,J,K,M,N,II,JJ,K1,IDEF,IERR,JDEF
      REAL D(N),E(N),E2(N),W(N),BD(N)
      REAL F,P,Q,R,S,EP,QP,ERR,TOT,EPS1,DELTA,EPSLON
      INTEGER IND(N)
      LOGICAL TYPE
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE RATQR,
C     NUM. MATH. 11, 264-272(1968) BY REINSCH AND BAUER.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 257-265(1971).
C
C     THIS SUBROUTINE FINDS THE ALGEBRAICALLY SMALLEST OR LARGEST
C     EIGENVALUES OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE
C     RATIONAL QR METHOD WITH NEWTON CORRECTIONS.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        EPS1 IS A THEORETICAL ABSOLUTE ERROR TOLERANCE FOR THE
C          COMPUTED EIGENVALUES.  IF THE INPUT EPS1 IS NON-POSITIVE,
C          OR INDEED SMALLER THAN ITS DEFAULT VALUE, IT IS RESET
C          AT EACH ITERATION TO THE RESPECTIVE DEFAULT VALUE,
C          NAMELY, THE PRODUCT OF THE RELATIVE MACHINE PRECISION
C          AND THE MAGNITUDE OF THE CURRENT EIGENVALUE ITERATE.
C          THE THEORETICAL ABSOLUTE ERROR IN THE K-TH EIGENVALUE
C          IS USUALLY NOT GREATER THAN K TIMES EPS1.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY.
C
C        M IS THE NUMBER OF EIGENVALUES TO BE FOUND.
C
C        IDEF SHOULD BE SET TO 1 IF THE INPUT MATRIX IS KNOWN TO BE
C          POSITIVE DEFINITE, TO -1 IF THE INPUT MATRIX IS KNOWN TO
C          BE NEGATIVE DEFINITE, AND TO 0 OTHERWISE.
C
C        TYPE SHOULD BE SET TO .TRUE. IF THE SMALLEST EIGENVALUES
C          ARE TO BE FOUND, AND TO .FALSE. IF THE LARGEST EIGENVALUES
C          ARE TO BE FOUND.
C
C     ON OUTPUT
C
C        EPS1 IS UNALTERED UNLESS IT HAS BEEN RESET TO ITS
C          (LAST) DEFAULT VALUE.
C
C        D AND E ARE UNALTERED (UNLESS W OVERWRITES D).
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS SET TO 0.0E0 IF THE SMALLEST EIGENVALUES HAVE BEEN
C          FOUND, AND TO 2.0E0 IF THE LARGEST EIGENVALUES HAVE BEEN
C          FOUND.  E2 IS OTHERWISE UNALTERED (UNLESS OVERWRITTEN BY BD).
C
C        W CONTAINS THE M ALGEBRAICALLY SMALLEST EIGENVALUES IN
C          ASCENDING ORDER, OR THE M LARGEST EIGENVALUES IN
C          DESCENDING ORDER.  IF AN ERROR EXIT IS MADE BECAUSE OF
C          AN INCORRECT SPECIFICATION OF IDEF, NO EIGENVALUES
C          ARE FOUND.  IF THE NEWTON ITERATES FOR A PARTICULAR
C          EIGENVALUE ARE NOT MONOTONE, THE BEST ESTIMATE OBTAINED
C          IS RETURNED AND IERR IS SET.  W MAY COINCIDE WITH D.
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC..
C
C        BD CONTAINS REFINED BOUNDS FOR THE THEORETICAL ERRORS OF THE
C          CORRESPONDING EIGENVALUES IN W.  THESE BOUNDS ARE USUALLY
C          WITHIN THE TOLERANCE SPECIFIED BY EPS1.  BD MAY COINCIDE
C          WITH E2.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          6*N+1      IF  IDEF  IS SET TO 1 AND  TYPE  TO .TRUE.
C                     WHEN THE MATRIX IS NOT POSITIVE DEFINITE, OR
C                     IF  IDEF  IS SET TO -1 AND  TYPE  TO .FALSE.
C                     WHEN THE MATRIX IS NOT NEGATIVE DEFINITE,
C          5*N+K      IF SUCCESSIVE ITERATES TO THE K-TH EIGENVALUE
C                     ARE NOT MONOTONE INCREASING, WHERE K REFERS
C                     TO THE LAST SUCH OCCURRENCE.
C
C     NOTE THAT SUBROUTINE TRIDIB IS GENERALLY FASTER AND MORE
C     ACCURATE THAN RATQR IF THE EIGENVALUES ARE CLUSTERED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      JDEF = IDEF
C     .......... COPY D ARRAY INTO W ..........
      DO 20 I = 1, N
   20 W(I) = D(I)
C
      IF (TYPE) GO TO 40
      J = 1
      GO TO 400
   40 ERR = 0.0E0
      S = 0.0E0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DEFINE
C                INITIAL SHIFT FROM LOWER GERSCHGORIN BOUND.
C                COPY E2 ARRAY INTO BD ..........
      TOT = W(1)
      Q = 0.0E0
      J = 0
C
      DO 100 I = 1, N
         P = Q
         IF (I .EQ. 1) GO TO 60
         IF (P .GT. EPSLON(ABS(D(I)) + ABS(D(I-1)))) GO TO 80
   60    E2(I) = 0.0E0
   80    BD(I) = E2(I)
C     .......... COUNT ALSO IF ELEMENT OF E2 HAS UNDERFLOWED ..........
         IF (E2(I) .EQ. 0.0E0) J = J + 1
         IND(I) = J
         Q = 0.0E0
         IF (I .NE. N) Q = ABS(E(I+1))
         TOT = AMIN1(W(I)-P-Q,TOT)
  100 CONTINUE
C
      IF (JDEF .EQ. 1 .AND. TOT .LT. 0.0E0) GO TO 140
C
      DO 110 I = 1, N
  110 W(I) = W(I) - TOT
C
      GO TO 160
  140 TOT = 0.0E0
C
  160 DO 360 K = 1, M
C     .......... NEXT QR TRANSFORMATION ..........
  180    TOT = TOT + S
         DELTA = W(N) - S
         I = N
         F = ABS(EPSLON(TOT))
         IF (EPS1 .LT. F) EPS1 = F
         IF (DELTA .GT. EPS1) GO TO 190
         IF (DELTA .LT. (-EPS1)) GO TO 1000
         GO TO 300
C     .......... REPLACE SMALL SUB-DIAGONAL SQUARES BY ZERO
C                TO REDUCE THE INCIDENCE OF UNDERFLOWS ..........
  190    IF (K .EQ. N) GO TO 210
         K1 = K + 1
         DO 200 J = K1, N
            IF (BD(J) .LE. (EPSLON(W(J)+W(J-1))) ** 2) BD(J) = 0.0E0
  200    CONTINUE
C
  210    F = BD(N) / DELTA
         QP = DELTA + F
         P = 1.0E0
         IF (K .EQ. N) GO TO 260
         K1 = N - K
C     .......... FOR I=N-1 STEP -1 UNTIL K DO -- ..........
         DO 240 II = 1, K1
            I = N - II
            Q = W(I) - S - F
            R = Q / QP
            P = P * R + 1.0E0
            EP = F * R
            W(I+1) = QP + EP
            DELTA = Q - EP
            IF (DELTA .GT. EPS1) GO TO 220
            IF (DELTA .LT. (-EPS1)) GO TO 1000
            GO TO 300
  220       F = BD(I) / Q
            QP = DELTA + F
            BD(I+1) = QP * EP
  240    CONTINUE
C
  260    W(K) = QP
         S = QP / P
         IF (TOT + S .GT. TOT) GO TO 180
C     .......... SET ERROR -- IRREGULAR END OF ITERATION.
C                DEFLATE MINIMUM DIAGONAL ELEMENT ..........
         IERR = 5 * N + K
         S = 0.0E0
         DELTA = QP
C
         DO 280 J = K, N
            IF (W(J) .GT. DELTA) GO TO 280
            I = J
            DELTA = W(J)
  280    CONTINUE
C     .......... CONVERGENCE ..........
  300    IF (I .LT. N) BD(I+1) = BD(I) * F / QP
         II = IND(I)
         IF (I .EQ. K) GO TO 340
         K1 = I - K
C     .......... FOR J=I-1 STEP -1 UNTIL K DO -- ..........
         DO 320 JJ = 1, K1
            J = I - JJ
            W(J+1) = W(J) - S
            BD(J+1) = BD(J)
            IND(J+1) = IND(J)
  320    CONTINUE
C
  340    W(K) = TOT
         ERR = ERR + ABS(DELTA)
         BD(K) = ERR
         IND(K) = II
  360 CONTINUE
C
      IF (TYPE) GO TO 1001
      F = BD(1)
      E2(1) = 2.0E0
      BD(1) = F
      J = 2
C     .......... NEGATE ELEMENTS OF W FOR LARGEST VALUES ..........
  400 DO 500 I = 1, N
  500 W(I) = -W(I)
C
      JDEF = -JDEF
      GO TO (40,1001), J
C     .......... SET ERROR -- IDEF SPECIFIED INCORRECTLY ..........
 1000 IERR = 6 * N + 1
 1001 RETURN
      END
      SUBROUTINE REBAK(NM,N,B,DL,M,Z)
C
      INTEGER I,J,K,M,N,I1,II,NM
      REAL B(NM,N),DL(N),Z(NM,M)
      REAL X
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE REBAKA,
C     NUM. MATH. 11, 99-110(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 303-314(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A GENERALIZED
C     SYMMETRIC EIGENSYSTEM BY BACK TRANSFORMING THOSE OF THE
C     DERIVED SYMMETRIC MATRIX DETERMINED BY  REDUC.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX SYSTEM.
C
C        B CONTAINS INFORMATION ABOUT THE SIMILARITY TRANSFORMATION
C          (CHOLESKY DECOMPOSITION) USED IN THE REDUCTION BY  REDUC
C          IN ITS STRICT LOWER TRIANGLE.
C
C        DL CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATION.
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
C
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS
C          IN ITS FIRST M COLUMNS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
C
      DO 100 J = 1, M
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
         DO 100 II = 1, N
            I = N + 1 - II
            I1 = I + 1
            X = Z(I,J)
            IF (I .EQ. N) GO TO 80
C
            DO 60 K = I1, N
   60       X = X - B(K,I) * Z(K,J)
C
   80       Z(I,J) = X / DL(I)
  100 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE REBAKB(NM,N,B,DL,M,Z)
C
      INTEGER I,J,K,M,N,I1,II,NM
      REAL B(NM,N),DL(N),Z(NM,M)
      REAL X
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE REBAKB,
C     NUM. MATH. 11, 99-110(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 303-314(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A GENERALIZED
C     SYMMETRIC EIGENSYSTEM BY BACK TRANSFORMING THOSE OF THE
C     DERIVED SYMMETRIC MATRIX DETERMINED BY  REDUC2.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX SYSTEM.
C
C        B CONTAINS INFORMATION ABOUT THE SIMILARITY TRANSFORMATION
C          (CHOLESKY DECOMPOSITION) USED IN THE REDUCTION BY  REDUC2
C          IN ITS STRICT LOWER TRIANGLE.
C
C        DL CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATION.
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
C
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS
C          IN ITS FIRST M COLUMNS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
C
      DO 100 J = 1, M
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
         DO 100 II = 1, N
            I1 = N - II
            I = I1 + 1
            X = DL(I) * Z(I,J)
            IF (I .EQ. 1) GO TO 80
C
            DO 60 K = 1, I1
   60       X = X + B(I,K) * Z(K,J)
C
   80       Z(I,J) = X
  100 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE REDUC(NM,N,A,B,DL,IERR)
C
      INTEGER I,J,K,N,I1,J1,NM,NN,IERR
      REAL A(NM,N),B(NM,N),DL(N)
      REAL X,Y
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE REDUC1,
C     NUM. MATH. 11, 99-110(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 303-314(1971).
C
C     THIS SUBROUTINE REDUCES THE GENERALIZED SYMMETRIC EIGENPROBLEM
C     AX=(LAMBDA)BX, WHERE B IS POSITIVE DEFINITE, TO THE STANDARD
C     SYMMETRIC EIGENPROBLEM USING THE CHOLESKY FACTORIZATION OF B.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRICES A AND B.  IF THE CHOLESKY
C          FACTOR L OF B IS ALREADY AVAILABLE, N SHOULD BE PREFIXED
C          WITH A MINUS SIGN.
C
C        A AND B CONTAIN THE REAL SYMMETRIC INPUT MATRICES.  ONLY THE
C          FULL UPPER TRIANGLES OF THE MATRICES NEED BE SUPPLIED.  IF
C          N IS NEGATIVE, THE STRICT LOWER TRIANGLE OF B CONTAINS,
C          INSTEAD, THE STRICT LOWER TRIANGLE OF ITS CHOLESKY FACTOR L.
C
C        DL CONTAINS, IF N IS NEGATIVE, THE DIAGONAL ELEMENTS OF L.
C
C     ON OUTPUT
C
C        A CONTAINS IN ITS FULL LOWER TRIANGLE THE FULL LOWER TRIANGLE
C          OF THE SYMMETRIC MATRIX DERIVED FROM THE REDUCTION TO THE
C          STANDARD FORM.  THE STRICT UPPER TRIANGLE OF A IS UNALTERED.
C
C        B CONTAINS IN ITS STRICT LOWER TRIANGLE THE STRICT LOWER
C          TRIANGLE OF ITS CHOLESKY FACTOR L.  THE FULL UPPER
C          TRIANGLE OF B IS UNALTERED.
C
C        DL CONTAINS THE DIAGONAL ELEMENTS OF L.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          7*N+1      IF B IS NOT POSITIVE DEFINITE.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      NN = IABS(N)
      IF (N .LT. 0) GO TO 100
C     .......... FORM L IN THE ARRAYS B AND DL ..........
      DO 80 I = 1, N
         I1 = I - 1
C
         DO 80 J = I, N
            X = B(I,J)
            IF (I .EQ. 1) GO TO 40
C
            DO 20 K = 1, I1
   20       X = X - B(I,K) * B(J,K)
C
   40       IF (J .NE. I) GO TO 60
            IF (X .LE. 0.0E0) GO TO 1000
            Y = SQRT(X)
            DL(I) = Y
            GO TO 80
   60       B(J,I) = X / Y
   80 CONTINUE
C     .......... FORM THE TRANSPOSE OF THE UPPER TRIANGLE OF INV(L)*A
C                IN THE LOWER TRIANGLE OF THE ARRAY A ..........
  100 DO 200 I = 1, NN
         I1 = I - 1
         Y = DL(I)
C
         DO 200 J = I, NN
            X = A(I,J)
            IF (I .EQ. 1) GO TO 180
C
            DO 160 K = 1, I1
  160       X = X - B(I,K) * A(J,K)
C
  180       A(J,I) = X / Y
  200 CONTINUE
C     .......... PRE-MULTIPLY BY INV(L) AND OVERWRITE ..........
      DO 300 J = 1, NN
         J1 = J - 1
C
         DO 300 I = J, NN
            X = A(I,J)
            IF (I .EQ. J) GO TO 240
            I1 = I - 1
C
            DO 220 K = J, I1
  220       X = X - A(K,J) * B(I,K)
C
  240       IF (J .EQ. 1) GO TO 280
C
            DO 260 K = 1, J1
  260       X = X - A(J,K) * B(I,K)
C
  280       A(I,J) = X / DL(I)
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- B IS NOT POSITIVE DEFINITE ..........
 1000 IERR = 7 * N + 1
 1001 RETURN
      END
      SUBROUTINE REDUC2(NM,N,A,B,DL,IERR)
C
      INTEGER I,J,K,N,I1,J1,NM,NN,IERR
      REAL A(NM,N),B(NM,N),DL(N)
      REAL X,Y
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE REDUC2,
C     NUM. MATH. 11, 99-110(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 303-314(1971).
C
C     THIS SUBROUTINE REDUCES THE GENERALIZED SYMMETRIC EIGENPROBLEMS
C     ABX=(LAMBDA)X OR BAY=(LAMBDA)Y, WHERE B IS POSITIVE DEFINITE,
C     TO THE STANDARD SYMMETRIC EIGENPROBLEM USING THE CHOLESKY
C     FACTORIZATION OF B.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRICES A AND B.  IF THE CHOLESKY
C          FACTOR L OF B IS ALREADY AVAILABLE, N SHOULD BE PREFIXED
C          WITH A MINUS SIGN.
C
C        A AND B CONTAIN THE REAL SYMMETRIC INPUT MATRICES.  ONLY THE
C          FULL UPPER TRIANGLES OF THE MATRICES NEED BE SUPPLIED.  IF
C          N IS NEGATIVE, THE STRICT LOWER TRIANGLE OF B CONTAINS,
C          INSTEAD, THE STRICT LOWER TRIANGLE OF ITS CHOLESKY FACTOR L.
C
C        DL CONTAINS, IF N IS NEGATIVE, THE DIAGONAL ELEMENTS OF L.
C
C     ON OUTPUT
C
C        A CONTAINS IN ITS FULL LOWER TRIANGLE THE FULL LOWER TRIANGLE
C          OF THE SYMMETRIC MATRIX DERIVED FROM THE REDUCTION TO THE
C          STANDARD FORM.  THE STRICT UPPER TRIANGLE OF A IS UNALTERED.
C
C        B CONTAINS IN ITS STRICT LOWER TRIANGLE THE STRICT LOWER
C          TRIANGLE OF ITS CHOLESKY FACTOR L.  THE FULL UPPER
C          TRIANGLE OF B IS UNALTERED.
C
C        DL CONTAINS THE DIAGONAL ELEMENTS OF L.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          7*N+1      IF B IS NOT POSITIVE DEFINITE.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      NN = IABS(N)
      IF (N .LT. 0) GO TO 100
C     .......... FORM L IN THE ARRAYS B AND DL ..........
      DO 80 I = 1, N
         I1 = I - 1
C
         DO 80 J = I, N
            X = B(I,J)
            IF (I .EQ. 1) GO TO 40
C
            DO 20 K = 1, I1
   20       X = X - B(I,K) * B(J,K)
C
   40       IF (J .NE. I) GO TO 60
            IF (X .LE. 0.0E0) GO TO 1000
            Y = SQRT(X)
            DL(I) = Y
            GO TO 80
   60       B(J,I) = X / Y
   80 CONTINUE
C     .......... FORM THE LOWER TRIANGLE OF A*L
C                IN THE LOWER TRIANGLE OF THE ARRAY A ..........
  100 DO 200 I = 1, NN
         I1 = I + 1
C
         DO 200 J = 1, I
            X = A(J,I) * DL(J)
            IF (J .EQ. I) GO TO 140
            J1 = J + 1
C
            DO 120 K = J1, I
  120       X = X + A(K,I) * B(K,J)
C
  140       IF (I .EQ. NN) GO TO 180
C
            DO 160 K = I1, NN
  160       X = X + A(I,K) * B(K,J)
C
  180       A(I,J) = X
  200 CONTINUE
C     .......... PRE-MULTIPLY BY TRANSPOSE(L) AND OVERWRITE ..........
      DO 300 I = 1, NN
         I1 = I + 1
         Y = DL(I)
C
         DO 300 J = 1, I
            X = Y * A(I,J)
            IF (I .EQ. NN) GO TO 280
C
            DO 260 K = I1, NN
  260       X = X + A(K,J) * B(K,I)
C
  280       A(I,J) = X
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- B IS NOT POSITIVE DEFINITE ..........
 1000 IERR = 7 * N + 1
 1001 RETURN
      END
      SUBROUTINE RG(NM,N,A,WR,WI,MATZ,Z,IV1,FV1,IERR)
C
      INTEGER N,NM,IS1,IS2,IERR,MATZ
      REAL A(NM,N),WR(N),WI(N),Z(NM,N),FV1(N)
      INTEGER IV1(N)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     OF A REAL GENERAL MATRIX.
C
C     ON INPUT
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT.
C
C        N  IS THE ORDER OF THE MATRIX  A.
C
C        A  CONTAINS THE REAL GENERAL MATRIX.
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT
C
C        WR  AND  WI  CONTAIN THE REAL AND IMAGINARY PARTS,
C        RESPECTIVELY, OF THE EIGENVALUES.  COMPLEX CONJUGATE
C        PAIRS OF EIGENVALUES APPEAR CONSECUTIVELY WITH THE
C        EIGENVALUE HAVING THE POSITIVE IMAGINARY PART FIRST.
C
C        Z  CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS
C        IF MATZ IS NOT ZERO.  IF THE J-TH EIGENVALUE IS REAL, THE
C        J-TH COLUMN OF  Z  CONTAINS ITS EIGENVECTOR.  IF THE J-TH
C        EIGENVALUE IS COMPLEX WITH POSITIVE IMAGINARY PART, THE
C        J-TH AND (J+1)-TH COLUMNS OF  Z  CONTAIN THE REAL AND
C        IMAGINARY PARTS OF ITS EIGENVECTOR.  THE CONJUGATE OF THIS
C        VECTOR IS THE EIGENVECTOR FOR THE CONJUGATE EIGENVALUE.
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR HQR
C           AND HQR2.  THE NORMAL COMPLETION CODE IS ZERO.
C
C        IV1  AND  FV1  ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 CALL  BALANC(NM,N,A,IS1,IS2,FV1)
      CALL  ELMHES(NM,N,IS1,IS2,A,IV1)
      IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  HQR(NM,N,IS1,IS2,A,WR,WI,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  ELTRAN(NM,N,IS1,IS2,A,IV1,Z)
      CALL  HQR2(NM,N,IS1,IS2,A,WR,WI,Z,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  BALBAK(NM,N,IS1,IS2,FV1,N,Z)
   50 RETURN
      END
      SUBROUTINE RGG(NM,N,A,B,ALFR,ALFI,BETA,MATZ,Z,IERR)
C
      INTEGER N,NM,IERR,MATZ
      REAL A(NM,N),B(NM,N),ALFR(N),ALFI(N),BETA(N),Z(NM,N)
      LOGICAL TF
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     FOR THE REAL GENERAL GENERALIZED EIGENPROBLEM  AX = (LAMBDA)BX.
C
C     ON INPUT
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT.
C
C        N  IS THE ORDER OF THE MATRICES  A  AND  B.
C
C        A  CONTAINS A REAL GENERAL MATRIX.
C
C        B  CONTAINS A REAL GENERAL MATRIX.
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT
C
C        ALFR  AND  ALFI  CONTAIN THE REAL AND IMAGINARY PARTS,
C        RESPECTIVELY, OF THE NUMERATORS OF THE EIGENVALUES.
C
C        BETA  CONTAINS THE DENOMINATORS OF THE EIGENVALUES,
C        WHICH ARE THUS GIVEN BY THE RATIOS  (ALFR+I*ALFI)/BETA.
C        COMPLEX CONJUGATE PAIRS OF EIGENVALUES APPEAR CONSECUTIVELY
C        WITH THE EIGENVALUE HAVING THE POSITIVE IMAGINARY PART FIRST.
C
C        Z  CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS
C        IF MATZ IS NOT ZERO.  IF THE J-TH EIGENVALUE IS REAL, THE
C        J-TH COLUMN OF  Z  CONTAINS ITS EIGENVECTOR.  IF THE J-TH
C        EIGENVALUE IS COMPLEX WITH POSITIVE IMAGINARY PART, THE
C        J-TH AND (J+1)-TH COLUMNS OF  Z  CONTAIN THE REAL AND
C        IMAGINARY PARTS OF ITS EIGENVECTOR.  THE CONJUGATE OF THIS
C        VECTOR IS THE EIGENVECTOR FOR THE CONJUGATE EIGENVALUE.
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR QZIT.
C           THE NORMAL COMPLETION CODE IS ZERO.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      TF = .FALSE.
      CALL  QZHES(NM,N,A,B,TF,Z)
      CALL  QZIT(NM,N,A,B,0.0E0,TF,Z,IERR)
      CALL  QZVAL(NM,N,A,B,ALFR,ALFI,BETA,TF,Z)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 TF = .TRUE.
      CALL  QZHES(NM,N,A,B,TF,Z)
      CALL  QZIT(NM,N,A,B,0.0E0,TF,Z,IERR)
      CALL  QZVAL(NM,N,A,B,ALFR,ALFI,BETA,TF,Z)
      IF (IERR .NE. 0) GO TO 50
      CALL  QZVEC(NM,N,A,B,ALFR,ALFI,BETA,Z)
   50 RETURN
      END
      SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
C
      INTEGER N,NM,IERR,MATZ
      REAL A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     OF A REAL SYMMETRIC MATRIX.
C
C     ON INPUT
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT.
C
C        N  IS THE ORDER OF THE MATRIX  A.
C
C        A  CONTAINS THE REAL SYMMETRIC MATRIX.
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT
C
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.
C
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT
C           AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO.
C
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  TRED1(NM,N,A,W,FV1,FV2)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
      END
      SUBROUTINE RSB(NM,N,MB,A,W,MATZ,Z,FV1,FV2,IERR)
C
      INTEGER N,MB,NM,IERR,MATZ
      REAL A(NM,MB),W(N),Z(NM,N),FV1(N),FV2(N)
      LOGICAL TF
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     OF A REAL SYMMETRIC BAND MATRIX.
C
C     ON INPUT
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT.
C
C        N  IS THE ORDER OF THE MATRIX  A.
C
C        MB  IS THE HALF BAND WIDTH OF THE MATRIX, DEFINED AS THE
C        NUMBER OF ADJACENT DIAGONALS, INCLUDING THE PRINCIPAL
C        DIAGONAL, REQUIRED TO SPECIFY THE NON-ZERO PORTION OF THE
C        LOWER TRIANGLE OF THE MATRIX.
C
C        A  CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
C        BAND MATRIX.  ITS LOWEST SUBDIAGONAL IS STORED IN THE
C        LAST  N+1-MB  POSITIONS OF THE FIRST COLUMN, ITS NEXT
C        SUBDIAGONAL IN THE LAST  N+2-MB  POSITIONS OF THE
C        SECOND COLUMN, FURTHER SUBDIAGONALS SIMILARLY, AND
C        FINALLY ITS PRINCIPAL DIAGONAL IN THE  N  POSITIONS
C        OF THE LAST COLUMN.  CONTENTS OF STORAGES NOT PART
C        OF THE MATRIX ARE ARBITRARY.
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT
C
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.
C
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT
C           AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO.
C
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (N .LE. NM) GO TO 5
      IERR = 10 * N
      GO TO 50
    5 IF (MB .GT. 0) GO TO 10
      IERR = 12 * N
      GO TO 50
   10 IF (MB .LE. N) GO TO 15
      IERR = 12 * N
      GO TO 50
C
   15 IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      TF = .FALSE.
      CALL  BANDR(NM,N,MB,A,W,FV1,FV2,TF,Z)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 TF = .TRUE.
      CALL  BANDR(NM,N,MB,A,W,FV1,FV1,TF,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
      END
      SUBROUTINE RSG(NM,N,A,B,W,MATZ,Z,FV1,FV2,IERR)
C
      INTEGER N,NM,IERR,MATZ
      REAL A(NM,N),B(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     FOR THE REAL SYMMETRIC GENERALIZED EIGENPROBLEM  AX = (LAMBDA)BX.
C
C     ON INPUT
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT.
C
C        N  IS THE ORDER OF THE MATRICES  A  AND  B.
C
C        A  CONTAINS A REAL SYMMETRIC MATRIX.
C
C        B  CONTAINS A POSITIVE DEFINITE REAL SYMMETRIC MATRIX.
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT
C
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.
C
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT
C           AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO.
C
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 CALL  REDUC(NM,N,A,B,FV2,IERR)
      IF (IERR .NE. 0) GO TO 50
      IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  TRED1(NM,N,A,W,FV1,FV2)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  REBAK(NM,N,B,FV2,N,Z)
   50 RETURN
      END
      SUBROUTINE RSGAB(NM,N,A,B,W,MATZ,Z,FV1,FV2,IERR)
C
      INTEGER N,NM,IERR,MATZ
      REAL A(NM,N),B(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     FOR THE REAL SYMMETRIC GENERALIZED EIGENPROBLEM  ABX = (LAMBDA)X.
C
C     ON INPUT
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT.
C
C        N  IS THE ORDER OF THE MATRICES  A  AND  B.
C
C        A  CONTAINS A REAL SYMMETRIC MATRIX.
C
C        B  CONTAINS A POSITIVE DEFINITE REAL SYMMETRIC MATRIX.
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT
C
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.
C
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT
C           AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO.
C
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 CALL  REDUC2(NM,N,A,B,FV2,IERR)
      IF (IERR .NE. 0) GO TO 50
      IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  TRED1(NM,N,A,W,FV1,FV2)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  REBAK(NM,N,B,FV2,N,Z)
   50 RETURN
      END
      SUBROUTINE RSGBA(NM,N,A,B,W,MATZ,Z,FV1,FV2,IERR)
C
      INTEGER N,NM,IERR,MATZ
      REAL A(NM,N),B(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     FOR THE REAL SYMMETRIC GENERALIZED EIGENPROBLEM  BAX = (LAMBDA)X.
C
C     ON INPUT
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT.
C
C        N  IS THE ORDER OF THE MATRICES  A  AND  B.
C
C        A  CONTAINS A REAL SYMMETRIC MATRIX.
C
C        B  CONTAINS A POSITIVE DEFINITE REAL SYMMETRIC MATRIX.
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT
C
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.
C
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT
C           AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO.
C
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 CALL  REDUC2(NM,N,A,B,FV2,IERR)
      IF (IERR .NE. 0) GO TO 50
      IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  TRED1(NM,N,A,W,FV1,FV2)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  REBAKB(NM,N,B,FV2,N,Z)
   50 RETURN
      END
      SUBROUTINE RSM(NM,N,A,W,M,Z,FWORK,IWORK,IERR)
C
      INTEGER N,NM,M,IWORK(N),IERR
      REAL A(NM,N),W(N),Z(NM,M),FWORK(1)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND ALL OF THE EIGENVALUES AND SOME OF THE EIGENVECTORS
C     OF A REAL SYMMETRIC MATRIX.
C
C     ON INPUT
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT.
C
C        N  IS THE ORDER OF THE MATRIX  A.
C
C        A  CONTAINS THE REAL SYMMETRIC MATRIX.
C
C        M  THE EIGENVECTORS CORRESPONDING TO THE FIRST M EIGENVALUES
C           ARE TO BE COMPUTED.
C           IF M = 0 THEN NO EIGENVECTORS ARE COMPUTED.
C           IF M = N THEN ALL OF THE EIGENVECTORS ARE COMPUTED.
C
C     ON OUTPUT
C
C        W  CONTAINS ALL N EIGENVALUES IN ASCENDING ORDER.
C
C        Z  CONTAINS THE ORTHONORMAL EIGENVECTORS ASSOCIATED WITH
C           THE FIRST M EIGENVALUES.
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT,
C           IMTQLV AND TINVIT.  THE NORMAL COMPLETION CODE IS ZERO.
C
C        FWORK  IS A TEMPORARY STORAGE ARRAY OF DIMENSION 8*N.
C
C        IWORK  IS AN INTEGER TEMPORARY STORAGE ARRAY OF DIMENSION N.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 10 * N
      IF (N .GT. NM .OR. M .GT. NM) GO TO 50
      K1 = 1
      K2 = K1 + N
      K3 = K2 + N
      K4 = K3 + N
      K5 = K4 + N
      K6 = K5 + N
      K7 = K6 + N
      K8 = K7 + N
      IF (M .GT. 0) GO TO 10
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  TRED1(NM,N,A,W,FWORK(K1),FWORK(K2))
      CALL  TQLRAT(N,W,FWORK(K2),IERR)
      GO TO 50
C     .......... FIND ALL EIGENVALUES AND M EIGENVECTORS ..........
   10 CALL  TRED1(NM,N,A,FWORK(K1),FWORK(K2),FWORK(K3))
      CALL  IMTQLV(N,FWORK(K1),FWORK(K2),FWORK(K3),W,IWORK,
     X             IERR,FWORK(K4))
      CALL  TINVIT(NM,N,FWORK(K1),FWORK(K2),FWORK(K3),M,W,IWORK,Z,IERR,
     X             FWORK(K4),FWORK(K5),FWORK(K6),FWORK(K7),FWORK(K8))
      CALL  TRBAK1(NM,N,A,FWORK(K2),M,Z)
   50 RETURN
      END
      SUBROUTINE RSP(NM,N,NV,A,W,MATZ,Z,FV1,FV2,IERR)
C
      INTEGER I,J,N,NM,NV,IERR,MATZ
      REAL A(NV),W(N),Z(NM,N),FV1(N),FV2(N)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     OF A REAL SYMMETRIC PACKED MATRIX.
C
C     ON INPUT
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT.
C
C        N  IS THE ORDER OF THE MATRIX  A.
C
C        NV  IS AN INTEGER VARIABLE SET EQUAL TO THE
C        DIMENSION OF THE ARRAY  A  AS SPECIFIED FOR
C        A  IN THE CALLING PROGRAM.  NV  MUST NOT BE
C        LESS THAN  N*(N+1)/2.
C
C        A  CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
C        PACKED MATRIX STORED ROW-WISE.
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT
C
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.
C
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT
C           AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO.
C
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (N .LE. NM) GO TO 5
      IERR = 10 * N
      GO TO 50
    5 IF (NV .GE. (N * (N + 1)) / 2) GO TO 10
      IERR = 20 * N
      GO TO 50
C
   10 CALL  TRED3(N,NV,A,W,FV1,FV2)
      IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 DO 40 I = 1, N
C
         DO 30 J = 1, N
            Z(J,I) = 0.0E0
   30    CONTINUE
C
         Z(I,I) = 1.0E0
   40 CONTINUE
C
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  TRBAK3(NM,N,NV,A,N,Z)
   50 RETURN
      END
      SUBROUTINE RST(NM,N,W,E,MATZ,Z,IERR)
C
      INTEGER I,J,N,NM,IERR,MATZ
      REAL W(N),E(N),Z(NM,N)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     OF A REAL SYMMETRIC TRIDIAGONAL MATRIX.
C
C     ON INPUT
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT.
C
C        N  IS THE ORDER OF THE MATRIX.
C
C        W  CONTAINS THE DIAGONAL ELEMENTS OF THE REAL
C        SYMMETRIC TRIDIAGONAL MATRIX.
C
C        E  CONTAINS THE SUBDIAGONAL ELEMENTS OF THE MATRIX IN
C        ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT
C
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.
C
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR IMTQL1
C           AND IMTQL2.  THE NORMAL COMPLETION CODE IS ZERO.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  IMTQL1(N,W,E,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 DO 40 I = 1, N
C
         DO 30 J = 1, N
            Z(J,I) = 0.0E0
   30    CONTINUE
C
         Z(I,I) = 1.0E0
   40 CONTINUE
C
      CALL  IMTQL2(NM,N,W,E,Z,IERR)
   50 RETURN
      END
      SUBROUTINE RT(NM,N,A,W,MATZ,Z,FV1,IERR)
C
      INTEGER N,NM,IERR,MATZ
      REAL A(NM,3),W(N),Z(NM,N),FV1(N)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     OF A SPECIAL REAL TRIDIAGONAL MATRIX.
C
C     ON INPUT
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT.
C
C        N  IS THE ORDER OF THE MATRIX  A.
C
C        A  CONTAINS THE SPECIAL REAL TRIDIAGONAL MATRIX IN ITS
C        FIRST THREE COLUMNS.  THE SUBDIAGONAL ELEMENTS ARE STORED
C        IN THE LAST  N-1  POSITIONS OF THE FIRST COLUMN, THE
C        DIAGONAL ELEMENTS IN THE SECOND COLUMN, AND THE SUPERDIAGONAL
C        ELEMENTS IN THE FIRST  N-1  POSITIONS OF THE THIRD COLUMN.
C        ELEMENTS  A(1,1)  AND  A(N,3)  ARE ARBITRARY.
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT
C
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.
C
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR IMTQL1
C           AND IMTQL2.  THE NORMAL COMPLETION CODE IS ZERO.
C
C        FV1  IS A TEMPORARY STORAGE ARRAY.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  FIGI(NM,N,A,W,FV1,FV1,IERR)
      IF (IERR .GT. 0) GO TO 50
      CALL  IMTQL1(N,W,FV1,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  FIGI2(NM,N,A,W,FV1,Z,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  IMTQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
      END
      SUBROUTINE SVD(NM,M,N,A,W,MATU,U,MATV,V,IERR,RV1)
C
      INTEGER I,J,K,L,M,N,II,I1,KK,K1,LL,L1,MN,NM,ITS,IERR
      REAL A(NM,N),W(N),U(NM,N),V(NM,N),RV1(N)
      REAL C,F,G,H,S,X,Y,Z,TST1,TST2,SCALE,PYTHAG
      LOGICAL MATU,MATV
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE SVD,
C     NUM. MATH. 14, 403-420(1970) BY GOLUB AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
C
C     THIS SUBROUTINE DETERMINES THE SINGULAR VALUE DECOMPOSITION
C          T
C     A=USV  OF A REAL M BY N RECTANGULAR MATRIX.  HOUSEHOLDER
C     BIDIAGONALIZATION AND A VARIANT OF THE QR ALGORITHM ARE USED.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.  NOTE THAT NM MUST BE AT LEAST
C          AS LARGE AS THE MAXIMUM OF M AND N.
C
C        M IS THE NUMBER OF ROWS OF A (AND U).
C
C        N IS THE NUMBER OF COLUMNS OF A (AND U) AND THE ORDER OF V.
C
C        A CONTAINS THE RECTANGULAR INPUT MATRIX TO BE DECOMPOSED.
C
C        MATU SHOULD BE SET TO .TRUE. IF THE U MATRIX IN THE
C          DECOMPOSITION IS DESIRED, AND TO .FALSE. OTHERWISE.
C
C        MATV SHOULD BE SET TO .TRUE. IF THE V MATRIX IN THE
C          DECOMPOSITION IS DESIRED, AND TO .FALSE. OTHERWISE.
C
C     ON OUTPUT
C
C        A IS UNALTERED (UNLESS OVERWRITTEN BY U OR V).
C
C        W CONTAINS THE N (NON-NEGATIVE) SINGULAR VALUES OF A (THE
C          DIAGONAL ELEMENTS OF S).  THEY ARE UNORDERED.  IF AN
C          ERROR EXIT IS MADE, THE SINGULAR VALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,IERR+2,...,N.
C
C        U CONTAINS THE MATRIX U (ORTHOGONAL COLUMN VECTORS) OF THE
C          DECOMPOSITION IF MATU HAS BEEN SET TO .TRUE.  OTHERWISE
C          U IS USED AS A TEMPORARY ARRAY.  U MAY COINCIDE WITH A.
C          IF AN ERROR EXIT IS MADE, THE COLUMNS OF U CORRESPONDING
C          TO INDICES OF CORRECT SINGULAR VALUES SHOULD BE CORRECT.
C
C        V CONTAINS THE MATRIX V (ORTHOGONAL) OF THE DECOMPOSITION IF
C          MATV HAS BEEN SET TO .TRUE.  OTHERWISE V IS NOT REFERENCED.
C          V MAY ALSO COINCIDE WITH A IF U IS NOT NEEDED.  IF AN ERROR
C          EXIT IS MADE, THE COLUMNS OF V CORRESPONDING TO INDICES OF
C          CORRECT SINGULAR VALUES SHOULD BE CORRECT.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          K          IF THE K-TH SINGULAR VALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C        RV1 IS A TEMPORARY STORAGE ARRAY.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
C
      DO 100 I = 1, M
C
         DO 100 J = 1, N
            U(I,J) = A(I,J)
  100 CONTINUE
C     .......... HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM ..........
      G = 0.0E0
      SCALE = 0.0E0
      X = 0.0E0
C
      DO 300 I = 1, N
         L = I + 1
         RV1(I) = SCALE * G
         G = 0.0E0
         S = 0.0E0
         SCALE = 0.0E0
         IF (I .GT. M) GO TO 210
C
         DO 120 K = I, M
  120    SCALE = SCALE + ABS(U(K,I))
C
         IF (SCALE .EQ. 0.0E0) GO TO 210
C
         DO 130 K = I, M
            U(K,I) = U(K,I) / SCALE
            S = S + U(K,I)**2
  130    CONTINUE
C
         F = U(I,I)
         G = -SIGN(SQRT(S),F)
         H = F * G - S
         U(I,I) = F - G
         IF (I .EQ. N) GO TO 190
C
         DO 150 J = L, N
            S = 0.0E0
C
            DO 140 K = I, M
  140       S = S + U(K,I) * U(K,J)
C
            F = S / H
C
            DO 150 K = I, M
               U(K,J) = U(K,J) + F * U(K,I)
  150    CONTINUE
C
  190    DO 200 K = I, M
  200    U(K,I) = SCALE * U(K,I)
C
  210    W(I) = SCALE * G
         G = 0.0E0
         S = 0.0E0
         SCALE = 0.0E0
         IF (I .GT. M .OR. I .EQ. N) GO TO 290
C
         DO 220 K = L, N
  220    SCALE = SCALE + ABS(U(I,K))
C
         IF (SCALE .EQ. 0.0E0) GO TO 290
C
         DO 230 K = L, N
            U(I,K) = U(I,K) / SCALE
            S = S + U(I,K)**2
  230    CONTINUE
C
         F = U(I,L)
         G = -SIGN(SQRT(S),F)
         H = F * G - S
         U(I,L) = F - G
C
         DO 240 K = L, N
  240    RV1(K) = U(I,K) / H
C
         IF (I .EQ. M) GO TO 270
C
         DO 260 J = L, M
            S = 0.0E0
C
            DO 250 K = L, N
  250       S = S + U(J,K) * U(I,K)
C
            DO 260 K = L, N
               U(J,K) = U(J,K) + S * RV1(K)
  260    CONTINUE
C
  270    DO 280 K = L, N
  280    U(I,K) = SCALE * U(I,K)
C
  290    X = AMAX1(X,ABS(W(I))+ABS(RV1(I)))
  300 CONTINUE
C     .......... ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS ..........
      IF (.NOT. MATV) GO TO 410
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 400 II = 1, N
         I = N + 1 - II
         IF (I .EQ. N) GO TO 390
         IF (G .EQ. 0.0E0) GO TO 360
C
         DO 320 J = L, N
C     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
  320    V(J,I) = (U(I,J) / U(I,L)) / G
C
         DO 350 J = L, N
            S = 0.0E0
C
            DO 340 K = L, N
  340       S = S + U(I,K) * V(K,J)
C
            DO 350 K = L, N
               V(K,J) = V(K,J) + S * V(K,I)
  350    CONTINUE
C
  360    DO 380 J = L, N
            V(I,J) = 0.0E0
            V(J,I) = 0.0E0
  380    CONTINUE
C
  390    V(I,I) = 1.0E0
         G = RV1(I)
         L = I
  400 CONTINUE
C     .......... ACCUMULATION OF LEFT-HAND TRANSFORMATIONS ..........
  410 IF (.NOT. MATU) GO TO 510
C     ..........FOR I=MIN(M,N) STEP -1 UNTIL 1 DO -- ..........
      MN = N
      IF (M .LT. N) MN = M
C
      DO 500 II = 1, MN
         I = MN + 1 - II
         L = I + 1
         G = W(I)
         IF (I .EQ. N) GO TO 430
C
         DO 420 J = L, N
  420    U(I,J) = 0.0E0
C
  430    IF (G .EQ. 0.0E0) GO TO 475
         IF (I .EQ. MN) GO TO 460
C
         DO 450 J = L, N
            S = 0.0E0
C
            DO 440 K = L, M
  440       S = S + U(K,I) * U(K,J)
C     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            F = (S / U(I,I)) / G
C
            DO 450 K = I, M
               U(K,J) = U(K,J) + F * U(K,I)
  450    CONTINUE
C
  460    DO 470 J = I, M
  470    U(J,I) = U(J,I) / G
C
         GO TO 490
C
  475    DO 480 J = I, M
  480    U(J,I) = 0.0E0
C
  490    U(I,I) = U(I,I) + 1.0E0
  500 CONTINUE
C     .......... DIAGONALIZATION OF THE BIDIAGONAL FORM ..........
  510 TST1 = X
C     .......... FOR K=N STEP -1 UNTIL 1 DO -- ..........
      DO 700 KK = 1, N
         K1 = N - KK
         K = K1 + 1
         ITS = 0
C     .......... TEST FOR SPLITTING.
C                FOR L=K STEP -1 UNTIL 1 DO -- ..........
  520    DO 530 LL = 1, K
            L1 = K - LL
            L = L1 + 1
            TST2 = TST1 + ABS(RV1(L))
            IF (TST2 .EQ. TST1) GO TO 565
C     .......... RV1(1) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
            TST2 = TST1 + ABS(W(L1))
            IF (TST2 .EQ. TST1) GO TO 540
  530    CONTINUE
C     .......... CANCELLATION OF RV1(L) IF L GREATER THAN 1 ..........
  540    C = 0.0E0
         S = 1.0E0
C
         DO 560 I = L, K
            F = S * RV1(I)
            RV1(I) = C * RV1(I)
            TST2 = TST1 + ABS(F)
            IF (TST2 .EQ. TST1) GO TO 565
            G = W(I)
            H = PYTHAG(F,G)
            W(I) = H
            C = G / H
            S = -F / H
            IF (.NOT. MATU) GO TO 560
C
            DO 550 J = 1, M
               Y = U(J,L1)
               Z = U(J,I)
               U(J,L1) = Y * C + Z * S
               U(J,I) = -Y * S + Z * C
  550       CONTINUE
C
  560    CONTINUE
C     .......... TEST FOR CONVERGENCE ..........
  565    Z = W(K)
         IF (L .EQ. K) GO TO 650
C     .......... SHIFT FROM BOTTOM 2 BY 2 MINOR ..........
         IF (ITS .EQ. 30) GO TO 1000
         ITS = ITS + 1
         X = W(L)
         Y = W(K1)
         G = RV1(K1)
         H = RV1(K)
         F = 0.5E0 * (((G + Z) / H) * ((G - Z) / Y) + Y / H - H / Y)
         G = PYTHAG(F,1.0E0)
         F = X - (Z / X) * Z + (H / X) * (Y / (F + SIGN(G,F)) - H)
C     .......... NEXT QR TRANSFORMATION ..........
         C = 1.0E0
         S = 1.0E0
C
         DO 600 I1 = L, K1
            I = I1 + 1
            G = RV1(I)
            Y = W(I)
            H = S * G
            G = C * G
            Z = PYTHAG(F,H)
            RV1(I1) = Z
            C = F / Z
            S = H / Z
            F = X * C + G * S
            G = -X * S + G * C
            H = Y * S
            Y = Y * C
            IF (.NOT. MATV) GO TO 575
C
            DO 570 J = 1, N
               X = V(J,I1)
               Z = V(J,I)
               V(J,I1) = X * C + Z * S
               V(J,I) = -X * S + Z * C
  570       CONTINUE
C
  575       Z = PYTHAG(F,H)
            W(I1) = Z
C     .......... ROTATION CAN BE ARBITRARY IF Z IS ZERO ..........
            IF (Z .EQ. 0.0E0) GO TO 580
            C = F / Z
            S = H / Z
  580       F = C * G + S * Y
            X = -S * G + C * Y
            IF (.NOT. MATU) GO TO 600
C
            DO 590 J = 1, M
               Y = U(J,I1)
               Z = U(J,I)
               U(J,I1) = Y * C + Z * S
               U(J,I) = -Y * S + Z * C
  590       CONTINUE
C
  600    CONTINUE
C
         RV1(L) = 0.0E0
         RV1(K) = F
         W(K) = X
         GO TO 520
C     .......... CONVERGENCE ..........
  650    IF (Z .GE. 0.0E0) GO TO 700
C     .......... W(K) IS MADE NON-NEGATIVE ..........
         W(K) = -Z
         IF (.NOT. MATV) GO TO 700
C
         DO 690 J = 1, N
  690    V(J,K) = -V(J,K)
C
  700 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO A
C                SINGULAR VALUE AFTER 30 ITERATIONS ..........
 1000 IERR = K
 1001 RETURN
      END
      SUBROUTINE TINVIT(NM,N,D,E,E2,M,W,IND,Z,
     X                  IERR,RV1,RV2,RV3,RV4,RV6)
C
      INTEGER I,J,M,N,P,Q,R,S,II,IP,JJ,NM,ITS,TAG,IERR,GROUP
      REAL D(N),E(N),E2(N),W(M),Z(NM,M),
     X       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
      REAL U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER,EPSLON,
     X       PYTHAG
      INTEGER IND(M)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE INVERSE ITERATION TECH-
C     NIQUE IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL
C     SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
C     USING INVERSE ITERATION.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E,
C          WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.
C          E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
C          THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE SUM
C          OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST CONTAIN
C          0.0E0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0E0
C          IF THE EIGENVALUES ARE IN DESCENDING ORDER.  IF  BISECT,
C          TRIDIB, OR  IMTQLV  HAS BEEN USED TO FIND THE EIGENVALUES,
C          THEIR OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE.
C
C        M IS THE NUMBER OF SPECIFIED EIGENVALUES.
C
C        W CONTAINS THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER.
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.
C
C     ON OUTPUT
C
C        ALL INPUT ARRAYS ARE UNALTERED.
C
C        Z CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS.
C          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          -R         IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
C                     EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS.
C
C        RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (M .EQ. 0) GO TO 1001
      TAG = 0
      ORDER = 1.0E0 - E2(1)
      Q = 0
C     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX ..........
  100 P = Q + 1
C
      DO 120 Q = P, N
         IF (Q .EQ. N) GO TO 140
         IF (E2(Q+1) .EQ. 0.0E0) GO TO 140
  120 CONTINUE
C     .......... FIND VECTORS BY INVERSE ITERATION ..........
  140 TAG = TAG + 1
      S = 0
C
      DO 920 R = 1, M
         IF (IND(R) .NE. TAG) GO TO 920
         ITS = 1
         X1 = W(R)
         IF (S .NE. 0) GO TO 510
C     .......... CHECK FOR ISOLATED ROOT ..........
         XU = 1.0E0
         IF (P .NE. Q) GO TO 490
         RV6(P) = 1.0E0
         GO TO 870
  490    NORM = ABS(D(P))
         IP = P + 1
C
         DO 500 I = IP, Q
  500    NORM = AMAX1(NORM, ABS(D(I))+ABS(E(I)))
C     .......... EPS2 IS THE CRITERION FOR GROUPING,
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                ROOTS ARE MODIFIED BY EPS3,
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ..........
         EPS2 = 1.0E-3 * NORM
         EPS3 = EPSLON(NORM)
         UK = Q - P + 1
         EPS4 = UK * EPS3
         UK = EPS4 / SQRT(UK)
         S = P
  505    GROUP = 0
         GO TO 520
C     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
  510    IF (ABS(X1-X0) .GE. EPS2) GO TO 505
         GROUP = GROUP + 1
         IF (ORDER * (X1 - X0) .LE. 0.0E0) X1 = X0 + ORDER * EPS3
C     .......... ELIMINATION WITH INTERCHANGES AND
C                INITIALIZATION OF VECTOR ..........
  520    V = 0.0E0
C
         DO 580 I = P, Q
            RV6(I) = UK
            IF (I .EQ. P) GO TO 560
            IF (ABS(E(I)) .LT. ABS(U)) GO TO 540
C     .......... WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
C                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY ..........
            XU = U / E(I)
            RV4(I) = XU
            RV1(I-1) = E(I)
            RV2(I-1) = D(I) - X1
            RV3(I-1) = 0.0E0
            IF (I .NE. Q) RV3(I-1) = E(I+1)
            U = V - XU * RV2(I-1)
            V = -XU * RV3(I-1)
            GO TO 580
  540       XU = E(I) / U
            RV4(I) = XU
            RV1(I-1) = U
            RV2(I-1) = V
            RV3(I-1) = 0.0E0
  560       U = D(I) - X1 - XU * V
            IF (I .NE. Q) V = E(I+1)
  580    CONTINUE
C
         IF (U .EQ. 0.0E0) U = EPS3
         RV1(Q) = U
         RV2(Q) = 0.0E0
         RV3(Q) = 0.0E0
C     .......... BACK SUBSTITUTION
C                FOR I=Q STEP -1 UNTIL P DO -- ..........
  600    DO 620 II = P, Q
            I = P + Q - II
            RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
            V = U
            U = RV6(I)
  620    CONTINUE
C     .......... ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                MEMBERS OF GROUP ..........
         IF (GROUP .EQ. 0) GO TO 700
         J = R
C
         DO 680 JJ = 1, GROUP
  630       J = J - 1
            IF (IND(J) .NE. TAG) GO TO 630
            XU = 0.0E0
C
            DO 640 I = P, Q
  640       XU = XU + RV6(I) * Z(I,J)
C
            DO 660 I = P, Q
  660       RV6(I) = RV6(I) - XU * Z(I,J)
C
  680    CONTINUE
C
  700    NORM = 0.0E0
C
         DO 720 I = P, Q
  720    NORM = NORM + ABS(RV6(I))
C
         IF (NORM .GE. 1.0E0) GO TO 840
C     .......... FORWARD SUBSTITUTION ..........
         IF (ITS .EQ. 5) GO TO 830
         IF (NORM .NE. 0.0E0) GO TO 740
         RV6(S) = EPS4
         S = S + 1
         IF (S .GT. Q) S = P
         GO TO 780
  740    XU = EPS4 / NORM
C
         DO 760 I = P, Q
  760    RV6(I) = RV6(I) * XU
C     .......... ELIMINATION OPERATIONS ON NEXT VECTOR
C                ITERATE ..........
  780    DO 820 I = IP, Q
            U = RV6(I)
C     .......... IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
C                WAS PERFORMED EARLIER IN THE
C                TRIANGULARIZATION PROCESS ..........
            IF (RV1(I-1) .NE. E(I)) GO TO 800
            U = RV6(I-1)
            RV6(I-1) = RV6(I)
  800       RV6(I) = U - RV4(I) * RV6(I-1)
  820    CONTINUE
C
         ITS = ITS + 1
         GO TO 600
C     .......... SET ERROR -- NON-CONVERGED EIGENVECTOR ..........
  830    IERR = -R
         XU = 0.0E0
         GO TO 870
C     .......... NORMALIZE SO THAT SUM OF SQUARES IS
C                1 AND EXPAND TO FULL ORDER ..........
  840    U = 0.0E0
C
         DO 860 I = P, Q
  860    U = PYTHAG(U,RV6(I))
C
         XU = 1.0E0 / U
C
  870    DO 880 I = 1, N
  880    Z(I,R) = 0.0E0
C
         DO 900 I = P, Q
  900    Z(I,R) = RV6(I) * XU
C
         X0 = X1
  920 CONTINUE
C
      IF (Q .LT. N) GO TO 100
 1001 RETURN
      END
      SUBROUTINE TQL1(N,D,E,IERR)
C
      INTEGER I,J,L,M,N,II,L1,L2,MML,IERR
      REAL D(N),E(N)
      REAL C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2,PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL1,
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
C     WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE QL METHOD.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES.
C
C        E HAS BEEN DESTROYED.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0E0
      TST1 = 0.0E0
      E(N) = 0.0E0
C
      DO 290 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (TST1 .LT. H) TST1 = H
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            TST2 = TST1 + ABS(E(M))
            IF (TST2 .EQ. TST1) GO TO 120
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0E0 * E(L))
         R = PYTHAG(P,1.0E0)
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
C
         DO 140 I = L2, N
  140    D(I) = D(I) - H
C
  145    F = F + H
C     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0E0
         C2 = C
         EL1 = E(L1)
         S = 0.0E0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            R = PYTHAG(P,E(I))
            E(I+1) = S * R
            S = E(I) / R
            C = P / R
            P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
  200    CONTINUE
C
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         TST2 = TST1 + ABS(E(L))
         IF (TST2 .GT. TST1) GO TO 130
  210    P = D(L) + F
C     .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
C
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      REAL D(N),E(N),Z(NM,N)
      REAL C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2,PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
C     WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1.
C
C        E HAS BEEN DESTROYED.
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0E0
      TST1 = 0.0E0
      E(N) = 0.0E0
C
      DO 240 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (TST1 .LT. H) TST1 = H
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            TST2 = TST1 + ABS(E(M))
            IF (TST2 .EQ. TST1) GO TO 120
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0E0 * E(L))
         R = PYTHAG(P,1.0E0)
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
C
         DO 140 I = L2, N
  140    D(I) = D(I) - H
C
  145    F = F + H
C     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0E0
         C2 = C
         EL1 = E(L1)
         S = 0.0E0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            R = PYTHAG(P,E(I))
            E(I+1) = S * R
            S = E(I) / R
            C = P / R
            P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
C
  200    CONTINUE
C
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         TST2 = TST1 + ABS(E(L))
         IF (TST2 .GT. TST1) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      SUBROUTINE TQLRAT(N,D,E2,IERR)
C
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      REAL D(N),E2(N)
      REAL B,C,F,G,H,P,R,S,T,EPSLON,PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES.
C
C        E2 HAS BEEN DESTROYED.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
C
      F = 0.0E0
      T = 0.0E0
      E2(N) = 0.0E0
C
      DO 290 L = 1, N
         J = 0
         H = ABS(D(L)) + SQRT(E2(L))
         IF (T .GT. H) GO TO 105
         T = H
         B = EPSLON(T)
         C = B * B
C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         S = SQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0E0 * S)
         R = PYTHAG(P,1.0E0)
         D(L) = S / (P + SIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
C     .......... RATIONAL QL TRANSFORMATION ..........
         G = D(M)
         IF (G .EQ. 0.0E0) G = B
         H = G
         S = 0.0E0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G .EQ. 0.0E0) G = B
            H = G * P / R
  200    CONTINUE
C
         E2(L) = S * G
         D(L) = H
C     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
         IF (H .EQ. 0.0E0) GO TO 210
         IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0E0) GO TO 130
  210    P = D(L) + F
C     .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      SUBROUTINE TRBAK1(NM,N,A,E,M,Z)
C
      INTEGER I,J,K,L,M,N,NM
      REAL A(NM,N),E(N),Z(NM,M)
      REAL S
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK1,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED1.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  TRED1
C          IN ITS STRICT LOWER TRIANGLE.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
C
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS
C          IN ITS FIRST M COLUMNS.
C
C     NOTE THAT TRBAK1 PRESERVES VECTOR EUCLIDEAN NORMS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
      IF (N .EQ. 1) GO TO 200
C
      DO 140 I = 2, N
         L = I - 1
         IF (E(I) .EQ. 0.0E0) GO TO 140
C
         DO 130 J = 1, M
            S = 0.0E0
C
            DO 110 K = 1, L
  110       S = S + A(I,K) * Z(K,J)
C     .......... DIVISOR BELOW IS NEGATIVE OF H FORMED IN TRED1.
C                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            S = (S / A(I,L)) / E(I)
C
            DO 120 K = 1, L
  120       Z(K,J) = Z(K,J) + S * A(I,K)
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE TRBAK3(NM,N,NV,A,M,Z)
C
      INTEGER I,J,K,L,M,N,IK,IZ,NM,NV
      REAL A(NV),Z(NM,M)
      REAL H,S
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED3.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT.
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANSFORMATIONS
C          USED IN THE REDUCTION BY  TRED3  IN ITS FIRST
C          N*(N+1)/2 POSITIONS.
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
C
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS
C          IN ITS FIRST M COLUMNS.
C
C     NOTE THAT TRBAK3 PRESERVES VECTOR EUCLIDEAN NORMS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
      IF (N .EQ. 1) GO TO 200
C
      DO 140 I = 2, N
         L = I - 1
         IZ = (I * L) / 2
         IK = IZ + I
         H = A(IK)
         IF (H .EQ. 0.0E0) GO TO 140
C
         DO 130 J = 1, M
            S = 0.0E0
            IK = IZ
C
            DO 110 K = 1, L
               IK = IK + 1
               S = S + A(IK) * Z(K,J)
  110       CONTINUE
C     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            S = (S / H) / H
            IK = IZ
C
            DO 120 K = 1, L
               IK = IK + 1
               Z(K,J) = Z(K,J) - S * A(IK)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE TRED1(NM,N,A,D,E,E2)
C
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL A(NM,N),D(N),E(N),E2(N)
      REAL F,G,H,SCALE
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX
C     TO A SYMMETRIC TRIDIAGONAL MATRIX USING
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
C          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER
C          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      DO 100 I = 1, N
         D(I) = A(N,I)
         A(N,I) = A(I,I)
  100 CONTINUE
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0E0
         SCALE = 0.0E0
         IF (L .LT. 1) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(D(K))
C
         IF (SCALE .NE. 0.0E0) GO TO 140
C
         DO 125 J = 1, L
            D(J) = A(L,J)
            A(L,J) = A(I,J)
            A(I,J) = 0.0E0
  125    CONTINUE
C
  130    E(I) = 0.0E0
         E2(I) = 0.0E0
         GO TO 300
C
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         F = D(L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         IF (L .EQ. 1) GO TO 285
C     .......... FORM A*U ..........
         DO 170 J = 1, L
  170    E(J) = 0.0E0
C
         DO 240 J = 1, L
            F = D(J)
            G = E(J) + A(J,J) * F
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
               G = G + A(K,J) * D(K)
               E(K) = E(K) + A(K,J) * F
  200       CONTINUE
C
  220       E(J) = G
  240    CONTINUE
C     .......... FORM P ..........
         F = 0.0E0
C
         DO 245 J = 1, L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
  245    CONTINUE
C
         H = F / (H + H)
C     .......... FORM Q ..........
         DO 250 J = 1, L
  250    E(J) = E(J) - H * D(J)
C     .......... FORM REDUCED A ..........
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
C
            DO 260 K = J, L
  260       A(K,J) = A(K,J) - F * E(K) - G * D(K)
C
  280    CONTINUE
C
  285    DO 290 J = 1, L
            F = D(J)
            D(J) = A(L,J)
            A(L,J) = A(I,J)
            A(I,J) = F * SCALE
  290    CONTINUE
C
  300 CONTINUE
C
      RETURN
      END
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
C
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL A(NM,N),D(N),E(N),Z(NM,N)
      REAL F,G,H,HH,SCALE
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A
C     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
C
C        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX
C          PRODUCED IN THE REDUCTION.
C
C        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      DO 100 I = 1, N
C
         DO 80 J = I, N
   80    Z(J,I) = A(J,I)
C
         D(I) = A(N,I)
  100 CONTINUE
C
      IF (N .EQ. 1) GO TO 510
C     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0E0
         SCALE = 0.0E0
         IF (L .LT. 2) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(D(K))
C
         IF (SCALE .NE. 0.0E0) GO TO 140
  130    E(I) = D(L)
C
         DO 135 J = 1, L
            D(J) = Z(L,J)
            Z(I,J) = 0.0E0
            Z(J,I) = 0.0E0
  135    CONTINUE
C
         GO TO 290
C
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
C
         F = D(L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
C     .......... FORM A*U ..........
         DO 170 J = 1, L
  170    E(J) = 0.0E0
C
         DO 240 J = 1, L
            F = D(J)
            Z(J,I) = F
            G = E(J) + Z(J,J) * F
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
               G = G + Z(K,J) * D(K)
               E(K) = E(K) + Z(K,J) * F
  200       CONTINUE
C
  220       E(J) = G
  240    CONTINUE
C     .......... FORM P ..........
         F = 0.0E0
C
         DO 245 J = 1, L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
  245    CONTINUE
C
         HH = F / (H + H)
C     .......... FORM Q ..........
         DO 250 J = 1, L
  250    E(J) = E(J) - HH * D(J)
C     .......... FORM REDUCED A ..........
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
C
            DO 260 K = J, L
  260       Z(K,J) = Z(K,J) - F * E(K) - G * D(K)
C
            D(J) = Z(L,J)
            Z(I,J) = 0.0E0
  280    CONTINUE
C
  290    D(I) = H
  300 CONTINUE
C     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
      DO 500 I = 2, N
         L = I - 1
         Z(N,L) = Z(L,L)
         Z(L,L) = 1.0E0
         H = D(I)
         IF (H .EQ. 0.0E0) GO TO 380
C
         DO 330 K = 1, L
  330    D(K) = Z(K,I) / H
C
         DO 360 J = 1, L
            G = 0.0E0
C
            DO 340 K = 1, L
  340       G = G + Z(K,I) * Z(K,J)
C
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * D(K)
  360    CONTINUE
C
  380    DO 400 K = 1, L
  400    Z(K,I) = 0.0E0
C
  500 CONTINUE
C
  510 DO 520 I = 1, N
         D(I) = Z(N,I)
         Z(N,I) = 0.0E0
  520 CONTINUE
C
      Z(N,N) = 1.0E0
      E(1) = 0.0E0
      RETURN
      END
      SUBROUTINE TRED3(N,NV,A,D,E,E2)
C
      INTEGER I,J,K,L,N,II,IZ,JK,NV,JM1
      REAL A(NV),D(N),E(N),E2(N)
      REAL F,G,H,HH,SCALE
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX, STORED AS
C     A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX
C     USING ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT.
C
C        A CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
C          INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL
C          ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.
C
C     ON OUTPUT
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL
C          TRANSFORMATIONS USED IN THE REDUCTION.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         IZ = (I * L) / 2
         H = 0.0E0
         SCALE = 0.0E0
         IF (L .LT. 1) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
            IZ = IZ + 1
            D(K) = A(IZ)
            SCALE = SCALE + ABS(D(K))
  120    CONTINUE
C
         IF (SCALE .NE. 0.0E0) GO TO 140
  130    E(I) = 0.0E0
         E2(I) = 0.0E0
         GO TO 290
C
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         F = D(L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         A(IZ) = SCALE * D(L)
         IF (L .EQ. 1) GO TO 290
         JK = 1
C
         DO 240 J = 1, L
            F = D(J)
            G = 0.0E0
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 220
C
            DO 200 K = 1, JM1
               G = G + A(JK) * D(K)
               E(K) = E(K) + A(JK) * F
               JK = JK + 1
  200       CONTINUE
C
  220       E(J) = G + A(JK) * F
            JK = JK + 1
  240    CONTINUE
C     .......... FORM P ..........
         F = 0.0E0
C
         DO 245 J = 1, L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
  245    CONTINUE
C
         HH = F / (H + H)
C     .......... FORM Q ..........
         DO 250 J = 1, L
  250    E(J) = E(J) - HH * D(J)
C
         JK = 1
C     .......... FORM REDUCED A ..........
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
C
            DO 260 K = 1, J
               A(JK) = A(JK) - F * E(K) - G * D(K)
               JK = JK + 1
  260       CONTINUE
C
  280    CONTINUE
C
  290    D(I) = A(IZ+1)
         A(IZ+1) = SCALE * SQRT(H)
  300 CONTINUE
C
      RETURN
      END
      SUBROUTINE TRIDIB(N,EPS1,D,E,E2,LB,UB,M11,M,W,IND,IERR,RV4,RV5)
C
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,M1,M2,M11,M22,TAG,IERR,ISTURM
      REAL D(N),E(N),E2(N),W(M),RV4(N),RV5(N)
      REAL U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,TST1,TST2,EPSLON
      INTEGER IND(M)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BISECT,
C     NUM. MATH. 9, 386-393(1967) BY BARTH, MARTIN, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 249-256(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVALUES OF A TRIDIAGONAL
C     SYMMETRIC MATRIX BETWEEN SPECIFIED BOUNDARY INDICES,
C     USING BISECTION.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        EPS1 IS AN ABSOLUTE ERROR TOLERANCE FOR THE COMPUTED
C          EIGENVALUES.  IF THE INPUT EPS1 IS NON-POSITIVE,
C          IT IS RESET FOR EACH SUBMATRIX TO A DEFAULT VALUE,
C          NAMELY, MINUS THE PRODUCT OF THE RELATIVE MACHINE
C          PRECISION AND THE 1-NORM OF THE SUBMATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY.
C
C        M11 SPECIFIES THE LOWER BOUNDARY INDEX FOR THE DESIRED
C          EIGENVALUES.
C
C        M SPECIFIES THE NUMBER OF EIGENVALUES DESIRED.  THE UPPER
C          BOUNDARY INDEX M22 IS THEN OBTAINED AS M22=M11+M-1.
C
C     ON OUTPUT
C
C        EPS1 IS UNALTERED UNLESS IT HAS BEEN RESET TO ITS
C          (LAST) DEFAULT VALUE.
C
C        D AND E ARE UNALTERED.
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS ALSO SET TO ZERO.
C
C        LB AND UB DEFINE AN INTERVAL CONTAINING EXACTLY THE DESIRED
C          EIGENVALUES.
C
C        W CONTAINS, IN ITS FIRST M POSITIONS, THE EIGENVALUES
C          BETWEEN INDICES M11 AND M22 IN ASCENDING ORDER.
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC..
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          3*N+1      IF MULTIPLE EIGENVALUES AT INDEX M11 MAKE
C                     UNIQUE SELECTION IMPOSSIBLE,
C          3*N+2      IF MULTIPLE EIGENVALUES AT INDEX M22 MAKE
C                     UNIQUE SELECTION IMPOSSIBLE.
C
C        RV4 AND RV5 ARE TEMPORARY STORAGE ARRAYS.
C
C     NOTE THAT SUBROUTINE TQL1, IMTQL1, OR TQLRAT IS GENERALLY FASTER
C     THAN TRIDIB, IF MORE THAN N/4 EIGENVALUES ARE TO BE FOUND.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      TAG = 0
      XU = D(1)
      X0 = D(1)
      U = 0.0E0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DETERMINE AN
C                INTERVAL CONTAINING ALL THE EIGENVALUES ..........
      DO 40 I = 1, N
         X1 = U
         U = 0.0E0
         IF (I .NE. N) U = ABS(E(I+1))
         XU = AMIN1(D(I)-(X1+U),XU)
         X0 = AMAX1(D(I)+(X1+U),X0)
         IF (I .EQ. 1) GO TO 20
         TST1 = ABS(D(I)) + ABS(D(I-1))
         TST2 = TST1 + ABS(E(I))
         IF (TST2 .GT. TST1) GO TO 40
   20    E2(I) = 0.0E0
   40 CONTINUE
C
      X1 = N
      X1 = X1 * EPSLON(AMAX1(ABS(XU),ABS(X0)))
      XU = XU - X1
      T1 = XU
      X0 = X0 + X1
      T2 = X0
C     .......... DETERMINE AN INTERVAL CONTAINING EXACTLY
C                THE DESIRED EIGENVALUES ..........
      P = 1
      Q = N
      M1 = M11 - 1
      IF (M1 .EQ. 0) GO TO 75
      ISTURM = 1
   50 V = X1
      X1 = XU + (X0 - XU) * 0.5E0
      IF (X1 .EQ. V) GO TO 980
      GO TO 320
   60 IF (S - M1) 65, 73, 70
   65 XU = X1
      GO TO 50
   70 X0 = X1
      GO TO 50
   73 XU = X1
      T1 = X1
   75 M22 = M1 + M
      IF (M22 .EQ. N) GO TO 90
      X0 = T2
      ISTURM = 2
      GO TO 50
   80 IF (S - M22) 65, 85, 70
   85 T2 = X1
   90 Q = 0
      R = 0
C     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
C                INTERVAL BY THE GERSCHGORIN BOUNDS ..........
  100 IF (R .EQ. M) GO TO 1001
      TAG = TAG + 1
      P = Q + 1
      XU = D(P)
      X0 = D(P)
      U = 0.0E0
C
      DO 120 Q = P, N
         X1 = U
         U = 0.0E0
         V = 0.0E0
         IF (Q .EQ. N) GO TO 110
         U = ABS(E(Q+1))
         V = E2(Q+1)
  110    XU = AMIN1(D(Q)-(X1+U),XU)
         X0 = AMAX1(D(Q)+(X1+U),X0)
         IF (V .EQ. 0.0E0) GO TO 140
  120 CONTINUE
C
  140 X1 = EPSLON(AMAX1(ABS(XU),ABS(X0)))
      IF (EPS1 .LE. 0.0E0) EPS1 = -X1
      IF (P .NE. Q) GO TO 180
C     .......... CHECK FOR ISOLATED ROOT WITHIN INTERVAL ..........
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 900
  180 X1 = X1 * (Q - P + 1)
      LB = AMAX1(T1,XU-X1)
      UB = AMIN1(T2,X0+X1)
      X1 = LB
      ISTURM = 3
      GO TO 320
  200 M1 = S + 1
      X1 = UB
      ISTURM = 4
      GO TO 320
  220 M2 = S
      IF (M1 .GT. M2) GO TO 940
C     .......... FIND ROOTS BY BISECTION ..........
      X0 = UB
      ISTURM = 5
C
      DO 240 I = M1, M2
         RV5(I) = UB
         RV4(I) = LB
  240 CONTINUE
C     .......... LOOP FOR K-TH EIGENVALUE
C                FOR K=M2 STEP -1 UNTIL M1 DO --
C                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) ..........
      K = M2
  250    XU = LB
C     .......... FOR I=K STEP -1 UNTIL M1 DO -- ..........
         DO 260 II = M1, K
            I = M1 + K - II
            IF (XU .GE. RV4(I)) GO TO 260
            XU = RV4(I)
            GO TO 280
  260    CONTINUE
C
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K)
C     .......... NEXT BISECTION STEP ..........
  300    X1 = (XU + X0) * 0.5E0
         IF ((X0 - XU) .LE. ABS(EPS1)) GO TO 420
         TST1 = 2.0E0 * (ABS(XU) + ABS(X0))
         TST2 = TST1 + (X0 - XU)
         IF (TST2 .EQ. TST1) GO TO 420
C     .......... IN-LINE PROCEDURE FOR STURM SEQUENCE ..........
  320    S = P - 1
         U = 1.0E0
C
         DO 340 I = P, Q
            IF (U .NE. 0.0E0) GO TO 325
            V = ABS(E(I)) / EPSLON(1.0E0)
            IF (E2(I) .EQ. 0.0E0) V = 0.0E0
            GO TO 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
            IF (U .LT. 0.0E0) S = S + 1
  340    CONTINUE
C
         GO TO (60,80,200,220,360), ISTURM
C     .......... REFINE INTERVALS ..........
  360    IF (S .GE. K) GO TO 400
         XU = X1
         IF (S .GE. M1) GO TO 380
         RV4(M1) = X1
         GO TO 300
  380    RV4(S+1) = X1
         IF (RV5(S) .GT. X1) RV5(S) = X1
         GO TO 300
  400    X0 = X1
         GO TO 300
C     .......... K-TH EIGENVALUE FOUND ..........
  420    RV5(K) = X1
      K = K - 1
      IF (K .GE. M1) GO TO 250
C     .......... ORDER EIGENVALUES TAGGED WITH THEIR
C                SUBMATRIX ASSOCIATIONS ..........
  900 S = R
      R = R + M2 - M1 + 1
      J = 1
      K = M1
C
      DO 920 L = 1, R
         IF (J .GT. S) GO TO 910
         IF (K .GT. M2) GO TO 940
         IF (RV5(K) .GE. W(L)) GO TO 915
C
         DO 905 II = J, S
            I = L + S - II
            W(I+1) = W(I)
            IND(I+1) = IND(I)
  905    CONTINUE
C
  910    W(L) = RV5(K)
         IND(L) = TAG
         K = K + 1
         GO TO 920
  915    J = J + 1
  920 CONTINUE
C
  940 IF (Q .LT. N) GO TO 100
      GO TO 1001
C     .......... SET ERROR -- INTERVAL CANNOT BE FOUND CONTAINING
C                EXACTLY THE DESIRED EIGENVALUES ..........
  980 IERR = 3 * N + ISTURM
 1001 LB = T1
      UB = T2
      RETURN
      END
      SUBROUTINE TSTURM(NM,N,EPS1,D,E,E2,LB,UB,MM,M,W,Z,
     X                  IERR,RV1,RV2,RV3,RV4,RV5,RV6)
C
      INTEGER I,J,K,M,N,P,Q,R,S,II,IP,JJ,MM,M1,M2,NM,ITS,
     X        IERR,GROUP,ISTURM
      REAL D(N),E(N),E2(N),W(MM),Z(NM,MM),
     X       RV1(N),RV2(N),RV3(N),RV4(N),RV5(N),RV6(N)
      REAL U,V,LB,T1,T2,UB,UK,XU,X0,X1,EPS1,EPS2,EPS3,EPS4,
     X       NORM,TST1,TST2,EPSLON,PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRISTURM
C     BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVALUES OF A TRIDIAGONAL
C     SYMMETRIC MATRIX WHICH LIE IN A SPECIFIED INTERVAL AND THEIR
C     ASSOCIATED EIGENVECTORS, USING BISECTION AND INVERSE ITERATION.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        EPS1 IS AN ABSOLUTE ERROR TOLERANCE FOR THE COMPUTED
C          EIGENVALUES.  IT SHOULD BE CHOSEN COMMENSURATE WITH
C          RELATIVE PERTURBATIONS IN THE MATRIX ELEMENTS OF THE
C          ORDER OF THE RELATIVE MACHINE PRECISION.  IF THE
C          INPUT EPS1 IS NON-POSITIVE, IT IS RESET FOR EACH
C          SUBMATRIX TO A DEFAULT VALUE, NAMELY, MINUS THE
C          PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE
C          1-NORM OF THE SUBMATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY.
C
C        LB AND UB DEFINE THE INTERVAL TO BE SEARCHED FOR EIGENVALUES.
C          IF LB IS NOT LESS THAN UB, NO EIGENVALUES WILL BE FOUND.
C
C        MM SHOULD BE SET TO AN UPPER BOUND FOR THE NUMBER OF
C          EIGENVALUES IN THE INTERVAL.  WARNING. IF MORE THAN
C          MM EIGENVALUES ARE DETERMINED TO LIE IN THE INTERVAL,
C          AN ERROR RETURN IS MADE WITH NO VALUES OR VECTORS FOUND.
C
C     ON OUTPUT
C
C        EPS1 IS UNALTERED UNLESS IT HAS BEEN RESET TO ITS
C          (LAST) DEFAULT VALUE.
C
C        D AND E ARE UNALTERED.
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS ALSO SET TO ZERO.
C
C        M IS THE NUMBER OF EIGENVALUES DETERMINED TO LIE IN (LB,UB).
C
C        W CONTAINS THE M EIGENVALUES IN ASCENDING ORDER IF THE MATRIX
C          DOES NOT SPLIT.  IF THE MATRIX SPLITS, THE EIGENVALUES ARE
C          IN ASCENDING ORDER FOR EACH SUBMATRIX.  IF A VECTOR ERROR
C          EXIT IS MADE, W CONTAINS THOSE VALUES ALREADY FOUND.
C
C        Z CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS.
C          IF AN ERROR EXIT IS MADE, Z CONTAINS THOSE VECTORS
C          ALREADY FOUND.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          3*N+1      IF M EXCEEDS MM.
C          4*N+R      IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
C                     EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS.
C
C        RV1, RV2, RV3, RV4, RV5, AND RV6 ARE TEMPORARY STORAGE ARRAYS.
C
C     THE ALGOL PROCEDURE STURMCNT CONTAINED IN TRISTURM
C     APPEARS IN TSTURM IN-LINE.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      T1 = LB
      T2 = UB
C     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES ..........
      DO 40 I = 1, N
         IF (I .EQ. 1) GO TO 20
         TST1 = ABS(D(I)) + ABS(D(I-1))
         TST2 = TST1 + ABS(E(I))
         IF (TST2 .GT. TST1) GO TO 40
   20    E2(I) = 0.0E0
   40 CONTINUE
C     .......... DETERMINE THE NUMBER OF EIGENVALUES
C                IN THE INTERVAL ..........
      P = 1
      Q = N
      X1 = UB
      ISTURM = 1
      GO TO 320
   60 M = S
      X1 = LB
      ISTURM = 2
      GO TO 320
   80 M = M - S
      IF (M .GT. MM) GO TO 980
      Q = 0
      R = 0
C     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
C                INTERVAL BY THE GERSCHGORIN BOUNDS ..........
  100 IF (R .EQ. M) GO TO 1001
      P = Q + 1
      XU = D(P)
      X0 = D(P)
      U = 0.0E0
C
      DO 120 Q = P, N
         X1 = U
         U = 0.0E0
         V = 0.0E0
         IF (Q .EQ. N) GO TO 110
         U = ABS(E(Q+1))
         V = E2(Q+1)
  110    XU = AMIN1(D(Q)-(X1+U),XU)
         X0 = AMAX1(D(Q)+(X1+U),X0)
         IF (V .EQ. 0.0E0) GO TO 140
  120 CONTINUE
C
  140 X1 = EPSLON(AMAX1(ABS(XU),ABS(X0)))
      IF (EPS1 .LE. 0.0E0) EPS1 = -X1
      IF (P .NE. Q) GO TO 180
C     .......... CHECK FOR ISOLATED ROOT WITHIN INTERVAL ..........
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940
      R = R + 1
C
      DO 160 I = 1, N
  160 Z(I,R) = 0.0E0
C
      W(R) = D(P)
      Z(P,R) = 1.0E0
      GO TO 940
  180 U = Q-P+1
      X1 = U * X1
      LB = AMAX1(T1,XU-X1)
      UB = AMIN1(T2,X0+X1)
      X1 = LB
      ISTURM = 3
      GO TO 320
  200 M1 = S + 1
      X1 = UB
      ISTURM = 4
      GO TO 320
  220 M2 = S
      IF (M1 .GT. M2) GO TO 940
C     .......... FIND ROOTS BY BISECTION ..........
      X0 = UB
      ISTURM = 5
C
      DO 240 I = M1, M2
         RV5(I) = UB
         RV4(I) = LB
  240 CONTINUE
C     .......... LOOP FOR K-TH EIGENVALUE
C                FOR K=M2 STEP -1 UNTIL M1 DO --
C                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) ..........
      K = M2
  250    XU = LB
C     .......... FOR I=K STEP -1 UNTIL M1 DO -- ..........
         DO 260 II = M1, K
            I = M1 + K - II
            IF (XU .GE. RV4(I)) GO TO 260
            XU = RV4(I)
            GO TO 280
  260    CONTINUE
C
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K)
C     .......... NEXT BISECTION STEP ..........
  300    X1 = (XU + X0) * 0.5E0
         IF ((X0 - XU) .LE. ABS(EPS1)) GO TO 420
         TST1 = 2.0E0 * (ABS(XU) + ABS(X0))
         TST2 = TST1 + (X0 - XU)
         IF (TST2 .EQ. TST1) GO TO 420
C     .......... IN-LINE PROCEDURE FOR STURM SEQUENCE ..........
  320    S = P - 1
         U = 1.0E0
C
         DO 340 I = P, Q
            IF (U .NE. 0.0E0) GO TO 325
            V = ABS(E(I)) / EPSLON(1.0E0)
            IF (E2(I) .EQ. 0.0E0) V = 0.0E0
            GO TO 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
            IF (U .LT. 0.0E0) S = S + 1
  340    CONTINUE
C
         GO TO (60,80,200,220,360), ISTURM
C     .......... REFINE INTERVALS ..........
  360    IF (S .GE. K) GO TO 400
         XU = X1
         IF (S .GE. M1) GO TO 380
         RV4(M1) = X1
         GO TO 300
  380    RV4(S+1) = X1
         IF (RV5(S) .GT. X1) RV5(S) = X1
         GO TO 300
  400    X0 = X1
         GO TO 300
C     .......... K-TH EIGENVALUE FOUND ..........
  420    RV5(K) = X1
      K = K - 1
      IF (K .GE. M1) GO TO 250
C     .......... FIND VECTORS BY INVERSE ITERATION ..........
      NORM = ABS(D(P))
      IP = P + 1
C
      DO 500 I = IP, Q
  500 NORM = AMAX1(NORM, ABS(D(I)) + ABS(E(I)))
C     .......... EPS2 IS THE CRITERION FOR GROUPING,
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                ROOTS ARE MODIFIED BY EPS3,
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ..........
      EPS2 = 1.0E-3 * NORM
      EPS3 = EPSLON(NORM)
      UK = Q - P + 1
      EPS4 = UK * EPS3
      UK = EPS4 / SQRT(UK)
      GROUP = 0
      S = P
C
      DO 920 K = M1, M2
         R = R + 1
         ITS = 1
         W(R) = RV5(K)
         X1 = RV5(K)
C     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
         IF (K .EQ. M1) GO TO 520
         IF (X1 - X0 .GE. EPS2) GROUP = -1
         GROUP = GROUP + 1
         IF (X1 .LE. X0) X1 = X0 + EPS3
C     .......... ELIMINATION WITH INTERCHANGES AND
C                INITIALIZATION OF VECTOR ..........
  520    V = 0.0E0
C
         DO 580 I = P, Q
            RV6(I) = UK
            IF (I .EQ. P) GO TO 560
            IF (ABS(E(I)) .LT. ABS(U)) GO TO 540
            XU = U / E(I)
            RV4(I) = XU
            RV1(I-1) = E(I)
            RV2(I-1) = D(I) - X1
            RV3(I-1) = 0.0E0
            IF (I .NE. Q) RV3(I-1) = E(I+1)
            U = V - XU * RV2(I-1)
            V = -XU * RV3(I-1)
            GO TO 580
  540       XU = E(I) / U
            RV4(I) = XU
            RV1(I-1) = U
            RV2(I-1) = V
            RV3(I-1) = 0.0E0
  560       U = D(I) - X1 - XU * V
            IF (I .NE. Q) V = E(I+1)
  580    CONTINUE
C
         IF (U .EQ. 0.0E0) U = EPS3
         RV1(Q) = U
         RV2(Q) = 0.0E0
         RV3(Q) = 0.0E0
C     .......... BACK SUBSTITUTION
C                FOR I=Q STEP -1 UNTIL P DO -- ..........
  600    DO 620 II = P, Q
            I = P + Q - II
            RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
            V = U
            U = RV6(I)
  620    CONTINUE
C     .......... ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                MEMBERS OF GROUP ..........
         IF (GROUP .EQ. 0) GO TO 700
C
         DO 680 JJ = 1, GROUP
            J = R - GROUP - 1 + JJ
            XU = 0.0E0
C
            DO 640 I = P, Q
  640       XU = XU + RV6(I) * Z(I,J)
C
            DO 660 I = P, Q
  660       RV6(I) = RV6(I) - XU * Z(I,J)
C
  680    CONTINUE
C
  700    NORM = 0.0E0
C
         DO 720 I = P, Q
  720    NORM = NORM + ABS(RV6(I))
C
         IF (NORM .GE. 1.0E0) GO TO 840
C     .......... FORWARD SUBSTITUTION ..........
         IF (ITS .EQ. 5) GO TO 960
         IF (NORM .NE. 0.0E0) GO TO 740
         RV6(S) = EPS4
         S = S + 1
         IF (S .GT. Q) S = P
         GO TO 780
  740    XU = EPS4 / NORM
C
         DO 760 I = P, Q
  760    RV6(I) = RV6(I) * XU
C     .......... ELIMINATION OPERATIONS ON NEXT VECTOR
C                ITERATE ..........
  780    DO 820 I = IP, Q
            U = RV6(I)
C     .......... IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
C                WAS PERFORMED EARLIER IN THE
C                TRIANGULARIZATION PROCESS ..........
            IF (RV1(I-1) .NE. E(I)) GO TO 800
            U = RV6(I-1)
            RV6(I-1) = RV6(I)
  800       RV6(I) = U - RV4(I) * RV6(I-1)
  820    CONTINUE
C
         ITS = ITS + 1
         GO TO 600
C     .......... NORMALIZE SO THAT SUM OF SQUARES IS
C                1 AND EXPAND TO FULL ORDER ..........
  840    U = 0.0E0
C
         DO 860 I = P, Q
  860    U = PYTHAG(U,RV6(I))
C
         XU = 1.0E0 / U
C
         DO 880 I = 1, N
  880    Z(I,R) = 0.0E0
C
         DO 900 I = P, Q
  900    Z(I,R) = RV6(I) * XU
C
         X0 = X1
  920 CONTINUE
C
  940 IF (Q .LT. N) GO TO 100
      GO TO 1001
C     .......... SET ERROR -- NON-CONVERGED EIGENVECTOR ..........
  960 IERR = 4 * N + R
      GO TO 1001
C     .......... SET ERROR -- UNDERESTIMATE OF NUMBER OF
C                EIGENVALUES IN INTERVAL ..........
  980 IERR = 3 * N + 1
 1001 LB = T1
      UB = T2
      RETURN
      END
