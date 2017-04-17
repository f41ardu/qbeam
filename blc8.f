C   23/11/92 211231441  MEMBER NAME  BLC8     (EISPACK.)    FVS
      SUBROUTINE CAXPY(N,CA,CX,INCX,CY,INCY)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      COMPLEX CX(1),CY(1),CA
      INTEGER I,INCX,INCY,IX,IY,N
C
      IF(N.LE.0)RETURN
      IF (ABS(REAL(CA)) + ABS(AIMAG(CA)) .EQ. 0.0 ) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        CY(IY) = CY(IY) + CA*CX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        CY(I) = CY(I) + CA*CX(I)
   30 CONTINUE
      RETURN
      END
      SUBROUTINE  CCOPY(N,CX,INCX,CY,INCY)
C
C     COPIES A VECTOR, X, TO A VECTOR, Y.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      COMPLEX CX(1),CY(1)
      INTEGER I,INCX,INCY,IX,IY,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        CY(IY) = CX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        CY(I) = CX(I)
   30 CONTINUE
      RETURN
      END
      COMPLEX FUNCTION CDOTC(N,CX,INCX,CY,INCY)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS, CONJUGATING THE FIRST
C     VECTOR.
C     JACK DONGARRA, LINPACK,  3/11/78.
C
      COMPLEX CX(1),CY(1),CTEMP
      INTEGER I,INCX,INCY,IX,IY,N
C
      CTEMP = (0.0,0.0)
      CDOTC = (0.0,0.0)
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        CTEMP = CTEMP + CONJG(CX(IX))*CY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      CDOTC = CTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        CTEMP = CTEMP + CONJG(CX(I))*CY(I)
   30 CONTINUE
      CDOTC = CTEMP
      RETURN
      END
      COMPLEX FUNCTION CDOTU(N,CX,INCX,CY,INCY)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      COMPLEX CX(1),CY(1),CTEMP
      INTEGER I,INCX,INCY,IX,IY,N
C
      CTEMP = (0.0,0.0)
      CDOTU = (0.0,0.0)
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        CTEMP = CTEMP + CX(IX)*CY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      CDOTU = CTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        CTEMP = CTEMP + CX(I)*CY(I)
   30 CONTINUE
      CDOTU = CTEMP
      RETURN
      END
      REAL FUNCTION CMACH(JOB)
      INTEGER JOB
C
C     SMACH COMPUTES MACHINE PARAMETERS OF FLOATING POINT
C     ARITHMETIC FOR USE IN TESTING ONLY.  NOT REQUIRED BY
C     LINPACK PROPER.
C
C     IF TROUBLE WITH AUTOMATIC COMPUTATION OF THESE QUANTITIES,
C     THEY CAN BE SET BY DIRECT ASSIGNMENT STATEMENTS.
C     ASSUME THE COMPUTER HAS
C
C        B = BASE OF ARITHMETIC
C        T = NUMBER OF BASE  B  DIGITS
C        L = SMALLEST POSSIBLE EXPONENT
C        U = LARGEST POSSIBLE EXPONENT
C
C     THEN
C
C        EPS = B**(1-T)
C        TINY = 100.0*B**(-L+T)
C        HUGE = 0.01*B**(U-T)
C
C     DMACH SAME AS SMACH EXCEPT T, L, U APPLY TO
C     DOUBLE PRECISION.
C
C     CMACH SAME AS SMACH EXCEPT IF COMPLEX DIVISION
C     IS DONE BY
C
C        1/(X+I*Y) = (X-I*Y)/(X**2+Y**2)
C
C     THEN
C
C        TINY = SQRT(TINY)
C        HUGE = SQRT(HUGE)
C
C
C     JOB IS 1, 2 OR 3 FOR EPSILON, TINY AND HUGE, RESPECTIVELY.
C
C
      REAL EPS,TINY,HUGE,S
C
      EPS = 1.0
   10 EPS = EPS/2.0
      S = 1.0 + EPS
      IF (S .GT. 1.0) GO TO 10
      EPS = 2.0*EPS
      CMACH =EPS
      IF( JOB .EQ. 1) RETURN
C
      S = 1.0
   20 TINY = S
      S = S/16.0
      IF (S*1.0 .NE. 0.0) GO TO 20
      TINY = (TINY/EPS)*100.
      S = REAL((1.0,0.0)/CMPLX(TINY,0.0))
      IF (S .NE. 1.0/TINY) TINY = SQRT(TINY)
      HUGE = 1.0/TINY
      IF (JOB .EQ. 1) CMACH = EPS
      IF (JOB .EQ. 2) CMACH = TINY
      IF (JOB .EQ. 3) CMACH = HUGE
      RETURN
      END
      SUBROUTINE CROTG(CA,CB,C,S)
      COMPLEX CA,CB,S
      REAL C
      REAL NORM,SCALE
      COMPLEX ALPHA
      IF (CABS(CA) .NE. 0.) GO TO 10
         C = 0.
         S = (1.,0.)
         CA = CB
         GO TO 20
   10 CONTINUE
         SCALE = CABS(CA) + CABS(CB)
         NORM = SCALE * SQRT((CABS(CA/SCALE))**2 + (CABS(CB/SCALE))**2)
         ALPHA = CA /CABS(CA)
         C = CABS(CA) / NORM
         S = ALPHA * CONJG(CB) / NORM
         CA = ALPHA * NORM
   20 CONTINUE
      RETURN
      END
      SUBROUTINE  CSCAL(N,CA,CX,INCX)
C
C     SCALES A VECTOR BY A CONSTANT.
C     JACK DONGARRA, LINPACK,  3/11/78.
C
      COMPLEX CA,CX(1)
      INTEGER I,INCX,N,NINCX
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        CX(I) = CA*CX(I)
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DO 30 I = 1,N
        CX(I) = CA*CX(I)
   30 CONTINUE
      RETURN
      END
      SUBROUTINE  CSROT (N,CX,INCX,CY,INCY,C,S)
C
C     APPLIES A PLANE ROTATION, WHERE THE COS AND SIN (C AND S) ARE REAL
C     AND THE VECTORS CX AND CY ARE COMPLEX.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      COMPLEX CX(1),CY(1),CTEMP
      REAL C,S
      INTEGER I,INCX,INCY,IX,IY,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        CTEMP = C*CX(IX) + S*CY(IY)
        CY(IY) = C*CY(IY) - S*CX(IX)
        CX(IX) = CTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        CTEMP = C*CX(I) + S*CY(I)
        CY(I) = C*CY(I) - S*CX(I)
        CX(I) = CTEMP
   30 CONTINUE
      RETURN
      END
      SUBROUTINE  CSSCAL(N,SA,CX,INCX)
C
C     SCALES A COMPLEX VECTOR BY A REAL CONSTANT.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      COMPLEX CX(1)
      REAL SA
      INTEGER I,INCX,N,NINCX
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        CX(I) = CMPLX(SA*REAL(CX(I)),SA*AIMAG(CX(I)))
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DO 30 I = 1,N
        CX(I) = CMPLX(SA*REAL(CX(I)),SA*AIMAG(CX(I)))
   30 CONTINUE
      RETURN
      END
      SUBROUTINE  CSWAP (N,CX,INCX,CY,INCY)
C
C     INTERCHANGES TWO VECTORS.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      COMPLEX CX(1),CY(1),CTEMP
      INTEGER I,INCX,INCY,IX,IY,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        CTEMP = CX(IX)
        CX(IX) = CY(IY)
        CY(IY) = CTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
   20 DO 30 I = 1,N
        CTEMP = CX(I)
        CX(I) = CY(I)
        CY(I) = CTEMP
   30 CONTINUE
      RETURN
      END
      INTEGER FUNCTION ICAMAX(N,CX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      COMPLEX CX(1)
      REAL SMAX
      INTEGER I,INCX,IX,N
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
      ICAMAX = 0
      IF( N .LT. 1 ) RETURN
      ICAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      SMAX = CABS1(CX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(CABS1(CX(IX)).LE.SMAX) GO TO 5
         ICAMAX = I
         SMAX = CABS1(CX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 SMAX = CABS1(CX(1))
      DO 30 I = 2,N
         IF(CABS1(CX(I)).LE.SMAX) GO TO 30
         ICAMAX = I
         SMAX = CABS1(CX(I))
   30 CONTINUE
      RETURN
      END
      REAL FUNCTION SCASUM(N,CX,INCX)
C
C     TAKES THE SUM OF THE ABSOLUTE VALUES OF A COMPLEX VECTOR AND
C     RETURNS A SINGLE PRECISION RESULT.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      COMPLEX CX(1)
      REAL STEMP
      INTEGER I,INCX,N,NINCX
C
      SCASUM = 0.0E0
      STEMP = 0.0E0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        STEMP = STEMP + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))
   10 CONTINUE
      SCASUM = STEMP
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DO 30 I = 1,N
        STEMP = STEMP + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))
   30 CONTINUE
      SCASUM = STEMP
      RETURN
      END
      REAL FUNCTION SCNRM2( N, CX, INCX)
      LOGICAL IMAG, SCALE
      INTEGER          NEXT
      REAL         CUTLO, CUTHI, HITEST, SUM, XMAX, ABSX, ZERO, ONE
      COMPLEX      CX(1)
      DATA         ZERO, ONE /0.0E0, 1.0E0/
C
C     UNITARY NORM OF THE COMPLEX N-VECTOR STORED IN CX() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON , 1978 JAN 08
C
C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C     HOPEFULLY APPLICABLE TO ALL MACHINES.
C         CUTLO = MAXIMUM OF  SQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C         CUTHI = MINIMUM OF  SQRT(V)      OVER ALL KNOWN MACHINES.
C     WHERE
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C     BRIEF OUTLINE OF ALGORITHM..
C
C     PHASE 1    SCANS ZERO COMPONENTS.
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C
C     VALUES FOR CUTLO AND CUTHI..
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                   UNIVAC AND DEC AT 2**(-103)
C                   THUS CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                   THUS CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
C
      IF(N .GT. 0) GO TO 10
         SCNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C                                                 BEGIN MAIN LOOP
      DO 210 I=1,NN,INCX
         ABSX = ABS(REAL(CX(I)))
         IMAG = .FALSE.
         GO TO NEXT,(30, 50, 70, 90, 110)
   30 IF( ABSX .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      SCALE = .FALSE.
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF( ABSX .EQ. ZERO) GO TO 200
      IF( ABSX .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 ASSIGN 110 TO NEXT
      SUM = (SUM / ABSX) / ABSX
  105 SCALE = .TRUE.
      XMAX = ABSX
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF( ABSX .GT. CUTLO ) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF( ABSX .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / ABSX)**2
         XMAX = ABSX
         GO TO 200
C
  115 SUM = SUM + (ABSX/XMAX)**2
      GO TO 200
C
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
   85 ASSIGN 90 TO NEXT
      SCALE = .FALSE.
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
      HITEST = CUTHI/FLOAT( N )
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
   90 IF(ABSX .GE. HITEST) GO TO 100
         SUM = SUM + ABSX**2
  200 CONTINUE
C                  CONTROL SELECTION OF REAL AND IMAGINARY PARTS.
C
      IF(IMAG) GO TO 210
         ABSX = ABS(AIMAG(CX(I)))
         IMAG = .TRUE.
      GO TO NEXT,(  50, 70, 90, 110 )
C
  210 CONTINUE
C
C              END OF MAIN LOOP.
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      SCNRM2 = SQRT(SUM)
      IF(SCALE) SCNRM2 = SCNRM2 * XMAX
  300 CONTINUE
      RETURN
      END
      SUBROUTINE SCOTG(SA,SB,C,S)
C
C     CONSTRUCT GIVENS PLANE ROTATION.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      REAL SA,SB,C,S,ROE,SCALE,R,Z
C
      ROE = SB
      IF( ABS(SA) .GT. ABS(SB) ) ROE = SA
      SCALE = ABS(SA) + ABS(SB)
      IF( SCALE .NE. 0.0 ) GO TO 10
         C = 1.0
         S = 0.0
         R = 0.0
         GO TO 20
   10 R = SCALE*SQRT((SA/SCALE)**2 + (SB/SCALE)**2)
      R = SIGN(1.0,ROE)*R
      C = SA/R
      S = SB/R
   20 Z = 1.0
      IF( ABS(SA) .GT. ABS(SB) ) Z = S
      IF( ABS(SB) .GE. ABS(SA) .AND. C .NE. 0.0 ) Z = 1.0/C
      SA = R
      SB = Z
      RETURN
      END
