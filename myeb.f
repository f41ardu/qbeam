      subroutine debwal (h,k,l,temp,dh)
      IMPLICIT NONE

      REAL*4 temp,atomff,s,d,dh,thb,rlambda,
     *       hkl,temp1,debtemp,DEBYETP,atomm,atomms
      real*4 h,k,l
      INTEGER n
C
C   Bestimmung von Wellenlaenge,Braggwinkel,s=SIN(thb)/lambda
C    Debye-W.-Faktor dh, Atommasse atomm [u], Gitterkonstante d [A] :
C
      d = 5.43102018
      atomm=28.086
      hkl = sqrt(h*h+k*k+l*l)
c germanium      atomm=72.59
C Atommasse in 1.E-27 [kg] :
      atomm=atomm*1.6605655
      debtemp=DEBYETP(14)
c     
      CALL DEBYEW(d/hkl,debtemp,atomm,temp,dh)
c
c
      return
      end
C
      SUBROUTINE DEBYEW(d,td,am,to,dh)

C---------------------------------------------------------------------C
C                                                                     C
C DEBYEWaller berechnet den D.-W.-Faktor dh, siehe auch PHI.SUB  .    C
C   d : Netzebenenabstand ;        [ Angstroem ] = 1.E-10 [ m ]       C
C  td : Debyetemperatur ;          [ K ]                              C
C  to : Temperatur ;               [ K ]                              C
C  am : Atommasse ;         1.E-27 [ kg ]                             C
C  dh : Debye-Waller-Faktor ;      [ dimensionslos ]                  C
C                                                                     C
C---------------------------------------------------------------------C
      REAL d,td,am,to,dh,h,kb,pi,sum,x,m
C   Festlegung der Konstanten :
         pi=4.e0*atan(1.e0)
         h=6.626196
C    Wirkungsquantum * 1.E34
         kb=1.380622
C    Boltzmann-Konst. * 1.E23
C  
        x=td/to
        sum=.25*x+PHI(x)
        m=6.*h*h*to/(am*kb*td*td)*sum
        dh=EXP(-m/(4.*d*d)*100.)
C Faktor 100 durch die Konstanten !
   20 CONTINUE
      RETURN
      END
C
      REAL FUNCTION DEBYETP(imat)

C---------------------------------------------------------------------C
c Vers. Th. Rautenstrauch                                             C
c Lit.: Salter,Adv.Phys. 14,1(1965) , Aldred & Hart (1973),           C
c       Batterman & Chipman,Phys.Rev.127,690(1962)                    C
c                                                                     C
C---------------------------------------------------------------------C
      IMPLICIT NONE

      INTEGER imat
      REAL*4    dt,DEBYETAB

      IF (imat.EQ.14) THEN
           dt=542.5
      ELSE IF (imat.EQ.32) THEN
         dt=295.
      ENDIF
C
      DEBYETP=dt
      RETURN
      END
C
      REAL FUNCTION PHI(z)

C------------------------------------------------------------------C
C Modifiziert Thomas Rautenstrauch 19/02/93                        C 
C    PHI(z) berechnet fuer den Debye-Waller Faktor die  Groesse    C
C    (1/z)*(Integral von 0...z ueber x/(EXP(x)-1)*dx)              C
C    mit z = Debyetemperatur/Temperatur                            C
C------------------------------------------------------------------C

      IMPLICIT NONE

      REAL y,z,pi,ph,z1
      real e,en,s,u,a

      pi=4.E0*ATAN(1.E0)
      IF (z.LE..08) THEN
        y=4.*SQRT(.6)*ATAN(z/SQRT(.6)/(z+4.))
       ELSE IF (z.GT.30.) THEN
        y=pi*pi/6.
       ELSE
        u=1.
        e=EXP(-z)
        en=e
        s=e
   99    u=u+1.
         en=en*e
         a=en/(u*u)
         s=s+a
         IF (a*1e8.GT.s) GOTO 99
        y=pi*pi/6.+z*LOG(1.-EXP(-z))-s
      ENDIF
  999 IF (z.NE.0) THEN
           y=y/z
        ELSE
           y=1.
      ENDIF
      PHI=y
      RETURN
      END

