c      
      subroutine strfac(hi,ki,li,alamda,psi)
c
c --- fhkl fuer si
c
      real*8 hi,ki,li
      real*4 h,k,l
      real*4 f,fp,fpp,dh
      real*4 astar,vc,rcl,alamda
      real*4 energy,gfak
      real*4 pi,zpi 
      complex*8 psi
      complex*8 fhkl,sk
c
      pi = 4.e0*atan(1.e0)
      zpi = 2.e0*pi
c      
cc      sk = 1.e0 +cexp((0.e0,1.e0)*2.d0*pi*(h+k+l))
c      
      h = real(hi)
      k = real(ki)
      l = real(li)
      sk = cexp((0.e0,1.e0)*zpi*(0.*h+0.*k+0.*l))+
     .     cexp((0.e0,1.e0)*zpi*(0.*h+1./2.*k+1./2.*l))+
     .     cexp((0.e0,1.e0)*zpi*(1./2.*h+0.*k+1./2.*l))+
     .     cexp((0.e0,1.e0)*zpi*(1./2.*h+1./2.*k+0.*l))+
     .     cexp((0.e0,1.e0)*zpi*(1./4.*h+1./4.*k+1./4.*l))+
     .     cexp((0.e0,1.e0)*zpi*(1./4.*h+3./4.*k+3./4.*l))+
     .     cexp((0.e0,1.e0)*zpi*(3./4.*h+1./4.*k+3./4.*l))+
     .     cexp((0.e0,1.e0)*zpi*(3./4.*h+3./4.*k+1./4.*l))

      gfak = cabs(sk)
      if((abs(h)+abs(k)+abs(l)).ne.0.) then 
        call debwal(h,k,l,293.50,dh)
      else 
        dh = 1.e0
      endif
c      write(6,*)'sk :',cabs(sk),' dh :',dh
      astar = 5.43102018
      vc = astar**3
      rcl=2.81774e-5   
c      
      call atomfor(h,k,l,alamda,f,fp,fpp)
c
      fhkl = cmplx((f+fp),fpp)*dh
c
      psi=-dcmplx(((rcl*alamda**2)/(pi*(vc)))*fhkl*gfak)
      return 
      end
c
      subroutine atomfor(h,k,l,alamda,f,fp,fpp)
C
C --- Interpolation Atomarer Streufaktoren f,fp,fpp
C     fuer Silizium
C     fp, fpp sind mit absorbtion auf der DesyVax ermittelt
C     Nach Cromer - Liebermann
C     Die Daten fuer f sind aus den International Tables of
C     X-Ray Diffraction Kap. 3.3 ff
C
C     Die fuer f' und f'' sind auf der DESY Vax mit       
C     den Routinen zum Packet Absorption bestimmt.
C
      implicit none
      real thovlm(20)
      real MAtScF(20)
      real y_spline(20)
      real x,f,fp,fpp,formfcal
      real a0,hkl
      real alamda
      real h,k,l
      integer ih,ik,il
      integer i
      logical start
      data start/.true./
      data a0/5.43102018e0/  !Gi Kn Silizium
      data thovlm/0.e0,0.05e0,0.10e0,0.15e0,0.20e0,0.25e0
     +           ,0.30e0,0.35e0,0.40e0,
     +    0.50e0,0.60e0,0.70e0,0.80e0,0.90e0,1.00e0,1.10e0,1.20e0,
     +                     1.30e0,1.40e0,1.50e0/
      data MAtScF/14.00e0,13.45e0,12.16e0,10.79e0,9.67e0
     +           ,8.84e0,8.22e0,7.70e0,
     +     7.20e0,6.24e0,5.31e0,4.47e0,3.75e0,3.16e0,2.96e0,2.35e0
     +    ,2.07e0,1.87e0,1.71e0,1.60e0/
c
c      write(6,*) ' atomfor ..'
        call spline(thovlm,MAtScF,20,2.e30,2.e30,y_spline)
c
      ih = int(h)
      ik = int(k)
      il = int(l)
      hkl = (h*h+k*k+l*l)
      x = sqrt(hkl)/(2.e0*a0)
c      write(6,*) ' x ' ,x 
c      
      call splint(thovlm,MAtScF,y_spline,20,x,f)
c
c      if((h+k+l).ne.0.) call fspline(12349./alamda) 
      if(start)then 
        call fspline(12349./alamda) 
        start=.false.
      endif
c      
      call fdata(12349./alamda,fp,fpp)
c
      return
      end
cc
      subroutine fspline(energy)
c
c --- berechnung von 100 spline koeffizienten in der neahe der 
c     gesuchten energie 
c     energy muss durch 5 teilbar sein
c     
      real energy, entest,testen,rdum
      real en(100),fp(100),fpp(100)
      real fps(100),fpps(100)
      integer ndata
      
      common/transfer/en,fp,fpp,fps,fpps  
      
      open(11,file='/dosc/fprime/fppp14.dat')
c 
c --- einlesen bis zu einer energie 50*5 eV unterhalb der 
c     gewuenschten energie
c
      energy = real(int(energy/5.)+1.)*5.      

      entest = energy - 50.*5. 
      rewind(11)
 5    continue
      read(11,*,end=100) testen,rdum,rdum
      if(testen.ne.entest) goto 5
c
c --- dann die naechsten hundert Werte fuer die Tabelle      
c     einlesen
c
      i = 0
 10   continue
      i = i + 1
      read(11,*,err=100) en(i),fp(i),fpp(i)
      if(i.ne.100) goto 10
      ndata = i
c
c --- spline koeffizienten bestimmen
c
      call SPLINE(en,fp,ndata,2e30,2e30,fps)
      call SPLINE(en,fpp,ndata,2e30,2e30,fpps)
c
      close(11)
      return
 100  write(6,*) 'Ein Fehler in fspline(energy) ...'
      return
      end
cc
       subroutine fdata(energy,fpo,fppo)
c     ermittelt mit splint werte fuer fp,fpp aus
c --- tabellierten werten von fp, fpp fuer si
c     ermittelt mit absorbtion auf vxdesy im bereich 250 ~ 10000 eV
c     
      real*4 energy,fpo,fppo
      real*4 en(100),fp(100),fpp(100)
      real*4 fps(100),fpps(100)
c      
      common/transfer/en,fp,fpp,fps,fpps
c
c --- fp und fpp mit splint berechnen
c
      call splint(en,fp,fps,100,energy,fpo)
      call splint(en,fpp,fpps,100,energy,fppo)
c
      return
      end
c
c --- aus Numerical Recipes  Press et all Cambridge University Press 1989 
c
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) then
         write(6,*) 'Bad XA input.'
         return
      endif
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END
C
C --- Spline
C
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      PARAMETER (NMAX=100)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF (YP1.GT..99d30) THEN
          Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.d0)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99e30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END
      
      

