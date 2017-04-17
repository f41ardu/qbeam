c
        function myexp(arg)
c
c--- Zerlegung der complexen Exponentialfunktion
c    Auswertung der trignometischen Funktionen ueber 
c    Hauptwerte  
c
        complex*16 myexp,arg
        real*8 rarg,iarg,pi
        real*8 resultrarg,resultiarg
        real*8 emuet
c
        pi = 4.d0*atan(1.d0)
c
        myexp = (0.d0,0.d0)
c        
        rarg = dble(arg)
        iarg = dimag(arg)
c
        iarg = mod(iarg,2.d0*pi)
c        
        if(rarg.ge.0.d0) then
           if(rarg.ge.15.d0) then 
             emuet = dexp(15.d0)
           else
             emuet=dexp(rarg)
           endif
        endif
c        
        if(rarg.lt.0.d0) then
           if(rarg.le.-15.d0) then
              emuet = dexp(-15.d0)
           else 
              emuet=dexp(rarg)
           endif
        endif
c        
        resultrarg=dcos(iarg)
        resultiarg=dsin(iarg)
c
        myexp = emuet*dcmplx(resultrarg,resultiarg)
        return
c        
        end
 
