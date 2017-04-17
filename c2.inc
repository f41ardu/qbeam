c
c --- neudeklaration aller common variablen
c     Copyright Th. Rautenstrauch
c        
        real*8 astar,bstar,cstar,e0,ph,pk,pl,t(3),xi(3),
     .         yi(3),zi(3),pi,alamda,thetao,dtheta,
     .         thetx,thetax,to,esse,sintb,eo,
     .         costb,thetab,kprp,kpar,alpha,zzeta,zz1,zz2,
     .         mh,mk,ml,phi,h(nsize),k(nsize),l(nsize),
     .         bhx(8),bhy(8),bhz(8),phio(10),bi,g0
        complex*16 psi(nsize,nsize),psio
        integer n,ipvec,iter,isinv,ncal,isprp,ispar,ik,jflag,kflag,
     .          izscan,ifac,iphi,ntheta,nz,nscan,nprint
c        
        include 'c2com.inc'
