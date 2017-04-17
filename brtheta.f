      subroutine brtheta                                                        
c
c --- geandert, es werden die chi werte aus unserem dyth program 
c     verwendet
c              16/02/93
                                                                                
      implicit none                                                  
                                                                                      
      include 'param.inc'
      integer i,j      
      real*8 ux(3),v(3),w(3)                             
      real*8 rcl,esse2,vcstar,vc,zeta,zeta1                            
      real*8 zero                                                               
      real*8 ch,ck,cl
      complex*16 fhkl                                     
      
      include 'c2.inc'                                                                          
      
      data zero/0.d0/                                                            
      pi=4.*datan(1.d0)                                                         
c      write(*,*)'      Begin brtheta...'                                       
      rcl=2.81774d-5                                                            
      ux(1)=astar                                                               
      ux(2)=0.d0                                                                
      ux(3)=0.d0                                                                
      v(1)=0.d0                                                                 
      v(2)=bstar                                                                
      v(3)=0.d0                                                                 
      call cross(ux,v,w)                                                        
      ux(1)=0.d0                                                                
      ux(2)=0.d0                                                                
      ux(3)=cstar                                                               
      call dot(w,ux,vcstar)                                                     
      vc=1.d0/vcstar                                                            
   74 format(3h  i,2h j,4x,5hh k l,5x,4hesse,4x,5x,6x,4hfhkl,8x,                
     1 13x,5x,3hpsi)                                                            
   60 format(1h*)                                                               
c     computation of psi(i,j)                                                   
      if(iphi.eq.1)print 74                                                     
      do 150 i=1,n                                                              
      do 150 j=1,n                                                                                                                                               
      if(i.eq.j)go to 150                                                       
      ch=h(i)-h(j)                                                              
      ck=k(i)-k(j)                                                              
      cl=l(i)-l(j)                                                                                                         
      call strfac(ch,ck,cl,real(alamda),fhkl)                                                
      psi(i,j) = fhkl
ccccc      write(6,*)' BRTHETA : ',ch,ck,cl,psi(i,j)                                                
 919  format('BRTHETA: ',3(f6.2,1x),2x,2f10.4)                                   
                                                                                
      ux(1)=astar*ch                                                            
      ux(2)=bstar*ck                                                            
      ux(3)=cstar*cl                                                            
      call dot(ux,ux,esse2)                                                     
      esse=0.5d0*dsqrt(esse2)                                                   
cc      psi(i,j)=-(rcl*(alamda)**2/(pi*vc))*fhkl*8.d0 
      write(6,*) psi(i,j)
cc      psi(i,j)=(-6.950058795054567E-06,-4.221856539822649E-07)
ccc      psi(i,j) =  (-.6752d-5,-.4158d-6)                                  
      if(iphi.eq.1)print 76,i,j,ch,ck,cl,esse,fhkl,psi(i,j)                     
   76 format(i3,i2,3(f4.1,1x),f9.4,2x,2f10.4,2x,2e13.5)                          
  150 continue               
      call strfac(0.,0.,0.,real(alamda),fhkl)                                          
      psio = fhkl
      write(6,*)zero,zero,zero,psio                                          
cc      psio=-(rcl*(alamda)**2/(pi*vc))*fhkl*8.d0 
ccc      psio = (-.1792d-4,-.4692d-6)
cc      psio = (-1.740014922413449E-05,-1.130706373945842E-07)
cc      write(6,*)zero,zero,zero,psio                                        
      ux(1)=h(2)*astar                                                          
      ux(2)=k(2)*bstar                                                          
      ux(3)=l(2)*cstar                                                          
      call dot(ux,ux,esse2)                                                     
      esse=0.5d0*dsqrt(esse2)                                                   
      sintb=esse*alamda                                                         
      costb=dsqrt(1.d0-sintb*sintb)                                             
      thetab=dasin(sintb)                                                       
      eo=-psio/(2.d0*sintb*costb)                                               
      thetao=thetab+eo                                                          
      dtheta=-dsqrt((dble(psi(2,1))/(2.d0*sintb*costb))**2                     
     &            + (dimag(psi(2,1))/(2.d0*sintb*costb))**2)                    
c                                                                               
c                                                                               
c      if(psi(2,1).eq.0.)dtheta=(pi/180.)/1000.                                 
c                                                                               
c                                                                               
      dtheta=dabs(dtheta)                                                       
      dtheta=ifac*dtheta                                                        
cc spezial spezial       
cc      dtheta = 9.69627e-5
c     values of eo and dtheta will be modified later for the                    
c     asymmetric case.                                                          
      if(iphi.eq.1)print 60                                                     
c     transformation of coordinates                                             
c     y - axis = normal inward                                                  
c     x - axis = parallel to ko tang                                            
c     z - axis = ux x uy                                                        
      bhx(1)=0.d0                                                               
      bhy(1)=0.d0                                                               
      bhz(1)=0.d0                                                               
c     define t-vector, component of m-vector perpendicular                      
c     to scattering vector h(2), k(2), l(2).                                    
      ux(1)=h(2)*astar                                                          
      ux(2)=k(2)*bstar                                                          
      ux(3)=l(2)*cstar                                                          
      v(1)=mh*astar                                                             
      v(2)=mk*bstar                                                             
      v(3)=ml*cstar                                                             
      call dot(ux,v,zeta)                                                       
      call dot(ux,ux,zeta1)                                                     
      t(1)=v(1)-ux(1)*zeta/zeta1                                                
      t(2)=v(2)-ux(2)*zeta/zeta1                                                
      t(3)=v(3)-ux(3)*zeta/zeta1                                                
c     define xi,yi,zi-vectors. these are orthogonal unit vectors                
c     corresponding to the x,y,z axes defined after statement 61.               
c     their contravariant components are referred to the                        
c     astar, bstar, cstar axes of the crystal.                                  
c     y-vector                                                                  
      ux(1)=ph*astar                                                            
      ux(2)=pk*bstar                                                            
      ux(3)=pl*cstar                                                            
      call dot(ux,ux,zeta)                                                      
      zeta=dsqrt(zeta)                                                          
      yi(1)=-ux(1)/zeta                                                         
      yi(2)=-ux(2)/zeta                                                         
      yi(3)=-ux(3)/zeta                                                         
c      write(*,*)'Calling angles(thetab)...',thetab                                    
      call angles(thetab)                                                       
      write(*,1121)thetax, 180.*thetax/pi                                       
 1121 format('thetax=',e13.5,'radians = ',f9.3,                                 
     & ' degrees ****************')                                             
 533  continue                                                                  
      dtheta=dtheta/dsqrt(dabs(bi))                                             
      eo=eo*(bi-1.d0)/(2.d0*bi)                                                 
      thetao=thetab+eo                                                          
c      write(*,*)'          End brtheta...'                                     
      return                                                                    
      end   
