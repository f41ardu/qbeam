      subroutine angles(theta)                                                  
                                                                                
      implicit none                                                  
      include 'param.inc'

      real*8 det(2),work(3)                                                     
      integer ipvt(3),job,job1,info                                             
                                                                                
      real*8 theta
      real*8 hmod,khtxww,khtzww,khy                        
      real*8 rkflag                             
      real*8 a1,a2,a3,b1,b2,b3,z,z2,opmod,tmod                         
      real*8 rc(3,3),rcxp(3),gammao,gammah                                      
      integer i,three                                                         
      
      real*8 ux(3),v(3),o(3),w(3),op(3)                  
      
      include 'c2.inc'
      data three /3/                                                      
      data job,job1 /01,0/                                                      
                                                                                
c     theta is the angle between o-vector and h(2),k(2),l(2) planes.            
c     thetax is the angle between o-vector and xtal surface.                    
c     determine contravariant components of x-unit vector. to do                
c     this, we first calculate the contravariant components of o-unit           
c     vector, opposite but otherwise coincident with ko, incident               
c     wavevector.                                                               
c     first we consider plane determined by h(2),k(2),l(2) vector               
c     and o(1),o(2),o(3) vector.                                                
c     o-unit vector is determined by three equations:                           
c                                                                               
c     (o x h) . t = tmod.hmod.sin(phi).cos(theta)                               
c                                                                               
c     o . t = tmod.cos(theta).cos(phi)                                          
c                                                                               
c     o . h = hmod.sin(theta)                                                   
c                                                                               
c     to understand above three equations, consider stereographic               
c     projection with h(2), k(2), l(2) vector (called simply h in               
c     above equations) at center of projection, o is a point on the             
c     equator, to the right of center. angle between h-vector and               
c     o-vector is pi/2 - theta. t-vector is on circumference,                   
c     rotated clockwise from o-side of equator by phi angle.                    
c      write(*,*)'        Beginning angles...'                                  
c      write(6,*) e0
      ux(1)=1.d0                                                                
      ux(2)=0.d0                                                                
      ux(3)=0.d0                                                                
      v(1)=h(2)*astar                                                           
      v(2)=k(2)*bstar                                                           
      v(3)=l(2)*cstar                                                           
      call dot(ux,t,a1)                                                         
      call dot(ux,v,b1)                                                         
      ux(1)=0.d0                                                                
      ux(2)=1.d0                                                                
      call dot(ux,t,a2)                                                         
      call dot(ux,v,b2)                                                         
      ux(2)=0.d0                                                                
      ux(3)=1.d0                                                                
      call dot(ux,t,a3)                                                         
      call dot(ux,v,b3)                                                         
c                                                                               
      rc(1,1)=v(2)*t(3)-v(3)*t(2)                                               
      rc(1,2)=v(3)*t(1)-v(1)*t(3)                                               
      rc(1,3)=v(1)*t(2)-v(2)*t(1)                                               
      call dot(v,v,z)                                                           
      hmod=dsqrt(z)                                                             
      call dot(t,t,z)                                                           
      tmod=dsqrt(z)                                                             
      z=(tmod*hmod/e0)*dsin(phi)*dcos(theta)                                    
      rcxp(1)=z                                                                 
      rc(2,1)=a1                                                                
      rc(2,2)=a2                                                                
      rc(2,3)=a3                                                                
      z=tmod*dcos(phi)*dcos(theta)                                              
      rcxp(2)=z                                                                 
      rc(3,1)=b1                                                                
      rc(3,2)=b2                                                                
      rc(3,3)=b3                                                                
      z=hmod*dsin(theta)                                                         
      rcxp(3)=z                                                                 
c      write(6,*) 'vor linpack '
c                                                                               
c     Use Linpack for matrix inversion, vectors...                              
c                                                                               
c      write(*,*)'Call dgefa, dgesl, dgedi...'                                  
      call sgefa(rc,three,three,ipvt,info)                                      
      call sgesl(rc,three,three,ipvt,rcxp,job1)                                 
      call sgedi(rc,three,three,ipvt,det,work,job)                              
                                                                                
c                                                                               
c     <<call lineq4(c,cxp,xxc,3,3,1,ipvec)>> (CDC routine)                      
c                                                                               
                                                                                
      o(1)=rcxp(1)                                                              
      o(2)=rcxp(2)                                                              
      o(3)=rcxp(3)                                                              
c     consider now plane determined by ph,pk,pl vector and o-vector.            
c     determine op(1), op(2), op(3) vector, component of o-vector               
c     normal to ph, pk, pl vector.                                              
      w(1)=ph*astar                                                             
      w(2)=pk*bstar                                                             
      w(3)=pl*cstar                                                             
      call dot(w,w,z2)                                                          
      call dot(o,w,z)                                                           
      op(1)=o(1)-w(1)*z/z2                                                      
      op(2)=o(2)-w(2)*z/z2                                                      
      op(3)=o(3)-w(3)*z/z2                                                      
      call dot(op,op,z)                                                         
      opmod=dsqrt(z)                                                            
      xi(1)=-op(1)/opmod                                                        
      xi(2)=-op(2)/opmod                                                        
      xi(3)=-op(3)/opmod                                                        
c     now determine z-unit vector (contravariant components).                   
c                                                                               
c     z = x cross y                                                             
c                                                                               
      call cross(xi,yi,zi)                                                      
      do 25 i=2,n                                                               
      ux(1)=h(i)*astar                                                          
      ux(2)=k(i)*bstar                                                          
      ux(3)=l(i)*cstar                                                          
      call dot(ux,xi,bhx(i))                                                    
      call dot(ux,yi,bhy(i))                                                    
   25 call dot(ux,zi,bhz(i))                                                    
c     calculate thetax                                                          
      call dot(o,xi,z)                                                          
      thetax=dacos(-z)                                                          
c      write(*,1121)thetax                                                      
 1121 format('thetax=',e13.5)                                                   
c     calculate bi = gammao/gammah (zac. 3.115)                                 
      gammao=dsin(thetax)                                                       
      g0=dabs(gammao)
      khtxww=(1.d0/alamda)*dcos(thetax)+bhx(2)                                  
      khtzww=bhz(2)                                                             
      kflag=0                                                                   
      rkflag=1.d0/alamda**2-khtxww**2-khtzww**2                                 
      if(rkflag.le.0.d0)go to 600                                                 
      khy=-dsqrt(1.d0/alamda**2-khtxww**2-khtzww**2)                            
c     khy = normal component of diffracted beam wavevector n.2                  
c     along y-axis (in vacuum).                                                 
      gammah=khy*alamda                                                         
      bi=gammao/gammah                                                          
      go to 700                                                                 
  600 continue                                                                  
      kflag=1                                                                   
      write(*,*)'kflag=1; negative exit angle!!'                                
  700 continue                                                                  
c      write(*,*)'              Ending angles...'                               
      return                                                                    
      end    
