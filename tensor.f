      subroutine tensor(alphax,betax,gammax,ao,e,ai,e0)                         
										
      implicit none                                                  
      include 'param.inc'
      integer i,j,ik
      real*8 alphax,betax,gammax
      real*8 ao(3,3),e(3,3,3),ai(3,3),e0
      real*8 det,cax,cbx,ccx                                                    
      real*8 pi
      pi=4.*datan(1.0d0)                                                        
      alphax=alphax*pi/180.d0                                                   
      betax=betax*pi/180.d0                                                     
      gammax=gammax*pi/180.d0                                                   
      cax=dcos(alphax)                                                          
      cbx=dcos(betax)                                                           
      ccx=dcos(gammax)                                                          
c     ai(3,3) is the inverse matrx of ao(3,3).                                  
c     construction of fundamental tensor (covariant components).                
      do 455 i=1,3                                                              
      do 455 j=1,3                                                              
  455 ao(i,j)=1.d0                                                              
      ao(1,2)=ccx                                                               
      ao(2,1)=ao(1,2)                                                           
      ao(1,3)=cbx                                                               
      ao(3,1)=ao(1,3)                                                           
      ao(2,3)=cax                                                               
      ao(3,2)=ao(2,3)                                                           
      det= 1.d0+2.d0*cax*cbx*ccx-cax**2-cbx**2-ccx**2                           
c     construction of ricci tensor (covariant components).                      
      e0=dsqrt(det)                                                             
      do 456 i=1,3                                                              
      do 456 j=1,3                                                              
      do 456 ik=1,3                                                             
      e(i,j,ik)=0.d0                                                            
      if(i.ne.j.and.j.ne.ik.and.ik.ne.i) e(i,j,ik)=e0                           
  456 continue                                                                  
      e(2,1,3)=-e0                                                              
      e(1,3,2)=-e0                                                              
      e(3,2,1)=-e0                                                              
c     contravariant components of fundamental tensor.                           
      ai(1,1)=(1.d0-cax*cax)/det                                                
      ai(1,2)=(cax*cbx-ccx)/det                                                 
      ai(1,3)=(cax*ccx-cbx)/det                                                 
      ai(2,3)=(cbx*ccx-cax)/det                                                 
      ai(2,2)=(1.d0-cbx*cbx)/det                                                
      ai(3,3)=(1.d0-ccx*ccx)/det                                                
      ai(2,1)=ai(1,2)                                                           
      ai(3,1)=ai(1,3)                                                           
      ai(3,2)=ai(2,3)                                                           
      return                                                                    
      end   
