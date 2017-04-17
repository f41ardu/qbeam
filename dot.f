      subroutine dot(u,v,zeta)                                                  
                                                                                
      implicit none                                                  

      integer i,j
c     u, v, are 3-dimensional vectors in general contravariant                  
c     coordinates.                                                              
c     ao(3,3) is the fundamental tensor.                                        
c     e(3,3,3) is the ricci tensor (not used in this routine).                  
c     at(3,3) is the inverse matrix of a(3,3), (not used in                     
c     this routine).                                                            
c     z = u . v                                                                 
      real*8 u(3),v(3),sum,zeta                        
                                                                                
      include 'c1.inc'

      sum=0.d0                                                                  
      do 30 i=1,3                                                               
      do 30 j=1,3                                                               
   30 sum=sum+ao(i,j)*u(i)*v(j)                                                 
      zeta=sum                                                                  
      return                                                                    
      end                                                                       
c           
