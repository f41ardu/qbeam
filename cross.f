      subroutine cross(u,v,w)                                                   
                                                                                
      implicit none                                                  

      integer i,j,k
c     u, v, w are 3-dimensional vectors in general contravariant                
c     coordinates.                                                              
c     w = u x v                                                                 
c     e(3,3,3) is the ricci tensor                                              
c     ao(3,3) is the fundamental tensor                                         
c     ao(3,3) and e(3,3,3) are in covariant coordinates.                        
c     ai(3,3) is the inverse matrix of ao.                                      
      real*8 u(3),v(3),w(3),wx(3)                      
      real*8 sum                                                                
      include 'c1.inc'                                                                          
      
      do 26 k=1,3                                                               
      sum=0.d0                                                                  
      do 25 i=1,3                                                               
      do 25 j=1,3                                                               
   25 sum=sum+e(i,j,k)*u(i)*v(j)                                                
   26 wx(k)=sum                                                                 
c     transformation from covariant into contravariant                          
c     coordinates.                                                              
      do 28 k=1,3                                                               
      sum=0.d0                                                                  
      do 27 i=1,3                                                               
   27 sum=sum+ai(i,k)*wx(i)                                                     
   28 w(k)=sum                                                                  
      return                                                                    
      end                                                                       
            
