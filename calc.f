      subroutine calc(theta,ro,lo,rh,lh)
										
      implicit none                                                  
      include 'param.inc'                                                                                
				       
      real*4 theta
      integer itr
      integer ipvt(nmax2),job,matz                                                   
      integer i,j,jj,n1,n2,n4,ji,n2p1,ierr                                                                                
      integer n41,info,in1,ix,job1,ikm1
      real*4 u(nsize,3)                        
      real*4 sigma(nsize,3),p(nsize,3)                                          
      real*4 cgamma(nsize),sgamma(nsize)
										
      real*4 aacr(maxsize,maxsize),aaci(maxsize,maxsize)                        
      real*4 lambdr(maxsize),lambdi(maxsize)                                    
      real*4 zr(maxsize,maxsize),zii(maxsize,maxsize)                           
      real*4 fv1(maxsize),fv2(maxsize),fv3(maxsize)                             
      real*4 roprp,ropar
      
      real*4 rhpar(nm1),rhprp(nm1),lhpar(nm1),lhprp(nm1)                        
      real*4 rh(nm1),lh(nm1)                                                    
      real*4 loprp,lopar,lo,ro                                                  
      real*4 rel(maxsize)                                
      real*4 gamma(nsize),khtx(nsize),khtz(nsize)                               
      real*4 sintx,costx,temp,fac,sintx2                               
      real*4 wmax,wmin

      complex*8 work(nmax2),det(2)
      complex*8 xc(nmax2,nmax2),cxp(maxsize)                   
      complex*8 lambda(maxsize),qqq,eoeh(nsize,maxsize)                        
      complex*8 esp(maxsize),aac(maxsize,maxsize),bbc(maxsize)                 
      complex*8 coi(maxsize),zaac(maxsize,maxsize)                             
      complex*8 s_w(maxsize),e(maxsize),v(maxsize,maxsize)
      complex*8 w(maxsize,maxsize)
      complex*8 a(nsize,nsize),b(nsize,nsize),c(nsize,nsize)                   
      complex*8 d(nsize,nsize),ac(nmax2,maxsize)                               
      complex*8 argn,q1,q2,q3,ui(nsize,3,maxsize)                              
      complex*8 uii(nsize,3,maxsize),sigmi(nsize,3,maxsize)                    
      complex*8 pii(nsize,3,maxsize),pihbar(nsize)                             
      complex*8 fracdh(nsize,maxsize),sinthi(maxsize)                          
      complex*8 ch1,ch2,dh1,dh2,xxc(maxsize),rosprp,ropprp                     
      complex*8 losprp,lopprp,lospar,loppar                                    
      complex*8 ch3,ch4,dh3,dh4,roppar,rospar                                  
      complex*8 rhsprp(nm1),rhpprp(nm1),rhspar(nm1),rhppar(nm1)                
      complex*8 lhsprp(nm1),lhpprp(nm1),lhspar(nm1),lhppar(nm1)                
      complex*8 myexp                                                                          
      
      
      include 'c2.inc'

      data job /01/                                                             
      n1=n-1                                                                    
      call angles(theta)                                                        
c      if(thetax.gt.0.e0) go to 430                                               
      ro=0.e0                                                                   
      lo=0.e0                                                                   
										
      do 2620 j=1,n1                                                            
      rh(j)=0.e0                                                                
      lh(j)=0.e0                                                                
 2620 continue                                                                  
c      go to 101                                                                
  430 continue                                                                  
  427 sintx=sin(thetax)                                                        
      costx=cos(thetax)                                                        
c      write(*,1110)thetax,sintx,costx                                          
 1110 format('thetax=',e13.5,' sin=',e13.5,' cos=',e13.5)                       
c     definition of u(j,i)                                                      
c     j = beam index (j = 1 corresponds to 0, 0, 0)                             
c     i = cartesian components                                                  
      u(1,1)=+costx                                                             
      u(1,2)=sintx                                                              
      u(1,3)=0.e0                                                               
      do 26 j=2,n                                                               
      u(j,1)=+costx/alamda+bhx(j)                                               
      u(j,2)=sintx/alamda+bhy(j)                                                
      u(j,3)=bhz(j)                                                             
      fac=sqrt(u(j,1)**2+u(j,2)**2+u(j,3)**2)                                  
      u(j,1)=u(j,1)/fac                                                         
      u(j,2)=u(j,2)/fac                                                         
      u(j,3)=u(j,3)/fac                                                         
c      write(*,1112)bhx(j),bhy(j),bhz(j)                                        
c      write(*,1111)j,fac,u(j,1),u(j,2),u(j,3)                                  
   26 continue                                                                  
 1112 format('bhx, bhy, bhz =',3(e13.5,2x))                                     
 1111 format('j=',i2,' fac=',e13.5,'u(j,i)=',3(e13.5,2x))                        
c     definition of sigma(j,i)                                                  
c     sigma(h) = u(h) x u(0)                                                    
      sigma(1,1)=0.e0                                                           
      sigma(1,2)=0.e0                                                           
      sigma(1,3)=1.e0                                                           
      do 27 j=2,n                                                               
	 sigma(j,1)=-u(j,3)*u(1,2)                                                 
	 sigma(j,2)=u(j,3)*u(1,1)                                                  
	 sigma(j,3)=u(j,1)*u(1,2)-u(j,2)*u(1,1)                                    
	 fac=sqrt(sigma(j,1)**2+sigma(j,2)**2+sigma(j,3)**2)                      
	 sigma(j,1)=sigma(j,1)/fac                                                 
	 sigma(j,2)=sigma(j,2)/fac                                                 
	 sigma(j,3)=sigma(j,3)/fac                                                 
   27 continue                                                                  
c     definition of p(j,i)                                                      
c     p(h) = u(h) x sigma(h)                                                    
      do 35 j=1,n                                                               
	 p(j,1)=-sigma(j,2)*u(j,3)+sigma(j,3)*u(j,2)                               
	 p(j,2)=-sigma(j,3)*u(j,1)+sigma(j,1)*u(j,3)                               
	 p(j,3)=-sigma(j,1)*u(j,2)+sigma(j,2)*u(j,1)                               
   35 continue                                                                  
      if(nprint.eq.0) go to 88                                                  
      do 28 j=1,n                                                               
	 print 60                                                                  
	 print 29,h(j),k(j),l(j),u(j,1),sigma(j,1),p(j,1)                          
	 write(1,29)h(j),k(j),l(j),u(j,1),sigma(j,1),p(j,1)                        
	 print 30,u(j,2),sigma(j,2),p(j,2)                                         
	 write(1,30)u(j,2),sigma(j,2),p(j,2)                                       
	 print 30,u(j,3),sigma(j,3),p(j,3)                                         
	 write(1,30)u(j,3),sigma(j,3),p(j,3)                                       
   28 continue                                                                  
   88 continue                                                                  
   29 format(i5,2i2,5x,2Hu=,e13.5,5x,6Hsigma=,e13.5,5x,3Hpi=,e13.5)                
   30 format(16x,e13.5,11x,e13.5,8x,e13.5)                                      
c     off diagonal terms of submatrices a,b,c,d                                 
      do 31 i=1,n                                                               
      do 31 j=1,n                                                               
	 if(i.eq.j)go to 31                                                        
      a(i,j)=psi(i,j)*(p(i,1)*p(j,1)+p(i,2)*p(j,2)+p(i,3)*p(j,3) )              
      b(i,j)=psi(i,j)*(p(i,1)*sigma(j,1)+p(i,2)*sigma(j,2)+p(i,3)*              
     1 sigma(j,3))                                                              
      c(i,j)=psi(i,j)*(sigma(i,1)*p(j,1)+sigma(i,2)*p(j,2)+sigma(i,3)*          
     1 p(j,3))                                                                  
      d(i,j)=psi(i,j)*(sigma(i,1)*sigma(j,1)+sigma(i,2)*sigma(j,2)+             
     1 sigma(i,3)*sigma(j,3))                                                   
   31 continue                                                                  
c     diagonal terms of submatrices a, b, c, d                                  
      do 32 i=1,n                                                               
      a(i,i)=psio-1.e0                                                          
      b(i,i)=(0.e0,0.e0)                                                        
      c(i,i)=(0.e0,0.e0)                                                        
      d(i,i)=psio-1.e0                                                          
   32 continue                                                                  
c     a matrix (2n*2n)                                                          
      do 33 i=1,n                                                               
      do 33 j=1,n                                                               
      ac(i,j)=a(i,j)                                                            
      ac(i,j+n)=b(i,j)                                                          
      ac(i+n,j)=c(i,j)                                                          
      ac(i+n,j+n)=d(i,j)                                                        
   33 continue                                                                  
      n2=2*n                                                                    
      n4=4*n                                                                    
      if(nprint.eq.0)go to 37                                                   
      print 60                                                                  
      do 36 i=1,n2                                                              
      do 36 j=1,n2                                                              
      print 34,i,j,ac(i,j)                                                      
      write(1,34)i,j,ac(i,j)                                                    
 36   continue                                                                  
   34 format(3h i=,i3,3x,2hj=,i3,5x,8hac(i,j)=,2e18.10)                            
   37 continue                                                                  
c     inversion                                                                 
      do 156 i=1,n2                                                             
      do 155 j=1,n2                                                             
	 xc(i,j)=ac(i,j)                                                        
  155 continue
  156 continue
c
c                                                                               
c     Matrix inversion with linpack routines....                                
c                                                                               
c      write(*,*)'Calling zgefa, zgedi...'                                      
      call cgefa(xc,nmax2,n2,ipvt,info)                                         
      call cgedi(xc,nmax2,n2,ipvt,det,work,job)                                 
c                                                                               
c     << call lineq4(ac,bc,xc,16,n2,n2,isinv) >> (CDC routine)                  
c     ac inverted is now in xc                                                  
c                                                                               
c     definition of aac(4n*4n) = q matrix in my notes                           
c     upper left (2n*2n)                                                        
c                                                                               
      do 158 i=1,n2                                                             
      do 158 j=1,n2                                                             
  158 aac(i,j)=(0.e0,0.e0)                                                      
      do 159 j=2,n                                                              
      aac(j,j)=-2.e0*bhy(j)*alamda                                              
  159 aac(j+n,j+n)=aac(j,j)                                                     
c     upper right (2n*2n)                                                       
      do 160 i=1,n2                                                             
      do 160 j=1,n2                                                             
  160 aac(i,j+n2)=-xc(i,j)                                                      
      do 161 j=1,n                                                              
      khtx(j)=costx/alamda+bhx(j)                                               
      khtz(j)=bhz(j)                                                            
      aac(j,n2+j)=aac(j,n2+j)-(bhy(j)**2+khtx(j)**2+khtz(j)**2)*alamda**        
     1 2                                                                        
  161 aac(j+n,n2+n+j)=aac(j,n2+j)                                               
c     lower left (2n*2n)                                                        
      do 162 i=1,n2                                                             
      do 163 j=1,n2                                                             
  163 aac(n2+i,j)=(0.e0,0.e0)                                                   
  162 aac(n2+i,i)=(1.e0,0.e0)                                                   
c     lower right (2n*2n)                                                       
      do 164 i=1,n2                                                             
      do 164 j=1,n2                                                             
  164 aac(i+n2,j+n2)=(0.e0,0.e0)                                                
c     the eigenvalues are the normal components of betao*alamda                 
      ncal=1                                                                    
      if(nprint.eq.0)go to 165                                                  
      print 60                                                                  
      write(1,60)                                                               
      do 166 i=1,n4                                                             
      do 166 j=1,n4                                                             
      print 167,i,j,aac(i,j)                                                    
      write(1,167)i,j,aac(i,j)                                                  
 166  continue                                                                  
 165  continue                                                                  
  167 format(3h i=,i3,2x,2hj=,i3,5x,9haac(i,j)=,2e18.10)                           
c                                                                               
c                                                                               
      do 780 i=1,n4                                                             
	 do 779 j=1,n4                                                          
	    aacr(i,j)=real(aac(i,j))                                           
	    aaci(i,j)=aimag(aac(i,j))                                           
 779        continue                                                            
 780     continue                                                               
c                                                                               
c                                                                               
	 matz=1                                                                 
	 ierr=0                                                                 
c         write(*,*)'     Calling cg...'                                        
	 call cg(maxsize,n4,aacr,aaci,lambdr,lambdi,matz,                       
     *    zr,zii,fv1,fv2,fv3,ierr)                                              
c         write(*,*)'     Return from cg...'                                    
	 if(ierr.ne.0)then                                                      
	    write(*,*)'    An error was found in cg...'                         
	    write(*,9876) ierr                                                  
 9876       format('     Error number: ',i5)                                    
	 endif                                                                  
c                                                                               
c                                                                               
	 do 782 i=1,n4                                                          
	    lambda(i)=cmplx(lambdr(i),lambdi(i))                               
c            write(*,1102)i,lambda(i)                                           
 782     continue                                                               
 1102    format(i2,' lambda=',2e13.5)                                           
										
c                                                                               
c    << 165 call cgeig(32,n4,aac,ncal,lambda) >> (CDC routine)                  
c     computation of eigenvectors                                               
										
      do 169 i=1,n4                                                             
      do 169 j=1,n2                                                             
      ac(j,i)=cmplx(zr(j,i),zii(j,i))                                          
 169  continue                                                                  
 916  format('j=',i2,'  i=',i2,'  ac(j,i)=',2e13.5)                             
c     the eigenvectors are now stored in ac(j,i)                                
c     j = beam index                                                            
c     i = eigenvalue                                                            
c     normalization of eigenvectors                                             
										
c      write(*,*)'Calculate normalizing factor, qqq....'                        
      do 85 i=1,n4                                                              
      qqq=(0.e0,0.e0)                                                           
      do 86 j=1,n2                                                              
										
   86 qqq=qqq+cabs(ac(j,i))**2                                                 
      qqq=csqrt(qqq)                                                           
c      write(*,915)i,qqq                                                        
 915  format('i=',i2,2x,'qqq=',2e13.5)                                          
      do 87 j=1,n2                                                              
   87 ac(j,i)=ac(j,i)/qqq                                                       
   85 continue                                                                  
										
      if(nprint.eq.0)go to 3800                                                 
      write(*,*)'List lambda(i), ac(ij)...'                                     
      print 60                                                                  
      write(1,60)                                                               
      do 38 i=1,n4                                                              
      print 41,lambda(i),ac(1,i),i                                              
      write(1,41)lambda(i),ac(1,i),1                                            
   41 format(8h deltao=,2e18.10,5x/,7hac( 1)=,2e13.5,' i=',i2)                    
      do 42 j=2,n2                                                              
      print 43,j,ac(j,i)                                                        
      write(1,43)j,ac(j,i)                                                      
 42   continue                                                                  
   43 format(3hac(,i2,i2,2h)=,2e13.5)                                                
   38 continue                                                                  
 3800 continue                                                                  
      do 80 i=1,n4                                                              
      esp(i)=-2.e0*pi*(0.e0,1.e0)*lambda(i)*to/(alamda*g0)                          
c      write(*,*)i,' esp(i)=',esp(i)                                            
 80   continue                                                                  
      if(nprint.eq.0)go to 81                                                   
      print 60                                                                  
      write(1,60)                                                               
      write(*,*)'List i, esp(i)...'                                             
      do 82 i=1,n4                                                              
      print 83,i,esp(i)                                                         
      write(1,83)i,esp(i)                                                       
 82   continue                                                                  
   83 format(i3,5x,7hesp(i)=,2e13.5)                                             
   81 continue                                                                  
										
c     ordering...                                                               
      do 90 i=1,n4                                                              
	rel(i)=real(esp(i))                                                      
 90   continue                                                                  
      n41 = n4-1 
      do 91 i=1,n41                                                             
      in1=i+1                                                                   
      do 91 j=in1,n4                                                            
	if(rel(i).le.rel(j))go to 91                                              
		temp=rel(i)                                                               
		rel(i)=rel(j)                                                             
		rel(j)=temp                                                               
		do 92 ji=1,n2                                                             
			qqq=ac(ji,i)                                                              
			ac(ji,i)=ac(ji,j)                                                         
			ac(ji,j)=qqq                                                              
 92             continue                                                                  
		qqq=esp(i)                                                                
		esp(i)=esp(j)                                                             
		esp(j)=qqq                                                                
		qqq=lambda(i)                                                             
		lambda(i)=lambda(j)                                                       
		lambda(j)=qqq                                                             
 91   continue                                                                  
      n2p1=n2+1                                                                 
      ik=n4                                                                     
      do 93 i=n2p1,n4                                                           
c********* TEST ************************************                            
c      if((rel(i)-rel(i-1)). gt. 15.e0) goto 94
c      if(rel(i) .gt. 15.e0) goto 94
c      if((rel(i)-rel(i-1))-10.e0)93,93,94
   93 continue                                                                  
      go to 113                                                                 
   94 ik=i-1                                                                    
  113 continue                                                                  
c      do i=1,ik                                                                           
c        write(6,*) 'i : ',i,' rel(i) : ',rel(i)
c      end do
										
										
c      do 9950 i=1,ik                                                           
c         write(6,*)i,' esp(i)= ',esp(i)                                        
c 9950 continue                                                                 
										
cthr      write(6,*)'ik=',ik                                                        
										
c     ik is the number of non zero psi's, that is, the order of                 
c     the final system                                                          
c     definition of eo/eh = eoeh(j,i)                                           
c     j = beam index (j = 1 corresponds to 000)                                 
c     i = eigenvalue                                                            
										
      do 54 j=1,n                                                               
      do 54 i=1,ik                                                              
      
      qqq=lambda(i)**2+                                                         
     .  2.e0*lambda(i)*bhy(j)*alamda+
     .  cmplx((khtx(j)**2+khtz(j)**2+bhy(j)**2)*alamda**2,0.e0)                                                   
   
   54 eoeh(j,i)=1.e0/qqq                                                        
      if (nprint.eq.0) go to 55                                                 
      print 60                                                                  
      write(1,60)                                                               
      do 56 j=1,n                                                               
      do 56 i=1,ik                                                              
      print 57,j,i,eoeh(j,i)                                                    
      write(1,57)j,i,eoeh(j,i)                                                  
 56   continue                                                                  
   57 format(3h j=,i3,3x,2hi=,i3,5x,10heoeh(j,i)=,2e18.10)                         
   55 continue                                                                  
										
c     definition of ui(1,ix,i)                                                  
c     1 = beam (000)                                                            
c     ix = 1, 2, 3, cartesian coordinates                                       
c     i = eigenvalue                                                            
      do 45 i=1,ik                                                              
      q1=cmplx(+costx,0.e0)                                                                 
      q2=lambda(i)                                                              
      q3=(0.e0,0.e0)                                                            
cthr
      fac=sqrt(cabs(q1)**2+cabs(q2)**2+cabs(q3)**2)                         
      ui(1,1,i)=q1/fac                                                          
      ui(1,2,i)=q2/fac                                                          
      ui(1,3,i)=q3/fac                                                          
   45 continue                                                                  
										
c     definition of ui(j,ix,i)                                                  
c     j = beam index                                                            
c     ix = cartesian component                                                  
c     i = eigenvalue                                                            
      do 46 i=1,ik                                                              
      do 46 j=2,n                                                               
      q1=cmplx(costx/alamda+bhx(j),0.e0)                                                    
      q2=lambda(i)/alamda+bhy(j)                                                
      q3=bhz(j)                                                                 
      fac=sqrt(cabs(q1)**2+cabs(q2)**2+cabs(q3)**2)                         
      ui(j,1,i)=q1/fac                                                          
      ui(j,2,i)=q2/fac                                                          
      ui(j,3,i)=q3/fac                                                          
   46 continue                                                                  
										
      if(nprint.eq.0)go to 69                                                   
      print 60                                                                  
      write(1,60)                                                               
      do 68 j=1,n                                                               
      print 70                                                                  
      write(1,70)                                                               
   70 format(127x)                                                              
      do 68 i=1,ik                                                              
      do 68 ix=1,3                                                              
      print 71,j,ix,i,ui(j,ix,i)                                                
      write(1,71)j,ix,i,ui(j,ix,i)                                              
   71 format(3h j=,i3,5x,3hix=,i3,5x,2hi=,i3,5x,11hui(j,ix,i)=,2e18.10)             
   68 continue                                                                  
   69 continue                                                                  
										
c     definition of gamma(h)                                                    
      gamma(1)=0.e0                                                             
      do 47 j=2,n                                                               
      if(khtx(j).gt.0.e0)gamma(j)=-atan(khtz(j)/khtx(j))                         
      if(khtx(j).lt.0.e0)gamma(j)=-atan(khtz(j)/khtx(j)) +pi                     
      if(khtx(j).eq.0.e0)gamma(j)=-(pi/2.e0)*(khtz(j)/abs(khtz(j)))              
      cgamma(j)=cos(gamma(j))                                                  
      sgamma(j)=sin(gamma(j))                                                  
   47 continue                                                                  
										
      if(nprint.eq.0)go to 7200                                                 
      print 60                                                                  
      write(1,60)                                                               
      do 72 j=1,n                                                               
      print 73,j,gamma(j)                                                       
      write(1,73)j,gamma(j)                                                     
   73 format(3h j=,i3,5x,9hgamma(j)=,f10.5)                                       
   72 continue                                                                  
 7200 continue                                                                  
										
c     change of reference                                                       
      do 48 i=1,ik                                                              
      do 48 j=2,n                                                               
      uii(j,1,i)=ui(j,1,i)*cgamma(j)-ui(j,3,i)*sgamma(j)                        
      uii(j,2,i)=ui(j,2,i)                                                      
      uii(j,3,i)=ui(j,1,i)*sgamma(j)+ui(j,3,i)*cgamma(j)                        
   48 continue                                                                  
										
      if(nprint.eq.0)go to 170                                                  
      print 60                                                                  
      write(1,60)                                                               
      do 77 j=2,n                                                               
      print 70                                                                  
      write(1,70)                                                               
      do 77 i=1,ik                                                              
      do 77 ix=1,3                                                              
      print 78,j,ix,i,uii(j,ix,i)                                               
      write(1,78)j,ix,i,uii(j,ix,i)                                             
   78 format(3h j=,i3,5x,3hix=,i3,5x,2hi=,i3,5x,12huii(j,ix,i)=,2e18.10)            
   77 continue                                                                  
  170 continue                                                                  
										
c     definition of sigmi(j,ix,i)                                               
c     sigmi(h) = uii(h) x uii(0)                                                
      do 49 i=1,ik                                                              
      do 49 j=2,n                                                               
      uii(1,1,i)=ui(1,1,i)*cgamma(j)-ui(1,3,i)*sgamma(j)                        
      uii(1,2,i)=ui(1,2,i)                                                      
      uii(1,3,i)=ui(1,1,i)*sgamma(j)+ui(1,3,i)*cgamma(j)                        
ccthr      
      q1=uii(j,2,i)*uii(1,3,i)-uii(j,3,i)*uii(1,2,i)                            
      q2=-uii(j,1,i)*uii(1,3,i)+uii(j,3,i)*uii(1,1,i)                           
      q3=uii(j,1,i)*uii(1,2,i)-uii(j,2,i)*uii(1,1,i)                            
      fac=sqrt(cabs(q1)**2+cabs(q2)**2+cabs(q3)**2)                         
      sigmi(j,1,i)=q1/fac                                                       
      sigmi(j,2,i)=q2/fac                                                       
      sigmi(j,3,i)=q3/fac                                                       
   49 continue                                                                  
										
      if(nprint.eq.0)go to 171                                                  
      print 60                                                                  
      write(1,60)                                                               
      do 79 j=2,n                                                               
      print 70                                                                  
      write(1,70)                                                               
      do 79 i=1,ik                                                              
      do 79 ix=1,3                                                              
      print 84,j,ix,i,sigmi(j,ix,i)                                             
      write(1,84)j,ix,i,sigmi(j,ix,i)                                           
   84 format(3h j=,i3,5x,3hix=,i3,5x,2hi=,i3,5x,14hsigmi(j,ix,i)=,
     .       2e18.10)          
   79 continue                                                                  
  171 continue                                                                  
										
c     definition of pii(j,ix,i)                                                 
c     pii(h)=uii(h) x sigmi(h)                                                  
      do 50 i=1,ik                                                              
      do 50 j=2,n                                                               
      pii(j,1,i)=-sigmi(j,2,i)*uii(j,3,i)+sigmi(j,3,i)*uii(j,2,i)               
      pii(j,2,i)=sigmi(j,1,i)*uii(j,3,i)-sigmi(j,3,i)*uii(j,1,i)                
      pii(j,3,i)=-sigmi(j,1,i)*uii(j,2,i)+sigmi(j,2,i)*uii(j,1,i)               
   50 continue                                                                  
										
      if(nprint.eq.0)go to 172                                                  
      print 60                                                                  
      write(1,60)                                                               
      do 131 j=2,n                                                              
      print 70                                                                  
      write(1,70)                                                               
      do 131 i=1,ik                                                             
      do 131 ix=1,3                                                             
      print 132,j,ix,i,pii(j,ix,i)                                              
      write(1,132)j,ix,i,pii(j,ix,i)                                            
  132 format(3h j=,i3,5x,3hix=,i3,5x,2hi=,i3,5x,12hpii(j,ix,i)=,
     .       2e18.10)            
  131 continue                                                                  
  172 continue                                                                  
										
c     definition of pihbarxprime=pihbar(j)                                      
      do 51 j=2,n                                                               
      qqq=1.e0/alamda**2-khtx(j)**2-khtz(j)**2                                    
      argn=csqrt(qqq)                                                          
      pihbar(j)=-argn*alamda                                                    
   51 continue                                                                  
										
      if(nprint.eq.0)go to 173                                                  
      print 60                                                                  
      write(1,60)                                                               
      do 133 j=2,n                                                              
      print 134,j,pihbar(j)                                                     
      write(1,134)j,pihbar(j)                                                   
  134 format(3h j=,i3,5x,13hpihbarxprime=,2e18.10)                                
  133 continue                                                                  
  173 continue                                                                  
										
c     definition of fracdh(j,i)=1/(1+delta(h,i))                                
      do 52 i=1,ik                                                              
   60 format(1h*)                                                               
      do 52 j=1,n                                                               
      fracdh(j,i)=csqrt(eoeh(j,i))                                             
   52 continue                                                                  
										
      if(nprint.eq.0)go to 174                                                  
      print 60                                                                  
      write(1,60)                                                               
      do 135 j=1,n                                                              
      print 70                                                                  
      write(1,70)                                                               
      do 135 i=1,ik                                                             
      print 136,j,i,fracdh(j,i)                                                 
      write(1,136)j,i,fracdh(j,i)                                               
  136 format(3h j=,i3,5x,2hi=,i3,5x,12hfracdh(j,i)=,2e18.10)                       
  135 continue                                                                  
  174 continue                                                                  
										
      sintx2=sintx**2                                                           
c     definition of sinthi(i)=sin(thetai)/sintx                                 
      do 53 i=1,ik                                                              
      qqq=lambda(i)**2+1.e0-sintx2                                              
      sinthi(i)=lambda(i)/(sintx*csqrt(qqq))                                   
   53 continue                                                                  
										
      if(nprint.eq.0)go to 175                                                  
      print 60                                                                  
      write(1,60)                                                               
      do 137 i=1,ik                                                             
      print 138,i,sinthi(i)                                                     
      write(1,138)i,sinthi(i)                                                   
  138 format(3h i=,i3,5x,10hsinthi(i)=,2e18.10)                                   
  137 continue                                                                  
  175 continue                                                                  
										
c     determination of the psi's                                                
c     perpendicular polarisation                                                
      do 58 i=1,ik                                                              
      aac(1,i)=ac(n+1,i)*(eoeh(1,i)+sinthi(i)*fracdh(1,i))                      
      do 58 j=2,n                                                               
      aac(2,i)=ac(1,i)*(fracdh(1,i)+sinthi(i)*eoeh(1,i))                        
      ch1=eoeh(j,i)*pii(j,3,i)+sigmi(j,1,i)*fracdh(j,i)/pihbar(j)               
      ch2=eoeh(j,i)*sigmi(j,3,i)-pii(j,1,i)*fracdh(j,i)/pihbar(j)               
      ch3=fracdh(j,i)*sigmi(j,3,i)-pii(j,1,i)*eoeh(j,i)/pihbar(j)               
      ch4=-fracdh(j,i)*pii(j,3,i)-sigmi(j,1,i)*eoeh(j,i)/pihbar(j)              
      aac(2*j-1,i)=ac(j,i)*ch1+ac(j+n,i)*ch2                                    
      aac(2*j,i)=ac(j,i)*ch3+ac(j+n,i)*ch4                                      
   58 continue                                                                  
										
c                                                                               
c      if(ik.eq.n2)go to 59                                                     
c                                                                               
										
c     calculate dynamical absorption coefficients for all                       
c     ik wavefields.                                                            
      do 191 i=1,ik                                                             
	rel(i)=real(esp(i))
	if(rel(i).lt.-15.e0) go to 192
	cxp(i)=myexp(esp(i))
c      write(*,*)'cxp(',i,')=',cxp(i)                                           
      go to 191
  192 cxp(i)=(0.e0,0.e0)
  191 continue                                                                  
										
      ikm1=ik-1                                                                 
      do 98 i=1,ikm1                                                            
      esp(i)=esp(i)-esp(ik)                                                     
      rel(i)=real(esp(i))                                                      
      if(rel(i).lt.-15.e0) go to 193
      coi(i)=myexp(esp(i))                                                      
      go to 98
  193 coi(i)=(0.e0,0.e0)
   98 continue                                                                  
      coi(ik)=(1.e0,0.e0)                                                       
										
      if(ik.eq.n2)go to 59                                                      
										
      do 63 i=1,ik                                                              
      aac(n2+1,i)=ac(n+1,i)*coi(i)*(eoeh(1,i)-sinthi(i)*fracdh(1,i))            
   63 aac(n2+2,i)=ac(1,i)*coi(i)*(fracdh(1,i)-sinthi(i)*eoeh(1,i))              
										
      if(ik.eq.n2+2)go to 59                                                    
										
      do 64 j=2,n                                                               
      do 64 i=1,ik                                                              
										
      if(n2+2*(j-1).gt.ik)go to 59                                              
										
      dh1=eoeh(j,i)*pii(j,3,i)-sigmi(j,1,i)*fracdh(j,i)/pihbar(j)               
      dh2=eoeh(j,i)*sigmi(j,3,i)+pii(j,1,i)*fracdh(j,i)/pihbar(j)               
      dh3=fracdh(j,i)*sigmi(j,3,i)+pii(j,1,i)*eoeh(j,i)/pihbar(j)               
      dh4=-fracdh(j,i)*pii(j,3,i)+sigmi(j,1,i)*eoeh(j,i)/pihbar(j)              
      aac(n2+2*j-1,i)=coi(i)*(ac(j,i)*dh1+ac(j+n,i)*dh2)                        
      aac(n2+2*j,i)=coi(i)*(ac(j,i)*dh3+ac(j+n,i)*dh4)                          
   64 continue                                                                  
   59 bbc(1)=(2.e0,0.e0)                                                        
      do 65 i=2,ik                                                              
   65 bbc(i)=(0.e0,0.e0)                                                        
      jflag=0                                                                   
   12 continue                                                                  
c                                                                               
c     Matrix inversion with Linpack routines...                                 
c     ...and calculation of vectors xxc...                                      
c                                                                               
      do 6660 i=1,maxsize
	 do 6659 j=1,ik
	    zaac(i,j)=aac(i,j)
 6659    continue                                                               
 6660 continue                                                                  
      job1=0
c      call csvdc(zaac,maxsize,ik,maxsize,s_w,e,u,maxsize,v,maxsize,
c     .            work,job1,info)
c      wmax = 0.
c      do j=1,ik
c        if(real(s_w(j)).gt.wmax)wmax=real(s_w(j))
c      end do
c      wmin=wmax*1.e-6
c      do j=1,ik
c        if(real(s_w(j)).lt.wmin)s_w(j)=(0.0,0.0)
c      end do
c
c      call svbksb(u,s_w,v,maxsize,ik,maxsize,maxsize,bbc,xxc)
c      write(*,*)'Calling zgefa, zgedi...'                                      
       call cgefa(zaac,maxsize,ik,ipvt,info)
c       if(info.ne.0) write(6,*) ' INFO : ',info
       if(info.ne.0) then 
		ik = info -1
	       goto 12
       endif
      call cgesl(zaac,maxsize,ik,ipvt,bbc,job1)
      call cgedi(zaac,maxsize,ik,ipvt,det,work,job)
c
      do 777 i=1,ik
	 xxc(i)=bbc(i)                                                          
777   continue                                                                  
										
c                                                                               
c     << call lineq4(aac,bbc,xxc,32,ik,1,isprp) >> (CDC routine)                
c                                                                               
c      if(isprp.eq.0.or.jflag.eq.1)go to 10                                     
c      do 11 i=1,7                                                              
c   11 aac(7,i)=aac(8,i)                                                        
c      jflag=1                                                                  
c      go to 12                                                                 
c     above statements only applicable when singularities occurr                
c     in the two-beam case and ik = 7.                                          
   10 continue                                                                  
										
      if(nprint.eq.0)go to 120                                                  
      print 60                                                                  
      do 20 j=1,ik                                                              
      do 20 i=1,ik                                                              
      print 21,j,i,aac(j,i)                                                     
      write(1,21)j,i,aac(j,i)                                                   
 20   continue                                                                  
   21 format(14h j,i,aac(j,i)=,2i3,2e20.6)                                       
      print 60                                                                  
      write(1,60)                                                               
      do 121 i=1,ik                                                             
      print 122,i,xxc(i)                                                        
      write(1,122)i,xxc(i)                                                      
 121  continue                                                                  
  122 format(3h i=,i3,5x,7hxxc(i)=,2e13.5)                                        
  120 continue                                                                  
										
c     xxc(i), i = 1, 2, ...ik are the psi=s                                     
c     in this version of nbeam, ro is the specular beam,                        
c     and lo, lh are the laue beams (lbeam, 12-23-85)                           
      ropprp=(0.e0,0.e0)                                                        
      rosprp=(0.e0,0.e0)                                                        
      lopprp=(0.e0,0.e0)                                                        
      losprp=(0.e0,0.e0)                                                        
c      write(*,*)'lopprp=',lopprp,' losprp=',losprp                             
										
      do 2621 j=1,n1                                                            
      rhpprp(j)=(0.e0,0.e0)                                                     
      rhsprp(j)=(0.e0,0.e0)                                                     
      lhpprp(j)=(0.e0,0.e0)                                                     
      lhsprp(j)=(0.e0,0.e0)                                                     
c      write(*,*)'j=',j,' lhsprp=',lhsprp,' lhpprp=',lhpprp                     
 2621 continue                                                                  
										
      do 66 i=1,ik                                                              
c         write(*,*)'i=',i,'cxp=',cxp(i)                                        
      ropprp=ropprp+xxc(i)*ac(1,i)*fracdh(1,i)                                  
      lopprp=lopprp+xxc(i)*ac(1,i)*fracdh(1,i)*cxp(i)                           
      rosprp=rosprp+.5d0*xxc(i)*ac(n+1,i)*(eoeh(1,i)-sinthi(i)*                 
     1 fracdh(1,i))                                                             
										
c*******test*******x2:                                                          
										
      losprp=losprp+2.e0*.5d0*xxc(i)*ac(n+1,i)*eoeh(1,i)*cxp(i)                   
c      write(*,*)'lopprp=',lopprp,' losprp=',losprp                             
										
      do 2622 j=1,n1                                                            
      jj=j+1                                                                    
      rhpprp(j)=rhpprp(j)+xxc(i)*fracdh(jj,i)*(ac(jj,i)*sigmi(jj,3,i)           
     1  -ac(n+jj,i)*pii(jj,3,i))                                                
      lhpprp(j)=lhpprp(j)+xxc(i)*fracdh(jj,i)*(ac(jj,i)*sigmi(jj,3,i)           
     1  -ac(n+jj,i)*pii(jj,3,i))*cxp(i)                                         
      rhsprp(j)=rhsprp(j)+xxc(i)*eoeh(jj,i)*(ac(jj,i)*pii(jj,3,i)               
     1  +ac(jj+n,i)*sigmi(jj,3,i))                                              
      lhsprp(j)=lhsprp(j)+xxc(i)*eoeh(jj,i)*(ac(jj,i)*pii(jj,3,i)               
     1  +ac(jj+n,i)*sigmi(jj,3,i))*cxp(i)                                       
c      write(*,*)'j=',j,' lhsprp=',lhsprp,' lhpprp=',lhpprp                     
 2622 continue                                                                  
   66 continue                                                                  
										
      if(nprint.eq.0)go to 139                                                  
      print 60                                                                  
      write(1,60)                                                               
      print 142,ropprp,rosprp                                                   
      write(1,142)ropprp,rosprp                                                 
  142 format(8h ropprp=,2e13.5,5x,7hrosprp=,2e13.5)                               
      do 2623 j=1,n1                                                            
      print 143,j+1,rhpprp(j),rhsprp(j)                                         
      write(1,143)j+1,rhpprp(j),rhsprp(j)                                       
 2623 continue                                                                  
  143 format('j=',i2,8h rhpprp=,2e13.5,5x,7hrhsprp=,2e13.5)                       
  139 continue                                                                  
										
      roprp=cabs(ropprp)**2+cabs(rosprp)**2                                   
      loprp=cabs(lopprp)**2+cabs(losprp)**2                                   
										
      do 2624 j=1,n1                                                            
      rhprp(j)=cabs(rhpprp(j))**2+cabs(rhsprp(j))**2                          
      lhprp(j)=cabs(lhpprp(j))**2+cabs(lhsprp(j))**2                          
 2624 continue                                                                  
c                                                                               
c                                                                               
c     parallel polarisation                                                     
c                                                                               
c                                                                               
      do 250 i=1,ik                                                             
      aac(1,i)=ac(1,i)*(fracdh(1,i)+sinthi(i)*eoeh(1,i))                        
      aac(2,i)=ac(n+1,i)*(eoeh(1,i)+sinthi(i)*fracdh(1,i))                      
  250 continue                                                                  
c                                                                               
c     Solution vectors with Linpack...                                          
c                                                                               
      bbc(1)=(2.e0,0.e0)                                                        
      do 767 i=2,ik                                                             
	 bbc(i)=(0.e0,0.e0)                                                     
 767  continue                                                                  
										
      do 6662 i=1,maxsize                                                       
	 do 6661 j=1,maxsize                                                    
	    zaac(i,j)=aac(i,j)                                                  
 6661    continue                                                               
 6662 continue                                                                  
										
      job1=0                                                                    
c      write(*,*)'Calling zgefa, zgesl, and zgedi...'                           
      call cgefa(zaac,maxsize,ik,ipvt,info)                                     
      call cgesl(zaac,maxsize,ik,ipvt,bbc,job1)                                 
      call cgedi(zaac,maxsize,ik,ipvt,det,work,job)                             
										
      do 770 i=1,ik                                                             
	 xxc(i)=bbc(i)                                                          
 770  continue                                                                  
										
c                                                                               
c                                                                               
c     << call lineq4(aac,bbc,xxc,32,ik,1,ispar) >> (CDC routine)                
										
      if(nprint.eq.0)go to 251                                                  
      print 60                                                                  
      write(1,60)                                                               
      do 252 i=1,ik                                                             
      print 122,i,xxc(i)                                                        
      write(1,122)i,xxc(i)                                                      
 252  continue                                                                  
  251 continue                                                                  
										
      roppar=(0.e0,0.e0)                                                        
      rospar=(0.e0,0.e0)                                                        
      loppar=(0.e0,0.e0)                                                        
      lospar=(0.e0,0.e0)                                                        
										
      do 2625 j=1,n1                                                            
      rhppar(j)=(0.e0,0.e0)                                                     
      rhspar(j)=(0.e0,0.e0)                                                     
      lhppar(j)=(0.e0,0.e0)                                                     
      lhspar(j)=(0.e0,0.e0)                                                     
 2625 continue                                                                  
										
      do 253 i=1,ik                                                             
      roppar=roppar                                                             
     *       -0.5d0*xxc(i)*ac(1,i)*(sinthi(i)*eoeh(1,i)-fracdh(1,i))            
      loppar=loppar+xxc(i)*ac(1,i)*fracdh(1,i)*cxp(i)                           
      rospar=rospar+xxc(i)*ac(1+n,i)*eoeh(1,i)                                  
      lospar=lospar+xxc(i)*ac(1+n,i)*eoeh(1,i)*cxp(i)                           
										
      do 2626 j=1,n1                                                            
      jj=j+1                                                                    
      rhppar(j)=rhppar(j)+xxc(i)*fracdh(jj,i)*(ac(jj,i)*sigmi(jj,3,i)           
     1  -ac(jj+n,i)*pii(jj,3,i))                                                
      lhppar(j)=lhppar(j)+xxc(i)*fracdh(jj,i)*(ac(jj,i)*sigmi(jj,3,i)           
     1  -ac(jj+n,i)*pii(jj,3,i))*cxp(i)                                         
      rhspar(j)=rhspar(j)+xxc(i)*eoeh(jj,i)*(ac(jj,i)*pii(jj,3,i)               
     1  +ac(jj+n,i)*sigmi(jj,3,i))                                              
      lhspar(j)=lhspar(j)+xxc(i)*eoeh(jj,i)*(ac(jj,i)*pii(jj,3,i)               
     1  +ac(jj+n,i)*sigmi(jj,3,i))*cxp(i)                                       
 2626 continue                                                                  
  253 continue                                                                  
										
      if(nprint.eq.0)go to 254                                                  
      print 60                                                                  
      write(1,60)                                                               
      print 255,roppar,rospar                                                   
      write(1,255)roppar,rospar                                                 
  255 format(8h roppar=,2e13.5,5x,7hrospar=,2e13.5)                               
										
      do 2627 j=1,n1                                                            
      print 256,j+1,rhppar(j),rhspar(j)                                         
      write(1,256)j+1,rhppar(j),rhspar(j)                                       
 2627 continue                                                                  
  256 format('j=',i2,8h rhppar=,2e13.5,5x,7hrhspar=,2e13.5)                       
  254 continue                                                                  
										
c      write(*,*)'loppar=',loppar,' lopprp=',lopprp                             
										
      ropar=cabs(roppar)**2+cabs(rospar)**2                                   
      lopar=cabs(loppar)**2+cabs(lospar)**2                                   
										
      do 2628 j=1,n1                                                            
      rhpar(j)=cabs(rhppar(j))**2+cabs(rhspar(j))**2                          
      lhpar(j)=cabs(lhppar(j))**2+cabs(lhspar(j))**2                          
 2628 continue                                                                  
c                                                                               
c                                                                               
      if(nprint.eq.0)go to 257                                                  
      print 60                                                                  
      write(1,60)                                                               
      print 258,roprp,ropar                                                     
      write(1,258)roprp,ropar                                                   
  258 format(7h roprp=,e13.5,5x,6hropar=,e13.5)                                   
      print 60                                                                  
      write(1,60)                                                               
  257 continue                                                                  
										
      ro=(kprp*roprp+kpar*ropar)/abs(bi)                                                  
      lo=(kprp*loprp+kpar*lopar)/abs(bi)                                                  
c      write(*,*)'ro=',ro,'  lo=',lo                                            
c
c      write(6,*) 'In calc ....', n1                                                                                 
      do 2729 j=1,n1                                                            
      rh(j)=(kprp*rhprp(j)+kpar*rhpar(j))/abs(bi)                                         
      lh(j)=(kprp*lhprp(j)+kpar*lhpar(j))/abs(bi)                                         
cc      write(6,*)j,' rh=',rh(j),
cc     .            ' lh=',lh(j)                                    
 2729 continue                                                                  
										
  101 continue                                                                  
c      write(*,*)'   Ending calc...'                                            
      return                                                                    
      end  
