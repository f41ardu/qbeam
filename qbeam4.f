       program qbeam4                                                           
c                                                                               
c      This version calculates all n beam intensities (15 Sep 89)               
c      This version has all real HKL (1 Sep 89)                                 
c                                                                               
c      This version incorporates some linpack routines.  13-Feb-89              
c                                                                               
c      qbeam2 incorporates an internal subroutine for matrix inversion,         
c      and an eispack routine.                                                  
c                                                                               
c      program qbeam1 is based on the CDC version of lbeam10.                   
c      Most of the CDC clutter has been removed. << SMD, 22-Jun-88 >>           
c                                                                               
c      program lbeam10 is constructed from lbeam7 by removing                   
c      the calculations of phases. (these statements have                       
c      been preserved as comments.) (11-may-87)                                 
c                                                                               
c
c ---  program qbeam4 
c      this version will be build on HP - Cluster at DESY 
c      non changes in running modules, only change for errors on 
c      this special version of f77 compiler 
c      date 1/02/92 f41ardu 
c      it run with data from steve at 3/02/93
c
      implicit none                                                 
                                                                                
      include 'param.inc'
c                                                                                
c       the size of the above plot arrays must                                  
c       always be at least ntheta+3.                                            

      real*8 x(1000),yro(1000),ylo(1000),yrh(1000),ylh(1000),
     .       yrl(1000),yll(1000),
     .       yh(8,1000),yl(8,1000),                                
     .       ro,lo,rh(nm1),lh(nm1),theta,                                        
     .       phiscan,phi1,phi2,                                 
     .       kipp,zx,zy,zz,alphax,betax,gammax,energy,salamda,
     .       lamda0,xstart,xend,xstep,mtheta
      real gikn,temp
      integer nphi,i,j,npts
      integer ichoise

      character*30 name1 !,name2,name3                                            
      character reply1,reply2,reply3,reply4,reply5,reply6,reply30                       
      character reply7,reply8,reply9,reply10 
c                                                                               
      include 'c1.inc'
      include 'c2.inc'
c
   99 format(4h bi=,e13.5,5x,15hthetax(thetab)=,e20.11)                           
  260 format(5h phi=,f10.7)                                                      
   61 format(4h eo=,e13.5,4x,7hthetao=,e16.9,5x,7hthetab=,
     .        e13.5,5x,/,              
     .        7hdtheta=,e13.5,5x,5hpsio=,2e13.5)                                         
   75 format(3h i=,i3,1x,2hh=,f7.4,5x,2hk=,f7.4,5x,                                 
     . 2hl=,f7.4,2x,/,4hbhx=,f10.6,5x,                                            
     . 4hbhy=,f10.6,5x,4hbhz=,f10.6)                                              
   24 format(7h ipvec=,i5)                                                       
   76 format(i3,i2,i5,2i2,f9.4,5x,2f10.4,5x,2e13.4)                             
   74 format(' i',2h j,4x,5hh k l,5x,4hesse,4x,5x,6x,4hfhkl,8x,                
     . 13x,5x,3hpsi)                                                            
   60 format(1h*)                                                               
  259 format(1h ,5hisinv,6h  ncal,7h  isprp,7h  ispar,                          
     . 6h  ik  ,8x/,2hro,19x,2hrh,12x,6hthetax,11x,5hnscan,
     .     14x,5htheta,         
     . 9x,5hjflag)                                                              
  113 format(2x,6hroell=,f8.4,2x,6hrophi=,f8.4,2x,6hloell=,                     
     . f8.4,                                                                    
     . 2x,6hlophi=,f8.4,2x,6hrhell=,f8.4,2x,6hrhphi=,f8.4,                      
     . 2x,6hlhell=,f8.4,2x,6hlhphi=,f8.4)                                       
  112 format(2i5,2i6,i3/,4e12.5,e17.11,i4,2x,e17.11,i4)                         
                                                                                
c      write(6,*) 'Enter filename for output data:'                             
c      read(5,2727) name1                                                       
      name1='out.dat'                                                             
      write(6,*)'filename= ',name1                                              
      open(1,file='out1.dat')
      open(30,file='inth.dat')
      open(40,file='intl.dat')                                    
      write(1,*)'filename= ',name1                                              
 2727 format(a)                                                                 
      rewind(1)                                                                 
c      write(6,*)'Enter a one-line descriptive title:'                          
c      read(5,2727)title                                                        
c      write(1,2727)title                                                       
                                                                                
      pi=4.*datan(1.0d0)                                                        
c                                                                               
c     general structures. bragg asymmetric case.                                
c     cork screw rule applies to all rotations                                  
c                                                                               
c      data astar,bstar,cstar/.18423,.18423,.18423/                              
c      data alphax,betax,gammax/90.d0,90.d0,90.d0/                               
cthr
c gitterkonstante aus der teworte arbeit nach becker      
c
c      astar = 1.d0/5.43102018d0
      write(6,*) ' temperatur [K]: '
      read(5,*) temp
      astar = 1./gikn(temp,60.*pi/180.)		
      bstar = astar
      cstar = bstar
      alphax = 90.d0
      betax =  90.d0
      gammax = 90.d0
c                                                                               
c     above two lines give reciprocal lattice constants.                        
c                                                                               
c      data ph,pk,pl/1.,1.,1./                                                   
c      data mh,mk,ml/3.,-3.,1./                                                  
c
c --- choose direction of kipp vector
c
      zx = -1.d0/dsqrt(6.d0)
      zy =  zx
      zz = -zx
c
      kipp = 0.00d0 
      write(6,*)'(kipp)',kipp               
      write(6,*)'Change any of these? (y/n)'                                    
      read(5,*)reply30                                                           
      if(reply30.eq.'y')then                                                     
         write(6,*)'Enter kipp:'                                                   
         read(5,*) kipp
      endif

cthr surface normal 
      ph = 4.d0+kipp*zx
      pk = 2.d0+kipp*zy
      pl = 2.d0+kipp*zz
c reference vektor, what does that mean ?? dont changed yet !
c 6,1,5 fuer einfallsrichtung antipar 2,1,1
c antipar -4,0,-4
c   =     -4,-4,0
c 
      write(6,*) 'Strahlteiler 1, Umlenk1 2, Umlenk2 3 '
      read(5,*) ichoise
c
      if(ichoise.eq.1) then 
      mh =  4.d0
      mk =  1.d0
      ml =  3.d0
c      
      h(1) =  0.d0
      k(1) =  0.d0 
      l(1) =  0.d0 
c     
      h(2) =  4.d0
      k(2) =  4.d0
      l(2) =  0.d0
c
      h(3) =  4.d0
      k(3) =  0.d0
      l(3) =  4.d0
c
      endif
c      
      if(ichoise.eq.2) then 
      mh =  4.d0 
      mk = -1.d0 
      ml =  5.d0 
c      
      h(1) =  0.d0
      k(1) =  0.d0 
      l(1) =  0.d0 
c     
      h(2) =  4.d0
      k(2) =  0.d0
      l(2) =  4.d0
c
      h(3) =  0.d0
      k(3) =  -4.d0
      l(3) =  4.d0 
c
      endif
c
      if(ichoise.eq.3) then 
      mh =  4. 
      mk =  5. 
      ml =  -1. 
c      
      h(1) =  0.d0
      k(1) =  0.d0 
      l(1) =  0.d0 
c     
      h(2) =  4.d0
      k(2) =  4.d0
      l(2) =  0.d0
c
      h(3) =  0.d0 
      k(3) =  4.d0 
      l(3) =  -4.d0  
c
      endif
c

c 
c --- choose a new reference vektor 
c       
c                                                                               
c     ph, pk, pl are miller indices of surface, not necessarily                 
c     equal to h(2),k(2),l(2). when they are different, we have                 
c     the non symmetric bragg case.                                             
c                                                                               
      write(6,2730)astar,bstar,cstar,alphax,betax,gammax                        
      write(6,2731)ph,pk,pl,mh,mk,ml                                            
 2730 format('(astar,bstar,cstar) = ',3(f9.5,1x),/                               
     & '(alphax,betax,gammax) = ',3(f9.4,1x))                                    
 2731 format('Surface indices (ph,pk,pl) = ',3(f5.2,1x),/                        
     & 'Phi reference direction (mh,mk,ml) = ',3(f5.2,1x))                       
      write(6,*)'Change any of these? (y/n)'                                    
      read(5,*)reply1                                                           
      if(reply1.eq.'y')then                                                     
 2732    write(6,*)'Enter astar, bstar, cstar:'                                 
         read(5,*)astar,bstar,cstar                                             
         write(6,*)'Enter alphax, betax, gammax:'                               
         read(5,*)alphax,betax,gammax                                           
         write(6,*)'Enter ph, pk, pl:'                                          
         read(5,*) ph, pk, pl                                                   
         write(6,*)'Enter mh, mk, ml:'                                          
         read(5,*) mh, mk, ml                                                   
         write(6,2730)astar,bstar,cstar,alphax,betax,gammax                     
         write(6,2731)ph,pk,pl,mh,mk,ml                                         
         write(6,*)'Change any of these again? (y/n)'                           
         read(5,*)reply2                                                        
         if(reply2.eq.'y') go to 2732                                           
      endif                                                                     
c                                                                               
c      data n,thetx,alamda                                                       
c     1/2,0.0,1.542/                                                             
cthr
      n = 3
cthr      alamda = 1.6629035279d0
       alamda = 1.6628925d0
c       alamda = 1.66276
c      alamda = 1.66289d0
c      data nphi,(phio(i),i=1,10  )                                              
c     1/1,10*0.61/
      nphi = 1
c      do i=1,10
         phio(1) = 0.d0
c      end do
      phi=phio(1)                                                               
                                                                                
      write(6,2740)n,alamda                                                     
      write(6,2741)nphi,(i,phio(i),i=1,nphi)                                    
 2740 format('Number of beams n = ',i2,';',5x,                                  
     & 'Wavelength alamda = ',f10.6)                                            
 2741 format('Number of scans, at different phi values (nphi) = ',i2,           
     & /('phio(',i2,')=',f15.10))                                               
                                                                                
      write(6,*)'Change any of these? (y/n)'                                    
      read(5,*)reply3                                                           
      if(reply3.eq.'y')then                                                     
         write(6,*)'(Maximum value of n is 8.)'                                 
 2743    write(6,*)'Enter n:'                                                   
         read(5,*)n                                                             
                                                                                
         write(6,*)'Do you want to enter energy instead of wavelength?          
     & (y or n)'                                                                
         read(5,*)reply3                                                        
         if(reply3.eq.'y')then                                                  
            write(6,*)'Enter energy (keV):'                                     
            read(5,*)energy                                                     
            alamda=12.394/energy                                                
         elseif(reply3.ne.'y')then                                              
            write(6,*)'Enter wavelength alamda (in Angstrom):'                  
            read(5,*)alamda                                                     
         endif                                                                  
                                                                                
         write(6,*)'Enter number of different phio values (nphi):'              
         read(5,*)nphi                                                          
         write(6,*)'Enter nphi values of phio (in radians):'                    
         read(5,2742)(phio(i),i=1,nphi)                                         
 2742    format(f15.10)                                                         
         write(6,2740)n,alamda                                                  
         write(6,2741)nphi,(i,phio(i),i=1,nphi)                                 
         write(6,*)'Change any of these again? (y/n)'                           
         read(5,*)reply4                                                        
         if(reply4.eq.'y')go to 2743                                            
      end if                                                                    
                                                                                
       nprint = 0                                                            
c      write(6,*)'Do you want all of the intermediate calculations'             
c      write(6,*)'printed (for debugging purposes)? (y/n)'                      
c      read(5,*)reply1                                                          
c      if(reply1.eq.'y')nprint=1                                                
c      write(6,2750)nprint                                                      
 2750 format('nprint =',i2)                                                      
                                                                                
       ifac = 5
       ntheta = 50                                                  
      write(6,2751)ifac,ntheta                                                  
 2751 format('Scan width in units of Darwin widths, ifac =',i4,/                
     & 'Number of data points per scan, ntheta =',i4)                           
      write(6,*)'Change either of these? (y/n)'                                 
      read(5,*)reply5                                                           
      if(reply5.eq.'y')then                                                     
 2752    write(6,*)'Enter ifac, ntheta:'                                        
         read(5,*)ifac, ntheta                                                  
         write(6,2751)ifac,ntheta                                               
         write(6,*)'Change again? (y/n)'                                        
         read(5,*)reply6                                                        
         if(reply6.eq.'y')go to 2752                                            
      endif                                                                     
                                                                                
c      data h(1),k(1),l(1)/0.,0.,0./                                             
c      data h(2),k(2),l(2)/2.,-2.,0./                                            
c      data h(3),k(3),l(3)/1.,1.,1./                                             
c      data h(4),k(4),l(4)/2.,0.,2./                                             
c      data h(5),k(5),l(5)/-1.,1.,1./                                            
c      data h(6),k(6),l(6)/1.,-1.,1./                                            
c      data h(7),k(7),l(7)/-1.,-1.,1./                                           
c      h(3) =  -4.d0
c      k(3) = 4.d0
c      l(3) =  0.d0
c
c      h(5) = -.1
c      k(5) = 1.
c      l(5) = 1.
c
c      h(6) = 1.
c      k(6) = -1.
c      l(7) = 1.
c
c      h(7) = -1.
c      k(7) = -1.
c      l(7) = 1.
c
      write(6,2760)(i,h(i),k(i),l(i),i=1,n)                                     
 2760 format(('hkl(',i2,')=',3(f7.4,1x)))                                        
      write(6,*)'Change these? (y/n)'                                           
      read(5,*)reply7                                                           
      if(reply7.eq.'y')then                                                     
 2763    continue                                                               
         do 2761 i=1,n                                                          
            write(6,2762)i                                                      
            read(5,*)h(i),k(i),l(i)                                             
 2761    continue                                                               
         write(6,2760)(i,h(i),k(i),l(i),i=1,n)                                  
         write(6,*)'Change these again? (y/n)'                                  
         read(5,*)reply8                                                        
         if(reply8.eq.'y')go to 2763                                            
      endif                                                                     
 2762 format('Enter hkl(',i2,'):')                                              
c                                                                               
c     alamda in angstroms                                                       
c     phi in radians                                                            
c     phi=0 corresponds to (mh,mk,ml) vector mostly antiparallel                
c     to k0. phi > 0 corresponds to anticlockwise rotation with respect         
c     to scattering vector.                                                     
c     mh, mk, ml are components of a reference vector, not                      
c     necessarily normal to h(2), k(2), l(2).                                   
c     note:  phi in nbeam = +beta in umweg.                                     
c     thetax in radians                                                         
c     all angles on output in radians                                           
c     vectors are always expressed in general contravariant                     
c     coordinates.                                                              
c     ao(3,3) is the fundamental tensor in covariant coordinates.               
c     e(3,3,3) is the ricci tensor in covariant coordinates.                    
c     ai(3,3) is the fundamental tensor in contravariant                        
c     coordinates.                                                              
c                                                                               
      kprp=1.0d0                                                                
      kpar=0.0d0                                                                
      to=5.d0                                                                 
                                                                                
      write(6,2770)kprp,kpar,to                                                 
 2770 format('Perpendicular polarization kprp =',f6.3,/                         
     & 'Parallel polarization kpar =',f6.3,/                                    
     & 'Crystal thickness =',f10.5,' mm')                                 
      write(6,*)'Change these? (y/n)'                                           
      read(5,*)reply9                                                           
      if(reply9.eq.'y')then                                                     
 2771    write(6,*)'Enter kprp, kpar:'                                          
         read(5,*)kprp,kpar                                                     
         write(6,*)'Enter thickness (to):'                                      
         read(5,*)to                                                            
         write(6,2770)kprp,kpar,to                                              
         write(6,*)'Change again? (y/n)'                                        
         read(5,*)reply10                                                       
         if(reply10.eq.'y')go to 2771                                           
      endif                                                                     
      to = to*1.d7
                                                                                
c                                                                               
c     kprp, kpar are the fractions of perpendicular and parallel                
c     polarizations in the incident beam. kprp+kpar=1.                          
c     to = thickness of crystal in angstroms.                                   
                                                                                
       phiscan = 0.
       phi1 = 0.01
       phi2 = .5
       npts = 101                               
                                                                                
c                                                                               
c      phiscan= 1.0 sets theta bragg angle, scans phi                           
c      phi1 = initial phi value                                                 
c      phi2 = final phi value                                                   
c      npts = number of data points in phi range                                
c                                                                               
      izscan = 0
      zz1    = .1e7
      zz2    = .1e8
      nz     = 100                                 
c                                                                               
c                                                                               
c      write(6,*)'Calling tensor...'                                            
      call tensor(alphax,betax,gammax,ao,e,ai,e0)                               
c                                                                               
c                                                                               
c                                                                               
c --- lamda scan lamda scan lamda scan lammda scan 
c
c 
c --- an other chnge for lamda scan alamda = 1.6629035279d0   
c     we will scan 31 points around 15 below, and 15 above
      lamda0 = alamda
      xstart =   alamda - alamda/200000.d0
      xend =     alamda + alamda/200000.d0
      xstep = dabs(xend - xstart)/5.d0
c      xstep = 2.638d-6 ! Das ist die nominelle Aufl”sung des 2-Ax MC
      call brtheta
ccc      do 11111 salamda=xstart,xend,xstep                                                  
      salamda = alamda
      alamda = salamda                     
      phi=phio(1)                                                            
c      write(6,*)'Calling brtheta...'                                           
cc      call brtheta                                                              
c                                                                               
c                                                                               
c     printing control statements...                                            
c                                                                               
c                                                                               
      if(iphi.eq.1)then                                                         
         print 24,ipvec                                                         
         write(1,24)ipvec                                                        
      do 25 i=2,n                                                               
         print 75,i,h(i),k(i),l(i),bhx(i),bhy(i),bhz(i)                         
         write(1,75)i,h(i),k(i),l(i),bhx(i),bhy(i),bhz(i)                       
         print 60                                                               
c         write(1,60)                                                           
         print 61,eo,thetao,thetab,dtheta,psio                                  
         write(1,61)eo,thetao,thetab,dtheta,psio                                
         print 60                                                               
c         write(1,60)                                                           
         print 260,phi                                                          
         write(1,260)phi                                                        
         print 60                                                               
c         write(1,60)                                                           
         print 99, bi,thetax                                                    
         write(1,99)bi,thetax                                                   
         print 60                                                               
c         write(1,60)                                                           
c         print 259                                                             
c         write(1,259)                                                          
         print 60                                                               
         write(1,60)                                                            
 25      continue                                                               
      endif                                                                     
 612  continue                                                                  
c                                                                               
c                                                                               
c                                                                               
c     note: theta = angle of incidence on lattice planes                        
c          thetax = angle of incidence on xtal surface.                         
c                                                                               
c 
c      theta0 = 60.d0                                                                              
       mtheta = (20./3600)*pi/180.
       do 189 nscan=1,ntheta                                                     
c         write(6,613)nscan                                                     
c 613     format('nscan=',i3)                                                   
c      theta=thetao-dtheta+(nscan-1)*2.d0*dtheta/(ntheta-1)                      
      theta=thetao-mtheta+(nscan-1)*2.d0*mtheta/(ntheta-1) 
      call calc(theta,ro,lo,rh,lh)                                              
c      write(6,*)'Returned to main program from calc...'                        
      if(kflag.eq.1)rh(1)=0.0d0                                                 
      yro(nscan)=ro                                                             
      ylo(nscan)=lo                                                             
      yrh(nscan)=rh(1)                                                          
      ylh(nscan)=lh(1)
      yrl(nscan)=rh(2)
      yll(nscan)=lh(2)
      
      yh(1,nscan) = ro
      yl(1,nscan) = lo
      do j=2,8
      yh(j,nscan) = rh(j-1)
      yl(j,nscan) = lh(j-1)
      end do                                                          
c                                                                               
c --- Ausgabe in Winkelsekunden !!
c
      
      x(nscan)=(theta-thetao)*(180.0/pi)*3600.0                                                     
      if(ichoise.eq.3) x(nscan)=-x(nscan)
cc      write(6,*) theta,thetao
      if(nscan.eq.1)then                                                        
      write(1,2730)astar,bstar,cstar,alphax,betax,gammax                        
      write(1,2731)ph,pk,pl,mh,mk,ml                                            
      write(1,2740)n,alamda                                                     
      write(1,2741)nphi,(i,phio(i),i=1,nphi)                                    
      write(1,2751)ifac,ntheta                                                  
      write(1,2760)(i,h(i),k(i),l(i),i=1,n)                                     
      write(1,2770)kprp,kpar,to                                                 
      write(6,*) ' Lambda : ', alamda 
c      write(1,928)                                                              
c      write(6,928)                                                              
      endif                                                                     
c      write(1,927)nscan,theta*180./pi,thetax*180./pi,                           
c     & ro,rh(1),rh(2),lh(1)                                                        
c      write(6,927)nscan,theta*180./pi,thetax*180./pi,                           
c     & ro,rh(1),rh(2),lh(1)                                                        
 927  format(i3,2(2x,f8.5,2x),4(2x,e10.5))                                      
 926  format(i3,' theta=',f7.4,' thetax=',f7.4,                                 
     & ' ro=',e8.3,' rh1=',e8.3,                                                 
     & ' rh2=',e8.3,' lh=',e8.3)                                                 
 928  format(' # ',4x,'theta',3x,                                               
     &3x,'thetax',3x,                                                           
     &5x,'ro',5x,                                                               
     &4x,'rh(1)',3x,                                                            
     &4x,'rh(2)',5x,                                                               
     &4x,'lh(1)',3x)                                                            
                                                                                
c      write(6,226)nscan,theta,ro,rh(1),rh(2),rh(3),rh(4)                       
c      write(6,229)rh(5),rh(6),rh(7)                                            
 229  format(3(e12.7,1x))                                                        
 226  format(i3,1x,e12.7,1x,5(e12.7,1x))                                           
ccc11111 format(2f16.8)
  26  format(e18.11,2x,7(e15.8,1x))                                              
c      print 112,isinv,ncal,isprp,ispar,ik,ro,rh,lo,lh,thetax,                  
c     1 nscan,theta,kflag                                                       
c      ncal = number of eigenvalues calculated by cgeig.                        
c      isinv,isprp,ispar = 0 when matrix is non-singular in lineq4              
c                                                                               
c                                                                               
  189 continue                                                                  
                                                                                
cc       write(6,*)'Plot the diffracted intensity? (y/n)'                          
cc       read(5,*)reply1                                                           
ccc         call ploth(ntheta,x,yh)
ccc         call plot2(ntheta,x,yl)
       call plot3d(bi,30,ntheta,x,(lamda0-salamda),yh)
       call plot3d(bi,40,ntheta,x,(lamda0-salamda),yl)
ccc 11111 continue
      close(40)
      close(30)
c lamda scan lamda scan lamda scan 
c
c                                 
  190 continue                                                                  
c                                                                               
c                                                                               
      close(1)                                                                  
  180 stop                                                                      
      end                                                                       
                                                                                
c
c --- normal theta scan 
c
      subroutine ploth(ndata,xdata,ydata1)
c
c --- datafile schreiben 
c  
      integer ndata,iunit
      real*8 xdata(ndata),ydata1(8,ndata)
      real*8 pi
100   format(2f16.8)
      open(15,file='c:\fortran\qbeam\output\int0h.dat') 
      open(16,file='c:\fortran\qbeam\output\int1h.dat')
      open(17,file='c:\fortran\qbeam\output\int2h.dat')
      open(18,file='c:\fortran\qbeam\output\int3h.dat')
      open(19,file='c:\fortran\qbeam\output\int4h.dat')
      open(20,file='c:\fortran\qbeam\output\int5h.dat')          
      open(21,file='c:\fortran\qbeam\output\int6h.dat')
      open(22,file='c:\fortran\qbeam\output\int7h.dat')          
      
c
c --- datenfile oeffnen 
c
      pi = 4.d0*datan(1.d0)
      do j = 1,8
       do i=1,ndata
          write(14+j,100) xdata(i),ydata1(j,i)
       end do
      close(14+j)
      end do
c
      
      return
      end
      subroutine plot2(ndata,xdata,ydata1)
c
c --- datafile schreiben 
c  
         
          
      integer ndata,iunit
      real*8 xdata(ndata),ydata1(8,ndata)
      real*8 pi
100   format(2f16.8)
      open(25,file='c:\fortran\qbeam\output\int0l.dat') 
      open(26,file='c:\fortran\qbeam\output\int1l.dat')
      open(27,file='c:\fortran\qbeam\output\int2l.dat')
      open(28,file='c:\fortran\qbeam\output\int3l.dat')
      open(29,file='c:\fortran\qbeam\output\int4l.dat')
      open(30,file='c:\fortran\qbeam\output\int5l.dat')          
      open(31,file='c:\fortran\qbeam\output\int6l.dat')
      open(32,file='c:\fortran\qbeam\output\int7l.dat')          

c
c --- datenfile oeffnen 
c
      pi = 4.d0*datan(1.d0)
      do j = 1,8
       do i=1,ndata
          write(24+j,100) xdata(i),ydata1(j,i)
       end do
      close(14+j)
      end do
c
      
      return
      end

c
c
c
c
c --- normal theta lamda scan 
c
      subroutine plot3d(bi,iunit,ndata,xdata,alamda,ydata)
c
c --- datafile schreiben 
c  
         
          
      integer ndata,iunit
      real*8 xdata(ndata),ydata(8,ndata),alamda,bi
      real*8 pi
100   format(10f16.8)  
c
      pi = 4.d0*datan(1.d0)
c      
      do i=1,ndata
          if(i .eq. 3) ydata(j,i) = ydata(j,i)*bi/0.5
          if(i .eq. 4) ydata(j,i) = ydata(j,i)*bi/2.
          write(iunit,100) alamda,xdata(i),(ydata(j,i),j=1,8)
      end do
      write(iunit,*)
c
      return
      end
