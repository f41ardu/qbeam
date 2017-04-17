	function gikn(T,ThetaBragg)
c 
c --- berechnet fuer einen Braggwinkel, die temperaturabhaengige
c     Gitterkonstante
c 8v7 miukbnmutttbuuuuuuuuuuuummmmmmuzcnvh ccccnnbzzd<essxd 
	real a0,T,ThetaBragg,gikn
        real deltaTheta
        real pi
	pi = 4.*atan(1.)
	a0 = 5.431002
        gikn = a0*(1. - 5.e-6*(deltaTheta(T)/
     .           tan(pi/2. - ThetaBragg)))
	end
c
	function deltaTheta(T)
c
	real a0,a1,a2,a3,a4,a5
	real deltaTheta,T
        real sum
c
c --- parameter nach bartonitz
c 
	a0 =  15.76
        a1 = -0.4246e-1
        a2 =  0.1353e-2
        a3 = -0.1006e-4
        a4 =  0.2433e-7
        a5 = -0.2120e-10
c
	sum = a0+a1*T+a2*T**2+a3*T**3+a4*T**4+a5*T**5
c
        deltaTheta = sum
	end
