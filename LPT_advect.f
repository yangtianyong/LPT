c======================================================================c
c------------ LAGRANGIAN PARTICLE TRACKING ADVECTION CODE -------------c
c======================================================================c
      subroutine lpt_advect

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      common /switches/ IPRESTART, IDRAG, ILIFT, IGRAV, IPRES, IMASS
     $                , ITWOWAY, IFOURWAY, IAGGLOM

      common /flow_prop/ brey, srey, bwall(ldim*2), bcwall(ldim*2)
     $                 , fvel(ldim,lpart)
     $              , fdvu(ldim,lpart),fdvv(ldim,lpart),fdvw(ldim,lpart)

      common /particles_prop/ diam_p, dt_p, dens_pf, rad_p, vol_p
     $                      , stkb_p, stkt_p, vfrac_p

      common /particles_info/ ppos(ldim,lpart), ppos0(ldim,lpart)
     $                      , pvel(ldim,lpart), pvel0(ldim,lpart)

      common /grads/  dur(lmaxv),dus(lmaxv),dut(lmaxv)
     $             ,  dvr(lmaxv),dvs(lmaxv),dvt(lmaxv)
     $             ,  dwr(lmaxv),dws(lmaxv),dwt(lmaxv)

      common /rk4/ rkpx1(ldim,lpart), rkpx2(ldim,lpart)
     $           , rkpx3(ldim,lpart), rkpx4(ldim,lpart)
     $           , rkpv1(ldim,lpart), rkpv2(ldim,lpart)
     $           , rkpv3(ldim,lpart), rkpv4(ldim,lpart)
     $           , rkpa1(ldim,lpart), rkpa2(ldim,lpart)
     $           , rkpa3(ldim,lpart), rkpa4(ldim,lpart)

 92   format (A,F9.5,F9.5,F9.5) 

      if(nid.EQ.0) write(*,*) 'Advecting particles.'

c     === STEP 1 ==========================================

c     === WRITE CURRENT PARTICLE POSITION AND VELOCITY FOR RK4 ===
      do jp = 1, lpart
      	do d=1,3
      		rkpx1(d,jp) = ppos0(d,jp)
      		rkpv1(d,jp) = pvel0(d,jp)
      	enddo
      enddo

c     === INTERPOLATE CURRENT VELOCITY FIELDS ===

      ntot  = nx1*ny1*nz1*nelv        ! Total number of cells

      call interp_v(fvel,ppos0,lpart) ! Fluid velocity at particle positions
      call gradm1(dur,dus,dut,vx)     ! Fluid du/dx, du/dy, du/dz at particle positions
      call gradm1(dvr,dvs,dvt,vy)     ! Fluid dv/dx, dv/dy, dv/dz at particle positions
      call gradm1(dwr,dws,dwt,vz)     ! Fluid dw/dx, dw/dy, dw/dz at particle positions
      call interp_v3v(fdvu,ppos0,dur,dus,dut,lpart) ! Fluid U-div field at particle positions
      call interp_v3v(fdvv,ppos0,dvr,dvs,dvt,lpart) ! Fluid V-div field at particle positions
      call interp_v3v(fdvw,ppos0,dwr,dws,dwt,lpart) ! Fluid W-div field at particle positions

      call lpt_force(rkpx1,rkpv1,rkpa1,1)

c     =====================================================

c     === STEP 2 ==========================================

c     === CALCULATE NEW POSITIONS AND VELOCITIES FOR K2 ===

      do jp = 1, lpart
      	do d=1,3
      		rkpx2(d,jp) = ppos0(d,jp) + 0.5D0*rkpv1(d,jp)*dt_p
      		rkpv2(d,jp) = pvel0(d,jp) + 0.5D0*rkpa1(d,jp)*dt_p
      	enddo
      enddo

      call lpt_force(rkpx2,rkpv2,rkpa2,0)

c     =====================================================

c     === STEP 3 ==========================================

c     === CALCULATE NEW POSITIONS AND VELOCITIES FOR K3 ===

      do jp = 1, lpart
      	do d=1,3
      		rkpx3(d,jp) = ppos0(d,jp) + 0.5D0*rkpv2(d,jp)*dt_p
      		rkpv3(d,jp) = pvel0(d,jp) + 0.5D0*rkpa2(d,jp)*dt_p
      	enddo
      enddo

      call lpt_force(rkpx3,rkpv3,rkpa3,0)

c     =====================================================

c     === STEP 4 ==========================================

c     === CALCULATE NEW POSITIONS AND VELOCITIES FOR K4 ===

      do jp = 1, lpart
      	do d=1,3
      		rkpx4(d,jp) = ppos0(d,jp) + rkpv3(d,jp)*dt_p
      		rkpv4(d,jp) = pvel0(d,jp) + rkpa3(d,jp)*dt_p
      	enddo
      enddo

      call lpt_force(rkpx4,rkpv4,rkpa4,0)

c     =====================================================

c     === UPDATE PARTICLE POSITIONS AND VELOCITIES ===

      do jp = 1, lpart
      	do d=1,3
      		ppos(d,jp) = ppos0(d,jp) + (dt_p/6.0D0)*(rkpv1(d,jp)
     $                                      +2.0D0*rkpv2(d,jp)
     $                                      +2.0D0*rkpv3(d,jp)
     $                                            +rkpv4(d,jp))

      		pvel(d,jp) = pvel0(d,jp) + (dt_p/6.0D0)*(rkpa1(d,jp)
     $                                      +2.0D0*rkpa2(d,jp)
     $                                      +2.0D0*rkpa3(d,jp)
     $                                            +rkpa4(d,jp))
      	enddo
      enddo

c     === CALL COLLISION SUBROUTINE TO CHECK FOR WALL COLLISIONS ===

      call LPT_collide

c     === CALL INTERPARTICLE COLLISION ===

      if(IFOURWAY .EQ. 1) then
      	if(nid.EQ.0) write(*,*) "ENTERING LPT_IP_COLLIDE."
      	call LPT_ip_collide
      endif

c     === CALL COLLISION AGAIN SUBROUTINE TO CHECK FOR WALL COLLISIONS ===

      call LPT_collide

c     === RESET PARTICLE START POSITIONS/VELOCITIES ===
      do jp = 1, lpart
      	do d=1,3
      		ppos0(d,jp) = ppos(d,jp)
      		pvel0(d,jp) = pvel(d,jp)
      	enddo
      enddo

      if(nid.EQ.0) write(*,*) 'Finished advecting particles.'       

      return
      end

c===========================================================================c
c------------ LAGRANGIAN PARTICLE TRACKING LINEAR FORCING CODE -------------c
c===========================================================================c
      subroutine lpt_force(rkx, rkv, rka, IFORCE)

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      common /switches/ IPRESTART, IDRAG, ILIFT, IGRAV, IPRES, IMASS
     $                , ITWOWAY, IFOURWAY, IAGGLOM

      common /flow_prop/ brey, srey, bwall(ldim*2), bcwall(ldim*2)
     $                 , fvel(ldim,lpart)
     $              , fdvu(ldim,lpart),fdvv(ldim,lpart),fdvw(ldim,lpart)

      common /particles_prop/ diam_p, dt_p, dens_pf, rad_p, vol_p
     $                      , stkb_p, stkt_p, vfrac_p

      common /particles_info/ ppos(ldim,lpart), ppos0(ldim,lpart)
     $                      , pvel(ldim,lpart), pvel0(ldim,lpart)

      common /grads/  dur(lmaxv),dus(lmaxv),dut(lmaxv)
     $             ,  dvr(lmaxv),dvs(lmaxv),dvt(lmaxv)
     $             ,  dwr(lmaxv),dws(lmaxv),dwt(lmaxv)

      common /forces/ pf1(ldim,lpart)
     $              , pf2(ldim,lpart)
     $              , pf3(ldim,lpart)
     $              , pf4(ldim,lpart)
     $              , pf5(ldim,lpart)

      integer IFORCE

      real rkx(ldim,lpart) ! Current position to calc deriv with
      real rkv(ldim,lpart) ! Current velocity to calc deriv with
      real rka(ldim,lpart) ! Current acceleration to calc deriv with

      real usl(ldim) ! Slip velocity for current particle
      real urel ! Magnitude of slip velocity for current particle
      real mvort ! Magnitude of vorticity for current particle
      real pbrey ! Particle bulk reynolds number

c     === DRAG VARIABLES ===
      real fre ! Schiller drag correlation
      real cd ! Drag coefficient

c     === LIFT VARIABLES ===
      real res
      real betal
      real cl ! Lift coefficient

c     === PRESSURE GRADIENT VARIABLES ===
      real cgp(ldim)

c     === ADDED MASS VARIABLES ===
      real ccp(ldim)
      real mad ! Added mass contribution

      real fx(ldim*2) ! Holds the x component of each force
      real fy(ldim*2) ! Holds the y component of each force
      real fz(ldim*2) ! Holds the z component of each force

 90   format (A,F7.3,F7.3,F7.3)
 91   format (A,F7.3)
 92   format (A,F9.5,F9.5,F9.5)

c     === LOOP OVER ALL PARTICLES ===
      do jp = 1,lpart

c     	=== GET SLIP VELOCITY ===
      	usl(1) = fvel(1,jp) - rkv(1,jp)
      	usl(2) = fvel(2,jp) - rkv(2,jp)
      	usl(3) = fvel(3,jp) - rkv(3,jp)

c     	=== IF IDENTICAL, MAKE NONZERO ===
      	if(abs(usl(1)) .LE. 1.0E-30) usl(1) = 1.0E-30
      	if(abs(usl(2)) .LE. 1.0E-30) usl(2) = 1.0E-30
      	if(abs(usl(3)) .LE. 1.0E-30) usl(3) = 1.0E-30

c       === GET MAGNITUDE OF SLIP VELOCITY ===
      	urel = sqrt(usl(1)**2
     $             +usl(2)**2
     $             +usl(3)**2)

c       === IF ZERO, MAKE NONZERO ===
      	if(urel .LE. 1.0E-30) urel = 1.0E-30

c       === GET MAGNITUDE OF VORTICITY ===
      	mvort = sqrt((fdvw(2,jp)-fdvv(3,jp))**2
     $              +(fdvu(3,jp)-fdvw(1,jp))**2
     $              +(fdvv(1,jp)-fdvu(2,jp))**2)

c       === IF ZERO, MAKE NONZERO ===
      	if(mvort .LE. 1.0E-30) mvort = 1.0E-30

c       === GET PARTICLE REYNOLDS NUMBER ===
      	pbrey = brey*diam_p*urel

c       =====================================
c       ------------ DRAG FORCE -------------
c       =====================================
      	if(IDRAG.EQ.1) then
c     		=== Schiller-Naumann drag force ===

c     	=== Drag correlations below ===
      		if(pbrey .GE. 1000) then
      			fre = 0.44D0*pbrey/24.0D0
      		else if(pbrey .LE. 0.5) then
      			fre = 1.0D0
      		else
      			fre = 1.0 + 0.15D0*pbrey**0.687D0
      		end if

c     	=== Calculate drag coefficient ===
      		cd = fre*24.0D0/pbrey

c         === CALCULATE DRAG FORCE ===
      		fx(1) = (urel*usl(1)*3*cd)/(4*diam_p*dens_pf)
      		fy(1) = (urel*usl(2)*3*cd)/(4*diam_p*dens_pf)
      		fz(1) = (urel*usl(3)*3*cd)/(4*diam_p*dens_pf)

c     		=== RECORD DRAG FORCE ===
      		if(IFORCE.EQ.1) then
      			pf1(1,jp) = fx(1)
      			pf1(2,jp) = fy(1)
      			pf1(3,jp) = fz(1)
      		endif

c     			=== OUTPUT SOME DRAG FORCES ===
c     		if(nid.EQ.0 .AND. jp .LE. 10) then
c     			write(*,92) 'FD: ',fx(1),fy(1),fz(1)
c     		endif

      	else

      		fx(1) = 0.0D0
      		fy(1) = 0.0D0
      		fz(1) = 0.0D0

c     		=== RECORD DRAG FORCE ===
      		if(IFORCE.EQ.1) then
      			pf1(1,jp) = fx(1)
      			pf1(2,jp) = fy(1)
      			pf1(3,jp) = fz(1)
      		endif

      	endif

c     	=====================================
c     	------------ LIFT FORCE -------------
c     	=====================================
      	if(ILIFT.EQ.1) then

c     		=== CALCULATE LIFT FORCE ===
      		res = (diam_p**2)*brey*mvort
      		betal = 0.5D0*(res/pbrey)

      		if(pbrey.LE.40) then
      			frel = (1-0.3314*sqrt(betal))*
     $		         exp(-1.0D0*(pbrey/10.0D0))+
     $		         0.3314*sqrt(betal)
      		else
      			frel = 0.0524*sqrt(betal*pbrey)
      		endif

      		cl = (4.1126D0*frel)/(sqrt(res))

      		fx(2) = ((3.0D0*cl)/(4.0D0*dens_pf))
     $		*(usl(2)*(fdvv(1,jp)-fdvu(2,jp))
     $		 -usl(3)*(fdvu(3,jp)-fdvw(1,jp)))

      		fy(2) = ((3.0D0*cl)/(4.0D0*dens_pf))
     $		*(usl(3)*(fdvw(2,jp)-fdvv(3,jp))
     $		 -usl(1)*(fdvv(1,jp)-fdvu(2,jp)))

      		fz(2) = ((3.0D0*cl)/(4.0D0*dens_pf))
     $		*(usl(1)*(fdvu(3,jp)-fdvw(1,jp))
     $		 -usl(2)*(fdvw(2,jp)-fdvv(3,jp)))

c     		=== RECORD LIFT FORCE ===
      		if(IFORCE.EQ.1) then
      			pf2(1,jp) = fx(2)
      			pf2(2,jp) = fy(2)
      			pf2(3,jp) = fz(2)
      		endif

c     		=== OUTPUT SOME LIFT FORCES ===
c     		if(nid.EQ.0 .AND. jp .LE. 10) then
c     			write(*,92) 'FL: ',fx(2),fy(2),fz(2)
c     		endif

      	else

      		fx(2) = 0.0D0
      		fy(2) = 0.0D0
      		fz(2) = 0.0D0

c     		=== RECORD LIFT FORCE ===
      		if(IFORCE.EQ.1) then
      			pf2(1,jp) = fx(2)
      			pf2(2,jp) = fy(2)
      			pf2(3,jp) = fz(2)
      		endif

      	endif

c     	========================================
c     	------------ GRAVITY FORCE -------------
c     	========================================
      	if(IGRAV.EQ.1) then

c     		=== CALCULATE GRAVITY FORCE ===
      		fx(3) = 0.0D0
      		fy(3) = 1.0D0*9.81D0*(1-(1/dens_pf))
      		fz(3) = 0.0D0

c     		=== RECORD GRAVITY FORCE ===
      		if(IFORCE.EQ.1) then
      			pf3(1,jp) = fx(3)
      			pf3(2,jp) = fy(3)
      			pf3(3,jp) = fz(3)
      		endif

c     		=== OUTPUT SOME GRAV FORCES ===
c     		if(nid.EQ.0 .AND. jp .LE. 10) then
c     			write(*,92) 'FG: ',fx(3),fy(3),fz(3)
c     		endif

      	else

      		fx(3) = 0.0D0
      		fy(3) = 0.0D0
      		fz(3) = 0.0D0

c     		=== RECORD GRAVITY FORCE ===
      		if(IFORCE.EQ.1) then
      			pf3(1,jp) = fx(3)
      			pf3(2,jp) = fy(3)
      			pf3(3,jp) = fz(3)
      		endif

      	endif

c     	==================================================
c     	------------ PRESSURE GRADIENT FORCE -------------
c     	==================================================
      	if(IPRES.EQ.1) then

      		cc1 = fvel(1,jp)*fdvu(1,jp)
      		cc2 = fvel(2,jp)*fdvu(2,jp)
      		cc3 = fvel(3,jp)*fdvu(3,jp)

      		cgp(1) = cc1+cc2+cc3

      		cc1 = fvel(1,jp)*fdvv(1,jp)
      		cc2 = fvel(2,jp)*fdvv(2,jp)
      		cc3 = fvel(3,jp)*fdvv(3,jp)

      		cgp(2) = cc1+cc2+cc3

      		cc1 = fvel(1,jp)*fdvw(1,jp)
      		cc2 = fvel(2,jp)*fdvw(2,jp)
      		cc3 = fvel(3,jp)*fdvw(3,jp)

      		cgp(3) = cc1+cc2+cc3

c     		=== CALCULATE PRESSURE FORCE ===
      		fx(4) = (1.0D0/dens_pf)*cgp(1)
      		fy(4) = (1.0D0/dens_pf)*cgp(2)
      		fz(4) = (1.0D0/dens_pf)*cgp(3)

c     		=== RECORD PRESSURE FORCE ===
      		if(IFORCE.EQ.1) then
      			pf4(1,jp) = fx(4)
      			pf4(2,jp) = fy(4)
      			pf4(3,jp) = fz(4)
      		endif

c     		=== OUTPUT SOME PRES FORCES ===
c     		if(nid.EQ.0 .AND. jp .LE. 10) then
c     			write(*,92) 'FP: ',fx(4),fy(4),fz(4)
c     		endif

      	else

      		fx(4) = 0.0D0
      		fy(4) = 0.0D0
      		fz(4) = 0.0D0

c     		=== RECORD PRESSURE FORCE ===
      		if(IFORCE.EQ.1) then
      			pf4(1,jp) = fx(4)
      			pf4(2,jp) = fy(4)
      			pf4(3,jp) = fz(4)
      		endif

      	endif

c     	=============================================
c     	------------ VIRTUAL MASS FORCE -------------
c     	=============================================
      	if(IMASS.EQ.1) then

      		cc1 = rkv(1,jp)*fdvu(1,jp)
      		cc2 = rkv(2,jp)*fdvu(2,jp)
      		cc3 = rkv(3,jp)*fdvu(3,jp)

      		ccp(1) = cc1+cc2+cc3

      		cc1 = rkv(1,jp)*fdvv(1,jp)
      		cc2 = rkv(2,jp)*fdvv(2,jp)
      		cc3 = rkv(3,jp)*fdvv(3,jp)

      		ccp(2) = cc1+cc2+cc3

      		cc1 = rkv(1,jp)*fdvw(1,jp)
      		cc2 = rkv(2,jp)*fdvw(2,jp)
      		cc3 = rkv(3,jp)*fdvw(3,jp)

      		ccp(3) = cc1+cc2+cc3

c     		=== CALCULATE VIRTUAL MASS FORCE ===
      		fx(5) = (1.0D0/(2*dens_pf))*ccp(1)
      		fy(5) = (1.0D0/(2*dens_pf))*ccp(2)
      		fz(5) = (1.0D0/(2*dens_pf))*ccp(3)

c     		=== RECORD VIRTUAL MASS FORCE ===
      		if(IFORCE.EQ.1) then
      			pf5(1,jp) = fx(5)
      			pf5(2,jp) = fy(5)
      			pf5(3,jp) = fz(5)
      		endif

      	else

      		fx(5) = 0.0D0
      		fy(5) = 0.0D0
      		fz(5) = 0.0D0

c     		=== RECORD VIRTUAL MASS FORCE ===
      		if(IFORCE.EQ.1) then
      			pf5(1,jp) = fx(5)
      			pf5(2,jp) = fy(5)
      			pf5(3,jp) = fz(5)
      		endif

      	endif

c     	=== ADD ALL THE FORCES TOGETHER FOR NEW ACCELERATIONS ===
      	if(IMASS.EQ.1) then
      		mad = (1.0D0+(1.0D0/(2.0D0*dens_pf)))
      		rka(1,jp) = (fx(1)+fx(2)+fx(3)+fx(4)+fx(5))/mad
      		rka(2,jp) = (fy(1)+fy(2)+fy(3)+fy(4)+fy(5))/mad
      		rka(3,jp) = (fz(1)+fz(2)+fz(3)+fz(4)+fz(5))/mad
      	else if(IMASS.EQ.0) then
      		rka(1,jp) = fx(1)+fx(2)+fx(3)+fx(4)+fx(5)
      		rka(2,jp) = fy(1)+fy(2)+fy(3)+fy(4)+fy(5)
      		rka(3,jp) = fz(1)+fz(2)+fz(3)+fz(4)+fz(5)
      	endif

      enddo

      return
      end
c-----------------------------------------------------------------------c
