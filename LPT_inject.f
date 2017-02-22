c======================================================================c
c------------ LAGRANGIAN PARTICLE TRACKING INJECTION CODE -------------c
c======================================================================c
      subroutine lpt_inject

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

 100  FORMAT (A9,I6.1,A,F7.3,F7.3,F7.3)
  99  FORMAT (I3.1,A)

      if(nid.EQ.0) write(*,*) "==================================="
      if(nid.EQ.0) write(*,*) "FIRST TIMESTEP. INJECTING PARTICLES"
      if(nid.EQ.0) write(*,*) "==================================="

c     === RANDOM INJECTION - LOOP OVER ALL PARTICLES ===
      do jp = 1,lpart
 10   	CONTINUE
c     	=== GET RANDOM X POSITION ===
      	ppos0(1,jp) = ((drand(0)*2.0D0)-1.0D0)*(bwall(2)-0.5D0)
      	if((ppos0(1,jp) + rad_p) .GE. bwall(2)) goto 10
      	if((ppos0(1,jp) - rad_p) .LE. bwall(1)) goto 10

 20   	CONTINUE
c     	=== GET RANDOM Y POSITION ===
      	ppos0(2,jp) = ((drand(0)*2.0D0)-1.0D0)*(bwall(4)-0.1D0)
      	if((ppos0(2,jp) + rad_p) .GE. bwall(4)) goto 20
      	if((ppos0(2,jp) - rad_p) .LE. bwall(3)) goto 20

 30   	CONTINUE
c     	=== GET RANDOM Z POSITION ===
      	ppos0(3,jp) = ((drand(0)*2.0D0)-1.0D0)*(bwall(6)-0.1D0)
      	if((ppos0(3,jp) + rad_p) .GE. bwall(6)) goto 30
      	if((ppos0(3,jp) - rad_p) .LE. bwall(5)) goto 30

c     	=== PARTICLE HAS POSITION, MAKE SURE IT'S NOT TOUCHING ANY OTHER ===
      	if(IFOURWAY.EQ.1) then
      		do jpp = 1,(jp-1)
      			if(sqrt(
     $             (ppos0(1,jpp)-ppos0(1,jp))**2
     $            +(ppos0(2,jpp)-ppos0(2,jp))**2
     $            +(ppos0(3,jpp)-ppos0(3,jp))**2
     $             ) .LE. diam_p) then
      				goto 10
      			endif
      		enddo
      	endif

      	if(nid.EQ.0 .AND. mod(jp,MAX((lpart/100),1)) .EQ. 0) then
      		write(*,99) 100*jp/lpart,'%'
      	endif        
      enddo
c     === FINISHED INJECTION LOOP OVER ALL PARTICLES ===

c     === GIVE PARTICLES VELOCITIES EQUAL TO CURRENT FLUID VELOCITY ===
      call interp_v(pvel0,ppos0,lpart)
      call interp_v(fvel,ppos0,lpart)

      if(nid.EQ.0) write(*,*) "FINISHED INJECTING"

      return
      end
c-----------------------------------------------------------------------c
