c======================================================================c
c------------ LAGRANGIAN PARTICLE TRACKING COLLISION CODE -------------c
c======================================================================c
      subroutine lpt_collide

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      common /flow_prop/ brey, srey, bwall(ldim*2), bcwall(ldim*2)
     $                 , fvel(ldim,lpart)
     $              , fdvu(ldim,lpart),fdvv(ldim,lpart),fdvw(ldim,lpart)

      common /particles_prop/ diam_p, dt_p, dens_pf, rad_p, vol_p
     $                      , stkb_p, stkt_p, vfrac_p

      common /particles_info/ ppos(ldim,lpart), ppos0(ldim,lpart)
     $                      , pvel(ldim,lpart), pvel0(ldim,lpart)

      integer dd ! Actual boundary direction, X=1, Y=2, Z=3
      real rebound ! Rebound distance

c     === LOOP OVER ALL PARTICLES TO CHECK FOR COLLISIONS ===
      do jp=1,lpart
c     	=== CHECK BOUNDARIES FOR COLLISION ===
      	do d=1,(2*ldim)

c         === GET DIRECTION REFS ===
      		if(d .EQ. 1 .OR. d .EQ. 2) dd = 1
      		if(d .EQ. 3 .OR. d .EQ. 4) dd = 2
      		if(d .EQ. 5 .OR. d .EQ. 6) dd = 3

c     		=== PERIODIC BOUNDARY ===
      		if(bcwall(d) .EQ. 1) then

c     			=== CHECK LOWER WALL ===
      			if(MOD(d,2.) .NE. 0) then 
      				if((ppos(dd,jp)) .LT. bwall(d)) then
      					ppos(dd,jp) = bwall(d+1) - (1.0D-6)
      				endif
c     			=== CHECK UPPER WALL ===
      			else if (MOD(d,2.) .EQ. 0) then
      				if((ppos(dd,jp)) .GT. bwall(d)) then
      					ppos(dd,jp) = bwall(d-1) + (1.0D-6)
      				endif
      			endif
      		endif

c     		=== WALL BOUNDARY ===
      		if(bcwall(d) .EQ. 0) then

c     			=== CHECK LOWER WALL ===
      			if(MOD(d,2.) .NE. 0) then 
      				if((ppos(dd,jp)-rad_p) .LE. bwall(d)) then
      					rebound = ppos(dd,jp)-rad_p-bwall(d)
      					ppos(dd,jp) = (bwall(d) + rad_p) - rebound
      					pvel(dd,jp) = -1.0D0*pvel(dd,jp)
      				endif
c     			=== CHECK UPPER WALL ===
      			else if (MOD(d,2.) .EQ. 0) then
      				if((ppos(dd,jp)+rad_p) .GE. bwall(d)) then
      					rebound = ppos(dd,jp)+rad_p-bwall(d)
      					ppos(dd,jp) = (bwall(d) - rad_p) - rebound
      					pvel(dd,jp) = -1.0D0*pvel(dd,jp)
      				endif
      			endif
      		endif

      	enddo
c     	=== END OF BOUNDARY CHECKING ===

      enddo
c     === END OF LOOP OVER ALL PARTICLES TO CHECK FOR COLLISIONS ===

      return
      end
c-----------------------------------------------------------------------c
