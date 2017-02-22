c======================================================================c
c--------------- LAGRANGIAN PARTICLE TRACKING MAIN CODE ---------------c
c======================================================================c
c     NOTE: VARIABLES STARTING WITH I,J,K,L,M,N ARE INTEGERS           c
c                                                                      c
c     SWITCHES:                                                        c
c     IPRESTART         -  Restart from file? (1-Yes,0-No)             c
c     IDRAG             -  Drag force (1-On,0-Off)                     c
c     ILIFT             -  Lift force (1-On,0-Off)                     c
c     IGRAV             -  Gravity/Buoyancy force (1-On,0-Off)         c
c     IPRES             -  Pressure Gradient force (1-On,0-Off)        c
c     IMASS             -  Virtual Mass force (1-On,0-Off)             c
c     ITWOWAY           -  Two-way coupling (1-On,0-Off)               c
c     IFOURWAY          -  Four-way coupling (1-On,0-Off)              c
c     IAGGLOM           -  DLVO Agglomeration (1-On,0-Off)             c
c                                                                      c
c     STATSWITCHES:                                                    c
c     I3DPLOT           -  Particle 3D plotting (X-Frequency,0-No)     c
c     ISTATGET          -  Particle stat recorder (X-Frequency,0-No)   c
c     ISTATFLUSH        -  Particle stat write (X-Frequency,0-No)      c
c                                                                      c
c     FLOW PROPERTIES:                                                 c
c     brey              -  Flow bulk reynolds number                   c
c     srey              -  Flow shear reynolds number                  c
c     bwall(ldim*2)     -  Domain boundaries                           c
c     bcwall(ldim*2)    -  Boundary condition                          c
c     fvel(ldim,lpart)  -  Fluid velocity at particle positions        c
c     fdvu(ldim,lpart)  -  Fluid u div field at particle positions     c
c     fdvv(ldim,lpart)  -  Fluid v div field at particle positions     c
c     fdvw(ldim,lpart)  -  Fluid w div field at particle positions     c
c                                                                      c
c     PARTICLE TRACKER VARIABLES:                                      c
c     diam_p            -  Particle ND diameter                        c
c     dt_p              -  Particle timestep                           c
c     dens_pf           -  Particle ND density                         c
c     rad_p             -  Particle ND radius                          c
c     stkb_p            -  Particle bulk stokes number                 c
c     stkt_p            -  Particle shear stokes number                c
c     vfrac_p           -  Particle volume fraction                    c
c     ppos(ldim,lpart)  -  Particle current position                   c
c     ppos0(ldim,lpart) -  Particle position last step                 c
c     pvel(ldim,lpart)  -  Particle current velocity                   c
c     pvel0(ldim,lpart) -  Particle velocity last step                 c
c                                                                      c
c     TWO-WAY COUPLING VARIABLES:                                      c
c     twowayforcex,y,z  -  Force on fluid element (i,j,k,v)            c
c     xmc, xxc          -  Minimum and maximum x value for cell        c
c     ymc, yxc          -  Minimum and maximum y value for cell        c
c     zmc, zxc          -  Minimum and maximum z value for cell        c
c     vol               -  Volume of cell (i,j,k,v)                    c
c======================================================================c
      subroutine lpt_main

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'  ! for nelx,nely,nelz

      common /switches/ IPRESTART, IDRAG, ILIFT, IGRAV, IPRES, IMASS
     $                , ITWOWAY, IFOURWAY, IAGGLOM

      common /statswitches/ I3DPLOT, ISTATGET, ISTATFLUSH, IFORCESTAT
     $                    , IWALLCONC

      if(nid.EQ.0) write(*,*) "ENTERING LPT_MAIN. ISTEP = ", istep

c     === CALL INITIALISATION IF FIRST TIMESTEP ===     
      if(istep.EQ.0) then
      	if(nid.EQ.0) write(*,*) "ENTERING LPT_INIT."
      	call lpt_init
      endif

c     === CALL INJECTION OR RESTART===
      if(istep.EQ.0 .AND. IPRESTART .EQ. 0) then
      	if(nid.EQ.0) write(*,*) "ENTERING LPT_INJECT."
      	call lpt_inject
      else if(istep.EQ.0 .AND. IPRESTART .EQ. 1) then
      	if(nid.EQ.0) write(*,*) "ENTERING LPT_RESTART."
      	call lpt_restart(0)
      endif

c     === CALL TWO-WAY SETUP ===
      if(istep.EQ.0 .AND. ITWOWAY .EQ. 1) then
      	if(nid.EQ.0) write(*,*) "ENTERING LPT_TWOWAY(0)."
      	call lpt_twoway(0)
      endif

c     === CALL TWO-WAY ===
      if(istep.GT.0 .AND. ITWOWAY .EQ. 1) then
      	if(nid.EQ.0) write(*,*) "ENTERING LPT_TWOWAY(1)."
      	call lpt_twoway(1)
      endif

c     === CALL PARTICLE FORCING FOR ADVECTION ===
      if(istep.GT.0) then
      	if(nid.EQ.0) write(*,*) "ENTERING LPT_ADVECT."
      	call lpt_advect
      endif

c     === CALL RESTART WRITING ===
      if(istep.GT.0 .AND. mod(ISTEP,IOSTEP) .EQ. 0) then
      	call lpt_restart(1)
      endif

c     === CALL 3D PLOTTING FUNCTION ===
      if(istep.GT.0 .AND. I3DPLOT .GT. 0) then
      	if(MOD(ISTEP,I3DPLOT) .EQ. 0) then
      		if(nid.eq.0) write(*,*) 'Calling LPT_3D.'
      		call lpt_3D
      	endif
      endif

c     === CALL STATS GATHERING FUNCTION ===
      if(istep.GT.0) then
      	if(MOD(ISTEP,ISTATGET) .EQ. 0) then
      		if(nid.eq.0) write(*,*) 'Calling LPT_STATGET.'
      		call lpt_statget
      	endif
      endif

c     === CALL WALLCONC FLUSHING FUNCTION ===
      if(istep.GT.0) then
      	if(MOD(ISTEP,ISTATFLUSH) .EQ. 0 
     $           .AND. IWALLCONC .EQ. 1) then
      		if(nid.eq.0) write(*,*) 'Calling LPT_WALLCONC.'
      		call lpt_wallconc
      	endif
      endif

c     === CALL STATS FLUSHING FUNCTION ===
      if(istep.GT.0) then
      	if(MOD(ISTEP,ISTATFLUSH) .EQ. 0) then
      		if(nid.eq.0) write(*,*) 'Calling LPT_STATFLUSH.'
      		call lpt_statflush
      	endif
      endif

      return
      end
c-----------------------------------------------------------------------c
