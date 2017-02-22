c======================================================================c
c------------- LAGRANGIAN PARTICLE TRACKING RESTART CODE --------------c
c======================================================================c
      subroutine lpt_restart (IWRITE) ! 1 = Write new file, 0 = Open

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      common /particles_info/ ppos(ldim,lpart), ppos0(ldim,lpart)
     $                      , pvel(ldim,lpart), pvel0(ldim,lpart)

      integer IWRITE

300   format(F9.5,F9.5,F9.5,F9.5,F9.5,F9.5)

c     === OPEN RESTART FILE ===
      if(IWRITE .EQ. 0) then

      	if(nid.eq.0) write(*,*) "Opening restart file."
      	open(unit = 80,file = "partrestart.dat")

      	do jp=1,lpart
      		read(80,300) ppos0(1,jp),ppos0(2,jp),ppos0(3,jp),
     $                 pvel0(1,jp),pvel0(2,jp),pvel0(3,jp)
      	enddo

      	close(80)
      endif

c     === WRITE TO RESTART FILE ===
      if(IWRITE .EQ. 1) then

      	if(nid.eq.0) write(*,*) "Writing restart file."
      	if(nid.eq.0) open(unit = 80,file = "partrestart.dat")
      	if(nid.eq.0) close(unit = 80, STATUS='DELETE') 
      	if(nid.eq.0) open(unit = 80,file = "partrestart.dat")

      	do jp=1,lpart
      		if(nid.eq.0) write(80,300)(ppos0(k,jp),k=1,3)
     $                             ,(pvel0(k,jp),k=1,3)
      	enddo

      	if(nid.eq.0) close(80)
      endif

      return
      end
c-----------------------------------------------------------------------c
