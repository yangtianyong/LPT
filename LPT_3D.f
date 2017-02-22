c======================================================================c
c----------- LAGRANGIAN PARTICLE TRACKING 3D PLOTTING CODE ------------c
c======================================================================c
      subroutine lpt_3D

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      common /particles_info/ ppos(ldim,lpart), ppos0(ldim,lpart)
     $                      , pvel(ldim,lpart), pvel0(ldim,lpart)

      common /statswitches/ I3DPLOT, ISTATGET, ISTATFLUSH, IFORCESTAT

      character*128 filename3D

208   format('partstat/3D/part',i5.5,'.3D')
2     format(F9.5,F9.5,F9.5)

c     === OPEN FILE PER TIMESTEP ===
      write(filename3D,208) istep
      open(unit = 73,file = filename3D)

c     === FILE FOR GNUPLOT TRACKING ===
      open(unit = 72,file = 'partstat/3D/particles.3D')

      write(73,*) 'X Y Z'
      do jp = 1, lpart
      	if(nid.eq.0) write(73,2)(ppos0(k,jp),k=1,3)
      enddo

      do jp = 1, lpart
      	if(nid.eq.0) write(72,2)(ppos0(k,jp),k=1,3)
      enddo

      close(73)
      close(72)

      return
      end
c-----------------------------------------------------------------------c
