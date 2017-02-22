c========================================================================c
c----------- LAGRANGIAN PARTICLE TRACKING WALL CONCENTRATION ------------c
c========================================================================c
      subroutine lpt_wallconc

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      common /statswitches/ I3DPLOT, ISTATGET, ISTATFLUSH, IFORCESTAT
     $                    , IWALLCONC

      common /particles_stat/ pvxc(ldim*40),pvyc(ldim*40),pvzc(ldim*40)
     $                    , pvx2c(ldim*40),pvy2c(ldim*40),pvz2c(ldim*40)
     $                    , pvxm(ldim*40),pvym(ldim*40),pvzm(ldim*40)
     $                    , pnc(ldim*40)
     $                    , pvxs(ldim*40),pvys(ldim*40),pvzs(ldim*40)
     $                    , pvx2s(ldim*40),pvy2s(ldim*40),pvz2s(ldim*40)
     $                    , pncs(ldim*40)

      common /wallconc/ wallcount, wallavg

4     format(1p16e17.9)

      wallcount = (pncs(1) + pncs(2)
     $          +  pncs(3) + pncs(4)
     $          +  pncs(5) + pncs(6)
     $          +  pncs(7) + pncs(8))/4.0D0

c     === OPEN FILE FOR WRITING STATISTICS ===
      if(nid.eq.0) open(unit = 178,file = 'partstat/partconc.dat'
     $                 ,position="append")

      wallavg = wallcount/real(ISTATFLUSH)

      if(nid.eq.0) write(178,4) time, wallavg

c     === ZERO STATS ARRAYS ===
      wallavg = 0.0D0

      if(nid.eq.0) close(178)

      return
      end
c-----------------------------------------------------------------------c
