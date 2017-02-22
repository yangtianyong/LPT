c===============================================================================c
c----------- LAGRANGIAN PARTICLE TRACKING STATISTICAL FLUSHING CODE ------------c
c===============================================================================c
      subroutine lpt_statflush

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      common /particles_stat/ pvxc(ldim*40),pvyc(ldim*40),pvzc(ldim*40)
     $                    , pvx2c(ldim*40),pvy2c(ldim*40),pvz2c(ldim*40)
     $                    , pvxm(ldim*40),pvym(ldim*40),pvzm(ldim*40)
     $                    , pnc(ldim*40)
     $                    , pvxs(ldim*40),pvys(ldim*40),pvzs(ldim*40)
     $                    , pvx2s(ldim*40),pvy2s(ldim*40),pvz2s(ldim*40)
     $                    , pncs(ldim*40)

      common /flow_prop/ brey, srey, bwall(ldim*2), bcwall(ldim*2)
     $                 , fvel(ldim,lpart)
     $              , fdvu(ldim,lpart),fdvv(ldim,lpart),fdvw(ldim,lpart)

      common /fourway/ avgcol(ldim*10), avgcolc(ldim*10) ! AVERAGE NUMBER OF COLLISIONS IN SLICE Y
     $               , pxcol(ldim*10), pxcolc(ldim*10) ! AVERAGE NET X MOMENTUM CHANGE IN SLICE Y
     $               , pycol(ldim*10), pycolc(ldim*10) ! AVERAGE NET Y MOMENTUM CHANGE IN SLICE Y
     $               , pzcol(ldim*10), pzcolc(ldim*10) ! AVERAGE NET Z MOMENTUM CHANGE IN SLICE Y
     $               , xmxcol, ymxcol, zmxcol
     $               , dxcol, dycol, dzcol

      common /statswitches/ I3DPLOT, ISTATGET, ISTATFLUSH, IFORCESTAT
     $                    , IWALLCONC

      common /forcestat/ pf1cx(ldim*40)
     $                 , pf1cy(ldim*40)
     $                 , pf1cz(ldim*40)
     $                 , pf2cx(ldim*40)
     $                 , pf2cy(ldim*40)
     $                 , pf2cz(ldim*40)
     $                 , pf3cx(ldim*40)
     $                 , pf3cy(ldim*40)
     $                 , pf3cz(ldim*40)
     $                 , pf4cx(ldim*40)
     $                 , pf4cy(ldim*40)
     $                 , pf4cz(ldim*40)
     $                 , pf5cx(ldim*40)
     $                 , pf5cy(ldim*40)
     $                 , pf5cz(ldim*40)
     $                 , pf1sx(ldim*40)
     $                 , pf1sy(ldim*40)
     $                 , pf1sz(ldim*40)
     $                 , pf2sx(ldim*40)
     $                 , pf2sy(ldim*40)
     $                 , pf2sz(ldim*40)
     $                 , pf3sx(ldim*40)
     $                 , pf3sy(ldim*40)
     $                 , pf3sz(ldim*40)
     $                 , pf4sx(ldim*40)
     $                 , pf4sy(ldim*40)
     $                 , pf4sz(ldim*40)
     $                 , pf5sx(ldim*40)
     $                 , pf5sy(ldim*40)
     $                 , pf5sz(ldim*40)

      real bl ! Bin length for stat gathering
      real cl ! Bin length for collision stat gathering

      character*128 filenamestat
      character*128 filenamecstat
      character*128 filenamefstat

      bl = (2.0*bwall(4))/real(ldim*40)
      cl = (2.0*bwall(4))/real(ldim*10)

c     === OPEN FILE FOR WRITING STATISTICS ===

      if(nid.eq.0) open(unit = 77,file = 'partstat/partstat.dat')
      if(nid.eq.0) write(77,'(A,1pe14.7)') '#time = ', time
      if(nid.eq.0) write(77,'(A)')
     $'# y  n  u  v  w  u2  v2  w2'

      if(nid.eq.0) open(unit = 78,file = 'partstat/partcolstat.dat')
      if(nid.eq.0) write(78,'(A,1pe14.7)') '#time = ', time
      if(nid.eq.0) write(78,'(A)')
     $'# y  ncol  pw  pv  pw'

      if(nid.eq.0) open(unit = 79,file = 'partstat/partfstat.dat')
      if(nid.eq.0) write(79,'(A,1pe14.7)') '#time = ', time
      if(nid.eq.0) write(79,'(A)')
     $'# y f1x f1y f1z f2x f2y f2z f3x f3y f3z f4x f4y f4z f5x f5y f5z'

      do s=1,ldim*40
      	pvxm(s) = pvxs(s)/real(istatflush)
      	pvym(s) = pvys(s)/real(istatflush)
      	pvzm(s) = pvzs(s)/real(istatflush)
      enddo

4     format(1p16e17.9)

      do s=1,ldim*40
      	if(nid.eq.0) write(77,4) (((bl*(s-1))-1.0D0)+0.5D0*bl)
     $  , pncs(s)/real(istatflush)
     $  , pvxs(s)/real(istatflush)
     $  , pvys(s)/real(istatflush)
     $  , pvzs(s)/real(istatflush)
     $  , pvx2s(s)/real(istatflush)-(pvxs(s)/real(istatflush))**2
     $  , pvy2s(s)/real(istatflush)-(pvys(s)/real(istatflush))**2
     $  , pvz2s(s)/real(istatflush)-(pvzs(s)/real(istatflush))**2

      	if(nid.eq.0) write(79,4) (((bl*(s-1))-1.0D0)+0.5D0*bl)
     $  , pf1sx(s)/real(istatflush)
     $  , pf1sy(s)/real(istatflush)
     $  , pf1sz(s)/real(istatflush)
     $  , pf2sx(s)/real(istatflush)
     $  , pf2sy(s)/real(istatflush)
     $  , pf2sz(s)/real(istatflush)
     $  , pf3sx(s)/real(istatflush)
     $  , pf3sy(s)/real(istatflush)
     $  , pf3sz(s)/real(istatflush)
     $  , pf4sx(s)/real(istatflush)
     $  , pf4sy(s)/real(istatflush)
     $  , pf4sz(s)/real(istatflush)
     $  , pf5sx(s)/real(istatflush)
     $  , pf5sy(s)/real(istatflush)
     $  , pf5sz(s)/real(istatflush)
      enddo

      do s=1,ldim*10
      	if(nid.eq.0) write(78,4) (((cl*(s-1))-1.0D0)+0.5D0*cl)
     $  , avgcol(s)/real(istatflush)
     $  , pxcol(s)/real(istatflush)
     $  , pycol(s)/real(istatflush)
     $  , pzcol(s)/real(istatflush)
      enddo

c     === ZERO STATS ARRAYS ===
      do d = 1,ldim*40
      	pvxc(d) = 0.0D0
      	pvyc(d) = 0.0D0
      	pvzc(d) = 0.0D0
      	pvxs(d) = 0.0D0
      	pvys(d) = 0.0D0
      	pvzs(d) = 0.0D0
      	pvx2c(d) = 0.0D0
      	pvy2c(d) = 0.0D0
      	pvz2c(d) = 0.0D0
      	pvx2s(d) = 0.0D0
      	pvy2s(d) = 0.0D0
      	pvz2s(d) = 0.0D0
      	pnc(d) = 0.0D0
      	pncs(d) = 0.0D0
      	pf1sx(d) = 0.0D0
      	pf1sy(d) = 0.0D0
      	pf1sz(d) = 0.0D0
      	pf2sx(d) = 0.0D0
      	pf2sy(d) = 0.0D0
      	pf2sz(d) = 0.0D0
      	pf3sx(d) = 0.0D0
      	pf3sy(d) = 0.0D0
      	pf3sz(d) = 0.0D0
      	pf4sx(d) = 0.0D0
      	pf4sy(d) = 0.0D0
      	pf4sz(d) = 0.0D0
      	pf5sx(d) = 0.0D0
      	pf5sy(d) = 0.0D0
      	pf5sz(d) = 0.0D0
      enddo

      do d = 1,ldim*10
      	avgcol(d) = 0.0D0
      	pxcol(d) = 0.0D0
      	pycol(d) = 0.0D0
      	pzcol(d) = 0.0D0
      enddo

      if(nid.eq.0) close(77)
      if(nid.eq.0) close(78)
      if(nid.eq.0) close(79)

      return
      end
c-----------------------------------------------------------------------c
