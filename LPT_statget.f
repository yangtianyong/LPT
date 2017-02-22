c================================================================================c
c----------- LAGRANGIAN PARTICLE TRACKING STATISTICAL GATHERING CODE ------------c
c================================================================================c
      subroutine lpt_statget

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      common /particles_info/ ppos(ldim,lpart), ppos0(ldim,lpart)
     $                      , pvel(ldim,lpart), pvel0(ldim,lpart)

      common /statswitches/ I3DPLOT, ISTATGET, ISTATFLUSH, IFORCESTAT
     $                    , IWALLCONC

      common /flow_prop/ brey, srey, bwall(ldim*2), bcwall(ldim*2)
     $                 , fvel(ldim,lpart)
     $              , fdvu(ldim,lpart),fdvv(ldim,lpart),fdvw(ldim,lpart)

      common /particles_stat/ pvxc(ldim*40),pvyc(ldim*40),pvzc(ldim*40)
     $                    , pvx2c(ldim*40),pvy2c(ldim*40),pvz2c(ldim*40)
     $                    , pvxm(ldim*40),pvym(ldim*40),pvzm(ldim*40)
     $                    , pnc(ldim*40)
     $                    , pvxs(ldim*40),pvys(ldim*40),pvzs(ldim*40)
     $                    , pvx2s(ldim*40),pvy2s(ldim*40),pvz2s(ldim*40)
     $                    , pncs(ldim*40)

      common /fourway/ avgcol(ldim*10), avgcolc(ldim*10) ! AVERAGE NUMBER OF COLLISIONS IN SLICE Y
     $               , pxcol(ldim*10), pxcolc(ldim*10) ! AVERAGE NET X MOMENTUM CHANGE IN SLICE Y
     $               , pycol(ldim*10), pycolc(ldim*10) ! AVERAGE NET Y MOMENTUM CHANGE IN SLICE Y
     $               , pzcol(ldim*10), pzcolc(ldim*10) ! AVERAGE NET Z MOMENTUM CHANGE IN SLICE Y
     $               , xmxcol, ymxcol, zmxcol
     $               , dxcol, dycol, dzcol

      common /forces/ pf1(ldim,lpart)
     $              , pf2(ldim,lpart)
     $              , pf3(ldim,lpart)
     $              , pf4(ldim,lpart)
     $              , pf5(ldim,lpart)

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

      bl = (2.0*bwall(4))/(ldim*40)

c     === STATS COUNTING ===
      do jp=1,lpart
      	do s=1,ldim*40
      		if((ppos0(2,jp) .GE. (bl*(s-1) - bwall(4))) .AND.
     $			 (ppos0(2,jp) .LT. (bl*(s)   - bwall(4)))) then

      			pnc(s) = pnc(s) + 1

      			pvxc(s) = pvxc(s) + pvel0(1,jp)
      			pvyc(s) = pvyc(s) + pvel0(2,jp)
      			pvzc(s) = pvzc(s) + pvel0(3,jp)

      			pvx2c(s) = pvx2c(s) + (pvel0(1,jp))**2
      			pvy2c(s) = pvy2c(s) + (pvel0(2,jp))**2
      			pvz2c(s) = pvz2c(s) + (pvel0(3,jp))**2

      			if(IFORCESTAT.EQ.1) then

      				pf1cx(s) = pf1cx(s) + pf1(1,jp)
      				pf1cy(s) = pf1cy(s) + pf1(2,jp)
      				pf1cz(s) = pf1cz(s) + pf1(3,jp)

      				pf2cx(s) = pf2cx(s) + pf2(1,jp)
      				pf2cy(s) = pf2cy(s) + pf2(2,jp)
      				pf2cz(s) = pf2cz(s) + pf2(3,jp)

      				pf3cx(s) = pf3cx(s) + pf3(1,jp)
      				pf3cy(s) = pf3cy(s) + pf3(2,jp)
      				pf3cz(s) = pf3cz(s) + pf3(3,jp)

      				pf4cx(s) = pf4cx(s) + pf4(1,jp)
      				pf4cy(s) = pf4cy(s) + pf4(2,jp)
      				pf4cz(s) = pf4cz(s) + pf4(3,jp)

      				pf5cx(s) = pf5cx(s) + pf5(1,jp)
      				pf5cy(s) = pf5cy(s) + pf5(2,jp)
      				pf5cz(s) = pf5cz(s) + pf5(3,jp)
      			endif

      		endif
      	enddo
      enddo

c     === STATS AVERAGING AND ZEROING OF COUNTING VARIABLES ===
      do s=1,ldim*40
      	pvxs(s) = pvxs(s)+(pvxc(s)/max(1e-30,pnc(s)))
      	pvxc(s) = 0
      	pvys(s) = pvys(s)+(pvyc(s)/max(1e-30,pnc(s)))
      	pvyc(s) = 0
      	pvzs(s) = pvzs(s)+(pvzc(s)/max(1e-30,pnc(s)))
      	pvzc(s) = 0

      	pvx2s(s) = pvx2s(s)+(pvx2c(s)/max(1e-30,pnc(s)))
      	pvx2c(s) = 0
      	pvy2s(s) = pvy2s(s)+(pvy2c(s)/max(1e-30,pnc(s)))
      	pvy2c(s) = 0
      	pvz2s(s) = pvz2s(s)+(pvz2c(s)/max(1e-30,pnc(s)))
      	pvz2c(s) = 0

      	if(IFORCESTAT.EQ.1) then

      		pf1sx(s) = pf1sx(s) + (pf1cx(s)/max(1e-30,pnc(s)))
      		pf1cx(s) = 0.0D0
      		pf1sy(s) = pf1sy(s) + (pf1cy(s)/max(1e-30,pnc(s)))
      		pf1cy(s) = 0.0D0
      		pf1sz(s) = pf1sz(s) + (pf1cz(s)/max(1e-30,pnc(s)))
      		pf1cz(s) = 0.0D0

      		pf2sx(s) = pf2sx(s) + (pf2cx(s)/max(1e-30,pnc(s)))
      		pf2cx(s) = 0.0D0
      		pf2sy(s) = pf2sy(s) + (pf2cy(s)/max(1e-30,pnc(s)))
      		pf2cy(s) = 0.0D0
      		pf2sz(s) = pf2sz(s) + (pf2cz(s)/max(1e-30,pnc(s)))
      		pf2cz(s) = 0.0D0

      		pf3sx(s) = pf3sx(s) + (pf3cx(s)/max(1e-30,pnc(s)))
      		pf3cx(s) = 0.0D0
      		pf3sy(s) = pf3sy(s) + (pf3cy(s)/max(1e-30,pnc(s)))
      		pf3cy(s) = 0.0D0
      		pf3sz(s) = pf3sz(s) + (pf3cz(s)/max(1e-30,pnc(s)))
      		pf3cz(s) = 0.0D0

      		pf4sx(s) = pf4sx(s) + (pf4cx(s)/max(1e-30,pnc(s)))
      		pf4cx(s) = 0.0D0
      		pf4sy(s) = pf4sy(s) + (pf4cy(s)/max(1e-30,pnc(s)))
      		pf4cy(s) = 0.0D0
      		pf4sz(s) = pf4sz(s) + (pf4cz(s)/max(1e-30,pnc(s)))
      		pf4cz(s) = 0.0D0

      		pf5sx(s) = pf5sx(s) + (pf5cx(s)/max(1e-30,pnc(s)))
      		pf5cx(s) = 0.0D0
      		pf5sy(s) = pf5sy(s) + (pf5cy(s)/max(1e-30,pnc(s)))
      		pf5cy(s) = 0.0D0
      		pf5sz(s) = pf5sz(s) + (pf5cz(s)/max(1e-30,pnc(s)))
      		pf5cz(s) = 0.0D0

      	endif

      	pncs(s) = pncs(s) + pnc(s)
      	pnc(s) = 0
      enddo

      if(IFOURWAY.EQ.1) then
c     === COLLISION STATS AVERAGING AND ZEROING OF COUNTING VARIABLES ===
      	do s=1,ldim*10
      		pxcol(s) = pxcol(s) + (pxcolc(s)/max(1e-30,avgcolc(s)))
      		pxcolc(s) = 0
      		pycol(s) = pycol(s) + (pycolc(s)/max(1e-30,avgcolc(s)))
      		pycolc(s) = 0
      		pzcol(s) = pzcol(s) + (pzcolc(s)/max(1e-30,avgcolc(s)))
      		pzcolc(s) = 0
      		avgcol(s) = avgcol(s) + avgcolc(s)
      		avgcolc(s) = 0
      	enddo
      endif

      return
      end
c-----------------------------------------------------------------------c
