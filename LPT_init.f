c======================================================================c
c---------- LAGRANGIAN PARTICLE TRACKING INITIALISATION CODE ----------c
c======================================================================c
      subroutine lpt_init

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

      common /particles_stat/ pvxc(ldim*40),pvyc(ldim*40),pvzc(ldim*40)
     $                    , pvx2c(ldim*40),pvy2c(ldim*40),pvz2c(ldim*40)
     $                    , pvxm(ldim*40),pvym(ldim*40),pvzm(ldim*40)
     $                    , pnc(ldim*40)
     $                    , pvxs(ldim*40),pvys(ldim*40),pvzs(ldim*40)
     $                    , pvx2s(ldim*40),pvy2s(ldim*40),pvz2s(ldim*40)
     $                    , pncs(ldim*40)

      common /particles_twoway/ twowayforcex(lx1,ly1,lz1,lelv)
     $                        , twowayforcey(lx1,ly1,lz1,lelv)
     $                        , twowayforcez(lx1,ly1,lz1,lelv)
     $                        ,          xmc(lx1,ly1,lz1,lelv)
     $                        ,          xxc(lx1,ly1,lz1,lelv)
     $                        ,          ymc(lx1,ly1,lz1,lelv)
     $                        ,          yxc(lx1,ly1,lz1,lelv)
     $                        ,          zmc(lx1,ly1,lz1,lelv)
     $                        ,          zxc(lx1,ly1,lz1,lelv)
     $                        ,          vol(lx1,ly1,lz1,lelv)
     $                        ,         inpc(lx1,ly1,lz1,lelv)
     $                        ,    ipartrefp(ldim*2,lpart)

      common /forces/ pf1(ldim,lpart)
     $              , pf2(ldim,lpart)
     $              , pf3(ldim,lpart)
     $              , pf4(ldim,lpart)
     $              , pf5(ldim,lpart)

      common /fourway/ avgcol(ldim*10), avgcolc(ldim*10) ! AVERAGE NUMBER OF COLLISIONS IN SLICE Y
     $               , pxcol(ldim*10), pxcolc(ldim*10) ! AVERAGE NET X MOMENTUM CHANGE IN SLICE Y
     $               , pycol(ldim*10), pycolc(ldim*10) ! AVERAGE NET Y MOMENTUM CHANGE IN SLICE Y
     $               , pzcol(ldim*10), pzcolc(ldim*10) ! AVERAGE NET Z MOMENTUM CHANGE IN SLICE Y
     $               , xmxcol, ymxcol, zmxcol
     $               , dxcol, dycol, dzcol

      common /statswitches/ I3DPLOT, ISTATGET, ISTATFLUSH, IFORCESTAT
     $                    , IWALLCONC

      common /wallconc/ wallcount, wallavg

      character*80 skipline
      character*128 filename_info

 101  FORMAT (A12,F9.3)
 102  FORMAT (A12,I5.1)
 103  FORMAT (A12,F9.3,F9.3)
 104  FORMAT (A12,E10.3)

      pi = 4.0D0*atan(1.0)

      if(nid.EQ.0) write(*,*) "================================"
      if(nid.EQ.0) write(*,*) "FIRST TIMESTEP. INITIALISING LPT"
      if(nid.EQ.0) write(*,*) "================================"

c     === ZERO PARTICLE ARRAYS ===
      do d = 1,ldim
      	do p = 1,lpart
      		ppos(d,p) = 0.0D0
      		ppos0(d,p) = 0.0D0
      		pvel(d,p) = 0.0D0
      		pvel0(d,p) = 0.0D0
      	enddo
      enddo
c     === END OF ZERO PARTICLE ARRAYS ===

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
      	pvxm(d) = 0.0D0
      	pvym(d) = 0.0D0
      	pvzm(d) = 0.0D0
      	pnc(d) = 0.0D0
      	pncs(d) = 0.0D0
      	avgcol(d) = 0.0D0
      	avgcolc(d) = 0.0D0
      	pxcol(d) = 0.0D0
      	pxcolc(d) = 0.0D0
      	pycol(d) = 0.0D0
      	pycolc(d) = 0.0D0
      	pzcol(d) = 0.0D0
      	pzcolc(d) = 0.0D0
      enddo
c     === END OF ZERO STATS ARRAYS ===

c     === ZERO TWO-WAY ARRAYS ===
      if(ITWOWAY.EQ.1) then
      	n = nx1*ny1*nz1*nelv
      	call izero(inpc,n)
      	call rzero(twowayforcex,n)
      	call rzero(twowayforcey,n)
      	call rzero(twowayforcez,n)
      	call rzero(xmc,n)
      	call rzero(xxc,n)
      	call rzero(ymc,n)
      	call rzero(yxc,n)
      	call rzero(zmc,n)
      	call rzero(zxc,n)
      	call rzero(vol,n)
      endif
c     === END OF ZEROING TWO WAY ARRAYS ===

c     === ZERO FORCE ARRAYS ===
      do d = 1,ldim
      	do p = 1,lpart
      		pf1(d,p) = 0.0D0
      		pf2(d,p) = 0.0D0
      		pf3(d,p) = 0.0D0
      		pf4(d,p) = 0.0D0
      		pf5(d,p) = 0.0D0
      	enddo
      enddo
c     === END OF ZEROING FORCE ARRAYS ===

c     === ZERO PARTICLE REFERENCE ARRAYS ===
      do d = 1,2*ldim
      	do p = 1,lpart
      		ipartrefp(d,p) = 0
      	enddo
      enddo
c     === END OF ZEROING PARTICLE REFERENCE ARRAYS ===

c     === WALL CONCENTRATION ARRAYS ===
      wallcount = 0.0D0
      wallavg = 0.0D0
c     === END OF ZEROING WALL CONCENTRATION ARRAYS ===

c     === READ IN INPUT FILE ===
      if(nid.EQ.0) write(*,*) "Reading in partinput.dat..."
      open(unit=108,file="partinput.dat")

      read(108,*) skipline 
      read(108,*) skipline
      read(108,*) IPRESTART
      if(nid.EQ.0) write(*,102) "IPRESTART: ", IPRESTART
      read(108,*) IDRAG
      if(nid.EQ.0) write(*,102) "IDRAG: ", IDRAG
      read(108,*) ILIFT
      if(nid.EQ.0) write(*,102) "ILIFT: ", ILIFT
      read(108,*) IGRAV
      if(nid.EQ.0) write(*,102) "IGRAV: ", IGRAV
      read(108,*) IPRES
      if(nid.EQ.0) write(*,102) "IPRES: ", IPRES
      read(108,*) IMASS
      if(nid.EQ.0) write(*,102) "IMASS: ", IMASS
      read(108,*) ITWOWAY
      if(nid.EQ.0) write(*,102) "ITWOWAY: ", ITWOWAY
      read(108,*) IFOURWAY
      if(nid.EQ.0) write(*,102) "IFOURWAY: ", IFOURWAY
      read(108,*) IAGGLOM
      if(nid.EQ.0) write(*,102) "IAGGLOM: ", IAGGLOM
      read(108,*) skipline
      read(108,*) srey
      if(nid.EQ.0) write(*,101) "srey: ", srey
      read(108,*) bwall(1), bcwall(1)
      read(108,*) bwall(2), bcwall(2)
      read(108,*) bwall(3), bcwall(3)
      read(108,*) bwall(4), bcwall(4)
      read(108,*) bwall(5), bcwall(5)
      read(108,*) bwall(6), bcwall(6)
      if(nid.EQ.0) write(*,103) "XMIN: ", bwall(1), bcwall(1)
      if(nid.EQ.0) write(*,103) "XMAX: ", bwall(2), bcwall(2) 
      if(nid.EQ.0) write(*,103) "YMIN: ", bwall(3), bcwall(3) 
      if(nid.EQ.0) write(*,103) "YMAX: ", bwall(4), bcwall(4) 
      if(nid.EQ.0) write(*,103) "ZMIN: ", bwall(5), bcwall(5) 
      if(nid.EQ.0) write(*,103) "ZMAX: ", bwall(6), bcwall(6)  
      read(108,*) skipline
      read(108,*) diam_p
      if(nid.EQ.0) write(*,101) "diam_p: ", diam_p
      read(108,*) dens_pf
      if(nid.EQ.0) write(*,101) "dens_pf: ", dens_pf
      read(108,*) skipline
      read(108,*) I3DPLOT
      if(nid.EQ.0) write(*,102) "I3DPLOT: ", I3DPLOT
      read(108,*) ISTATGET
      if(nid.EQ.0) write(*,102) "ISTATGET: ", ISTATGET
      read(108,*) ISTATFLUSH
      if(nid.EQ.0) write(*,102) "ISTATFLUSH: ", ISTATFLUSH
      read(108,*) IFORCESTAT
      if(nid.EQ.0) write(*,102) "IFORCESTAT: ", IFORCESTAT
      read(108,*) IWALLCONC
      if(nid.EQ.0) write(*,102) "IWALLCONC: ", IWALLCONC

      close(108)
c     === END OF READING IN INPUT FILE ===

c     === CALCULATE OTHER VALUES FROM INPUT FILES ===
      brey = param(2)**-1
      if(nid.EQ.0) write(*,101) "brey: ", brey
      dt_p = -1.0D0*param(12)
      if(nid.EQ.0) write(*,101) "dt_p: ", dt_p
      rad_p = 0.5D0*diam_p
      if(nid.EQ.0) write(*,104) "rad_p: ", rad_p
      vol_p = (4.0/3.0)*pi*rad_p**3
      if(nid.EQ.0) write(*,104) "vol_p: ", vol_p
      stkb_p = (diam_p*diam_p*brey*dens_pf)/(18.0D0)
      if(nid.EQ.0) write(*,101) "stkb_p: ", stkb_p
      stkt_p = stkb_p*srey*srey/brey
      if(nid.EQ.0) write(*,101) "stkt_p: ", stkt_p
      vfrac_p = (lpart*vol_p)/((2.0D0*bwall(2))
     $                        *(2.0D0*bwall(4))
     $                        *(2.0D0*bwall(6)))
      if(nid.EQ.0) write(*,104) "vfrac_p: ", vfrac_p
c     === END OF CALCULATING OTHER VALUES ===

c     === GEOMETRY PARAMETERS FOR FOURWAY ===
      if(IFOURWAY.EQ.1) then
      	xmxcol = 2.0D0*bwall(2)
      	ymxcol = 2.0D0*bwall(4)
      	zmxcol = 2.0D0*bwall(6)
      	dxcol = xmxcol/(10.0D0*ldim)
      	dycol = ymxcol/(10.0D0*ldim)
      	dzcol = zmxcol/(10.0D0*ldim)
      	if(nid.EQ.0) write(*,101) "dxcol: ", dxcol
      	if(nid.EQ.0) write(*,101) "dycol: ", dycol
      	if(nid.EQ.0) write(*,101) "dzcol: ", dzcol
      endif
c     === END OF GEOMETRY PARAMETERS FOR FOURWAY ===

c     === MAKE A PARTICLE INFO FILE ===
      open(unit = 74,file = 'partstat/partinfo.dat')

      if(nid.EQ.0) write(74,*) "SWITCHES"
      if(nid.EQ.0) write(74,*) "========"
      if(nid.EQ.0) write(74,*) "  "

      if(nid.EQ.0) write(74,102) "IPRESTART: ", IPRESTART
      if(nid.EQ.0) write(74,102) "IDRAG: ", IDRAG
      if(nid.EQ.0) write(74,102) "ILIFT: ", ILIFT
      if(nid.EQ.0) write(74,102) "IGRAV: ", IGRAV
      if(nid.EQ.0) write(74,102) "IPRES: ", IPRES
      if(nid.EQ.0) write(74,102) "IMASS: ", IMASS
      if(nid.EQ.0) write(74,102) "ITWOWAY: ", ITWOWAY
      if(nid.EQ.0) write(74,102) "IFOURWAY: ", IFOURWAY
      if(nid.EQ.0) write(74,102) "IAGGLOM: ", IAGGLOM

      if(nid.EQ.0) write(74,*) "  "
      if(nid.EQ.0) write(74,*) "FLUID DATA"
      if(nid.EQ.0) write(74,*) "=========="
      if(nid.EQ.0) write(74,*) "  "

      if(nid.EQ.0) write(74,101) "srey: ", srey
      if(nid.EQ.0) write(74,101) "brey: ", brey
      if(nid.EQ.0) write(74,103) "XMIN: ", bwall(1), bcwall(1)
      if(nid.EQ.0) write(74,103) "XMAX: ", bwall(2), bcwall(2)
      if(nid.EQ.0) write(74,103) "YMIN: ", bwall(3), bcwall(3)
      if(nid.EQ.0) write(74,103) "YMAX: ", bwall(4), bcwall(4)
      if(nid.EQ.0) write(74,103) "ZMIN: ", bwall(5), bcwall(5)
      if(nid.EQ.0) write(74,103) "ZMAX: ", bwall(6), bcwall(6)

      if(nid.EQ.0) write(74,*) "  "
      if(nid.EQ.0) write(74,*) "PARTICLE DATA"
      if(nid.EQ.0) write(74,*) "============="
      if(nid.EQ.0) write(74,*) "  "

      if(nid.EQ.0) write(74,101) "diam_p: ", diam_p
      if(nid.EQ.0) write(74,101) "dens_pf: ", dens_pf
      if(nid.EQ.0) write(74,101) "dt_p: ", dt_p
      if(nid.EQ.0) write(74,104) "rad_p: ", rad_p
      if(nid.EQ.0) write(74,104) "vol_p: ", vol_p
      if(nid.EQ.0) write(74,101) "stkb_p: ", stkb_p
      if(nid.EQ.0) write(74,101) "stkt_p: ", stkt_p
      if(nid.EQ.0) write(74,104) "vfrac_p: ", vfrac_p
      close(74)

      return
      end
c-----------------------------------------------------------------------c
