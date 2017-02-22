c======================================================================c
c----- LAGRANGIAN PARTICLE TRACKING INTERPARTICLE COLLISION CODE ------c
c======================================================================c
      subroutine lpt_ip_collide

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      common /particles_prop/ diam_p, dt_p, dens_pf, rad_p, vol_p
     $                      , stkb_p, stkt_p, vfrac_p

      common /particles_info/ ppos(ldim,lpart), ppos0(ldim,lpart)
     $                      , pvel(ldim,lpart), pvel0(ldim,lpart)

      common /fourway/ avgcol(ldim*10), avgcolc(ldim*10) ! AVERAGE NUMBER OF COLLISIONS IN SLICE Y
     $               , pxcol(ldim*10), pxcolc(ldim*10) ! AVERAGE NET X MOMENTUM CHANGE IN SLICE Y
     $               , pycol(ldim*10), pycolc(ldim*10) ! AVERAGE NET Y MOMENTUM CHANGE IN SLICE Y
     $               , pzcol(ldim*10), pzcolc(ldim*10) ! AVERAGE NET Z MOMENTUM CHANGE IN SLICE Y
     $               , xmxcol, ymxcol, zmxcol
     $               , dxcol, dycol, dzcol

      common /flow_prop/ brey, srey, bwall(ldim*2), bcwall(ldim*2)
     $                 , fvel(ldim,lpart)
     $              , fdvu(ldim,lpart),fdvv(ldim,lpart),fdvw(ldim,lpart)

      integer npref(ldim*10,ldim*10,ldim*10,ldim*20) ! REFERENCE TO NPSth PARTICLE IN SLICE X,Y,Z
      integer nps(ldim*10,ldim*10,ldim*10) ! NUMBER OF PARTICLES IN SLICE X,Y,Z
      integer npsx(lpart), npsy(lpart), npsz(lpart) ! SLICE PARTICLE J IS CONTAINED IN
      integer ncc1, ncc2 ! REFERENCE TO COLLIDING PARTICLE 
      integer ipx, ipy, ipz, nxtp ! PARTICLE REFERENCES
      real nfull(ldim), nhat(ldim) ! RELATIVE VELOCITY VECTOR (AND UNIT VECTOR)
      real vfull(ldim), vhat(ldim) ! RELATIVE VELOCITY VECTOR (AND UNIT VECTOR)
      real vmag1, vmag2, tcol1, tcol2, dpart ! VELOCITY MAGNITUDE, COLLISION TIME

c     === RESET NPS ===
      do sx = 1, 10*ldim
      	do sy = 1, 10*ldim
      		do sz = 1, 10*ldim
      			nps(sx,sy,sz) = 0
      		enddo
      	enddo
      enddo
c     === END OF RESET NPS ===

c     === COUNT THE NUMBER OF PARTICLES IN EACH SLICE ===
      do jp = 1, lpart
      	ipx = int((ppos(1,jp)+bwall(2))/dxcol) + 1
      	ipy = int((ppos(2,jp)+bwall(4))/dycol) + 1
      	ipz = int((ppos(3,jp)+bwall(6))/dzcol) + 1

      	if(ipx.GT.10*ldim) ipx = 10*ldim
      	if(ipy.GT.10*ldim) ipy = 10*ldim
      	if(ipz.GT.10*ldim) ipz = 10*ldim

      	nps(ipx,ipy,ipz) = nps(ipx,ipy,ipz) + 1

      	npsx(jp) = ipx
      	npsy(jp) = ipy
      	npsz(jp) = ipz

      	npref(ipx,ipy,ipz,nps(ipx,ipy,ipz)) = jp ! Reference of NPSth particle in slice I,J,K
      enddo
c     === END OF CHECKING ALLOCATED PARTICLES ===

c     === CHECK ALL SLICES FOR COLLISIONS ===
      do sx=1,10*ldim
      	do sy=1,10*ldim
      		do sz=1,10*ldim
c     			=== CHECK SLICE FOR COLLISION
      			do c1=1,nps(sx,sy,sz) ! FIRST PARTICLE C1
      				nxtp = (c1+1)
      				do c2 = nxtp,nps(sx,sy,sz) ! CHECKING WITH SECOND PARTICLE C2
      					cc1 = npref(sx,sy,sz,c1)
      					cc2 = npref(sx,sy,sz,c2)
c     					=== CHECK CC1 AND CC2 FOR POTENTIAL COLLISION ===
      					dpart = sqrt((ppos(1,cc1)-ppos(1,cc2))**2
     $                      +(ppos(2,cc1)-ppos(2,cc2))**2
     $                      +(ppos(3,cc1)-ppos(3,cc2))**2) ! CALCULATE SEPARATION

      					if(dpart .LE. diam_p) then 

c     						=== PARTICLE HAS COLLIDED ===
      						do d=1,3
      							nfull(d) = ppos(d,cc1)-ppos(d,cc2)
      							vfull(d) = pvel(d,cc1)-pvel(d,cc2)
      						enddo

      						do d=1,3
      							nhat(d) = nfull(d)/
     $							sqrt(nfull(1)**2+nfull(2)**2+nfull(3)**2)
      							vhat(d) = vfull(d)/
     $							sqrt(vfull(1)**2+vfull(2)**2+vfull(3)**2)
      						enddo

      						vmag1 = sqrt(pvel(1,cc1)**2
     $												+pvel(2,cc1)**2
     $												+pvel(3,cc1)**2)

      						vmag2 = sqrt(pvel(1,cc2)**2
     $												+pvel(2,cc2)**2
     $												+pvel(3,cc2)**2)

      						tcol1 = (diam_p-dpart)/vmag1
      						tcol2 = (diam_p-dpart)/vmag2

      						do d=1,3
      							ppos(d,cc1)=ppos(d,cc1)+((diam_p-dpart)+1e-30)*nhat(d)
      						enddo

      						a1 = (pvel(1,cc1)*nhat(1)+
     $									pvel(2,cc1)*nhat(2)+
     $									pvel(3,cc1)*nhat(3))

      						a2 = (pvel(1,cc2)*nhat(1)+
     $									pvel(2,cc2)*nhat(2)+
     $									pvel(3,cc2)*nhat(3))

      						do d=1,3
      							pvel(d,cc1)=(pvel(d,cc1)-(a1-a2)*nhat(d))
      							pvel(d,cc2)=(pvel(d,cc2)+(a1-a2)*nhat(d))
      							ppos(d,cc1)=(ppos(d,cc1)+pvel(d,cc1)*tcol1)
      							ppos(d,cc2)=(ppos(d,cc2)+pvel(d,cc2)*tcol1)
      						enddo

      						avgcolc(sy) = avgcolc(sy) + 1.0D0
      						pxcolc(sy) = pxcolc(sy) + (a1-a2)*nhat(1)
      						pycolc(sy) = pycolc(sy) + (a1-a2)*nhat(2)
      						pzcolc(sy) = pzcolc(sz) + (a1-a2)*nhat(3)

c     						=== END OF PARTICLE COLLISION ===
      					endif
      				enddo
      			enddo
c     			=== END OF CHECKING ONE SLICE FOR COLLISION
      		enddo
      	enddo
      enddo
c     === END OF CHECKING ALL SLICES FOR COLLISIONS ===

      return
      end
c-----------------------------------------------------------------------c
