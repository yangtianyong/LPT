c=======================================================================c
c---------- LAGRANGIAN PARTICLE TRACKING TWOWAY COUPLING CODE ----------c
c=======================================================================c
      subroutine lpt_twoway(ITWO)

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

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

      common /particles_info/ ppos(ldim,lpart), ppos0(ldim,lpart)
     $                      , pvel(ldim,lpart), pvel0(ldim,lpart)

      common /particles_prop/ diam_p, dt_p, dens_pf, rad_p, vol_p
     $                      , stkb_p, stkt_p, vfrac_p

      integer ITWO
      integer n
      integer nparts
      real gvol

      n = nx1*ny1*nz1*nelv

c     === SET UP TWO-WAY COUPLING ===
      if(ITWO.EQ.0) then
c     	=== FIND CELL BOUNDARIES AND VOLUMES ===
      	do e = 1, nelv
      		do i = 1, nx1
      			do j = 1, ny1
      				do k = 1, nz1

      				if(i.GE.2) then
      					xmc(i,j,k,e) = xm1(i,j,k,e) - 
     $         (xm1(i,j,k,e)-xm1(i-1,j,k,e))*0.5D0
      				else
      					xmc(i,j,k,e)  = xm1(i,j,k,e)
      				endif

      				if(i.LE.7) then
      					xxc(i,j,k,e) = xm1(i,j,k,e) +
     $         (xm1(i+1,j,k,e)-xm1(i,j,k,e))*0.5D0
      				else
      					xxc(i,j,k,e) = xm1(i,j,k,e)
      				endif

      				if(j.GE.2) then
      					ymc(i,j,k,e) = ym1(i,j,k,e) - 
     $         (ym1(i,j,k,e)-ym1(i,j-1,k,e))*0.5D0
      				else
      					ymc(i,j,k,e) = ym1(i,j,k,e)
      				endif

      				if(j.LE.7) then
      					yxc(i,j,k,e) = ym1(i,j,k,e) +
     $         (ym1(i,j+1,k,e)-ym1(i,j,k,e))*0.5D0
      				else
      					yxc(i,j,k,e) = ym1(i,j,k,e)
      				endif

      				if(k.GE.2) then
      					zmc(i,j,k,e) = zm1(i,j,k,e) - 
     $     (zm1(i,j,k,e)-zm1(i,j,k-1,e))*0.5D0
      				else
      					zmc(i,j,k,e) = zm1(i,j,k,e)
      				endif

      				if(k.LE.7) then
      					zxc(i,j,k,e) = zm1(i,j,k,e) +
     $         (zm1(i,j,k+1,e)-zm1(i,j,k,e))*0.5D0
      				else
      					zxc(i,j,k,e) = zm1(i,j,k,e)
      				endif

      				vol(i,j,k,e) =
     $        (xxc(i,j,k,e) - xmc(i,j,k,e))
     $       *(yxc(i,j,k,e) - ymc(i,j,k,e))
     $       *(zxc(i,j,k,e) - zmc(i,j,k,e))
  
      				enddo
      			enddo
      		enddo
      	enddo

      	gvol = glsum(vol,n)
      	if(nid.eq.0) write(*,*) 'Real volume: ', volvm1
      	if(nid.eq.0) write(*,*) 'Calculated volume: ', gvol
      endif

c     === TWO-WAY COUPLING ===
      if(ITWO.EQ.1) then
c       === RE-ZERO FORCE ARRAYS ===
      	call rzero(twowayforcex,n)
      	call rzero(twowayforcey,n)
      	call rzero(twowayforcez,n)
      	call izero(inpc,n)

c       === LOOP OVER ALL THE PARTICLES ===
      	do jp=1,lpart
      		ix = ipartrefp(1,jp)
      		iy = ipartrefp(2,jp)
      		iz = ipartrefp(3,jp)
      		iel = ipartrefp(4,jp)

c         === CHECK TO SEE IF PARTICLE HAS BEEN FOUND PREVIOUSLY
      		if(ipartrefp(5,jp) .EQ. 1 .AND.
     $       ppos0(1,jp) .GE. xmc(ix,iy,iz,iel) .AND.
     $       ppos0(1,jp) .LT. xxc(ix,iy,iz,iel) .AND.
     $       ppos0(2,jp) .GE. ymc(ix,iy,iz,iel) .AND.
     $       ppos0(2,jp) .LT. yxc(ix,iy,iz,iel) .AND.
     $       ppos0(3,jp) .GE. zmc(ix,iy,iz,iel) .AND.
     $       ppos0(3,jp) .LT. zxc(ix,iy,iz,iel)) then

c           === PARTICLE HAS BEEN FOUND ===
      			inpc(ix,iy,iz,iel) = inpc(ix,iy,iz,iel) + 1 ! Update number of particles in cell
      			twowayforcex(ix,iy,iz,iel) =
     $			twowayforcex(ix,iy,iz,iel) + pf1(1,jp) ! Update x force for particle
      			twowayforcey(ix,iy,iz,iel) =
     $			twowayforcey(ix,iy,iz,iel) + pf1(2,jp) ! Update y force for particle
      			twowayforcez(ix,iy,iz,iel) =
     $			twowayforcez(ix,iy,iz,iel) + pf1(3,jp) ! Update z force for particle
          else
            do iel = 1,nelv ! Loop over all elements
      				do ix = 1,nx1 ! Loop over all x cells
      					if(ppos0(1,jp) .GE. xmc(ix,iy,iz,iel) .AND.
     $					ppos0(1,jp) .LT. xxc(ix,iy,iz,iel)) then
      						ipartrefp(1,jp) = ix
      						do iy = 1,ny1 ! Loop over all y cells
      							if(ppos0(2,jp) .GE. ymc(ix,iy,iz,iel) .AND.
     $							ppos0(2,jp) .LT. yxc(ix,iy,iz,iel)) then
      								ipartrefp(2,jp) = iy
      								do iz = 1,nz1 ! Loop over all z cells
      									if(ppos0(3,jp) .GE. zmc(ix,iy,iz,iel) .AND.
     $									ppos0(3,jp) .LT. zxc(ix,iy,iz,iel)) then
      										ipartrefp(3,jp) = iz
      										ipartrefp(4,jp) = iel
      										ipartrefp(5,jp) = 1
      										inpc(ix,iy,iz,iel) = inpc(ix,iy,iz,iel) + 1 ! Update number of particles in cell
      										twowayforcex(ix,iy,iz,iel) =
     $										twowayforcex(ix,iy,iz,iel) + pf1(1,jp) ! Update x force for particle
      										twowayforcey(ix,iy,iz,iel) =
     $										twowayforcey(ix,iy,iz,iel) + pf1(2,jp) ! Update y force for particle
      										twowayforcez(ix,iy,iz,iel) =
     $										twowayforcez(ix,iy,iz,iel) + pf1(3,jp) ! Update z force for particle
      									endif
      								enddo
      							endif
      						enddo
      					endif
      				enddo
      			enddo
      		endif
c         === END OF CHECKING 1 PARTICLE ===
      	enddo
c       === END OF LOOPING OVER ALL PARTICLES ===
      	nparts = iglsum(inpc,n)

      	if(nid.eq.0) write(*,*) 'Number of particles: ', lpart
      	if(nid.eq.0) write(*,*) 'Clc number of particles: ', nparts

      	do iel = 1,nelv
      		do ix = 1,nx1
      			do iy = 1,ny1
      				do iz = 1,nz1
      					if(inpc(ix,iy,iz,iel).GT.0) then
      						twowayforcex(ix,iy,iz,iel) = 
     $						-1.0D0*twowayforcex(ix,iy,iz,iel)
     $						*vol_p/vol(ix,iy,iz,iel)
      						twowayforcey(ix,iy,iz,iel) = 
     $						-1.0D0*twowayforcey(ix,iy,iz,iel)
     $						*vol_p/vol(ix,iy,iz,iel)
      						twowayforcez(ix,iy,iz,iel) = 
     $						-1.0D0*twowayforcez(ix,iy,iz,iel)
     $						*vol_p/vol(ix,iy,iz,iel)
      					endif
      				enddo
      			enddo
      		enddo
      	enddo

      endif

      return
      end
c-----------------------------------------------------------------------c
