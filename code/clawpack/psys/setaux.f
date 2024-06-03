!     ============================================
      subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)
!     ============================================
!
!     # set auxiliary arrays 
!     # variable coefficient acoustics
!     #  aux(i,1) = kappa in ith cell
!
      implicit double precision (a-h,o-z)
      dimension aux(1-mbc:maxmx+mbc, 1)
      common /commat/ dkappaA, dkappaB

      open(unit=31,file='fort.aux',status='unknown',form='formatted')
!

       xupper = xlower + mx*dx
       do 100 i=1-mbc,mx+mbc
          xcell = xlower + (i-0.5d0)*dx
!
!         # discrete layers:
!         # layer A between j and j+frac1
!         # layer B between j+frac1 and j+1

          frac1 = 1.d0/2.d0

          ix = xcell
          xfrac = xcell - ix
          if (xfrac .lt. frac1) then
!             # layer A:
              aux(i,1) = dkappaA
            else
!             # layer B:
              aux(i,1) = dkappaB
            endif

  100     continue

!      # make material uniform in ghost cells:
       do ibc=1,mbc
          aux(1-ibc,1) = aux(1,1)
          aux(mx+ibc,1) = aux(mx,1)
          enddo

	do i=1-mbc,mx+mbc
          write(31,701) aux(i,1)
  701     format(3e16.6)
          enddo

       close(unit=31)
!
       return
       end
