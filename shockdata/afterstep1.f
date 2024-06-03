!     ============================================
      subroutine afterstep1(maxmx,mbc,mx,meqn,q,
     &            xlower,dx,t,dt,maux,aux)
!     ============================================
!
!     # called from claw1 after successful step.
!     # use to set time-dependent aux arrays or perform other tasks
!     # which must be done every time step.
!
!
!     
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, *)
      common /comchar/ xchar(1000),x0char,t0char,dtchar,nchar


      !----------------------------------------------------
      ! advance ODEs for characteristic curves:

      ! number of characteristics being actively tracked:
      nchar_track = min(max(0, floor((t-t0char)/dtchar)), nchar)
      ! write(66,*) t,nchar_track

!     write(66,*) "In afterstep at t = ",t
      do k=1,nchar_track
          xc = xchar(k)
          i0 = floor((xc - xlower)/dx - 0.5d0)
          xi0 = xlower + (i0-0.5d0)*dx

          ! linear interpolation of characteristic velocity to xc:
          w1 = (xc-xi0)/dx
          w2 = 1.d0 - w1

          bulk1 = sigmap(q(i0+1,1),i0+1,aux(i0+1,2),aux(i0+1,3))
          rho1 = aux(i0+1,1)
          c1 = dsqrt(bulk1/rho1)

          bulk2 = sigmap(q(i0,1),i0,aux(i0,2),aux(i0,3))
          rho2 = aux(i0,1)
          c2 = dsqrt(bulk2/rho2)

          vc = w1*c1 + w2*c2
!         write(66,666) k,i0,xc,c1,c2
!         write(66,667) rho1,bulk1,rho2,bulk2
! 666     format(/,'#',i3,i5,3e16.6)
! 667     format(10x,4e16.6)

          ! advance location xc by Euler:
          xchar(k) = xc + dt*vc

          enddo
      write(31,*) t,xchar(1:nchar)

      sigmax = 0.d0
      dsigmax = 0.d0
      imax = 0
      idmax = 0
      sigold = sigma(q(0,1),0,aux(0,2),aux(0,3))
      entropy = 0.d0
      do i=1,mx
         sig = sigma(q(i,1),i,aux(i,2),aux(i,3))
         dsig = sig - sigold
         sigold = sig
         if (sig .gt. sigmax) then
             sigmax = sig
             imax = i
         endif
         if (dsig .gt. dsigmax) then
             dsigmax = dsig
             idmax = i
         endif
         entropy = entropy + 0.5d0*q(i,2)**2 / aux(i,1) 
     &               + sig/aux(i,2) - q(i,1)
          
      enddo
      xmax = xlower + (imax-0.5d0)*dx
      xdmax = xlower + (idmax-0.5d0)*dx
      write(68,681) t,xmax,xdmax
      entropy = entropy * dx
      write(69,681) t,entropy
  681 format(3e16.6)
!     write(66,*) 'imax = ',imax

          

      return
      end
