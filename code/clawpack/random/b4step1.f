c     ============================================
      subroutine b4step1(maxmx,mbc,mx,meqn,q,
     &            xlower,dx,t,dt,maux,aux)
c     ============================================
c
c     # called from claw1 before each call to step1.
c     # use to set time-dependent aux arrays or perform other tasks
c     # which must be done every time step.
c
c
c     
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, *)
      common /comtimereverse/ trtime, trdone

      if (t.ge.trtime .and. trdone.eq.0.d0) then
        do i=1-mbc,mx+mbc
          q(i,2)=-q(i,2)
        enddo
        trdone=1.d0
      endif

      i = 1
      sig = sigma(q(i,1),i,aux(i,2),aux(i,3))
      u = q(i,2)/aux(i,1)
      write(17,1005) t,q(i,1),q(i,2),sig,u
 1005 format(e16.8,4e16.8)

      i = mx
      sig = sigma(q(i,1),i,aux(i,2),aux(i,3))
      u = q(i,2)/aux(i,1)
      write(18,1005) t,q(i,1),q(i,2),sig,u

      i = 12
      sig = sigma(q(i,1),i,aux(i,2),aux(i,3))
      u = q(i,2)/aux(i,1)
      write(28,1005) t,q(i,1),q(i,2),sig,u

c
      return
      end
