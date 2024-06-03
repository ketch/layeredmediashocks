c
c
c =========================================================
       subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, *)
c
!     Initialize with uniform pressure everywhere:
!     p = kappa(x)/v^2 = 1
!     v(x,0) = 1./sqrt(kappa(x))
      do 150 i=1,mx
	 xcell = xlower + (i-0.5d0)*dx
         q(i,1) = dsqrt(aux(i,1))
         q(i,2) = 0.d0
         p = pressure(q(i,1),aux(i,1))
  150    continue
c
      return
      end
