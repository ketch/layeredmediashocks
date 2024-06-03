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
      common /comic/ ic, a2
      common /commat/ dKA, dKB, rhoA, rhoB
c
c
      rhomean = 0.5d0*(rhoA+rhoB)
      dKmean  = 2.d0/(1.d0/dKA + 1.d0/dKB)
      Zmean = dsqrt(rhomean*dKmean)

      do 150 i=1,mx
        xcell = xlower + (i-0.5d0)*dx

        if(ic==1) then ! Zero IC
            q(i,1) = 0.d0
            q(i,2) = 0.d0

        elseif(ic==2) then !Gaussian IC
            q(i,2) = 0.d0!-a2*dexp(-((xcell-75.d0)/10.d0)**2.d0)*aux(i,1)
            sig = a2*dexp(-((xcell-75.d0)/30.d0)**2.d0)
            q(i,1) = dlog(sig+1)/aux(i,2)

        elseif(ic==3) then ! entropy tests:
            q(i,2) = -a2*dexp(-((xcell-75.d0)/10.d0)**2.d0)*aux(i,1)
            sig = (1.d0-0.5d0*Zmean*q(i,2)/aux(i,1))**2.d0-1.d0
            q(i,1) = dlog(sig+1)/aux(i,2)
        endif
  150    continue
c
      return
      end
