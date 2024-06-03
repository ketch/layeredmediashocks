      subroutine rp1(maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                  fwave,s,amdq,apdq)
c     =====================================================
c
c     # Riemann solver for the p-system equations in 1d,
c     #  with spatially varying pressure function
c     #   v_t - u_x = 0
c     #   u_t + p(v,x)_x =0
c     # where v=specific volume, u=velocity
c
c     # Lagrangian isothermal flow: p(v) = a(x)/v^2
c     # Lagrangian isentropic flow: p(v) = kappa(x)/v^gamma
c     # aux(i,1) = kappa(i)
c
c     # function pressure(v,i) gives pressure relation in ith cell
c     # function pressurep(v,i) gives dp/dv
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # On output, fwave contains the waves as jumps in f,
c     #            s the speeds,
c     #
c     #            amdq = A^- Delta q, 
c     #            apdq = A^+ Delta q,
c     #                   the decomposition of the flux difference
c     #                       f(qr(i-1)) - f(ql(i))
c     #                   into leftgoing and rightgoing parts respectively.
c     #
c
c     # Note that the ith Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c
      implicit double precision (a-h,o-z)
c
      dimension auxl(1-mbc:maxm+mbc, 1)
      dimension auxr(1-mbc:maxm+mbc, 1)
      dimension fwave(1-mbc:maxm+mbc, meqn, mwaves)
      dimension    s(1-mbc:maxm+mbc, mwaves)
      dimension   ql(1-mbc:maxm+mbc, meqn)
      dimension   qr(1-mbc:maxm+mbc, meqn)
      dimension apdq(1-mbc:maxm+mbc, meqn)
      dimension amdq(1-mbc:maxm+mbc, meqn)
c
c
c     # split the jump in q at each interface into waves
c
      do 20 i = 2-mbc, mx+mbc
         vi  = ql(i,1)
         vim = qr(i-1,1)
         ui  = ql(i,2)
         uim = qr(i-1,2)

c        #linearize on each side:

         ppi = -pressurep(vi,auxl(i,1))
         ppim = -pressurep(vim,auxr(i-1,1))
         ci = dsqrt(ppi)
         cim = dsqrt(ppim)
         zi = ci
         zim = cim

         du = ui- uim
         dp = pressure(vi,auxl(i,1))
     &          - pressure(vim,auxr(i-1,1))
         b1 = -(zi*du - dp) / (zim + zi)
         b2 = -(zim*du + dp) / (zim + zi)
c
c        # Compute the waves.
c
         fwave(i,1,1) = b1
         fwave(i,2,1) = b1 * zim
         s(i,1) = -cim
c
         fwave(i,1,2) = b2
         fwave(i,2,2) = b2*(-zi)
         s(i,2) = ci

   20    continue
c
c     # compute the leftgoing and rightgoing fluctuations:
c     # Note s(i,1) < 0   and   s(i,2) > 0.
c
      do 220 m=1,meqn
         do 220 i = 2-mbc, mx+mbc
            amdq(i,m) = fwave(i,m,1)
            apdq(i,m) = fwave(i,m,2)
  220       continue
c
      return
      end


c     --------------------------------------------
      double precision function pressure(v,dkappa)
c     --------------------------------------------
      implicit double precision (a-h,o-z)

c     # pressure relation in ith cell
c     # p=kappa/v**2

      pressure = dkappa / v**2

      return
      end


c     --------------------------------------------
      double precision function pressurep(v,dkappa)
c     --------------------------------------------
      implicit double precision (a-h,o-z)

c     # derivative of pressure relation in ith cell
c     # dp/dv = -2*kappa/v**3

      pressurep = -2*dkappa/v**3

      return
      end
