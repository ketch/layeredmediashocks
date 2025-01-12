c =========================================================
      subroutine output(meqn,mbc,mx,xlower,dx,q,t,iframe,aux,maux)
c =========================================================
c
c     # Output the results for a general system of conservation laws
c     # in 1 dimension
c
      implicit double precision (a-h,o-z)
      parameter (ndim=1)
      dimension q(1-mbc:mx+mbc, meqn)
      dimension aux(1-mbc:mx+mbc, 3)
      character*10 fname1, fname2, fname3
      logical outaux

      outaux = .false.
c
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw1.m
c     # The same format is used by the amrclaw package.  
c     # Here it's adapted to output just the single grid.
c
c     # first create the file name and open file
c
         fname1 = 'fort.qxxxx'
         fname2 = 'fort.txxxx'
         fname3 = 'fort.axxxx'
         nstp = iframe
         do 55 ipos = 10, 7, -1
            idigit = mod(nstp,10)
            fname1(ipos:ipos) = char(ichar('0') + idigit)
            fname2(ipos:ipos) = char(ichar('0') + idigit)
            fname3(ipos:ipos) = char(ichar('0') + idigit)
            nstp = nstp / 10
 55      continue

         open(unit=50,file=fname1,status='unknown',form='formatted')
         open(unit=60,file=fname2,status='unknown',form='formatted')

c
c     # the following parameters are used in amrclaw where there are
c     # multiple grids.  Here they are all set to 1:
      ngrids = 1
      mptr = 1
      level = 1

      write(50,1001) mptr,level,mx
 1001 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 mx')

      write(50,1002) xlower,dx
 1002 format(e18.8,'    xlow', /,
     &       e18.8,'    dx', /)
c
        do 10 i=1-mbc,mx+mbc
          do m=1,meqn
c            # exponents with more than 2 digits cause problems reading
c            # into matlab... reset tiny values to zero:
             if (dabs(q(i,m)) .lt. 1d-90) q(i,m) = 0.d0
             enddo
c
          sig = sigma(q(i,1),i,aux(i,2),aux(i,3))
          u = q(i,2)/aux(i,1)
c          write(50,1005) (q(i,m), m=1,meqn),sig,u
          write(50,1005) sig,u
 1005     format(4e32.16)
c
 10       continue
        write(50,*) ' '
 20     continue
      write(50,*) ' '

c      write(60,1000) t,meqn+2,ngrids,ndim,maux
      write(60,1000) t,meqn,ngrids,maux,1
 1000 format(e18.8,'    time', /, 
     &       i5,'                 meqn'/,
     &       i5,'                 ngrids'/,
     &       i5,'                 maux'/,
     &       i5,'                 ndim'/,/)
c

      close(unit=50)
      close(unit=60)

      return
      end
