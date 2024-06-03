      subroutine setprob
      implicit double precision (a-h,o-z)
      character*25 fname

      common /combc/ omega
      common /comwall/ pi,t1,a1,tw1
      common /combcu/ ubc(1:1000)
      common /combcn/ nt,ntmax
      common /comtimereverse/ trtime, trdone
      common /commat/ dKA, dKB, rhoA, rhoB
      common /comic/ ic, a2
      common /comchar/ xchar(1000),x0char,t0char,dtchar,nchar

c
c     # Set the material parameters for the acoustic equations
c
      iunit = 7
      fname = 'setprob.data'
      
      call opendatafile(iunit, fname)

      pi = 4.d0*datan(1.d0)

c     # parameters for wall
      read(7,*) t1
      read(7,*) a1
      read(7,*) tw1
      read(7,*) dKA
      read(7,*) dKB
      read(7,*) rhoA
      read(7,*) rhoB
      read(7,*) trtime
      read(7,*) ic
      read(7,*) a2
      read(7,*) nchar
      read(7,*) x0char
      read(7,*) t0char
      read(7,*) dtchar

      if (nchar.gt.1000) then
         write(6,*) 'nchar limited to 1000'
         stop
         endif

      do k=1,nchar
         xchar(k) = x0char
         enddo

      return
      end
