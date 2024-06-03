      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname

      common /comwall/ pi,t1,a1,tw1
      common /commat/ dkappaA, dkappaB

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
      read(7,*) dkappaA
      read(7,*) dkappaB

      return
      end
