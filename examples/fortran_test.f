!     Example of fortran code calling the c++ library
!     compile with: gfortran fortran_test.f `qcdloop-config --ldflags`
      
      program main
      implicit none
      double complex qli1
      double precision m, mu2

!     init qcdloop objects
      call qlinit()

      m = 5d0
      mu2 = 1.7d0**2     
!     call some test function
      write(*,*) qli1(m, mu2, 0)
      write(*,*) qli1(m, mu2, 1)
      write(*,*) qli1(m, mu2, 2)
      
      end
