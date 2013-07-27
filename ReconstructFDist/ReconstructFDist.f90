!  ReconstructFDist.f90 
!
!  FUNCTIONS:
!  ReconstructFDist - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: ReconstructFDist
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program ReconstructFDist
    use MathUtils, ONLY : logspace_s
    implicit none
    integer, parameter :: N=100
    real*8 R(N), r1/2.0/, r2/100.0/
    R=0
    call logspace_s(r1, r2, R, N-10)
  
    print *,R
    end program ReconstructFDist

