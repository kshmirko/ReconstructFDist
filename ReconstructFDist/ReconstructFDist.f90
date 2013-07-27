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
    use MathUtils, ONLY : linspace_s, basisMatrix, calcCoefs
    

    implicit none
    integer, parameter :: N=100, M=5
    integer I
    real*8 R(N), Y(N), r1/1.0/, r2/6.0/, KNTS(M), B(N,M), YY(N,1), C(M)


   
    
    KNTS = 0
    R=0
    call linspace_s(r1, r2, R, N)
    call linspace_s(r1, r2-1.0, KNTS,M)
    call basisMatrix(B, R, KNTS, N,M)
    

    Y = sin(R*2*3.1415/180.0)
    do I=1,N
        write(*,'(7F8.3)') R(I),Y(I),B(I,:)
    enddo

    
    YY = RESHAPE(Y, (/N,1/))
    
    call calcCoefs(B,YY,C,N,M)
    
    do I=1,M
        print*, knts(I),C(I)
    enddo
    end program ReconstructFDist

