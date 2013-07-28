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
    use MathUtils, ONLY : DECOMPOSE_VECTOR_S, LINSPACE_S, DECOMPOSE_MATRIX_S
    use mieScatter, only  : calcScatt
    

    implicit none
    integer, parameter :: N=100, M=5
    integer I
    real*8 R0(N), X(N), Y(N,M), r1/0.05/, r2/2.0/, KNTS(M), B(N,M), YY(N,1), C(M,M), K(N,M)
    real*8 QEXT(N),QSCA(N),qabs(N),qbk(N),qpr(N),alb(N),g(N)
    complex*16 RI/(1.33,-0.005)/
       
    KNTS = 0
    R0=0
    call linspace_s(r1, r2, R0, N)
    
    X = 6.28*R0/0.355

    call calcScatt(RI, X,QEXT,QSCA,qabs,qbk,qpr,alb,g,n)
    K(:,1) = QEXT
    K(:,2) = Qbk
    K(:,3) = Qbk
    K(:,4) = Qbk
    K(:,5) = Qbk

    CALL DECOMPOSE_MATRIX_S(R0, K, R1, R2-1.5, 'LOG', KNTS, C, N,M, .FALSE.)

    DO I=1,N
    WRITE(*,'(6F8.3)')R0(I), K(I,:)
    ENDDO
    WRITE(*,'(A)') '–≈«”À‹“¿“ –¿«ÀŒ∆≈Õ»ﬂ œŒ ¡¿«»—Õ€Ã ‘”Õ ÷»ﬂÃ:'
    DO I=1,M
    WRITE(*,'(6F8.3)')KNTS(I),C(I,:)
    ENDDO
    end program ReconstructFDist

