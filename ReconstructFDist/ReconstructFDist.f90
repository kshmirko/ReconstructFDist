!  RECONSTRUCTFDIST.F90 
!
!  FUNCTIONS:
!  RECONSTRUCTFDIST - ENTRY POINT OF CONSOLE APPLICATION.
!

!****************************************************************************
!
!  PROGRAM: RECONSTRUCTFDIST
!
!  PURPOSE:  ENTRY POINT FOR THE CONSOLE APPLICATION.
!
!****************************************************************************

    PROGRAM RECONSTRUCTFDIST
    USE MATHUTILS, ONLY : DECOMPOSE_VECTOR_S, LINSPACE_S, LOGSPACE_S, DECOMPOSE_MATRIX_S, TRAPZ
    USE MIESCATTER, ONLY  : CALCSCATT
    USE IOMODULE
    

    IMPLICIT NONE
    INTEGER, PARAMETER  :: NRRI = 100, NIRI=30, NR1=30, NR2=20, NN=200, M=5
    REAL*8  RRI(NRRI), IRI(NIRI)
    REAL*8  R1(NR1), R2(NR2), C(M,1)
    REAL*8  R0(NN), X(NN), WL/1.064/, KNTS(M)
    REAL*8  QEXT(NN,1),QSCA(NN),QABS(NN),QBK(NN,1),QPR(NN),ALB(NN),G(NN)
    REAL*8 R01/0.05/, R02/0.5/
    REAL*8 R11/0.2/, R12/2.0/, R13/5.0/, S1, S2
    CHARACTER*80 FNAME
    COMPLEX*16 RI/(1.33,-0.005)/
    INTEGER I,J 

    CALL LINSPACE_S(R01, R02, R1, NR1)
    CALL LINSPACE_S(R11, R12, R2, NR2)
   

    DO I=1,NR1
        DO J=1, NR2
            IF((R2(J)-R1(I))<0.01) CYCLE
            CALL LINSPACE_S(R1(I), R2(J), R0, NN)
            CALL LOGSPACE_S(R1(I), R2(J), KNTS, M)
            X = 6.28*R0/WL
            CALL CALCSCATT(RI, X,QEXT(:,1),QSCA,QABS,QBK(:,1),QPR,ALB,G,NN)
            CALL DECOMPOSE_MATRIX_S(R0, QBK, R1(I), R2(J), 'LOG', KNTS, C, NN,M, .FALSE.)
            S1 = TRAPZ(R0, QBK(:,1), NN)
            S2 = TRAPZ(KNTS, C, M)
            WRITE(*,'(4F8.4, F8.2)') R1(I), R2(J), S1, S2, (S1-S2)/S1*100
        END DO
    END DO
    WRITE(*,'(A)') 'END OF PROGRAM'
!    CALL PREPAREH5FILE(FNAME, RRI, NRRI, IRI, NIRI, R1, NR1, R2, NR2)
    CALL PREPAREH5FILE()
!    INTEGER, PARAMETER :: N=100, M=5
!    INTEGER I
!    REAL*8 R0(N), X(N), Y(N,M), R1/0.05/, R2/2.0/, KNTS(M), B(N,M), YY(N,1), C(M,M), K(N,M)
!    REAL*8 QEXT(N),QSCA(N),QABS(N),QBK(N),QPR(N),ALB(N),G(N)
!    COMPLEX*16 RI/(1.33,-0.005)/
!       
!    KNTS = 0
!    R0=0
!    CALL LINSPACE_S(R1, R2, R0, N)
!    
!    X = 6.28*R0/0.355
!
!    CALL CALCSCATT(RI, X,QEXT,QSCA,QABS,QBK,QPR,ALB,G,N)
!    K(:,1) = QEXT
!    K(:,2) = QBK
!    K(:,3) = QBK
!    K(:,4) = QBK
!    K(:,5) = QBK
!
!    CALL DECOMPOSE_MATRIX_S(R0, K, R1, R2-1.5, 'LOG', KNTS, C, N,M, .FALSE.)
!
!    DO I=1,N
!    WRITE(*,'(6F8.3)')R0(I), K(I,:)
!    ENDDO
!    WRITE(*,'(A)') 'ÐÅÇÓËÜÒÀÒ ÐÀÇËÎÆÅÍÈß ÏÎ ÁÀÇÈÑÍÛÌ ÔÓÍÊÖÈßÌ:'
!    DO I=1,M
!    WRITE(*,'(6F8.3)')KNTS(I),C(I,:)
!    ENDDO
    END PROGRAM RECONSTRUCTFDIST

