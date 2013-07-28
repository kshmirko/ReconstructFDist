MODULE MathUtils
USE mkl95_PRECISION, ONLY: WP => DP
USE mkl95_LAPACK, ONLY: GELSS
IMPLICIT NONE

CONTAINS

!ccccccccccccccccccccccccccccccccccccccccccccc
!������� ��������������-�������������� 
!������������������
! R=10**[lR1,lR2,...,lRn], dlR = (lRn-lR1)/(N-1)
! lR1 = LOG10(R1)
!ccccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE LOGSPACE_S(r1, r2, R, N)
REAL*8, INTENT(IN)  ::  r1, r2
INTEGER*4, INTENT(IN)   :: N
INTEGER I
REAL*8, INTENT(inout) :: R(N)
REAL*8  lr1, lr2, dr
! ������� ��������� ������������ ������������������
R = (/(i,i=0,N-1)/)
! �������� ����� � ������ ������� ��������� � ����������
lr1 = LOG10(r1)
lr2 = LOG10(r2)
! ������������ ��� ����������������� � ����������
dr = (lr2-lr1)/(N-1)

! ����������� �� ����������, ���������� ���������� ������
R = 10**(R(1:N)*dr+lr1)
END SUBROUTINE LOGSPACE_S

!cccccccccccccccccccccccccccccccccccccccccccc
!������� �������������� ������������������
!R = [R1,R2,...Rn], dR = (Rn-R1)/(N-1)
!cccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE LINSPACE_S(r1, r2, R, N)
REAL*8, INTENT(IN)  ::  r1, r2
INTEGER*4, INTENT(IN)   :: N
INTEGER I
REAL*8 :: R(N), dr

! ������� ��������� ������������ ������������������ � ����� 1 ������� � ����
R = (/(i,i=0,N-1)/)

! ������������ ��� ����������������� 
dr = (r2-r1)/(N-1)

!��������� ������������������ � �������� ��������
R = R(1:N)*dr+r1
END SUBROUTINE linspace_s


!ccccccccccccccccccccccccccccccccccccccccccccc
!������� ��������������-�������������� 
!������������������
! R=10**[lR1,lR2,...,lRn], dlR = (lRn-lR1)/(N-1)
! lR1 = LOG10(R1)
!ccccccccccccccccccccccccccccccccccccccccccccc
FUNCTION LOGSPACE_F(r1, r2, N) RESULT(R)
REAL*8, INTENT(IN)  ::  r1, r2
INTEGER*4, INTENT(IN)   :: N
INTEGER I
REAL*8 :: R(N), lr1, lr2, dr
! ������� ��������� ������������ ������������������
R = (/(i,i=0,N-1)/)
! �������� ����� � ������ ������� ��������� � ����������
lr1 = LOG10(r1)
lr2 = LOG10(r2)
! ������������ ��� ����������������� � ����������
dr = (lr2-lr1)/(N-1)

! ����������� �� ����������, ���������� ���������� ������
R = 10**(R*dr+lr1)
END FUNCTION LOGSPACE_F


!cccccccccccccccccccccccccccccccccccccccccccc
!������� �������������� ������������������
!R = [R1,R2,...Rn], dR = (Rn-R1)/(N-1)
!cccccccccccccccccccccccccccccccccccccccccccc
FUNCTION LINSPACE_F(r1, r2, N) RESULT(R)
REAL*8, INTENT(IN)  ::  r1, r2
INTEGER*4, INTENT(IN)   :: N
INTEGER I
REAL*8 :: R(N), dr

! ������� ��������� ������������ ������������������ � ����� 1 ������� � ����
R = (/(i,i=0,N-1)/)

! ������������ ��� ����������������� 
dr = (r2-r1)/(N-1)

!��������� ������������������ � �������� ��������
R = R*dr+r1
END FUNCTION LINSPACE_F

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ������� 1 ����������� �������� �������
! tri(Xc,t)=max(0, 1-|Xc-t|)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE TRI(Xc, t, B, N)
REAL*8, INTENT(IN)  :: Xc
INTEGER, INTENT(IN) :: N
REAL*8, INTENT(IN)  :: t(N)
REAL*8, INTENT(OUT) :: B(N)
INTEGER I

B = MAX(0.0,1.0-ABS(Xc-t))

END SUBROUTINE TRI

SUBROUTINE BASISMATRIX(B,R,knts,N,M)
INTEGER, INTENT(IN) :: N,M
REAL*8, INTENT(IN)  :: R(N), knts(M)
REAL*8, INTENT(OUT) :: B(N,M)
INTEGER I

! ������ ������� ������� - ����������� �������� ������ � ������� � ����� knts
do I=1,M
    CALL TRI(knts(I), R, B(:,I), N)
enddo

!� ������� � ���������� �������� �������� ������� 
!����� � ������ ����� ��������������
where (R .LT. knts(1))
    B(:,1)=0.0
endwhere

where (R .GT. knts(M))
    B(:,M)=0.0
endwhere
END SUBROUTINE BASISMATRIX

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C ��������� ������������ ���������� 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE CALCCOEFS(B,F,C,N,M)
INTEGER, INTENT(IN) ::  N,M
REAL*8, INTENT(IN)  ::  B(N,M), F(N,1)
REAL*8, INTENT(OUT) ::  C(M)
INTEGER RANK, INFO
REAL*8 S(M), BB(N,M), FF(N,1)
BB(:,:) = B(:,:)
FF(:,:) = F(:,:)
CALL GELSS( BB, FF, RANK, S, 0.00000001_WP, INFO=INFO )
C = FF(1:M,1)
END SUBROUTINE CALCCOEFS


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C ��������� ��������� ������������ ���������� �� ������������ ������ 
!C �������� �������� ������� (R,F) � ������� ������ {KNTS}.
!C ��������� ���������� ������������ � ���� ������������� {C}
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE DECOMPOSE_VECTOR_S(R, F, R1, R2, FLAG, KNTS, C, N,M)
INTEGER*4, INTENT(IN)   ::  N,M
REAL*8, INTENT(IN)      ::  R(N), F(N,1), R1, R2
REAL*8, INTENT(OUT)     ::  C(M,1), KNTS(M)
CHARACTER*3, INTENT(IN) ::  FLAG
REAL*8  :: TR1, TR2, TDR, RR(N), B(N,M), LKNTS(M), S(M), FF(N,1)
INTEGER I, RANK, INFO

!��������� ���, � ����� ����� � ������ �������
IF (FLAG .EQ. 'LOG') THEN
    CALL LOGSPACE_S(R1, R2, KNTS, M)
    LKNTS = LOG10(KNTS)
    TDR = LKNTS(2)-LKNTS(1)
    RR = LOG10(R)
    TR1=LKNTS(1)
ELSE
    CALL LINSPACE_S(R1, R2, KNTS, M)
    LKNTS=KNTS
    TDR = KNTS(2)-KNTS(1)
    TR1 = LKNTS(1)
    RR=R
ENDIF
LKNTS=(LKNTS-TR1)/TDR
RR=(RR-TR1)/TDR

CALL BASISMATRIX(B,RR,LKNTS,N,M)
FF=F
CALL GELSS( B, FF, RANK, S, 0.00000001_WP, INFO=INFO )
C(1:M,1) = FF(1:M,1)

END SUBROUTINE DECOMPOSE_VECTOR_S

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C ��������� ��������� ������������ ���������� �� ������������ ������ 
!C �������� �������� ������� (R,Fi), ��� I=1..M - ����� ������� 
!C � ������� ������ {KNTS}.
!C ��������� ���������� ������������ � ���� ������������� {C}
!C R1, R2 - ������� ��������� ����������
!C FLAG - ������ ���������� 'LOG' - ���-���������������, 'LIN' - ���������������
!C N,M - ����������� ������� ����������
!C R(N), F(N,M) - ��� ������������
!C KNTS(M) - ������ � ��������� [R1, ... RM]
!C TRANS - .TRUE. - ��������� ����������������, .FALSE. - ���
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE DECOMPOSE_MATRIX_S(R, F, R1, R2, FLAG, KNTS, C, N,M, TRANS)
INTEGER*4, INTENT(IN)   ::  N,M
REAL*8, INTENT(IN)      ::  R(N), F(N,M), R1, R2
REAL*8, INTENT(OUT)     ::  C(M,M), KNTS(M)
CHARACTER*3, INTENT(IN) ::  FLAG
LOGICAL TRANS !IF .TRUE. THAN PERFORM TRANSPOSE OPERATOR.
REAL*8  :: TR1, TR2, TDR, RR(N), B(N,M), LKNTS(M), S(M), FF(N,M)
INTEGER I, RANK, INFO

!��������� ���, � ����� ����� � ������ �������
IF (FLAG .EQ. 'LOG') THEN
    CALL LOGSPACE_S(R1, R2, KNTS, M)
    LKNTS = LOG10(KNTS)
    TDR = LKNTS(2)-LKNTS(1)
    RR = LOG10(R)
    TR1=LKNTS(1)
ELSE
    CALL LINSPACE_S(R1, R2, KNTS, M)
    LKNTS=KNTS
    TDR = KNTS(2)-KNTS(1)
    TR1 = LKNTS(1)
    RR=R
ENDIF
LKNTS=(LKNTS-TR1)/TDR
RR=(RR-TR1)/TDR

CALL BASISMATRIX(B,RR,LKNTS,N,M)
FF=F
CALL GELSS( B, FF, RANK, S, 0.00000001_WP, INFO=INFO )

C(1:M,1:M) = FF(1:M,1:M)
IF(TRANS .EQ. .TRUE.) C = TRANSPOSE(C)


END SUBROUTINE DECOMPOSE_MATRIX_S



END MODULE MathUtils