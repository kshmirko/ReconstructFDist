module MathUtils
USE mkl95_PRECISION, ONLY: WP => DP
USE mkl95_LAPACK, ONLY: GELSS
implicit none

contains

subroutine logspace_s(r1, r2, R, N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!������� ��������������-�������������� 
!������������������
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8, intent(in)  ::  r1, r2
integer*4, intent(in)   :: N
integer I
real*8, intent(inout) :: R(N)
real*8  lr1, lr2, dr
! ������� ��������� ������������ ������������������
R = (/(i,i=0,N-1)/)
! �������� ����� � ������ ������� ��������� � ����������
lr1 = log10(r1)
lr2 = log10(r2)
! ������������ ��� ����������������� � ����������
dr = (lr2-lr1)/(N-1)

! ����������� �� ����������, ���������� ���������� ������
R = 10**(R(1:N)*dr+lr1)

end subroutine logspace_s

subroutine linspace_s(r1, r2, R, N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!������� �������������� ������������������
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8, intent(in)  ::  r1, r2
integer*4, intent(in)   :: N
integer I
real*8 :: R(N), dr

! ������� ��������� ������������ ������������������ � ����� 1 ������� � ����
R = (/(i,i=0,N-1)/)

! ������������ ��� ����������������� 
dr = (r2-r1)/(N-1)

!��������� ������������������ � �������� ��������
R = R(1:N)*dr+r1
end subroutine linspace_s


function logspace_f(r1, r2, N) result(R)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!������� ��������������-�������������� 
!������������������
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8, intent(in)  ::  r1, r2
integer*4, intent(in)   :: N
integer I
real*8 :: R(N), lr1, lr2, dr
! ������� ��������� ������������ ������������������
R = (/(i,i=0,N-1)/)
! �������� ����� � ������ ������� ��������� � ����������
lr1 = log10(r1)
lr2 = log10(r2)
! ������������ ��� ����������������� � ����������
dr = (lr2-lr1)/(N-1)

! ����������� �� ����������, ���������� ���������� ������
R = 10**(R*dr+lr1)

end function logspace_f

function linspace_f(r1, r2, N) result(R)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!������� �������������� ������������������
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8, intent(in)  ::  r1, r2
integer*4, intent(in)   :: N
integer I
real*8 :: R(N), dr

! ������� ��������� ������������ ������������������ � ����� 1 ������� � ����
R = (/(i,i=0,N-1)/)

! ������������ ��� ����������������� 
dr = (r2-r1)/(N-1)

!��������� ������������������ � �������� ��������
R = R*dr+r1
end function linspace_f


subroutine tri(Xc, t, B, N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ������� 1 ����������� �������� �������
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8, intent(in)  :: Xc
integer, intent(in) :: N
real*8, intent(in)  :: t(N)
real*8, intent(out) :: B(N)
integer I

B = max(0.0,1.0-abs(Xc-t))

end subroutine tri

subroutine basisMatrix(B,R,knts,N,M)
integer, intent(IN) :: N,M
real*8, intent(IN)  :: R(N), knts(M)
real*8, intent(OUT) :: B(N,M)
integer I

! ������ ������� ������� - ����������� �������� ������ � ������� � ����� knts
do I=1,M
    call tri(knts(I), R, B(:,I), N)
enddo

where (R .LT. knts(1))
    B(:,1)=0.0
endwhere

where (R .GT. knts(M))
    B(:,M)=0.0
endwhere
end subroutine basisMatrix


subroutine calcCoefs(B,F,C,N,M)
integer, intent(IN) ::  N,M
real*8, intent(IN)  ::  B(N,M), F(N,1)
real*8, intent(OUT) ::  C(M)
INTEGER RANK, INFO
REAL*8 S(M), BB(N,M), FF(N,1)
BB(:,:) = B(:,:)
FF(:,:) = F(:,:)
CALL GELSS( BB, FF, RANK, S, 0.00000001_WP, INFO=INFO )
C = FF(1:M,1)
end subroutine calcCoefs

endmodule MathUtils