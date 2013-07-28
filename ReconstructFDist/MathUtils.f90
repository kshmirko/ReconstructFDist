MODULE MathUtils
USE mkl95_PRECISION, ONLY: WP => DP
USE mkl95_LAPACK, ONLY: GELSS
IMPLICIT NONE

CONTAINS

!ccccccccccccccccccccccccccccccccccccccccccccc
!создаем логарифмически-эквидистантную 
!последовательность
! R=10**[lR1,lR2,...,lRn], dlR = (lRn-lR1)/(N-1)
! lR1 = LOG10(R1)
!ccccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE LOGSPACE_S(r1, r2, R, N)
REAL*8, INTENT(IN)  ::  r1, r2
INTEGER*4, INTENT(IN)   :: N
INTEGER I
REAL*8, INTENT(inout) :: R(N)
REAL*8  lr1, lr2, dr
! Создаем монотонно возрастающую последовательность
R = (/(i,i=0,N-1)/)
! получаем лесую и правую границы интервала в логарифмах
lr1 = LOG10(r1)
lr2 = LOG10(r2)
! рассчитываем шаг последовательноси в логарифмах
dr = (lr2-lr1)/(N-1)

! избавляемся от логарифмов, потенцируя полученный вектор
R = 10**(R(1:N)*dr+lr1)
END SUBROUTINE LOGSPACE_S

!cccccccccccccccccccccccccccccccccccccccccccc
!создаем эквидистантную последовательность
!R = [R1,R2,...Rn], dR = (Rn-R1)/(N-1)
!cccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE LINSPACE_S(r1, r2, R, N)
REAL*8, INTENT(IN)  ::  r1, r2
INTEGER*4, INTENT(IN)   :: N
INTEGER I
REAL*8 :: R(N), dr

! Создаем монотонно возрастающую последовательность с шагом 1 начиная с нуля
R = (/(i,i=0,N-1)/)

! рассчитываем шаг последовательноси 
dr = (r2-r1)/(N-1)

!Вормируем последовательность в заданных границах
R = R(1:N)*dr+r1
END SUBROUTINE linspace_s


!ccccccccccccccccccccccccccccccccccccccccccccc
!создаем логарифмически-эквидистантную 
!последовательность
! R=10**[lR1,lR2,...,lRn], dlR = (lRn-lR1)/(N-1)
! lR1 = LOG10(R1)
!ccccccccccccccccccccccccccccccccccccccccccccc
FUNCTION LOGSPACE_F(r1, r2, N) RESULT(R)
REAL*8, INTENT(IN)  ::  r1, r2
INTEGER*4, INTENT(IN)   :: N
INTEGER I
REAL*8 :: R(N), lr1, lr2, dr
! Создаем монотонно возрастающую последовательность
R = (/(i,i=0,N-1)/)
! получаем лесую и правую границы интервала в логарифмах
lr1 = LOG10(r1)
lr2 = LOG10(r2)
! рассчитываем шаг последовательноси в логарифмах
dr = (lr2-lr1)/(N-1)

! избавляемся от логарифмов, потенцируя полученный вектор
R = 10**(R*dr+lr1)
END FUNCTION LOGSPACE_F


!cccccccccccccccccccccccccccccccccccccccccccc
!создаем эквидистантную последовательность
!R = [R1,R2,...Rn], dR = (Rn-R1)/(N-1)
!cccccccccccccccccccccccccccccccccccccccccccc
FUNCTION LINSPACE_F(r1, r2, N) RESULT(R)
REAL*8, INTENT(IN)  ::  r1, r2
INTEGER*4, INTENT(IN)   :: N
INTEGER I
REAL*8 :: R(N), dr

! Создаем монотонно возрастающую последовательность с шагом 1 начиная с нуля
R = (/(i,i=0,N-1)/)

! рассчитываем шаг последовательноси 
dr = (r2-r1)/(N-1)

!Вормируем последовательность в заданных границах
R = R*dr+r1
END FUNCTION LINSPACE_F

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! создаем 1 треугольную базисную функцию
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

! Каждый столбец матрицы - треугольный базисный вектор с центром в узлах knts
do I=1,M
    CALL TRI(knts(I), R, B(:,I), N)
enddo

!У первого и последнего базисных векторов убираем 
!левую и правую части соответственно
where (R .LT. knts(1))
    B(:,1)=0.0
endwhere

where (R .GT. knts(M))
    B(:,M)=0.0
endwhere
END SUBROUTINE BASISMATRIX

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C ВЫЧИСЛЯЕТ КОЭФФИЦИЕНТЫ РАЗЛОЖЕНИЯ 
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
!C СЛЕДУЮЩАЯ ПРОЦЕДУРА ОСУЩЕСТВЛЯЕТ РАЗЛОЖЕНИЕ ПО ТРЕУГОЛЬНОМУ БАЗИСУ 
!C ТАБЛИЧНО ЗАДАННЦЮ ФУНКЦИЮ (R,F) В УЗЛОВЫХ ТОЧКАХ {KNTS}.
!C РЕЗУЛЬТАТ РАЗЛОЖЕНИЯ ВОЗВРАЩАЕТСЯ В ВИДЕ КОЭФФИЦИЕНТОВ {C}
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE DECOMPOSE_VECTOR_S(R, F, R1, R2, FLAG, KNTS, C, N,M)
INTEGER*4, INTENT(IN)   ::  N,M
REAL*8, INTENT(IN)      ::  R(N), F(N,1), R1, R2
REAL*8, INTENT(OUT)     ::  C(M,1), KNTS(M)
CHARACTER*3, INTENT(IN) ::  FLAG
REAL*8  :: TR1, TR2, TDR, RR(N), B(N,M), LKNTS(M), S(M), FF(N,1)
INTEGER I, RANK, INFO

!ФОРМИРУЕМ ШАГ, А ТАКЖЕ ЛЕВУЮ И ПРАВУЮ ГРАНИЦЫ
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
!C СЛЕДУЮЩАЯ ПРОЦЕДУРА ОСУЩЕСТВЛЯЕТ РАЗЛОЖЕНИЕ ПО ТРЕУГОЛЬНОМУ БАЗИСУ 
!C ТАБЛИЧНО ЗАДАННЦЮ ФУНКЦИЮ (R,Fi), ГДЕ I=1..M - ЧИСЛО ФУНКЦИЙ 
!C В УЗЛОВЫХ ТОЧКАХ {KNTS}.
!C РЕЗУЛЬТАТ РАЗЛОЖЕНИЯ ВОЗВРАЩАЕТСЯ В ВИДЕ КОЭФФИЦИЕНТОВ {C}
!C R1, R2 - ГРАНИЦЫ ИНТЕРВАЛА РАЗЛОЖЕНИЯ
!C FLAG - СПОСОБ РАЗЛОЖЕНИЯ 'LOG' - ЛОГ-ЭКВИДИСТАНТИНЫЙ, 'LIN' - ЭКВИДИСТАНТИНЫЙ
!C N,M - РАЗМЕРНОСТИ ВХОДНЫХ ПАРАМЕТРОВ
!C R(N), F(N,M) - ЧТО РАСКЛАДЫВАЕМ
!C KNTS(M) - ВЕКТОР С ОТСЧЕТАМИ [R1, ... RM]
!C TRANS - .TRUE. - ВЫПОЛНЯЕМ ТРАНСПОНИРОВАНИЕ, .FALSE. - НЕТ
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE DECOMPOSE_MATRIX_S(R, F, R1, R2, FLAG, KNTS, C, N,M, TRANS)
INTEGER*4, INTENT(IN)   ::  N,M
REAL*8, INTENT(IN)      ::  R(N), F(N,M), R1, R2
REAL*8, INTENT(OUT)     ::  C(M,M), KNTS(M)
CHARACTER*3, INTENT(IN) ::  FLAG
LOGICAL TRANS !IF .TRUE. THAN PERFORM TRANSPOSE OPERATOR.
REAL*8  :: TR1, TR2, TDR, RR(N), B(N,M), LKNTS(M), S(M), FF(N,M)
INTEGER I, RANK, INFO

!ФОРМИРУЕМ ШАГ, А ТАКЖЕ ЛЕВУЮ И ПРАВУЮ ГРАНИЦЫ
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