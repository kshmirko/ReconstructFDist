MODULE MATHUTILS
USE MKL95_PRECISION, ONLY: WP => DP
USE MKL95_LAPACK, ONLY: GELSS
IMPLICIT NONE

CONTAINS

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!янгдюел кнцюпхтлхвеяйх-щйбхдхярюмрмсч 
!онякеднбюрекэмнярэ
! R=10**[LR1,LR2,...,LRN], DLR = (LRN-LR1)/(N-1)
! LR1 = LOG10(R1)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE LOGSPACE_S(R1, R2, R, N)
REAL*8, INTENT(IN)  ::  R1, R2
INTEGER*4, INTENT(IN)   :: N
INTEGER I
REAL*8, INTENT(INOUT) :: R(N)
REAL*8  LR1, LR2, DR
! янгдюел лнмнрнммн бнгпюярючысч онякеднбюрекэмнярэ
R = (/(I,I=0,N-1)/)
! онксвюел кеясч х опюбсч цпюмхжш хмрепбюкю б кнцюпхтлюу
LR1 = LOG10(R1)
LR2 = LOG10(R2)
! пюяявхршбюел ьюц онякеднбюрекэмнях б кнцюпхтлюу
DR = (LR2-LR1)/(N-1)

! хгаюбкъеляъ нр кнцюпхтлнб, онремжхпсъ онксвеммши бейрнп
R = 10**(R(1:N)*DR+LR1)
END SUBROUTINE LOGSPACE_S

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!янгдюел щйбхдхярюмрмсч онякеднбюрекэмнярэ
!R = [R1,R2,...RN], DR = (RN-R1)/(N-1)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE LINSPACE_S(R1, R2, R, N)
REAL*8, INTENT(IN)  ::  R1, R2
INTEGER*4, INTENT(IN)   :: N
INTEGER I
REAL*8 :: R(N), DR

! янгдюел лнмнрнммн бнгпюярючысч онякеднбюрекэмнярэ я ьюцнл 1 мювхмюъ я мскъ
R = (/(I,I=0,N-1)/)

! пюяявхршбюел ьюц онякеднбюрекэмнях 
DR = (R2-R1)/(N-1)

!бнплхпсел онякеднбюрекэмнярэ б гюдюммшу цпюмхжюу
R = R(1:N)*DR+R1
END SUBROUTINE LINSPACE_S


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!янгдюел кнцюпхтлхвеяйх-щйбхдхярюмрмсч 
!онякеднбюрекэмнярэ
! R=10**[LR1,LR2,...,LRN], DLR = (LRN-LR1)/(N-1)
! LR1 = LOG10(R1)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
FUNCTION LOGSPACE_F(R1, R2, N) RESULT(R)
REAL*8, INTENT(IN)  ::  R1, R2
INTEGER*4, INTENT(IN)   :: N
INTEGER I
REAL*8 :: R(N), LR1, LR2, DR
! янгдюел лнмнрнммн бнгпюярючысч онякеднбюрекэмнярэ
R = (/(I,I=0,N-1)/)
! онксвюел кеясч х опюбсч цпюмхжш хмрепбюкю б кнцюпхтлюу
LR1 = LOG10(R1)
LR2 = LOG10(R2)
! пюяявхршбюел ьюц онякеднбюрекэмнях б кнцюпхтлюу
DR = (LR2-LR1)/(N-1)

! хгаюбкъеляъ нр кнцюпхтлнб, онремжхпсъ онксвеммши бейрнп
R = 10**(R*DR+LR1)
END FUNCTION LOGSPACE_F


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!янгдюел щйбхдхярюмрмсч онякеднбюрекэмнярэ
!R = [R1,R2,...RN], DR = (RN-R1)/(N-1)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
FUNCTION LINSPACE_F(R1, R2, N) RESULT(R)
REAL*8, INTENT(IN)  ::  R1, R2
INTEGER*4, INTENT(IN)   :: N
INTEGER I
REAL*8 :: R(N), DR

! янгдюел лнмнрнммн бнгпюярючысч онякеднбюрекэмнярэ я ьюцнл 1 мювхмюъ я мскъ
R = (/(I,I=0,N-1)/)

! пюяявхршбюел ьюц онякеднбюрекэмнях 
DR = (R2-R1)/(N-1)

!бнплхпсел онякеднбюрекэмнярэ б гюдюммшу цпюмхжюу
R = R*DR+R1
END FUNCTION LINSPACE_F

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! янгдюел 1 рпесцнкэмсч аюгхямсч тсмйжхч
! TRI(XC,T)=MAX(0, 1-|XC-T|)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE TRI(XC, T, B, N)
REAL*8, INTENT(IN)  :: XC
INTEGER, INTENT(IN) :: N
REAL*8, INTENT(IN)  :: T(N)
REAL*8, INTENT(OUT) :: B(N)
INTEGER I

B = MAX(0.0,1.0-ABS(XC-T))

END SUBROUTINE TRI

SUBROUTINE BASISMATRIX(B,R,KNTS,N,M)
INTEGER, INTENT(IN) :: N,M
REAL*8, INTENT(IN)  :: R(N), KNTS(M)
REAL*8, INTENT(OUT) :: B(N,M)
INTEGER I

! йюфдши ярнкаеж люрпхжш - рпесцнкэмши аюгхямши бейрнп я жемрпнл б сгкюу KNTS
DO I=1,M
    CALL TRI(KNTS(I), R, B(:,I), N)
ENDDO

!с оепбнцн х онякедмецн аюгхямшу бейрнпнб сахпюел 
!кебсч х опюбсч вюярх яннрберярбеммн
WHERE (R .LT. KNTS(1))
    B(:,1)=0.0
ENDWHERE

WHERE (R .GT. KNTS(M))
    B(:,M)=0.0
ENDWHERE
END SUBROUTINE BASISMATRIX

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C бшвхякъер йнщттхжхемрш пюгкнфемхъ 
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
!C якедсчыюъ опнжедспю нясыеярбкъер пюгкнфемхе он рпесцнкэмнлс аюгхяс 
!C рюакхвмн гюдюммжч тсмйжхч (R,F) б сгкнбшу рнвйюу {KNTS}.
!C пегскэрюр пюгкнфемхъ бнгбпюыюеряъ б бхде йнщттхжхемрнб {C}
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE DECOMPOSE_VECTOR_S(R, F, R1, R2, FLAG, KNTS, C, N,M)
INTEGER*4, INTENT(IN)   ::  N,M
REAL*8, INTENT(IN)      ::  R(N), F(N,1), R1, R2
REAL*8, INTENT(OUT)     ::  C(M,1), KNTS(M)
CHARACTER*3, INTENT(IN) ::  FLAG
REAL*8  :: TR1, TR2, TDR, RR(N), B(N,M), LKNTS(M), S(M), FF(N,1)
INTEGER I, RANK, INFO

!тнплхпсел ьюц, ю рюйфе кебсч х опюбсч цпюмхжш
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
!C якедсчыюъ опнжедспю нясыеярбкъер пюгкнфемхе он рпесцнкэмнлс аюгхяс 
!C рюакхвмн гюдюммжч тсмйжхч (R,FI), цде I=1..M - вхякн тсмйжхи 
!C б сгкнбшу рнвйюу {KNTS}.
!C пегскэрюр пюгкнфемхъ бнгбпюыюеряъ б бхде йнщттхжхемрнб {C}
!C R1, R2 - цпюмхжш хмрепбюкю пюгкнфемхъ
!C FLAG - яоняна пюгкнфемхъ 'LOG' - кнц-щйбхдхярюмрхмши, 'LIN' - щйбхдхярюмрхмши
!C N,M - пюглепмнярх бундмшу оюпюлерпнб
!C R(N), F(N,M) - врн пюяйкюдшбюел
!C KNTS(M) - бейрнп я нряверюлх [R1, ... RM]
!C TRANS - .TRUE. - бшонкмъел рпюмяонмхпнбюмхе, .FALSE. - мер
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE DECOMPOSE_MATRIX_S(R, F, R1, R2, FLAG, KNTS, C, N,M, TRANS)
INTEGER*4, INTENT(IN)   ::  N,M
REAL*8, INTENT(IN)      ::  R(N), F(N,M), R1, R2
REAL*8, INTENT(OUT)     ::  C(M,M), KNTS(M)
CHARACTER*3, INTENT(IN) ::  FLAG
LOGICAL TRANS !IF .TRUE. THAN PERFORM TRANSPOSE OPERATOR.
REAL*8  :: TR1, TR2, TDR, RR(N), B(N,M), LKNTS(M), S(M), FF(N,M)
INTEGER I, RANK, INFO

!тнплхпсел ьюц, ю рюйфе кебсч х опюбсч цпюмхжш
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

FUNCTION TRAPZ(X,Y,N) RESULT(RES)
    INTEGER, INTENT(IN) ::  N
    INTEGER I
    REAL*8, INTENT(IN)  ::  X(N), Y(N)
    REAL*8  RES
    RES=0

    DO I=1, N-1
        RES=RES+0.5*(Y(I)+Y(I+1))*(X(I+1)-X(I))
    END DO
    return
END FUNCTION TRAPZ

END MODULE MATHUTILS