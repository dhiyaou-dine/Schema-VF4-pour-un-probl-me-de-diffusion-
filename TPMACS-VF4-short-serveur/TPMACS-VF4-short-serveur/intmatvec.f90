module intmatvec
  use longr 
  USE imprime
  USE parmmage
  implicit none

  interface operator (*)
     module procedure matvec
  end interface
contains
  FUNCTION matvec(A, U)
    !     *--------------------------------------
    !     * Ce sous programme calcule le produit d'une matrice creuse par un vecteur
    ! Centrale Nantes, M. SAAD

    IMPLICIT NONE

    !--------------------------
    ! Declaration des arguments
    !--------------------------
    TYPE(MatCreux), intent(in)                        :: A
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 ),intent(in) :: U
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 )            :: matvec
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)      :: oldprf
    INTEGER               :: i, j, ndim

    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'MATVEC'

    !------
    ! Corps
    !------
    matvec = 0.D0
    ndim = SIZE(A%IndPL) -1 
    DO i = 1, ndim
       DO j = A%IndPL(i) ,  A%IndPL(i+1) -1 
          matvec(i) =  matvec(i) + A%TMat( j )*U(A%Indc(j))
       END DO
    ENDDO

    !------------
    ! Impressions
    !------------
!!$  IF (iprint >= 2) THEN
!!$     WRITE(uprint,*)' SOLUTION DE MATVEC'
!!$     WRITE(uprint,100) (matvec(i), i=1,ndim)
!!$  END IF
100 FORMAT (10(E10.3, 2X))
    !-----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf

    RETURN
  END FUNCTION matvec


end module intmatvec
