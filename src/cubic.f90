! *****************************************************************************
! *                                                                           *
! * Purpose: Find the roots of a cubic equation by first applying Newton's    *
! *          Method to find a real root and then using synthetic division to  *
! *          find a quadratic equation whose roots can be found by using the  *
! *          quadratic formula.                                               *
! *****************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Program Execution !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PROGRAM cubicSolver
  IMPLICIT NONE
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Type Declarations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  TYPE QUADRATIC
    REAL :: A, B, C
    COMPLEX :: R1, R2
  END TYPE QUADRATIC
  
  TYPE CUBIC
    REAL :: A, B, C, D
    REAL :: R1
    TYPE(QUADRATIC), POINTER :: quadraticEq
  END TYPE CUBIC
  
  TYPE(CUBIC), POINTER :: cubicEq
  
  NULLIFY(cubicEq)
  ALLOCATE(cubicEq)
  NULLIFY(cubicEq%quadraticEq)
  
  WRITE(*,*) "Enter the coefficients for a cubic equation to be solved."
  WRITE(*, '(A)', advance = "no") "A: "
  READ(*,*) cubicEq%A
  WRITE(*, '(A)', advance = "no") "B: "
  READ(*,*) cubicEq%B
  WRITE(*, '(A)', advance = "no") "C: "
  READ(*,*) cubicEq%C
  WRITE(*, '(A)', advance = "no") "D: "
  READ(*,*) cubicEq%D
  
  WRITE(*,*) "Cubic equation to be solved:", cubicEq%A, "x^3 +", cubicEq%B, &
    "x^2 +", cubicEq%C, "x +", cubicEq%D
  
  ! Seed differently maybe?
  cubicEq%R1 = newtonsRoot(cubicEq, 0.)
  
  WRITE(*,*) "R1 =", cubicEq%R1
  
  ALLOCATE(cubicEq%quadraticEq)
  cubicEq%quadraticEq = syntheticDivision(cubicEq)
  
  WRITE(*,*) "Cubic reduces to: "
  WRITE(*,*) cubicEq%QuadraticEq%A, "x^2 +", cubicEq%QuadraticEq%B, "x +", &
                                                          cubicEq%QuadraticEq%C
  
  cubicEq%quadraticEq = quadraticRoots(cubicEq)
  
  WRITE(*,*) "R2 =", cubicEq%quadraticEq%R1
  WRITE(*,*) "R3 =", cubicEq%quadraticEq%R2
  
  DEALLOCATE(cubicEq%quadraticEq)
  DEALLOCATE(cubicEq)

  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! newtonsRoot                                         
  !     Uses Newtons Method to find a single            
  !     real root of a cubic equation                   
  !                                                     
  ! Param:  cubicEq - Cubic equation to be solved       
  !               X - X value to begin finding the root 
  !                   This should be less than the first
  !                   max/min if leading coefficient is 
  !                   positive, or after the second     
  !                   max/min if the leading coefficient
  !                   is negative.                      
  ! Return: newtonsRoot - A single real root            
  !            What if the approximation fails?         
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION newtonsRoot(cubicEq, Xin)
  TYPE(CUBIC), POINTER :: cubicEq
  REAL :: newtonsRoot
  REAL :: M, Xin, X, Y, prevX
  INTEGER :: matching
  
  matching = 0
  X = Xin
  DO
    prevX = X
    Y = cubicEq%A * X**3 + cubicEq%B * X**2 + cubicEq%C * X + cubicEq%D
    M = 3 * cubicEq%A * X**2 + 2 * cubicEq%B * X + cubicEq%C
    IF (M == 0) THEN
      IF (X == 0) X = X + 1
      IF (cubicEq%A < 0) THEN
        X = ABS(X) * 2
      ELSE
        X = -ABS(X) * 2
      ENDIF
      WRITE(*,*) "X =", X
    ENDIF
    X = prevX - Y / M
    IF(X == prevX) matching = matching + 1
    IF(matching > 2) EXIT
  END DO
  
  newtonsRoot = X
  
  END FUNCTION newtonsRoot
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! syntheticDivision
  !     Does something
  !
  ! Param:  cubicEq -
  !            root - 
  ! Return: syntheticDivision - 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION syntheticDivision(cubicEq)
  TYPE(CUBIC), POINTER, intent(in) :: cubicEq
  TYPE(QUADRATIC) :: syntheticDivision
  
  syntheticDivision%A = cubicEq%A
  syntheticDivision%B = syntheticDivision%A * cubicEq%R1 + cubicEq%B
  syntheticDivision%C = syntheticDivision%B * cubicEq%R1 + cubicEq%C
  
  END FUNCTION syntheticDivision
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Sub-Routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! quadraticRoots
  !     
  ! Param:  cubicEq - cubicEq containing quadratic
  !                   to be solved
  ! Return: quadraticRoots - quadratic with roots solved
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION quadraticRoots(cubicEq)
  TYPE(CUBIC), POINTER, intent(in) :: cubicEq
  TYPE(QUADRATIC) :: quadraticRoots
  
  quadraticRoots = cubicEq%quadraticEq
  
  quadraticRoots%R1 = (-(quadraticRoots%B) - sqrt(CMPLX((quadraticRoots%B)**2 &
           - 4* quadraticRoots%A * quadraticRoots%C))) / (2 * quadraticRoots%A)
  quadraticRoots%R2 = (-(quadraticRoots%B) + sqrt(CMPLX((quadraticRoots%B)**2 &
           - 4* quadraticRoots%A * quadraticRoots%C))) / (2 * quadraticRoots%A)
  
  END FUNCTION quadraticRoots
  
  SUBROUTINE printQuadratic(cubicEq)
    TYPE(CUBIC), POINTER, intent(in) :: cubicEq
    
    write(*,*) "Cubic reduces to: "
    WRITE(*,*) cubicEq%QuadraticEq%A, "x^2 +", cubicEq%QuadraticEq%B, "x +", cubicEq%QuadraticEq%C
  END SUBROUTINE printQuadratic
  
  END PROGRAM cubicSolver
