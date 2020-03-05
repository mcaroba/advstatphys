program LJ_B2

! This program calculates the B2(T) coefficient according to some
! parametrization of the Lennard-Jones potential. The user is expected
! to pass epsilon, sigma and temperature values. The code also expects
! the discretization along the r axis. A simple trapezoid rule is used
! for numerical integration. Copyright by Miguel A. Caro (2018).
!
! Input info
!
! The code reads in this order and units:
! 
! epsilon(eV) sigma(Angst.) T(K) dr(Angst.) 

  implicit none

  real*8 :: epsilon, sigma, T, dr, pi, kB, r, tol, S, dS, dS_prev
  integer :: i

! Boltzmann's constant in eV/K
  kB = 8.6173303d-5
  pi = dacos(-1.d0)

  read(*,*) epsilon, sigma, T, dr

! Tolerance is by default 1.d-10 times dr. The integration will stop once
! each step contributes less than this amount to the integral
  tol = 1.d-10 * dr

! Initialize
  S = 0.d0
  r = dr/2.d0
  dS_prev = 1.d10
  do
    dS = -2.d0*pi * dr * r**2 * (dexp( 4.d0*epsilon/kB/T * &
         (sigma**6/r**6 - sigma**12/r**12) ) - 1.d0)
!   Break loop if below tolerance
    if( dabs(dS) < tol .and. dabs(dS_prev) < tol ) then
      S = S + dS
      exit
    else
      S = S + dS
      dS_prev = dS
    end if
    r = r + dr
  end do

  write(*,'(A,F10.4,A)') "B2(T) = ", S, " A^3/particle"
  write(*,'(A,F10.4,A)') "      = ", 0.602214d0*S, " cm^3/mol"

end program

