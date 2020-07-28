! Copyright by Miguel A. Caro (2018).

module potentials

  implicit none

  contains

!**********************************************************************************************
! This subroutine returns the distance between ri and rj under
! certain boundary conditions
  subroutine get_distance(posi, posj, L, PBC, dist, d)

    implicit none

    real*8, intent(in) :: posi(1:3), posj(1:3), L(1:3)
    logical, intent(in) :: PBC(1:3)
    real*8, intent(out) :: d
    real*8, intent(out) :: dist(1:3)
    real*8 :: d2
    integer :: i

    d2 = 0.d0
    do i = 1, 3
      if( PBC(i) )then
        dist(i) = modulo(posj(i) - posi(i), L(i))
        if( dist(i) > L(i)/2.d0 )then
          dist(i) = dist(i) - L(i)
        end if
      else
        dist(i) = posj(i) - posi(i)
      end if
      d2 = d2 + dist(i)**2
    end do
    d = dsqrt(d2)

  end subroutine get_distance
!**********************************************************************************************




!**********************************************************************************************
! This returns potential energy and force for a 1/r type potential. G is a
! constant prefactor; it could be the gravity constant or e^2/(4 pi eps0), etc.
  subroutine pairwise_electrostatic_potential(posi, posj, Zi, Zj, G, L, PBC, &
                                              Epot, fi)

    implicit none

    real*8, intent(in) :: posi(1:3), posj(1:3), Zi, Zj, G, L(1:3)
    logical, intent(in) :: PBC(1:3)
    real*8, intent(out) :: Epot, fi(1:3)
    real*8 :: d, dist(1:3)

    call get_distance(posi, posj, L, PBC, dist, d)

    Epot = G * Zi*Zj / d

!   The force on i is calculated assuming the convention that dist(1) = xj - xi
    fi(1:3) = - G * Zi*Zj * dist(1:3)/d**3.d0

  end subroutine
!**********************************************************************************************




!**********************************************************************************************
! This subroutine returns the potential energy due to 2 hard spheres interacting
  subroutine hs_potential(posi, posj, Ri, Rj, E0, alpha, L, PBC, Epot, fi)

    implicit none

    real*8, intent(in) :: posi(1:3), posj(1:3), Ri, Rj, E0, alpha, L(1:3)
    logical, intent(in) :: PBC(1:3)
    real*8, intent(out) :: Epot, fi(1:3)
    real*8 :: d, decay_length, pi, dist(1:3)

    pi = dacos(-1.d0)

    decay_length = alpha * (Ri + Rj)
    
    call get_distance(posi, posj, L, PBC, dist, d)

    Epot = E0/2. * erfc((d - Ri - Rj) / decay_length)

!   The force on i is calculated assuming the convention that dist(1) = xj - xi
    fi(1:3) = - E0/dsqrt(pi) * dexp(-((d-Ri-Rj)/decay_length)**2) * dist(1:3) &
              / d / decay_length

  end subroutine hs_potential
!**********************************************************************************************




!**********************************************************************************************
! This subroutine returns the potential energy due to a hard sphere interacting
! with the walls of a container
  subroutine hs_container_potential(posi, Ri, E0, alpha, L, Epot, fi)

    implicit none

    real*8, intent(in) :: posi(1:3), Ri, E0, alpha, L(1:3)
    real*8, intent(out) :: Epot, fi(1:3)
    real*8 :: d(1:6), decay_length, pi, dist(1:3)
    integer :: i

    pi = dacos(-1.d0)

    decay_length = alpha * Ri
    
    d(1) = -posi(1)
    d(2) = -posi(2)
    d(3) = -posi(3)
    d(4) = -posi(1) + L(1)
    d(5) = -posi(2) + L(2)
    d(6) = -posi(3) + L(3)

    Epot = 0.d0
    do i = 1, 6
      Epot = Epot + E0/2. * erfc((dabs(d(i)) - Ri) / decay_length)
    end do

!   The force on i is calculated assuming the convention that dist(1) = xj - xi
    fi = 0.d0
    do i = 1, 3
      fi(i) = - E0/dsqrt(pi) * dexp(-((dabs(d(i))-Ri)/decay_length)**2) * d(i) &
                / dabs(d(i)) / decay_length
      fi(i) = fi(i) - E0/dsqrt(pi) * dexp(-((dabs(d(i+3))-Ri)/decay_length)**2) &
              * d(i+3) / dabs(d(i+3)) / decay_length
    end do

  end subroutine hs_container_potential
!**********************************************************************************************




!**********************************************************************************************
! This subroutine returns the interaction energy between two particles according to the
! Lennard-Jones potential
  subroutine lj_potential(posi, posj, sigmai, sigmaj, epsiloni, epsilonj, &
                          Rcut, L, PBC, Epot, fi)

    implicit none

    real*8, intent(in) :: posi(1:3), posj(1:3), sigmai, sigmaj, &
                          epsiloni, epsilonj, L(1:3), Rcut
    logical, intent(in) :: PBC(1:3)
    real*8, intent(out) :: Epot, fi(1:3)
    real*8 :: d, epsilon, sigma, pi, dist(1:3), E0

    pi = dacos(-1.d0)

    sigma = 0.5d0 * (sigmai + sigmaj)
    epsilon = dsqrt(epsiloni*epsilonj)

    E0 = 4.d0*epsilon * ( (sigma/Rcut)**12 - (sigma/Rcut)**6 )
    
    call get_distance(posi, posj, L, PBC, dist, d)

    if( d < Rcut )then
      Epot = 4.d0*epsilon * ( (sigma/d)**12 - (sigma/d)**6 ) - E0
!     The force on i is calculated assuming the convention that dist(1) = xj - xi
      fi(1:3) = 4.d0*epsilon/d**2 * (-12.d0*(sigma/d)**12 + &
                6.d0*(sigma/d)**6) * dist(1:3)
    else
!     There is no discontinuity of the potential at Rcut, but there
!     is a discontinuity of the force
      Epot = 0.d0
      fi(1:3) = 0.d0
    end if

  end subroutine lj_potential
!**********************************************************************************************



end module
