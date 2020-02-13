module potentials

  implicit none

  contains

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

end module
