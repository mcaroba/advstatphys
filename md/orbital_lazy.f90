program orbital

  use potentials
  use integrators

  implicit none

  real*8 :: pos(1:3, 1:3), vel(1:3, 1:3), m(1:3), dt, L(1:3), E_sum, Eij
  real*8 :: fi_sum(1:3, 1:3), fi(1:3), fi_prev(1:3, 1:3), fi_array(1:2, 1:3)
  real*8 :: xi_prev(1:3, 1:3), xi_array(1:2, 1:3), new_pos(1:3), new_vel(1:3)
  integer :: step, n, i, j
  logical :: PBC(1:3)

! Sun initial conditions
  pos(1,1:3) = (/ 0.d0, 0.d0, 0.d0 /)
  vel(1,1:3) = (/ 0.d0, -0.335d0, 0.d0 /)
  m(1) = 10.d0

! Earth initial conditions
  pos(2,1:3) = (/ -1.d0, 0.d0, 0.d0 /)
  vel(2,1:3) = (/ 0.d0, 4.d0, 0.d0 /)
  m(2) = 1.d0

! Moon initial conditions
  pos(3,1:3) = (/ -0.95d0, 0.d0, 0.d0 /)
  vel(3,1:3) = (/ 0.d0, -1.d0, 0.d0 /)
  m(3) = 0.1d0


! Time step
  dt = 1.d-3
  n = 10000

! BC
  PBC = .false.
  L = (/ 1.d0, 1.d0, 1.d0 /)


  open(unit=10, file="orbit_lazy", status="unknown")
! Run for n time steps
  do step = 0, n
    E_sum = 0.d0
    fi_sum = 0.d0
    do i = 1, 3
      do j = i+1, 3
        call pairwise_electrostatic_potential(pos(i, 1:3), pos(j, 1:3), m(i), &
                                              m(j), -1.d0, L, PBC, Eij, fi(1:3))
        E_sum = E_sum + Eij
        fi_sum(i, 1:3) = fi_sum(i, 1:3) + fi(1:3)
        fi_sum(j, 1:3) = fi_sum(j, 1:3) - fi(1:3)
      end do
    end do
    write(10,*) dfloat(step)*dt, E_sum, pos(1, 1:2), pos(2, 1:2), pos(3, 1:2)
    do i = 1, 3
      call lazy_man(pos(i, 1:3), vel(i, 1:3), fi_sum(i, 1:3), m(i), dt, &
                    new_pos, new_vel)
      pos(i, 1:3) = new_pos(1:3)
      vel(i, 1:3) = new_vel(1:3)
    end do  
  end do
  close(10)
  
end program
