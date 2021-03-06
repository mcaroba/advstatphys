module integrators

  implicit none

  contains

! Lazy Man's approach
  subroutine lazy_man(x_in, v_in, F, m, dt, x_out, v_out)

    implicit none

    real*8, intent(in) :: x_in(1:3), v_in(1:3)
    real*8, intent(out) :: x_out(1:3), v_out(1:3)
    real*8, intent(in) :: F(1:3), m, dt

    x_out(1:3) = x_in(1:3) + v_in(1:3)*dt + F(1:3)/m/2.d0*dt**2
    v_out(1:3) = v_in(1:3) + F(1:3)/m*dt

  end subroutine


! Regular Verlet
  subroutine verlet(x_in, F, m, dt, x_out)

    implicit none

    real*8, intent(in) :: x_in(1:2,1:3)
    real*8, intent(out) :: x_out(1:3)
    real*8, intent(in) :: F(1:3), m, dt

    x_out(1:3) = 2.d0*x_in(2,1:3) - x_in(1,1:3) + F(1:3)/m*dt**2

  end subroutine


! Velocity Verlet is two subroutines
  subroutine velocity_verlet_vel(v_in, F, m, dt, v_out)

    implicit none

    real*8, intent(in) :: v_in(1:3)
    real*8, intent(out) :: v_out(1:3)
    real*8, intent(in) :: F(1:2,1:3), m, dt

    v_out(1:3) = v_in(1:3) + (F(2,1:3) + F(1,1:3))/2.d0/m*dt

  end subroutine
  subroutine velocity_verlet_pos(x_in, v_in, F, m, dt, x_out)

    implicit none

    real*8, intent(in) :: x_in(1:3), v_in(1:3)
    real*8, intent(out) :: x_out(1:3)
    real*8, intent(in) :: F(1:3), m, dt

    x_out(1:3) = x_in(1:3) + v_in(1:3)*dt + F(1:3)/m/2.d0*dt**2

  end subroutine




  subroutine berendsen_thermostat(vel, T0, T, tau, dt)

    implicit none

    real*8, intent(inout) :: vel(:,:)
    real*8, intent(in) :: T0, T, tau, dt

    vel = vel * dsqrt(1.d0 + dt/tau * (T0/T - 1.d0))

  end subroutine




  subroutine berendsen_barostat(vector_3d, P0, P, tau, dt)
!   Berendsen barostat that uses the compressibility of water
!   and takes P in bar
    implicit none

    real*8, intent(inout) :: vector_3d(:)
    real*8, intent(in) :: P0, P, tau, dt

    vector_3d = vector_3d * (1.d0 + dt/tau * 4.5d-5 * (P - P0))**(1.d0/3.d0)

  end subroutine





  subroutine remove_cm_vel(vel, M)

!   I should adapt this code to mixed boundary conditions, where
!   the CM velocity can be removed per Cartesian dimension independently

    implicit none

    real*8, intent(inout) :: vel(:,:)
    real*8, intent(in) :: M(:)
    real*8 :: cm_pos(1:3), cm_vel(1:3), total_mass
    integer :: Np, i

    Np = size(vel, 2)

    cm_vel = 0.d0
    total_mass = 0.d0
    do i = 1, Np
      cm_vel(1:3) = cm_vel(1:3) + M(i)*vel(1:3,i)
      total_mass = total_mass + M(i)
    end do
    cm_vel = cm_vel / total_mass
    do i = 1, Np
      vel(1:3,i) = vel(1:3,i) - cm_vel(1:3)
    end do

  end subroutine




  subroutine scale_box(pos, L, init_L, final_L, t, init_t, final_t)

    implicit none

    real*8, intent(inout) :: pos(:,:), L(:)
    real*8, intent(in) :: init_L(:), final_L(:)
    integer, intent(in) :: t, init_t, final_t
    real*8 :: f(1:3)
    integer :: natoms, i

!   We scale the box and positions only if within the rescaling interval
    if( t > init_t .and. t <= final_t )then
      f(1:3) = ( init_L(1:3) + (final_L(1:3) - init_L(1:3)) / dfloat(final_t - init_t) * &
                 dfloat(t - init_t) ) / L(1:3)
      L(1:3) = L(1:3) * f(1:3)
      natoms = size(pos,2)
      do i = 1, natoms
        pos(1:3, i) = pos(1:3, i) * f(1:3)
      end do
    end if

  end subroutine

end module
