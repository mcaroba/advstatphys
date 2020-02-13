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

end module
