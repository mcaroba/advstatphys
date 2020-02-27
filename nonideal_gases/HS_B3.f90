! Copyright by Miguel A. Caro (2018).

program integrate_MC_spherical_B3

  implicit none

  real*8, allocatable :: random(:,:)
  real*8 :: r12, r13, r23, phi12, phi13, t12, t13, int, pi, weight_in, weight_tot
  integer :: N, N_in, i, j

  pi = dacos(-1.d0)

  N = 1d6
  N_in = 0

  int = 0.d0

  do i = 1, 1000

    weight_in = 0.d0
    weight_tot = 0.d0

    allocate( random(1:6, 1:N) )
    call init_random_seed()
    call random_number(random)

    do j = 1, N
      r12 = random(1, j)
      r13 = random(2, j)
      phi12 = random(3,j) * pi
      phi13 = random(4,j) * pi
      t12 = random(5,j) * 2.d0*pi
      t13 = random(6,j) * 2.d0*pi

      r23 = dsqrt( r12**2 + r13**2 - 2.d0*r12*r13*(dsin(phi12)* &
            dsin(phi13)*dcos(t12-t13) + dcos(phi12)*dcos(phi13)) )

      if( r23 < 1.d0 )then
        N_in = N_in + 1
        weight_in = weight_in + dsin(phi12)*dsin(phi13)*r12**2*r13**2
      end if
      weight_tot = weight_tot + dsin(phi12)*dsin(phi13)*r12**2*r13**2
    end do

    int = int + 16.d0*pi**2/9.d0 * weight_in/weight_tot /3.d0

    write(*,*) N*i, N_in, int/dfloat(i), int/dfloat(i)-5.d0/18.d0*pi**2

    deallocate( random )
  end do

end program


! The code below is NOT copyright by Miguel A. Caro
subroutine init_random_seed()

!*******************************************************************************
! This subroutine was copy-pasted from the gfortran online documentation. It
! initializes the random number generator.
!*******************************************************************************

  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t

   call random_seed(size = n)
   allocate(seed(n))
   ! First try if the OS provides a random number generator
   open(newunit=un, file="/dev/urandom", access="stream", &
     form="unformatted", action="read", status="old", iostat=istat)
   if (istat == 0) then
     read(un) seed
     close(un)
   else
! Fallback to XOR:ing the current time and pid. The PID is
! useful in case one launches multiple instances of the same
! program in parallel.
     call system_clock(t)
     if (t == 0) then
       call date_and_time(values=dt)
       t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
           + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
           + dt(3) * 24_int64 * 60 * 60 * 1000 &
           + dt(5) * 60 * 60 * 1000 &
           + dt(6) * 60 * 1000 + dt(7) * 1000 &
           + dt(8)
     end if
     pid = getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
       seed(i) = lcg(t)
     end do
   end if
   call random_seed(put=seed)
   contains
   ! This simple PRNG might not be good enough for real work, but is
   ! sufficient for seeding a better PRNG.
     function lcg(s)
       integer :: lcg
       integer(int64) :: s
       if (s == 0) then
         s = 104729
       else
         s = mod(s, 4294967296_int64)
       end if
         s = mod(s * 279470273_int64, 4294967291_int64)
         lcg = int(mod(s, int(huge(0), int64)), kind(0))
       end function lcg
end subroutine init_random_seed
