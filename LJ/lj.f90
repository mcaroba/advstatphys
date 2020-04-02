program lj

! Copyright by Miguel A. Caro (2018)

  use potentials
  use integrators
  use neighbors

  implicit none

  real*8 :: L(1:3), d, Epot, fi(1:3), new_pos(1:3), new_vel(1:3), surf, conv_to_bar
  real*8 :: dt, xi_array(1:2, 1:3), T, kB, Ep, Ek, d0, dist(1:3), v_mod, pt, P, R_neighbors
  real*8, allocatable :: pos(:, :), vel(:, :), f(:,:), xi_prev(:, :), pts(:), M(:)
  logical :: PBC(1:3), container
  integer :: i, Nt, N, j, step, Np, k, p_cycle, write_xyz, pcount, nbuild, write_log
  integer, allocatable :: list(:,:), nn(:)

  real*8 :: epsilon, sigma, Rcut, dist2(1:3), f_img(-1:1, -1:1, -1:1), tau
  integer :: i2, j2, k2

  container = .false.

!  L = (/ 100.d0, 100.d0, 100.d0 /)
  read(*,*) Np, L(1), T
  L(2) = L(1)
  L(3) = L(1)

  kB = 8.6173303d-5
  conv_to_bar = 1.6021766208d6

  Nt = 1000000
!  Np = 500
!  T = 300.d0
  dt = 2.d0
  tau = 100.d0

  nbuild = 100
  R_neighbors = 15.d0
  write_xyz = 100
  write_log = 100
  p_cycle = 1000
  allocate( pts(1:p_cycle) )

  if( container )then
    PBC = .false.
  else
    PBC = .true.
  end if
  surf = 2.d0*L(1)*L(2) + 2.d0*L(1)*L(3) + 2.d0*L(2)*L(3)

  epsilon = 0.0117d0
  sigma = 3.235d0
  Rcut = 12.d0
  allocate( M(1:Np) )
  M = 39.95d0 * 103.62d0

! Safety distance (particles must be further away than this distance during initialization)
  d0 = 1.5d0*sigma

  allocate( pos(1:3, 1:Np) )
  allocate( vel(1:3, 1:Np) )
  allocate( f(1:3, 1:Np) )
  allocate( xi_prev(1:3, 1:Np) )
  allocate( nn(1:Np) )
  allocate( list(1:Np, 1:Np) )

! Initialize random number generator
  call init_random_seed()

! Make sure that hard spheres are initially far away from each other
  do i = 1, Np
!   Give particle i an initial postion
    call random_number(new_pos)
    do k = 1, 3
      pos(k, i) = d0/2.d0 + (L(k)-d0)*new_pos(k)
    end do
!   Check that this initial position is at least d0 away from any other sphere
    j = 1
    do while(j < i)
      call get_distance(pos(1:3, i), pos(1:3, j), L, PBC, dist, d)
      if( d > d0 )then
        j = j + 1
      else
        j = 1
        call random_number(new_pos)
        do k = 1, 3
          pos(k, i) = d0/2.d0 + (L(k)-d0)*new_pos(k)
        end do
      end if
    end do
  end do

! Initialize velocities by randomizing the direction
  call random_number(vel)
  vel = 2.d0*vel-1.d0
  do i = 1, Np
    v_mod = dsqrt(dot_product(vel(1:3, i), vel(1:3, i)))
    vel(1:3, i) = vel(1:3, i) / v_mod
  end do
  if( .not. container )then
    do i = 1,3
      vel(i,:) = dsqrt(3.d0*kB*T / M(:) * dfloat(Np-1)/dfloat(Np) ) * vel(i,:)
    end do
    call remove_cm_vel(vel, M)
  else
    do i = 1, 3
      vel(i,:) = dsqrt(3.d0*kB*T / M(:)) * vel(i,:)
    end do
  end if


  open(unit=10, file="log.dat", status="unknown")
  open(unit=20, file="trj.xyz", status="unknown")





  do step = 0, Nt
    f = 0.d0
    Ep = 0.d0
    Ek = 0.d0
    pt = 0.d0
!   Rebuild list every now and then
    if( int(step/nbuild)*nbuild == step )then
      call build_neighbors(pos, Np, L, PBC, R_neighbors, nn, list)
    end if
    do i = 1, Np
!==================================================================================
!     If the gas is in a container, the pressure is computed from the force exerted on the
!     walls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Figure out a way to do LJ with the container
!      if( container )then
!        call hs_container_potential(pos(1:3,i), R, E0, alpha, L, PBC, Epot, fi)
!        f(1:3, i) = f(1:3, i) + fi(1:3)
!        Ep = Ep + Epot
!        pt = pt + dsqrt(dot_product(fi,fi))*dt
!      end if
!==================================================================================
      do k = 1, nn(i)
        j = list(k,i)
        if( j > i )then
!          call hs_potential(pos(1:3, i), pos(1:3, j), R, R, E0, alpha, L, PBC, Epot, fi)
          call lj_potential(pos(1:3, i), pos(1:3, j), sigma, sigma, epsilon, epsilon, Rcut, L, PBC, Epot, fi)
          f(1:3, i) = f(1:3, i) + fi(1:3)
          f(1:3, j) = f(1:3, j) - fi(1:3)
          Ep = Ep + Epot
        end if
      end do
      if( step == 0 )then
        xi_prev(1:3, i) = pos(1:3, i)
        call lazy_man(pos(1:3, i), vel(1:3, i), f(1:3, i), M(i), dt, new_pos, new_vel)
      else
        xi_array(1, 1:3) = xi_prev(1:3, i)
        xi_array(2, 1:3) = pos(1:3, i)
        call verlet(xi_array, f(1:3, i), M(i), dt, new_pos)
      end if
      xi_prev(1:3, i) = pos(1:3, i)
      pos(1:3, i) = new_pos(1:3)
!     Interpolate the velocity to estimate the temperature
      vel(1:3, i) = (pos(1:3, i) - xi_prev(1:3, i)) / dt
      Ek = Ek + 0.5d0 * M(i) * dot_product(vel(1:3, i), vel(1:3, i))

!==================================================================================
!     I can't figure out how to correctly implement the pressure calculation
!     under periodic boundary conditions. I'll leave this for the future
!     If there is no container, the pressure is computed like this
!      if( .not. container )then
!       Position of image of particle i which is closest to middle of the cell
!        call get_distance(0.5d0*L, pos(1:3, i), L, PBC, dist, d)
!       The absolute position of copy i is
!        new_pos = 0.5d0*L(1:3) + dist(1:3)
!        pt = pt + dot_product(new_pos(1:3), f(1:3,i))
!        do k = 1, nn(i)
!          j = list(k,i)
!          if( j > i )then
!           Position of image of particle j which is closest to 1, wrt 1
!            call get_distance(pos(1:3, i), pos(1:3, j), L, PBC, dist2, d)
!           The absolute position of copy j is
!            new_pos = 0.5d0*L(1:3) + dist(1:3) + dist2(1:3)
!           Add this one only if copy is outside the local cell
!            i2 = floor(new_pos(1)/L(1))
!            j2 = floor(new_pos(2)/L(2))
!            k2 = floor(new_pos(3)/L(3))
!            if( i2 /= 0 .or. j2 /= 0 .or. k2 /= 0 )then
!              call lj_potential(pos(1:3, i), pos(1:3, j), sigma, sigma, epsilon, epsilon, Rcut, L, PBC, Epot, fi)
!              new_pos(1:3) = (/ dfloat(i2)*L(1), dfloat(j2)*L(2), dfloat(k2)*L(3) /)
!              pt = pt + dot_product(new_pos(1:3), -fi(1:3))
!            end if
!          end if
!        end do
!      end if
!==================================================================================
    end do

!   Remove CM velocity (should be close to zero since we initialized properly)
    call remove_cm_vel(vel, M)
!   Rescale the velocities to control the temperature
    call berendsen_thermostat(vel, T, 2.d0/3.d0/dfloat(Np-1)/kB*Ek, tau, dt)
!   Readjust the positions according to new velocities
    pos = xi_prev + vel*dt

!   Write trajectory
    do i = 1, Np
      if( int(step/write_xyz)*write_xyz == step )then
        if(i == 1)then
          write(20,*) Np
          write(20,*) 'Lattice="', L(1), 0.d0, 0.d0, 0.d0, L(2), 0.d0, 0.d0, 0.d0, L(3), '"'
        end if
        call get_distance(0.5d0*L, pos(1:3, i), L, PBC, new_pos, d)
        write(20,*) "Ar", 0.5d0*L+new_pos(1:3)
      end if
    end do


!   Estimate pressure within cycle
    if( step == 0 )then
      if( container )then
        P = pt/dt/surf
      else
        P = (2.d0*Ek + pt) / 3.d0 / (L(1)*L(2)*L(3))
      end if
    else if( step <= p_cycle )then
      if( container )then
        pts(step) = pt
      else
        pts(step) = pt + 2.d0*Ek
      end if
      P = 0.d0
      do k = 1, step
        P = P + pts(k)
      end do
      if( container )then
        P = P/dfloat(step)/dt/surf
      else
        P = P/dfloat(step) / 3.d0 / (L(1)*L(2)*L(3)) 
      end if
    else
      P = 0.d0
      do k = 2, p_cycle
        pts(k-1) = pts(k)
        P = P + pts(k)
      end do
      if( container )then
        pts(p_cycle) = pt
      else
        pts(p_cycle) = pt + 2.d0*Ek
      end if
      P = P + pts(p_cycle)
      if( container )then
        P = P/dfloat(p_cycle)/dt/surf
      else
        P = P/dfloat(p_cycle) / 3.d0 / (L(1)*L(2)*L(3))
      end if
    end if
!   Write time, potential energy, kinetic energy, instantaneous temperature, pressure within cycle
    if( int(step/write_log)*write_log == step )then
      write(10,*) dfloat(step)*dt, Ep, Ek, 2.d0/3.d0/dfloat(Np-1)/kB*Ek, P*conv_to_bar
    end if
  end do

  close(10)
  close(20)

end program

























!=================================================================================================
!=================================================================================================
subroutine init_random_seed()

!***********************************************************************************
! This subroutine was copy-pasted from the gfortran online documentation. It
! initializes the random number generator.
!***********************************************************************************

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
!=================================================================================================
!=================================================================================================
