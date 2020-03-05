program hs
! Copyright by Miguel A. Caro (2018)

  use potentials
  use integrators
  use neighbors

  implicit none

  real*8 :: L(1:3), d, Epot, alpha, R, E0, fi(1:3), M, new_pos(1:3), &
            new_vel(1:3), surf, conv_to_bar
  real*8 :: dt, xi_array(1:2, 1:3), T, kB, Ep, Ek, v_int(1:3), d0, dist(1:3), &
            v_mod, pt, P, R_neighbors
  real*8, allocatable :: pos(:, :), vel(:, :), f(:,:), xi_prev(:, :), pts(:)
  logical :: PBC(1:3), container
  integer :: i, Nt, N, j, step, Np, k, p_cycle, write_xyz, pcount, nbuild
  integer, allocatable :: list(:,:), nn(:)

! Our particles are inside a container
  container = .true.

! We read in number of particles, side length and temperature, everything else
! is hard-coded. Writing a better interface to pass input parameters should be
! straightforward.
  read(*,*) Np, L(1), T
  L(2) = L(1)
  L(3) = L(1)

! Some constants
  kB = 8.6173303d-5
  conv_to_bar = 1.6021766208d6

! Number of particles and time step (in fs)
  Nt = 100000
  dt = 1.d0

  nbuild = 100
  R_neighbors = 7.d0
  write_xyz = 100000
  p_cycle = 1000
  allocate( pts(1:p_cycle) )

  if( container )then
    PBC = .false.
  else
    PBC = .true.
  end if
  surf = 2.d0*L(1)*L(2) + 2.d0*L(1)*L(3) + 2.d0*L(2)*L(3)

  R = 2.d0
  alpha = 0.05d0
  E0 = 1.d0
  M = 39.95d0 * 103.62d0

! Safety distance (HS must be further away than this during initialization)
  d0 = 2.d0*R*(1.d0+10.d0*alpha)

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
  vel = dsqrt(3.d0*kB*T / M) * vel

  open(unit=10, file="log.dat", status="unknown")
  open(unit=20, file="trj.xyz", status="unknown")

! Run the dynamics
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
!     If the gas is in a container, the pressure is computed from the force
!     exerted on the walls
      if( container )then
        call hs_container_potential(pos(1:3,i), R, E0, alpha, L, PBC, Epot, fi)
        f(1:3, i) = f(1:3, i) + fi(1:3)
        Ep = Ep + Epot
        pt = pt + dsqrt(dot_product(fi,fi))*dt
      end if
      do k = 1, nn(i)
        j = list(k,i)
        if( j > i )then
          call hs_potential(pos(1:3, i), pos(1:3, j), R, R, E0, alpha, L, PBC, &
                            Epot, fi)
          f(1:3, i) = f(1:3, i) + fi(1:3)
          f(1:3, j) = f(1:3, j) - fi(1:3)
          Ep = Ep + Epot
        end if
      end do
      if( step == 0 )then
        xi_prev(1:3, i) = pos(1:3, i)
        call lazy_man(pos(1:3, i), vel(1:3, i), f(1:3, i), M, dt, new_pos, &
                      new_vel)
      else
        xi_array(1, 1:3) = xi_prev(1:3, i)
        xi_array(2, 1:3) = pos(1:3, i)
        call verlet(xi_array, f(1:3, i), M, dt, new_pos)
      end if
      xi_prev(1:3, i) = pos(1:3, i)
      pos(1:3, i) = new_pos(1:3)
!     Interpolate the velocity to estimate the temperature
      v_int(1:3) = (pos(1:3, i) - xi_prev(1:3, i)) / dt
      Ek = Ek + 0.5d0 * M * dot_product(v_int, v_int)

!     Write trajectory
      if( int(step/write_xyz)*write_xyz == step )then
        if(i == 1)then
          write(20,*) Np
          write(20,*) 'Lattice="', L(1), 0.d0, 0.d0, 0.d0, L(2), 0.d0, 0.d0, &
                      0.d0, L(3), '"'
        end if
        call get_distance(0.5d0*L, pos(1:3, i), L, PBC, new_pos, d)
        write(20,*) "Ar", 0.5d0*L+new_pos(1:3)
      end if
    end do

!   Estimate pressure within cycle
    if( step == 0 )then
      if( container )then
        P = pt/dt/surf
      end if
    else if( step <= p_cycle )then
      if( container )then
        pts(step) = pt
      end if
      P = 0.d0
      do k = 1, step
        P = P + pts(k)
      end do
      if( container )then
        P = P/dfloat(step)/dt/surf
      end if
    else
      P = 0.d0
      do k = 2, p_cycle
        pts(k-1) = pts(k)
        P = P + pts(k)
      end do
      if( container )then
        pts(p_cycle) = pt
      end if
      P = P + pts(p_cycle)
      if( container )then
        P = P/dfloat(p_cycle)/dt/surf
      end if
    end if
!   Write time, potential energy, kinetic energy, instantaneous temperature,
!   pressure within cycle
    write(10,*) dfloat(step)*dt, Ep, Ek, 2.d0/3.d0/dfloat(Np)/kB*Ek, P*conv_to_bar
  end do

  close(10)
  close(20)

end program
