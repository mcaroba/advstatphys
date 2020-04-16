program viscosity_HS

! Copyright by Miguel A. Caro (2018)
!
! Note: this program contains a lot of repeated code from previous sections.
! In future editions of these notes I will simplify the code by separating
! independent chunks into more self-contained modules

  use potentials
  use integrators
  use neighbors
  use random

  implicit none

  real*8 :: L(1:3), d, Epot, alpha, R, E0, fi(1:3), M, new_pos(1:3), new_vel(1:3)
  real*8 :: dt, xi_array(1:2, 1:3), T, kB, Ep, Ek, v_int(1:3), d0, dist(1:3), &
            v_mod, px, R_neighbors, Lpx
  real*8, allocatable :: pos(:, :), vel(:, :), f(:,:), xi_prev(:, :), pts(:), &
                         vx_av(:)
  logical :: PBC(1:3)
  integer :: i, Nt, N, j, step, Np, k, write_xyz, nbuild, Npx, nlayers, &
             n_vx_max, n_vx_min
  integer, allocatable :: list(:,:), nn(:), n_av(:)
  real*8 :: vx_max, vx_min, conv, d1, d2
  integer :: step0_flux, write_log

! Reads in number of particles, side length of box (Angstrom),
! temperature (K) and hard-sphere radius (Angstrom)
  read(*,*) Np, L(1), T, R
  L(2) = L(1)
  L(3) = L(1)

! We divide the simulation box into 20 layers
  nlayers = 20
  Lpx = L(3)/dfloat(nlayers)
  allocate( vx_av(1:nlayers) )
  allocate( n_av(1:nlayers) )
! This is the number of steps between momentum exchange events
! Increase to decrease flow and vice versa
  Npx = 100
  px = 0.d0
! This is the time step when we start monitoring momentum exchange
  step0_flux = 1000000
  vx_av = 0.d0
  n_av = 0

  kB = 8.6173303d-5
! From eV/Angst^3 to kg / m / s^2
  conv = 1.6021766208d-19 / 1.d-10**3

! Total number of time steps and time step (in fs)
  Nt = 2000000
  dt = 1.d0

! Build neighbor list every 100 steps and keep a 7 Angstrom
! radius for neighbor builds
  nbuild = 100
  R_neighbors = 7.d0
  write_xyz = 2000000

  PBC = .true.

  write_log = 100

!  R = 2.d0
  alpha = 0.05d0
  E0 = 1.d0
! This mass is in eV*fs^2/A^2
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
  open(unit=30, file="vx_field.dat", status="unknown")


  do step = 0, Nt
    f = 0.d0
    Ep = 0.d0
    Ek = 0.d0
!   Rebuild list every now and then
    if( int(step/nbuild)*nbuild == step )then
      call build_neighbors(pos, Np, L, PBC, R_neighbors, nn, list)
    end if

!   Transfer momentum
!   NOTE: the code below is only valid for particles of the same mass!!!!!!!!
    if( int(step/Npx)*Npx == step )then
      vx_max = -1.d10
      vx_min = 1.d10
      n_vx_max = 0
      n_vx_min = 0
!     Find particles in bottom and central layers
      do i = 1, Np
        call get_distance((/0.d0, 0.d0, 0.d0 /), pos(1:3, i), L, PBC, dist, d)
!       Find particle in bottom layer with highest vx
        if( dist(3) > 0.d0 .and. dist(3) <= Lpx )then
          v_int(1:3) = (pos(1:3, i) - xi_prev(1:3, i)) / dt
          if( v_int(1) > vx_max )then
            vx_max = v_int(1)
            n_vx_max = i
          end if
!       Find particle in central layer with lowest vx
        else if( dist(3) > -L(3)/2.d0 .and. dist(3) <= -L(3)/2.d0+Lpx )then
          v_int(1:3) = (pos(1:3, i) - xi_prev(1:3, i)) / dt
          if( v_int(1) < vx_min )then
            vx_min = v_int(1)
            n_vx_min = i
          end if
        end if
      end do
!     Exchange momenta between particle in bottom and central layers
      if( n_vx_max /= 0 .and. n_vx_min /= 0 )then
        i = n_vx_max
        j = n_vx_min

!       d1 is the momentum of the particle in the bottom layer
        d1 = pos(1,i) - xi_prev(1,i)
!       d2 is the momentum of the particle in the middle layer
        d2 = pos(1,j) - xi_prev(1,j)

!       We subtract d1 and add d2 to the momentum of i
        xi_prev(1, i) = pos(1, i) - d2
!       We subtract d2 and add d1 to the momentum of j
        xi_prev(1, j) = pos(1, j) - d1

!       We do this only after step0_flux-1 steps, before that we are
!       equilibrating the system
        if( step >= step0_flux )then
!         Flux is measured towards the middle layer (positive z direction),
!         therefore we add d1 and subtract d2 twice to obtain the momentum
!         we have transferred. Because we are reversing the process, in the
!         physical simulation the momentum would flow in the opposite
!         direction. Therefore we change the sign
          px = px - M*2.d0*(d1-d2)/dt
        end if
      end if
    end if

!   Force evaluation and integration
    do i = 1, Np
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

!   Calculate vx field
    if( step >= step0_flux )then
      do i = 1, Np
        call get_distance(0.5d0*L, pos(1:3, i), L, PBC, dist, d)
        d = dist(3) + L(3)/2.d0
        j = 1 + int(d/Lpx)
        vx_av(j) = vx_av(j) + (pos(1, i) - xi_prev(1, i)) / dt
        n_av(j) = n_av(j) + 1
      end do
    end if

!   Write time, potential energy, kinetic energy, instantaneous temperature
    if( int(step/write_log)*write_log == step )then
      write(10,*) dfloat(step)*dt, Ep, Ek, 2.d0/3.d0/dfloat(Np)/kB*Ek
    end if
  end do

  close(10)
  close(20)


! Write vx field
  vx_av = vx_av / dfloat(n_av)
  do i = 1, nlayers
    write(30, *) Lpx*(dfloat(i)-0.5d0), vx_av(i)
  end do
  write(30, *) Lpx*(dfloat(nlayers+1)-0.5d0), vx_av(1)
  close(30)


! We divide the momentum by two because half is transferred up and half down
! (because of the PBC)
  write(*,*) "Momentum flux = ", px/2.d0/dfloat(Nt-step0_flux)/dt/L(1)/L(2), &
             "eV / Angst^3"
  write(*,*) "              = ", px/2.d0/dfloat(Nt-step0_flux)/dt/L(1)/L(2)* & 
             conv, "kg / m / s^2"

end program
