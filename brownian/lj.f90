program lj

! Copyright by Miguel A. Caro (2020)

  use potentials
  use integrators
  use neighbors
  use analyze

  implicit none

  real*8 :: L(1:3), d, Epot, fi(1:3), new_pos(1:3), new_vel(1:3), conv_to_bar
  real*8 :: dt, xi_array(1:2, 1:3), T, kB, Ep, Ek, d0, dist(1:3), v_mod, R_neighbors
  real*8, allocatable :: pos(:, :), vel(:, :), f(:,:), xi_prev(:, :), M(:)
  logical :: PBC(1:3)
  integer :: i, Nt, N, j, step, Np, k, write_xyz, nbuild, write_log
  integer, allocatable :: list(:,:), nn(:)

  real*8 :: dist2(1:3), f_img(-1:1, -1:1, -1:1), tau = 100.d0
  integer :: i2, j2, k2

  character*128 :: keyword, crap
  character*1 :: keyword_first
  integer :: iostatus, l_max = -1, i_ub
  real*8 :: crap_v3(1:3), rdf_min = 4.5d0
  real*8 :: V, virial, P, update_bar, tau_p = 100.d0, P0
  logical :: rescale = .false., restart = .false.
  real*8, allocatable :: Ql(:)

  integer :: Nspecies = 1
  real*8, allocatable :: M0(:), epsilon(:), sigma(:), Rcut(:), pos0(:,:)
  character*128, allocatable :: species(:)
  logical :: planar = .false., wrap = .true.
  integer, allocatable :: which_species(:), formula(:)
  real*8 :: msd, msd_std


  write(*,*)
  write(*,'(1X,A)') "*************************************"
  write(*,'(1X,A)') "  Welcome to the CBUF MD simulator   "
  write(*,*)        "   (CBUF = Crappy BUt Functional)    "
  write(*,*)

  PBC = .true.

! Fundamental physical constants
  kB = 8.6173303d-5
  conv_to_bar = 1.6021766208d6

! Options for the neighbor list interface and writing out
! Ideally there would be a runtime check on whether the
! lists are built too often or, especially, not often enough
  nbuild = 100
  R_neighbors = 15.d0
  write_xyz = 100
  write_log = 100


! Read input file
  open( unit = 10, file = "input", status = "old" )
  iostatus = 0
  do while( iostatus == 0 )
    read(10, *, iostat=iostatus) keyword_first
    if( keyword_first == "!" .or. keyword_first == "#" )then
      cycle
    else
      backspace(10)
      read(10, *, iostat=iostatus) keyword
      if( keyword == "T" )then
        backspace(10)
        read(10, *) keyword, keyword, T
      else if( keyword == "P" )then
        backspace(10)
        read(10, *) keyword, keyword, P0
      else if( keyword == "Np" )then
        backspace(10)
        read(10, *) keyword, keyword, Np
      else if( keyword == "Nt" )then
        backspace(10)
        read(10, *) keyword, keyword, Nt
      else if( keyword == "dt" )then
        backspace(10)
        read(10, *) keyword, keyword, dt
      else if( keyword == "restart" )then
        backspace(10)
        read(10, *) keyword, keyword, restart
      else if( keyword == "L" )then
        backspace(10)
        read(10, *) keyword, keyword, L(1)
        L(2:3) = L(1)
      else if( keyword == "tau" )then
        backspace(10)
        read(10, *) keyword, keyword, tau
      else if( keyword == "tau_p" )then
        backspace(10)
        read(10, *) keyword, keyword, tau_p
      else if( keyword == "rescale" )then
        backspace(10)
        read(10, *) keyword, keyword, rescale
      else if( keyword == "rdf_min" )then
        backspace(10)
        read(10, *) keyword, keyword, rdf_min
      else if( keyword == "write_log" )then
        backspace(10)
        read(10, *) keyword, keyword, write_log
      else if( keyword == "write_xyz" )then
        backspace(10)
        read(10, *) keyword, keyword, write_xyz
      else if( keyword == "Nspecies" )then
        backspace(10)
        read(10, *) keyword, keyword, Nspecies
      else if( keyword == "planar" )then
        backspace(10)
        read(10, *) keyword, keyword, planar
      else if( keyword == "wrap" )then
        backspace(10)
        read(10, *) keyword, keyword, wrap
      else if( keyword == "l_max" )then
        backspace(10)
        read(10, *) keyword, keyword, l_max
        if( l_max >= 0 )then
          allocate( Ql(0:l_max) )
        end if
      end if
    end if
  end do
  close(10)


  allocate( formula(1:Nspecies) )
  formula(1) = Np
! Read input file for formula only
  open( unit = 10, file = "input", status = "old" )
  iostatus = 0
  do while( iostatus == 0 )
    read(10, *, iostat=iostatus) keyword_first
    if( keyword_first == "!" .or. keyword_first == "#" )then
      cycle
    else
      backspace(10)
      read(10, *, iostat=iostatus) keyword
      if( keyword == "formula" )then
        backspace(10)
        read(10, *) keyword, keyword, formula(1:Nspecies)
      end if
    end if
  end do
  close(10)



! Read potential parameter file
  allocate( species(1:Nspecies) )
  allocate( epsilon(1:Nspecies) )
  allocate( sigma(1:Nspecies) )
  allocate( Rcut(1:Nspecies) )
  allocate( M0(1:Nspecies) )
  open( unit = 10, file = "pot.param", status = "old" )
! Read in the values
  iostatus = 0
  do while( iostatus == 0 )
    read(10, *, iostat=iostatus) keyword_first
    if( keyword_first == "!" .or. keyword_first == "#" )then
      cycle
    else
      backspace(10)
      read(10, *, iostat=iostatus) keyword
      if( keyword == "epsilon" )then
        backspace(10)
        read(10, *) keyword, keyword, epsilon(1:Nspecies)
      else if( keyword == "sigma" )then
        backspace(10)
        read(10, *) keyword, keyword, sigma(1:Nspecies)
      else if( keyword == "Rcut" )then
        backspace(10)
        read(10, *) keyword, keyword, Rcut(1:Nspecies)
      else if( keyword == "M0" )then
        backspace(10)
        read(10, *) keyword, keyword, M0(1:Nspecies)
      else if( keyword == "species" )then
        backspace(10)
        read(10, *) keyword, keyword, species(1:Nspecies)
      end if
    end if
  end do
  close(10)



! If this is a restart calculation, we handle it here
! We only read the last frame
  if( restart )then
    open(unit=10, file="trj.xyz", status="old")
    iostatus = 0
    do while( iostatus == 0 )
      read(10, *, iostat=iostatus) Np
      do i = 1, Np+1
        read(10, *, iostat=iostatus)
      end do
    end do
  end if


! Allocate all arrays
  allocate( pos(1:3, 1:Np) )
  allocate( pos0(1:3, 1:Np) )
  allocate( vel(1:3, 1:Np) )
  allocate( f(1:3, 1:Np) )
  allocate( xi_prev(1:3, 1:Np) )
  allocate( nn(1:Np) )
  allocate( list(1:Np, 1:Np) )
  allocate( which_species(1:Np) )
! The mass M0 needs to be read from the pot.param file
  allocate( M(1:Np) )


  if( restart )then
    do i = 1, Np+3
      backspace(10)
    end do
    read(10, *) Np
    read(10, *) crap, L(1), crap, crap, crap, L(2), crap, crap, crap, L(3)
    do i = 1, Np
      read(10,*) crap, pos(1:3, i), vel(1:3, i)
      do j = 1, Nspecies
        if( crap == species(j) )then
          which_species(i) = j
          M(i) = M0(j) * 103.62d0
          exit
        end if
      end do
    end do
    close(10)
  else
    k = 1
    do i = 1, Nspecies
      do j = 1, formula(i)
        which_species(k) = i
        M(k) = M0(i) * 103.62d0
        k = k + 1
      end do
    end do
    write(*,*)
    write(*,'(1X,A)') "Generating starting configuration:"
    write(*,*)
    write(*,'(1X,A)') ">...................................<"
    write(*,'(1X,A)',advance='no') '['
    update_bar = dfloat(Np)/36.d0
    i_ub = 1
!   If this is a new configuration, randomize
!   Safety distance (particles must be further away than this distance during initialization)
    d0 = 2.d0
!   Initialize random number generator
    call init_random_seed()
!   Make sure that atoms are initially far away from each other
    do i = 1, Np
!     Update progress bar every Np/36 iterations
      if(dfloat(i) > dfloat(i_ub)*update_bar)then
        write(*,'(A)', advance='no') '='
        i_ub = i_ub + 1
      end if
!     Give particle i an initial postion
      call random_number(new_pos)
      do k = 1, 3
        pos(k, i) = d0/2.d0 + (L(k)-d0)*new_pos(k)
        if( planar )then
          pos(3, i) = L(3) / 2.d0
        end if
      end do
!     Check that this initial position is at least d0 away from any other sphere
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
            if( planar )then
              pos(3, i) = L(3) / 2.d0
            end if
          end do
        end if
      end do
    end do
!   Initialize velocities by randomizing the direction
    call random_number(vel)
    vel = 2.d0*vel-1.d0
    if( planar )then
      vel(3, 1:Np) = 0.d0
    end if
    do i = 1, Np
      v_mod = dsqrt(dot_product(vel(1:3, i), vel(1:3, i)))
      vel(1:3, i) = vel(1:3, i) / v_mod
    end do
    do i = 1,3
      vel(i,:) = dsqrt(3.d0*kB*T / M(:) * dfloat(Np-1)/dfloat(Np) ) * vel(i,:)
    end do
    call remove_cm_vel(vel, M)
    write(*,'(A)') ']'
    write(*,*)
  end if

! Save initial configuration for statistics
  pos0 = pos

! Open log and trajectory files
  open(unit=10, file="log.dat", status="unknown")
  open(unit=20, file="trj.xyz", status="unknown")


  write(*,'(1X,A)') "Running MD simulation"
  write(*,*)
  write(*,'(1X,A,1X,F8.1,1X,A)') "T  =", T, "K"
  write(*,'(1X,A,1X,F8.4,1X,A)') "L  =", L(1), "Angst."
  write(*,'(1X,A,1X,I8,1X,A)') "Np =", Np, "particles"
  write(*,*)
  write(*,'(1X,A)') "Progress:"
  write(*,*)
  write(*,'(1X,A)') ">...................................<"
  write(*,'(1X,A)',advance='no') '['
  update_bar = dfloat(Nt)/36.d0
  i_ub = 1


! Do the molecular dynamics
  do step = 0, Nt
!   Update progress bar every Nt/36 iterations
    if(dfloat(step) > dfloat(i_ub)*update_bar)then
      write(*,'(A)', advance='no') '='
      i_ub = i_ub + 1
    end if
    f = 0.d0
    Ep = 0.d0
    Ek = 0.d0
    virial = 0.d0
!   Rebuild list every now and then
    if( int(step/nbuild)*nbuild == step )then
      call build_neighbors(pos, Np, L, PBC, R_neighbors, nn, list)
    end if
    do i = 1, Np
      do k = 1, nn(i)
        j = list(k,i)
        if( j > i )then
          call lj_potential(pos(1:3, i), pos(1:3, j), sigma(which_species(i)), sigma(which_species(j)), &
                            epsilon(which_species(i)), epsilon(which_species(j)), &
                            0.5d0*(Rcut(which_species(i))+Rcut(which_species(j))), L, PBC, Epot, fi)
          f(1:3, i) = f(1:3, i) + fi(1:3)
          f(1:3, j) = f(1:3, j) - fi(1:3)
          Ep = Ep + Epot
!         Compute the virial
          call get_distance( pos(1:3, i), pos(1:3, j), L, PBC, dist, d )
          virial = virial - dot_product(fi,dist)
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
    end do

!   Get pressure from the virial
    V = L(1)*L(2)*L(3)
    P = conv_to_bar* ( kB*dfloat(Np)*T/V + virial/V/3.d0 )
!   Scale box if necessary
    if( rescale )then
      call berendsen_barostat(L, P0, P, tau_p, dt)
      do i = 1, Np
        call berendsen_barostat(xi_prev(1:3,i), P0, P, tau_p, dt)
      end do
    end if
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
          write(20,*) 'Lattice="', L(1), 0.d0, 0.d0, 0.d0, L(2), 0.d0, 0.d0, 0.d0, L(3), '" ', &
                      'Properties=species:S:1:pos:R:3:vel:R:3'
        end if
        if( wrap )then
          call get_distance(0.5d0*L, pos(1:3, i), L, PBC, new_pos, d)
          write(20,*) trim(species(which_species(i))), 0.5d0*L+new_pos(1:3), vel(1:3, i)
        else
          write(20,*) trim(species(which_species(i))), pos(1:3, i), vel(1:3, i)
        end if
      end if
    end do


!   Write time, potential energy, kinetic energy, instantaneous temperature, pressure,
!   and order parameters. The pressure only makes sense if we're rescaling the box
    if( int(step/write_log)*write_log == step )then
!     Compute mean squared displacement and its std
      msd = sum( (pos-pos0)**2 )/dfloat(Np)
      msd_std = dsqrt( sum( ((pos-pos0)**2 - msd)**2 / dfloat(Np) ) )
      if( l_max >= 0 )then
        do i = 0, l_max
          call get_Ql(pos, L, PBC, nn, list, rdf_min, i, Ql(i))
        end do
        write(10,*) dfloat(step)*dt, Ep, Ek, 2.d0/3.d0/dfloat(Np-1)/kB*Ek, P, &
                    msd, msd_std, Ql(0:l_max)
      else
        write(10,*) dfloat(step)*dt, Ep, Ek, 2.d0/3.d0/dfloat(Np-1)/kB*Ek, P, &
                    msd, msd_std
      end if
    end if
  end do

  close(10)
  close(20)

  write(*,'(A)') ']'

  write(*,*)
  write(*,'(1X,A)') "*************************************"
  write(*,*)


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
