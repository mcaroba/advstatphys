module analyze

! Copyright (c) 2020 Miguel A. Caro

  use potentials

  contains


!**************************************************************************
  subroutine get_Qlm(pos, Lv, PBC, nn, nn_list, rmax, l, m, Qlm)

    implicit none

!   Input variables
    real*8, intent(in) :: pos(:,:), Lv(:), rmax
    integer, intent(in) :: nn(:), nn_list(:,:), l, m
    logical, intent(in) :: PBC(:)
!   Output variables
    complex*16, intent(out) :: Qlm
!   Internal variables
    real*8 :: fc, rmin, dist(1:3), d, theta, phi, pi, Ncoord
    integer :: natoms, i, j, k
    complex*16 :: sum

    natoms = size(nn,1)

    pi = dacos(-1.d0)

!   Hardcode the soft cutoff
    rmin = 0.9d0*rmax

    Qlm = 0.d0
    do i = 1, natoms
      Ncoord = 0.d0
      sum = 0.d0
      do k = 1, nn(i) 
        j = nn_list(k, i)
        call get_distance(pos(1:3,i), pos(1:3,j), Lv, PBC, dist, d)
        if( d < rmin )then
          theta = dacos( dist(3) / d )
          phi = datan2( dist(2), dist(1) )
          sum = sum + ylm_double(l, m, theta, phi)
          Ncoord = Ncoord + 1.d0
        else if( d < rmax )then
          theta = dacos( dist(3) / d )
          phi = datan2( dist(2), dist(1) )
          fc = 0.5d0 + 0.5d0*dcos((d-rmin)/(rmax-rmin)*pi)
          sum = sum + fc * ylm_double(l, m, theta, phi)
          Ncoord = Ncoord + fc
        end if
      end do
      if( Ncoord > 0.d0 )then
!        Qlm = Qlm + sum/Ncoord
        Qlm = Qlm + sum/6.d0
      end if
    end do
    Qlm = Qlm/dfloat(natoms)


  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine get_Ql(pos, Lv, PBC, nn, nn_list, rmax, l, Ql)

    implicit none

!   Input variables
    real*8, intent(in) :: pos(:,:), Lv(:), rmax
    integer, intent(in) :: nn(:), nn_list(:,:), l
    logical, intent(in) :: PBC(:)
!   Output variables
    real*8, intent(out) :: Ql
!   Internal variables
    real*8 :: pi
    complex*16 :: Qlm
    integer :: m

    pi = dacos(-1.d0)

    Ql = 0.d0
    do m = -l, l
      call get_Qlm(pos, Lv, PBC, nn, nn_list, rmax, l, m, Qlm)
      Ql = Ql + realpart(Qlm)**2 + imagpart(Qlm)**2
    end do
    Ql = dsqrt(4.d0*pi/dfloat(2*l+1) * Ql)

  end subroutine
!**************************************************************************




!**************************************************************************
  function plm_double(l, m, x)
!   Associated Legendre polynomial Plm(x).
!   l is > 0, m can take values from 0 to l and x is within the [-1, 1] domain
    implicit none
    integer :: l, m, i, l2
    real*8 :: plm_double, x, dfact_arg, pl2m, pmm, pmp1m, sq
    if(l < 0 .or. m < 0 .or. m > l .or. dabs(x) > 1.d0)then
        write(*,*) "Bad arguments for associated Legendre polynomial"
    end if
    pmm = 1.d0
    if(m > 0) then
      sq = dsqrt( (1.d0 - x**2) )
      dfact_arg = -1.d0
      do i = 1, m
        dfact_arg = dfact_arg + 2.d0
        pmm =  -pmm * dfact_arg * sq
      end do
    end if
    if(l == m)then
      plm_double = pmm
    else
      pmp1m = x * (2.d0*dfloat(m) + 1.d0) * pmm
      if(l == m+1) then
        plm_double = pmp1m
      else
        do l2 = m+2, l
          pl2m = ( x*dfloat(2*l2-1)*pmp1m - dfloat(l2+m-1)*pmm ) / dfloat(l2-m)
          pmm = pmp1m
          pmp1m = pl2m
        end do
        plm_double = pl2m
      end if
    end if
  end function
!**************************************************************************




!**************************************************************************
  function ylm_double(l, m, theta, phi)
    implicit none
    real*8 :: plm, theta, phi, pref, pi, ylm_r, ylm_i, fact1, fact2
    complex*16 :: ylm_double
    integer :: l, m, i, sgn
    pi = dacos(-1.d0)
!   Factorials for |m|
    fact1 = 1.d0
    do i = 1, l-abs(m)
      fact1 = fact1 * dfloat(i)
    end do
    fact2 = 1.d0
    do i = 1, l+abs(m)
      fact2 = fact2 * dfloat(i)
    end do
!   Sign
    sgn = 1
    do i = 1, m
      sgn = -1 * sgn
    end do
!   Prefactor
    pref = dsqrt( dfloat(2*l+1)/(4.d0*pi) * fact1 / fact2 )
!   This is pl|m| 
    plm = plm_double(l, abs(m), dcos(theta))
!   Real part
    ylm_r = pref * plm * dcos( dfloat(abs(m))*phi )
!   Imaginary part
    ylm_i = pref * plm * dsin( dfloat(abs(m))*phi )
!   For m >= 0
    if( m >= 0 )then
      ylm_double = ylm_r * (1.d0, 0.d0) + ylm_i * (0.d0, 1.d0)
!   For m < 0
    else
      ylm_double = dfloat(sgn) * ( ylm_r * (1.d0, 0.d0) - ylm_i * (0.d0, 1.d0) )
    end if
  end function
!**************************************************************************

end module
