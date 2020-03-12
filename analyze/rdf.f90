! Copyright by Miguel Caro (2018)

program rdf

  use potentials

  implicit none

  character*256 :: filename, junk
  real*8 :: L(1:3), dr, rmax, dist(1:3), d, pi
  real*8, allocatable :: pos(:,:), g(:)
  integer :: natoms, i, iostatus = 0, nsnaps, j, k
  logical :: PBC(1:3)

  pi = dacos(-1.d0)

  read(*,*) filename, L(1:3)

  rmax = 1.d10
  do i = 1, 3
    if( L(i)/2.d0 < rmax )then
      rmax = L(i)/2.d0
    end if
  end do

  PBC = .true.

  dr = 0.1d0
  allocate( g(0:int(rmax/dr)-1) )
  g = 0.d0

  open(unit=10, file=trim(adjustl(filename)), status="old")

  read(10, *, iostat=iostatus) natoms
  allocate( pos(1:natoms, 1:3) )

  nsnaps = 0
  do while( iostatus == 0 )
    nsnaps = nsnaps + 1

    read(10, *)
    do i = 1, natoms
      read(10, *, iostat=iostatus) junk, pos(i, 1:3)
    end do
    do i = 1, natoms
      do j = i+1, natoms
        call get_distance(pos(i,1:3), pos(j,1:3), L, PBC, dist, d)
        if( d < dr*int(rmax/dr) )then
          k = int(d/dr)
          g(k) = g(k) + 2.d0/(dfloat(k)+0.5d0)**2/dr**3
        end if
      end do
    end do
    read(10, *, iostat=iostatus) natoms
  end do
  close(10)

  g = g/4.d0/pi/dfloat(natoms)**2*L(1)*L(2)*L(3)/dfloat(nsnaps)

  open(unit=10, file="rdf.dat", status="unknown")
  do i = 0, int(rmax/dr)-1
    write(10,*) dfloat(i)*dr, g(i)
  end do
  close(10)

end program
