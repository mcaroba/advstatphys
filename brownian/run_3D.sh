mkdir -p 3D/
#for T in $(seq 150 50 400); do
for T in 175 225 275 325 375; do

echo "Doing T = $T K"

# Do the equilibration
cat>input<<eof
Np = 100
L = 25.
T = $T
Nt = 50000
dt = 2.
tau = 100.
restart = .false.
write_xyz = 50000
write_log = 50000
planar = .false.
Nspecies = 1
formula = 100
eof

./lj.ex > /dev/null




# Now do production
cat>input<<eof
Np = 100
L = 25.
T = $T
Nt = 50000
dt = 2.
tau = 100.
restart = .true.
write_xyz = 50000
write_log = 100
planar = .false.
Nspecies = 1
formula = 100
eof

./lj.ex > /dev/null

mv log.dat 3D/log_$T.dat


done
