Np=32
tau=20.
points=5001

T_list="10 20 50 70 100 150 200 350 500 1000"
P_list="0 1 2 5 10 20 50 100 200 500 1000 2000 5000"

rm -f entropies.dat

for T in $T_list; do
for P in $P_list; do

echo "Doing T = $T; P = $P"

python groups.py $Np > groups
echo "1-$Np" > supergroups

ln -sf ../results_cv_${Np}/trj_${T}_${P}.xyz traj.xyz
cell=$(awk 'NR==2{print $2/10.; exit}' traj.xyz)

cat>input<<eof
points = $points
tau = $tau
cell = $cell $cell $cell
temperature = $T
format = xyz
estimate_velocities = .false.
renormalize_dos = .false.
eof

DoSPT > /dev/null

S=$(tail -1 entropy | awk '{print $6/32.}')
dof=$(tail -1 fluidicity | awk '{print $2}')

echo "$T $P $S $dof" >> entropies.dat

done
echo "" >> entropies.dat
done
