Np=32
L=20.
Nt=50000
dt=2.

T_list="10 20 50 70 100 150 200 350 500 1000"
P_list="0 1 2 5 10 20 50 100 200 500 1000 2000 5000"

for T in $T_list; do
for P in $P_list; do

echo "Doing T = $T; P = $P"

# For the initial structure in each T series we initialize
if [ "$P" == 0 ]; then
restart=".false."
else
restart=".true."
fi

# This is to equilibrate
cat>input<<eof
Np = $Np
L = $L
T = $T
P = $P
Nt = $Nt
dt = $dt
tau = 100.
tau_p = 100.
restart = $restart
rescale = .true.
eof

./lj.ex > /dev/null


# This is production
cat>input<<eof
Np = $Np
L = $L
T = $T
P = $P
Nt = $Nt
dt = $dt
tau = 100.
tau_p = 100.
l_max = 6
restart = .true.
rescale = .true.
eof

./lj.ex > /dev/null

mkdir -p results_$Np
cp log.dat results_$Np/log_${T}_${P}.dat
cp trj.xyz results_$Np/trj_${T}_${P}.xyz

done
done



# Exctract the results
rm results_$Np/Qs.dat
for T in $T_list; do
for P in $P_list; do
awk -v T0=$T -v P0=$P '{T += $4; P += $5; Q0 += $6; Q2 += $8; Q4 += $10; Q6 += $12; n += 1} \
END {print T0, P0, T/n, P/n, Q0/n, Q2/n, Q4/n, Q6/n}' results_$Np/log_${T}_${P}.dat >> results_$Np/Qs.dat
done
echo "" >> results_$Np/Qs.dat
done
