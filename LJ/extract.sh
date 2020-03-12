n=100
for T0 in $(seq 80 10 130); do
for i in $(seq 1 1 10); do
awk -v n=$n -v i=$i -v T0=$T0 '{if($1 > 10000){count += 1; p += $5; t += $4; \
v += $2}} END {print n, (n/i/0.0001)**(1./3.), T0, p/count, t/count, v/count}' \
log_${n}_${i}_${T}.dat
done
echo ""; echo ""
done
