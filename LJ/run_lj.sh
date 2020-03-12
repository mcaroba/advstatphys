n=100
for T in $(seq 80 10 130); do
for i in $(seq 1 1 10); do
echo "Doing $n $T $i"
L=$(echo "print (float($n)/float($i)/0.0001)**(1./3.)" | python)
echo "$n $L $T" | ./lj.ex
mv log.dat log_${n}_${i}_${T}.dat
done
echo ""
done
