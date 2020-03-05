#!/bin/bash -l

Ts="84.791 88.336 92.303 95.058 101.398 102.014 \
    105.513 108.146 113.318 117.501 123.990"

Bs="249.34 229.89 211.79 200.87 178.73 177.65 \
    166.06 160.27 149.58 140.58 127.99"

dr=0.01

n=$(echo $Ts | wc -w)

for epsilon in 0.0115 0.0116 0.0117 0.0118 0.0119; do
for sigma in 3.21 3.22 3.23 3.24 3.25; do
#epsilon=0.0117; sigma=3.235
error=0.
for i in $(seq 1 1 $n); do
T=$(echo $Ts | awk -v col=$i '{print $col}')
B=$(echo $Bs | awk -v col=$i '{print $col}')
Bc=$(echo "$epsilon $sigma $T $dr" | ./LJ_B2.ex | awk 'NR==2{print -$2}')
error=$(echo $error $B $Bc | awk '{print $1+($2-$3)**2}')
#echo $T $B $Bc
done
echo $epsilon $sigma $error
done
done

