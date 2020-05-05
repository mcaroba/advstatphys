# This code builds the lj.f90 code including the external module dependencies

mkdir -p temp
cd temp
for module in potentials neighbors integrators; do
#if [ ! -f "$module.o" ]; then
gfortran -O3 -c ../../md/$module.f90
#fi
done

gfortran -O3 -o ../lj.ex ../lj.f90 *.o

cd ..
rm -rf temp
