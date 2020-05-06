# This code builds the lj.f90 code including the external module dependencies

# Create temp dir
mkdir -p temp
cd temp

####################################################
# Compile modules
for module in potentials neighbors integrators; do
gfortran -O3 -c ../../md/$module.f90
done
gfortran -O3 -c ../../analyze/analyze.f90
####################################################

# Compile program
gfortran -O3 -o ../lj.ex ../lj.f90 *.o

# Clean up
cd ..
rm -rf temp
