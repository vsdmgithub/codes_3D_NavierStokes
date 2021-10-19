#!/bin/sh
sed -i 's/ECCODE/4.0D0/g' NSE_code.f90 
sed -i 's/run_EC/run_V4/g' modules-primary/system_initialcondition.f90 
sed -i 's/run_EC/run_V4/g' modules-primary/system_basicoutput.f90 
sed -i 's/ex_ECCODE/ex_V4CODE/g' makefile
make
sed -i 's/ex_V4CODE/ex_ECCODE/g' makefile
sed -i 's/run_V4/run_EC/g' modules-primary/system_initialcondition.f90 
sed -i 's/run_V4/run_EC/g' modules-primary/system_basicoutput.f90 
sed -i 's/4.0D0/ECCODE/g' NSE_code.f90 
cp run.sub run_V4.sub
sed -i 's/ECCODE/V4CODE/g' run_V4.sub 
sbatch run_V4.sub









