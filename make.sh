#/bin/bash

#code="ptmc_SiO4_TAA_nanoparticles.f90"
#module2="variables_SiO4_TAA_nanoparticles.f90"
#module1="input_SiO4_TAA_nanoparticles.f90"

code="ptmc_SiO4_SDA.f90"
module2="variables_SiO4_SDA.f90"
module1="input_SiO4_SDA.f90"

output="xcutable"
flag0="-O5"
flag1="-fbounds-check -pg -C"

compiler1="mpif90"
compiler2="/opt/intel/bin/ifort"

#---------------------------------------
#rm $output
#compile the program
$compiler1 $flag0 $module2 $module1 $code -o $output

