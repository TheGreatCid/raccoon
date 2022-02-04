#!/bin/sh

for i in  1e3 1e4 1e5
do 
	sed -i "22s/.*/sr=$i/" ../genprojects/kalthoff/HomCube/mechanicaltherm_Damage.i
       mpirun -n 16 ./raccoon-opt -i ../genprojects/kalthoff/HomCube/mechanicaltherm_Damage.i
done

