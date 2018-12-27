#!/bin/bash

filename="klist.elec.irr"
while read -r line
do
   echo "$line" > tmp
   k_id=$(awk '{print $1}' tmp)
   kx=$(awk '{print $2}' tmp)
   ky=$(awk '{print $3}' tmp)
   kz=$(awk '{print $4}' tmp)
   rm tmp
   mkdir k_$k_id
   cd k_$k_id
   mkdir GaAs.save
   cp ../charge-density.dat GaAs.save
   cp ../data-file.xml GaAs.save
   cp ../nscf.in .
   sed -i 's/kx/'"$kx"'/g' nscf.in
   sed -i 's/ky/'"$ky"'/g' nscf.in
   sed -i 's/kz/'"$kz"'/g' nscf.in
   cp ../pw2wannier.in .
   cp ../GaAs.nnkp .
   sed -i 's/kx/'"$kx"'/g' GaAs.nnkp
   sed -i 's/ky/'"$ky"'/g' GaAs.nnkp
   sed -i 's/kz/'"$kz"'/g' GaAs.nnkp
   cd ..
done < "$filename"
