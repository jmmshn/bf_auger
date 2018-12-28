#!/bin/bash

for((i=1;i<=6;i++))
do
cd phonon_$i
grep freq GaAs.dyn1 | cut -b 45-57 > omega
cd ..
done
