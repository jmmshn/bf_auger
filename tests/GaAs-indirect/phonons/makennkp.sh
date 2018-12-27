#!/bin/bash
FILE='GaAs.nnkp'
echo "writing $FILE"

cat > $FILE <<EOF
calc_only_A  :  F

begin real_lattice
-2.805452740     0.000000000     2.805452740
 0.000000000     2.805452740     2.805452740
-2.805452740     2.805452740     0.000000000 
end real_lattice

begin recip_lattice
-1.119816638    -1.119816638    1.119816638
 1.119816638     1.119816638    1.119816638
-1.119816638     1.119816638   -1.119816638
end recip_lattice

begin kpoints
EOF

NQNSCF=`grep 'number of k points' ph2.out | tail -1 | awk '{print $5}'`
echo $NQNSCF >> $FILE
grep 'k(' ph2.out | tail -$NQNSCF | cut -b 22-56 >> $FILE

cat >> $FILE <<EOF
end kpoints

begin projections
   0
end projections

begin nnkpts
   1
EOF

J=0
while [[ $J -lt $NQNSCF ]] ; do
    ((J++))
    echo $J $J " 0 0 0 " >> $FILE
done

cat >> $FILE <<EOF
end nnkpts

begin exclude_bands
   0
end exclude_bands
EOF
