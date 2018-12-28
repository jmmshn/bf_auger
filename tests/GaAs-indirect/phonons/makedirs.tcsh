#!/bin/tcsh
alias calc 'awk "BEGIN{ print \!* }" '

set QLIST = kirr_cartesian.dat
set NQ = `wc $QLIST | awk '{print $1}'`
echo "found $NQ phonons"

#-- loop over the phonon q-vectors
set IQ = 0
while ( $IQ < $NQ )
    set IQ = `calc $IQ + 1`
    set LINE = `sed -ne "${IQ}p" $QLIST`
    set QX   = $LINE[1]
    set QY   = $LINE[2]
    set QZ   = $LINE[3]
    #-- Make subdirectories
    echo "phonon $IQ"
    mkdir -pv phonon_$IQ
    cd phonon_$IQ
    sed -e "s/QX/$QX/g" -e "s/QY/$QY/g" -e "s/QZ/$QZ/g" ../ph1.in > ph1.in
    sed -e "s/QX/$QX/g" -e "s/QY/$QY/g" -e "s/QZ/$QZ/g" ../ph2.in > ph2.in
    cp ../nscf.in .
    cp -r ../GaAs.save .
    #sbatch ../submit.sh
    cd ..
end
