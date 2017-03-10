##---------------------------------------------------------------------------
## Author:      Simone Conti
## Copyright:   (C) 2015 Simone Conti 
## License:     GPLv3+
##---------------------------------------------------------------------------

include("${SRCDIR}/../wrdmacro.cmake")

# Set the command line
set(TEST_ARGS -imol ${SRCDIR}/../4COF.pdb -itrj ${SRCDIR}/../4COF.pdb -ia HOLE --POINT \(3.5:-1.1:139.9\) --VECT  \(0.34:0.01:0.94\) --SELE "*" --TITLE 4COF --SEED 1989 )

# Run wordom
wrd_run(${WORDOM} "hole1" "${TEST_ARGS}")

# Compare output files
wrd_compare_file("4COF-prof.dat")
wrd_compare_file("4COF-traj.vtf")

