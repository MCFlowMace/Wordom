##---------------------------------------------------------------------------
## Author:      Simone Conti
## Copyright:   (C) 2015 Simone Conti 
## License:     GPLv3+
##---------------------------------------------------------------------------

include("${SRCDIR}/../wrdmacro.cmake")

# Set the command line
set(TEST_ARGS -imol ${SRCDIR}/../2VL0.pdb -itrj ${SRCDIR}/../2VL0.pdb -iA ${SRCDIR}/hole.inp )

# Run wordom
wrd_run(${WORDOM} "hole2" "${TEST_ARGS}")

# Compare output files
wrd_compare_file("2VL0a-prof.dat")
wrd_compare_file("2VL0a-traj.vtf")

