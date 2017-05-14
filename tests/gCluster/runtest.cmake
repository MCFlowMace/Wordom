##---------------------------------------------------------------------------
## Author:      Simone Conti
## Copyright:   (C) 2014 Simone Conti 
## License:     GPLv3+
##---------------------------------------------------------------------------

include("${SRCDIR}/../wrdmacro.cmake")

# Set the command line
set(TEST_ARGS -imol ${SRCDIR}/../loop1.pdb -itrj ${SRCDIR}/../loop1.dcd -iA ${SRCDIR}/gCluster1.inp)

# Run wordom
wrd_run(${WORDOM} "gCluster1" "${TEST_ARGS}")

# Compare output files
wrd_compare_file("g_RMSD.out")
wrd_compare_file("g_RMSD")

# Set the command line
set(TEST_ARGS2 -imol ${SRCDIR}/../loop1.pdb -itrj ${SRCDIR}/../loop1.dcd -iA ${SRCDIR}/gCluster2.inp)

# Run wordom
wrd_run(${WORDOM} "gCluster2" "${TEST_ARGS2}")

# Compare output files
wrd_compare_file("g_DRMS.out")
wrd_compare_file("g_DRMS")

