##---------------------------------------------------------------------------
## Author:      Simone Conti
## Copyright:   (C) 2014 Simone Conti 
## License:     GPLv3+
##---------------------------------------------------------------------------

include("${SRCDIR}/../wrdmacro.cmake")

# Set the command line
set(TEST_ARGS -imol ${SRCDIR}/../loop1.pdb -itrj ${SRCDIR}/../loop1.dcd -iA ${SRCDIR}/cluster1.inp -otxt cluster1.out)

# Run wordom
wrd_run(${WORDOM} "cluster1" "${TEST_ARGS}")

# Compare output files
wrd_compare_file("cluster1.out")
wrd_compare_file("c_leader")
wrd_compare_file("c_qt")
wrd_compare_file("matrix_cluster.dat")

