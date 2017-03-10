##---------------------------------------------------------------------------
## Author:      Simone Conti
## Copyright:   (C) 2014 Simone Conti 
## License:     GPLv3+
##---------------------------------------------------------------------------

include("${SRCDIR}/../wrdmacro.cmake")

# Set the command line
set(TEST_ARGS -imol ${SRCDIR}/../loop1.pdb -itrj ${SRCDIR}/../loop1.dcd -iA ${SRCDIR}/cluster2.inp -otxt cluster2.out)

# Copy input files called inside the wordom inp file
execute_process( COMMAND ${CMAKE_COMMAND} -E copy ${SRCDIR}/matrix_cluster.dat matrix_cluster.dat)
execute_process( COMMAND ${CMAKE_COMMAND} -E copy ${SRCDIR}/c_qt c_qt)

# Run wordom
wrd_run(${WORDOM} "proj" "${TEST_ARGS}")

# Compare output files
wrd_compare_file("cluster2.out")
wrd_compare_file("c_hiero")

