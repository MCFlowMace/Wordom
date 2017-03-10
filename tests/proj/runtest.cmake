##---------------------------------------------------------------------------
## Author:      Simone Conti
## Copyright:   (C) 2014 Simone Conti 
## License:     GPLv3+
##---------------------------------------------------------------------------

include("${SRCDIR}/../wrdmacro.cmake")

# Set the command line
set(TEST_ARGS -imol ${SRCDIR}/../loop1.pdb -itrj ${SRCDIR}/../loop1.dcd -iA ${SRCDIR}/proj.inp -otxt proj.out)

# Copy input files called inside the wordom inp file
execute_process( COMMAND ${CMAKE_COMMAND} -E copy ${SRCDIR}/eigvec.txt eigvec.txt)

# Run wordom
wrd_run(${WORDOM} "proj" "${TEST_ARGS}")

# Compare output files
wrd_compare_file("proj.out")
wrd_compare_file("proj.dcd")

