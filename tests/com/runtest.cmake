##---------------------------------------------------------------------------
## Author:      Simone Conti
## Copyright:   (C) 2014 Simone Conti 
## License:     GPLv3+
##---------------------------------------------------------------------------

include("${SRCDIR}/../wrdmacro.cmake")

# Set the command line
set(TEST_ARGS -imol ${SRCDIR}/../loop1.pdb -itrj ${SRCDIR}/../loop1.dcd -iA ${SRCDIR}/com.inp -otxt com.out)

# Run wordom
wrd_run(${WORDOM} "com" "${TEST_ARGS}")

# Compare output files
wrd_compare_file("com.out")

