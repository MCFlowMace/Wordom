##---------------------------------------------------------------------------
## Author:      Simone Conti
## Copyright:   (C) 2014 Simone Conti 
## License:     GPLv3+
##---------------------------------------------------------------------------

include("${SRCDIR}/../wrdmacro.cmake")

# Set the command line
set(TEST_ARGS -imol ${SRCDIR}/looprad.pdb -itrj ${SRCDIR}/../loop1.dcd -iA ${SRCDIR}/surface.inp -otxt surface.out -end 5)

# Run wordom
wrd_run(${WORDOM} "surface" "${TEST_ARGS}")

# Compare output files
wrd_compare_file("surface.out")

