##---------------------------------------------------------------------------
## Author:      Simone Conti
## Copyright:   (C) 2014 Simone Conti 
## License:     GPLv3+
##---------------------------------------------------------------------------

include("${SRCDIR}/../wrdmacro.cmake")

# Set the command line
set(TEST_ARGS -imol ${SRCDIR}/../pentamer2.pdb -itrj ${SRCDIR}/../pentamer2.dcd -iA ${SRCDIR}/flux.inp -otxt flux.out)

# Run wordom
wrd_run(${WORDOM} "flux" "${TEST_ARGS}")

# Compare output files
wrd_compare_file("flux.out")
wrd_compare_file("pent-trans.dat")

