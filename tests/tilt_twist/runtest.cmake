##---------------------------------------------------------------------------
## Author:      Simone Conti
## Copyright:   (C) 2014 Simone Conti 
## License:     GPLv3+
##---------------------------------------------------------------------------

include("${SRCDIR}/../wrdmacro.cmake")

# Set the command line
set(TEST_ARGS -imol ${SRCDIR}/../pentamer.pdb -itrj ${SRCDIR}/../pentamer.dcd -iA ${SRCDIR}/tilt_twist.inp -otxt tilt_twist.out)

# Run wordom
wrd_run(${WORDOM} "tilt_twist" "${TEST_ARGS}")

# Compare output files
wrd_compare_file("tilt_twist.out")

