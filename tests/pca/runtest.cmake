##---------------------------------------------------------------------------
## Author:      Simone Conti
## Copyright:   (C) 2014 Simone Conti 
## License:     GPLv3+
##---------------------------------------------------------------------------

include("${SRCDIR}/../wrdmacro.cmake")

# Set the command line
set(TEST_ARGS -imol ${SRCDIR}/../loop1.pdb -itrj ${SRCDIR}/../loop1.dcd -iA ${SRCDIR}/pca.inp -otxt pca.out)

# Run wordom
wrd_run(${WORDOM} "pca" "${TEST_ARGS}")

# Compare output files
wrd_compare_file("pca.out")
wrd_compare_file("pca-eigval.txt")
wrd_compare_file("pca-eigvec.txt")
wrd_compare_file("pca-matrix.txt")

