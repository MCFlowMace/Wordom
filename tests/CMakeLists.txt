##---------------------------------------------------------------------------
## Author:      Simone Conti
## Copyright:   (C) 2014 Simone Conti 
## License:     GPLv3+
##---------------------------------------------------------------------------

include("${PROJECT_SOURCE_DIR}/tests/wrdmacro.cmake")

# Distances (DISTANCE), contacts (CONSTACTS), angles (ANGLE), Dihedral angles (DIHEDRAL), 
#   within (WITHIN), Radius of gyration (RGYR) and Hydrogen bond detection (HBOND)
wrd_add_test("geometry")

# RMSD, DRMS, RMSF
wrd_add_test("moldiff")

# Principal Component Analysis (PCA)
wrd_add_test("pca")

# Principal Component Analysis projection (PROJECT)
wrd_add_test("proj")

# RMSD- and DRMS-based Clustering (CLUSTER)
wrd_add_test("cluster1")
wrd_add_test("cluster2")

# RMSD- and DRMS-based Clustering on GPU (GCLUSTER)
if(CUDA_FOUND)
	wrd_add_test("gCluster")
endif(CUDA_FOUND)

# Molecular Surface (SURF)
wrd_add_test("surface")

# Protein Structure Network (PSN)
wrd_add_test("psn")

# Secondary Structure Assignment (SSA) 
wrd_add_test("ssa")

# Cross-Correlation (CORR)
wrd_add_test("corr")

# Tilt (TILT) and twist (TWIST) angles
wrd_add_test("tilt_twist")

# HOLE module
wrd_add_test("hole1")
wrd_add_test("hole2")
wrd_add_test("hole3")

# COM module
wrd_add_test("com")

# FLUX module
wrd_add_test("flux")
