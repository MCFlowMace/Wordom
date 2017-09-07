// ------------------------------------------------------------------
// Copyright (C) 2003  University of Zurich and Univ. of Modena & Reggio Emilia
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ------------------------------------------------------------------
/*! \file tools.h
 \brief Mathematical tools and more
 
 Mathematical and other basic functions used by several files and
 functions in the wordom program. When adding a function that might be
 useful to many, place it here
*/
#ifndef TOOLS
#define TOOLS

/*!
 Diagonalization function imported from lapack.
 Diagonalizes a matrix and extracts eigenvectors and eigenvalues
 Not to be called directly, rather through DiagMatrix
*/
#ifdef LAPACK
  extern void ssyev_(char *, char *, int *, float *, int *, float *, float *, float *, int *);
#endif

// ------------------------------------------------------------------
//! Computes Root Mean Square Deviation - QCP algorithm
/*!
 RmsdCalc computes the Root Mean Square Deviation between two sets of
 coordinates, either with or without previous superposition. The
 movcoor set is modified if superposition to refcoor is required.
 Coor set are in coor[3][nato] format. The QCP algorithm is used; faster
 than Kabsch, still in testing
*/
float RmsdCalcQCP(float **refcoor, float **movcoor, int nato, int super );

#ifndef TSNE_INCLUDE
// ------------------------------------------------------------------
//! Vectorial Distance structure
/*!
  If a set of distances is to be computed with vectorialization enabled,
  this structure is used. Arrays of coordinates have to be allocated
  and aligned.
*/
typedef struct distvcts_
{
   int ndata;
   float *dist_t, *disty_t, *distz_t, *x0_t, *y0_t, *z0_t, *x1_t, *y1_t, *z1_t;
   float *dist, *disty, *distz, *x0, *y0, *z0, *x1, *y1, *z1;
} _distvcts;
// ------------------------------------------------------------------
//! General matrix structure
/*!
  general matrix structure with linear float arrangment and 
*/
typedef struct matrix_
{
   int      size;
   float  **inpmat;
   float   *matrix;
   float   *eigval;
   float  **eigvec;
} wrd_Matrix;
// ------------------------------------------------------------------
//! Matrix diagonalization structure
/*!
  If a (dense) matrix is to be diagonalized and its
  eigenvectors/eigenvalues extracted, it should be handled as 
  an _eigen structure
*/
typedef struct eigen_
{
   int      size;
   float  **inpmat;
   float   *eigval;
   float  **eigvec;
} _eigen;
// ------------------------------------------------------------------
//! Periodic Boundary Conditions structure
/*!
  If PBC are present, its parameters should be stored in a
  _pbcdata structure
*/
typedef struct pbcdata
{
   int          pbc;
   float        pbcbox[3];
   float        xsize;
   float        ysize;
   float        zsize;
} _pbcdata;
// ------------------------------------------------------------------
//! Computes Root Mean Square Deviation - better algorithm
/*!
 RmsdCalc computes the Root Mean Square Deviation between two sets of
 coordinates, either with or without previous superposition. The
 movcoor set is modified if superposition to refcoor is required.
 First attempt with QCP algorithm, Kabsch-like in case of failure
 Coor set are in coor[3][nato] format
*/
float RmsdCalc( float **refcoor,        /*!< reference coordinate set */
                float **movcoor,        /*!< moving coordinate set */
                int n_ato,              /*!< number of atoms for both sets */
                int super_flag          /*!< whether superposition is required */
              );
// ------------------------------------------------------------------
//! Computes Root Mean Square Deviation - preserves original data
/*!
 RmsdCalc_preserve computes the Root Mean Square Deviation between two
 sets of coordinates, either with or without previous superposition.
 By using local arrays before calling RmsdCalc it preserves the original
 data. Coor set are in coor[nato][3] format.
*/
float RmsdCalc_preserve(float **refcoor, float **movcoor, int nato, int super );
// ------------------------------------------------------------------
//! Computes Root Mean Square Deviation - no superposition
/*!
 RmsdCalc_nosup computes the Root Mean Square Deviation between two
 sets of coordinates, without previous superposition. Coor sets are in 
 coor[3][nato] format
*/
float RmsdCalc_nosup(float **refcoor, float **movcoor, int nato);
// ------------------------------------------------------------------
//! Computes Root Mean Square Deviation - Kabsch algorithm
/*!
 RmsdCalcKabsch computes the Root Mean Square Deviation between two sets
 of coordinates, either with or without previous superposition. The
 movcoor set is modified if superposition to refcoor is required.
 Coor set are in coor[nato][3] format. The "traditional" Kabsch
 is used.
*/
float RmsdCalcKabsch( float **refcoor,        /*!< reference coordinate set */
                      float **movcoor,        /*!< moving coordinate set */
                      int n_ato,              /*!< number of atoms for both sets */
                      int super_flag          /*!< whether superposition is required */
                    );
// ------------------------------------------------------------------
//! Computes Root Mean Square Deviation - Kabsch algorithm
/*!
 RmsdCalcKabsch3n computes the Root Mean Square Deviation between two
 sets of coordinates, either with or without previous superposition. The
 movcoor set is modified if superposition to refcoor is required.
 Coor set are in coor[3][nato] format. The "traditional" Kabsch
 is used.
*/
float RmsdCalcKabsch3n( float **refcoor,        /*!< reference coordinate set */
                        float **movcoor,        /*!< moving coordinate set */
                        int n_ato,              /*!< number of atoms for both sets */
                        int super_flag          /*!< whether superposition is required */
                      );
// ------------------------------------------------------------------
//! Computes Root Mean Square Deviation among Molecules
/*!
 Applies RmsdCalc between the provided molecules with the provided 
 selection string
*/
float MolRmsd(Molecule *mol1, Molecule *mol2, char *selestring );
// ------------------------------------------------------------------
//! 
/*!
  
*/
// ------------------------------------------------------------------
//! 
/*!
 */
void CalcRMSDRotTrans(int npoints, float **referenceSet, float **tomoveSet, float rotMatrix[3][3], float transVec[3]);
// ------------------------------------------------------------------
//! Superposition function
/*!
 This superimposes two coordinate sets to minimize their RMSD
 according to the Kabsch algorithm. RMSD is minimized if the 
 rototraslation according to rot_matrix and vector is performed
*/
void superimpose ( int n_ato,                   /*!< number of atoms for both sets */
                   float **refcoor,             /*!< reference coordinate set */
                   float **movcoor,             /*!< moving coordinate set */
                   float rot_matrix[3][3],      /*!< rotation matrix*/
                   float vector[3]              /*!< traslation vector */
                 );
// ------------------------------------------------------------------
//! Superposition function
/*!
 This returns a new molecule, mol3, which is mol2 superimposed on mol1
 Selection must refer to the same number of atom in both input molecules
 (mol1 and mol2), or NULL is returned
*/
Molecule * molSuperimpose(  Molecule *mol1,     /*!< reference molecule */
                            Molecule *mol2,     /*!< moving molecule */
                            Selection *sele1,   /*!< selection on mol1 */
                            Selection *sele2    /*!< selection on mol2 */
                         );
// ------------------------------------------------------------------
//! Computes Distance Matrix for DRMS calculations
/*!
 
*/
void DistMtxCalc( int     pbcflag,      /*!< Periodic Boundary Conditions flag */
                  float  *pbcbox,       /*!< PBC box (x, y and z size of box) */
                  int     n_ato,        /*!< number of atoms */
                  float **coor,         /*!< coordinates' set */
                  float  *dist_mtx      /*!< internal distances matrix */
                );
// ------------------------------------------------------------------
//! Computes Distance Root Mean Square deviation (2 matrices)
/*!
 DRMS is computed between two distance matrices
*/
void Drms( int k , float  * distance , float * dist_mtx , float * drms );
//void Drms2( int k , float  * distance , float * dist_mtx , float * drms );
// ------------------------------------------------------------------
//!  Computes Distance Root Mean Square deviation (2 coor sets)
/*!
 DRMS is computed between two coordinates' sets. DistMtxCalc is used.
*/
float DrmsCalc( float **refcoor,        /*!< reference coordinate set */
                float **movcoor,        /*!< moving coordinate set */
                int     n_ato,          /*!< number of atoms for both sets */
                int     pbcflag,        /*!< Periodic Boundary Conditions flag */
                float  *pbcbox );       /*!< PBC box (x, y and z size of box) */
// ------------------------------------------------------------------
//! Matrix diagonalization + eigenvectors/eigenvalues
/*!
 
*/
void DiagMatrix ( _eigen *);
//! int comparison for qsort 
/*!
 
*/
//int compare_ints( const void* a, const void* b );
// ------------------------------------------------------------------

// ------------------------------------------------------------------
//! Computes a vector of distances
/*!
 
*/
/*void distance( int , float *, float * ,float *,float *,float *,float *,float *);*/
// ------------------------------------------------------------------
//! Dimension distvcts structure - header
/*!
 
*/
void FillDistvcts ( _distvcts * );
// ------------------------------------------------------------------
//! Compute Distance
/*!
 
*/
void vDistance ( int, float *, float *, float *, float *, float *, float *, float * );
// ------------------------------------------------------------------
//! Compute Distance with Periodic Boundary Conditions
/*!
 
*/
void vDistancePBC ( int, float *, float *, float *, float *, float *, float *, float * , float * );
// ------------------------------------------------------------------
//! Computes Distance
/*!
 
*/
//void vdistance ( _distvcts * );
// ------------------------------------------------------------------
//! Computes Distance with Periodic Boundary Conditions
/*!
 
*/
void vdistancePBC ( float *, _distvcts * );
// ------------------------------------------------------------------
// ------------------------------------------------------------------
//! More tools to compute distances with periodic boundary conditions
/*!
 
*/
float PBC1( CoorSet *coorset1, int atom1, CoorSet *coorset2, int atom2);
void DistanceAxis( float *distances, float xcoor1, float ycoor1, float zcoor1, float xcoor2, float ycoor2, float zcoor2, Pbc *pbc );
void DistanceAxes( float *distances, float xcoor1, float ycoor1, float zcoor1, float xcoor2, float ycoor2, float zcoor2, Pbc *pbc );
inline float DistanceCoor( float xcoor1, float ycoor1, float zcoor1, float xcoor2, float ycoor2, float zcoor2, Pbc *pbc );
float * DistanceVectCoor( CoorSet *coorset1, int *piSele1, CoorSet *coorset2, int *piSele2, int iNumOfAtomPairs );
float * DistanceSelCoor( CoorSet *coorset1, Selection *sele1, CoorSet *coorset2, Selection *sele2 );
//int DistanceArray( );
// ------------------------------------------------------------------
//! compute angle 
/*!
 compute angle formed by segment atom1-atom2 and segment atom2-atom3
*/
double AngleCalc(double *a, double *b, double *c);
// ------------------------------------------------------------------
//! compute dihedral angle 
/*!
 Compute the dihedral angle between four points in 3D
*/
double DiheCalc(double *a, double *b, double *c, double *d);
// ------------------------------------------------------------------
//! local version of scopy
/*!
 
*/
void  scopy( int * , float * , int   * , float * , int   * );
//! local version of saxpy
/*!
 
*/
void  saxpy( int * , float * , float * , int   * , float *, int * );
//! local version of snrm2
/*!
 
*/
// ------------------------------------------------------------------
float snrm2( int * , float * , int   * );
// ------------------------------------------------------------------
//! bubble sort 
/*!
 
*/
void DoubleSort_int( int , int *, int * );
// ------------------------------------------------------------------
//! checks if int list is sorted 
/*!
 checks if all elements in listA->nframe are sorted and returns 1 if true, else 0
*/
int checkSorted( struct intlist * listA );
//----------------------------------------------------------
// --- String Functions ---
// ------------------------------------------------------------------
//! Expands a selection string
/*!
 parses a selection string and adds relevant "/"s and "*"s until 
 the /chain/segment/resid/atom structure sequence is 
 complete. Missing sections are added to the left (ie a single 
 section is taken as atom identifier) to be compatible with the
 old notation /segment/resid/atom
*/
char * parseSlash( char *datastring );

int wrd_isnumber( char *nstring );

//----------------------------------------------------------
// --- Math Functions ---
// ------------------------------------------------------------------
//! Computes the scalar product (a,b)
/*!
 
*/
float dot_prod(float * , float * ) ;
//! Computes the scalar product (a,b)
/*!
 
*/
double wrd_vect_dotprod(double *a, double *b, size_t n);
// ------------------------------------------------------------------
//! Computes the vector product res = a x b
/*!
 
*/
void cross_prod(float * , float * , float * ) ;
//! Computes the vector product res = a x b
/*!
 
*/
void wrd_vect_crossprod(double *a, double *b, double *res) ;
// ------------------------------------------------------------------
//! Computes the triple product -> a.(b X c)
/*!
 
*/
//float triple_prod(float * , float * , float * ) ;
double wrd_vect_tripleprod(double *a, double *b, double *c);
// ------------------------------------------------------------------
//! Normalizes a vector (vec)
/*!
 
*/
void w_norm(float * ) ;
//! Normalizes a vector (vec)
/*!
 
*/
void wrd_vect_norm(double *vec, size_t n);
// ------------------------------------------------------------------
//! Determines the sign of a float variable
/*!
 
*/
float wsign(float ) ;
// ------------------------------------------------------------------
//! 
/*!
 
*/
float *DistanceVectCoor( CoorSet *coorset1, int *piSele1, CoorSet *coorset2, int *piSele2, int iNumOfAtomPairs );
// ------------------------------------------------------------------
//! Coverts 3 char Res Name to 1 char Res Name
/*!
 
*/
void Res3ToRes1(char *cResCode3, char *cResCode1);
// ------------------------------------------------------------------
//======================================================================
#endif
#endif
