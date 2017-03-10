// -------------------------------------------------------------------------
// Copyright (C) 2009  Francesca Fanelli, Angelo Felline
// University of Modena and Reggio Emilia - ITALY
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
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
// -------------------------------------------------------------------------

#ifndef SURF
#define SURF

struct inp_surf
{
  // === Common to ARVO and GEPOL Algorithms ===========================
  int         iWriteSurfFlag;                    // If 1 write surface in output file/standard out
  int         iAlgoFlag;                         // If 0 use ARVO Algorithm, else use GEPOL Algorithm  
  
  int         iNumOfAtoms;                       // Number of atoms used in calculation
  int         iNumOfRealAtoms;                   // Number of atoms in PDB file
  
  float       dSolventRadius;                    // Solvent molecule radius  
  char        cTitle[64];                        // Output column header string    
  Selection   sele1;                             // Wordom atom selection facility

  // === Used by GEPOL Algorithm =======================================
  int         iMaxSize;                          // Vectors and arrays length
  int         memsize;                           // permitted memsize
  int         iSurfType;                         // 0=vdW; 1=ASA; 2=ESA
  int         iGhostFlag;                        // 0 if there are no ghost atoms
  
  int         iSolvRadAddedFalg;                 // Used to add solvent radius only one time
  
  int         iNDIV;                             // Diviosion Level
  float       fOFAC;                             // Overlapping Factor
  float       fRMIN;                             // Radius of smallest sphere
  
  int         ppiJVT1[3][60];                    // This has the information about the vertices
  int         ppiJVT2[3][4];                     // This has the information about the vertices
  
  int        *piAtomType;                        /* Vector of length iNumOfRealAtoms
                                                  * 
                                                  * 1 Initial center with radius == 0
                                                  * 2 Sphere engulfed by another
                                                  * 3 Ghost sphere
                                                  * 4 Sphere with final area == 0
                                                  * 5 Semi ghost sphere
                                                  * 6 Real sphere
                                                  */
                                                  
  int        *piIndexVect;                       // Used in GeoCav
  
  float     **ppfAtomCoord;                      // Atom coordinates
  float      *pfDistVect;                        // Used in Bulk function
  int        *piEngulfedAtoms;                   // Engulfed atoms list, used in Bulk
  int        *piZeroSurf;                        // Atoms with accessible surface area zero, used in Bulk
  int        *piZeroSurfMZ5;                     // Atoms with accessible surface area zero, used in MZero5
  int        *piZeroSurfC5;                      // Atoms with accessible surface area zero, used in Clean5
  int        *piZeroCheck;                       // Atoms with accessible surface area zero, used in MZero5

  double      ppdCV[32][3];                      // ?
  double      pdTriaCenter[3];                   // Center coordinates of a spherical triangle, calculated by CalcTriaCenter
  
  int         iNumOfTesserae;                    // = 60*pow(4, (iNDIV-1))
  int         iNumOfPoints;                      // Number of Points
  
  float      *pfXCoord1;                         // Tesserae X Coordinates
  float      *pfYCoord1;                         // Tesserae Y Coordinates
  float      *pfZCoord1;                         // Tesserae Z Coordinates
  
  int        *piUsefulAtoms;                     // Used in Clean5
  
  float     **ppfPCoord;                         // Points coordinates, used in GeoCav
  float      *pfAPVect;                          // Used in GeoCav
  
  // === Used by ARVO Algorithm ========================================
  int         iCirclesPlaneForAtom;              // Maximum number of local circles in the plane for one atom
  int         iArcsAnglesCircleIntersect;        // Maximum number of arcs and angles which arise from the local circles intersections
  int         iNumOfNeighbors;                   // Maximum number of neighbors
  
  double     *pdAtomsRadii;                      // Vector with all atom radii 
  double    **ppdCircles;                        // Used in MakeCircles
  double    **ppdArcs;                           // Used in MakeCircles
  double    **ppdNewArcs;                        // Used in CircleToArcs
  double     *pdAngles;                          // Used in MakeNewArcs
  double     *pdNewAngles;                       // Used in DeleteEqual
  double    **ppdAtomNeighbors;                  // List of stom neighbors
  
  double      dIntPointAngleA1;                  // 
  double      dIntPointAngleA2;                  // Angles of two intersection points
  double      dIntPointAngleB1;                  // 
  double      dIntPointAngleB2;                  // 
    
  int        *piNeighborsNumber;                 // Number of Neighbors for all atoms
  int        *piNeighborsStartIndex;             // Start of neighbors indices for ith atom in vector piNeighborsIndices
  int        *piNeighborsIndices;                // Vector of neighbors indices for each atom
  int        *piIndicesVector;                   // Used in MakeNeighbors
  int        *piIndicesVector2;                  // Used in GetSurf
  
  int        *piWithinAtoms;                     // List of atoms within 10 A from selected atoms
  float     **ppfWithinAtomsCoord;               // Coordinates of atoms within 10 A from selected atoms
   
  // === Used by ARVO and GEPOL Algorithms in SURFCORR Module ==========
  int         iNumOfFrames;                      // Number of frames in passed trajectory
  int         iFrameNum;                         // Number of used frames
  
  int         iWriteRawFlag;                     // if 1 writes sufaces at the end of output file
  
  float      *pfSele1Surf;                       // Used if iCorrType == 0, stores sele1 surface values
  float      *pfSele2Surf;                       // Used if iCorrType == 0, stores sele2 surface values
  
  Selection    sele2;                             // Wordom atom selection facility
  
  // === Used by ARVO and GEPOL Algorithms in SURFCLUST Modeule ========
  float       fClustBin;                         // Clust Bin
  float       fSmallestSurf;                     // Smallest Surface
  float       fLargestSurf;                      // Biggest Surface

};
// =====================================================================

// === Function prototypes =============================================

// === Common Functions ================================================
int    Read_SURF(char **input, int inp_index, struct inp_surf *inp_surf, Molecule *molecule, char *outstring);
int    Compute_SURF(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring);
// === GEPOL Functions =================================================
double ComputeGEPOLSurf(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring);
void   FillJVTArrays(struct inp_surf *inp_surf);
void   CompTriaVertexCoord(struct inp_surf *inp_surf);
void   Divide(struct inp_surf *inp_surf);
void   CalcTriaCenter(struct inp_surf *inp_surf, double dCoord1X, double dCoord1Y, double dCoord1Z, double dCoord2X, double dCoord2Y, double dCoord2Z, double dCoord3X, double dCoord3Y, double dCoord3Z);
void   TriaToFour(double ppdTmpCV[6][3]);
void   Shell(struct inp_surf *inp_surf);
void   Bulk(struct inp_surf *inp_surf);
int    MZero5(struct inp_surf *inp_surf, int iNumOfStartingAtoms);
void   Clean5(struct inp_surf *inp_surf);
void   GenAtoms(struct inp_surf *inp_surf);
void   AddSolvRad(struct inp_surf *inp_surf);
void   GeoCav(struct inp_surf *inp_surf);
float  AreaSum(struct inp_surf *inp_surf);
// === ARVO Functions ==================================================
double ComputeARVOSurf(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring);
void   MakeNeighbors(struct inp_surf *inp_surf);
int    Neighbors(int iIndex, struct inp_surf *inp_surf);
int    NorthPoleTest(struct inp_surf *inp_surf);
void   MolRotation(struct inp_surf *inp_surf);
double GetSurf(int iAtomNumber, struct inp_surf *inp_surf);
void   MakeCircles(struct inp_surf *inp_surf, int iNumOfLocalAtoms);
int    CircleToArcs(struct inp_surf *inp_surf, int iNumOfLocalAtoms);
int    MakeNewArcs(int iIndex, int iNumOfLocalAtoms, struct inp_surf *inp_surf);
void   CirclesIntersection(int iIndex1, int iIndex2, struct inp_surf *inp_surf);
int    CircleInCircle(int iIndex1, int iIndex2, struct inp_surf *inp_surf);
void   DirectSort(int iNumOfAngles, struct inp_surf *inp_surf);
void   InvertedSort(int iNumOfAngles, struct inp_surf *inp_surf);
int    DeleteEqual(int iNumOfAngles, struct inp_surf *inp_surf);
int    PointInCircle(double dCompT, double dCompS, int iJJIndex, struct inp_surf *inp_surf);
double AvgIntegral(int iNumOfArcs, double dAtomCoordZ, double dAtomRadius, struct inp_surf *inp_surf);
double Fract(double dA, double dB, double dC, double dSinPhi, double dCosPhi, int iExp);
// =====================================================================
// === Surface-based clustering ==================================================================================
int    Read_SURFCLUST(char **input, int inp_index, struct inp_surf *inp_surf, Molecule *molecule, char *outstring, int iNumOfFrames);
int    Compute_SURFCLUST(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring);
void   Compute_ARVO_SURFCLUST(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring);
void   Compute_GEPOL_SURFCLUST(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring);
int    Post_SURFCLUST(struct inp_surf *inp_surf, struct sopt *OPT);
// ===============================================================================================================
// === Surface Correlation =======================================================================================
int    Read_SURFCORR(char **input, int inp_index, struct inp_surf *inp_surf, Molecule *molecule, char *outstring, int iNumOfFrames);
int    Compute_SURFCORR(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring);
void   Compute_ARVO_SURFCORR(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring);
void   Compute_GEPOL_SURFCORR(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring);
int    Post_SURFCORR(struct inp_surf *inp_surf, struct sopt *OPT);
// ===============================================================================================================
#endif
