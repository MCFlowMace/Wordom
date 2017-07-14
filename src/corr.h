/* =========================================================================
   Copyright (C) 2009  Francesca Fanelli, Angelo Felline
   University of Modena and Reggio Emilia - ITALY
  
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
   ========================================================================= */
   
   
/*
  The algorithms implemented in this module were written from scratch
  following the methods described in these papers:
    DCC: 
    * J.A. McCammon and S.C. Harvey, Dynamics of proteins and nucleic acids; Cambridge Univ Pr, 1988.
    
    LMI
    * A. Kraskov, H. Stoegbauer, and P. Grassberger, Physical Review E, 2004, 69(6), 66138.
    * O.F. Lange and H. Grubmuller, PROTEINS-NEW YORK-, 2006, 62(4), 1053.
    
    DCOR:
    * A. Roy and C. B. Post, JCTC. 2012, 8, 3009−3014
*/

struct inp_corr
{
  // Almost General Variables
  int         iPDBFlag;                                                 // 1 if PBC
  int         iNumOfFrames;                                             // Number of frames
  int         iNumOfAtoms;                                              // Number of atoms
  int         iNumOfRes;                                                // Number of residues
  int         iNumOfSelAtoms;                                           // Number of selected atoms
  int         iNumOfSelRes;                                             // Number of selected residues
  int         iMassFlag;                                                // 1 = Mass Centre, 0 = Geometrical Centre
  int         iResFlag;                                                 // 1 = Correlates by residues, 0 = Correlates by atoms
  int         iMultiAtomFlag;                                           // 1 = More than 1 atom per residue
  int         iMatrixDim;                                               // Correlation Matrix dimension
  int        *piSelResLen;                                              // Number of atoms per selected residue 
  int         iZeroMassWarning;                                         // 1 if one or more atoms have mass = 0.0
  int        *piProgNum;                                                // Progressive res numbers, used if iMultiAtomFlag == 1
  int         iCorrType;                                                // Correlation Type: 0: DCCM; 1: LMI, 2: DCOR
  int         iFrameNum;                                                // Frame number
  int         iNumOfThreads;                                            // Number of threads
  int         iNumOfFrmPerThreads;                                      // Used in multithreading, number of frames per thread
  int       **ppiFrmListPerThreads;                                     // Frames list of each threads
  
  int         iNumOfPairsPerThreads;                                    // Used in multithreading, number of atom/residue pairs per thread
  int       **ppiPairsListFrm1;                                         // Frame pairs list of each threads
  int       **ppiPairsListFrm2;                                         // Frame pairs list of each threads
  
  char        cTitle[512];                                              // Output column header string
  char      **pcLabels;                                                 // Labels (segname:restype+resnum), used if iMultiAtomFlag == 1
  
  float       pfAxesDist[3];                                            // Used to store distances if PBC
  float      *pfMasses;                                                 // Used to store atom masses, used if iMultiAtomFlag and iMassFlag == 1 
  float     **ppfVirtAtomCoord;                                         // Frame Virtual Atom Coordinates, used if iMultiAtomFlag == 1
  float     **ppfVirtRefAtomCoord;                                      // Reference Virtual Atom Coordinates, used if iMultiAtomFlag == 1
  
  FILE       *FOutFile;                                                 // Output File Handler
  FILE       *FVerbOutFile;                                             // Verbose Output File Handler
  
  Selection   sele;                                                     // Selection Structure
  
  pthread_t  *corr_threads;                                             // Used in multithreading computations (DCor only at the moment)
  
  Pbc        *pPBCInfo;
  
  // DCC Variables
  float     **ppfMeanDeviation;                                         // Atom coord mean deviations
  float      *pfSquareMeanDeviation;                                    // Atom coord square mean deviations

  
  // LMI Variables
  double     *pdCovMatrix;                                              // Covariance Matrix
  double     *pdOutput;

  // DCOR Variables
  int       iThreadResNum;                                              // Used by threaded computations
  double  **ppfResFrameAvgDist;                                         // Average dist of each frame of each selected atoms
  double ***pppfThreadResFrameAvgDist;                                  // Average dist of each frame of each selected atoms, one for each thread
  double   *pfResOverallAvgDist;                                        // Overall average dist of each residue
  double  **ppfThreadResOverallAvgDist;                                 // Overall average dist of each residue, one for each thread
  double ***pppfThreadDistCov;                                          // Distance covariance matrix, one for each thread
  double  **ppfDistCov;                                                 // Final Distance covariance matrix
  double  **ppfDCorr;                                                   // Distance correlation matrix
  
  // FLUCT Variables
  int        iFirstRoundFlag;                                           // Flag used to rewind the trj
  int        iFirstFrame;                                               // Flag used to save min and max distances
  int        iVerboseFlag;                                              // Flag used to print some useful (?) info
  int        iNumOFSubSele;                                             // Number of sub selections for overall fluct calculations
  int        iMatchSubSele;                                             // If 1, all sub-seles will be pairwise matched to calculated cross-sub-sele-fluct
  int        iNumOfSubMatch;                                            // Number of sub-sale matches
  int        iStdDevFlag;                                               // If 1, ppfFluctMatrix will be filled with standard-dev instead of variance
  int        iGetFramesFlag;                                            // If 1 print the avg distances of each [sub]seles and matches
  int       *piMasterSeleIndexes;                                       // Used to map master-sele to molecule's residues/atoms
  double   **ppfMeanDistance;                                           // Stores the mean dist between two residues
  double   **ppfFluctMatrix;                                            // Fluctuation matrix
  double   **ppfFrameMatrix;                                            // Fluctuation matrix of a single frame
  double   **ppfMinDistances;                                           // Stores the min dist between two residues
  double   **ppfMaxDistances;                                           // Stores the max dist between two residues
  double     fOverallFluct;                                             // Overall fluctation calculated over the main selection
  double     fOverallAvgDist;                                           // Average Overall distance over the main selection
  double    *pfSubSeleFluct;                                            // Overall fluctation calculated over sub-selection
  double    *pfSubSeleAvgDist;                                          // Overall fluctation calculated over sub-selection
  double    *pfMatchSubSeleFluct;                                       // Pairwise sub-sele Overall fluctations
  double    *pfMatchSubSeleAvgDist;                                     // Pairwise sub-sele Overall fluctations
  Selection *pSubSele;
  
  // LMI & DCOR                                                         
  double   ***pppdCoord;                                                // Holds selected atoms coordinates
  
};

// General Functions
int    Read_CORR(char **input, int input_index, struct inp_corr *inp_corr, Molecule *molecule, char *outstring, int iNumOfFrames);
int    Compute_CORR(struct inp_corr *inp_corr, Molecule *molecule, CoorSet *trj_crd, char *outstring);
int    Post_CORR (struct inp_corr *inp_corr, Molecule *molecule, struct sopt *OPT, CoorSet *trj_crd);

// LMI Functions
void     CovarMatrix(struct inp_corr *inp_corr);
void     GaussCorrMatrix(struct inp_corr *inp_corr);
void     SubMat(double *pdDestMatrix, double *pdSourceMatrix, int iSourceMatrixDim, int iStartRow, int iEndRow, int iStartCol, int iEndCol);
void     PushMat(double *pdDestMatrix, double *pdSourceMatrix, int iStartRow, int iStartCol, int iTransposeFlag);
double   CalcEntropy(double *pdMatrix, int iMatrixDim);
double   MatDet(double **ppdMatrix, int iMatrixSize);
double **VectToMat(double *pdVector, int iVectorSize);

// DCor Functions
void    *DCorStep1(void *pviThreadData);                                // Multithreading function
void    *DCorStep2(void *pviThreadData);                                // Multithreading function
