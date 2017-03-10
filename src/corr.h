/* =====================================================================
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
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
   ===================================================================== */

struct inp_corr
{
  // Almost General Variables
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
  int         iCorrType;                                                // Correlation Type: 0: DCCM; 1: LMI
  int         iFrameNum;                                                // Frame number
  
  char        cTitle[512];                                               // Output column header string
  char      **pcLabels;                                                 // Labels (segname:restype+resnum), used if iMultiAtomFlag == 1
  
  float       pfAxesDist[3];                                            // Used to store distances if PBC
  
  FILE       *FOutFile;                                                 // Output File Handler
  
  Selection   sele;                                                     // Selection Structure
  
  // DCC Variables
  float     **ppfMeanDeviation;                                         // Atom coord mean deviations
  float      *pfSquareMeanDeviation;                                    // Atom coord square mean deviations
  float      *pfMasses;                                                 // Used to store atom masses, used if iMultiAtomFlag and iMassFlag == 1 
  float     **ppfVirtAtomCoord;                                         // Frame Virtual Atom Coordinates, used if iMultiAtomFlag == 1
  float     **ppfVirtRefAtomCoord;                                      // Reference Virtual Atom Coordinates, used if iMultiAtomFlag == 1
  
  // LMI Variables
  double     *pdCovMatrix;                                              // Covariance Matrix
  double     *pdOutput;
  
  
  double   ***pppdCoord;                                                // Holds selected atoms coordinates
};

// General Functions
int     Read_CORR(char **input, int input_index, struct inp_corr *inp_corr, Molecule *molecule, char *outstring, int iNumOfFrames);
int     Compute_CORR(struct inp_corr *inp_corr, Molecule *molecule, CoorSet *trj_crd, char *outstring);
int    Post_CORR (struct inp_corr *inp_corr, Molecule *molecule, struct sopt *OPT, CoorSet *trj_crd);

// MI & LMI Functions
void     CovarMatrix(struct inp_corr *inp_corr);
void     GaussCorrMatrix(struct inp_corr *inp_corr);
void     SubMat(double *pdDestMatrix, double *pdSourceMatrix, int iSourceMatrixDim, int iStartRow, int iEndRow, int iStartCol, int iEndCol);
void     PushMat(double *pdDestMatrix, double *pdSourceMatrix, int iStartRow, int iStartCol, int iTransposeFlag);
double   CalcEntropy(double *pdMatrix, int iMatrixDim);
double   MatDet(double **ppdMatrix, int iMatrixSize);
double **VectToMat(double *pdVector, int iVectorSize);
