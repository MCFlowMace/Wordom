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
#include "wordom.h"
#include "tools.h"
#include "fileio.h"
#include "datahandler.h"
#include "time.h"
#include "corr.h"

#define  VERSTRING "0.3b"
#define  PI        3.14159265358979323846

void CovarMatrix(struct inp_corr *inp_corr)
{
  int     ii, jj;
  int     iAtom1, iAtom2;
  int     iCoord1, iCoord2;
  int     iFrameNum;
  
  for(iAtom1=0, ii=0; iAtom1<inp_corr->iNumOfSelAtoms; iAtom1++)
  {
    for(iCoord1=0; iCoord1<3; iCoord1++, ii++)
    {
      for(iAtom2=0, jj=0; iAtom2<inp_corr->iNumOfSelAtoms; iAtom2++)
      {
        for(iCoord2=0; iCoord2<3; iCoord2++, jj++)
        {
          inp_corr->pdCovMatrix[ii*inp_corr->iNumOfSelAtoms*3+jj] = 0.0;
          for(iFrameNum=0; iFrameNum<inp_corr->iNumOfFrames; iFrameNum++)
          {
            inp_corr->pdCovMatrix[ii*inp_corr->iNumOfSelAtoms*3+jj] += (inp_corr->pppdCoord[iAtom1][iFrameNum][iCoord1]*
                                                                        inp_corr->pppdCoord[iAtom2][iFrameNum][iCoord2]);
          }
          inp_corr->pdCovMatrix[ii*inp_corr->iNumOfSelAtoms*3+jj] /= inp_corr->iNumOfFrames;
        }
      }
    }
  }
  
}

void SubMat(double *pdDestMatrix, double *pdSourceMatrix, int iSourceMatrixDim, int iStartRow, int iEndRow, int iStartCol, int iEndCol)
{
  /*
   * This function copies elements from iStartRow to iEndRow-1 and
   * from iStartCol to iEndCol-1 from the iSourceMatrixDim^2 dimension
   * pdSourceMatrixDim matrix into pdDestMatrix matrix
   */
  
  int     ii, jj;
  int     iRowNum, iColNum;
  
  for(iRowNum=iStartRow, ii=0; iRowNum<iEndRow; iRowNum++, ii++)
  {
    for(iColNum=iStartCol, jj=0; iColNum<iEndCol; iColNum++, jj++)
    {
      pdDestMatrix[ii*(iEndCol-iStartCol)+jj] = pdSourceMatrix[iRowNum*iSourceMatrixDim+iColNum];
    }
  }
}

void PushMat(double *pdDestMatrix, double *pdSourceMatrix, int iStartRow, int iStartCol, int iTransposeFlag)
{
  /*
   * This function copies all elements of 3x3 matrix SourceMatrix
   * into the 3x3 submatrix of the 6x6 matrix DestMatrix which
   * starts at row iStartRow and col iStartCol
   */

  int     ii, jj;
  int     iRowNum, iColNum;
  
  for(iRowNum=iStartRow, ii=0; iRowNum<(iStartRow+3); iRowNum++, ii++)
  {
    for(iColNum=iStartCol, jj=0; iColNum<(iStartCol+3); iColNum++, jj++)
    {
      if(iTransposeFlag == 1)
        pdDestMatrix[iRowNum*3*2+iColNum] = pdSourceMatrix[jj*3+ii];
        
      else
        pdDestMatrix[iRowNum*3*2+iColNum] = pdSourceMatrix[ii*3+jj];
    }
  }
}

double **VectToMat(double *pdVector, int iVectorSize)
{
  /*
   * This function transforms vector pdVector
   * of size iVectorSize in a matrix of
   * size sqrt(iVectorSize) x sqrt(iVectorSize)
   */
  
  int	     ii;
  int      iRowNum, iColNum, iPos;
  int      iMatrixSize;
  
  double **ppdMatrix;
    
  iMatrixSize = (int) sqrt(iVectorSize);
  ppdMatrix = (double **) malloc(iMatrixSize*sizeof(double *));
  for(ii=0; ii<iMatrixSize; ii++)
    ppdMatrix[ii] = (double *) malloc(iMatrixSize*sizeof(double));


  iRowNum =  0;
  iColNum = -1;
  iPos    = -1;
  
  for(ii=0; ii<iVectorSize; ii++)
  {
    iPos++;
    if(iPos == iMatrixSize)
    {
      iRowNum++;
      iColNum = 0;
      iPos    = 0;
    }
    
    else
      iColNum++;
    
    ppdMatrix[iRowNum][iColNum] = pdVector[ii];
  }

  return ppdMatrix;
}

double MatDet(double **ppdMatrix, int iMatrixSize)
{
  /*
   * An __unefficient__ algorithm to
   * calculate  matrix  determinants
   */
  
  int	     ii, jj, kk, mm;
  double   dDet;
  double **ppdSubMatrix;
  
  if(iMatrixSize == 2)
    dDet = (ppdMatrix[0][0] * ppdMatrix[1][1] - ppdMatrix[1][0] * ppdMatrix[0][1]);
  
  else
  {
    dDet = 0.0;
    for (kk=0; kk<iMatrixSize; kk++)
    {
      ppdSubMatrix = (double **) malloc((iMatrixSize-1)*sizeof(double *));
      for(ii=0; ii<iMatrixSize-1; ii++)
        ppdSubMatrix[ii] = (double *) malloc((iMatrixSize-1)*sizeof(double));
      
      for(ii=1; ii<iMatrixSize; ii++)
      {
        mm = 0;
        for(jj=0; jj<iMatrixSize; jj++)
        {
          if(jj == kk)
            continue;

          ppdSubMatrix[ii-1][mm] = ppdMatrix[ii][jj];
          mm++;
        }
      }
      
      dDet += pow(-1.0,1.0+kk+1.0) * ppdMatrix[0][kk] * MatDet(ppdSubMatrix,iMatrixSize-1);
      
      for (ii=0; ii<iMatrixSize-1; ii++)
        free(ppdSubMatrix[ii]);
      free(ppdSubMatrix);
    }
  }
  
  return dDet;
}

double  CalcEntropy(double *pdMatrix, int iMatrixDim)
{
  // Calculates entropy of iMatrixDim-variate Gaussian with Covariance matrix pdMatrix
  int      ii;
  int      iTmpMatrixSize;
  
  double   dDet;
  double **ppdTmpMatrix;
  
  iTmpMatrixSize = iMatrixDim*iMatrixDim;
  ppdTmpMatrix = VectToMat(pdMatrix, iTmpMatrixSize);
  
  dDet = MatDet(ppdTmpMatrix, iMatrixDim);
  
  for(ii=0; ii<iMatrixDim; ii++)
    free(ppdTmpMatrix[ii]);
  free(ppdTmpMatrix);

  return 0.5*(iMatrixDim*(1.0+log(2.0*PI))+log(dDet));
}

void GaussCorrMatrix(struct inp_corr *inp_corr)
{
  int     iAtom1, iAtom2;
  int     iDimension;
  
  double  pdSubMatDim3[9];
  double  pdSubMatDim6[36];
  double  dEntropy1, dEntropy2, dEntropy3;
  
  iDimension = inp_corr->iNumOfSelAtoms * 3;
  
  for(iAtom1=0; iAtom1<inp_corr->iNumOfSelAtoms; iAtom1++)
  {
    inp_corr->pdOutput[iAtom1*inp_corr->iNumOfSelAtoms+iAtom1] = 2000.0;
    for(iAtom2=iAtom1+1; iAtom2<inp_corr->iNumOfSelAtoms; iAtom2++)
    {
      SubMat(pdSubMatDim3, inp_corr->pdCovMatrix, iDimension, iAtom1*3, iAtom1*3+3, iAtom1*3, iAtom1*3+3);
      PushMat(pdSubMatDim6, pdSubMatDim3, 0, 0, 0);
      dEntropy1 = CalcEntropy(pdSubMatDim3, 3);
      
      SubMat(pdSubMatDim3, inp_corr->pdCovMatrix, iDimension, iAtom2*3, iAtom2*3+3, iAtom2*3, iAtom2*3+3);
      PushMat(pdSubMatDim6, pdSubMatDim3, 3, 3, 0);
      dEntropy2 = CalcEntropy(pdSubMatDim3, 3);
      
      SubMat(pdSubMatDim3, inp_corr->pdCovMatrix, iDimension, iAtom1*3, iAtom1*3+3, iAtom2*3, iAtom2*3+3);
      PushMat(pdSubMatDim6, pdSubMatDim3, 0, 3, 0);
      PushMat(pdSubMatDim6, pdSubMatDim3, 3, 0, 1);
      dEntropy3 = CalcEntropy(pdSubMatDim6, 6);
      
      inp_corr->pdOutput[iAtom2*inp_corr->iNumOfSelAtoms+iAtom1] = dEntropy1 + dEntropy2 - dEntropy3;
      inp_corr->pdOutput[iAtom1*inp_corr->iNumOfSelAtoms+iAtom2] = dEntropy1 + dEntropy2 - dEntropy3;
    }
  }
}

int Read_CORR(char **input, int input_index, struct inp_corr *inp_corr, Molecule *molecule, char *outstring, int iNumOfFrames)
{
  int     ii, jj;
  int     iWinpOptFlag=0;
  int     iResNum, iResCont, iVirtAtom;
  
  char    cWordomInpBuffer[1024], cTmpString[10];
  char    cLastRes[50], cOutPutFileName[512];
  char    cResCode1[3];
  
  float   fTotalWeight;
  float   fVirtAtomXCoord, fVirtAtomYCoord, fVirtAtomZCoord;

  extern short int     no_frame_par;
  no_frame_par = 1;
  
  memset ( cWordomInpBuffer, '\0', sizeof(cWordomInpBuffer));
  
  // === example winp file ===
  //BEGIN
  //--TITLE CORR
  //--TYPE DCCM
  //--SELE /*/*/*
  //--LEVEL RES
  //--TYPE DCC
  //--MASS 0
  //END
  // =========================
  
  // === Default Values ======
  sprintf(inp_corr->sele.selestring, "/*/*/CA");
  inp_corr->iMassFlag = 0;
  inp_corr->iResFlag  = 1;
  inp_corr->iCorrType = 0;
  // =========================
  
  while( strncmp (cWordomInpBuffer, "END", 3))
  {
    iWinpOptFlag = 0;
    sprintf( cWordomInpBuffer, "%s", input[input_index]);
    if( !strncmp(cWordomInpBuffer, "BEGIN", 5) || !strncmp(cWordomInpBuffer, "END", 3) || cWordomInpBuffer[0] == '#')
      iWinpOptFlag = 1;
    else if ( !strncmp(cWordomInpBuffer, "--TITLE", 7))
    {
      sscanf( cWordomInpBuffer, "--TITLE %s", inp_corr->cTitle);
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--TYPE", 6))
    {
      sscanf( cWordomInpBuffer, "--TYPE %s", cTmpString);
      if(strcmp(cTmpString, "DCC")==0)
        inp_corr->iCorrType=0;
      else if(strcmp(cTmpString, "LMI")==0)
        inp_corr->iCorrType=1;
      else
      {
        fprintf( stderr, "CORR module: Could NOT understand --TYPE option: %s\n", cTmpString);
        exit(5);
      }
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--LEVEL", 7))
    {
      sscanf( cWordomInpBuffer, "--LEVEL %s", cTmpString);
      if(strcmp(cTmpString, "RES")==0)
        inp_corr->iResFlag=1;
      else if(strcmp(cTmpString, "ATM")==0)
        inp_corr->iResFlag=0;
      else
      {
        fprintf( stderr, "CORR module: Could NOT understand --LEVEL option: %s\n", cTmpString);
        exit(5);
      }
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--MASS", 6))
    {
      sscanf(cWordomInpBuffer, "--MASS %d", &inp_corr->iMassFlag);
      if(inp_corr->iMassFlag==1)
        inp_corr->iMassFlag=1;
      else
        inp_corr->iMassFlag=0;
        
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--SELE", 6))
    {
      sscanf(cWordomInpBuffer, "--SELE %[^\n]%*c ", inp_corr->sele.selestring);
      iWinpOptFlag = 1;
    }
    if( iWinpOptFlag==0 )
    {
      fprintf( stderr, "CORR module: Could NOT understand option: %s\n", cWordomInpBuffer);
      exit(5);
    }
    input_index++;
  }
  
  GetSele(inp_corr->sele.selestring, &inp_corr->sele, molecule);
  if(inp_corr->sele.nselatm == 0)
  {
    fprintf( stderr, "CORR module: Empty selection: %s\n", inp_corr->sele.selestring);
    exit(5);
  }
  
  // === General Settings ==============================================
  inp_corr->iNumOfFrames=iNumOfFrames;
  
  inp_corr->iNumOfRes=molecule->nRes;
  inp_corr->iNumOfAtoms=molecule->nato;

  inp_corr->iNumOfSelRes=0;
  inp_corr->iNumOfSelAtoms=inp_corr->sele.nselatm;
  inp_corr->iMultiAtomFlag=0;
  
  cTmpString[0]='\0';
  cLastRes[0]='\0';
  
  if(inp_corr->iResFlag==1)                                             // Calculates Correlation by Residues
  {                                                                     // Calculates the number of selected residues
    for(ii=0; ii<inp_corr->sele.nselatm; ii++)
    {
      sprintf(cTmpString, "%s%d", molecule->rawmol.segId[inp_corr->sele.selatm[ii]-1], 
                                  molecule->rawmol.resn[inp_corr->sele.selatm[ii]-1]);
      if(strcmp(cTmpString, cLastRes)!=0)
      {
        inp_corr->iNumOfSelRes++;
        strcpy(cLastRes, cTmpString);
      }
      else
        inp_corr->iMultiAtomFlag=1;                                     // More than 1 atom per residue
    }
    inp_corr->iMatrixDim=inp_corr->iNumOfSelRes;
  }
  
  else                                                                  // Calculates Correlation by Atoms
    inp_corr->iMatrixDim=inp_corr->iNumOfSelAtoms;

  strcpy(cOutPutFileName, inp_corr->cTitle);
  inp_corr->FOutFile=fopen(cOutPutFileName, "w");
  
  if(inp_corr->iMassFlag==1)
  {
    inp_corr->pfMasses = (float *) calloc(inp_corr->sele.nselatm, sizeof(float));
    for(ii=0; ii<inp_corr->sele.nselatm; ii++)
      inp_corr->pfMasses[ii] = molecule->rawmol.bFac[inp_corr->sele.selatm[ii]-1];
      
    inp_corr->iZeroMassWarning=0;
    for(ii=0; ii<inp_corr->sele.nselatm; ii++)
    {
      inp_corr->pfMasses[ii] = molecule->rawmol.bFac[inp_corr->sele.selatm[ii]-1];
      if(inp_corr->pfMasses[ii]==0.0)
        inp_corr->iZeroMassWarning=1;
    }
  }
  
  if(inp_corr->iMultiAtomFlag==1)                                       // Calculates the number of atoms per selected residue
  {
    inp_corr->piSelResLen = (int   *) calloc(inp_corr->iMatrixDim, sizeof(int   ));
    inp_corr->piProgNum   = (int   *) calloc(inp_corr->iMatrixDim, sizeof(int   ));
    inp_corr->pcLabels    = (char **) calloc(inp_corr->iMatrixDim, sizeof(char *));
    for(ii=0; ii<inp_corr->iMatrixDim; ii++)
      inp_corr->pcLabels[ii] = (char *) calloc(20, sizeof(char));
    
    cTmpString[0]='\0';
    cLastRes[0]='\0';
    sprintf(cLastRes, "%s%d", molecule->rawmol.segId[inp_corr->sele.selatm[0]-1], 
                              molecule->rawmol.resn[inp_corr->sele.selatm[0]-1]);
    
    iResNum=-1;
    iResCont=-1;

    Res3ToRes1(molecule->rawmol.restype[inp_corr->sele.selatm[0]-1], cResCode1);
    sprintf(inp_corr->pcLabels[0], "%s:%s%d", molecule->rawmol.segId[inp_corr->sele.selatm[0]-1], cResCode1,
                                              molecule->rawmol.resn[inp_corr->sele.selatm[0]-1]);
                                                        
    inp_corr->piProgNum[0]=molecule->rawmol.presn[inp_corr->sele.selatm[0]-1];
    
    for(ii=0; ii<inp_corr->sele.nselatm; ii++)
    {
      iResCont++;
      sprintf(cTmpString, "%s%d", molecule->rawmol.segId[inp_corr->sele.selatm[ii]-1], 
                                  molecule->rawmol.resn[inp_corr->sele.selatm[ii]-1]);
      
      if(strcmp(cTmpString, cLastRes)!=0)
      {
        iResNum++;
        inp_corr->piSelResLen[iResNum]=iResCont;
        strcpy(cLastRes, cTmpString);
        Res3ToRes1(molecule->rawmol.restype[inp_corr->sele.selatm[ii]-1], cResCode1);
        sprintf(inp_corr->pcLabels[iResNum+1], "%s:%s%d", molecule->rawmol.segId[inp_corr->sele.selatm[ii]-1], cResCode1,
                                                          molecule->rawmol.resn[inp_corr->sele.selatm[ii]-1]);
                                                        
        inp_corr->piProgNum[iResNum+1]=molecule->rawmol.presn[inp_corr->sele.selatm[ii]-1];
        iResCont=0;
      }
    }
    iResNum++;
    inp_corr->piSelResLen[iResNum]=iResCont+1;
    
    inp_corr->ppfVirtAtomCoord    = (float **) calloc(inp_corr->iMatrixDim, sizeof(float *));
    inp_corr->ppfVirtRefAtomCoord = (float **) calloc(inp_corr->iMatrixDim, sizeof(float *));
    for(ii=0; ii<inp_corr->iMatrixDim; ii++)                            // Stores Geometrical/Mass centre coordinates
    {
      inp_corr->ppfVirtAtomCoord[ii]    = (float *) calloc(3, sizeof(float));
      inp_corr->ppfVirtRefAtomCoord[ii] = (float *) calloc(3, sizeof(float));
    }
    
    // === Calculates geo/mass centres of reference molecule res =======
    iResNum=0;
    iResCont=0;
    iVirtAtom=-1;
    fTotalWeight=0.0;
    
    fVirtAtomXCoord = 0.0;
    fVirtAtomYCoord = 0.0;
    fVirtAtomZCoord = 0.0;
    
    for(ii=0; ii<inp_corr->sele.nselatm; ii++)
    {
      if(iResCont==inp_corr->piSelResLen[iResNum])
      {
        if(inp_corr->iMassFlag==0)
        {
          fVirtAtomXCoord = (fVirtAtomXCoord/(float)inp_corr->piSelResLen[iResNum]);
          fVirtAtomYCoord = (fVirtAtomYCoord/(float)inp_corr->piSelResLen[iResNum]);
          fVirtAtomZCoord = (fVirtAtomZCoord/(float)inp_corr->piSelResLen[iResNum]);
        }
        else
        {
          fVirtAtomXCoord = (fVirtAtomXCoord/fTotalWeight);
          fVirtAtomYCoord = (fVirtAtomYCoord/fTotalWeight);
          fVirtAtomZCoord = (fVirtAtomZCoord/fTotalWeight);
          
          fTotalWeight = 0.0;
        }
        
        iResNum++;
        iResCont=0;
        iVirtAtom++;
        
        inp_corr->ppfVirtRefAtomCoord[iVirtAtom][0]=fVirtAtomXCoord;
        inp_corr->ppfVirtRefAtomCoord[iVirtAtom][1]=fVirtAtomYCoord;
        inp_corr->ppfVirtRefAtomCoord[iVirtAtom][2]=fVirtAtomZCoord;
        
        fVirtAtomXCoord = 0.0;
        fVirtAtomYCoord = 0.0;
        fVirtAtomZCoord = 0.0;
      }
      
      if(inp_corr->iMassFlag==0)
      {
        // Geometrical centre
        fVirtAtomXCoord += molecule->coor.xcoor[inp_corr->sele.selatm[ii]-1];
        fVirtAtomYCoord += molecule->coor.ycoor[inp_corr->sele.selatm[ii]-1];
        fVirtAtomZCoord += molecule->coor.zcoor[inp_corr->sele.selatm[ii]-1];
      }
      else
      {
        // Mass centre
        fVirtAtomXCoord += (molecule->coor.xcoor[inp_corr->sele.selatm[ii]-1] * inp_corr->pfMasses[ii]);
        fVirtAtomYCoord += (molecule->coor.ycoor[inp_corr->sele.selatm[ii]-1] * inp_corr->pfMasses[ii]);
        fVirtAtomZCoord += (molecule->coor.zcoor[inp_corr->sele.selatm[ii]-1] * inp_corr->pfMasses[ii]);
        
        fTotalWeight += inp_corr->pfMasses[ii];
      }
      iResCont++;
    }
    
    // Process the last residue
    if(inp_corr->iMassFlag==0)
    {
      fVirtAtomXCoord = (fVirtAtomXCoord/(float)inp_corr->piSelResLen[iResNum]);
      fVirtAtomYCoord = (fVirtAtomYCoord/(float)inp_corr->piSelResLen[iResNum]);
      fVirtAtomZCoord = (fVirtAtomZCoord/(float)inp_corr->piSelResLen[iResNum]);
    }
    else
    {
      fVirtAtomXCoord = (fVirtAtomXCoord/fTotalWeight);
      fVirtAtomYCoord = (fVirtAtomYCoord/fTotalWeight);
      fVirtAtomZCoord = (fVirtAtomZCoord/fTotalWeight);
      
      fTotalWeight = 0.0;
    }
        
    iVirtAtom++;
    inp_corr->ppfVirtRefAtomCoord[iVirtAtom][0]=fVirtAtomXCoord;
    inp_corr->ppfVirtRefAtomCoord[iVirtAtom][1]=fVirtAtomYCoord;
    inp_corr->ppfVirtRefAtomCoord[iVirtAtom][2]=fVirtAtomZCoord;
    // =================================================================
    
  }
  
  // ===================================================================
  
  if(inp_corr->iCorrType == 0)
  {
    // --TYPE DCC
    
    inp_corr->pfSquareMeanDeviation     = (float  *) calloc(inp_corr->iMatrixDim, sizeof(float  ));
    
    inp_corr->ppfMeanDeviation          = (float **) calloc(inp_corr->iMatrixDim, sizeof(float *));
    for(ii=0; ii<inp_corr->iMatrixDim; ii++)
      inp_corr->ppfMeanDeviation[ii]    = (float  *) calloc(inp_corr->iMatrixDim, sizeof(float  ));
  }
  
  else if(inp_corr->iCorrType == 1)
  {
    // --TYPE LMI
    
    if(inp_corr->iResFlag == 1)
      inp_corr->iNumOfSelAtoms = inp_corr->iNumOfSelRes;
    
    inp_corr->iFrameNum=-1;
    
    inp_corr->pppdCoord = (double ***) malloc(inp_corr->iNumOfSelAtoms * sizeof(double **));
    for(ii=0; ii<inp_corr->iNumOfSelAtoms; ii++)
    {
      inp_corr->pppdCoord[ii] = (double **) malloc(inp_corr->iNumOfFrames * sizeof(double *));
      for(jj=0; jj<inp_corr->iNumOfFrames; jj++)
      {
        inp_corr->pppdCoord[ii][jj] = (double *) malloc(3 * sizeof(double));
      }
    }
    
    inp_corr->pdCovMatrix = (double *) calloc(inp_corr->iNumOfSelAtoms*3*inp_corr->iNumOfSelAtoms*3, sizeof(double));
    inp_corr->pdOutput    = (double *) calloc(inp_corr->iNumOfSelAtoms*inp_corr->iNumOfSelAtoms,     sizeof(double));
    
    for(ii=0; ii<inp_corr->iNumOfSelAtoms*inp_corr->iNumOfSelAtoms; ii++)
      inp_corr->pdOutput[ii] = 1.0;
  }
  
  sprintf( outstring, " %10s ", inp_corr->cTitle);
  return 12;
}

int Compute_CORR(struct inp_corr *inp_corr, Molecule *molecule, CoorSet *trj_crd, char *outstring)
{
  int     ii, jj;
  int     iAtomNumA, iAtomNumB;
  int     iResNum, iResCont, iVirtAtom;
  int     iAtomNum;
  
  float   fXDiffA, fYDiffA, fZDiffA;
  float   fXDiffB, fYDiffB, fZDiffB;
  float   fVirtAtomXCoord, fVirtAtomYCoord, fVirtAtomZCoord;
  float   fTotalWeight;
  
  inp_corr->iFrameNum++;                                                // Update frame number
  
  if(inp_corr->iMultiAtomFlag == 1)
  {
    // More than one atom per residue
    
    // === Calculates geo/mass centres of current frame res ==========
    iResNum=0;
    iResCont=0;
    iVirtAtom=-1;
    fTotalWeight=0.0;
    
    fVirtAtomXCoord = 0.0;
    fVirtAtomYCoord = 0.0;
    fVirtAtomZCoord = 0.0;
    
    for(ii=0; ii<inp_corr->sele.nselatm; ii++)
    {
      if(iResCont==inp_corr->piSelResLen[iResNum])
      {
        if(inp_corr->iMassFlag==0)
        {
          fVirtAtomXCoord = (fVirtAtomXCoord/(float)inp_corr->piSelResLen[iResNum]);
          fVirtAtomYCoord = (fVirtAtomYCoord/(float)inp_corr->piSelResLen[iResNum]);
          fVirtAtomZCoord = (fVirtAtomZCoord/(float)inp_corr->piSelResLen[iResNum]);
        }
        else
        {
          fVirtAtomXCoord = (fVirtAtomXCoord/fTotalWeight);
          fVirtAtomYCoord = (fVirtAtomYCoord/fTotalWeight);
          fVirtAtomZCoord = (fVirtAtomZCoord/fTotalWeight);
          
          fTotalWeight = 0.0;
        }
        
        iResNum++;
        iResCont=0;
        iVirtAtom++;
        
        inp_corr->ppfVirtAtomCoord[iVirtAtom][0]=fVirtAtomXCoord;
        inp_corr->ppfVirtAtomCoord[iVirtAtom][1]=fVirtAtomYCoord;
        inp_corr->ppfVirtAtomCoord[iVirtAtom][2]=fVirtAtomZCoord;
        
        fVirtAtomXCoord = 0.0;
        fVirtAtomYCoord = 0.0;
        fVirtAtomZCoord = 0.0;
      }
      
      if(inp_corr->iMassFlag==0)
      {
        // Geometrical centre
        fVirtAtomXCoord += trj_crd->xcoor[inp_corr->sele.selatm[ii]-1];
        fVirtAtomYCoord += trj_crd->ycoor[inp_corr->sele.selatm[ii]-1];
        fVirtAtomZCoord += trj_crd->zcoor[inp_corr->sele.selatm[ii]-1];
      }
      else
      {
        // Mass centre
        fVirtAtomXCoord += (trj_crd->xcoor[inp_corr->sele.selatm[ii]-1] * inp_corr->pfMasses[ii]);
        fVirtAtomYCoord += (trj_crd->ycoor[inp_corr->sele.selatm[ii]-1] * inp_corr->pfMasses[ii]);
        fVirtAtomZCoord += (trj_crd->zcoor[inp_corr->sele.selatm[ii]-1] * inp_corr->pfMasses[ii]);
        
        fTotalWeight += inp_corr->pfMasses[ii];
      }
      iResCont++;
    }
    
    // Process the last residue
    if(inp_corr->iMassFlag==0)
    {
      fVirtAtomXCoord = (fVirtAtomXCoord/(float)inp_corr->piSelResLen[iResNum]);
      fVirtAtomYCoord = (fVirtAtomYCoord/(float)inp_corr->piSelResLen[iResNum]);
      fVirtAtomZCoord = (fVirtAtomZCoord/(float)inp_corr->piSelResLen[iResNum]);
    }
    else
    {
      fVirtAtomXCoord = (fVirtAtomXCoord/fTotalWeight);
      fVirtAtomYCoord = (fVirtAtomYCoord/fTotalWeight);
      fVirtAtomZCoord = (fVirtAtomZCoord/fTotalWeight);
      
      fTotalWeight = 0.0;
    }
        
    iVirtAtom++;
    inp_corr->ppfVirtAtomCoord[iVirtAtom][0]=fVirtAtomXCoord;
    inp_corr->ppfVirtAtomCoord[iVirtAtom][1]=fVirtAtomYCoord;
    inp_corr->ppfVirtAtomCoord[iVirtAtom][2]=fVirtAtomZCoord;
  }
  

  if(inp_corr->iCorrType == 0)
  {
    // --TYPE DCC
    if(inp_corr->iMultiAtomFlag==0)
    {
      // One atom per residue or correlation by atoms
      if(trj_crd->pbc_flag == 0)
      {
        // No PBC 
        for(ii=0; ii<inp_corr->iMatrixDim; ii++)
        {
          iAtomNumA=inp_corr->sele.selatm[ii]-1;
          
          fXDiffA=(trj_crd->xcoor[iAtomNumA]-molecule->coor.xcoor[iAtomNumA]);
          fYDiffA=(trj_crd->ycoor[iAtomNumA]-molecule->coor.ycoor[iAtomNumA]);
          fZDiffA=(trj_crd->zcoor[iAtomNumA]-molecule->coor.zcoor[iAtomNumA]);
          
          inp_corr->pfSquareMeanDeviation[ii] += ((fXDiffA*fXDiffA)+(fYDiffA*fYDiffA)+(fZDiffA*fZDiffA));
          
          for(jj=0; jj<inp_corr->iMatrixDim; jj++)
          {
            iAtomNumB=inp_corr->sele.selatm[jj]-1;
            
            fXDiffB=(trj_crd->xcoor[iAtomNumB]-molecule->coor.xcoor[iAtomNumB]);
            fYDiffB=(trj_crd->ycoor[iAtomNumB]-molecule->coor.ycoor[iAtomNumB]);
            fZDiffB=(trj_crd->zcoor[iAtomNumB]-molecule->coor.zcoor[iAtomNumB]);
            
            inp_corr->ppfMeanDeviation[ii][jj] += ((fXDiffA*fXDiffB)+(fYDiffA*fYDiffB)+(fZDiffA*fZDiffB));
          }
        }
      }
      
      else
      {
        // PBC present !
        for(ii=0; ii<inp_corr->iMatrixDim; ii++)
        {
          iAtomNumA=inp_corr->sele.selatm[ii]-1;          
          DistanceAxes(inp_corr->pfAxesDist,
                       trj_crd->xcoor[iAtomNumA],
                       trj_crd->ycoor[iAtomNumA],
                       trj_crd->zcoor[iAtomNumA], 
                       molecule->coor.xcoor[iAtomNumA],
                       molecule->coor.ycoor[iAtomNumA],
                       molecule->coor.zcoor[iAtomNumA],
                       trj_crd->pbc);
          
          fXDiffA=inp_corr->pfAxesDist[0];
          fYDiffA=inp_corr->pfAxesDist[1];
          fZDiffA=inp_corr->pfAxesDist[2];
          
          inp_corr->pfSquareMeanDeviation[ii] += ((fXDiffA*fXDiffA)+(fYDiffA*fYDiffA)+(fZDiffA*fZDiffA));
          
          for(jj=0; jj<inp_corr->iMatrixDim; jj++)
          {
            iAtomNumB=inp_corr->sele.selatm[jj]-1;

            DistanceAxes(inp_corr->pfAxesDist,
                         trj_crd->xcoor[iAtomNumB],
                         trj_crd->ycoor[iAtomNumB],
                         trj_crd->zcoor[iAtomNumB],
                         molecule->coor.xcoor[iAtomNumB],
                         molecule->coor.ycoor[iAtomNumB],
                         molecule->coor.zcoor[iAtomNumB],
                         trj_crd->pbc);
            
            fXDiffB=inp_corr->pfAxesDist[0];
            fYDiffB=inp_corr->pfAxesDist[1];
            fZDiffB=inp_corr->pfAxesDist[2];
            
            inp_corr->ppfMeanDeviation[ii][jj] += ((fXDiffA*fXDiffB)+(fYDiffA*fYDiffB)+(fZDiffA*fZDiffB));
          }
        }
      }
    }
    
    else
    {
      // More than one atom per residue
      if(trj_crd->pbc_flag == 0)
      {
        // No PBC
        for(ii=0; ii<inp_corr->iMatrixDim; ii++)
        {
          fXDiffA=(inp_corr->ppfVirtAtomCoord[ii][0]-inp_corr->ppfVirtRefAtomCoord[ii][0]);
          fYDiffA=(inp_corr->ppfVirtAtomCoord[ii][1]-inp_corr->ppfVirtRefAtomCoord[ii][1]);
          fZDiffA=(inp_corr->ppfVirtAtomCoord[ii][2]-inp_corr->ppfVirtRefAtomCoord[ii][2]);
          
          inp_corr->pfSquareMeanDeviation[ii] += ((fXDiffA*fXDiffA)+(fYDiffA*fYDiffA)+(fZDiffA*fZDiffA));
          
          for(jj=0; jj<inp_corr->iMatrixDim; jj++)
          {
            fXDiffB=(inp_corr->ppfVirtAtomCoord[jj][0]-inp_corr->ppfVirtRefAtomCoord[jj][0]);
            fYDiffB=(inp_corr->ppfVirtAtomCoord[jj][1]-inp_corr->ppfVirtRefAtomCoord[jj][1]);
            fZDiffB=(inp_corr->ppfVirtAtomCoord[jj][2]-inp_corr->ppfVirtRefAtomCoord[jj][2]);
            
            inp_corr->ppfMeanDeviation[ii][jj] += ((fXDiffA*fXDiffB)+(fYDiffA*fYDiffB)+(fZDiffA*fZDiffB));
          }
        }
      }
      
      else
      {
        // PBC Present
        for(ii=0; ii<inp_corr->iMatrixDim; ii++)
        {
          DistanceAxes(inp_corr->pfAxesDist,
                       inp_corr->ppfVirtAtomCoord[ii][0],
                       inp_corr->ppfVirtAtomCoord[ii][1],
                       inp_corr->ppfVirtAtomCoord[ii][2],
                       inp_corr->ppfVirtRefAtomCoord[ii][0],
                       inp_corr->ppfVirtRefAtomCoord[ii][1],
                       inp_corr->ppfVirtRefAtomCoord[ii][2],
                       trj_crd->pbc);
          
          fXDiffA=inp_corr->pfAxesDist[0];
          fYDiffA=inp_corr->pfAxesDist[1];
          fZDiffA=inp_corr->pfAxesDist[2];
          
          inp_corr->pfSquareMeanDeviation[ii] += ((fXDiffA*fXDiffA)+(fYDiffA*fYDiffA)+(fZDiffA*fZDiffA));
          
          for(jj=0; jj<inp_corr->iMatrixDim; jj++)
          {
            DistanceAxes(inp_corr->pfAxesDist,
                         inp_corr->ppfVirtAtomCoord[jj][0],
                         inp_corr->ppfVirtAtomCoord[jj][1],
                         inp_corr->ppfVirtAtomCoord[jj][2],
                         inp_corr->ppfVirtRefAtomCoord[jj][0],
                         inp_corr->ppfVirtRefAtomCoord[jj][1],
                         inp_corr->ppfVirtRefAtomCoord[jj][2],
                         trj_crd->pbc);
            
            fXDiffB=inp_corr->pfAxesDist[0];
            fYDiffB=inp_corr->pfAxesDist[1];
            fZDiffB=inp_corr->pfAxesDist[2];
            
            inp_corr->ppfMeanDeviation[ii][jj] += ((fXDiffA*fXDiffB)+(fYDiffA*fYDiffB)+(fZDiffA*fZDiffB));
          }
        }
      }
    }
  }
  
  else if(inp_corr->iCorrType == 1)
  {
    // --TYPE LMI
    if(inp_corr->iMultiAtomFlag == 0)
    {
      // One atom per residue or correlation by atoms
      if(trj_crd->pbc_flag == 0)
      {
        // No PBC
        for(ii=0; ii<inp_corr->sele.nselatm; ii++)
        {
          iAtomNum = inp_corr->sele.selatm[ii] - 1;
          
          // Set X, Y & Z coordinates for all selected atoms for all processed frames
          inp_corr->pppdCoord[ii][inp_corr->iFrameNum][0] = (trj_crd->xcoor[iAtomNum]-molecule->coor.xcoor[iAtomNum]) / 10.0;
          inp_corr->pppdCoord[ii][inp_corr->iFrameNum][1] = (trj_crd->ycoor[iAtomNum]-molecule->coor.ycoor[iAtomNum]) / 10.0;
          inp_corr->pppdCoord[ii][inp_corr->iFrameNum][2] = (trj_crd->zcoor[iAtomNum]-molecule->coor.zcoor[iAtomNum]) / 10.0;
        }
      }
      
      else
      {
        // PBC Present !
        for(ii=0; ii<inp_corr->sele.nselatm; ii++)
        {
          iAtomNum = inp_corr->sele.selatm[ii] - 1;
          
          // Set X, Y & Z coordinates for all selected atoms for all processed frames
          DistanceAxes(inp_corr->pfAxesDist,
                       trj_crd->xcoor[iAtomNum],
                       trj_crd->ycoor[iAtomNum],
                       trj_crd->zcoor[iAtomNum],
                       molecule->coor.xcoor[iAtomNum],
                       molecule->coor.ycoor[iAtomNum],
                       molecule->coor.zcoor[iAtomNum],
                       trj_crd->pbc);
          
          inp_corr->pppdCoord[ii][inp_corr->iFrameNum][0] = (inp_corr->pfAxesDist[0] / 10.0);
          inp_corr->pppdCoord[ii][inp_corr->iFrameNum][1] = (inp_corr->pfAxesDist[1] / 10.0);
          inp_corr->pppdCoord[ii][inp_corr->iFrameNum][2] = (inp_corr->pfAxesDist[2] / 10.0);
        }
      }
    }
    
    else
    {
      // More than one atom per residue
      
      if(trj_crd->pbc_flag == 0)
      {
        // No PBC
        for(ii=0; ii<inp_corr->iMatrixDim; ii++)
        {
          inp_corr->pppdCoord[ii][inp_corr->iFrameNum][0] = (inp_corr->ppfVirtAtomCoord[ii][0]-inp_corr->ppfVirtRefAtomCoord[ii][0]) / 10.0;
          inp_corr->pppdCoord[ii][inp_corr->iFrameNum][1] = (inp_corr->ppfVirtAtomCoord[ii][1]-inp_corr->ppfVirtRefAtomCoord[ii][1]) / 10.0;
          inp_corr->pppdCoord[ii][inp_corr->iFrameNum][2] = (inp_corr->ppfVirtAtomCoord[ii][2]-inp_corr->ppfVirtRefAtomCoord[ii][2]) / 10.0;
        }
      }
      
      else
      {
        // PBC Present !
        for(ii=0; ii<inp_corr->iMatrixDim; ii++)
        {
          DistanceAxes(inp_corr->pfAxesDist,
                       inp_corr->ppfVirtAtomCoord[ii][0],
                       inp_corr->ppfVirtAtomCoord[ii][1],
                       inp_corr->ppfVirtAtomCoord[ii][2],
                       inp_corr->ppfVirtRefAtomCoord[ii][0],
                       inp_corr->ppfVirtRefAtomCoord[ii][1],
                       inp_corr->ppfVirtRefAtomCoord[ii][2],
                       trj_crd->pbc);
          
          inp_corr->pppdCoord[ii][inp_corr->iFrameNum][0] = (inp_corr->pfAxesDist[0] / 10.0);
          inp_corr->pppdCoord[ii][inp_corr->iFrameNum][1] = (inp_corr->pfAxesDist[1] / 10.0);
          inp_corr->pppdCoord[ii][inp_corr->iFrameNum][2] = (inp_corr->pfAxesDist[2] / 10.0);
        }
      }
    }
  }
  
  sprintf( outstring, "            ");
  return 12;
}

int Post_CORR(struct inp_corr *inp_corr, Molecule *molecule, struct sopt *OPT, CoorSet *trj_crd)
{
  int      ii, jj;
  int      iProgNumA, iProgNumB;
  int      iRowNum, iColNum, iPos;
  int      iMatrixSize, iVectorSize;
  
  double   dCorrValue=0.0;
  double **ppdMatrix;
  
  char     cLabelA[20], cLabelB[20], cResCode1[3];
  
  time_t   time_Today;
  
  // === Some data post-processing =====================================
  if(inp_corr->iCorrType == 0)
  {
    // --TYPE DCC
    
    for(ii=0; ii<inp_corr->iMatrixDim; ii++)
    {
      inp_corr->pfSquareMeanDeviation[ii] = sqrt((inp_corr->pfSquareMeanDeviation[ii]/(float)inp_corr->iNumOfFrames));
      for(jj=0; jj<inp_corr->iMatrixDim; jj++)
      {
        inp_corr->ppfMeanDeviation[ii][jj] = (inp_corr->ppfMeanDeviation[ii][jj]/(float)inp_corr->iNumOfFrames);
      }
    }
  }
  
  else if(inp_corr->iCorrType == 1)
  {
    // --TYPE LMI
    CovarMatrix(inp_corr);
    GaussCorrMatrix(inp_corr);
    
    // Transforms mutual information into generalized correlation coeff
    for(ii=0; ii<inp_corr->iNumOfSelAtoms*inp_corr->iNumOfSelAtoms; ii++)
    {
      if(inp_corr->pdOutput[ii] > 0.0)
        inp_corr->pdOutput[ii] = sqrt(1.0-exp(-2.0/3.0*inp_corr->pdOutput[ii]));
      
      else
        inp_corr->pdOutput[ii] = 0.0;
    }
    
    iVectorSize = inp_corr->iNumOfSelAtoms*inp_corr->iNumOfSelAtoms;
    iMatrixSize = (int) sqrt(iVectorSize);
    
    ppdMatrix = (double **) malloc(iMatrixSize*sizeof(double *));
    for(ii=0; ii<iMatrixSize; ii++)
      ppdMatrix[ii] = (double *) malloc(iMatrixSize*sizeof(double));
    
    iRowNum =  0;
    iColNum = -1;
    iPos    = -1;
    
    for(ii=0; ii<iVectorSize; ii++)
    {
      iPos++;
      if(iPos == iMatrixSize)
      {
        iRowNum++;
        iColNum = 0;
        iPos    = 0;
      }
      
      else
        iColNum++;
      
      ppdMatrix[iRowNum][iColNum] = inp_corr->pdOutput[ii];
    }
  }
  // ===================================================================
  
  // === It's time for results !!! =====================================
  time(&time_Today);
  fprintf(inp_corr->FOutFile, "# ==============================================================\n");
  fprintf(inp_corr->FOutFile, "# ***                   WORDOM CORR MODULE                   ***\n");
  fprintf(inp_corr->FOutFile, "# ==============================================================\n");
  fprintf(inp_corr->FOutFile, "#\n");
  fprintf(inp_corr->FOutFile, "# Version   : %s\n", VERSTRING);
  fprintf(inp_corr->FOutFile, "# License   : GPL 3\n");
  fprintf(inp_corr->FOutFile, "# Copyright : Fanelli, Felline\n");
  fprintf(inp_corr->FOutFile, "#             University of Modena\n");
  fprintf(inp_corr->FOutFile, "#             Modena - Italy\n");
  fprintf(inp_corr->FOutFile, "#\n");
  //fprintf(inp_corr->FOutFile, "# Date      : %s", asctime(localtime(&time_Today)));
  //fprintf(inp_corr->FOutFile, "#\n");
  //fprintf(inp_corr->FOutFile, "# Mol File  : %s\n", OPT->IMOL_FILE);
  fprintf(inp_corr->FOutFile, "# Res Num   : %d\n", inp_corr->iNumOfRes);
  //fprintf(inp_corr->FOutFile, "# Traj File : %s\n", OPT->ITRJ_FILE);
  fprintf(inp_corr->FOutFile, "# Frame Num : %d\n", inp_corr->iNumOfFrames);
  
  if(trj_crd->pbc_flag == 0)
    fprintf(inp_corr->FOutFile, "# PBC       : No\n");
  else
    fprintf(inp_corr->FOutFile, "# PBC       : Yes\n");
  
  fprintf(inp_corr->FOutFile, "#\n");
  fprintf(inp_corr->FOutFile, "# Title     : %s\n", inp_corr->cTitle);
  
  if(inp_corr->iCorrType == 0)
    fprintf(inp_corr->FOutFile, "# Type      : DCC\n");
  else if(inp_corr->iCorrType == 1)
    fprintf(inp_corr->FOutFile, "# Type      : LMI\n");
  
  fprintf(inp_corr->FOutFile, "# Sele      : %s\n", inp_corr->sele.selestring);
  
  if(inp_corr->iResFlag==1)
    fprintf(inp_corr->FOutFile, "# Sele Res  : %d\n", inp_corr->iMatrixDim);
  else
    fprintf(inp_corr->FOutFile, "# Sele Atm  : %d\n", inp_corr->iMatrixDim);
    
  if(inp_corr->iResFlag==1)
    fprintf(inp_corr->FOutFile, "# Level     : RES\n");
  else
    fprintf(inp_corr->FOutFile, "# Level     : ATM\n");
  
  if(inp_corr->iMassFlag==1)
    fprintf(inp_corr->FOutFile, "# Mass      : Yes\n");
  else
    fprintf(inp_corr->FOutFile, "# Mass      : No\n");
  
  fprintf(inp_corr->FOutFile, "#\n");
  
  if(inp_corr->iMultiAtomFlag == 1 && trj_crd->pbc_flag == 1)
    fprintf(inp_corr->FOutFile, "# Warning!  : PBC not taken into account in the calculation of geo/mass centers\n#\n");
    
  if(inp_corr->iZeroMassWarning == 1)
    fprintf(inp_corr->FOutFile, "# Warning!  : Atom(s) with zero mass\n#\n");
  
  fprintf(inp_corr->FOutFile, "# ==============================================================\n");
  
  fprintf(inp_corr->FOutFile, "#%8s   %8s   %15s   %15s   %5s\n", "ProgNumA", "ProgNumB", "LabelA", "LabelB", "CORR");
  
  for(ii=0; ii<inp_corr->iMatrixDim; ii++)
  {
    for(jj=0; jj<inp_corr->iMatrixDim; jj++)
    {
      if(inp_corr->iResFlag==1)
      {
        if(inp_corr->iMultiAtomFlag==0)
        {
          Res3ToRes1(molecule->rawmol.restype[inp_corr->sele.selatm[ii]-1], cResCode1);
          sprintf(cLabelA, "%s:%s%d", molecule->rawmol.segId[inp_corr->sele.selatm[ii]-1], cResCode1,
                                      molecule->rawmol.resn[inp_corr->sele.selatm[ii]-1]);
        
          Res3ToRes1(molecule->rawmol.restype[inp_corr->sele.selatm[jj]-1], cResCode1);
          sprintf(cLabelB, "%s:%s%d", molecule->rawmol.segId[inp_corr->sele.selatm[jj]-1], cResCode1,
                                      molecule->rawmol.resn[inp_corr->sele.selatm[jj]-1]);
          
          iProgNumA=molecule->rawmol.presn[inp_corr->sele.selatm[ii]-1];
          iProgNumB=molecule->rawmol.presn[inp_corr->sele.selatm[jj]-1];
        }
        
        else
        {
          strcpy(cLabelA, inp_corr->pcLabels[ii]);
          strcpy(cLabelB, inp_corr->pcLabels[jj]);
          iProgNumA=inp_corr->piProgNum[ii];
          iProgNumB=inp_corr->piProgNum[jj];
        }
        
      }
      else
      {
        Res3ToRes1(molecule->rawmol.restype[inp_corr->sele.selatm[ii]-1], cResCode1);
        sprintf(cLabelA, "%s:%s%d:%s", molecule->rawmol.segId[inp_corr->sele.selatm[ii]-1], cResCode1,
                                       molecule->rawmol.resn[inp_corr->sele.selatm[ii]-1],
                                       molecule->rawmol.atmtype[inp_corr->sele.selatm[ii]-1]);
        
        Res3ToRes1(molecule->rawmol.restype[inp_corr->sele.selatm[jj]-1], cResCode1);
        sprintf(cLabelB, "%s:%s%d:%s", molecule->rawmol.segId[inp_corr->sele.selatm[jj]-1], cResCode1,
                                       molecule->rawmol.resn[inp_corr->sele.selatm[jj]-1],
                                       molecule->rawmol.atmtype[inp_corr->sele.selatm[jj]-1]);

        iProgNumA=atoi(molecule->rawmol.atmId[inp_corr->sele.selatm[ii]-1]);
        iProgNumB=atoi(molecule->rawmol.atmId[inp_corr->sele.selatm[jj]-1]);
	
      }
      
      if(inp_corr->iCorrType == 0)
        dCorrValue = (inp_corr->ppfMeanDeviation[ii][jj]/(inp_corr->pfSquareMeanDeviation[ii]*inp_corr->pfSquareMeanDeviation[jj]));
        
      else if(inp_corr->iCorrType == 1)
        dCorrValue = ppdMatrix[ii][jj];
      
      fprintf(inp_corr->FOutFile, " %8d   %8d   %15s   %15s   %5.2f\n", iProgNumA, iProgNumB, cLabelA, cLabelB, dCorrValue);
    }
  }
  
  return 0;
  // ===================================================================
}
