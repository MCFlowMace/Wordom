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

    FLUCT
    * ???
*/


#include "wordom.h"
#include "tools.h"
#include "fileio.h"
#include "datahandler.h"
#include "time.h"
#include "corr.h"

#define  VERSTRING "0.4a"
#define  PI        3.14159265358979323846

// a pointer to inp_corr, used in multithreadning
struct inp_corr  *pInpCorr;
                  
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
  int     ii, jj, kk, mm, xx, yy;
  int     iWinpOptFlag=0;
  int     iResNum, iResCont, iVirtAtom;
  int     iInputBegIndex;
  int     iSeleIdx;
  int     iMoleculeItems;
  
  char    cWordomInpBuffer[1024], cTmpString[99];
  char    cLastRes[50], cOutPutFileName[512];
  char    cResCode1[3];
  
  float   fTotalWeight;
  float   fVirtAtomXCoord, fVirtAtomYCoord, fVirtAtomZCoord;

  extern short int     no_frame_par;
  no_frame_par = 1;
  
  memset ( cWordomInpBuffer, '\0', sizeof(cWordomInpBuffer));
  
  // === example winp file ===
  //BEGIN CORR
  //--TITLE CORR
  //--TYPE DCC
  //--SELE /*/*/*
  //--LEVEL RES
  //--MASS 0
  //END
  // =========================
  
  // === Default Values ================================================
  inp_corr->iPDBFlag         =  0;
  sprintf(inp_corr->sele.selestring, "/*/*/CA");
  inp_corr->iMassFlag        =  0;
  inp_corr->iResFlag         =  1;
  inp_corr->iCorrType        =  0;
  inp_corr->iNumOfThreads    =  0;
  inp_corr->iFirstRoundFlag  =  0;
  inp_corr->iVerboseFlag     =  0;
  inp_corr->iMatchSubSele    =  0;
  inp_corr->iNumOfSubMatch   =  0;
  inp_corr->iStdDevFlag      =  0;
  inp_corr->iGetFramesFlag   = -1;
  // ===================================================================
  
  iInputBegIndex = input_index;
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
    else if ( !strncmp(cWordomInpBuffer, "--VERBOSE", 9))
    {
      inp_corr->iVerboseFlag = 1;
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--STDEV", 7))
    {
      inp_corr->iStdDevFlag = 1;
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--MATCHSUBSELE", 10))
    {
      inp_corr->iMatchSubSele = 1;
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--TYPE", 6))
    {
      sscanf( cWordomInpBuffer, "--TYPE %s", cTmpString);
      if(strcmp(cTmpString, "DCC")==0)
        inp_corr->iCorrType=0;
      else if(strcmp(cTmpString, "LMI")==0)
        inp_corr->iCorrType=1;
      else if(strcmp(cTmpString, "DCOR")==0)
        inp_corr->iCorrType=2;
      else if(strcmp(cTmpString, "FLUCT")==0)
        inp_corr->iCorrType=3;
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
    else if ( !strncmp(cWordomInpBuffer, "--NT ", 5))
    {
      sscanf(cWordomInpBuffer, "--NT %d", &inp_corr->iNumOfThreads);
      if(inp_corr->iNumOfThreads<0)
      {
        fprintf( stderr, "CORR module: --NT option needs a number >= 0, passed %d\n", inp_corr->iNumOfThreads);
        exit(5);
      }

      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--SELE", 6))
    {
      sscanf(cWordomInpBuffer, "--SELE %[^\n]%*c ", inp_corr->sele.selestring);
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--SUBSELE", 9))
    {
      inp_corr->iNumOFSubSele+=1;
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--GETFRAMES ", 12))
    {
      sscanf( cWordomInpBuffer, "--GETFRAMES %s", cTmpString);
      
      if(strcmp(cTmpString, "YES")==0)
        inp_corr->iGetFramesFlag = 1;
      else if(strcmp(cTmpString, "Yes")==0)
        inp_corr->iGetFramesFlag = 1;
      else if(strcmp(cTmpString, "yes")==0)
        inp_corr->iGetFramesFlag = 1;
      else if(strcmp(cTmpString, "1")==0)
        inp_corr->iGetFramesFlag = 1;
      else if(strcmp(cTmpString, "NO")==0)
        inp_corr->iGetFramesFlag = 0;
      else if(strcmp(cTmpString, "No")==0)
        inp_corr->iGetFramesFlag = 0;
      else if(strcmp(cTmpString, "no")==0)
        inp_corr->iGetFramesFlag = 0;
      else if(strcmp(cTmpString, "0")==0)
        inp_corr->iGetFramesFlag = 0;
      else
      {
        fprintf( stderr, "CORR module: invalid value passed to --GETFRAMES option: %s\n", cTmpString);
        exit(5);
      }
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
  
  if(inp_corr->iNumOfThreads > 0 && inp_corr->iCorrType != 2)
    fprintf( stderr, "CORR module: Sorry, multi-threading is only available in DCor calculation, this option will be ignored\n");
  
  
  if(inp_corr->iStdDevFlag == 1 && inp_corr->iCorrType != 3)
    fprintf( stderr, "CORR module: --STDEV option is only available with --TYPE FLUCT, this option will be ignored\n");
  
  if(inp_corr->iNumOFSubSele > 0 && inp_corr->iCorrType != 3)
    fprintf( stderr, "CORR module: --SUBSELE option is only available with --TYPE FLUCT, this option will be ignored\n");
  
  if(inp_corr->pfMatchSubSeleFluct == 1)
  { 
    if(inp_corr->iCorrType != 3)
      fprintf( stderr, "CORR module: --MATCHSUBSELE option is only available with --TYPE FLUCT, this option will be ignored\n");
      
    if(inp_corr->iNumOFSubSele < 2)
      fprintf( stderr, "CORR module: --MATCHSUBSELE is only available if you define at least 2 sub-selections, this option will be ignored\n");
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
  
  if(inp_corr->iNumOfThreads == 0 && inp_corr->iCorrType)
    inp_corr->iNumOfThreads = 1;
  
  if(inp_corr->iCorrType != 3 && inp_corr->iGetFramesFlag != -1)
    fprintf( stderr, "CORR module: --GETFRAMES option is aveilable only when --TYPE is FLUCT, it will be ignored\n");
  
  if(inp_corr->iCorrType == 3 && inp_corr->iGetFramesFlag == -1)
    inp_corr->iGetFramesFlag = 0;
  
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
  
  else if(inp_corr->iCorrType == 2)
  {
    // --TYPE DCOR *-*

    inp_corr->iFrameNum=-1;
    // allocated threads structure
    inp_corr->corr_threads = (pthread_t *) calloc(inp_corr->iNumOfThreads, sizeof(pthread_t));
    
    // number of frames per threads
    inp_corr->iNumOfFrmPerThreads = (int) (ceil(iNumOfFrames / inp_corr->iNumOfThreads) + 1);

    //inp_corr->iNumOfPairsPerThreads = (int) (((iNumOfFrames * iNumOfFrames) + iNumOfFrames) / 2.0);
    inp_corr->iNumOfPairsPerThreads = (int) (iNumOfFrames * iNumOfFrames);
    inp_corr->iNumOfPairsPerThreads = (int) (ceil(inp_corr->iNumOfPairsPerThreads / inp_corr->iNumOfThreads) + 1);

    //inp_corr->ppfThreadResOverallAvgDist = (double **) calloc(inp_corr->iNumOfThreads, sizeof(double *));
    //for(ii=0; ii<inp_corr->iNumOfThreads; ++ii)
    //  inp_corr->ppfThreadResOverallAvgDist[ii] = (double *) calloc(inp_corr->iMatrixDim, sizeof(double));
    //
    //inp_corr->pppfThreadResFrameAvgDist = (double ***) calloc(inp_corr->iNumOfThreads, sizeof(double **));
    //for(ii=0; ii<inp_corr->iNumOfThreads; ++ii)
    //{
    //  inp_corr->pppfThreadResFrameAvgDist[ii] = (double **) calloc(inp_corr->iMatrixDim, sizeof(double *));
    //  for(jj=0; jj<inp_corr->iMatrixDim; jj++)
    //    inp_corr->pppfThreadResFrameAvgDist[ii][jj] = (double *) calloc(inp_corr->iNumOfPairsPerThreads, sizeof(double));
    //}

    inp_corr->pfResOverallAvgDist = (double *) calloc(inp_corr->iMatrixDim, sizeof(double));
    inp_corr->ppfResFrameAvgDist = (double **) calloc(inp_corr->iMatrixDim, sizeof(double *));
    for(ii=0; ii<inp_corr->iMatrixDim; ii++)
      inp_corr->ppfResFrameAvgDist[ii] = (double *) calloc(inp_corr->iNumOfFrames, sizeof(double));
    
    inp_corr->pppdCoord = (double ***) calloc(iNumOfFrames, sizeof(double **));
    for(ii=0; ii<iNumOfFrames; ii++)
    {
      inp_corr->pppdCoord[ii] = (double **) calloc(inp_corr->iMatrixDim, sizeof(double *));
      for(jj=0; jj<inp_corr->iMatrixDim; jj++)
      {
        inp_corr->pppdCoord[ii][jj] = (double *) calloc(3, sizeof(double));
      }
    }

    inp_corr->pppfThreadDistCov = (double ***) calloc(inp_corr->iNumOfThreads, sizeof(double **));
    for(ii=0; ii<inp_corr->iNumOfThreads; ++ii)
    {
      inp_corr->pppfThreadDistCov[ii] = (double **) calloc(inp_corr->iMatrixDim, sizeof(double *));
      for(jj=0; jj<inp_corr->iMatrixDim; ++jj)
        inp_corr->pppfThreadDistCov[ii][jj] = (double *) calloc(inp_corr->iMatrixDim, sizeof(double));
    }

    inp_corr->ppfDCorr   = (double **) calloc(inp_corr->iMatrixDim, sizeof(double *));
    inp_corr->ppfDistCov = (double **) calloc(inp_corr->iMatrixDim, sizeof(double *));
    for(ii=0; ii<inp_corr->iMatrixDim; ii++)
    {
      inp_corr->ppfDCorr[ii]   = (double *) calloc(inp_corr->iMatrixDim, sizeof(double));
      inp_corr->ppfDistCov[ii] = (double *) calloc(inp_corr->iMatrixDim, sizeof(double));
    }


    
    // ppiFrmListPerThreads holds the list of frames that each threds will process
    // after allocation all position is set to -1
    inp_corr->ppiFrmListPerThreads = (int **) calloc(inp_corr->iNumOfThreads, sizeof(int *));
    for(ii=0; ii<inp_corr->iNumOfThreads; ++ii)
    {
      inp_corr->ppiFrmListPerThreads[ii] = (int *) calloc(inp_corr->iNumOfFrmPerThreads, sizeof(int));
      for(jj=0; jj<inp_corr->iNumOfFrmPerThreads; ++jj)
        inp_corr->ppiFrmListPerThreads[ii][jj] = -1;
    }
    
    // copy the list of franes that each threads will process
    // the last list will have a list of trailing -1
    
    jj =  0; // thread id
    kk = -1; // obj id
    for(ii=0; ii<iNumOfFrames; ++ii)
    {
      kk += 1;
      if(kk == inp_corr->iNumOfFrmPerThreads)
      {
        kk = 0;
        jj += 1;
      }
      inp_corr->ppiFrmListPerThreads[jj][kk] = ii;
    }

    //printf("\n");
    //printf("THR  NUM  FRM\n");
    //for(ii=0; ii<inp_corr->iNumOfThreads; ++ii)
    //  for(jj=0; jj<inp_corr->iNumOfFrmPerThreads; ++jj)
    //    printf("%3d  %3d  %3d\n", ii, jj, inp_corr->ppiFrmListPerThreads[ii][jj]);
    //exit(5);
    
    // inp_corr->ppiPairsListFrm1 and inp_corr->ppiPairsListFrm2 hold the first and the last
    // frame of each pair. All items of both arrays are set to -1
    inp_corr->ppiPairsListFrm1 = (int **) calloc(inp_corr->iNumOfThreads, sizeof(int *));
    inp_corr->ppiPairsListFrm2 = (int **) calloc(inp_corr->iNumOfThreads, sizeof(int *));
    for(ii=0; ii<inp_corr->iNumOfThreads; ++ii)
    {
      inp_corr->ppiPairsListFrm1[ii] = (int *) calloc(inp_corr->iNumOfPairsPerThreads, sizeof(int));
      inp_corr->ppiPairsListFrm2[ii] = (int *) calloc(inp_corr->iNumOfPairsPerThreads, sizeof(int));
      for(jj=0; jj<inp_corr->iNumOfPairsPerThreads; ++jj)
      {
        inp_corr->ppiPairsListFrm1[ii][jj] = -1;
        inp_corr->ppiPairsListFrm2[ii][jj] = -1;
      }
    }

    // populates ppiPairsListFrm1 and ppiPairsListFrm2 with all frame pairs
    kk =  0; // thread id
    mm = -1; // pair id
    for(ii=0; ii<iNumOfFrames; ++ii)
    {
      for(jj=0; jj<iNumOfFrames; ++jj)
      {
        mm += 1;
        if(mm == inp_corr->iNumOfPairsPerThreads)
        {
          mm = 0;
          kk += 1;
        }
        inp_corr->ppiPairsListFrm1[kk][mm] = ii;
        inp_corr->ppiPairsListFrm2[kk][mm] = jj;
      }
    }
    
    //printf("\n");
    //printf("THR  NUM  FN1  FN2\n");
    //for(ii=0; ii<inp_corr->iNumOfThreads; ++ii)
    //  for(jj=0; jj<inp_corr->iNumOfPairsPerThreads; ++jj)
    //    printf("%3d  %3d  %3d  %3d\n", ii, jj, inp_corr->ppiPairsListFrm1[ii][jj], inp_corr->ppiPairsListFrm2[ii][jj]);
    //exit(0);
    
    
    // copy inp_corr pointer to global pInpCorr
    pInpCorr = inp_corr;
  }
  
  else if(inp_corr->iCorrType == 3)
  {

    // --TYPE FLUCT

    inp_corr->iFirstRoundFlag = 1; 
    inp_corr->iFirstFrame     = 1;

    inp_corr->ppfMeanDistance  = (double **) calloc(inp_corr->iMatrixDim, sizeof(double *));
    inp_corr->ppfFluctMatrix   = (double **) calloc(inp_corr->iMatrixDim, sizeof(double *));
    inp_corr->ppfMinDistances  = (double **) calloc(inp_corr->iMatrixDim, sizeof(double *));
    inp_corr->ppfMaxDistances  = (double **) calloc(inp_corr->iMatrixDim, sizeof(double *));
    for(ii=0; ii<inp_corr->iMatrixDim; ++ii)
    {
      inp_corr->ppfMeanDistance[ii] = (double  *) calloc(inp_corr->iMatrixDim, sizeof(double));
      inp_corr->ppfFluctMatrix[ii]  = (double  *) calloc(inp_corr->iMatrixDim, sizeof(double));
      inp_corr->ppfMinDistances[ii] = (double  *) calloc(inp_corr->iMatrixDim, sizeof(double));
      inp_corr->ppfMaxDistances[ii] = (double  *) calloc(inp_corr->iMatrixDim, sizeof(double));
    }

    if(inp_corr->iGetFramesFlag == 1)
    {
      inp_corr->ppfFrameMatrix  = (double **) calloc(inp_corr->iMatrixDim, sizeof(double *));
      for(ii=0; ii<inp_corr->iMatrixDim; ++ii)
        inp_corr->ppfFrameMatrix[ii]  = (double  *) calloc(inp_corr->iMatrixDim, sizeof(double));
    }

    if(inp_corr->iNumOFSubSele > 0)
    {

      if(inp_corr->iMultiAtomFlag==1)
      {
        fprintf( stderr, "CORR module: sorry, but sub-selections with --LEVEL RES and more than 1 atom per residue is not yet implemented\n");
        exit(5);
      }
            
      // this vecotor maps of master-sele res/atm indexes to molecule indexes
      
      if(inp_corr->iResFlag==1)
        iMoleculeItems = molecule->nRes;
      else
        iMoleculeItems = molecule->nato;
      
      inp_corr->piMasterSeleIndexes = (int *) calloc(iMoleculeItems, sizeof(int));
      
      // set all position to -1, which means 'not selected in master-sele'
      
      for(ii=0; ii<iMoleculeItems; ++ii)
        inp_corr->piMasterSeleIndexes[ii] = -1;

      for(ii=0; ii<inp_corr->iMatrixDim; ii++)
      {
        if(inp_corr->iResFlag==1)
        {
          if(inp_corr->iMultiAtomFlag==0)
            iSeleIdx=molecule->rawmol.presn[inp_corr->sele.selatm[ii]-1];
          
          else
            iSeleIdx=inp_corr->piProgNum[ii];
        }
        
        else
          iSeleIdx=atoi(molecule->rawmol.atmId[inp_corr->sele.selatm[ii]-1]);
        
        // so, now we know that the iSeleIdx-th res/atm in passed molecule is the ii-th res/atm in master-sele
        inp_corr->piMasterSeleIndexes[iSeleIdx] = ii;
      }

      inp_corr->pSubSele         = (Selection *) calloc(inp_corr->iNumOFSubSele, sizeof(Selection));
      inp_corr->pfSubSeleFluct   = (double    *) calloc(inp_corr->iNumOFSubSele, sizeof(double));
      inp_corr->pfSubSeleAvgDist = (double    *) calloc(inp_corr->iNumOFSubSele, sizeof(double));
      kk = -1;
      input_index = iInputBegIndex;
      memset(cWordomInpBuffer, '\0', sizeof(cWordomInpBuffer));
      while( strncmp (cWordomInpBuffer, "END", 3))
      {
        sprintf(cWordomInpBuffer, "%s", input[input_index]);
        if ( !strncmp(cWordomInpBuffer, "--SUBSELE", 9))
        {
          // @@@ HERE

          kk += 1;
          sscanf(cWordomInpBuffer, "--SUBSELE %[^\n]%*c ", inp_corr->pSubSele[kk].selestring);
          GetSele(inp_corr->pSubSele[kk].selestring, &inp_corr->pSubSele[kk], molecule);
          if(inp_corr->pSubSele[kk].nselatm == 0)
          {
            fprintf( stderr, "CORR module: sub-selection #%d: %s, is empty\n", kk+1, inp_corr->pSubSele[kk].selestring);
            exit(5);
          }
          
          for(ii=0; ii<inp_corr->pSubSele[kk].nselatm; ++ii)
          {
            if(inp_corr->iResFlag==1)
              iSeleIdx=molecule->rawmol.presn[inp_corr->pSubSele[kk].selatm[ii]-1];
            else
              iSeleIdx=atoi(molecule->rawmol.atmId[inp_corr->pSubSele[kk].selatm[ii]-1]);
            
            if(inp_corr->piMasterSeleIndexes[iSeleIdx] == -1)
            {
              
              if(inp_corr->iResFlag == 1)
              {
                Res3ToRes1(molecule->rawmol.restype[inp_corr->pSubSele[kk].selatm[0]-1], cResCode1);
                sprintf(cTmpString, "%s:%s%d", molecule->rawmol.segId[inp_corr->pSubSele[kk].selatm[0]-1], cResCode1, molecule->rawmol.resn[inp_corr->pSubSele[kk].selatm[0]-1]);
              }
              
              else
              {
                Res3ToRes1(molecule->rawmol.restype[inp_corr->sele.selatm[ii]-1], cResCode1);
                sprintf(cTmpString, "%s:%s%d:%s", molecule->rawmol.segId[inp_corr->sele.selatm[ii]-1], cResCode1,
                                                  molecule->rawmol.resn[inp_corr->sele.selatm[ii]-1],
                                                  molecule->rawmol.atmtype[inp_corr->sele.selatm[ii]-1]);
              }
              
              fprintf( stderr, "CORR module: atm/res %s selected with sub-selection #%d: %s, must be also selected in master-sele\n", cTmpString, kk+1, inp_corr->pSubSele[kk].selestring);
              exit(5);
            }
          }
          
        }
        input_index++;
      }
      
      inp_corr->iNumOfSubMatch        = (int) ((inp_corr->iNumOFSubSele * inp_corr->iNumOFSubSele) - inp_corr->iNumOFSubSele) / 2.0;
      inp_corr->pfMatchSubSeleFluct   = (double *) calloc(inp_corr->iNumOfSubMatch, sizeof(double));
      inp_corr->pfMatchSubSeleAvgDist = (double *) calloc(inp_corr->iNumOfSubMatch, sizeof(double));
    }
    
    if(inp_corr->iVerboseFlag == 1)
      printf("{Corr::Fluct} End of Read_CORR\n");
      
  }

  if(inp_corr->iGetFramesFlag == 1)
  {
    
    sprintf(cOutPutFileName, "%s_frames.dat", inp_corr->cTitle);
    inp_corr->FVerbOutFile=fopen(cOutPutFileName, "w");
    
    fprintf(inp_corr->FVerbOutFile, "%12s   ", "nFr");
    fprintf(inp_corr->FVerbOutFile, "%12s   ", "MasterSele");
    if(inp_corr->iNumOFSubSele > 0)
    {
      for(xx=0; xx<inp_corr->iNumOFSubSele; ++xx)
      {
        sprintf(cTmpString, "Sele%d", xx+1);
        fprintf(inp_corr->FVerbOutFile, "%12s   ", cTmpString);
      }
    }
    
    if(inp_corr->iNumOfSubMatch > 0)
    {
      mm = 0;
      for(xx=0; xx<inp_corr->iNumOFSubSele; ++xx)
      {
        for(yy=xx+1; yy<inp_corr->iNumOFSubSele; ++yy)
        {
          mm += 1;
          sprintf(cTmpString, "Match%d", mm);
          fprintf(inp_corr->FVerbOutFile, "%12s   ", cTmpString);
          
        }
      }
    }
    fprintf(inp_corr->FVerbOutFile, "\n");
  }

  inp_corr->iFrameNum = 0;
  sprintf( outstring, " %10s ", inp_corr->cTitle);
  return 12;
}

int Compute_CORR(struct inp_corr *inp_corr, Molecule *molecule, CoorSet *trj_crd, char *outstring)
{
  int     ii, jj, mm, xx, yy;
  int     iAtomNumA, iAtomNumB;
  int     iResNum, iResCont, iVirtAtom;
  int     iAtomNum;
  int     iTotPairs, iIndexA, iIndexB, iProgNumA, iProgNumB;
  
  float   fXDiffA, fYDiffA, fZDiffA;
  float   fXDiffB, fYDiffB, fZDiffB;
  float   fPairDist, fThisDist;
  float   fVirtAtomXCoord, fVirtAtomYCoord, fVirtAtomZCoord;
  float   fTotalWeight;
  
  inp_corr->iFrameNum++;                                                // Update frame number
  inp_corr->iPDBFlag = trj_crd->pbc_flag;
  
  inp_corr->pPBCInfo = trj_crd->pbc;
  
  if(inp_corr->iFirstFrame == 1)
    inp_corr->iFirstFrame = 0;
  
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
  
  else if(inp_corr->iCorrType == 2)
  {
    // --TYPE DCOR *-*
    if(inp_corr->iMultiAtomFlag == 0)
    {
      // One atom per residue or correlation by atoms
      for(ii=0; ii<inp_corr->sele.nselatm; ii++)
      {
        iAtomNum = inp_corr->sele.selatm[ii] - 1;
        
        inp_corr->pppdCoord[inp_corr->iFrameNum][ii][0] = trj_crd->xcoor[iAtomNum];
        inp_corr->pppdCoord[inp_corr->iFrameNum][ii][1] = trj_crd->ycoor[iAtomNum];
        inp_corr->pppdCoord[inp_corr->iFrameNum][ii][2] = trj_crd->zcoor[iAtomNum];
      }
    }
    
    else
    {
      // More than one atom per residue
      
      for(ii=0; ii<inp_corr->iMatrixDim; ii++)
      {
        inp_corr->pppdCoord[inp_corr->iFrameNum][ii][0] = inp_corr->ppfVirtAtomCoord[ii][0];
        inp_corr->pppdCoord[inp_corr->iFrameNum][ii][1] = inp_corr->ppfVirtAtomCoord[ii][1];
        inp_corr->pppdCoord[inp_corr->iFrameNum][ii][2] = inp_corr->ppfVirtAtomCoord[ii][2];
      }
    }
  }

  else if(inp_corr->iCorrType == 3)
  {
    // --TYPE FLUCT @@@
    
    if(inp_corr->iGetFramesFlag == 1 && inp_corr->iFirstRoundFlag == 0)
    {
      // initializes inp_corr->ppfFrameMatrix
      for(ii=0; ii<inp_corr->iMatrixDim; ++ii)
      {
        inp_corr->ppfFrameMatrix[ii][ii] = 0.0;
        for(jj=ii+1; jj<inp_corr->iMatrixDim; ++jj)
        {
          inp_corr->ppfFrameMatrix[ii][jj] = 0.0;
          inp_corr->ppfFrameMatrix[jj][ii] = 0.0;
        }
      }
    }
    
    for(ii=0; ii<inp_corr->sele.nselatm; ++ii)
    {
      for(jj=ii+1; jj<inp_corr->sele.nselatm; ++jj)
      {
        if(inp_corr->iMultiAtomFlag == 0)
        {
          // One atom per residue or correlation by atoms
          
          iAtomNumA = inp_corr->sele.selatm[ii] - 1;
          iAtomNumB = inp_corr->sele.selatm[jj] - 1;
          fThisDist = DistanceCoor(trj_crd->xcoor[iAtomNumA],
                                   trj_crd->ycoor[iAtomNumA],
                                   trj_crd->zcoor[iAtomNumA],
                                   trj_crd->xcoor[iAtomNumB],
                                   trj_crd->ycoor[iAtomNumB],
                                   trj_crd->zcoor[iAtomNumB],
                                   trj_crd->pbc);
        }
        
        else
        {
          // More than one atom per residue
          fThisDist = DistanceCoor(inp_corr->ppfVirtAtomCoord[ii][0],
                                   inp_corr->ppfVirtAtomCoord[ii][1],
                                   inp_corr->ppfVirtAtomCoord[ii][2],
                                   inp_corr->ppfVirtAtomCoord[jj][0],
                                   inp_corr->ppfVirtAtomCoord[jj][1],
                                   inp_corr->ppfVirtAtomCoord[jj][2],
                                   trj_crd->pbc);
        }

        if(inp_corr->iFirstRoundFlag == 1)
        {
          inp_corr->ppfMeanDistance[ii][jj] += fThisDist;
          if(inp_corr->iFirstFrame == 1)
          {
            inp_corr->ppfMinDistances[ii][jj] = fThisDist;
            inp_corr->ppfMaxDistances[ii][jj] = fThisDist;
          }
          
          else
          {
            if(fThisDist < inp_corr->ppfMinDistances[ii][jj])
              inp_corr->ppfMinDistances[ii][jj] = fThisDist;
            if(fThisDist > inp_corr->ppfMaxDistances[ii][jj])
              inp_corr->ppfMaxDistances[ii][jj] = fThisDist;
          }
        }

        else
        {
          if(inp_corr->iGetFramesFlag == 1)
            inp_corr->ppfFrameMatrix[ii][jj] = fThisDist;

          fThisDist = (fThisDist - inp_corr->ppfMeanDistance[ii][jj]);
          fThisDist = fThisDist * fThisDist;
          inp_corr->ppfFluctMatrix[ii][jj] += fThisDist;

        }
      }
    }
    
    if(inp_corr->iGetFramesFlag == 1 && inp_corr->iFirstRoundFlag == 0)
    {
      // calculates the overall fluctuation over master selection
      fprintf(inp_corr->FVerbOutFile, "%12d   ", inp_corr->iFrameNum);
      fThisDist = 0.0;
      iTotPairs = 0;
      for(ii=0; ii<inp_corr->iMatrixDim; ++ii)
      {
        for(jj=ii+1; jj<inp_corr->iMatrixDim; ++jj)
        {
          iTotPairs += 1;
          fThisDist += inp_corr->ppfFrameMatrix[ii][jj];
        }
      }
      fThisDist = fThisDist / (float) iTotPairs;
      fprintf(inp_corr->FVerbOutFile, "%12.3f   ", fThisDist);
      
      if(inp_corr->iNumOFSubSele > 0)
      {
        // calculates the sub-sele fluctuations
        for(xx=0; xx<inp_corr->iNumOFSubSele; ++xx)
        {
          iTotPairs = 0;
          fThisDist = 0.0;
          for(ii=0; ii<inp_corr->pSubSele[xx].nselatm; ++ii)
          {
            for(jj=ii+1; jj<inp_corr->pSubSele[xx].nselatm; ++jj)
            {
              if(inp_corr->iResFlag==1)
              {
                iIndexA = molecule->rawmol.presn[inp_corr->pSubSele[xx].selatm[ii]-1];
                iIndexB = molecule->rawmol.presn[inp_corr->pSubSele[xx].selatm[jj]-1];
              }
              
              else
              {
                iIndexA = atoi(molecule->rawmol.atmId[inp_corr->pSubSele[xx].selatm[ii]-1]);
                iIndexB = atoi(molecule->rawmol.atmId[inp_corr->pSubSele[xx].selatm[jj]-1]);
              }
              
              iProgNumA = inp_corr->piMasterSeleIndexes[iIndexA];
              iProgNumB = inp_corr->piMasterSeleIndexes[iIndexB];
              
              fThisDist += inp_corr->ppfFrameMatrix[iProgNumA][iProgNumB];
              iTotPairs += 1;
            }
          }
          
          fThisDist = fThisDist / (float) iTotPairs;
          fprintf(inp_corr->FVerbOutFile, "%12.3f   ", fThisDist);
        }
        
        
        if(inp_corr->iNumOfSubMatch > 0)
        {
          // calculates matcthed sub-sele fluctuations
          // i.e. the overall-fluctuations between all res/atm of one sub-sele vs the res/atm of onother sele
          mm = -1;
          for(xx=0; xx<inp_corr->iNumOFSubSele; ++xx)
          {
            for(yy=xx+1; yy<inp_corr->iNumOFSubSele; ++yy)
            {
              mm += 1;
              fThisDist = 0.0;
              iTotPairs = 0;
      
              for(ii=0; ii<inp_corr->pSubSele[xx].nselatm; ++ii)
              {
                for(jj=0; jj<inp_corr->pSubSele[yy].nselatm; ++jj)
                {
                  if(inp_corr->iResFlag==1)
                  {
                    iIndexA = molecule->rawmol.presn[inp_corr->pSubSele[xx].selatm[ii]-1];
                    iIndexB = molecule->rawmol.presn[inp_corr->pSubSele[yy].selatm[jj]-1];
                  }
                  
                  else
                  {
                    iIndexA = atoi(molecule->rawmol.atmId[inp_corr->pSubSele[xx].selatm[ii]-1]);
                    iIndexB = atoi(molecule->rawmol.atmId[inp_corr->pSubSele[yy].selatm[jj]-1]);
                  }
                  
                  iProgNumA = inp_corr->piMasterSeleIndexes[iIndexA];
                  iProgNumB = inp_corr->piMasterSeleIndexes[iIndexB];
                  
                  fThisDist += inp_corr->ppfFrameMatrix[iProgNumA][iProgNumB];
                  iTotPairs += 1;
                }
              }
              fThisDist = fThisDist / (float) iTotPairs;
              fprintf(inp_corr->FVerbOutFile, "%12.3f   ", fThisDist);
            }
          }
        }
        fprintf(inp_corr->FVerbOutFile, "\n");
      }
    }
  }
  
  sprintf( outstring, "            ");
  return 12;
}

void *DCorStep1(void *pviThreadNum)
{

  int         ii, jj;
  int         tt, oo, rr, qq;
  double fThisDist;
  
  tt = (int) pviThreadNum;

  for(oo=0; oo<pInpCorr->iNumOfPairsPerThreads; ++oo)
  {
    // rr and qq are frames indexes
    rr = pInpCorr->ppiPairsListFrm1[tt][oo];
    qq = pInpCorr->ppiPairsListFrm2[tt][oo];
    
    if(rr == -1 || qq == -1)
      break;
    
    if(rr != qq)
    {
      // ii is a residue index
      for(ii=0; ii<pInpCorr->iMatrixDim; ++ii)
      {
        
        fThisDist = DistanceCoor(pInpCorr->pppdCoord[rr][ii][0],
                                 pInpCorr->pppdCoord[rr][ii][1],
                                 pInpCorr->pppdCoord[rr][ii][2],
                                 pInpCorr->pppdCoord[qq][ii][0],
                                 pInpCorr->pppdCoord[qq][ii][1],
                                 pInpCorr->pppdCoord[qq][ii][2],
                                 pInpCorr->pPBCInfo);

        pInpCorr->ppfResFrameAvgDist[ii][rr] += fThisDist;
        pInpCorr->pfResOverallAvgDist[ii]    += fThisDist;
        
        //pInpCorr->pppfThreadResFrameAvgDist[tt][ii][oo] += fThisDist;
        //pInpCorr->ppfThreadResOverallAvgDist[tt][ii]    += fThisDist;
      }
    }
  }

  return NULL;
}



void *DCorStep2(void *pviThreadNum)
{

  int         ii, jj;
  int         tt, oo, rr, qq;
  double fTempValA, fTempValB, fTempValC;
  
  tt = (int) pviThreadNum;

  for(oo=0; oo<pInpCorr->iNumOfPairsPerThreads; ++oo)
  {
    // rr and qq are frames indexes
    rr = pInpCorr->ppiPairsListFrm1[tt][oo];
    qq = pInpCorr->ppiPairsListFrm2[tt][oo];
    
    if(rr == -1 || qq == -1)
      break;
    
    // ii and jj are residue indexes
    for(ii=0; ii<pInpCorr->iMatrixDim; ++ii)
    {
      for(jj=ii; jj<pInpCorr->iMatrixDim; ++jj)
      {
        if(rr == qq)
        {
          fTempValA = (pInpCorr->ppfResFrameAvgDist[ii][rr] * -2) + pInpCorr->pfResOverallAvgDist[ii];
          
          if(ii == jj)
            fTempValB = fTempValA;
          
          else
            fTempValB = (pInpCorr->ppfResFrameAvgDist[jj][rr] * -2) + pInpCorr->pfResOverallAvgDist[jj];
        }
        
        else
        {
          fTempValA = DistanceCoor(pInpCorr->pppdCoord[rr][ii][0],
                                   pInpCorr->pppdCoord[rr][ii][1],
                                   pInpCorr->pppdCoord[rr][ii][2],
                                   pInpCorr->pppdCoord[qq][ii][0],
                                   pInpCorr->pppdCoord[qq][ii][1],
                                   pInpCorr->pppdCoord[qq][ii][2],
                                   pInpCorr->pPBCInfo);
          fTempValA = fTempValA - pInpCorr->ppfResFrameAvgDist[ii][rr] - pInpCorr->ppfResFrameAvgDist[ii][qq] + pInpCorr->pfResOverallAvgDist[ii];
          
          if(ii == jj)
            fTempValB = fTempValA;
          
          else
          {
            fTempValB = DistanceCoor(pInpCorr->pppdCoord[rr][jj][0],
                                     pInpCorr->pppdCoord[rr][jj][1],
                                     pInpCorr->pppdCoord[rr][jj][2],
                                     pInpCorr->pppdCoord[qq][jj][0],
                                     pInpCorr->pppdCoord[qq][jj][1],
                                     pInpCorr->pppdCoord[qq][jj][2],
                                     pInpCorr->pPBCInfo);
            fTempValB = fTempValB - pInpCorr->ppfResFrameAvgDist[jj][rr] - pInpCorr->ppfResFrameAvgDist[jj][qq] + pInpCorr->pfResOverallAvgDist[jj];
          }
        }

        pInpCorr->pppfThreadDistCov[tt][ii][jj] += (fTempValA * fTempValB);
      }
    }
  }

  return NULL;
}

int Post_CORR(struct inp_corr *inp_corr, Molecule *molecule, struct sopt *OPT, CoorSet *trj_crd)
{
  int           ii, jj, ff, rr, qq, tt, xx, yy, mm;
  int           iProgNumA, iProgNumB;
  int           iIndexA, iIndexB;
  int           iRowNum, iColNum, iPos;
  int           iTotPairs;
  int           iMatrixSize, iVectorSize;
  int           iRunningThreads, iRunnedThreads, iExitStatus;

  float         fThisDist, fThisFluct, fThisFlex;
  float         fTempValA=0.0, fTempValB=0.0;
  double        fSqrNumOfFrames;
  
  double        dCorrValue=0.0;
  double      **ppdMatrix;
  
  char          cLabelA[20], cLabelB[20], cResCode1[3], cTemp[150];
  
  time_t        time_Today;
  
  
  if(inp_corr->iFirstRoundFlag == 1)
  {
    // fluct calculation needs a second run over trj!
    
    if(inp_corr->iVerboseFlag == 1)
      printf("{Corr::Fluct} End of first run; calculating average distances and rewind trj\n");
    
    inp_corr->iFirstRoundFlag = 0;
    inp_corr->iFrameNum = 0;
    
    for(ii=0; ii<inp_corr->iMatrixDim; ++ii)
      for(jj=ii+1; jj<inp_corr->iMatrixDim; ++jj)
        inp_corr->ppfMeanDistance[ii][jj] = inp_corr->ppfMeanDistance[ii][jj] / inp_corr->iNumOfFrames;

    return 999;
  }
  
  if(inp_corr->iVerboseFlag == 1)
    printf("{Corr::Fluct} Entering in Post_CORR\n");

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
  
  else if(inp_corr->iCorrType == 2)
  {
    // --TYPE DCOR *-*
    inp_corr->iFrameNum = inp_corr->iNumOfFrames;
    
    fSqrNumOfFrames = (double) (inp_corr->iNumOfFrames * inp_corr->iNumOfFrames);

    // === 1st Step =============================================== //
    for(tt=0; tt<inp_corr->iNumOfThreads; ++tt)
      iExitStatus = pthread_create(&inp_corr->corr_threads[tt], NULL, DCorStep1, (void *) tt);

    for(tt=0; tt<inp_corr->iNumOfThreads; ++tt)
    {
      if(pthread_join(inp_corr->corr_threads[tt], NULL))
      {
        fprintf(stderr, "Corr Module: something bad happened to thread #%d in the 1st calculation step\n", tt);
        exit(5);
      }
    }
    
    // === merging togheter DCorStep1 results from different threads
    //for(tt=0; tt<inp_corr->iNumOfThreads; ++tt)
    //{
    //  for(ii=0; ii<inp_corr->iMatrixDim; ++ii)
    //  {
    //    inp_corr->pfResOverallAvgDist[ii] += inp_corr->ppfThreadResOverallAvgDist[tt][ii];
    //
    //    for(jj=0; jj<inp_corr->iNumOfPairsPerThreads; ++jj)
    //    {
    //      
    //      ff = inp_corr->ppiPairsListFrm1[tt][jj];
    //      if(ff != -1 && inp_corr->pppfThreadResFrameAvgDist[tt][ii][jj] != 0.0)
    //      {
    //        inp_corr->ppfResFrameAvgDist[ii][ff] += inp_corr->pppfThreadResFrameAvgDist[tt][ii][jj];
    //        inp_corr->pppfThreadResFrameAvgDist[tt][ii][jj] = 0.0;
    //      }
    //    }
    //  }
    //}

    // calculating ai.|a.j and a..
    for(ii=0; ii<inp_corr->iMatrixDim; ++ii)
    {
      inp_corr->pfResOverallAvgDist[ii] = inp_corr->pfResOverallAvgDist[ii] / fSqrNumOfFrames;
      for(tt=0; tt<inp_corr->iNumOfFrames; ++tt)
        inp_corr->ppfResFrameAvgDist[ii][tt] = inp_corr->ppfResFrameAvgDist[ii][tt] /  ((double) inp_corr->iNumOfFrames);
    }
    
    // === 2nd Step ================================================= //
    for(tt=0; tt<inp_corr->iNumOfThreads; ++tt)
      iExitStatus = pthread_create(&inp_corr->corr_threads[tt], NULL, DCorStep2, (void *) tt);

    for(tt=0; tt<inp_corr->iNumOfThreads; ++tt)
    {
      if(pthread_join(inp_corr->corr_threads[tt], NULL))
      {
        fprintf(stderr, "Corr Module: something bad happened to thread #%d in the 2nd calculation step\n", tt);
        exit(5);
      }
    }
    
    // === merging all pppfThreadDistCov data togheter
    for(tt=0; tt<inp_corr->iNumOfThreads; ++tt)
      for(ii=0; ii<inp_corr->iMatrixDim; ++ii)
        for(jj=ii; jj<inp_corr->iMatrixDim; ++jj)
          inp_corr->ppfDistCov[ii][jj] += inp_corr->pppfThreadDistCov[tt][ii][jj];

    // === finalizing =============================================== //
    for(rr=0; rr<inp_corr->iMatrixDim; ++rr)
      for(qq=rr; qq<inp_corr->iMatrixDim; ++qq)
        inp_corr->ppfDistCov[rr][qq] = sqrt(inp_corr->ppfDistCov[rr][qq] / fSqrNumOfFrames);

    for(rr=0; rr<inp_corr->iMatrixDim; ++rr)
    {
      inp_corr->ppfDCorr[rr][rr] = 1.0;
      for(qq=rr+1; qq<inp_corr->iMatrixDim; ++qq)
        inp_corr->ppfDCorr[rr][qq] = inp_corr->ppfDistCov[rr][qq] / (sqrt(inp_corr->ppfDistCov[rr][rr]*inp_corr->ppfDistCov[qq][qq]));
    }
  }
  
  else if(inp_corr->iCorrType == 3)
  {
    // --TYPE FLUCT @@@
    
    // finalizes fluctuation matrix calculation
    for(ii=0; ii<inp_corr->iMatrixDim; ++ii)
    {
      for(jj=ii+1; jj<inp_corr->iMatrixDim; ++jj)
      {
        inp_corr->ppfFluctMatrix[ii][jj] = inp_corr->ppfFluctMatrix[ii][jj] / (float) inp_corr->iNumOfFrames;

        // calculates standard deviation
        if(inp_corr->iStdDevFlag == 1)
          inp_corr->ppfFluctMatrix[ii][jj] = sqrt(inp_corr->ppfFluctMatrix[ii][jj]);
      }
    }

    // calculates the overall fluctuation over master selection
    inp_corr->fOverallFluct   = 0.0;
    inp_corr->fOverallAvgDist = 0.0;
    for(ii=0; ii<inp_corr->iMatrixDim; ++ii)
    {
      for(jj=0; jj<inp_corr->iMatrixDim; ++jj)
      {
        if(ii < jj)
        {
          inp_corr->fOverallFluct   += inp_corr->ppfFluctMatrix[ii][jj];
          inp_corr->fOverallAvgDist += inp_corr->ppfMeanDistance[ii][jj];
        }
        if(ii > jj)
        {
          inp_corr->fOverallFluct   += inp_corr->ppfFluctMatrix[jj][ii];
          inp_corr->fOverallAvgDist += inp_corr->ppfMeanDistance[jj][ii];
        }
      }
    }
    inp_corr->fOverallFluct   = sqrt(inp_corr->fOverallFluct / (float) (inp_corr->iMatrixDim*inp_corr->iMatrixDim));
    inp_corr->fOverallAvgDist = inp_corr->fOverallAvgDist / (float) (inp_corr->iMatrixDim*inp_corr->iMatrixDim);
    
    
    if(inp_corr->iNumOFSubSele > 0)
    {
      // calculates the sub-sele fluctuations
      for(xx=0; xx<inp_corr->iNumOFSubSele; ++xx)
      {
        iTotPairs = 0;
        inp_corr->pfSubSeleFluct[xx]   = 0.0;
        inp_corr->pfSubSeleAvgDist[xx] = 0.0;
        for(ii=0; ii<inp_corr->pSubSele[xx].nselatm; ++ii)
        {
          for(jj=0; jj<inp_corr->pSubSele[xx].nselatm; ++jj)
          {
            if(inp_corr->iResFlag==1)
            {
              iIndexA = molecule->rawmol.presn[inp_corr->pSubSele[xx].selatm[ii]-1];
              iIndexB = molecule->rawmol.presn[inp_corr->pSubSele[xx].selatm[jj]-1];
            }
            
            else
            {
              iIndexA = atoi(molecule->rawmol.atmId[inp_corr->pSubSele[xx].selatm[ii]-1]);
              iIndexB = atoi(molecule->rawmol.atmId[inp_corr->pSubSele[xx].selatm[jj]-1]);
            }
            
            iProgNumA = inp_corr->piMasterSeleIndexes[iIndexA];
            iProgNumB = inp_corr->piMasterSeleIndexes[iIndexB];
            
            if(iProgNumA < iProgNumB)
            {
              inp_corr->pfSubSeleFluct[xx]   += inp_corr->ppfFluctMatrix[iProgNumA][iProgNumB];
              inp_corr->pfSubSeleAvgDist[xx] += inp_corr->ppfMeanDistance[iProgNumA][iProgNumB];
              
            }
            
            if(iProgNumB < iProgNumA)
            {
              inp_corr->pfSubSeleFluct[xx]   += inp_corr->ppfFluctMatrix[iProgNumB][iProgNumA];
              inp_corr->pfSubSeleAvgDist[xx] += inp_corr->ppfMeanDistance[iProgNumB][iProgNumA];
            }
            
            iTotPairs += 1;
          }
        }
        
        inp_corr->pfSubSeleFluct[xx]   = sqrt(inp_corr->pfSubSeleFluct[xx] / (float) iTotPairs);
        inp_corr->pfSubSeleAvgDist[xx] = inp_corr->pfSubSeleAvgDist[xx] / (float) iTotPairs;
      }
      
      
      if(inp_corr->iNumOfSubMatch > 0)
      {
        // calculates matcthed sub-sele fluctuations
        // i.e. the overall-fluctuations between all res/atm of one sub-sele vs the res/atm of onother sele
        mm = -1;
        for(xx=0; xx<inp_corr->iNumOFSubSele; ++xx)
        {
          for(yy=xx+1; yy<inp_corr->iNumOFSubSele; ++yy)
          {
            mm += 1;
            inp_corr->pfMatchSubSeleFluct[mm] = 0.0;
            iTotPairs = 0;

            for(ii=0; ii<inp_corr->pSubSele[xx].nselatm; ++ii)
            {
              for(jj=0; jj<inp_corr->pSubSele[yy].nselatm; ++jj)
              {
                if(inp_corr->iResFlag==1)
                {
                  iIndexA = molecule->rawmol.presn[inp_corr->pSubSele[xx].selatm[ii]-1];
                  iIndexB = molecule->rawmol.presn[inp_corr->pSubSele[yy].selatm[jj]-1];
                }
                
                else
                {
                  iIndexA = atoi(molecule->rawmol.atmId[inp_corr->pSubSele[xx].selatm[ii]-1]);
                  iIndexB = atoi(molecule->rawmol.atmId[inp_corr->pSubSele[yy].selatm[jj]-1]);
                }
                
                iProgNumA = inp_corr->piMasterSeleIndexes[iIndexA];
                iProgNumB = inp_corr->piMasterSeleIndexes[iIndexB];
                
                if(iProgNumA < iProgNumB)
                {
                  inp_corr->pfMatchSubSeleFluct[mm]   += inp_corr->ppfFluctMatrix[iProgNumA][iProgNumB];
                  inp_corr->pfMatchSubSeleAvgDist[mm] += inp_corr->ppfMeanDistance[iProgNumA][iProgNumB];
                }
                if(iProgNumB < iProgNumA)
                {
                  inp_corr->pfMatchSubSeleFluct[mm]   += inp_corr->ppfFluctMatrix[iProgNumB][iProgNumA];
                  inp_corr->pfMatchSubSeleAvgDist[mm] += inp_corr->ppfMeanDistance[iProgNumB][iProgNumA];
                }
                
                iTotPairs += 1;
              }
            }
            inp_corr->pfMatchSubSeleFluct[mm]   = sqrt(inp_corr->pfMatchSubSeleFluct[mm] / (float) iTotPairs);
            inp_corr->pfMatchSubSeleAvgDist[mm] = inp_corr->pfMatchSubSeleAvgDist[mm] / (float) iTotPairs;
          }
        }
      }
    }
  }
  
  // ===================================================================
  
  // === It's time for results !!! =====================================
  time(&time_Today);
  fprintf(inp_corr->FOutFile, "# ==============================================================\n");
  fprintf(inp_corr->FOutFile, "# ***                   WORDOM CORR MODULE                   ***\n");
  fprintf(inp_corr->FOutFile, "# ==============================================================\n");
  fprintf(inp_corr->FOutFile, "#\n");
  fprintf(inp_corr->FOutFile, "# Version        : %s\n", VERSTRING);
  fprintf(inp_corr->FOutFile, "# License        : GPL 3\n");
  fprintf(inp_corr->FOutFile, "# Copyright      : Fanelli, Felline\n");
  fprintf(inp_corr->FOutFile, "#                  University of Modena\n");
  fprintf(inp_corr->FOutFile, "#                  Modena - Italy\n");
  fprintf(inp_corr->FOutFile, "#\n");
  fprintf(inp_corr->FOutFile, "# Date           : %s", asctime(localtime(&time_Today)));
  fprintf(inp_corr->FOutFile, "#\n");
  fprintf(inp_corr->FOutFile, "# Mol File       : %s\n", OPT->IMOL_FILE);
  fprintf(inp_corr->FOutFile, "# Res Num        : %d\n", inp_corr->iNumOfRes);
  fprintf(inp_corr->FOutFile, "# Traj File      : %s\n", OPT->ITRJ_FILE);
  fprintf(inp_corr->FOutFile, "# Frame Num      : %d\n", inp_corr->iNumOfFrames);
  
  if(trj_crd->pbc_flag == 0)
    fprintf(inp_corr->FOutFile, "# PBC            : No\n");
  else
    fprintf(inp_corr->FOutFile, "# PBC            : Yes\n");
  
  fprintf(inp_corr->FOutFile, "#\n");
  fprintf(inp_corr->FOutFile, "# Title          : %s\n", inp_corr->cTitle);
  
  if(inp_corr->iCorrType == 0)
    fprintf(inp_corr->FOutFile, "# Type           : DCC\n");
  else if(inp_corr->iCorrType == 1)
    fprintf(inp_corr->FOutFile, "# Type           : LMI\n");
  else if(inp_corr->iCorrType == 2)
    fprintf(inp_corr->FOutFile, "# Type           : DCOR\n");
  else if(inp_corr->iCorrType == 3)
    fprintf(inp_corr->FOutFile, "# Type           : FLUCT\n");
  
  fprintf(inp_corr->FOutFile, "# Sele           : %s\n", inp_corr->sele.selestring);
  
  if(inp_corr->iResFlag==1)
    fprintf(inp_corr->FOutFile, "# Sele Res       : %d\n", inp_corr->iMatrixDim);
  else
    fprintf(inp_corr->FOutFile, "# Sele Atm       : %d\n", inp_corr->iMatrixDim);
    
  if(inp_corr->iResFlag==1)
    fprintf(inp_corr->FOutFile, "# Level          : RES\n");
  else
    fprintf(inp_corr->FOutFile, "# Level          : ATM\n");
  
  if(inp_corr->iMassFlag==1)
    fprintf(inp_corr->FOutFile, "# Mass           : Yes\n");
  else
    fprintf(inp_corr->FOutFile, "# Mass           : No\n");
  
  if(inp_corr->iCorrType == 3)
  {
    if(inp_corr->iStdDevFlag == 1)
      fprintf(inp_corr->FOutFile, "# StDev Flag     : Yes\n");
    else
      fprintf(inp_corr->FOutFile, "# StDev Flag     : No\n");
  }

  fprintf(inp_corr->FOutFile, "# Num of Threads : %d\n", inp_corr->iNumOfThreads);
  fprintf(inp_corr->FOutFile, "#\n");
  
  if(inp_corr->iMultiAtomFlag == 1 && trj_crd->pbc_flag == 1)
    fprintf(inp_corr->FOutFile, "# Warning!     : PBC not taken into account in the calculation of geo/mass centers\n#\n");
    
  if(inp_corr->iZeroMassWarning == 1)
    fprintf(inp_corr->FOutFile, "# Warning!     : Atom(s) with zero mass\n#\n");
  
  if(inp_corr->iCorrType == 3)
  {
    fprintf(inp_corr->FOutFile, "# Overall Fluct  : %s %.3f\n", inp_corr->sele.selestring, inp_corr->fOverallFluct);
    fprintf(inp_corr->FOutFile, "# Overall Dist   : %s %.3f\n", inp_corr->sele.selestring, inp_corr->fOverallAvgDist);
  }
  
  if(inp_corr->iNumOFSubSele > 0)
  {
    fprintf(inp_corr->FOutFile, "#\n");
    fprintf(inp_corr->FOutFile, "# %3s   %-50s   %8s   %8s\n", "Num", "SubSele", "Fluct", "Dist");
    for(xx=0; xx<inp_corr->iNumOFSubSele; ++xx)
      fprintf(inp_corr->FOutFile, "# %3d   %-50s   %8.3f   %8.3f\n", xx+1, inp_corr->pSubSele[xx].selestring, inp_corr->pfSubSeleFluct[xx], inp_corr->pfSubSeleAvgDist[xx]);
  }
  
  if(inp_corr->iNumOfSubMatch > 0)
  {
    fprintf(inp_corr->FOutFile, "#\n");
    fprintf(inp_corr->FOutFile, "# %3s   %-50s   %-50s   %8s   %8s\n", "Num", "SubSele1", "SubSele2", "Fluct", "Dist");
    mm = -1;
    for(xx=0; xx<inp_corr->iNumOFSubSele; ++xx)
    {
      for(yy=xx+1; yy<inp_corr->iNumOFSubSele; ++yy)
      {
        mm += 1;
        fprintf(inp_corr->FOutFile, "# %3d   %-50s   %-50s   %8.3f   %8.3f\n", mm+1, inp_corr->pSubSele[xx].selestring, inp_corr->pSubSele[yy].selestring, inp_corr->pfMatchSubSeleFluct[mm], inp_corr->pfMatchSubSeleAvgDist[mm]);
      }
    }
  }
  
  
  if(inp_corr->iCorrType == 3)
  {
    fprintf(inp_corr->FOutFile, "#\n");
    fprintf(inp_corr->FOutFile, "# Info: Fluct describes the extent of the pairwise fluctations\n");
    fprintf(inp_corr->FOutFile, "# Info: Flex is the difference between the largest and the smallest dist\n");
    fprintf(inp_corr->FOutFile, "# Info: AvgDist is the average distance\n");
    fprintf(inp_corr->FOutFile, "# ====================================================================================\n");
    fprintf(inp_corr->FOutFile, "#%8s   %8s   %15s   %15s   %7s   %7s   %7s\n", "ProgNumA", "ProgNumB", "LabelA", "LabelB", "Fluct", "Flex", "AvgDist");
  }
  
  else
  {
    fprintf(inp_corr->FOutFile, "# ==============================================================\n");
    fprintf(inp_corr->FOutFile, "#%8s   %8s   %15s   %15s   %5s\n", "ProgNumA", "ProgNumB", "LabelA", "LabelB", "CORR");
  }
  
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
      
      else if(inp_corr->iCorrType == 2)
      {
        if(ii == jj)
          dCorrValue = 1.0;
        else if(ii < jj)
          dCorrValue = inp_corr->ppfDCorr[ii][jj];
        else if(ii > jj)
          dCorrValue = inp_corr->ppfDCorr[jj][ii];
      }
      
      else if(inp_corr->iCorrType == 3)
      {
        if(ii == jj)
        {
          fThisFluct = 0.0;
          fThisFlex  = 0.0;
          fThisDist  = 0.0;
        }
        
        else if(ii < jj)
        {
          fThisFluct = inp_corr->ppfFluctMatrix[ii][jj];
          fThisFlex  = inp_corr->ppfMaxDistances[ii][jj] - inp_corr->ppfMinDistances[ii][jj];
          fThisDist  = inp_corr->ppfMeanDistance[ii][jj];
        }
        
        else if(ii > jj)
        {
          fThisFluct = inp_corr->ppfFluctMatrix[jj][ii];
          fThisFlex  = inp_corr->ppfMaxDistances[jj][ii] - inp_corr->ppfMinDistances[jj][ii];
          fThisDist  = inp_corr->ppfMeanDistance[jj][ii];
        }

      }
      
      if(inp_corr->iCorrType == 3)
        fprintf(inp_corr->FOutFile, " %8d   %8d   %15s   %15s   %7.3f   %7.3f   %7.3f\n", iProgNumA, iProgNumB, cLabelA, cLabelB, fThisFluct, fThisFlex, fThisDist);
      else
        fprintf(inp_corr->FOutFile, " %8d   %8d   %15s   %15s   %5.3f\n", iProgNumA, iProgNumB, cLabelA, cLabelB, dCorrValue);
    }
  }
  
  close(inp_corr->FOutFile);
  return 0;
  // ===================================================================
}
