// -------------------------------------------------------------------------
// Copyright (C) 2009  Francesca Fanelli, Angelo Felline, Michele Seeber
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
#include <ctype.h>
#include "wordom.h"
#include "tools.h"
#include "fileio.h"
#include "datahandler.h"
#include "psn.h"
#include <signal.h>

#define TRUE  1
#define FALSE 0

int linkwalk(int nres, int **linklist, int *nlinks, struct simplecluster *cluster)
{
  int          ii, jj, kk;
  int          nclusters, thiscluster, *checked_links, *checked_res;
  int          *cluspop, **clusters;
  
  checked_links = calloc( nres, sizeof(int));        // whether the link to a res has been checked/added
  checked_res = calloc( nres, sizeof(int));          // whether a res has been checked/added
  
  for(ii=0; ii<nres; ii++)
    cluster->cluspop[ii]=0;
    
  for(ii=0; ii<nres; ii++)
    for(jj=0; jj<nres; jj++)
      cluster->clusters[ii][jj]=0;
    
  cluspop = cluster->cluspop;
  clusters = cluster->clusters;
  
  
  nclusters = 0;
  for( ii=0; ii<nres; ii++)
    if( checked_res[ii] == 0 && nlinks[ii]>0 )
    {
      // this residue is not included in any previous cluster! a new cluster is born
      thiscluster = nclusters;
      nclusters++;
      clusters[thiscluster][cluspop[thiscluster]] = ii;
      cluspop[thiscluster]++;
      checked_links[ii] = 1;
      
      for(jj=0; jj<nlinks[ii]; jj++)
      {
        clusters[thiscluster][cluspop[thiscluster]] = linklist[ii][jj];
        cluspop[thiscluster]++;
        checked_links[linklist[ii][jj]] = 1;
      }
      
      checked_res[ii] = 1;
      // now let's follow its links and gather all the followers
      for( jj=0; jj<cluspop[thiscluster]; jj++)
        if( checked_res[clusters[thiscluster][jj]] == 0)
        {
          for( kk=0; kk<nlinks[clusters[thiscluster][jj]]; kk++)
          {
            if( checked_links[linklist[clusters[thiscluster][jj]][kk]] == 0)
            {
              clusters[thiscluster][cluspop[thiscluster]]  = linklist[clusters[thiscluster][jj]][kk];
              cluspop[thiscluster]++;
              checked_links[linklist[clusters[thiscluster][jj]][kk]] = 1;
            }
          }
          checked_res[clusters[thiscluster][jj]] = 1;
        }
    }
  
  free(checked_links);
  free(checked_res);
  
  cluster->nclusters = nclusters;
  
  return nclusters;
}

// ------------------------------------------------------------------

float fGetNormFactor(char *resType)
{
  // this returns the normalization value for each residue according to
  // Kannan and Vishveshwara, JMB,(1999)292,441-464
  
  int          ii;
  float        norm;
  char         thisres[5] = { '\0', '\0', '\0', '\0', '\0' };

  for(ii=0; ii<strlen(resType); ii++)
    thisres[ii] = toupper(resType[ii]);
  
  if(!strcmp(      thisres, "ALA"))
    norm = 55.7551;

  else if(!strcmp( thisres, "ARG"))
    norm = 93.7891;

  else if(!strcmp( thisres, "ARGN"))
    norm = 93.7891;

  else if(!strcmp( thisres, "ASN"))
    norm = 73.4097;

  else if(!strcmp( thisres, "ASP"))
    norm = 75.1507;

  else if(!strcmp( thisres, "ASPH"))
    norm = 75.1507;

  else if(!strcmp( thisres, "CYS"))
    norm = 54.9528;

  else if(!strcmp( thisres, "CYN"))
    norm = 54.9528;

  else if(!strcmp( thisres, "GLN"))
    norm = 78.1301;

  else if(!strcmp( thisres, "GLU"))
    norm = 78.8288;

  else if(!strcmp( thisres, "GLUH"))
    norm = 78.8288;

  else if(!strcmp( thisres, "GLY"))
    norm = 47.3129;

  else if(!strcmp( thisres, "HIS"))
    norm = 83.7357;

  else if(!strcmp( thisres, "HID"))
    norm = 83.7357;

  else if(!strcmp( thisres, "HSD"))
    norm = 83.7357;

  else if(!strcmp( thisres, "HSC"))
    norm = 83.7357;
    
  else if(!strcmp( thisres, "HSP"))
    norm = 83.7357;
    
  else if(!strcmp( thisres, "HSE"))
    norm = 83.7357;
  
  else if(!strcmp( thisres, "ILE"))
    norm = 67.9452;
    
  else if(!strcmp( thisres, "LEU"))
    norm = 72.2517;
    
  else if(!strcmp( thisres, "LYS"))
    norm = 69.6096;
    
  else if(!strcmp( thisres, "LYP"))
    norm = 69.6096;
    
  else if(!strcmp( thisres, "LYSN"))
    norm = 69.6096;
    
  else if(!strcmp( thisres, "MET"))
    norm = 69.2569;
    
  else if(!strcmp( thisres, "PHE"))
    norm = 93.3082;
    
  else if(!strcmp( thisres, "PRO"))
    norm = 51.331;
    
  else if(!strcmp( thisres, "SER"))
    norm = 61.3946;
    
  else if(!strcmp( thisres, "THR"))
    norm = 63.7075;
    
  else if(!strcmp( thisres, "TRP"))
    norm = 106.703;
    
  else if(!strcmp( thisres, "TYR"))
    norm = 100.719;
    
  else if(!strcmp( thisres, "VAL"))
    norm = 62.3673;
  
  /* === Normalization Factors added by Fanelli and Co-Workers === */
  
  // 11-cis-retinal
  else if(!strcmp( thisres, "RET"))
    norm = 170.1355;
  
  // GDP
  else if(!strcmp( thisres, "GDP"))
    norm = 220.1921;
   
  // GTP
  else if(!strcmp( thisres, "GTP"))
    norm = 274.7802;
    
  // GTP in Galpha (PDB : 1CIP, 1CUL and 1TND)
  else if(!strcmp( thisres, "GTA"))
    norm = 361.3333;
   
  // GDP in Galpha (PDB : 1TAG)
  else if(!strcmp( thisres, "GDA"))
    norm = 293.0;

  // Mg ion in GTP-Bound forms
  else if(!strcmp( thisres, "MGT"))
    norm = 22.0147;
  
  // Mg ion in GDP-Bound forms
  else if(!strcmp( thisres, "MGD"))
    norm = 14.6585;
  
  // Mg ion in GTP- & GDP- Bound forms
  else if(!strcmp( thisres, "MGX"))
    norm = 19.2477;

  // Mg ion in Galpha GTP- & GDP- Bound forms (PDB : 1CIP, 1CUL, 1TAG and 1TND)
  else if(!strcmp( thisres, "MGA"))
    norm = 23.83333;
  
  // Structural water in 1CIP, 1CUL, 1TAG, 1TND and 1GZM
  else if(!strcmp( thisres, "H2O"))
    norm = 27.0;
  // generic water in solvent
  else if(!strcmp( thisres, "HOH") || !strcmp( thisres, "WAT") || !strcmp( thisres, "SOL"))
    norm = 27.0;
  
  // ZM241385 - R.C. Stevens et al, Science vol 322, 21 Nov 2008
  else if(!strcmp( thisres, "ZMA"))
    norm = 139.000;
  
  // Lys296 covalently linked to 11-cis-retinal in rhodopsins
  else if(!strcmp( thisres, "LYR"))
    norm = 109.9498;

  // Lys296 + 11-cis-retinal in rhodopsins
  else if(!strcmp( thisres, "KRT"))
    norm = 262.5612;
  
  // MN 
  else if(!strcmp( thisres, "MNG"))
    norm = 23.5;
  
  else
  {
    fprintf( stderr, "Could not recognize residue >%s<: will take norm as 999999.99 (low interaction propensity)\n\n", resType);
    norm = 999999.99;
  }
  
  return norm;
}

float GetUserNormFactor(struct inp_psn *inp_psn, char *resType)
{
  int              ii;

  for(ii=0; ii<inp_psn->iNumOfParam; ii++)
    if(strcmp(resType, inp_psn->pcParamResVect[ii]) == 0)
      return inp_psn->pfParamValVect[ii];
  
  return -999.999;
}

// === GetImin Section =================================================
void GetImin(struct inp_psn *inp_psn, struct sopt *OPT, Molecule *molecule)
{
  int           ii, jj, kk, mm;
  int           iIterNum=0;
  int          *piNodeClusters;
  int           iNumOfNodes;
  int           iNumOfFrames;
  int          *piClustSize;
  int         **ppiIntMatrix;
  int           iDecNum=2;
  int           DEBUG=0;
  
  char          cResName[6];
  
  double        fI0Size, fIcSize;
  double        fLowerImin=0.0, fHigherImin=10.0, fHalfImin;
  double        fHalfIminBigClsSize;
  double        fEpsilon=0.01;
  double        fLastHalfImin;
  double        fIminTest, fIminTestSize;
  double        fLastIminTest, fLastIminTestSize;
  double        fPCNSize=50.0;
  double        fPSNDeltaA, fPSNDeltaB;

  FILE         *FRawFile;

  // default values //
  fLowerImin  =  0.0;
  fHigherImin = 10.0;
  // ============== //
  
  // close "old" raw_file
  fclose(inp_psn->outfile);
  
  fPCNSize  = inp_psn->fPCNSize;
  iDecNum   = inp_psn->iDecNum;
  fEpsilon  = (1.0 / (float) (pow(10.0, (float) iDecNum)));
  fHalfImin = (((fHigherImin - fLowerImin) / 2.0) + fLowerImin);
  //fIcSize   = (fIcSize / fPCNSize); It is uninitialized!
  
  FRawFile = fopen(inp_psn->cRawFileName, "r");
  if(FRawFile == NULL)
  {
    printf("PSN Module : Unable to open file : %s\n", inp_psn->cRawFileName);
    exit(1);
  }
  
  iNumOfNodes  = inp_psn->iNumOfNodes;
  iNumOfFrames = inp_psn->iNumOfFrames;
  
  // === Some Allocations === //
  piNodeClusters = (int  *) calloc(iNumOfNodes, sizeof(int  ));
  //piClustSize    = (int  *) calloc(1000,        sizeof(int  ));
  piClustSize    = (int  *) calloc(iNumOfNodes, sizeof(int  ));
  ppiIntMatrix   = (int **) calloc(iNumOfNodes, sizeof(int *));
  for(ii=0; ii<iNumOfNodes; ++ii)
    ppiIntMatrix[ii] = (int *) calloc(iNumOfNodes, sizeof(int));
  // ======================== //

  fI0Size = GetBigClsSize(FRawFile, piNodeClusters, piClustSize, 0.0, iNumOfNodes, ppiIntMatrix);
  fPCNSize = (float) (100.0 / fPCNSize);
  fIcSize = (fI0Size / fPCNSize);
  
  if(DEBUG == 1)
  {
    printf("File      : %s\n", inp_psn->cRawFileName);
    printf("Nodes     : %d\n", iNumOfNodes);
    printf("Frames    : %d\n", iNumOfFrames);
    printf("Min Imin  : %f\n", fLowerImin);
    printf("Max Imin  : %f\n", fHigherImin);
    printf("Dec Num   : %d\n", iDecNum);
    printf("Epsilon   : %f\n", fEpsilon);
    printf("LCS Fract : %f\n", (100.0 / fPCNSize));
    printf("I0Size    : %f\n", fI0Size);
    printf("IcSize    : %f\n", fIcSize);
    printf("\n");
  }

  if(DEBUG == 1)
    printf(" #      Imin        Size         Delta\n");

  while(1)
  {
    iIterNum++;
    
    fHalfIminBigClsSize = GetBigClsSize(FRawFile, piNodeClusters, piClustSize, fHalfImin, iNumOfNodes, ppiIntMatrix);
    
    if(DEBUG == 1)
      printf("%2d      %f    %f    %+f\n", iIterNum, fHalfImin, fHalfIminBigClsSize, (fHalfIminBigClsSize - fIcSize));
    
    if(fHalfIminBigClsSize > fIcSize)
    {
      fLowerImin    = fHalfImin;
      fHigherImin   = fHigherImin;
      fLastHalfImin = fHalfImin;
      fHalfImin     = (((fHigherImin - fLowerImin) / 2.0) + fLowerImin);
    }
    
    else if(fHalfIminBigClsSize <= fIcSize)
    {
      fLowerImin    = fLowerImin;
      fHigherImin   = fHalfImin;
      fLastHalfImin = fHalfImin;
      fHalfImin     = (((fHigherImin - fLowerImin) / 2.0) + fLowerImin);
    }
    
    /*
    else
    {
      printf("*-* ==================================================== *-*\n");
      break;
    }*/
    
    if(fabs(fLastHalfImin - fHalfImin) <= fEpsilon)
      break;
  }

  if(DEBUG == 1)
    printf("\n*** convergence reached in %d step(s) with value %g, now refinement ***\n\n", iIterNum, fLastHalfImin);
  
  fIminTest = (floorf((fHalfImin * (1/fEpsilon)) + 0.5) / (1/fEpsilon)) - (fEpsilon*3);
  
  while(TRUE)
  {
    iIterNum++;
    fIminTestSize = GetBigClsSize(FRawFile, piNodeClusters, piClustSize, fIminTest, iNumOfNodes, ppiIntMatrix);
    
    if(DEBUG == 1)
      printf("%2d      %f    %f    %+f\n", iIterNum, fIminTest, fIminTestSize, (fIminTestSize - fIcSize));
    
    if(fIminTestSize <= fIcSize)
      break;
    
    else
    {
      fLastIminTest     = fIminTest;
      fLastIminTestSize = fIminTestSize;
      fIminTest         = fIminTest + fEpsilon;
    }
  }
  
  if(DEBUG == 1)
  {
    printf("\n");
    printf("=== Pre Ic Value ===\n");
  }
  
  fLastIminTestSize = GetBigClsSize(FRawFile, piNodeClusters, piClustSize, fLastIminTest, iNumOfNodes, ppiIntMatrix);
  fPSNDeltaA = (fLastIminTestSize - fIcSize);
  if(fPSNDeltaA < 0.0)
    printf("Warning: something goes wrong in Pre Ic Imin Calculation!! Call Angelo\n");

  fIminTestSize = GetBigClsSize(FRawFile, piNodeClusters, piClustSize, fIminTest, iNumOfNodes, ppiIntMatrix);
  fPSNDeltaB = (fIminTestSize - fIcSize);
  if(fPSNDeltaB > 0.0)
    printf("Warning: something goes wrong in Post Ic Imin Calculation!! Call Angelo\n");
    
  if(DEBUG == 1)
  {
    printf("Ic Value : %f\n", fLastIminTest);
    printf("PSN LCS  : %.2f%% (Delta: %+.2f)\n", fLastIminTestSize, fPSNDeltaA);
    printf("\n");
    printf("=== Post Ic Value ===\n");
    printf("Ic Value : %f\n", fIminTest);
    printf("PSN LCS  : %.2f%% (Delta: %+.2f)\n", fIminTestSize, fPSNDeltaB);
    printf("\n");
    printf("Best Ic Value         Imin              Size              Delta   Best\n");
    printf("----------------------------------------------------------------------\n");
  }
  
  inp_psn->fPreIcValue       = fLastIminTest;
  inp_psn->fPostIcValue      = fIminTest;
  inp_psn->fPreIcValueDelta  = fPSNDeltaA;
  inp_psn->fPostIcValueDelta = fPSNDeltaB;
  
  inp_psn->fIntMinStart = inp_psn->fPreIcValue;
  inp_psn->fIntMinStop  = inp_psn->fPostIcValue;
  inp_psn->fIntMinStep  = fEpsilon;
  
  if(fabs(fPSNDeltaA) < fabs(fPSNDeltaB))
  {
    inp_psn->iBestIcValue = 0;
    
    if(DEBUG == 1)
    {
      printf("Pre  Ic Value:      %12.8f       %12.8f       %+3.2f   *\n", fLastIminTest, fLastIminTestSize, fPSNDeltaA);
      printf("Post Ic Value:      %12.8f       %12.8f       %+3.2f    \n", fIminTest, fIminTestSize, fPSNDeltaB);
    }
  }
  
  else if(fabs(fPSNDeltaA) > fabs(fPSNDeltaB))
  {
    inp_psn->iBestIcValue = 1;
    
    if(DEBUG == 1)
    {
      printf("Pre  Ic Value:      %12.8f       %12.8f       %+3.2f    \n", fLastIminTest, fLastIminTestSize, fPSNDeltaA);
      printf("Post Ic Value:      %12.8f       %12.8f       %+3.2f   *\n", fIminTest, fIminTestSize, fPSNDeltaB);
    }
  }
  
  else
  {
    inp_psn->iBestIcValue = 2;
    
    if(DEBUG == 1)
    {
      printf("Pre  Ic Value:      %12.8f       %12.8f       %+3.2f   *\n", fLastIminTest, fLastIminTestSize, fPSNDeltaA);
      printf("Post Ic Value:      %12.8f       %12.8f       %+3.2f   *\n", fIminTest, fIminTestSize, fPSNDeltaB);
    }
  }

  fclose(FRawFile);
  
  // === De-Allocations === //
  free(piNodeClusters);
  free(piClustSize);
  for(ii=0; ii<iNumOfNodes; ++ii)
    free(ppiIntMatrix[ii]);
  // ======================== //

  inp_psn->outfile = fopen(inp_psn->cRawFileName, "w");
  
  fprintf( inp_psn->outfile, "# ======================================================================\n");
  fprintf( inp_psn->outfile, "#                       --- Protein Structure Network Analysis ---\n");
  fprintf( inp_psn->outfile, "#\n");
  fprintf( inp_psn->outfile, "# Copyright     : Francesca Fanelli\n");
  fprintf( inp_psn->outfile, "#                 Angelo Felline,\n");
  fprintf( inp_psn->outfile, "#                 Michele Seeber (2009)\n");
  fprintf( inp_psn->outfile, "# License       : GPL v 3\n");
  fprintf( inp_psn->outfile, "#\n");
  fprintf( inp_psn->outfile, "# PSN Ver       : 0.8c\n");
  fprintf( inp_psn->outfile, "#\n");
  fprintf( inp_psn->outfile, "# Date          : %s", asctime(localtime(&inp_psn->time_tToday)));
  fprintf( inp_psn->outfile, "#\n");
  fprintf( inp_psn->outfile, "# Mol Name      : %s\n", OPT->IMOL_FILE);
  fprintf( inp_psn->outfile, "# Seg Num       : %d\n", molecule->nSeg);
  fprintf( inp_psn->outfile, "# Res Num       : %d\n", inp_psn->tot_nresidues);
  fprintf( inp_psn->outfile, "# Trj Name      : %s\n", OPT->ITRJ_FILE);
  fprintf( inp_psn->outfile, "# Frame Num     : %d\n", inp_psn->iNumOfFrames);
  if(inp_psn->iPBCFlag == 0)
    fprintf(inp_psn->outfile, "# PBC           : No\n");
  else
    fprintf(inp_psn->outfile, "# PBC           : Yes\n");
  fprintf( inp_psn->outfile, "#\n");
  fprintf( inp_psn->outfile, "# Title         : %s\n", inp_psn->title);
  fprintf( inp_psn->outfile, "# Selection     : %s\n", inp_psn->sele.selestring);
  fprintf( inp_psn->outfile, "# Sele Res Num  : %d\n", inp_psn->iNumOfSelRes);
  fprintf( inp_psn->outfile, "# DistCutOff    : %f\n", inp_psn->fDistCutoff);
  fprintf( inp_psn->outfile, "# IntMin        : %s\n", inp_psn->cIminRange);
  fprintf( inp_psn->outfile, "# LCS Fract     : %f\n", inp_psn->fPCNSize);
  fprintf( inp_psn->outfile, "# Dec Pos       : %d\n", inp_psn->iDecNum);
  fprintf( inp_psn->outfile, "# Pre  Ic Val   : %f (delta %3.2f)\n", inp_psn->fPreIcValue, fabs(inp_psn->fPreIcValueDelta));
  fprintf( inp_psn->outfile, "# Post Ic Val   : %f (delta %3.2f)\n", inp_psn->fPostIcValue, fabs(inp_psn->fPostIcValueDelta));
  if(inp_psn->iBestIcValue == 0)
    fprintf( inp_psn->outfile, "# Best Ic Val   : PRE\n");
  else if(inp_psn->iBestIcValue == 1)
    fprintf( inp_psn->outfile, "# Best Ic Val   : POST\n");
  else if(inp_psn->iBestIcValue == 2)
    fprintf( inp_psn->outfile, "# Best Ic Val   : BOTH\n");
  
  fprintf( inp_psn->outfile, "# StableCutOff  : %f\n", inp_psn->fStableCutoff);
  
  if(inp_psn->iHubEqFlag == 0)
    fprintf( inp_psn->outfile, "# HubEquation   : NO\n");
  else
    fprintf( inp_psn->outfile, "# HubEquation   : YES\n");
  
  fprintf( inp_psn->outfile, "# HubContCutOff : %i\n", inp_psn->iHubContCutoff);
  fprintf( inp_psn->outfile, "# Termini       : %i\n", inp_psn->iTerminiFlag);
  fprintf( inp_psn->outfile, "# Proximity     : %i\n", inp_psn->iProxCutOff);
  fprintf( inp_psn->outfile, "#\n");
  fprintf( inp_psn->outfile, "# Num Of NoLink : %d\n", inp_psn->iNumOfNoLinkPairs);
  for(ii=0; ii<inp_psn->iNumOfNoLinkPairs; ii++)
    fprintf( inp_psn->outfile, "# No Link       : %d %d   %d\n", ii+1, inp_psn->ppiNoLinkPairs[ii][0]+1, inp_psn->ppiNoLinkPairs[ii][1]+1);
  
  fprintf(inp_psn->outfile, "# Num Of Merge  : %d\n", inp_psn->iNumOfMerge);
  for(ii=0; ii<inp_psn->iNumOfMerge; ii++)
    fprintf(inp_psn->outfile, "# Merge         : %d %s\n", ii+1, inp_psn->ppcMergedResidues[ii]);
  
  fprintf(inp_psn->outfile, "# Num Of Param  : %d\n", inp_psn->iNumOfParam);
  for(ii=0; ii<inp_psn->iNumOfParam; ii++)
    fprintf(inp_psn->outfile, "# Param         : %d %5s %f\n", ii+1, inp_psn->pcParamResVect[ii], inp_psn->pfParamValVect[ii]);
  
  fprintf( inp_psn->outfile, "# ======================================================================\n");
  fprintf( inp_psn->outfile, "\n");
  
  fprintf( inp_psn->outfile, "*** Seg Info ***\n");
  jj=1;
  for(ii=0; ii<molecule->nSeg; ii++)
  {
    fprintf( inp_psn->outfile, "%2d   %5s   %6d   %10d\n", ii+1, molecule->segment[ii].segName, molecule->segment[ii].nRpS, jj);
    jj=jj+molecule->segment[ii].nRpS;
  }
  fprintf( inp_psn->outfile, "========================\n");
  
  fprintf( inp_psn->outfile, "*** Seq ***\n");
  
  kk=0;
  for(ii=0; ii<molecule->nSeg; ii++)
  {
    for(jj=0; jj<molecule->segment[ii].nRpS; jj++)
    {
      kk++;
      if     (strncmp(molecule->segment[ii].pRes[jj].resType, "ALA", 3)==0) strcpy(cResName, "A\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "ARG", 3)==0) strcpy(cResName, "R\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "ASN", 3)==0) strcpy(cResName, "N\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "ASP", 3)==0) strcpy(cResName, "D\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "CYS", 3)==0) strcpy(cResName, "C\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "CYN", 3)==0) strcpy(cResName, "C\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "GLU", 3)==0) strcpy(cResName, "E\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "GLN", 3)==0) strcpy(cResName, "Q\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "GLY", 3)==0) strcpy(cResName, "G\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "HIS", 3)==0) strcpy(cResName, "H\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "HSC", 3)==0) strcpy(cResName, "H\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "HSD", 3)==0) strcpy(cResName, "H\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "HSP", 3)==0) strcpy(cResName, "H\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "HID", 3)==0) strcpy(cResName, "H\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "ILE", 3)==0) strcpy(cResName, "I\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "LEU", 3)==0) strcpy(cResName, "L\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "LYS", 3)==0) strcpy(cResName, "K\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "LYP", 3)==0) strcpy(cResName, "K\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "MET", 3)==0) strcpy(cResName, "M\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "PHE", 3)==0) strcpy(cResName, "F\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "PRO", 3)==0) strcpy(cResName, "P\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "SER", 3)==0) strcpy(cResName, "S\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "THR", 3)==0) strcpy(cResName, "T\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "TRP", 3)==0) strcpy(cResName, "W\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "TYR", 3)==0) strcpy(cResName, "Y\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "VAL", 3)==0) strcpy(cResName, "V\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "RET", 3)==0) strcpy(cResName, "X\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "GDP", 3)==0) strcpy(cResName, "d\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "GTP", 3)==0) strcpy(cResName, "t\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "Mg+", 3)==0) strcpy(cResName, "m\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "ZMA", 3)==0) strcpy(cResName, "Z\0");
      else                                                                  strcpy(cResName, "U\0");
      
      fprintf(inp_psn->outfile, "%10d   %s:%s%d\n", kk, molecule->segment[ii].segName, cResName, molecule->segment[ii].pRes[jj].resn);
    }
  }
  fprintf( inp_psn->outfile, "========================\n");

  // === re-initialize all vectors === //
  inp_psn->iFrameNum= -1;
  for(ii=0; ii<molecule->nRes; ++ii)
  {
    inp_psn->piNumResInteractions[ii] = 0;
    inp_psn->piNodeDegree[ii]         = 0;
    
    for(jj=0; jj<molecule->nRes; ++jj)
    {
      inp_psn->ppiIntAtomPairs[ii][jj]    = 0;
      inp_psn->ppfIntStrength[ii][jj]     = 0.0;
      inp_psn->ppfHubsIntStrength[ii][jj] = 0.0;
      inp_psn->ppfAvgResInt[ii][jj]       = 0.0;
    }
    
    for(jj=0; jj<inp_psn->iMaxResInteractions; ++jj)
      inp_psn->ppiInteractions[ii][jj] = 0;
  }
  
  for(ii=0; ii<inp_psn->iNumOfIntMinStep; ++ii)
  {
    for(jj=0; jj<molecule->nRes; ++jj)
    {
      inp_psn->ppiTmpStableHubs[ii][jj] = 0;
      inp_psn->ppiStableHubs[ii][jj]    = 0;
      
      for(kk=0; kk<molecule->nRes; ++kk)
      {
        inp_psn->pppiStableResInt[ii][jj][kk]         = 0;
        inp_psn->pppfStableResIntStrength[ii][jj][kk] = 0.0;
        
        for(mm=0; mm<3; ++mm)
          inp_psn->ppppiHubCorr[ii][jj][kk][mm] = 0;
      }
    }
    
    for(jj=0; jj<inp_psn->iMaxClustNum; ++jj)
      for(kk=0; kk<molecule->nRes; ++kk)
        inp_psn->pppiClusters[ii][jj][kk] = 0;
    
    for(jj=0; jj<2; ++jj)
      inp_psn->ppiLargestClusterSize[ii][jj] = 0;
  }
  
  for(ii=0; ii<molecule->nRes; ++ii)
    for(jj=0; jj<molecule->nRes; ++jj)
      inp_psn->ppiResResIntFreq[ii][jj] = 0;
  
  // ================================= //

  return;
}

double GetBigClsSize(FILE *FRawFile, int *piNodeClusters, int *piClustSize, float fImin, int iNumOfNodes, int **ppiIntMatrix)
{
  int           ii, jj;
  int           iNumOfFrames=0;
  int           iTmpRes1, iTmpRes2;
  int           iBigSize=0;
  int           iLastCluster=0;
  int           iAtmNum=0;
  
  float         fTmpRIS, fTmpHIS;
  double        fAvgBigSize=0.0;
  
  char          cLine[5000];
  char          cMagic[8];
  
  GetIminClearAll(iNumOfNodes, piNodeClusters, piClustSize, iLastCluster);
  iLastCluster=0;
  
  for(ii=0; ii<iNumOfNodes; ++ii)
    for(jj=0; jj<iNumOfNodes; ++jj)
      ppiIntMatrix[ii][jj] = 0;
  
  rewind(FRawFile);
  
  while(fgets(cLine, 5000, FRawFile) != NULL)
  {
    if(strncmp(cLine, ">INT", 4) == 0)
    {
      // skip header
     if(fgets(cLine, 5000, FRawFile)==NULL && ( !feof(FRawFile) || ferror(FRawFile) ))
     {
        fprintf(stderr, "Warning! Premature end of file reached!\n");
     }
      
      iNumOfFrames++;
      GetIminClearAll(iNumOfNodes, piNodeClusters, piClustSize, iLastCluster);
      iLastCluster=0;
      iBigSize = 0;
      
      while(1)
      {
        if(fgets(cLine, 5000, FRawFile)==NULL && ( !feof(FRawFile) || ferror(FRawFile) ))
        {
            fprintf(stderr, "Warning! Premature end of file reached!\n");
        }
        sscanf(cLine, "%s %d %d %f %f %d", cMagic, &iTmpRes1, &iTmpRes2, &fTmpRIS, &fTmpHIS, &iAtmNum);

        if(strcmp(cMagic, "&") != 0)
          break;
          
        iTmpRes1--;
        iTmpRes2--;

        if(fTmpRIS >= fImin)
        {
          iLastCluster = UpdateClusters(iTmpRes1, iTmpRes2, piNodeClusters, iNumOfNodes, iLastCluster);
          ppiIntMatrix[iTmpRes1][iTmpRes2]++;
        }
      }
      
      for(ii=0; ii<iNumOfNodes; ++ii)
        piClustSize[piNodeClusters[ii]]++;
      
      for(ii=1; ii<iLastCluster+1; ++ii)
      {
        if(iBigSize < piClustSize[ii])
          iBigSize = piClustSize[ii];
      }
      
      fAvgBigSize = fAvgBigSize + iBigSize;
    }
  }

  fAvgBigSize = (fAvgBigSize / (float) iNumOfFrames);
  fAvgBigSize = ((fAvgBigSize * 100) / (float) iNumOfNodes);

  return fAvgBigSize;
}

int UpdateClusters(int iTmpRes1, int iTmpRes2, int *piNodeClusters, int iNumOfNodes, int iLastCluster)
{
  int           ii;
  int           iTmpRes1Cls, iTmpRes2Cls;

  iTmpRes1Cls = piNodeClusters[iTmpRes1];
  iTmpRes2Cls = piNodeClusters[iTmpRes2];

  if(iTmpRes1Cls == iTmpRes2Cls && iTmpRes1Cls == 0)
  {
    iLastCluster++;
    piNodeClusters[iTmpRes1] = iLastCluster;
    piNodeClusters[iTmpRes2] = iLastCluster;
    return iLastCluster;
  }
  
  else if(iTmpRes1Cls == iTmpRes2Cls && iTmpRes1Cls != 0)
    return iLastCluster;
  
  else if(iTmpRes1Cls != iTmpRes2Cls)
  {
    if(iTmpRes1Cls == 0)
    {
      piNodeClusters[iTmpRes1] = iTmpRes2Cls;
      return iLastCluster;
    }
    
    if(iTmpRes2Cls == 0)
    {
      piNodeClusters[iTmpRes2] = iTmpRes1Cls;
      return iLastCluster;
    }
    
    for(ii=0; ii<iNumOfNodes; ii++)
      if(piNodeClusters[ii] == iTmpRes2Cls)
        piNodeClusters[ii] = iTmpRes1Cls;
    
    return iLastCluster;
  }
  
  return iLastCluster;
}

void GetIminClearAll(int iNumOfNodes, int *piNodeClusters, int *piClustSize, int iLastCluster)
{
  int           ii;

  for(ii=0; ii<iNumOfNodes; ++ii)
    piNodeClusters[ii] = 0;
  
  //for(ii=0; ii<1000; ++ii)
  for(ii=0; ii<iNumOfNodes; ++ii)
    piClustSize[ii] = 0;
}
// =====================================================================

int Read_iPSG ( char **input, int inp_index, struct inp_psn *inp_psn, char *printout, Molecule *molecule, int iNumOfFrames, struct sopt *OPT)
{
  int              ii, jj, kk, ww, oo;
  int              rescounter, atm_idx, found, seleflag=0;
  int              iTempRes1=0, iTempRes2=0, iInputFileStartingLine=0;
  int              iMergeResNum=0, iMergeDestRes=0;
  int              iTempInt;
  int              natmhere, gotit;
  
  char             bbatoms[50], buffer[10240], cIntMinRange[50];
  char             cCOOTerAtoms[50], cNHTerAtoms[50], cResName[6];
  char             cTmpString1[50], cTmpString2[50];
  char            *cToken, cMergeResType[5], cTempChArray[50];
  
  float            fTempFloat;
  
  extern short int     no_frame_par;
  no_frame_par = 1;
  
  memset ( buffer, '\0', sizeof(buffer));
  memset ( cTempChArray, '\0', sizeof(cTempChArray));
  
  iInputFileStartingLine = inp_index;
  
  // === Default Options ===
  sprintf(inp_psn->title, "PSNANALYSIS");
  inp_psn->fDistCutoff       = 4.5;
  inp_psn->fIntMinStart      = 0.0;
  inp_psn->fIntMinStop       = 0.0;
  inp_psn->fIntMinStep       = 1.0;
  inp_psn->iHubContCutoff    = 3;
  inp_psn->fStableCutoff     = 0.5;
  inp_psn->iWarningFlag      = 0;
  inp_psn->iVerboseFlag      = 0;
  inp_psn->iTerminiFlag      = 1;
  inp_psn->iProxCutOff       = 3;
  inp_psn->iNumOfNoLinkPairs = 0;
  inp_psn->iNumOfMerge       = 0;
  inp_psn->iMergeMaxResNum   = 0;
  inp_psn->iNumOfParam       = 0;
  inp_psn->iNumOfFrames      = iNumOfFrames;
  inp_psn->iIntType          = 0;
  inp_psn->iHubEqFlag        = 1;
  
  // === getimin section ===
  inp_psn->iGetIminFlag      =  0;
  inp_psn->iDecNum           =  2;
  inp_psn->fPCNSize          = 50.0;
  inp_psn->iSecondRound      =  0;
  // =======================
  
  while( strncmp (buffer, "END", 3))
  {
    gotit = 0;
    sprintf( buffer, "%s", input[inp_index]);
    if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
      gotit = 1;
    else if ( !strncmp(buffer, "--TITLE", 7))
    {
      sscanf( buffer, "--TITLE %s", inp_psn->title);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--SELE", 6))
    {
      sscanf(buffer, "--SELE %[^\n]%*c ", inp_psn->sele.selestring);
      seleflag = 1;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--DISTCUTOFF", 12))
    {
      sscanf( buffer, "--DISTCUTOFF %f", &inp_psn->fDistCutoff);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--STABLECUTOFF", 14))
    {
      sscanf( buffer, "--STABLECUTOFF %f", &inp_psn->fStableCutoff);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--HUBCONTCUTOFF", 15))
    {
      sscanf( buffer, "--HUBCONTCUTOFF %d", &inp_psn->iHubContCutoff);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--TERMINI", 9))
    {
      sscanf( buffer, "--TERMINI %d", &inp_psn->iTerminiFlag);
      if(inp_psn->iTerminiFlag!=1)
        inp_psn->iTerminiFlag=0;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--PROXIMITY", 11))
    {
      sscanf( buffer, "--PROXIMITY %d", &inp_psn->iProxCutOff);
      if(inp_psn->iProxCutOff<1)
      {
        fprintf( stderr, "PSN module: --PROXIMITY must be an integer >= 1\n");
        exit(5);
      }
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--VERBOSE", 9))
    {
      sscanf( buffer, "--VERBOSE %d", &inp_psn->iVerboseFlag);
      if(inp_psn->iVerboseFlag!=1)
        inp_psn->iVerboseFlag=0;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--INTTYPE", 9))
    {
      sscanf( buffer, "--INTTYPE %s", cTempChArray);
      
      if(strcmp(cTempChArray, "SC") == 0)
        inp_psn->iIntType = 0;

      else if(strcmp(cTempChArray, "SC+CA") == 0)
        inp_psn->iIntType = 1;

      else if(strcmp(cTempChArray, "ALL") == 0)
        inp_psn->iIntType = 2;
        
      else
      {
        fprintf( stderr, "PSN module: invalid --INTTYPE value ''%s'', valid values are: SC, SC+CA, ALL\n", cTempChArray);
        exit(5);
      }
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--HUBEQ", 7))
    {
      sscanf(buffer, "--HUBEQ %s", cTempChArray);
      
      if(strcmp(cTempChArray, "YES") == 0)
        inp_psn->iHubEqFlag = 1;

      else if(strcmp(cTempChArray, "NO") == 0)
        inp_psn->iHubEqFlag = 0;

      else
      {
        fprintf( stderr, "PSN module: invalid --HUBEQ value ''%s'', valid values are: YES, NO\n", cTempChArray);
        exit(5);
      }
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--NOLINK", 8))
    {
      inp_psn->iNumOfNoLinkPairs++;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--PARAM", 7))
    {
      inp_psn->iNumOfParam++;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--MERGE", 7))
    {
      inp_psn->iNumOfMerge++;
      buffer[strlen(buffer)-1] = '\0';
      
      // This is --MERGE
      cToken = strtok(buffer, " ");
      
      // This is New Residu Type
      cToken = strtok(NULL, " ");
      
      // This is Destination Residue
      cToken = strtok(NULL, " ");
      
      cToken = strtok(NULL, " ");
      
      iMergeResNum = 0;
      
      while(cToken != NULL)
      {

        if(strcmp(cToken, "\n") != 0)
          iMergeResNum++;
        
        cToken = strtok(NULL, " ");
      }
      
      if(iMergeResNum > inp_psn->iMergeMaxResNum)
        inp_psn->iMergeMaxResNum = iMergeResNum;
      gotit = 1;
    }
    
    else if ( !strncmp(buffer, "--INTMIN", 8))
    {
      sscanf( buffer, "--INTMIN %[^\n]%*c", cIntMinRange);
      if(strchr(cIntMinRange, ':')!=NULL)
      {
        sscanf(cIntMinRange, "%[a-z-A-Z]:%f:%d", cTmpString1, &fTempFloat, &iTempInt);
        if(strcmp(cTmpString1, "AUTO") == 0)
        {
          // automagic imin identification (aka GetImin)
          inp_psn->iGetIminFlag = 1;
          inp_psn->iSecondRound = 0;
          inp_psn->fPCNSize     = fTempFloat;
          inp_psn->iDecNum      = iTempInt;
          inp_psn->fIntMinStart = 0.0;
          inp_psn->fIntMinStop  = 1.0;
          inp_psn->fIntMinStep  = 1.0;
        }
        
        else
          sscanf(cIntMinRange, "%f:%f:%f", &inp_psn->fIntMinStart, &inp_psn->fIntMinStop, &inp_psn->fIntMinStep);
      }
      
      else
      {

        if(strcmp(cIntMinRange, "AUTO") == 0)
        {
          inp_psn->iGetIminFlag = 1;
          inp_psn->iSecondRound = 0;
          inp_psn->fPCNSize     = 50.0;
          inp_psn->iDecNum      = 2;
          inp_psn->fIntMinStart = 0.0;
          inp_psn->fIntMinStop  = 1.0;
          inp_psn->fIntMinStep  = 1.0;
        }
        
        else
        {
          inp_psn->fIntMinStart = atof(cIntMinRange);
          inp_psn->fIntMinStop=inp_psn->fIntMinStart;
          inp_psn->fIntMinStep=1.0;
        }
      }
      
      strcpy(inp_psn->cIminRange, cIntMinRange);
      gotit = 1;
    }
    
    if( gotit==0 )
    {
      fprintf( stderr, "PSN module: Could NOT understand option: %s\n", buffer);
      exit(5);
    }
    inp_index++;
  }
  
  if (!seleflag)
  {
    fprintf( stderr, "PSN module: You did not supply a selection\n");
    exit(99);
  }
  
  // === GetImin Check ===
  if(inp_psn->iGetIminFlag == 1 && inp_psn->iVerboseFlag == 0)
    inp_psn->iVerboseFlag = 1;
  // =====================
  
  //tot_residues = molecule->nRes;
  inp_psn->tot_nresidues = 0;
  for( ii=0; ii<molecule->nSeg ; ii++ )
   for( jj=0; jj<molecule->segment[ii].nRpS ; jj++ )
    inp_psn->tot_nresidues++;
  inp_psn->tot_nresidues = molecule->nRes;
  inp_psn->iNumOfNodes   = molecule->nRes;
  
  // === NoLink Section ================================================
  if(inp_psn->iNumOfNoLinkPairs != 0)
  {
    inp_psn->ppiNoLinkPairs = calloc(inp_psn->iNumOfNoLinkPairs, sizeof(int *));
    for(ii=0; ii<inp_psn->iNumOfNoLinkPairs; ii++)
      inp_psn->ppiNoLinkPairs[ii] = calloc(2, sizeof(int));
    
    inp_index = iInputFileStartingLine;
    sprintf( buffer, "%s", input[inp_index]);
    ii = -1;
    while( strncmp (buffer, "END", 3))
    {
      sprintf( buffer, "%s", input[inp_index]);
      if ( !strncmp(buffer, "--NOLINK", 8))
      {
        ii++;
        sscanf(buffer, "--NOLINK %d %d", &iTempRes1, &iTempRes2);
        
        if(iTempRes1 > inp_psn->tot_nresidues)
        {
          fprintf( stderr, "PSN module: Invalid residue number in --NOLINK option : %d\n", iTempRes1);
          exit(99);
        }

        if(iTempRes2 > inp_psn->tot_nresidues)
        {
          fprintf( stderr, "PSN module: Invalid residue number in --NOLINK option : %d\n", iTempRes2);
          exit(99);
        }
        
        inp_psn->ppiNoLinkPairs[ii][0] = iTempRes1 - 1;
        inp_psn->ppiNoLinkPairs[ii][1] = iTempRes2 - 1;
      }
      inp_index++;
    }
  }
  // ===================================================================

  // === Param Section =================================================
  if(inp_psn->iNumOfParam != 0)
  {
    inp_psn->pcParamResVect = (char **) calloc(inp_psn->iNumOfParam, sizeof(char *));
    for(ii=0; ii<inp_psn->iNumOfParam; ii++)
      inp_psn->pcParamResVect[ii] = (char *) calloc(50, sizeof(char));

    inp_psn->pfParamValVect = (float *) calloc(inp_psn->iNumOfParam, sizeof(float));
    
    inp_index = iInputFileStartingLine;
    sprintf( buffer, "%s", input[inp_index]);
    ii = -1;
    while( strncmp (buffer, "END", 3))
    {
      sprintf( buffer, "%s", input[inp_index]);
      if ( !strncmp(buffer, "--PARAM", 7))
      {
        ii++;
        sscanf(buffer, "--PARAM %s %f", cTmpString1, &fTempFloat);
        
        if(fTempFloat <= 0.0)
        {
          fprintf( stderr, "PSN module: Normalization factors in --PARAM option must be >= 0.0, given %f\n", fTempFloat);
          exit(99);
        }
        
        strcpy(inp_psn->pcParamResVect[ii], cTmpString1);
        inp_psn->pfParamValVect[ii] = fTempFloat;
      }
      inp_index++;
    }
  }
  // ===================================================================
  
  // take note of atoms to consider for interaction calculations
  
  if(inp_psn->iIntType == 0)
  {
    // SC atoms only
    sprintf( bbatoms, "@(CA|C|N|O|O*T?|H*)");
    if(inp_psn->iTerminiFlag==1)
    {
      sprintf(cCOOTerAtoms, "@(C|O*T?)");
      sprintf(cNHTerAtoms,  "N");
    }
  }

  else if(inp_psn->iIntType == 1)
  {
    // SC + CA atoms
    sprintf( bbatoms, "@(C|N|O|O*T?|H*)");
    if(inp_psn->iTerminiFlag==1)
    {
      sprintf(cCOOTerAtoms, "@(C|O*T?)");
      sprintf(cNHTerAtoms,  "N");
    }
  }
  
  else if(inp_psn->iIntType == 2)
  {
    // ALL atoms
    sprintf( bbatoms, "@(H*)");
    if(inp_psn->iTerminiFlag==1)
    {
      sprintf(cCOOTerAtoms, "@(C|O*T?)");
      sprintf(cNHTerAtoms,  "N");
    }
  }
  
  // sprintf( bbatoms, "@(CA|C|N|O|H|OT1|OT2|OXT|H*)");
  
  inp_psn->reslength = calloc( molecule->nRes, sizeof(int));
  inp_psn->reslength2 = calloc( molecule->nRes, sizeof(int));
  
  inp_psn->atmList = calloc( molecule->nRes, sizeof(int *));
  for( ii=0; ii<molecule->nRes; ii++ )
    inp_psn->atmList[ii] = calloc(500, sizeof(int)); // 500 as maximum number of atoms in a residue
  
  inp_psn->res_norm = calloc( molecule->nRes, sizeof(float));
  atm_idx = 0;
  rescounter = 0;
  natmhere = 0;
  for( ii=0; ii<molecule->nRes ; ii++ )
  {
    fTempFloat = GetUserNormFactor(inp_psn, molecule->rawmol.restype[atm_idx]);
    if(fTempFloat < 0)
      fTempFloat = fGetNormFactor(molecule->rawmol.restype[atm_idx]);
    inp_psn->res_norm[ii] = fTempFloat;

    if(inp_psn->res_norm[ii]>999)
      inp_psn->iWarningFlag=1;
    
    for( kk=0; kk<molecule->nApR[ii] ; kk++ )  // run along atoms     
    {
      if(fnmatch(bbatoms, molecule->rawmol.atmtype[atm_idx+kk], FNM_EXTMATCH))
      {
        inp_psn->atmList[rescounter][natmhere] = molecule->rawmol.atomn[atm_idx+kk];
        natmhere++;
      }
      
      if(inp_psn->iTerminiFlag==1)
      {
        if( molecule->rawmol.segend[atm_idx+molecule->nApR[ii]-1] == 1 ) 
          if(!fnmatch(cCOOTerAtoms, molecule->rawmol.atmtype[atm_idx+kk], FNM_EXTMATCH))
          {
            inp_psn->atmList[rescounter][natmhere] = molecule->rawmol.atomn[atm_idx+kk];
            natmhere++;
          }

        if( molecule->rawmol.segbeg[atm_idx] == 1 )
          if(!fnmatch(cNHTerAtoms, molecule->rawmol.atmtype[atm_idx+kk], FNM_EXTMATCH))
          {
            inp_psn->atmList[rescounter][natmhere] = molecule->rawmol.atomn[atm_idx+kk];
            natmhere++;
          }
      }
    }       
    
    inp_psn->reslength[rescounter] = natmhere;
    inp_psn->reslength2[rescounter] = molecule->nApR[ii];
    rescounter++;
    natmhere = 0;
    atm_idx += molecule->nApR[ii];
  }
  
  // if SELE, then not selected residues'first atoms are set to -1 and residues neglected in computation
  GetSele( inp_psn->sele.selestring, &inp_psn->sele, molecule );
  if( inp_psn->sele.nselatm == 0)
  {
    fprintf( stderr, "Selected Atoms = 0!\n");
    exit(1);
  }
  for( ii=0; ii<inp_psn->tot_nresidues; ii++)
  {
    if( inp_psn->reslength[ii] == 0 )
      continue;
    found = 0;
    for( jj=0; jj<inp_psn->sele.nselatm; jj++ )
      if( inp_psn->atmList[ii][0] == inp_psn->sele.selatm[jj] )
      {
        found = 1;
        break;
      }
    if (!found)
      inp_psn->atmList[ii][0] = -1;
  }
  
  if(inp_psn->iVerboseFlag==1)
  {
    sprintf(inp_psn->cRawFileName, "%s%s", "raw", inp_psn->title);
    inp_psn->outfile = O_File(inp_psn->cRawFileName, "w");
  }
  
  inp_psn->cluster.cluspop = calloc( molecule->nRes, sizeof(int));
  inp_psn->cluster.clusters = calloc( molecule->nRes, sizeof(int *));
  for( ii=0; ii<molecule->nRes; ii++)
    inp_psn->cluster.clusters[ii] = calloc( molecule->nRes, sizeof(int));

  inp_psn->iMaxNumOfLink=0;
  
  // cast uses floor: adding 1.5 to actually be sure to be adding 1
  inp_psn->iNumOfIntMinStep=(int)(((inp_psn->fIntMinStop - inp_psn->fIntMinStart)/inp_psn->fIntMinStep)+1.5);

  /*
    ppiIntAtomPairs allocation
    This is a NumOfRes x NumOfRes matrix that stores the
    the number of interacting atom pairs at a given frame
  */
  inp_psn->ppiIntAtomPairs=(int **) calloc(molecule->nRes, sizeof (int *));
  for(ii=0; ii<molecule->nRes; ii++)
    inp_psn->ppiIntAtomPairs[ii] = (int *) calloc(molecule->nRes, sizeof(int));

  /*
    ppcIntAtomNamePairs allocation
    This is a NumOfRes x NumOfRes matrix that stores the
    the names of interacting atoms at a given frame
  */
  inp_psn->pppcIntAtomNamePairs=(char ***) calloc(molecule->nRes, sizeof (char **));
  for(ii=0; ii<molecule->nRes; ii++)
  {
    inp_psn->pppcIntAtomNamePairs[ii] = (char **) calloc(molecule->nRes, sizeof(char *));
    for(jj=0; jj<molecule->nRes; jj++)
    {
      inp_psn->pppcIntAtomNamePairs[ii][jj] = (char *) calloc(1000, sizeof(char));
    }
  }

  /*
    ppfIntStrength allocation
    This is a NumOfRes x NumOfRes matrix that stores the
    Interaction Strengths of all Res-Res at a given frame
  */  
  inp_psn->ppfIntStrength=(float **) calloc(molecule->nRes, sizeof (float *));
  for(ii=0; ii<molecule->nRes; ii++)
    inp_psn->ppfIntStrength[ii] = (float *) calloc(molecule->nRes, sizeof(float));


  /*
   ppfHubsIntStrength allocation
   This is a nRes x nRes matrix that stores Res-Res I.S.
   used for hubs identification at a given frame
  */
  
  inp_psn->ppfHubsIntStrength=(float **) calloc(molecule->nRes, sizeof (float *));
  for(ii=0; ii<molecule->nRes; ii++)
    inp_psn->ppfHubsIntStrength[ii] = (float *) calloc(molecule->nRes, sizeof(float));
  
  inp_psn->iFrameNum= -1;


  // The highest number of residues interactions
  inp_psn->iMaxResInteractions=100*(int)inp_psn->fDistCutoff;
  
  /*
    ppiInteractions allocation
    The Map of Residue interaction at a given frame at a given Imin
  
  */
  inp_psn->ppiInteractions = (int **) calloc(molecule->nRes, sizeof (int *));
  for(ii=0; ii<molecule->nRes; ii++)
    inp_psn->ppiInteractions[ii] = (int *) calloc(inp_psn->iMaxResInteractions, sizeof(int));

  /*
    piNumResInteractions allocation
    Number of interactions of all residues at a given frame at a given Imin
  */
  inp_psn->piNumResInteractions=(int *) calloc(molecule->nRes, sizeof(int));
  
  /* ======================================================================= */
  
  /*
    ppiTmpStableHubs allocation
    This is a iNumOfRes x iNumOfRes matrix
    used to store averaged res-res interaction strength
  */
  inp_psn->ppiTmpStableHubs=(int **) calloc(inp_psn->iNumOfIntMinStep, sizeof(int *));
  for(ii=0; ii<inp_psn->iNumOfIntMinStep; ii++)
  {
    inp_psn->ppiTmpStableHubs[ii]=(int *) calloc(molecule->nRes, sizeof(int));
  }
  
  /*
    ppiStableHubs allocation
    This is a iNumOfRes x iNumOfRes matrix
    used to store averaged res-res interaction strength
  */
  inp_psn->ppiStableHubs=(int **) calloc(inp_psn->iNumOfIntMinStep, sizeof(int *));
  for(ii=0; ii<inp_psn->iNumOfIntMinStep; ii++)
  {
    inp_psn->ppiStableHubs[ii]=(int *) calloc(molecule->nRes, sizeof(int));
  }

  /*
    ppfAvgResInt allocation
    This is a iNumOfRes x iNumOfRes matrix
    used to store averaged res-res interaction strength
  */
  inp_psn->ppfAvgResInt=(float **) calloc(molecule->nRes, sizeof(float *));
  for(ii=0; ii<molecule->nRes; ii++)
  {
    inp_psn->ppfAvgResInt[ii]=(float *) calloc(molecule->nRes, sizeof(float));
  }

  /*
    ppiResResIntFreq allocation
    This is a iNumOfRes x iNumOfRes matrix
    used to store the number of frames in which two res have
    interaction strength != 0.0
  */
  inp_psn->ppiResResIntFreq=(int **) calloc(molecule->nRes, sizeof(int *));
  for(ii=0; ii<molecule->nRes; ii++)
    inp_psn->ppiResResIntFreq[ii]=(int *) calloc(molecule->nRes, sizeof(int));
  
  for(ii=0; ii<molecule->nRes; ii++)
    for(jj=0; jj<molecule->nRes; jj++)
      inp_psn->ppiResResIntFreq[ii][jj] = 0;

  /*
     pppiStableResInt allocation
     used to store the stable Residue Interactions
  */
  inp_psn->pppiStableResInt=(int ***) calloc(inp_psn->iNumOfIntMinStep, sizeof(int **));
  for(ii=0; ii<inp_psn->iNumOfIntMinStep; ii++)
  {
    inp_psn->pppiStableResInt[ii]=(int **) calloc(molecule->nRes, sizeof(int *));
    for(jj=0; jj<molecule->nRes; jj++)
    {
      inp_psn->pppiStableResInt[ii][jj]=(int *) calloc(molecule->nRes, sizeof(int));
    }
  }

  /*
     pppfStableResIntStrength allocation
     used to store stable Residue Interaction Strength
  */
  inp_psn->pppfStableResIntStrength=(float ***) calloc(inp_psn->iNumOfIntMinStep, sizeof(float **));
  for(ii=0; ii<inp_psn->iNumOfIntMinStep; ii++)
  {
    inp_psn->pppfStableResIntStrength[ii]=(float **) calloc(molecule->nRes, sizeof(float *));
    for(jj=0; jj<molecule->nRes; jj++)
    {
      inp_psn->pppfStableResIntStrength[ii][jj]=(float *) calloc(molecule->nRes, sizeof(float));
    }
  }

  /*
    ppppiHubCorr allocation
    used to store Hub Correlations
  */
  
  inp_psn->ppppiHubCorr=(int ****) calloc(inp_psn->iNumOfIntMinStep, sizeof(int ***));
  for(ii=0; ii<inp_psn->iNumOfIntMinStep; ii++)
  {
    inp_psn->ppppiHubCorr[ii]=(int ***) calloc(molecule->nRes, sizeof(int **));
    for(jj=0; jj<molecule->nRes; jj++)
    {
      inp_psn->ppppiHubCorr[ii][jj]=(int **) calloc(molecule->nRes, sizeof(int *));
      for(kk=0; kk<molecule->nRes; kk++)
      {
        inp_psn->ppppiHubCorr[ii][jj][kk]=(int *) calloc(3, sizeof(int));
      }
    }
  }

  /*
     pppiClusters allocation
     used to store the Stable Cluster Compositions
  */
  inp_psn->iMaxClustNum=(int)(((float)molecule->nRes/2.0)+1);
  inp_psn->pppiClusters=(int ***) calloc(inp_psn->iNumOfIntMinStep, sizeof(int **));
  for(ii=0; ii<inp_psn->iNumOfIntMinStep; ii++)
  {
    inp_psn->pppiClusters[ii]=(int **) calloc(inp_psn->iMaxClustNum, sizeof(int *));
    for(jj=0; jj<inp_psn->iMaxClustNum; jj++)
    {
      inp_psn->pppiClusters[ii][jj]=(int *) calloc(molecule->nRes, sizeof(int));
    }
  }
  
  /*
     piLargestClusterSize Allocation
     This vector stores the averaged size of largest cluster at each Imin step
  */
  
  inp_psn->ppiLargestClusterSize=(int **) calloc(inp_psn->iNumOfIntMinStep, sizeof(int*));
  for(ii=0; ii<inp_psn->iNumOfIntMinStep; ii++)
  {
    inp_psn->ppiLargestClusterSize[ii]=(int *) calloc(2, sizeof(int));
  }
  
  /*
     piNodeDegree Allocation
     This vector stores the degree of each node
  */
  
  inp_psn->piNodeDegree = (int *) calloc(molecule->nRes, sizeof(int));
  
  // ===================================================================
  
  
  // === Merge Section =================================================
  if(inp_psn->iNumOfMerge != 0)
  {
    inp_psn->iMergeMaxResNum++;
    
    inp_psn->ppcMergedResidues = calloc(inp_psn->iNumOfMerge, sizeof(char *));
    for(ii=0; ii<inp_psn->iNumOfMerge; ii++)
      inp_psn->ppcMergedResidues[ii] = calloc(500, sizeof(char));

    //inp_psn->ppiMergedResidues = calloc(inp_psn->iNumOfMerge, sizeof(int *));
    //for(ii=0; ii<inp_psn->iNumOfMerge; ii++)
    //{
      //inp_psn->ppiMergedResidues[ii] = calloc(inp_psn->iMergeMaxResNum, sizeof(int));
      //for(jj=0; jj<inp_psn->iMergeMaxResNum; jj++)
      //{
        //inp_psn->ppiMergedResidues[ii][jj] = -1;
      //}
    //}
    
    inp_index = iInputFileStartingLine;
    sprintf( buffer, "%s", input[inp_index]);
    
    ww = -1;
    while( strncmp (buffer, "END", 3))
    {
      sprintf( buffer, "%s", input[inp_index]);
      if ( !strncmp(buffer, "--MERGE", 7))
      {
        ww++;
        buffer[strlen(buffer)-1] = '\0';
        strncpy(inp_psn->ppcMergedResidues[ww], buffer+8, strlen(buffer));
        
        // This is --MERGE
        cToken = strtok(buffer, " ");
        
        // This is New Residu Type
        cToken = strtok(NULL, " ");
        sscanf(cToken, "%d", &iMergeDestRes);
        iMergeDestRes--;
        
        // This is Destination Residue
        cToken = strtok(NULL, " ");
        sscanf(cToken, "%s", cMergeResType);
        
        fTempFloat = GetUserNormFactor(inp_psn, cMergeResType);
        if(fTempFloat < 0)
          fTempFloat = fGetNormFactor(cMergeResType);
        inp_psn->res_norm[iMergeDestRes] = fTempFloat;

        cToken = strtok(NULL, " ");
        
        //kk = 0;
        while(cToken != NULL)
        {
          if(strcmp(cToken, "\n") != 0)
          {
            sscanf(cToken, "%d", &iMergeResNum);
            iMergeResNum--;
            if(iMergeResNum != iMergeDestRes)
            {
              for(oo=0; oo<inp_psn->reslength[iMergeResNum]; oo++)
              {
                inp_psn->atmList[iMergeDestRes][inp_psn->reslength[iMergeDestRes]+oo] = inp_psn->atmList[iMergeResNum][oo];
                inp_psn->atmList[iMergeResNum][oo] = -1;
              }
              
              inp_psn->reslength[iMergeDestRes] = inp_psn->reslength[iMergeDestRes] + inp_psn->reslength[iMergeResNum];  
              inp_psn->reslength[iMergeResNum] = 0;
            }
          }
          cToken = strtok(NULL, " ");
        }
      }
      inp_index++;
    }
  }
  // ===================================================================
  
  cTmpString1[0]='\0';
  cTmpString2[0]='\0';
  inp_psn->iNumOfSelRes = 0;
  for(ii=0; ii<inp_psn->sele.nselatm; ii++)
  {
    sprintf(cTmpString1, "%s%d", molecule->rawmol.segId[inp_psn->sele.selatm[ii]-1], molecule->rawmol.resn[inp_psn->sele.selatm[ii]-1]);
    if(strcmp(cTmpString1, cTmpString2)!=0)
    {
      inp_psn->iNumOfSelRes++;
      strcpy(cTmpString2, cTmpString1);
    }
  }
  
  time(&inp_psn->time_tToday);
  
  
  if(inp_psn->iVerboseFlag==1)
  {
    fprintf( inp_psn->outfile, "# ======================================================================\n");
    fprintf( inp_psn->outfile, "#                       --- Protein Structure Network Analysis ---\n");
    fprintf( inp_psn->outfile, "#\n");
    fprintf( inp_psn->outfile, "# Copyright     : Francesca Fanelli\n");
    fprintf( inp_psn->outfile, "#                 Angelo Felline,\n");
    fprintf( inp_psn->outfile, "#                 Michele Seeber (2009)\n");
    fprintf( inp_psn->outfile, "# License       : GPL v 3\n");
    fprintf( inp_psn->outfile, "#\n");
    fprintf( inp_psn->outfile, "# PSN Ver       : 0.8c\n");
    fprintf( inp_psn->outfile, "#\n");
    fprintf( inp_psn->outfile, "# Date          : %s", asctime(localtime(&inp_psn->time_tToday)));
    fprintf( inp_psn->outfile, "#\n");
    fprintf( inp_psn->outfile, "# Mol Name      : %s\n", OPT->IMOL_FILE);
    fprintf( inp_psn->outfile, "# Seg Num       : %d\n", molecule->nSeg);
    fprintf( inp_psn->outfile, "# Res Num       : %d\n", inp_psn->tot_nresidues);
    fprintf( inp_psn->outfile, "# Trj Name      : %s\n", OPT->ITRJ_FILE);
    fprintf( inp_psn->outfile, "# Frame Num     : %d\n", iNumOfFrames);
    if(inp_psn->iPBCFlag == 0)
      fprintf(inp_psn->outfile, "# PBC           : No\n");
    else
      fprintf(inp_psn->outfile, "# PBC           : Yes\n");
    fprintf( inp_psn->outfile, "#\n");
    fprintf( inp_psn->outfile, "# Title         : %s\n", inp_psn->title);
    fprintf( inp_psn->outfile, "# Selection     : %s\n", inp_psn->sele.selestring);
    fprintf( inp_psn->outfile, "# Sele Res Num  : %d\n", inp_psn->iNumOfSelRes);
    fprintf( inp_psn->outfile, "# DistCutOff    : %f\n", inp_psn->fDistCutoff);
    fprintf( inp_psn->outfile, "# IntMin        : %s\n", inp_psn->cIminRange);
    fprintf( inp_psn->outfile, "# IntMinStart   : %f\n", inp_psn->fIntMinStart);
    fprintf( inp_psn->outfile, "# IntMinStop    : %f\n", inp_psn->fIntMinStop);
    fprintf( inp_psn->outfile, "# IntMinStep    : %f\n", inp_psn->fIntMinStep);
    fprintf( inp_psn->outfile, "# StableCutOff  : %f\n", inp_psn->fStableCutoff);
    if(inp_psn->iHubEqFlag == 0)
      fprintf( inp_psn->outfile, "# HubEquation   : NO\n");
    else
      fprintf( inp_psn->outfile, "# HubEquation   : YES\n");
    fprintf( inp_psn->outfile, "# HubContCutOff : %i\n", inp_psn->iHubContCutoff);
    fprintf( inp_psn->outfile, "# Termini       : %i\n", inp_psn->iTerminiFlag);
    fprintf( inp_psn->outfile, "# Proximity     : %i\n", inp_psn->iProxCutOff);
    
    if(inp_psn->iIntType == 0)
      fprintf( inp_psn->outfile, "# Int Type      : SC\n");
    else if(inp_psn->iIntType == 1)
      fprintf( inp_psn->outfile, "# Int Type      : SC+CA\n");
    else if(inp_psn->iIntType == 2)
      fprintf( inp_psn->outfile, "# Int Type      : ALL\n");
    fprintf( inp_psn->outfile, "# Int Type      : %d\n", inp_psn->iIntType);
    
    fprintf( inp_psn->outfile, "#\n");
    fprintf( inp_psn->outfile, "# Num Of NoLink : %d\n", inp_psn->iNumOfNoLinkPairs);
    for(ii=0; ii<inp_psn->iNumOfNoLinkPairs; ii++)
      fprintf( inp_psn->outfile, "# No Link       : %d %d   %d\n", ii+1, inp_psn->ppiNoLinkPairs[ii][0]+1, inp_psn->ppiNoLinkPairs[ii][1]+1);
    
    fprintf(inp_psn->outfile, "# Num Of Merge  : %d\n", inp_psn->iNumOfMerge);
    for(ii=0; ii<inp_psn->iNumOfMerge; ii++)
      fprintf(inp_psn->outfile, "# Merge         : %d %s\n", ii+1, inp_psn->ppcMergedResidues[ii]);
    
    fprintf(inp_psn->outfile, "# Num Of Param  : %d\n", inp_psn->iNumOfParam);
    for(ii=0; ii<inp_psn->iNumOfParam; ii++)
      fprintf(inp_psn->outfile, "# Param         : %d %5s %f\n", ii+1, inp_psn->pcParamResVect[ii], inp_psn->pfParamValVect[ii]);
    
    fprintf( inp_psn->outfile, "# ======================================================================\n");
    fprintf( inp_psn->outfile, "\n");
    
    fprintf( inp_psn->outfile, "*** Seg Info ***\n");
    jj=1;
    for(ii=0; ii<molecule->nSeg; ii++)
    {
      fprintf( inp_psn->outfile, "%2d   %5s   %6d   %10d\n", ii+1, molecule->segment[ii].segName, molecule->segment[ii].nRpS, jj);
      jj=jj+molecule->segment[ii].nRpS;
    }
    fprintf( inp_psn->outfile, "========================\n");
    
    fprintf( inp_psn->outfile, "*** Seq ***\n");
    
    kk=0;
    for(ii=0; ii<molecule->nSeg; ii++)
    {
      for(jj=0; jj<molecule->segment[ii].nRpS; jj++)
      {
        kk++;
        if     (strncmp(molecule->segment[ii].pRes[jj].resType, "ALA", 3)==0) strcpy(cResName, "A\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "ARG", 3)==0) strcpy(cResName, "R\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "ASN", 3)==0) strcpy(cResName, "N\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "ASP", 3)==0) strcpy(cResName, "D\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "CYS", 3)==0) strcpy(cResName, "C\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "CYN", 3)==0) strcpy(cResName, "C\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "GLU", 3)==0) strcpy(cResName, "E\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "GLN", 3)==0) strcpy(cResName, "Q\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "GLY", 3)==0) strcpy(cResName, "G\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "HIS", 3)==0) strcpy(cResName, "H\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "HSC", 3)==0) strcpy(cResName, "H\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "HSD", 3)==0) strcpy(cResName, "H\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "HSP", 3)==0) strcpy(cResName, "H\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "HID", 3)==0) strcpy(cResName, "H\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "ILE", 3)==0) strcpy(cResName, "I\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "LEU", 3)==0) strcpy(cResName, "L\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "LYS", 3)==0) strcpy(cResName, "K\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "LYP", 3)==0) strcpy(cResName, "K\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "MET", 3)==0) strcpy(cResName, "M\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "PHE", 3)==0) strcpy(cResName, "F\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "PRO", 3)==0) strcpy(cResName, "P\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "SER", 3)==0) strcpy(cResName, "S\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "THR", 3)==0) strcpy(cResName, "T\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "TRP", 3)==0) strcpy(cResName, "W\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "TYR", 3)==0) strcpy(cResName, "Y\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "VAL", 3)==0) strcpy(cResName, "V\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "RET", 3)==0) strcpy(cResName, "X\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "GDP", 3)==0) strcpy(cResName, "d\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "GTP", 3)==0) strcpy(cResName, "t\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "Mg+", 3)==0) strcpy(cResName, "m\0");
        else if(strncmp(molecule->segment[ii].pRes[jj].resType, "ZMA", 3)==0) strcpy(cResName, "Z\0");
        else                                                                  strcpy(cResName, "U\0");
        
        fprintf(inp_psn->outfile, "%10d   %5s:%s%-9d   %.5f\n", kk, molecule->segment[ii].segName, cResName, molecule->segment[ii].pRes[jj].resn, inp_psn->res_norm[kk-1]);
      }
    }
    fprintf( inp_psn->outfile, "========================\n");
  }
  return 0;
}
// ------------------------------------------------------------------
int Compute_PSG ( struct inp_psn *inp_psn, struct sopt *OPT, CoorSet *trj_crd, char *outprint, Molecule *molecule )
{
  int                   ii, jj, kk, ll, ww, yy;
  int                   iTempRes1, iTempRes2;
  int                   resnumber1, resnumber2, interacting;
  int                   iClustSize;
  long double           dist;   
  float                 fIntMinIterator;
  float                 fIntStrength;
  int                   iIntMinIterationNum, iContactsCounter;
  int                   res1, res2;
  char                  res1_atmname[10], res2_atmname[10], cDistRep[100];

  inp_psn->iFrameNum++;  // update the number of actual frame

  if(inp_psn->iVerboseFlag==1)
    fprintf( inp_psn->outfile, "nFr: %d\n", inp_psn->iFrameNum+1);
  
  for(ii=0; ii<inp_psn->tot_nresidues; ii++)
  {
    for(jj=0; jj<inp_psn->tot_nresidues; jj++)
    {
      inp_psn->ppfIntStrength[ii][jj]          = 0.0;
      inp_psn->ppfHubsIntStrength[ii][jj]      = 0.0;
      inp_psn->ppiIntAtomPairs[ii][jj]         =   0;
      inp_psn->pppcIntAtomNamePairs[ii][jj][0] = '\0';
    }
  }
  
  /*
  for(ii=0; ii<inp_psn->tot_nresidues; ii++)
  {
    for(jj=0; jj<inp_psn->tot_nresidues; jj++)
    {
      inp_psn->ppfHubsIntStrength[ii][jj]=0.0;
    }
  }
  */

  // compute I.S. for each residue pair and place it in inp_psn->ppfIntStrength
  for(ii=0; ii<inp_psn->tot_nresidues; ii++)
  {
    if(inp_psn->atmList[ii][0] == -1 )
      continue;
    
    for( jj=ii+inp_psn->iProxCutOff; jj<inp_psn->tot_nresidues; jj++)
    {
      if(inp_psn->atmList[jj][0] == -1 )
        continue;
        
      interacting = 0;
      //toofar = 0;
      
      if(trj_crd->pbc_flag == 0 || trj_crd->pbc_flag == -1)
      {
        // No PBC 
        for( kk=0; kk<inp_psn->reslength[ii]; kk++)
        {
          for( ll=0; ll<inp_psn->reslength[jj]; ll++)
          {
            
            dist = sqrt( (trj_crd->xcoor[inp_psn->atmList[ii][kk] - 1] - trj_crd->xcoor[inp_psn->atmList[jj][ll] - 1])*(trj_crd->xcoor[inp_psn->atmList[ii][kk] - 1] - trj_crd->xcoor[inp_psn->atmList[jj][ll] - 1]) +
                         (trj_crd->ycoor[inp_psn->atmList[ii][kk] - 1] - trj_crd->ycoor[inp_psn->atmList[jj][ll] - 1])*(trj_crd->ycoor[inp_psn->atmList[ii][kk] - 1] - trj_crd->ycoor[inp_psn->atmList[jj][ll] - 1]) +
                         (trj_crd->zcoor[inp_psn->atmList[ii][kk] - 1] - trj_crd->zcoor[inp_psn->atmList[jj][ll] - 1])*(trj_crd->zcoor[inp_psn->atmList[ii][kk] - 1] - trj_crd->zcoor[inp_psn->atmList[jj][ll] - 1]) );
            
            if( dist <= inp_psn->fDistCutoff )
            {
              interacting++;
              
              res1_atmname[0] = '\0';
              res2_atmname[0] = '\0';
              cDistRep    [0] = '\0';
              
              strcat(res1_atmname, molecule->rawmol.atmtype[inp_psn->atmList[ii][kk] - 1]);
              strcat(res2_atmname, molecule->rawmol.atmtype[inp_psn->atmList[jj][ll] - 1]);
              sprintf(cDistRep, "%.3Lf", dist);
              
              strcat(inp_psn->pppcIntAtomNamePairs[ii][jj], res1_atmname);
              strcat(inp_psn->pppcIntAtomNamePairs[ii][jj], ":");
              strcat(inp_psn->pppcIntAtomNamePairs[ii][jj], res2_atmname);
              strcat(inp_psn->pppcIntAtomNamePairs[ii][jj], ":");
              strcat(inp_psn->pppcIntAtomNamePairs[ii][jj], cDistRep);
              strcat(inp_psn->pppcIntAtomNamePairs[ii][jj], ",");
              
              strcat(inp_psn->pppcIntAtomNamePairs[jj][ii], res1_atmname);
              strcat(inp_psn->pppcIntAtomNamePairs[jj][ii], ":");
              strcat(inp_psn->pppcIntAtomNamePairs[jj][ii], res2_atmname);
              strcat(inp_psn->pppcIntAtomNamePairs[jj][ii], ":");
              strcat(inp_psn->pppcIntAtomNamePairs[jj][ii], cDistRep);
              strcat(inp_psn->pppcIntAtomNamePairs[jj][ii], ",");
            }
          }
        }
      }
      
      else
      {
        // PBC !
        for( kk=0; kk<inp_psn->reslength[ii]; kk++)
        {
          for( ll=0; ll<inp_psn->reslength[jj]; ll++)
          {
            dist = DistanceCoor(trj_crd->xcoor[inp_psn->atmList[ii][kk] - 1], trj_crd->ycoor[inp_psn->atmList[ii][kk] - 1], trj_crd->zcoor[inp_psn->atmList[ii][kk] - 1],
                                trj_crd->xcoor[inp_psn->atmList[jj][ll] - 1], trj_crd->ycoor[inp_psn->atmList[jj][ll] - 1], trj_crd->zcoor[inp_psn->atmList[jj][ll] - 1],
                                trj_crd->pbc);
            
            if( dist<=inp_psn->fDistCutoff )
              interacting++;
          }
        }
      }
      
      // Update ppiIntAtomPairs
      inp_psn->ppiIntAtomPairs[ii][jj] = interacting;
      inp_psn->ppiIntAtomPairs[jj][ii] = interacting;
      //strcat(inp_psn->pppcIntAtomNamePairs[ii][jj], "?:?,");


      // Update ppfIntStrength
      fIntStrength = (interacting/sqrt(inp_psn->res_norm[ii]*inp_psn->res_norm[jj]))*100;
      inp_psn->ppfIntStrength[ii][jj] = fIntStrength;
      inp_psn->ppfIntStrength[jj][ii] = fIntStrength;
      
      // Update ppiResResIntFreq
      if(fIntStrength != 0.0)
      {
        inp_psn->ppiResResIntFreq[ii][jj]++;
        inp_psn->ppiResResIntFreq[jj][ii]++;
      }
      
      // Update IntStrength for Hubs identification
      if(inp_psn->iHubEqFlag == 1)
      {
        // use modified equation for hubs
        inp_psn->ppfHubsIntStrength[ii][jj] = (interacting/inp_psn->res_norm[ii])*100;
        inp_psn->ppfHubsIntStrength[jj][ii] = (interacting/inp_psn->res_norm[jj])*100;
      }
      
      else
      {
        // use normal equation
        inp_psn->ppfHubsIntStrength[ii][jj] = fIntStrength;
        inp_psn->ppfHubsIntStrength[jj][ii] = fIntStrength;
      }
   }
  }
  
  // === NoLink Section ================================================
  if(inp_psn->iNumOfNoLinkPairs != 0)
  {
    // Set to 0.0 Interaction Strength
    for(ii=0; ii<inp_psn->iNumOfNoLinkPairs; ii++)
    {
      iTempRes1 = inp_psn->ppiNoLinkPairs[ii][0];
      iTempRes2 = inp_psn->ppiNoLinkPairs[ii][1];
      
      inp_psn->ppfIntStrength[iTempRes1][iTempRes2] = 0.0;
      inp_psn->ppfIntStrength[iTempRes2][iTempRes1] = 0.0;
      
      inp_psn->ppfHubsIntStrength[iTempRes1][iTempRes2] = 0.0;
      inp_psn->ppfHubsIntStrength[iTempRes2][iTempRes1] = 0.0;
    }
  }
  // ===================================================================
  
  if(inp_psn->iVerboseFlag==1)
  {
    fprintf(inp_psn->outfile, ">INT\n");
    fprintf(inp_psn->outfile, "  Res1 Res2    RIS    HIS    ATM   Pairs\n");
  }
  for(ii=0; ii<inp_psn->tot_nresidues; ii++)
  {
    for(jj=0; jj<inp_psn->tot_nresidues; jj++)
    {
      /* === post-proc === */
      inp_psn->ppfAvgResInt[ii][jj]=inp_psn->ppfAvgResInt[ii][jj]+inp_psn->ppfIntStrength[ii][jj];
      /* === post-proc === */
      
      if(inp_psn->iVerboseFlag==1)
      {
        if(inp_psn->ppfIntStrength[ii][jj]!=0.0 || inp_psn->ppfHubsIntStrength[ii][jj]!=0.0)
        {
          
          inp_psn->pppcIntAtomNamePairs[ii][jj][strlen(inp_psn->pppcIntAtomNamePairs[ii][jj]) - 1] = '\0';
          
          fprintf( inp_psn->outfile, "& %4d %4d %6.3f %6.3f    %3d   %s\n", molecule->pRes[ii]->presn, molecule->pRes[jj]->presn, inp_psn->ppfIntStrength[ii][jj], inp_psn->ppfHubsIntStrength[ii][jj], inp_psn->ppiIntAtomPairs[ii][jj], inp_psn->pppcIntAtomNamePairs[ii][jj]);
        }
      }
    }
  }
  
  if(inp_psn->iVerboseFlag==1)
    fprintf(inp_psn->outfile, ">CLS\n");
  
  fIntMinIterator=inp_psn->fIntMinStart;  
  iIntMinIterationNum=0;
  for(iIntMinIterationNum=0; iIntMinIterationNum<inp_psn->iNumOfIntMinStep; iIntMinIterationNum++)
  {

    for(ii=0; ii<inp_psn->tot_nresidues; ii++)
      inp_psn->piNumResInteractions[ii]=0;
    
    for(ii=0; ii<inp_psn->tot_nresidues; ii++)
      for(jj=0; jj<inp_psn->iMaxResInteractions; jj++)
        inp_psn->ppiInteractions[ii][jj]=0;
      
    for(resnumber1=1; resnumber1<=molecule->nRes; resnumber1++)
    {
      iContactsCounter=0;
      for(resnumber2=1; resnumber2<=molecule->nRes; resnumber2++)
      {
        if(inp_psn->ppfIntStrength[resnumber1-1][resnumber2-1]>fIntMinIterator)
        {
          inp_psn->ppiInteractions[resnumber1-1][inp_psn->piNumResInteractions[resnumber1-1]]=resnumber2-1;
          inp_psn->piNumResInteractions[resnumber1-1]++;
          iContactsCounter++;
          
          /* === post-proc === */
          inp_psn->pppfStableResIntStrength[iIntMinIterationNum][resnumber1-1][resnumber2-1] = inp_psn->pppfStableResIntStrength[iIntMinIterationNum][resnumber1-1][resnumber2-1] + inp_psn->ppfIntStrength[resnumber1-1][resnumber2-1];
          inp_psn->pppiStableResInt[iIntMinIterationNum][resnumber1-1][resnumber2-1]++;
          if(inp_psn->ppfHubsIntStrength[resnumber1-1][resnumber2-1]>fIntMinIterator)
          {
            inp_psn->ppiTmpStableHubs[iIntMinIterationNum][resnumber1-1]++;
          }
          /* === post-proc === */
        }
      }
      
      if(inp_psn->iMaxNumOfLink<iContactsCounter)
        inp_psn->iMaxNumOfLink=iContactsCounter;
          
      iContactsCounter=0;
      
    }
    
    linkwalk( molecule->nRes, inp_psn->ppiInteractions, inp_psn->piNumResInteractions, &inp_psn->cluster);
    iClustSize=0;
    // run over the number of found clusters
    if(inp_psn->iVerboseFlag==1)
      fprintf(inp_psn->outfile, "Imin: %5.3f\n", fIntMinIterator);
    for(ww=0; ww<inp_psn->cluster.nclusters; ww++)
    {
      if(inp_psn->iVerboseFlag==1)
        fprintf(inp_psn->outfile, "C %d: ", ww+1);
      // run over the residues in the ww-th cluster
      for(yy=0; yy<inp_psn->cluster.cluspop[ww]; yy++)
      {
        if(inp_psn->cluster.cluspop[ww]>iClustSize)
        {
          iClustSize=inp_psn->cluster.cluspop[ww];
        }
        /* === post-proc === */
        inp_psn->pppiClusters[iIntMinIterationNum][ww][inp_psn->cluster.clusters[ww][yy]]++;
        /* === post-proc === */
        
        if(inp_psn->iVerboseFlag==1)
          fprintf(inp_psn->outfile, "%d ", molecule->pRes[inp_psn->cluster.clusters[ww][yy]]->presn );
      }
      if(inp_psn->iVerboseFlag==1)
        fprintf(inp_psn->outfile, "\n");
    }
    
    /* === post-proc === */
    inp_psn->ppiLargestClusterSize[iIntMinIterationNum][0]=inp_psn->ppiLargestClusterSize[iIntMinIterationNum][0]+iClustSize;
    /* === post-proc === */
    
    // *-* this line must be the last of this for cycle *-* //
    fIntMinIterator += inp_psn->fIntMinStep;
  }

  // === post-proc === 
  for(ii=0; ii<inp_psn->iNumOfIntMinStep; ii++)
  {
    res1=res2=0;
    for(jj=0; jj<molecule->nRes; jj++)
    {
      if(inp_psn->ppiTmpStableHubs[ii][jj]>inp_psn->iHubContCutoff)
      {
        inp_psn->ppiStableHubs[ii][jj]++;
        res1=1;
      }
      else
      {
        res1=0;
      }
    
      for(kk=0; kk<molecule->nRes; kk++)
      {
        if(inp_psn->ppiTmpStableHubs[ii][kk]>inp_psn->iHubContCutoff)
        {
          res2=1;
        }
        else
        {
          res2=0;
        }
        
        if(res1==res2)
        {
          inp_psn->ppppiHubCorr[ii][jj][kk][0]++;
          if(res1==1)
          {
            inp_psn->ppppiHubCorr[ii][jj][kk][1]++;
          }
          else
          {
            inp_psn->ppppiHubCorr[ii][jj][kk][2]++;
          }
        }
      }
    }
  }
  
  for(ii=0; ii<inp_psn->iNumOfIntMinStep; ii++)
  {
    for(jj=0; jj<molecule->nRes; jj++)
    {
      inp_psn->ppiTmpStableHubs[ii][jj]=0;
    }
  }
  
  // === post-proc === 
  inp_psn->iPBCFlag = trj_crd->pbc_flag;
  return 0;
}
// ------------------------------------------------------------------
int Post_PSG(struct inp_psn *inp_psn, int iNumOfFrames, Molecule *molecule, struct sopt *OPT)
{
  
  int                   ii, jj, kk, mm, ww, xx, yy, zz;
  int                   atmCounter1, atmCounter2, resindex;
  int                   iHubFlag, iClustNum, iClustSize, iResIntStabFlag=0;
  char                  cFileName[80], cPDBFileName[80], cResName[6];
  float                 fStable, fFreq, fImin, fStableAvgIntStrength;
  float                 fHubCorr, fHubTimes, fNotHubTimes;
  FILE                  *DistilledOutput;  

  if(inp_psn->iGetIminFlag == 1)
  {
    if(inp_psn->iSecondRound == 0)
    {
      inp_psn->iSecondRound = 1;
      GetImin(inp_psn, OPT, molecule);
      
      // rewind trj and let's start again from the first frame!!
      return 999;
    }
  }

  sprintf(cFileName, "%s%s", "avg", inp_psn->title);
  
  DistilledOutput=O_File(cFileName, "w");
  
  time(&inp_psn->time_tToday);

  fprintf(DistilledOutput, "# ======================================================================\n");
  fprintf(DistilledOutput, "#                       --- Protein Structure Network Analysis ---\n");
  fprintf(DistilledOutput, "#\n");
  fprintf(DistilledOutput, "# Copyright     : Francesca Fanelli\n");
  fprintf(DistilledOutput, "#                 Angelo Felline,\n");
  fprintf(DistilledOutput, "#                 Michele Seeber (2009)\n");
  fprintf(DistilledOutput, "# License       : GPL v 3\n");
  fprintf(DistilledOutput, "#\n");
  fprintf(DistilledOutput, "# PSN Ver       : 0.8c\n");
  fprintf(DistilledOutput, "#\n");
  fprintf(DistilledOutput, "# Date          : %s", asctime(localtime(&inp_psn->time_tToday)));
  fprintf(DistilledOutput, "#\n"); 
  fprintf(DistilledOutput, "# Mol Name      : %s\n", OPT->IMOL_FILE);
  fprintf(DistilledOutput, "# Seg Num       : %d\n", molecule->nSeg);
  fprintf(DistilledOutput, "# Res Num       : %d\n", inp_psn->tot_nresidues);
  fprintf(DistilledOutput, "# Trj Name      : %s\n", OPT->ITRJ_FILE);
  fprintf(DistilledOutput, "# Frame Num     : %d\n", iNumOfFrames);
  
  if(inp_psn->iPBCFlag == 0)
    fprintf(DistilledOutput, "# PBC           : No\n");
  else if(inp_psn->iPBCFlag == -1)
    fprintf(DistilledOutput, "# PBC           : No (forced)\n");
  else
    fprintf(DistilledOutput, "# PBC           : Yes\n");
  
  fprintf(DistilledOutput, "#\n");
  fprintf(DistilledOutput, "# Title         : %s\n", inp_psn->title);
  fprintf(DistilledOutput, "# Selection     : %s\n", inp_psn->sele.selestring);
  fprintf(DistilledOutput, "# Sele Res Num  : %d\n", inp_psn->iNumOfSelRes);
  fprintf(DistilledOutput, "# DistCutOff    : %f\n", inp_psn->fDistCutoff);
  
  fprintf(DistilledOutput, "# IntMin        : %s\n", inp_psn->cIminRange);
  if(inp_psn->iGetIminFlag == 1)
  {
    fprintf(DistilledOutput, "# LCS Fract     : %f\n", inp_psn->fPCNSize);
    fprintf(DistilledOutput, "# Dec Pos       : %d\n", inp_psn->iDecNum);
    fprintf(DistilledOutput, "# Pre  Ic Val   : %f (delta %3.2f)\n", inp_psn->fPreIcValue, fabs(inp_psn->fPreIcValueDelta));
    fprintf(DistilledOutput, "# Post Ic Val   : %f (delta %3.2f)\n", inp_psn->fPostIcValue, fabs(inp_psn->fPostIcValueDelta));
    
    if(inp_psn->iBestIcValue == 0)
      fprintf(DistilledOutput, "# Best Ic Val   : PRE\n");
    else if(inp_psn->iBestIcValue == 1)
      fprintf(DistilledOutput, "# Best Ic Val   : POST\n");
    else if(inp_psn->iBestIcValue == 2)
      fprintf(DistilledOutput, "# Best Ic Val   : BOTH\n");
  }

  else
  {
    fprintf(DistilledOutput, "# IntMinStart   : %f\n", inp_psn->fIntMinStart);
    fprintf(DistilledOutput, "# IntMinStop    : %f\n", inp_psn->fIntMinStop);
    fprintf(DistilledOutput, "# IntMinStep    : %f\n", inp_psn->fIntMinStep);
  }
  
  fprintf(DistilledOutput, "# StableCutOff  : %f\n", inp_psn->fStableCutoff);
  
  if(inp_psn->iHubEqFlag == 0)
    fprintf(DistilledOutput, "# HubEquation   : NO\n");
  else
    fprintf(DistilledOutput, "# HubEquation   : YES\n");
  
  fprintf(DistilledOutput, "# HubContCutOff : %i\n", inp_psn->iHubContCutoff);
  fprintf(DistilledOutput, "# Termini       : %i\n", inp_psn->iTerminiFlag);
  fprintf(DistilledOutput, "# Proximity     : %i\n", inp_psn->iProxCutOff);
  
  if(inp_psn->iIntType == 0)
    fprintf(DistilledOutput, "# Int Type      : SC\n");
  else if(inp_psn->iIntType == 1)
    fprintf(DistilledOutput, "# Int Type      : SC+CA\n");
  else if(inp_psn->iIntType == 2)
    fprintf(DistilledOutput, "# Int Type      : ALL\n");
  
  fprintf(DistilledOutput, "#\n");
  fprintf(DistilledOutput, "# Num Of NoLink : %d\n", inp_psn->iNumOfNoLinkPairs);
  for(ii=0; ii<inp_psn->iNumOfNoLinkPairs; ii++)
    fprintf(DistilledOutput, "# No Link       : #%d %d   %d\n", ii+1, inp_psn->ppiNoLinkPairs[ii][0]+1, inp_psn->ppiNoLinkPairs[ii][1]+1);
  
  fprintf(DistilledOutput, "# Num Of Merge  : %d\n", inp_psn->iNumOfMerge);
  for(ii=0; ii<inp_psn->iNumOfMerge; ii++)
    fprintf(DistilledOutput, "# Merge         : #%d %s\n", ii+1, inp_psn->ppcMergedResidues[ii]);
  
  fprintf(DistilledOutput, "# Num Of Param  : %d\n", inp_psn->iNumOfParam);
  for(ii=0; ii<inp_psn->iNumOfParam; ii++)
    fprintf(DistilledOutput, "# Param         : #%d %s %f\n", ii+1, inp_psn->pcParamResVect[ii], inp_psn->pfParamValVect[ii]);
  
  fprintf(DistilledOutput, "# ======================================================================\n");

  fprintf(DistilledOutput, "\n");

  fprintf( DistilledOutput, "*** Seg Info ***\n");
  jj=1;
  for(ii=0; ii<molecule->nSeg; ii++)
  {
    fprintf( DistilledOutput, "%2d   %5s   %6d   %10d\n", ii+1, molecule->segment[ii].segName, molecule->segment[ii].nRpS, jj);
    jj=jj+molecule->segment[ii].nRpS;
  }
  fprintf( DistilledOutput, "========================\n");
  
  fprintf( DistilledOutput, "*** Seq ***\n");
  
  kk=0;
  for(ii=0; ii<molecule->nSeg; ii++)
  {
    for(jj=0; jj<molecule->segment[ii].nRpS; jj++)
    {
      kk++;

      if     (strncmp(molecule->segment[ii].pRes[jj].resType, "ALA", 3)==0) strcpy(cResName, "A\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "ARG", 3)==0) strcpy(cResName, "R\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "ASN", 3)==0) strcpy(cResName, "N\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "ASP", 3)==0) strcpy(cResName, "D\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "CYS", 3)==0) strcpy(cResName, "C\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "CYN", 3)==0) strcpy(cResName, "C\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "GLU", 3)==0) strcpy(cResName, "E\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "GLN", 3)==0) strcpy(cResName, "Q\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "GLY", 3)==0) strcpy(cResName, "G\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "HIS", 3)==0) strcpy(cResName, "H\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "HSC", 3)==0) strcpy(cResName, "H\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "HSD", 3)==0) strcpy(cResName, "H\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "HSP", 3)==0) strcpy(cResName, "H\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "HID", 3)==0) strcpy(cResName, "H\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "ILE", 3)==0) strcpy(cResName, "I\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "LEU", 3)==0) strcpy(cResName, "L\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "LYS", 3)==0) strcpy(cResName, "K\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "LYP", 3)==0) strcpy(cResName, "K\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "MET", 3)==0) strcpy(cResName, "M\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "PHE", 3)==0) strcpy(cResName, "F\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "PRO", 3)==0) strcpy(cResName, "P\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "SER", 3)==0) strcpy(cResName, "S\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "THR", 3)==0) strcpy(cResName, "T\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "TRP", 3)==0) strcpy(cResName, "W\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "TYR", 3)==0) strcpy(cResName, "Y\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "VAL", 3)==0) strcpy(cResName, "V\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "RET", 3)==0) strcpy(cResName, "X\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "GDP", 3)==0) strcpy(cResName, "d\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "GTP", 3)==0) strcpy(cResName, "t\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "Mg+", 3)==0) strcpy(cResName, "m\0");
      else if(strncmp(molecule->segment[ii].pRes[jj].resType, "ZMA", 3)==0) strcpy(cResName, "Z\0");
      else                                                                  strcpy(cResName, "U\0");
      
      fprintf(DistilledOutput, "%-10d   %5s:%-10d   %5s:%1s%-10d   %.5f\n", kk, molecule->segment[ii].segName, molecule->segment[ii].pRes[jj].resn, molecule->segment[ii].segName, cResName, molecule->segment[ii].pRes[jj].resn, inp_psn->res_norm[kk - 1]);
    }
  }
  fprintf( DistilledOutput, "========================\n");

  fStable=(inp_psn->fStableCutoff*iNumOfFrames);
    
  fprintf(DistilledOutput, "*** Averaged Interaction Strength ***\n");
  fprintf(DistilledOutput, "%-10s   %-10s   %-10s   %-s\n", "Res1", "Res2", "I.S.", "Freq");
  atmCounter1 = 0;
  
  for(ii=0; ii<molecule->nRes; ii++)
  {
    atmCounter2 = 0;
    for( jj=0; jj<=ii; jj++)
      atmCounter2 += inp_psn->reslength2[jj];
      
    for(jj=ii+1; jj<molecule->nRes; jj++)
    {
      if(inp_psn->ppfAvgResInt[ii][jj]>0.0)
      {
        //fprintf(DistilledOutput, " %4s:%-4d   %4s:%-4d   %5.3f   %.3f\n", molecule->rawmol.segId[atmCounter1], molecule->rawmol.resn[atmCounter1], molecule->rawmol.segId[atmCounter2], molecule->rawmol.resn[atmCounter2], (inp_psn->ppfAvgResInt[ii][jj]/(float)iNumOfFrames), ((inp_psn->ppiResResIntFreq[ii][jj]/(float)iNumOfFrames)*100.0));
        fprintf(DistilledOutput, "%4s:%-4d   %4s:%-4d   %7.3f   %.3f\n", molecule->rawmol.segId[atmCounter1], molecule->rawmol.resn[atmCounter1], molecule->rawmol.segId[atmCounter2], molecule->rawmol.resn[atmCounter2], (inp_psn->ppfAvgResInt[ii][jj]/(float)iNumOfFrames), ((inp_psn->ppiResResIntFreq[ii][jj]/(float)iNumOfFrames)*100.0));
      }
      atmCounter2 += inp_psn->reslength2[jj];
    }
    atmCounter1 += inp_psn->reslength2[ii];
  }

  fprintf(DistilledOutput, "========================\n");

  fprintf(DistilledOutput, "*** Stable Residue Interactions ***\n");
  //fprintf(DistilledOutput, " Imin     Res1        Res2 Stable   Freq\n");
  fprintf(DistilledOutput, "%5s   %-9s   %-9s   %6s   %6s   %s\n", "Imin", "  Res1", "  Res2", "Stable", "Freq", "AvgInt");
  
  for(ii=0; ii<inp_psn->iNumOfIntMinStep; ii++)
  {
    fImin=(inp_psn->fIntMinStart+(inp_psn->fIntMinStep*ii));
    atmCounter1 = 0;
    for(jj=0; jj<molecule->nRes; jj++)
    {
      atmCounter2 = 0;
      for( kk=0; kk<=jj; kk++)
        atmCounter2 += inp_psn->reslength2[kk];
        
      for(kk=jj+1; kk<molecule->nRes; kk++)
      {
        if(inp_psn->pppiStableResInt[ii][jj][kk] >= fStable)
          iResIntStabFlag = 1;
        else
          iResIntStabFlag = 0;
        
        fFreq=((inp_psn->pppiStableResInt[ii][jj][kk]*100)/(float)iNumOfFrames);
        
        if(inp_psn->pppiStableResInt[ii][jj][kk] != 0.0)
          fStableAvgIntStrength = (inp_psn->pppfStableResIntStrength[ii][jj][kk] / (float) inp_psn->pppiStableResInt[ii][jj][kk]);
        else
          fStableAvgIntStrength = 0.0;
        
        fprintf(DistilledOutput, "%5.3f   %4s:%-4d   %4s:%-4d   %6d   %6.2f   %7.3f\n",
                                 fImin,
                                 molecule->rawmol.segId[atmCounter1], molecule->rawmol.resn[atmCounter1],
                                 molecule->rawmol.segId[atmCounter2], molecule->rawmol.resn[atmCounter2],
                                 iResIntStabFlag, fFreq, fStableAvgIntStrength);

        atmCounter2 += inp_psn->reslength2[kk];
      }
      atmCounter1 += inp_psn->reslength2[jj];
    }
  }
  fprintf(DistilledOutput, "========================\n");  

  fprintf(DistilledOutput, "*** Hubs Frequencies ***\n");
  fprintf(DistilledOutput, " Imin     Res     Stable      Frq\n");  
  for(ii=0; ii<inp_psn->iNumOfIntMinStep; ii++)
  {
    // === clears original BetaFactors  ==========================
    for( ww=0; ww< molecule->nSeg; ww++)
      for( xx=0; xx< molecule->segment[ww].nRpS; xx++)
       for( yy=0; yy< molecule->segment[ww].pRes[xx].nApR; yy++)
       {
         molecule->segment[ww].pRes[xx].pAto[yy].occup = 0.0;
         molecule->segment[ww].pRes[xx].pAto[yy].bFac  = 0.0;
       }
       
    for( ww=0; ww< molecule->nato; ww++)
    {
      molecule->rawmol.occup[ww] = 0.0;
      molecule->rawmol.bFac[ww]  = 0.0;
    }
    // ===========================================================

    fImin=(inp_psn->fIntMinStart+(inp_psn->fIntMinStep*ii));
    atmCounter1 = 0;
    for(jj=0; jj<molecule->nRes; jj++)
    {
      fFreq=((inp_psn->ppiStableHubs[ii][jj]*100)/(float)iNumOfFrames);
      if(inp_psn->ppiStableHubs[ii][jj]>=fStable)
      {
        iHubFlag=1;
      }
      else
      {
        iHubFlag=0;
      }
      fprintf(DistilledOutput, "%5.3f   %4s:%-4d   %4d   %6.2f\n", fImin, molecule->rawmol.segId[atmCounter1], molecule->rawmol.resn[atmCounter1], iHubFlag, fFreq);
      
      for( mm=0; mm<inp_psn->reslength2[jj]; mm++)
      {
        molecule->rawmol.occup[atmCounter1+mm] = fFreq/100.0;
        molecule->rawmol.bFac[atmCounter1+mm]  = iHubFlag;
      }
      atmCounter1 += inp_psn->reslength2[jj];
    }
    // write down a pdb with cluster-index on b-factor field
    sprintf( cPDBFileName, "hub%s-%4.3f.pdb", inp_psn->title, (inp_psn->fIntMinStart+(inp_psn->fIntMinStep*ii)));
    WritePdb_unstr( cPDBFileName, molecule );
  }
  fprintf(DistilledOutput, "========================\n");  
  
  fprintf(DistilledOutput, "*** Hub Correlations ***\n");
  fprintf(DistilledOutput, " Imin    Res1        Res2         Corr     HubCorr  NotHubCorr \n");
  for(ii=0; ii<inp_psn->iNumOfIntMinStep; ii++)
  {
    fImin=(inp_psn->fIntMinStart+(inp_psn->fIntMinStep*ii));
    atmCounter1 = 0;
    for(jj=0; jj<molecule->nRes; jj++)
    {
      atmCounter2 = 0;
      for(kk=0; kk<molecule->nRes; kk++)
      {
        fHubCorr=((float)inp_psn->ppppiHubCorr[ii][jj][kk][0]/(float)iNumOfFrames);
        fHubTimes=((float)inp_psn->ppppiHubCorr[ii][jj][kk][1]/(float)iNumOfFrames);
        fNotHubTimes=((float)inp_psn->ppppiHubCorr[ii][jj][kk][2]/(float)iNumOfFrames);
        if(inp_psn->ppppiHubCorr[ii][jj][kk][2]>inp_psn->ppppiHubCorr[ii][jj][kk][1])
        {
          fHubCorr=fHubCorr*-1.0;
        }
        fprintf(DistilledOutput, "%5.3f   %4s:%-4d   %4s:%-4d   %6.2f      %6.2f      %6.2f\n", fImin, molecule->rawmol.segId[atmCounter1], molecule->rawmol.resn[atmCounter1], molecule->rawmol.segId[atmCounter2], molecule->rawmol.resn[atmCounter2], fHubCorr, fHubTimes, fNotHubTimes);
        atmCounter2 += inp_psn->reslength2[kk];
      }
      atmCounter1 += inp_psn->reslength2[jj];
    }
  }
  fprintf(DistilledOutput, "========================\n");  

  fprintf(DistilledOutput, "*** Stable Cluster Compositions ***\n");
  for(ii=0; ii<inp_psn->iNumOfIntMinStep; ii++)
  {
    // === clears original BetaFactors  ==========================
    for( ww=0; ww< molecule->nSeg; ww++)
      for( xx=0; xx< molecule->segment[ww].nRpS; xx++)
       for( yy=0; yy< molecule->segment[ww].pRes[xx].nApR; yy++)
         molecule->segment[ww].pRes[xx].pAto[yy].bFac = 0.0;
    for( ww=0; ww< molecule->nato; ww++)
      molecule->rawmol.bFac[ww] = 0.0;
    // ===========================================================
  
    for(jj=0; jj<inp_psn->tot_nresidues; jj++)
      inp_psn->piNumResInteractions[jj]=0;
     
    for(jj=0; jj<inp_psn->tot_nresidues; jj++)
      for(kk=0; kk<inp_psn->iMaxResInteractions; kk++)
        inp_psn->ppiInteractions[jj][kk]=0;
    
    fImin=(inp_psn->fIntMinStart+(inp_psn->fIntMinStep*ii));
    for(jj=0; jj<inp_psn->tot_nresidues; jj++)
    {
      for(kk=0; kk<inp_psn->tot_nresidues; kk++)
      {
        if(inp_psn->pppiStableResInt[ii][jj][kk]>fStable)
        {
          inp_psn->ppiInteractions[jj][inp_psn->piNumResInteractions[jj]]=kk;
          inp_psn->piNumResInteractions[jj]++;
        }
      }
    }
    linkwalk( molecule->nRes, inp_psn->ppiInteractions, inp_psn->piNumResInteractions, &inp_psn->cluster);
    fprintf(DistilledOutput, "Imin: %5.3f\n", fImin);

    iClustNum=0;
    iClustSize=0;
    for(ww=0; ww<inp_psn->cluster.nclusters; ww++)
    {
      if(inp_psn->cluster.cluspop[ww]>2)
      {
        if(iClustSize<inp_psn->cluster.cluspop[ww])
        {
          iClustSize=inp_psn->cluster.cluspop[ww];
        }
        iClustNum++;
        fprintf(DistilledOutput, "C %2d: ", iClustNum);
        for(yy=0; yy<inp_psn->cluster.cluspop[ww]; yy++)
        {
          atmCounter2=0;
          for(zz=0; zz<inp_psn->cluster.clusters[ww][yy]; zz++)
              atmCounter2 += inp_psn->reslength2[zz];
          
          //fprintf(DistilledOutput, "%4s:%-d ", molecule->rawmol.segId[atmCounter2], inp_psn->cluster.clusters[ww][yy]+1); *-*
          fprintf(DistilledOutput, "%4s:%-4d ", molecule->rawmol.segId[atmCounter2], molecule->rawmol.resn[atmCounter2]);
          atmCounter1 = 0;
          resindex = inp_psn->cluster.clusters[ww][yy];
          for( mm=0; mm < resindex; mm++ )
            atmCounter1 += inp_psn->reslength2[mm];
          for( mm=0; mm<inp_psn->reslength2[resindex]; mm++ )
            molecule->rawmol.bFac[atmCounter1+mm] = iClustNum;
        }
        fprintf(DistilledOutput, "\n");
      }
    }

    // write down a pdb with cluster-index on b-factor field
    sprintf( cPDBFileName, "cls%s-%4.3f.pdb", inp_psn->title, fImin);
    WritePdb_unstr( cPDBFileName, molecule );
    inp_psn->ppiLargestClusterSize[ii][1]=iClustSize;
  }
  fprintf(DistilledOutput, "========================\n");

  fprintf(DistilledOutput, "*** Largest Cluster Size ***\n");
  fprintf(DistilledOutput, "  Imin       Sim%%        Sim       Avg%%       Avg\n");
  for(ii=0; ii<inp_psn->iNumOfIntMinStep; ii++)
  {
    fImin=(inp_psn->fIntMinStart+(inp_psn->fIntMinStep*ii));
    fprintf(DistilledOutput, "%6.3f     %6.2f     %6.2f     %6.2f     %5d\n",
                                                            fImin,
                                                            (((inp_psn->ppiLargestClusterSize[ii][0]/(float)iNumOfFrames)*100)/(float)inp_psn->tot_nresidues),
                                                            (inp_psn->ppiLargestClusterSize[ii][0]/(float)iNumOfFrames),
                                                            ((inp_psn->ppiLargestClusterSize[ii][1]*100)/(float)inp_psn->tot_nresidues),
                                                            inp_psn->ppiLargestClusterSize[ii][1]);
  }  
  fprintf(DistilledOutput, "========================\n");
  
  // === Node Degrees =============================================== //
  for(ii=0; ii<inp_psn->iNumOfIntMinStep; ii++)
  {
    fImin=(inp_psn->fIntMinStart+(inp_psn->fIntMinStep*ii));

    for(ww=0; ww<molecule->nRes; ww++)
      inp_psn->piNodeDegree[ww] = 0;

    // === clears original BetaFactors  ==========================
    for( ww=0; ww< molecule->nSeg; ww++)
      for( xx=0; xx< molecule->segment[ww].nRpS; xx++)
       for( yy=0; yy< molecule->segment[ww].pRes[xx].nApR; yy++)
         molecule->segment[ww].pRes[xx].pAto[yy].bFac  = 0.0;
       
    for( ww=0; ww< molecule->nato; ww++)
      molecule->rawmol.bFac[ww]  = 0.0;
    // ===========================================================
    
    for(jj=0; jj<molecule->nRes-1; jj++)
    {
      for(kk=jj+1; kk<molecule->nRes; kk++)
      {
        if(inp_psn->pppiStableResInt[ii][jj][kk] >= fStable)
        {
          inp_psn->piNodeDegree[jj]++;
          inp_psn->piNodeDegree[kk]++;
        }
      }
    }
    
    atmCounter1 = 0;
    for(jj=0; jj<molecule->nRes; jj++)
    {
      for( mm=0; mm<inp_psn->reslength2[jj]; mm++)
        molecule->rawmol.bFac[atmCounter1+mm] = (float)inp_psn->piNodeDegree[jj];

      atmCounter1 += inp_psn->reslength2[jj];
    }
    
    // write down a pdb with degree values on b-factor field
    sprintf(cPDBFileName, "deg%s-%4.3f.pdb", inp_psn->title, fImin);
    WritePdb_unstr( cPDBFileName, molecule );
  }
  // === Node Degrees =============================================== //
  
  return 0;  
}
// ------------------------------------------------------------------

/*
 * ========================== *
 * === PATHFINDER SECTION === *
 * ========================== *
 */

// ------------------------------------------------------------------
// === Macro Definitions ===============================================
#define PSNPATHVERSTRING     "0.7a"                                     // Version string
// =====================================================================

// === Global Variables ================================================
char    cRawPSNFileName[512];                                           // RAWPSN file name
char    cAvgPSNFileName[512];                                           // AVGPSN file name
char    cCorrFileName[512];                                             // Corr file name
char    cOutPutFileName[512];                                           // Output file name
char    cLogFileName[512];                                              // Log file name
char  **pcSequence;                                                     // Used to store sequence in the form segname:resname+resnum
char  **pcRawSequence;                                                  // Used in SequenceCheck
char  **pcPathList;                                                     // List of paths
char  **pcSegNameList;                                                  // List of segment names
char  **pcTmpPaths;                                                     // Stores tmp paths
char    cStatFileName[512];                                             // Stat file name
char    cFrameFileName[512];                                            // Frame file name
char  **ppcRes1;                                                        // First residue | pairs file name
char  **ppcRes2;                                                        // Last residue | pairs file name
char    cPathTitle[512];                                                // Used only as log file name
char   *pcPath;                                                         // Temp path
char  **ppcNodeLabel;                                                   // List of node labels, used with piNodeIndex

int     MAXNUMOFPATHS=5000;                                             // Max num of paths
int     GAIN=10;                                                        // Multiplier
int    *piRes1List;                                                     // The first residue number list
int    *piRes2List;                                                     // The last residue number list
int     iRes1, iRes2;                                                   // The first and the last residue in the path
int     iNumOfRes1;                                                     // Nummber of Res1 residues
int     iNumOfRes2;                                                     // Nummber of Res2 residues
int     iNumOfRes;                                                      // Number of Residues
int     iNumOfSeg;                                                      // Number of Segments
int     iNumOfFrames;                                                   // Number of Frames
int     iNumOfEdges;                                                    // Number of Edges
int     iNumOfPaths=0;                                                  // Number of Paths
int   **ppiPrevNode;                                                    // The node that comes right before "ii"
int    *piVisited;                                                      // Visited nodes
int    *piPathFreq;                                                     // Path Frequencies
int    *piNodesInPath;                                                  // Number of nodes in selected path
int    *piSegOffSet;                                                    // Used to manage segment offset
int     iLastDepth=0;                                                   // Used by GenPath
int   **ppiPathMatrix;                                                  // Used by GenPath
int    *piTmpPathsLen;                                                  // Tmp Path Length
int     iNumOfGenPaths=0;                                               // Number Of Path(s) generated by GenPath
int     iFrameCont;                                                     // Frame number
int     iStatFlag=0;                                                    // Used to manage --STAT option
int     iModeFlag=0;                                                    // Used to manage --MODE option
int    *piCorrVector;                                                   // Used to store correlated residues in --MODE cp
int     iNumOfCorrRes=0;                                                // Number of Correlated Residues, used in --MODE cp
int     iFrameFlag=0;                                                   // Used by --FRAME option
int    *piBadFrames;                                                    // Used to select productive frames
int     iNumOfBadFrames=0;                                              // Total number of bad frames
int     iLogFlag=0;                                                     // If 1 a log file will be created
int     iPairCont=0;                                                    // Res Pairs counter
int     iMatchMode=0;                                                   // Res Pairs Match Mode : 0 = Cross; 1 = Linear
int     iRes1FileFlag=0;                                                // 1 if the first residue is a file
int     iRes2FileFlag=0;                                                // 1 if the last  residue is a file
int     iMinLen=0;                                                      // Path minimum length
int     iFrameBlockNum;                                                 // Actul frame-block number
int     iFramePerBlock=100;                                             // Num of frames per block file
int     iAbsCorrFlag=0;                                                 // Absolute Correlation Values Flag
int     iWeightFlag=0;                                                  // Use weighted links
int    *piClusters;                                                     // The cluster of each node
int     iCheckClustFlag=1;                                              // If 0 does not check for bad frames
int     iStartFrame=1;                                                  // First frame to process
int     iStartPair=1;                                                   // First Pair to process
int     iCalcStep=0;                                                    // Used for debugging purpose
int     iPSNTypeFlag=0;                                                 // Type of PSN file: 0 = raw, 1 = avg
int    *piNodeIndex;                                                    // List of node indexes, used with ppcNodeLabel
int     iNumOfPoxInt=0;                                                 // Num of possible interactions among nodes

float  *pfTmpScore;                                                     // Tmp Path score
float  *pfTmpWeights;                                                   // Tmp Path weights
float   fIntMinCutOff=-1.0;                                             // Interaction Strength to use
float   fCorrCutOff=-2.0;                                               // Lowest Correlation cutoff allowed
float  *pfPathLength;                                                   // The length of the shortest path between iRes1 and iRes2
float  *pfPathWeights;                                                  // The weights of all paths
float **ppfIntMatrix;                                                   // PSN Interaction Matrix
float **ppfRealIntMatrix;                                               // PSN Interaction Matrix with real Interaction Strengths
float **ppfCorMatrix;                                                   // Cirr Matrix
float  *pfPathScores;                                                   // Scores of all path
float   fFreqCutOff=100.0;                                              // Good frames cutoff
float **ppfTmpCorrStat;                                                 // Tmp Corr Stat
float **ppfCorrStatList;                                                // Min, Max and Avg Correlation Statistics of all paths
float   fMinFreq=0.0;                                                   // Path minimum frequency
float   fIntMinFreq=0.0;                                                // Minimum interaction frequency
float **ppfAvgIntStrengths;                                             // Avg Res-Res Interaction Strengths
float **ppfAvgIntFrequencies;                                           // Avg Res-Res Interaction Frequencies
float **ppfStableLink;                                                  // Stable links: 1 stable, 0 unstable
float  *pfIntFreqValList;                                               // List of all interaction frequency values
float **ppfIntMatrixBackup;                                             // A copy of ppfIntMatrix

FILE   *FRawFile;                                                       // rawPSN file handler
FILE   *FAvgFile;                                                       // avgPSN file handler
FILE   *FCorFile;                                                       // Corr file handler
FILE   *FOutFile;                                                       // Output file handler
FILE   *FStatFile;                                                      // Stat file handler
FILE   *FFrameFile;                                                     // Frame file handler
FILE   *FLogFile;                                                       // File with error logs

time_t  time_Start;                                                     // Used to compute calculation time
time_t  time_Finish;                                                    // Used to compute calculation time
// =====================================================================

void PSNPathError(int iErrNum)
{
  
  /*
    PSNPATH         <--->  0
    GetPSNParam     <--->  1
    SequenceCheck   <--->  2
    GetCorrData     <--->  3
    SetIntVar       <--->  4
    ClearAll        <--->  5
    GetBadFrames    <--->  6
    GetCorrRes      <--->  7
    GetPaths        <--->  8
    FindPath        <--->  9
    GenPath         <---> 10
    SavePath        <---> 11
    SaveData        <---> 12
    ClearIntMatrix  <---> 13
    ClearAllLite    <---> 14
    GetFrameInt     <---> 15
    SavePathLite    <---> 16
    GetPSGData      <---> 17
    LoadAvgLabIdx   <---> 18
    LoadStableLinks <---> 19
    IntMatrixFilter <---> 20
  */
  
  char  cErrorFileName[512];
  FILE *FErrorFile;
  
  time(&time_Finish);
  
  fprintf(stderr, "PSNPATH Error - Something terrible happened to system memory ...\n\n");
  fprintf(stderr, "   Date         : %s\n", asctime(localtime(&time_Finish)));
  fprintf(stderr, "   Title        : %s\n", cPathTitle);
  fprintf(stderr, "   Frame        : %d\n", iFrameCont);
  fprintf(stderr, "   Pair         : %d\n", iPairCont);
  fprintf(stderr, "   First Res    : %s\n", pcSequence[iRes1]);
  fprintf(stderr, "   Last  Res    : %s\n", pcSequence[iRes2]);
  if(iModeFlag == 0)
    fprintf(stderr, "   Out File     : %s\n", cOutPutFileName);
  
  else if(iModeFlag == 1)
    fprintf(stderr, "   Out File     : %s\n", cFrameFileName);

  fprintf(stderr, "   Step         : %d\n", iCalcStep);
  fprintf(stderr, "   Error        : %d\n", iErrNum);
  fprintf(stderr, "   Num Of Res   : %d\n", iNumOfRes);
  fprintf(stderr, "   Num Of GenP  : %d\n", iNumOfGenPaths);
  fprintf(stderr, "   Num Of Paths : %d\n", iNumOfPaths);
  fprintf(stderr, "   MEMPATH      : %d\n", MAXNUMOFPATHS);
  
  if(iModeFlag == 0)
    fprintf(stderr, "   MEMGAIN      : %d\n", GAIN);
  
  fprintf(stderr, "\n");
  fprintf(stderr, "To restart the calculation from this point add/modify these lines in your wordom PSNPATH input file and rerun the job.\n\n");
  
  if(iModeFlag == 0)
  {
    // --MODE FULL
    fprintf(stderr, "--STARTPAIR %d\n", iPairCont);
    fprintf(stderr, "--MEMPATH   %d\n", (MAXNUMOFPATHS + (int)(MAXNUMOFPATHS / 3.0)));
    fprintf(stderr, "--MEMGAIN   %d\n", (GAIN + (int)(GAIN / 3.0)));
  }
  
  else if(iModeFlag == 1)
  {
    // --MODE LITE
    fprintf(stderr, "--STARTFRAME %d\n", iFrameCont);
    fprintf(stderr, "--MEMPATH    %d\n", (MAXNUMOFPATHS + (int)(MAXNUMOFPATHS / 3.0)));
    fprintf(stderr, "\n");
    fprintf(stderr, "Open ``%s'' file and delete all lines below the string ``@f%d'' (included).\n", cFrameFileName, iFrameCont);
  }
  
  fprintf(stderr, "\n");
  fprintf(stderr, "Now Exit\n\n");
  
  // ================================================================ //
  
  sprintf(cErrorFileName, "PSNPath-Error-%s.txt", cPathTitle);
  FErrorFile = fopen(cErrorFileName, "a+");
  
  if(FErrorFile == NULL)
  {
    fprintf( stderr, "PSNPATH module: Unable to create file: %s\n", cErrorFileName);
    exit(1);
  }
  
  fprintf(FErrorFile, "PSNPATH Error - Something terrible happened to system memory ...\n\n");
  fprintf(FErrorFile, "   Date         : %s\n", asctime(localtime(&time_Finish)));
  fprintf(FErrorFile, "   Title        : %s\n", cPathTitle);
  fprintf(FErrorFile, "   Frame        : %d\n", iFrameCont);
  fprintf(FErrorFile, "   Pair         : %d\n", iPairCont);
  fprintf(FErrorFile, "   First Res    : %s\n", pcSequence[iRes1]);
  fprintf(FErrorFile, "   Last  Res    : %s\n", pcSequence[iRes2]);
  if(iModeFlag == 0)
    fprintf(FErrorFile, "   Out File     : %s\n", cOutPutFileName);
  
  else if(iModeFlag == 1)
    fprintf(FErrorFile, "   Out File     : %s\n", cFrameFileName);

  fprintf(FErrorFile, "   Step         : %d\n", iCalcStep);
  fprintf(FErrorFile, "   Error        : %d\n", iErrNum);
  fprintf(FErrorFile, "   Num Of Res   : %d\n", iNumOfRes);
  fprintf(FErrorFile, "   Num Of GenP  : %d\n", iNumOfGenPaths);
  fprintf(FErrorFile, "   Num Of Paths : %d\n", iNumOfPaths);
  fprintf(FErrorFile, "   MEMPATH      : %d\n", MAXNUMOFPATHS);
  
  if(iModeFlag == 0)
    fprintf(FErrorFile, "   MEMGAIN      : %d\n", GAIN);
  
  fprintf(FErrorFile, "\n");
  fprintf(FErrorFile, "To restart the calculation from this point add/modify these lines in your wordom PSNPATH input file and rerun the job.\n\n");
  
  if(iModeFlag == 0)
  {
    // --MODE FULL
    fprintf(FErrorFile, "--STARTPAIR %d\n", iPairCont);
    fprintf(FErrorFile, "--MEMPATH   %d\n", (MAXNUMOFPATHS + (int)(MAXNUMOFPATHS / 3.0)));
    fprintf(FErrorFile, "--MEMGAIN   %d\n", (GAIN + (int)(GAIN / 3.0)));
  }
  
  else if(iModeFlag == 1)
  {
    // --MODE LITE
    fprintf(FErrorFile, "--STARTFRAME %d\n", iFrameCont);
    fprintf(FErrorFile, "--MEMPATH    %d\n", (MAXNUMOFPATHS + (int)(MAXNUMOFPATHS / 3.0)));
    fprintf(FErrorFile, "\n");
    fprintf(FErrorFile, "Open ``%s'' file and delete all lines below the string ``@f%d'' (included).\n", cFrameFileName, iFrameCont);
  }
  
  fprintf(FErrorFile, "\n");
  fprintf(FErrorFile, "Now Exit\n\n");
  
  fclose(FErrorFile);
  exit(1);
}

int PSNPATH(char **ppcInput, int iInputLineNum)
{
  int     ii, jj;                                                       // Some iterators
  int     iOriginalInputLineNum;
  int     iOptionFlag;                                                  // Used to catch invalid options
  int     iNumOfResInFile;                                              // Number of residues in resfile
  int     iTmpResNum;                                                   // Used to generate output file names
  int     iBlockCont;                                                   // Save frame in current block
  int     iCalcFlag;                                                    // Used with iCheckClustFlag
  
  int     iAVGFILEOptFlag = 0;
  int     iMININTFREQOptFlag = 0;
  
  char    cTmpSegName[5], cTmpResCode[2];                               // Used to generate output file names
  char    pcInputLine[10240];                                            // Input lines
  char   *pcTest1, *pcTest2;                                            // Used to check if a residue is a file name
  char    cTmpRes1[10240], cTmpRes2[10240];                             // Used to get --PAIR options
  char    cMatchMode[10];                                               // Used to get --MATCH mode
  char    cLine[10240];                                                   // Residue file lines
  
  float   fBadFramesVal=0.0;                                            // Temp Bad frames freq value
  
  FILE   *FResFileHandler;                                              // Used to handle residue files

  strcpy(cPathTitle, "-");
  
  memset ( pcInputLine, '\0', sizeof(pcInputLine));
  
  iOriginalInputLineNum = iInputLineNum;
  
  signal(SIGSEGV, PSNPathError);
  iCalcStep = 0;
  
  // Process input file directives
  while(strncmp(pcInputLine, "END", 3) != 0)
  {
    iOptionFlag = 0;
    sprintf(pcInputLine, "%s", ppcInput[iInputLineNum]);
    
    if(!strncmp(pcInputLine, "BEGIN", 5) || !strncmp(pcInputLine, "END", 3) || pcInputLine[0] == '#')
      iOptionFlag = 1;
    
    else if(strncmp(pcInputLine, "--TITLE", 7) == 0)
    {
      sscanf(pcInputLine, "--TITLE %s", cPathTitle);
      iOptionFlag = 1;
    }
    
    else if (strncmp(pcInputLine, "--PSN", 5) == 0)
    {
      sscanf(pcInputLine, "--PSN %s", cRawPSNFileName);
      
      FRawFile = fopen(cRawPSNFileName, "r");
      if(FRawFile == NULL)
      {
        fprintf( stderr, "PSNPATH module: Unable to open file: %s\n", cRawPSNFileName);
        return 1;
      }
      
      iOptionFlag = 1;
    }
    
    else if (strncmp(pcInputLine, "--AVGFILE", 9) == 0)
    {
      sscanf(pcInputLine, "--AVGFILE %s", cAvgPSNFileName);
      
      FAvgFile = fopen(cAvgPSNFileName, "r");
      if(FAvgFile == NULL)
      {
        fprintf( stderr, "PSNPATH module: Unable to open file: %s\n", cAvgPSNFileName);
        return 1;
      }
      
      iAVGFILEOptFlag = 1;
      iOptionFlag = 1;
    }
    
    else if (strncmp(pcInputLine, "--CORR", 6) == 0)
    {
      sscanf(pcInputLine, "--CORR %s", cCorrFileName);
      
      FCorFile = fopen(cCorrFileName, "r");
      if(FCorFile == NULL)
      {
        fprintf( stderr, "PSNPATH module: Unable to open file: %s\n", cCorrFileName);
        return 1;
      }
      
      iOptionFlag = 1;
    }

    else if (strncmp(pcInputLine, "--IMIN", 6) == 0)
    {
      sscanf(pcInputLine, "--IMIN %f", &fIntMinCutOff);
      iOptionFlag = 1;
    }

    else if (strncmp(pcInputLine, "--CUTOFF", 8) == 0)
    {
      sscanf(pcInputLine, "--CUTOFF %f", &fCorrCutOff);
      iOptionFlag = 1;
    }

    else if (strncmp(pcInputLine, "--MININTFREQ", 12) == 0)
    {
      sscanf(pcInputLine, "--MININTFREQ %f", &fIntMinFreq);
      iOptionFlag = 1;
      iMININTFREQOptFlag = 1;
    }

    else if (strncmp(pcInputLine, "--PAIR", 6) == 0)
    {
      sscanf(pcInputLine, "--PAIR %s %s", cTmpRes1, cTmpRes2);
      iOptionFlag = 1;
    }

    else if (strncmp(pcInputLine, "--MAXBAD", 8) == 0)
    {
      sscanf(pcInputLine, "--MAXBAD %f", &fFreqCutOff);
      iOptionFlag = 1;
    }
    
    else if (strncmp(pcInputLine, "--FBLOCK", 8) == 0)
    {
      sscanf(pcInputLine, "--FBLOCK %d", &iFramePerBlock);
      iOptionFlag = 1;
    }

    else if (strncmp(pcInputLine, "--ABSCORR", 9) == 0)
    {
      sscanf(pcInputLine, "--ABSCORR %d", &iAbsCorrFlag);
      if(iAbsCorrFlag != 0 && iAbsCorrFlag != 1)
      {
        fprintf( stderr, "PSNPATH module: Invalid --ABSCORR value: %i\n", iAbsCorrFlag);
        return 1;
      }

      iOptionFlag = 1;
    }

    else if (strncmp(pcInputLine, "--WEIGHT", 8) == 0)
    {
      sscanf(pcInputLine, "--WEIGHT %d", &iWeightFlag);
      if(iWeightFlag != 0 && iWeightFlag != 1 && iWeightFlag != 2)
      {
        fprintf( stderr, "PSNPATH module: Invalid --WEIGHT value: %i\n", iWeightFlag);
        return 1;
      }

      iOptionFlag = 1;
    }

    else if (strncmp(pcInputLine, "--CHECKCLUST", 12) == 0)
    {
      sscanf(pcInputLine, "--CHECKCLUST %d", &iCheckClustFlag);
      if(iCheckClustFlag != 0 && iCheckClustFlag != 1)
      {
        fprintf( stderr, "PSNPATH module: Invalid --CHECKCLUST value: %i\n", iCheckClustFlag);
        return 1;
      }

      iOptionFlag = 1;
    }

    else if (strncmp(pcInputLine, "--FRAME", 7) == 0)
    {
      sscanf(pcInputLine, "--FRAME %d", &iFrameFlag);
      if(iFrameFlag != 0 && iFrameFlag != 1)
      {
        fprintf( stderr, "PSNPATH module: Invalid --FRAME value: %i\n", iFrameFlag);
        return 1;
      }

      iOptionFlag = 1;
    }

    else if (strncmp(pcInputLine, "--STAT", 6) == 0)
    {
      sscanf(pcInputLine, "--STAT %d", &iStatFlag);
      if(iStatFlag != 0 && iStatFlag != 1)
      {
        fprintf( stderr, "PSNPATH module: Invalid --STAT value: %i\n", iStatFlag);
        return 1;
      }

      iOptionFlag = 1;
    }
    
    else if (strncmp(pcInputLine, "--MODE", 6) == 0)
    {
      sscanf(pcInputLine, "--MODE %s", cMatchMode);
      if(strcmp(cMatchMode, "FULL") == 0)
        iModeFlag = 0;
        
      else if(strcmp(cMatchMode, "LITE") == 0)
        iModeFlag = 1;
        
      else
      {
        fprintf( stderr, "PSNPATH module: Invalid --MODE option: %s; Valid values are : FULL (default), LITE\n", cMatchMode);
        exit(5);
      }
      
      iOptionFlag = 1;
    }

    else if (strncmp(pcInputLine, "--TYPE", 6) == 0)
    {
      sscanf(pcInputLine, "--TYPE %s", cMatchMode);
      if(strcmp(cMatchMode, "RAW") == 0)
        iPSNTypeFlag = 0;
        
      else if(strcmp(cMatchMode, "AVG") == 0)
        iPSNTypeFlag = 1;
        
      else if(strcmp(cMatchMode, "MIX") == 0)
        iPSNTypeFlag = 2;
        
      else
      {
        fprintf( stderr, "PSNPATH module: Invalid --PSNTYPE option: %s; Valid values are : RAW (default), AVG or MIX\n", cMatchMode);
        exit(5);
      }
      
      iOptionFlag = 1;
    }

    else if (strncmp(pcInputLine, "--LOG", 5) == 0)
    {
      sscanf(pcInputLine, "--LOG %d", &iLogFlag);
      if(iLogFlag != 0 && iLogFlag != 1)
      {
        fprintf( stderr, "PSNPATH module: Invalid --LOG value: %i\n", iLogFlag);
        return 1;
      }
      
      iOptionFlag = 1;
    }

    else if (strncmp(pcInputLine, "--OFFSET", 8) == 0)
    {
      iOptionFlag = 1;
    }
    
    else if (strncmp(pcInputLine, "--MATCH", 7) == 0)
    {
      sscanf(pcInputLine, "--MATCH %s", cMatchMode);
      if(strcmp(cMatchMode, "CROSS") == 0)
        iMatchMode = 0;
        
      else if(strcmp(cMatchMode, "LINEAR") == 0)
        iMatchMode = 1;
        
      else
      {
        fprintf( stderr, "PSNPATH module: Invalid --MATCH option: %s; Valid values are : CROSS (default), LINEAR\n", cMatchMode);
        exit(5);
      }
      
      iOptionFlag = 1;
    }
    
    else if (strncmp(pcInputLine, "--MINLEN", 8) == 0)
    {
      sscanf(pcInputLine, "--MINLEN %d", &iMinLen);
      iOptionFlag = 1;
    }

    else if (strncmp(pcInputLine, "--MINFREQ", 9) == 0)
    {
      sscanf(pcInputLine, "--MINFREQ %f", &fMinFreq);
      iOptionFlag = 1;
    }

    else if (strncmp(pcInputLine, "--MEMPATH", 9) == 0)
    {
      sscanf(pcInputLine, "--MEMPATH %d", &MAXNUMOFPATHS);
      iOptionFlag = 1;
    }
    
    else if (strncmp(pcInputLine, "--MEMGAIN", 9) == 0)
    {
      sscanf(pcInputLine, "--MEMGAIN %d", &GAIN);
      iOptionFlag = 1;
    }
    
    else if (strncmp(pcInputLine, "--STARTFRAME", 12) == 0)
    {
      sscanf(pcInputLine, "--STARTFRAME %d", &iStartFrame);
      iOptionFlag = 1;
    }
    
    else if (strncmp(pcInputLine, "--STARTPAIR", 11) == 0)
    {
      sscanf(pcInputLine, "--STARTPAIR %d", &iStartPair);
      iOptionFlag = 1;
    }
    
    if(iOptionFlag == 0)
    {
      fprintf( stderr, "PSNPATH module: Could NOT understand option: %s\n", pcInputLine);
      return 1;
    }
    
    iInputLineNum++;
  }

  // Check some option values
  if(fIntMinCutOff < 0.0)
  {
    fprintf( stderr, "PSNPATH module: Invalid --IMIN value\n");
    return 1;
  }

  if(fCorrCutOff < -1.0 || fCorrCutOff > 1.0)
  {
    fprintf( stderr, "PSNPATH module: Invalid --CUTOFF value\n");
    return 1;
  }
  
  if(iModeFlag == 1)
  {
    // --MODE LITE
    if(iStatFlag == 1)
    {
      fprintf( stderr, "PSNPATH module: With --MODE LITE, --STAT option is ignored\n");
      iStatFlag = 0;
    }
    
    if(iLogFlag == 1)
    {
      fprintf( stderr, "PSNPATH module: With --MODE LITE, --LOG option is ignored\n");
      iLogFlag = 0;
    }
    
    if(iFrameFlag == 1)
    {
      fprintf( stderr, "PSNPATH module: With --MODE LITE, --FRAME option is ignored\n");
      iFrameFlag = 0;
    }
  }
  
  if(iMININTFREQOptFlag == 1)
    if(iPSNTypeFlag == 0)
      fprintf( stderr, "PSNPATH module: With --TYPE RAW, --MININTFREQ option is ignored\n");
  
  if(iAVGFILEOptFlag == 0)
  {
    if(iPSNTypeFlag == 2)
    {
      fprintf( stderr, "PSNPATH module: With --TYPE MIX, --AVGFILE option is needed\n");
      return 1;
    }
  }

  else if(iAVGFILEOptFlag == 1)
  {
    if(iPSNTypeFlag != 2)
    {
      fprintf( stderr, "PSNPATH module: With --TYPE RAW or AVG, --AVGFILE option is ignored\n");
    }
  }

  if(iLogFlag == 1)
  {
    sprintf(cLogFileName, "%s.log", cPathTitle);
    FLogFile = fopen(cLogFileName, "w");
    
    if(FLogFile == NULL)
    {
      fprintf( stderr, "PSNPATH module: Unable to open file: %s\n", cLogFileName);
      return 1;
    }
    
    fprintf(FLogFile, "#  Pair    Res1      Res2     BadFrames     NullFrames   TotPaths   HighFreqVal   HighFreqLen   HighFreqCor   HighFreqPath\n");
  }
  
  
  if(iPSNTypeFlag == 1)
    iCheckClustFlag = 0;
  
  // Check if only one residue is submitted as
  // first residue or if a file if specified
  pcTest1 = strstr(cTmpRes1, ".txt");
  pcTest2 = strstr(cTmpRes1, ".TXT");
  if(pcTest1 == NULL && pcTest2 == NULL)
  {
    iNumOfRes1 = 1;
    ppcRes1    = (char **) calloc( 1, sizeof(char *));
    piRes1List = (int   *) calloc( 1, sizeof(int   ));
    ppcRes1[0] = (char  *) calloc(50, sizeof(char  ));
    strcpy(ppcRes1[0], cTmpRes1);
  }
  else
  {
    // First residue is a file
    iRes1FileFlag = 1;
    FResFileHandler = fopen(cTmpRes1, "r");
    if(FResFileHandler == NULL)
    {
      fprintf( stderr, "PSNPATH module: Unable to open file: %s\n", cTmpRes1);
      return 1;
    }
    
    iNumOfResInFile = 0;
    while(fgets(cLine, 500, FResFileHandler) != NULL)
    {
      if(strlen(cLine) > 1)
        iNumOfResInFile++;
    }
    
    iNumOfRes1    = iNumOfResInFile;
    piRes1List    = (int   *) calloc(iNumOfResInFile, sizeof(int   ));
    ppcRes1       = (char **) calloc(iNumOfResInFile, sizeof(char *));
    for(ii=0; ii<iNumOfResInFile; ii++)
      ppcRes1[ii] = (char  *) calloc(50,              sizeof(char  ));
    
    rewind(FResFileHandler);
    ii=-1;
    while(fgets(cLine, 500, FResFileHandler) != NULL)
    {
      if(strcmp(cLine, "\n") != 0)
      {
        ii++;
        cLine[strlen(cLine)-1] = '\0';
        strcpy(ppcRes1[ii], cLine);
      }
    }
  }
  
  // Check if only one residue is submitted as
  // last residue  or if a file if specified
  pcTest1 = strstr(cTmpRes2, ".txt");
  pcTest2 = strstr(cTmpRes2, ".TXT");
  if(pcTest1 == NULL && pcTest2 == NULL)
  {
    // Last residue is a residue :)
    iNumOfRes2 = 1;
    piRes2List = (int   *) calloc( 1, sizeof(int   ));
    ppcRes2    = (char **) calloc( 1, sizeof(char *));
    ppcRes2[0] = (char  *) calloc(50, sizeof(char  ));
    strcpy(ppcRes2[0], cTmpRes2);
  }
  else
  {
    // Last residue is a file
    iRes2FileFlag = 1;
    FResFileHandler = fopen(cTmpRes2, "r");
    if(FResFileHandler == NULL)
    {
      fprintf( stderr, "PSNPATH module: Unable to open file: %s\n", cTmpRes2);
      return 1;
    }
    
    iNumOfResInFile = 0;
    while(fgets(cLine, 500, FResFileHandler) != NULL)
    {
      if(strlen(cLine) > 1)
        iNumOfResInFile++;
    }
    
    iNumOfRes2    = iNumOfResInFile;
    piRes2List    = (int   *) calloc(iNumOfResInFile, sizeof(int   ));
    ppcRes2       = (char **) calloc(iNumOfResInFile, sizeof(char *));
    for(ii=0; ii<iNumOfResInFile; ii++)
      ppcRes2[ii] = (char  *) calloc(50,              sizeof(char  ));
    
    rewind(FResFileHandler);
    
    ii=-1;
    while(fgets(cLine, 500, FResFileHandler) != NULL)
    {
      if(strcmp(cLine, "\n") != 0)
      {
        ii++;
        cLine[strlen(cLine)-1] = '\0';
        strcpy(ppcRes2[ii], cLine);
      }
    }
  }
  
  if(GetPSNParam(ppcInput, iOriginalInputLineNum) == 1)                 // Gets some psn parameters
    return 1;
  
  iNumOfPoxInt = (int) ((((float) iNumOfRes * (float) iNumOfRes) - (float) iNumOfRes) / 2.0);
  
  // *** Allocations ***
  
  if(iPSNTypeFlag == 1)
  {
    // allocates matrices for avg paths search
    
    ppfAvgIntStrengths   = (float **) calloc(iNumOfRes, sizeof(float *));
    ppfAvgIntFrequencies = (float **) calloc(iNumOfRes, sizeof(float *));
    
    for(ii=0; ii<iNumOfRes; ii++)
    {
      ppfAvgIntStrengths[ii]   = (float *) calloc(iNumOfRes, sizeof(float));
      ppfAvgIntFrequencies[ii] = (float *) calloc(iNumOfRes, sizeof(float));
    }
  
    ppcNodeLabel = (char **) calloc(iNumOfRes, sizeof(char *));
    for(ii=0; ii<iNumOfRes; ii++)
      ppcNodeLabel[ii] = (char *) calloc(50, sizeof(char));
  }
  
  if(iPSNTypeFlag == 2)
  {
    // mix paths search
    
    ppcNodeLabel             = (char  **) calloc(iNumOfRes, sizeof(char  *));
    ppfStableLink            = (float **) calloc(iNumOfRes, sizeof(float *));
    
    if(iPSNTypeFlag == 2 && fIntMinFreq == 0.0)
      ppfIntMatrixBackup       = (float **) calloc(iNumOfRes, sizeof(float *));
    
    for(ii=0; ii<iNumOfRes; ii++)
    {
      ppcNodeLabel[ii]       = (char  *) calloc(50,         sizeof(char   ));
      ppfStableLink[ii]      = (float *) calloc(iNumOfRes,  sizeof(float  ));
      if(iPSNTypeFlag == 2 && fIntMinFreq == 0.0)
        ppfIntMatrixBackup[ii] = (float *) calloc(iNumOfRes,  sizeof(float  ));
    }
    
    pfIntFreqValList = (float *) calloc(iNumOfPoxInt, sizeof(float));
  }
  
  if(iModeFlag == 0 || iModeFlag == 1)
  {
    // --MODE FULL & LITE
    piVisited     = (int    *) calloc(iNumOfRes, sizeof(int       ));   // Visited nodes
    
    // *-* SIZE *-*
    //pfPathLength  = (float  *) calloc(MAXNUMOFPATHS, sizeof(float )); // The length of the shortest path between iRes1 and iRes2
    pfPathLength  = (float  *) calloc(iNumOfRes, sizeof(float ));       // The length of the shortest path between iRes1 and iRes2
    
    ppiPathMatrix = (int   **) calloc(MAXNUMOFPATHS, sizeof(int  *));   // Used by GenPath
    for(ii=0; ii<MAXNUMOFPATHS; ii++)
      ppiPathMatrix[ii] = (int    *) calloc(iNumOfRes, sizeof(int ));
  }
  
  if(iModeFlag == 1)
  {
    // --MODE LITE
    pcPath          = (char   *) calloc(MAXNUMOFPATHS, sizeof(char  )); // Stores tmp path in lite mode
    
    if(iCheckClustFlag == 1)
      piClusters    = (int    *) calloc(iNumOfRes,     sizeof(int   )); // Stores the cluster of each node
  }

  if(iModeFlag == 0)
  {
    // --MODE FULL
    piTmpPathsLen  = (int    *) calloc(MAXNUMOFPATHS, sizeof(int    )); // Tmp Path Length
    pfTmpScore     = (float  *) calloc(MAXNUMOFPATHS, sizeof(float  )); // Tmp Path score
    ppfTmpCorrStat = (float **) calloc(MAXNUMOFPATHS, sizeof(float *)); // Tmp Corr Stat
    pfTmpWeights   = (float  *) calloc(MAXNUMOFPATHS, sizeof(float  )); // Tmp Weigths
    pcTmpPaths     = (char  **) calloc(MAXNUMOFPATHS, sizeof(char  *)); // Stores tmp paths
    
    for(ii=0; ii<MAXNUMOFPATHS; ii++)
    {
      ppfTmpCorrStat[ii] = (float  *) calloc(6,             sizeof(float ));
      pcTmpPaths[ii]     = (char   *) calloc(MAXNUMOFPATHS, sizeof(char));
    }
    
    piPathFreq    = (int   *) calloc(GAIN*MAXNUMOFPATHS, sizeof(int  ));// Path Frequencies
    piNodesInPath = (int   *) calloc(GAIN*MAXNUMOFPATHS, sizeof(int  ));// Number of nodes in selected path
    pfPathScores  = (float *) calloc(GAIN*MAXNUMOFPATHS, sizeof(float));// Scores of all path
    pfPathWeights = (float *) calloc(GAIN*MAXNUMOFPATHS, sizeof(float));// Scores of all path
  
    ppfCorrStatList = (float **) calloc(GAIN*MAXNUMOFPATHS, sizeof(float *));// Min, Max and Avg Correlation Statistics of all paths
    for(ii=0; ii<GAIN*MAXNUMOFPATHS; ii++)
    ppfCorrStatList[ii] = (float *) calloc(6, sizeof(float));
  
    pcPathList = (char **) calloc(GAIN*MAXNUMOFPATHS, sizeof(char *));    // List of paths
    for(ii=0; ii<GAIN*MAXNUMOFPATHS; ii++)
      pcPathList[ii] = (char *) calloc(MAXNUMOFPATHS, sizeof(char  ));
  
    // Memory Check
    if(pcPathList == NULL)
    {
      fprintf( stderr, "PSNPATH module: Not enough memory\nTry to decrease --MEMPATH and/or--MEMGAIN values\n");
      fprintf( stderr, "Actual values:\n--MEMPATH : %d\n--MEMGAIN : %d\n", MAXNUMOFPATHS, GAIN);
      return 1;
    }
  }
  // *******************
  
  if(iMatchMode == 1)
  {
    // Some sanity checks for --MATCH LINEAR option

    // Checks if both residues in --PAIR option are files
    if((iRes1FileFlag + iRes2FileFlag) != 2)
    {
      fprintf( stderr, "PSNPATH module: --MATCH LINEAR needs two files as --PAIR option\n");
      return 1;
    }
    
    // Checks if both files have the same number of lines
    if(iNumOfRes1 != iNumOfRes2)
    {
      fprintf( stderr, "PSNPATH module: --MATCH LINEAR needs two files with the same number of lines\n");
      return 1;
    }
  }

  if(SequenceCheck() == 1)                                              // Check seq in PSN & Corr files
    return 1;

  GetCorrData();                                                        // Get Corr data
  
  if(iPSNTypeFlag == 2)
  {
    LoadAvgLabIdx();                                                    // Loads node labels in psn avg file
    LoadStableLinks();                                                  // Loads link frequencies
  }
  
  if(iModeFlag == 0)
  {
    // --MODE FULL

    iPairCont = 0;
    
    if(iPSNTypeFlag == 1)
    {
      ClearIntMatrix();
      LoadAvgLabIdx();
      GetPSGData();
      iNumOfFrames = 1;
    }
    
    // Loops over all selected pairs
    if(iMatchMode == 0)
    {
      // --MATCH CROSS
      for(ii=0; ii<iNumOfRes1; ii++)
      {
        iRes1 = piRes1List[ii];
        
        for(jj=0; jj<iNumOfRes2; jj++)
        {

          iRes2 = piRes2List[jj];
          
          iPairCont++;
          
          if(iPairCont < iStartPair)
            continue;

          time(&time_Start);
          
          SetIntVar();                                                  // Set internal variables
          ClearAll();
          
          if(iCheckClustFlag == 1)
          {
            fBadFramesVal = GetBadFrames();
            if(fBadFramesVal > fFreqCutOff)                             // Check Bad frames
            {
              if(iLogFlag == 1)
                fprintf(FLogFile, "%7d   %7s   %7s      %6.2f     %s\n",
                iPairCont, pcSequence[iRes1], pcSequence[iRes2], fBadFramesVal, "***SKIPPED***");
              
              continue;
            }
          }
          
          // Generate output file name components
          sscanf(pcSequence[iRes1], "%[^:]:%1s%d", cTmpSegName, cTmpResCode, &iTmpResNum);
          sprintf(cTmpRes1, "%s_%s%d", cTmpSegName, cTmpResCode, iTmpResNum);
          
          sscanf(pcSequence[iRes2], "%[^:]:%1s%d", cTmpSegName, cTmpResCode, &iTmpResNum);
          sprintf(cTmpRes2, "%s_%s%d", cTmpSegName, cTmpResCode, iTmpResNum);
          
          sprintf(cOutPutFileName, "PSNPath-%s-%s.path", cTmpRes1, cTmpRes2);
          
          FOutFile = fopen(cOutPutFileName, "w");
          if(FOutFile == NULL)
          {
            fprintf( stderr, "PSNPATH module: Unable to open file: %s\n", cOutPutFileName);
            return 1;
          }

          if(iFrameFlag == 1)
          {
            sprintf(cFrameFileName, "PSNPath-%s-%s.frame", cTmpRes1, cTmpRes2);
            FFrameFile = fopen(cFrameFileName, "w");
            if(FFrameFile == NULL)
            {
              fprintf( stderr, "PSNPATH module: Unable to open file: %s\n", cFrameFileName);
              return 1;
            }
            else
              fprintf(FFrameFile, "#   Frame   Paths\n");
          }
          
          if(iStatFlag == 1)
          {
            sprintf(cStatFileName, "PSNPath-%s-%s.stat", cTmpRes1, cTmpRes2);
            FStatFile = fopen(cStatFileName, "w");
            if(FStatFile == NULL)
            {
              fprintf( stderr, "PSNPATH module: Unable to open file: %s\n", cStatFileName);
              return 1;
            }
            else
              fprintf(FStatFile, "#   Frame        Path      Length   AvgScore   ScoreMin   ScoreMax DeltaScore\n");
          }
          
          GetCorrRes();                                                 // Gets correlated residues
          rewind(FRawFile);
          
          if(iPSNTypeFlag == 0 || iPSNTypeFlag == 2)
          {
            if(iPSNTypeFlag == 2 && fIntMinFreq == 0.0)
              IntMatrixFilter();
            
            GetPaths();                                                 // Calculates paths
          }
          
          else if(iPSNTypeFlag == 1)
          {
            FindPath();
            GenPath(iRes1, 0);
            SavePath();
          }

          SaveData();

          fclose(FOutFile);
          
          if(iFrameFlag == 1)
            fclose(FFrameFile);

          if(iStatFlag == 1)
            fclose(FStatFile);
          
          rewind(FRawFile);
          
          ClearAll();
        }
      }
    }

    else
    {
      // --MATCH LINEAR
      for(ii=0; ii<iNumOfRes1; ii++)
      {
        iRes1 = piRes1List[ii];
        iRes2 = piRes2List[ii];
        iPairCont++;
        
        if(iPairCont < iStartPair)
          continue;
        
        time(&time_Start);
        
        SetIntVar();                                                    // Set internal variables
        ClearAll();
        
        if(iCheckClustFlag == 1)
        {
          fBadFramesVal = GetBadFrames();
          
          if(fBadFramesVal > fFreqCutOff)                               // Check Bad frames
          {
            if(iLogFlag == 1)
              fprintf(FLogFile, "%7d   %7s   %7s      %6.2f     %s\n",
              iPairCont, pcSequence[iRes1], pcSequence[iRes2], fBadFramesVal, "***SKIPPED***");
              
            continue;
          }
        }
        
        // Generate output file name components
        sscanf(pcSequence[iRes1], "%[^:]:%1s%d", cTmpSegName, cTmpResCode, &iTmpResNum);
        sprintf(cTmpRes1, "%s_%s%d", cTmpSegName, cTmpResCode, iTmpResNum);
        
        sscanf(pcSequence[iRes2], "%[^:]:%1s%d", cTmpSegName, cTmpResCode, &iTmpResNum);
        sprintf(cTmpRes2, "%s_%s%d", cTmpSegName, cTmpResCode, iTmpResNum);
        
        sprintf(cOutPutFileName, "PSNPath-%s-%s.path", cTmpRes1, cTmpRes2);
        
        FOutFile = fopen(cOutPutFileName, "w");
        if(FOutFile == NULL)
        {
          fprintf( stderr, "PSNPATH module: Unable to open file: %s\n", cOutPutFileName);
          return 1;
        }
        
        if(iFrameFlag == 1)
        {
          sprintf(cFrameFileName, "PSNPath-%s-%s.frame", cTmpRes1, cTmpRes2);
          FFrameFile = fopen(cFrameFileName, "w");
          if(FFrameFile == NULL)
          {
            fprintf( stderr, "PSNPATH module: Unable to open file: %s\n", cFrameFileName);
            return 1;
          }
          else
            fprintf(FFrameFile, "#   Frame   Paths\n");
        }
        
        if(iStatFlag == 1)
        {
          sprintf(cStatFileName, "PSNPath-%s-%s.stat", cTmpRes1, cTmpRes2);
          FStatFile = fopen(cStatFileName, "w");
          if(FStatFile == NULL)
          {
            fprintf( stderr, "PSNPATH module: Unable to open file: %s\n", cStatFileName);
            return 1;
          }
          else
            fprintf(FStatFile, "#   Frame        Path      Length   AvgScore   ScoreMin   ScoreMax DeltaScore\n");
        }
        
        GetCorrRes();                                                   // Gets correlated residues
        rewind(FRawFile);
        
        if(iPSNTypeFlag == 0 || iPSNTypeFlag == 2)
        {
          if(iPSNTypeFlag == 2 && fIntMinFreq == 0.0)
            IntMatrixFilter();
          
          GetPaths();                                                   // Calculates paths
        }
        
        else if(iPSNTypeFlag == 1)
        {
          FindPath();
          GenPath(iRes1, 0);
          SavePath();
        }
        
        SaveData();                                                     // Save paths to file
        fclose(FOutFile);
        
        if(iFrameFlag == 1)
          fclose(FFrameFile);

        if(iStatFlag == 1)
          fclose(FStatFile);
        
        rewind(FRawFile);
        ClearAll();
      }
    }
    
    if(iLogFlag == 1)
      fclose(FLogFile);
  }
  
  else if(iModeFlag == 1)
  {
    // --MODE LITE
    
    if(iPSNTypeFlag == 1)
    {
      ClearIntMatrix();
      LoadAvgLabIdx();
      GetPSGData();
      iNumOfFrames = 1;
    }
    
    // Open the first output file
    iFrameBlockNum = 1;
    sprintf(cFrameFileName, "PSNPath-%s_fb-%d.fblock", cPathTitle , iFrameBlockNum);
    FFrameFile = fopen(cFrameFileName, "a+");
    if(FFrameFile == NULL)
    {
      fprintf( stderr, "PSNPATH module: Unable to create file: %s\n", cFrameFileName);
      return 1;
    }
    
    iBlockCont = 0;
    rewind(FRawFile);
    SetIntVar();
    
    if(iPSNTypeFlag == 0 || iPSNTypeFlag == 2)
      ClearIntMatrix();
    
    iNumOfPaths = iNumOfRes;
    ClearAllLite();
    iFrameCont = 0;
    
    iPairCont = 0;
    
    while(1)
    {
      iFrameCont++;      
      if(iFrameCont > iNumOfFrames)
        break;

      if(iFrameCont >= iStartFrame)
      {
        if(iPSNTypeFlag == 0 || iPSNTypeFlag == 2)
          ClearIntMatrix();
        
        ClearAllLite();
      }
      
      if(iPSNTypeFlag == 0 || iPSNTypeFlag == 2)
        GetFrameInt();

      iBlockCont++;
      
      if(iBlockCont > iFramePerBlock)
      {
        fclose(FFrameFile);
        iBlockCont = 1;
        iFrameBlockNum++;
        
        sprintf(cFrameFileName, "PSNPath-%s_fb-%d.fblock", cPathTitle , iFrameBlockNum);
        FFrameFile = fopen(cFrameFileName, "a+");
        if(FFrameFile == NULL)
        {
          fprintf( stderr, "PSNPATH module: Unable to create file: %s\n", cFrameFileName);
          return 1;
        }
      }
      
      if(iFrameCont < iStartFrame)
        continue;
      
      fprintf(FFrameFile, "@f%d\n", iFrameCont);

      if(iMatchMode == 0)
      {
        // --MATCH CROSS
        iPairCont = 0;
        for(ii=0; ii<iNumOfRes1; ii++)
        {
          iRes1 = piRes1List[ii];
          
          for(jj=0; jj<iNumOfRes2; jj++)
          {
            iRes2 = piRes2List[jj];
            
            iPairCont++;
            
            // skip non-productive pairs
            iCalcFlag = 1;
            if(iCheckClustFlag == 1) {
              if(piClusters[iRes1] == piClusters[iRes2] && piClusters[iRes1] != 0)
                iCalcFlag = 1;
              else
                iCalcFlag = 0;
            }

            if(iCalcFlag == 1)
            {
              if(iPSNTypeFlag == 2 && fIntMinFreq == 0.0)
                IntMatrixFilter();
              
              ClearAllLite();
              SetIntVar();                                              // Set internal variables
              FindPath();
              GenPath(iRes1, 0);
              SavePathLite();
            }
          }
        }
      }
      
      else
      {
        // --MATCH LINEAR
        iPairCont = 0;
        for(ii=0; ii<iNumOfRes1; ii++)
        {
          iRes1 = piRes1List[ii];
          iRes2 = piRes2List[ii];
          
          iPairCont++;

          // skip non-productive pairs
          iCalcFlag = 1;
          if(iCheckClustFlag == 1) {
            if(piClusters[iRes1] == piClusters[iRes2] && piClusters[iRes1] != 0)
              iCalcFlag = 1;
            else
              iCalcFlag = 0;
          }

          if(iCalcFlag == 1)
          {
            if(iPSNTypeFlag == 2 && fIntMinFreq == 0.0)
              IntMatrixFilter();
            
            ClearAllLite();
            SetIntVar();                                                // Set internal variables
            FindPath();
            GenPath(iRes1, 0);
            SavePathLite();
          }
        }
      }
    }
    

    fclose(FFrameFile);
  }
  
  fclose(FRawFile);
  return 0;                                                             // All ok
}

int GetPSNParam(char **ppcInput, int iInputLineNum)
{
  /*
     Gets some psn parameters
     number of residues
     number of segments
     number of franes
     applies segment offsets
     check passed residues
     etc
  */

  int     ii, jj;                                                       // Some iterators
  int     iSegNameFlag=0;                                               // Flag used to check segname passed using -o option
  int     iFirstResNum, iMaxResNum;                                     // the first and the last residue number in a segment
  int     iSeleRes1, iSeleRes2;                                         // selected resnums
  int     iProgResNum;                                                  // used to load sequence
  int     iTmpResNum;                                                   // Res num used to add offest using -o option
  int     iSegNum, iResNum;
  int    *piSegFirstResNum;
  int    *piSegLastResNum;
  int     iLineNum;                                                     // Used to parse input file and find --OFFSET directive(s)
  int     iTmpSegOffSet;
  
  char    cTmpSegName[10], cTmpResType[5];                              // Res seg and type used to add offest using -o option
  char    cResCode[15];                                                 // used to read sequence in the form segname:resname+resnum
  char    cSeleSeg1[10], cSeleSeg2[10];                                 // selected segnames
  char    cLine[100], cJunk[10], cSeg[10];                              // string buffers
  
  iCalcStep = 1;
  
  
  // Reads number of segments, residues and frames from psn file header
  /*
  for(ii=0; ii<15; ii++)
    fgets(cLine, 100, FRawFile);
  
  printf("*-*|%s|\n", cLine);
  sscanf(cLine, "%s %s %s %s %d", cJunk, cJunk, cJunk, cJunk, &iNumOfSeg);
  
  fgets(cLine, 100, FRawFile);
  sscanf(cLine, "%s %s %s %s %d", cJunk, cJunk, cJunk, cJunk, &iNumOfRes);
  printf("*-*|%s|\n", cLine);
  
  fgets(cLine, 100, FRawFile);
  fgets(cLine, 100, FRawFile);
  sscanf(cLine, "%s %s %s %s %d", cJunk, cJunk, cJunk, cJunk, &iNumOfFrames);
  printf("*-*|%s|\n", cLine);
  */
  
  // Reads number of segments, residues and frames from psn file header
  while(1)
  {
    if(fgets(cLine, 100, FRawFile)==NULL && ( !feof(FRawFile) || ferror(FRawFile) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
    
    if(strncmp(cLine, "# Seg Num       : ", 18) == 0)
    {
      sscanf(cLine, "%s %s %s %s %d", cJunk, cJunk, cJunk, cJunk, &iNumOfSeg);
    }
    
    if(strncmp(cLine, "# Res Num       : ", 18) == 0)
    {
      sscanf(cLine, "%s %s %s %s %d", cJunk, cJunk, cJunk, cJunk, &iNumOfRes);
    }
    
    if(strncmp(cLine, "# Frame Num     : ", 18) == 0)
    {
      sscanf(cLine, "%s %s %s %s %d", cJunk, cJunk, cJunk, cJunk, &iNumOfFrames);
      break;
    }
  }
  
  while(1)
  {
    if(fgets(cLine, 100, FRawFile)==NULL)
    {
      fprintf( stderr, "PSNPATH module: Unable to read Segment names in file: %s\n", cRawPSNFileName);
      return 1;
    }
    
    if(strcmp(cLine, "*** Seg Info ***\n")==0)
      break;
  }
  
  piSegOffSet   = (int   *) calloc(iNumOfSeg,  sizeof(int   ));
  pcSegNameList = (char **) calloc(iNumOfSeg,  sizeof(char *));
  for(ii=0; ii<iNumOfSeg; ii++)
    pcSegNameList[ii] = (char *) calloc(10,    sizeof(char  ));
  
  piSegFirstResNum = (int *) calloc(iNumOfSeg, sizeof(int   ));
  piSegLastResNum  = (int *) calloc(iNumOfSeg, sizeof(int   ));

  for(ii=0; ii<iNumOfSeg; ii++)
  {
    if(fgets(cLine, 100, FRawFile)==NULL && ( !feof(FRawFile) || ferror(FRawFile) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
    sscanf(cLine, "%s %s %d %d", cJunk, cSeg, &iMaxResNum, &iFirstResNum);
    iSegNum=-1;
    for(jj=0; jj<iNumOfSeg; jj++)
    {
      if(strcmp(pcSegNameList[jj], cSeg)==0)
        iSegNum=jj;
    }
    if(iSegNum==-1)
      iSegNum=ii;
    
    strcpy(pcSegNameList[iSegNum], cSeg);
    piSegFirstResNum[ii]  =  iFirstResNum;
    piSegLastResNum[ii]  +=  iMaxResNum;
  }
  
  // Applies offset(s)
  iLineNum = iInputLineNum;
  while(strncmp(ppcInput[iLineNum], "END", 3) != 0)
  {
    if(strncmp(ppcInput[iLineNum], "--OFFSET", 8) == 0)
    {
      iSegNameFlag=0;
      sscanf(ppcInput[iLineNum], "--OFFSET %s %d", cTmpSegName, &iTmpSegOffSet);
      for(jj=0; jj<iNumOfSeg; jj++)
      {
        if(strcmp(pcSegNameList[jj], cTmpSegName)==0)
        {
          piSegOffSet[jj] = iTmpSegOffSet;
          iSegNameFlag=1;
          break;
        }
      }
      if(iSegNameFlag==0)
      {
        fprintf( stderr, "PSNPATH module: Invalid Segment names: %s\n", ppcInput[iLineNum]);
        return 1;
      }
    }
    iLineNum++;
  }
  
  
  if(iCheckClustFlag == 1)
    piBadFrames    = (int    *)  calloc(iNumOfFrames, sizeof(int    )); // Used to select productive frames
  
  ppfIntMatrix     = (float **) calloc(iNumOfRes,     sizeof(float *)); // PSN Interaction Map
  ppfRealIntMatrix = (float **) calloc(iNumOfRes,     sizeof(float *)); // PSN Interaction Map, real values
  pcSequence       = (char  **) calloc(iNumOfRes,     sizeof(char  *)); // Molecule Sequence
  pcRawSequence    = (char  **) calloc(iNumOfRes,     sizeof(char  *)); // Molecule Sequence
  ppiPrevNode      = (int   **) calloc(iNumOfRes,     sizeof(int   *)); // Previous res in a path
  
  ppfCorMatrix  = (float **) calloc(iNumOfRes,     sizeof(float *));    // Used to store Corr data
  for(ii=0; ii<iNumOfRes; ii++)
    ppfCorMatrix[ii] = (float *) calloc(iNumOfRes, sizeof(float  ));
    
  piCorrVector = (int *)  calloc(iNumOfRes,        sizeof(int    ));    // Used to store correlated res in --mode cp
  
  for(ii=0; ii<iNumOfRes; ii++)
  {
    ppfIntMatrix[ii]     = (float *) calloc(iNumOfRes, sizeof(float ));
    ppfRealIntMatrix[ii] = (float *) calloc(iNumOfRes, sizeof(float ));
    pcSequence[ii]       = (char  *) calloc(15,        sizeof(char  ));
    pcRawSequence[ii]    = (char  *) calloc(15,        sizeof(char  ));
    ppiPrevNode[ii]      = (int   *) calloc(iNumOfRes, sizeof(int   ));
  }
  
  while(1)
  {
    if(fgets(cLine, 100, FRawFile)==NULL)
    {
      fprintf( stderr, "PSNPATH module: Unable to laod sequence from file: %s\n", cRawPSNFileName);
      return 1;
    }
    
    if(strcmp(cLine, "*** Seq ***\n")==0)
      break;
  }
  
  for(ii=0; ii<iNumOfRes; ii++)
  {
    if(fgets(cLine, 100, FRawFile)==NULL && ( !feof(FRawFile) || ferror(FRawFile) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
    
    // reading a raw file
    if(iPSNTypeFlag == 0 || iPSNTypeFlag == 2)
      sscanf(cLine, "%d %s", &iProgResNum, cResCode);
    
    // reading an avg file
    else if(iPSNTypeFlag == 1)
      sscanf(cLine, "%d %s %s", &iProgResNum, cJunk, cResCode);
    
    strcpy(pcRawSequence[iProgResNum-1], cResCode);
    sscanf(cResCode, "%[^:]:%1s%d", cTmpSegName, cTmpResType, &iTmpResNum);
    
    for(jj=0; jj<iNumOfSeg; jj++)
    {
      if(strcmp(pcSegNameList[jj], cTmpSegName)==0)
      {
        iTmpResNum=iTmpResNum+piSegOffSet[jj];
        cResCode[0]='\0';
        sprintf(cResCode, "%s:%s%d", cTmpSegName, cTmpResType, iTmpResNum);
      }
    }
    
    strcpy(pcSequence[iProgResNum-1], cResCode);
  }

  // Checking passed residues
  for(ii=0; ii<iNumOfRes1; ii++)
  {
    sscanf(ppcRes1[ii], "%[^:]:%d", cSeleSeg1, &iSeleRes1);
    iTmpResNum=0;
    for(jj=0; jj<iNumOfRes; jj++)
    {
      sscanf(pcSequence[jj], "%[^:]:%1s%d", cSeg, cJunk, &iResNum);
      if(strcmp(cSeg, cSeleSeg1) == 0 && iResNum == iSeleRes1)
      {
        piRes1List[ii] = jj;
        iTmpResNum = 1;
        break;
      }
    }
    
    if(iTmpResNum == 0)
    {
      fprintf( stderr, "PSNPATH module: Invalid residue number: %s\n", ppcRes1[ii]);
      return 1;
    }
  }

  for(ii=0; ii<iNumOfRes2; ii++)
  {
    sscanf(ppcRes2[ii], "%[^:]:%d", cSeleSeg2, &iSeleRes2);
    iTmpResNum = 0;
    for(jj=0; jj<iNumOfRes; jj++)
    {
      sscanf(pcSequence[jj], "%[^:]:%1s%d", cSeg, cJunk, &iResNum);
      if(strcmp(cSeg, cSeleSeg2) == 0 && iResNum == iSeleRes2)
      {
        piRes2List[ii] = jj;
        iTmpResNum = 1;
        break;
      }
    }
    
    if(iTmpResNum == 0)
    {
      fprintf( stderr, "PSNPATH module: Invalid residue number: %s\n", ppcRes2[ii]);
      return 1;
    }
  }
  
  iCalcStep = 0;
  
  return 0;
}

int SequenceCheck()
{
  /*
     Check Sequence in PSN and Corr files
  */
  
  int     iTmpRes1=0, iTmpRes2=0, iLastRes=0;
  
  char    cLine[200], cJunk[20], cResName[20];
  
  float   fCorr=0.0;
  
  iCalcStep = 2;
  
  cLine[0]='\0';
  cJunk[0]='\0';
  
  iLastRes = -1;
  
  while(fgets(cLine, 200, FCorFile)!=NULL)                              // skip header
  {
    if(strncmp(cLine, "#", 1)!=0)
    {
      sscanf(cLine, " %d %d %s %s %f", &iTmpRes1, &iTmpRes2, cJunk, cResName, &fCorr);
      
      iTmpRes2--;
      
      if(iLastRes > iTmpRes2)
        break;
      
      else
      {
        iLastRes = iTmpRes2;
        if(strcmp(pcRawSequence[iTmpRes2], cResName) != 0)
        {
          
          fprintf(stderr, "PSNPATH module: %s and %s sequences are different\n", cRawPSNFileName, cCorrFileName);
          fprintf(stderr, "%s != %s\n", pcRawSequence[iTmpRes2], cResName);
          return 1;
        }
      }
    }
  }
  
  iCalcStep = 0;
  
  rewind(FCorFile);
  
  return 0;
}

void GetCorrData()
{
  /*
     Reads Corr data from Corr file
  */
  
  int     ii, jj;                                                       // Some iterators
  int     iTmpRes1=0, iTmpRes2=0;
  
  char    cLine[200], cJunk[200];
  
  float   fCorr=0.0;
  
  iCalcStep = 3;
  
  cLine[0] = '\0';
  cJunk[0] = '\0';
  
  for(ii=0; ii<iNumOfRes; ii++)
    for(jj=0; jj<iNumOfRes; jj++)
      ppfCorMatrix[ii][jj]=0.0;
  
  while(fgets(cLine, 200, FCorFile)!=NULL)                              // skip header
  {
    if(strncmp(cLine, "#", 1)!=0)
    {
      sscanf(cLine, " %d %d %s %s %f", &iTmpRes1, &iTmpRes2, cJunk, cJunk, &fCorr);
      
      if(iAbsCorrFlag == 1)
        fCorr = fabs(fCorr);

      ppfCorMatrix[iTmpRes1-1][iTmpRes2-1]=fCorr;
    }
  }
  
  fclose(FCorFile);
  
  iCalcStep = 0;
}

void SetIntVar()
{
  /*
     Set internal variables used in paths computations
  */
  
  int     ii;
  
  iCalcStep = 4;
  
  if(iModeFlag == 0)
  {
    for(ii=0; ii<iNumOfPaths+1; ii++)
    {
      piNodesInPath[ii]      = 0;
      piPathFreq[ii]         = 0;
      pfPathScores[ii]       = 0.0;
      pcPathList[ii][0]      = '\0';
      ppfCorrStatList[ii][0] = 0.0;
      ppfCorrStatList[ii][1] = 0.0;
      ppfCorrStatList[ii][2] = 0.0;
      ppfCorrStatList[ii][3] = 0.0;
      ppfCorrStatList[ii][4] = 0.0;
      ppfCorrStatList[ii][5] = 0.0;
    }
  }
  
  iNumOfPaths     = 0;
  iNumOfGenPaths  = 0;
  iNumOfCorrRes   = 0;
  iNumOfBadFrames = 0;
  iLastDepth      = 0;
  iNumOfEdges     = 0;
  
  iCalcStep = 0;
  
}

float GetBadFrames()
{
  /* 
     Get non-productive frames
     
     Check if in a frame iRes1 and iRes2 are in the same cluster
      
     iCheckSum is 1 if iRes1 and iRes2 are in different clusters
     and than the frame is not productive and so that frame will
     not be processed by FindPath
     
  */

  int     ii, jj;
  int     iResNum=0, iResNumLen=-1, iFrmNum=-1, iTmpFrame=0;
  int     iCheckSum=0, iLastFrame=-1;
  
  float   fImin=-9.9;
  
  char    cLine[5000], cJunk[100], cBuffer[10];
  
  iCalcStep = 6;
  
  strcpy(cLine, "\0");
  
  rewind(FRawFile);
  
  iNumOfBadFrames = 0;
  
  while(1)
  {
    if(fgets(cLine, 5000, FRawFile)==NULL)
      break;
      
    if(strncmp(cLine, "nFr: ", 5)==0)
    {
      sscanf(cLine, "%s%d", cJunk, &iTmpFrame);
      iCheckSum=0;
      iFrmNum++;
    }

    else if(strncmp(cLine, "Imin:", 5)==0)
      sscanf(cLine, "%s%f", cJunk, &fImin);
    
    else if(strncmp(cLine, "C", 1)==0 && fImin==fIntMinCutOff && iLastFrame!=iFrmNum)
    {
      for(ii=0; ii<strlen(cLine); ii++)
      {
        if(cLine[ii]==':')
        {
          iResNumLen=-1;
          for(jj=ii+2; jj<strlen(cLine); jj++)
          {
            if(cLine[jj]!=' ')
            {
              iResNumLen++;
              cBuffer[iResNumLen]=cLine[jj];
            }
            
            else
            {
              iResNumLen++;
              cBuffer[iResNumLen]='\0';
              iResNumLen=-1;
              sscanf(cBuffer, "%d", &iResNum);
              iResNum--;
              
              if(iResNum == iRes1)
                iCheckSum++;
                
              if(iResNum == iRes2)
                iCheckSum++;
                
            }
          }
        }
      }
      
      if(iCheckSum == 1)
      {
        iCheckSum = 0;
        piBadFrames[iFrmNum] = 1;
        iLastFrame = iFrmNum;
        iNumOfBadFrames++;
      }
      
      else
      {
        iCheckSum = 0;
        piBadFrames[iFrmNum] = 0;
      }
    }
    
  }
  rewind(FRawFile);
  
  if(iNumOfBadFrames>0 && iModeFlag == 0)
  {
    iNumOfPaths=1;
    strcpy(pcPathList[0], "");
    strcpy(pcPathList[0], "[NULL_PATH]");
    piPathFreq[0]=iNumOfBadFrames;
    piNodesInPath[0]=0;
    pfPathScores[0]=0.0;
  }
  
  iCalcStep = 0;
  
  return ((iNumOfBadFrames*100)/(float)iNumOfFrames);
}

void GetCorrRes()
{
  /*
      Gets Correlated Residues
  */
  
  int     ii;                                                           // Some iterators

  iCalcStep = 7;

  iNumOfCorrRes = -1;
  
  for(ii=0; ii<iNumOfRes; ii++)
  {
    if(ppfCorMatrix[iRes1][ii] >= fCorrCutOff || ppfCorMatrix[iRes2][ii] >= fCorrCutOff)
    {
      iNumOfCorrRes++;
      piCorrVector[ii] = 1;
    }
    
    else
      piCorrVector[ii] = 0;
  }
  
  iNumOfCorrRes++;
  
  iCalcStep = 0;
}

void LoadAvgLabIdx()
{
  // Loads Node labels and indexes
  
  int     iNodeIdx;
  
  char    cLine[1000];
  char    cNodeLab1[30], cNodeLab2[30];
  FILE   *FileHandler;

  iCalcStep = 18;
  
  if(iPSNTypeFlag == 2)
    FileHandler = FAvgFile;
  else  
    FileHandler = FRawFile;
  
  rewind(FileHandler);
  
  while(1)
  {
    if(fgets(cLine, 1000, FileHandler)==NULL && ( !feof(FileHandler) || ferror(FileHandler) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
    if(strncmp(cLine, "*** Seq ***", 11) == 0)
      break;
  }
  
  while(1)
  {
    if(fgets(cLine, 1000, FileHandler)==NULL && ( !feof(FileHandler) || ferror(FileHandler) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
    
    if(strncmp(cLine, "========================", 24) == 0)
      break;
    
    else
    {
      sscanf(cLine, "%d %s %s", &iNodeIdx, cNodeLab1, cNodeLab2);
      strcat(ppcNodeLabel[iNodeIdx - 1], cNodeLab1);
    }
  }
  
  iCalcStep = 0;
}

int NodeLabToIndex(char *cNodeLab)
{
  
  int     ii;
  
  for(ii=0; ii<iNumOfRes; ++ii)
  {
    if(strcmp(cNodeLab, ppcNodeLabel[ii]) == 0)
    {
      return ii;
    }
  }
  
  printf("PSNPATH Error: unmapped node: %s\n", cNodeLab);
  exit(1);
}

void GetPSGData()
{
  // Get PSG data from an avg file
  
  int     ii, jj;
  int     iTmpRes1, iTmpRes2, iJunk;
  
  char    cLine[1000];
  char    cTmpRes1[20], cTmpRes2[20];
  
  float   fTmpImin, fTmpFreq, fTmpAvgInt;
  
  iCalcStep = 17;
  

  for(ii=0; ii<iNumOfRes; ++ii)
  {
    for(jj=0; jj<iNumOfRes; ++jj)
    {
      ppfIntMatrix[ii][jj]     = -1.0;
    }
  }
  
  while(1)
  {
    if(fgets(cLine, 1000, FRawFile)==NULL && ( !feof(FRawFile) || ferror(FRawFile) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
    if(strncmp(cLine, "*** Stable Residue Interactions ***", 35) == 0)
    {
      // skip header
      if(fgets(cLine, 1000, FRawFile)==NULL && ( !feof(FRawFile) || ferror(FRawFile) ))
      {
        fprintf(stderr, "Warning! Premature end of file reached!\n");
      }
      break;
    }
  }
  
  while(1)
  {
    if(fgets(cLine, 1000, FRawFile)==NULL && ( !feof(FRawFile) || ferror(FRawFile) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
    
    sscanf(cLine, "%f %s %s %d %f %f", &fTmpImin, cTmpRes1, cTmpRes2, &iJunk, &fTmpFreq, &fTmpAvgInt);
    
    if(fTmpImin == fIntMinCutOff)
    {
      if(fTmpFreq >= fIntMinFreq)
      {
        iTmpRes1 = NodeLabToIndex(cTmpRes1);
        iTmpRes2 = NodeLabToIndex(cTmpRes2);

        ppfAvgIntFrequencies[iTmpRes1][iTmpRes2] = fTmpFreq;
        ppfAvgIntFrequencies[iTmpRes2][iTmpRes1] = fTmpFreq;
        
        ppfAvgIntStrengths[iTmpRes1][iTmpRes2]   = fTmpAvgInt;
        ppfAvgIntStrengths[iTmpRes2][iTmpRes1]   = fTmpAvgInt;
        
        ppfRealIntMatrix[iTmpRes1][iTmpRes2]     = fTmpAvgInt;
        ppfRealIntMatrix[iTmpRes2][iTmpRes1]     = fTmpAvgInt;
        
        if(iWeightFlag == 0)
        {
          ppfIntMatrix[iTmpRes1][iTmpRes2] = 1.0;
          ppfIntMatrix[iTmpRes2][iTmpRes1] = 1.0;
        }
        
        else
        {
          ppfIntMatrix[iTmpRes1][iTmpRes2] = (100 - ppfAvgIntStrengths[iTmpRes1][iTmpRes2]);
          ppfIntMatrix[iTmpRes2][iTmpRes1] = (100 - ppfAvgIntStrengths[iTmpRes2][iTmpRes1]);
        }
      }
      
      break;
    }
  }
  
  while(1)
  {
    if(fgets(cLine, 1000, FRawFile)==NULL && ( !feof(FRawFile) || ferror(FRawFile) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
    
    if(strncmp(cLine, "========================", 24) == 0)
      break;
    
    sscanf(cLine, "%f %s %s %d %f %f", &fTmpImin, cTmpRes1, cTmpRes2, &iJunk, &fTmpFreq, &fTmpAvgInt);
    
    if(fTmpImin == fIntMinCutOff)
    {
      if(fTmpFreq >= fIntMinFreq)
      {
        iTmpRes1 = NodeLabToIndex(cTmpRes1);
        iTmpRes2 = NodeLabToIndex(cTmpRes2);

        ppfAvgIntFrequencies[iTmpRes1][iTmpRes2] = fTmpFreq;
        ppfAvgIntFrequencies[iTmpRes2][iTmpRes1] = fTmpFreq;
        
        ppfAvgIntStrengths[iTmpRes1][iTmpRes2]   = fTmpAvgInt;
        ppfAvgIntStrengths[iTmpRes2][iTmpRes1]   = fTmpAvgInt;
        
        ppfRealIntMatrix[iTmpRes1][iTmpRes2]     = fTmpAvgInt;
        ppfRealIntMatrix[iTmpRes2][iTmpRes1]     = fTmpAvgInt;
        
        if(iWeightFlag == 0)
        {
          ppfIntMatrix[iTmpRes1][iTmpRes2] = 1.0;
          ppfIntMatrix[iTmpRes2][iTmpRes1] = 1.0;
        }
        
        else
        {
          ppfIntMatrix[iTmpRes1][iTmpRes2] = (100 - ppfAvgIntStrengths[iTmpRes1][iTmpRes2]);
          ppfIntMatrix[iTmpRes2][iTmpRes1] = (100 - ppfAvgIntStrengths[iTmpRes2][iTmpRes1]);
        }
      }
    }
    
    else
      break;
  }
  
  iCalcStep = 0;

}

void GetUniqIntFreqList()
{
  int     ii, jj, kk;
  int     iAddValueFlag;
  
  for(ii=0; ii<iNumOfRes; ++ii)
  {
    for(jj=ii+1; jj<iNumOfRes; ++jj)
    {
      iAddValueFlag = 1;
      for(kk=0; kk<iNumOfPoxInt; ++kk)
      {
        if(pfIntFreqValList[kk] == 0)
          break;
          
        if(ppfStableLink[ii][jj] == pfIntFreqValList[kk])
        {
          iAddValueFlag = 0;
          break;
        }
      }
      
      if(iAddValueFlag == 1)
      {
        pfIntFreqValList[kk] = ppfStableLink[ii][jj];
      }
    }
  }

  qsort(pfIntFreqValList, iNumOfPoxInt, sizeof(float), ArraySorter);
}

int ArraySorter(const void *valueA, const void *valueB)
{
  //return (*(float *) valueB - *(float *) valueA);
  return (*(float *) valueB >= *(float *) valueA) ? 1 : -1;
}

void LoadStableLinks()
{
  // Get stable links from psn avg file
  
  int     ii, jj;
  int     iTmpRes1, iTmpRes2, iJunk;
  
  char    cLine[1000];
  char    cTmpRes1[20], cTmpRes2[20];
  
  float   fTmpImin, fTmpFreq, fTmpAvgInt;
  
  iCalcStep = 19;
  
  for(ii=0; ii<iNumOfRes; ++ii)
  {
    for(jj=0; jj<iNumOfRes; ++jj)
    {
      ppfIntMatrix[ii][jj]     = -1.0;
    }
  }
  
  while(1)
  {
    if(fgets(cLine, 1000, FAvgFile)==NULL && ( !feof(FAvgFile) || ferror(FAvgFile) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
    if(strncmp(cLine, "*** Stable Residue Interactions ***", 35) == 0)
    {
      // skip header
      if(fgets(cLine, 1000, FAvgFile)==NULL && ( !feof(FAvgFile) || ferror(FAvgFile) ))
      {
        fprintf(stderr, "Warning! Premature end of file reached!\n");
      }
      break;
    }
  }
  
  while(1)
  {
    if(fgets(cLine, 1000, FAvgFile)==NULL && ( !feof(FAvgFile) || ferror(FAvgFile) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
    sscanf(cLine, "%f %s %s %d %f %f", &fTmpImin, cTmpRes1, cTmpRes2, &iJunk, &fTmpFreq, &fTmpAvgInt);
    if(fTmpImin == fIntMinCutOff)
    {

      iTmpRes1 = NodeLabToIndex(cTmpRes1);
      iTmpRes2 = NodeLabToIndex(cTmpRes2);
      
      ppfStableLink[iTmpRes1][iTmpRes2] = fTmpFreq;
      ppfStableLink[iTmpRes2][iTmpRes1] = fTmpFreq;
      
      break;
    }
  }
  
  while(1)
  {
    if(fgets(cLine, 1000, FAvgFile)==NULL && ( !feof(FAvgFile) || ferror(FAvgFile) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
    
    if(strncmp(cLine, "========================", 24) == 0)
      break;
    
    sscanf(cLine, "%f %s %s %d %f %f", &fTmpImin, cTmpRes1, cTmpRes2, &iJunk, &fTmpFreq, &fTmpAvgInt);
    
    if(fTmpImin == fIntMinCutOff)
    {
      iTmpRes1 = NodeLabToIndex(cTmpRes1);
      iTmpRes2 = NodeLabToIndex(cTmpRes2);
      ppfStableLink[iTmpRes1][iTmpRes2] = fTmpFreq;
      ppfStableLink[iTmpRes2][iTmpRes1] = fTmpFreq;
    }
    
    else
      break;
  }
  
  //GetUniqIntFreqList();
  iCalcStep = 0;
}

void GetPaths()
{
  /*
     Reads res-res interactions from raw_psn_file and computes paths
  */

  int     iTmpRes1, iTmpRes2, iTmpFrame;
  char    cMagic[20];
  char    cLine[1000];
  float   fTmpRIS, fTmpHIS;
  
  iCalcStep = 8;
  
  iFrameCont=0;
  while(iFrameCont<iNumOfFrames)
  {
    while(1)
    {
      if(fgets(cLine, 1000, FRawFile)==NULL && ( !feof(FRawFile) || ferror(FRawFile) ))
      {
        fprintf(stderr, "Warning! Premature end of file reached!\n");
      }
      if(strncmp(cLine, "nFr:", 4) == 0)
      {
        sscanf(cLine, "%s%d", cMagic, &iTmpFrame);
        iFrameCont++;
        break;
      }
    }
    
    if(iCheckClustFlag == 1)
    {
      if(piBadFrames[iFrameCont-1] == 1)
      {
        if(iFrameFlag == 1)
          fprintf(FFrameFile, "%9d   [NULL_PATH]\n", iFrameCont);
          
        continue;
      }
    }
    
    if(fgets(cLine, 1000, FRawFile)==NULL && ( !feof(FRawFile) || ferror(FRawFile) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
    if(fgets(cLine, 1000, FRawFile)==NULL && ( !feof(FRawFile) || ferror(FRawFile) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
    
    iNumOfEdges=0;
    ClearAll();
    iCalcStep = 8;

    while(1)
    {
      if(fgets(cLine, 100, FRawFile)==NULL && ( !feof(FRawFile) || ferror(FRawFile) ))
      {
        fprintf(stderr, "Warning! Premature end of file reached!\n");
      }
      sscanf(cLine, "%s %d %d %f %f", cMagic, &iTmpRes1, &iTmpRes2, &fTmpRIS, &fTmpHIS);
      
      if(strcmp(cMagic, "&")!=0)
        break;
      
      iTmpRes1--;
      iTmpRes2--;
      
      if(fTmpRIS>=fIntMinCutOff)
      {
        //printf("%5d   %5d   %5d   %5d   %7.3f   %7.3f   %7.3f   %7.3f\n", iRes1, iRes2, iTmpRes1, iTmpRes2, ppfCorMatrix[iRes1][iTmpRes1], ppfCorMatrix[iRes2][iTmpRes1], ppfCorMatrix[iRes1][iTmpRes2], ppfCorMatrix[iRes2][iTmpRes2]);
        
        if(iPSNTypeFlag == 2)
        {
          if(ppfStableLink[iTmpRes1][iTmpRes2] < fIntMinFreq)
            continue;
        }
        
        ppfRealIntMatrix[iTmpRes1][iTmpRes2] = fTmpRIS;
        
        if(iWeightFlag == 0)
        {
          ppfIntMatrix[iTmpRes1][iTmpRes2] = 1.0;
          if(iPSNTypeFlag == 2 && fIntMinFreq == 0.0)
            ppfIntMatrixBackup[iTmpRes1][iTmpRes2] = ppfIntMatrix[iTmpRes1][iTmpRes2];
        }
        
        else if(iWeightFlag == 1)
        {
          ppfIntMatrix[iTmpRes1][iTmpRes2] = (100.0 - fTmpRIS);
          if(iPSNTypeFlag == 2 && fIntMinFreq == 0.0)
            ppfIntMatrixBackup[iTmpRes1][iTmpRes2] = ppfIntMatrix[iTmpRes1][iTmpRes2];
        }
        
        else if(iWeightFlag == 2)
        {
          ppfIntMatrix[iTmpRes1][iTmpRes2] = (100.0/fTmpRIS) * (100.0/fTmpRIS);
          if(iPSNTypeFlag == 2 && fIntMinFreq == 0.0)
            ppfIntMatrixBackup[iTmpRes1][iTmpRes2] = ppfIntMatrix[iTmpRes1][iTmpRes2];
        }
      }
    }
    
    if(iPSNTypeFlag == 2 && fIntMinFreq == 0.0)
      IntMatrixFilter();
      
    FindPath();
    GenPath(iRes1, 0);
    SavePath();
    iCalcStep = 0;
  }
}

void GetFrameInt()
{
  /*
     Reads res-res interactions from raw_psn_file
  */
  
  int     ii, jj;
  int     iTmpRes1=0, iTmpRes2=0;
  int     iResNum=0;
  int     iClsNum=0;
  int     iResNumLen=0;
  int     iAtm;
  
  char    cMagic[20];
  char    cLine[5000];
  char    cJunk[100], cBuffer[10];
  
  float   fTmpRIS=0.0, fTmpHIS=0.0, fImin=-9.9;
  
  iCalcStep = 15;
  
  // Reach interaction section
  while(1)
  {
    if(fgets(cLine, 5000, FRawFile)==NULL && ( !feof(FRawFile) || ferror(FRawFile) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
    if(strncmp(cLine, ">INT", 4) == 0)
      break;
  }
  if(fgets(cLine, 5000, FRawFile)==NULL && ( !feof(FRawFile) || ferror(FRawFile) ))
  {
    fprintf(stderr, "Warning! Premature end of file reached!\n");
  }
  // Load interactions
  iNumOfEdges=0;
  while(1)
  {
    if(fgets(cLine, 5000, FRawFile)==NULL && ( !feof(FRawFile) || ferror(FRawFile) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
    sscanf(cLine, "%s %d %d %f %f %d", cMagic, &iTmpRes1, &iTmpRes2, &fTmpRIS, &fTmpHIS, &iAtm);
    
    if(strcmp(cMagic, "&")!=0)
      break;
    
    iTmpRes1--;
    iTmpRes2--;
    
    //printf("*-* [%d] {%d : %d} = |%.2f| <-> |%d|\n", iFrameCont, iTmpRes1, iTmpRes2, fTmpRIS, ppfStableLink[iTmpRes1][iTmpRes2]);
    if(fTmpRIS>=fIntMinCutOff)
    {

      if(iPSNTypeFlag == 2)
      {
        if(ppfStableLink[iTmpRes1][iTmpRes2] < fIntMinFreq)
          continue;
      }

      ppfRealIntMatrix[iTmpRes1][iTmpRes2] = fTmpRIS;
      
      if(iWeightFlag == 0)
      {
        ppfIntMatrix[iTmpRes1][iTmpRes2] = 1.0;
        if(iPSNTypeFlag == 2 && fIntMinFreq == 0.0)
          ppfIntMatrixBackup[iTmpRes1][iTmpRes2] = ppfIntMatrix[iTmpRes1][iTmpRes2];
      }
      
      else
      {
        if(iPSNTypeFlag != 2)
        {
          if(iWeightFlag == 1)
            ppfIntMatrix[iTmpRes1][iTmpRes2] = (100.0 - fTmpRIS);
          
          else if(iWeightFlag == 2)
            ppfIntMatrix[iTmpRes1][iTmpRes2] = (100.0/fTmpRIS) * (100.0/fTmpRIS);
        }
        
        if(iPSNTypeFlag == 2)
          ppfIntMatrix[iTmpRes1][iTmpRes2] = (100.0/ppfStableLink[iTmpRes1][iTmpRes2]) * (100.0/ppfStableLink[iTmpRes1][iTmpRes2]);
        
        if(iPSNTypeFlag == 2 && fIntMinFreq == 0.0)
          ppfIntMatrixBackup[iTmpRes1][iTmpRes2] = ppfIntMatrix[iTmpRes1][iTmpRes2];
      }
    }
  }
  
  if(iCheckClustFlag == 1)
  {
    // Load cluster id of each node
    while(1)
    {
      if(feof(FRawFile))
        break;
      
      if(fgets(cLine, 5000, FRawFile)==NULL && ( !feof(FRawFile) || ferror(FRawFile) ))
      {
        fprintf(stderr, "Warning! Premature end of file reached!\n");
      }
      if(strncmp(cLine, "nFr: ", 5) == 0)
        break;

      else if(strncmp(cLine, "Imin:", 5)==0)
        sscanf(cLine, "%s%f", cJunk, &fImin);
      
      else if(strncmp(cLine, "C", 1) == 0 && fImin == fIntMinCutOff)
      {
        sscanf(cLine, "C %d: ", &iClsNum);
        for(ii=0; ii<strlen(cLine); ii++)
        {
          if(cLine[ii]==':')
          {
            iResNumLen=-1;
            for(jj=ii+2; jj<strlen(cLine); jj++)
            {
              if(cLine[jj]!=' ')
              {
                iResNumLen++;
                cBuffer[iResNumLen]=cLine[jj];
              }
              
              else
              {
                iResNumLen++;
                cBuffer[iResNumLen]='\0';
                iResNumLen=-1;
                sscanf(cBuffer, "%d", &iResNum);
                iResNum--;
                piClusters[iResNum] = iClsNum;
              }
            }
          }
        }
      }
    }
  }
  
  iCalcStep = 0;
}

void IntMatrixFilter()
{
  int     ii, jj, kk;                                                   // Some iterators
  int     iLastCalcStep;
  
  iLastCalcStep = iCalcStep;
  iCalcStep     = 20;
  
  for(kk=90; kk>=0; kk=kk-10)
  {
    for(ii=0; ii<iNumOfRes; ++ii)
    {
      for(jj=ii; jj<iNumOfRes; ++jj)
      {
        if(ppfIntMatrixBackup[ii][jj] != -1.0 && ppfStableLink[ii][jj] >= kk)
        {
          ppfIntMatrix[ii][jj] = ppfIntMatrixBackup[ii][jj];
          ppfIntMatrix[jj][ii] = ppfIntMatrixBackup[jj][ii];
        }
        
        else
        {
          ppfIntMatrix[ii][jj] = -1.0;
          ppfIntMatrix[jj][ii] = -1.0;
        }
      }
    }
    
    if(ClustCheck() == 1)
    {
      //printf("*-* %d %d %d %d", iFrameCont, iRes1+1, iRes2+1, kk);
      break;
    }
  }
  
  iCalcStep = iLastCalcStep;
}

int ClustCheck()
{
  int     ii, jj, kk;                                                   // Some iterators
  int     iLastClust=0;
  int    *piNodeClust;
  int     iOldClust;
  int     iLastCalcStep;
  int     iSameClustFlag;
  
  iLastCalcStep = iCalcStep;
  iCalcStep     = 21;
  
  piNodeClust = (int *) calloc(iNumOfRes, sizeof(int));
  
  for(ii=0; ii<iNumOfRes; ++ii)
  {
    for(jj=0; jj<iNumOfRes; ++jj)
    {
      if(ppfIntMatrix[ii][jj] != -1.0)
      {
        if(piNodeClust[ii] == piNodeClust[jj] && piNodeClust[ii] != 0)
        {
          // they are already in the same cluster, nothing to do
        }
        
        else if(piNodeClust[ii] == piNodeClust[jj] && piNodeClust[ii] == 0)
        {
          // they are both not yet clusterized
          iLastClust++;
          piNodeClust[ii] = iLastClust;
          piNodeClust[jj] = iLastClust;
        }
        
        else if(piNodeClust[ii] != piNodeClust[jj] && piNodeClust[ii] == 0 && piNodeClust[jj] != 0)
        {
          piNodeClust[ii] = piNodeClust[jj];
        }

        else if(piNodeClust[ii] != piNodeClust[jj] && piNodeClust[ii] != 0 && piNodeClust[jj] == 0)
        {
          piNodeClust[jj] = piNodeClust[ii];
        }
        
        else if(piNodeClust[ii] != piNodeClust[jj] && piNodeClust[ii] != 0 && piNodeClust[jj] != 0)
        {
          iOldClust = piNodeClust[jj];
          for(kk=0; kk<iNumOfRes; ++kk)
          {
            if(piNodeClust[kk] == iOldClust)
              piNodeClust[kk] = piNodeClust[ii];
          }
        }
      }
    }
  }
  
  if(piNodeClust[iRes1] == piNodeClust[iRes2] && piNodeClust[iRes1] != 0)
    iSameClustFlag = 1;
  else
    iSameClustFlag = 0;

  free(piNodeClust);
  iCalcStep = iLastCalcStep;
  
  return iSameClustFlag;
}

void ClearIntMatrix()
{
  int     ii, jj;                                                       // Some iterators
  
  iCalcStep = 13;
  
  if(iCheckClustFlag == 1)
  {
    for(ii=0; ii<iNumOfRes; ii++)
    {
      piClusters[ii] = 0;
      for(jj=0; jj<iNumOfRes; jj++)
      {
        ppfIntMatrix[ii][jj]         = -1.0;
        if(iPSNTypeFlag == 2 && fIntMinFreq == 0.0)
          ppfIntMatrixBackup[ii][jj] = -1.0;
      }
    }
  }
  
  else
  {
    for(ii=0; ii<iNumOfRes; ii++)
      for(jj=0; jj<iNumOfRes; jj++)
        ppfIntMatrix[ii][jj] = -1.0;
  }
  
  iCalcStep = 0;
}

void ClearAll()
{
  // Clears some matrices and vectors
  int     ii, jj;                                                       // Some iterators
  int     iStopFlag = 0;
  
  iCalcStep = 5;
  
  if(iPSNTypeFlag == 0 || iPSNTypeFlag == 2)
  {
    for(ii=0; ii<iNumOfRes; ii++)
    {
      for(jj=0; jj<iNumOfRes; jj++)
      {
        ppfIntMatrix[ii][jj]         = -1.0;
        if(iPSNTypeFlag == 2 && fIntMinFreq == 0.0)
          ppfIntMatrixBackup[ii][jj] = -1.0;
        
        ppiPrevNode[ii][jj]          =    0;
      }
    }
  }
  
  if(iPSNTypeFlag == 1)
  {
    for(ii=0; ii<iNumOfRes; ii++)
    {
      for(jj=0; jj<iNumOfRes; jj++)
      {
        ppiPrevNode[ii][jj]  =    0;
      }
    }
  }
  
  for(ii=0; ii<iNumOfRes; ii++)
  {
    if(pfPathLength[ii] != 999999.9999)
      pfPathLength[ii]   = 999999.9999;
      
    else
      break;
  }
  
  for(ii=0; ii<iNumOfRes; ii++)
  {
    if(piVisited[ii] != 0)
      piVisited[ii]   = 0;
    
    else
      break;
  }

  iStopFlag = 0;
  
  for(ii=0; ii<MAXNUMOFPATHS; ii++)
  {
    
    if(iStopFlag == 1)
      break;
    
    for(jj=0; jj<iNumOfRes; jj++)
    {

      if(ppiPathMatrix[ii][jj] != -1)
        ppiPathMatrix[ii][jj]   = -1;
        
      else
      {
        if(jj == 0)
          iStopFlag = 1;
          
        break;
      }
    }
  }
  
  iNumOfGenPaths  = 0;
  
  iCalcStep = 0;
}

void ClearAllLite()
{
  // Clears some matrices and vectors
  int     ii, jj;                                                       // Some iterators
  int     iStopFlag = 0;
  
  iCalcStep = 14;
  
  for(ii=0; ii<iNumOfRes; ii++)
    for(jj=0; jj<iNumOfRes; jj++)
      ppiPrevNode[ii][jj] = 0;
    
  for(ii=0; ii<iNumOfRes; ii++)
  {
    if(pfPathLength[ii] != 999999.9999)
      pfPathLength[ii]   = 999999.9999;
      
    else
      break;
  }

  for(ii=0; ii<iNumOfRes; ii++)
  {
    if(piVisited[ii] != 0)
      piVisited[ii]   = 0;
    
    else
      break;
  }
  
  iStopFlag = 0;
  
  for(ii=0; ii<MAXNUMOFPATHS; ii++)
  {
    
    if(iStopFlag == 1)
      break;
    
    for(jj=0; jj<iNumOfRes; jj++)
    {

      if(ppiPathMatrix[ii][jj] != -1)
        ppiPathMatrix[ii][jj]   = -1;
        
      else
      {
        if(jj == 0)
          iStopFlag = 1;
          
        break;
      }
    }
  }
  
  iNumOfGenPaths = 0;
  
  iCalcStep = 0;
}

void FindPath()
{
  // Finds shortest path(s)
  int     ii, jj, kk;                                                   // Some iterators
  int     iTmp;
  
  float  fTmpSum;

  iCalcStep = 9;

  pfPathLength[iRes2]=0.0;
  
  iLastDepth=0;
    
  for(kk=0; kk<iNumOfRes; kk++)
  {
    iTmp=-1;
    for(jj=0; jj<iNumOfRes; jj++)
    {
      if(!piVisited[jj] && ((iTmp == -1) || (pfPathLength[jj] < pfPathLength[iTmp])))
      {
        iTmp=jj;
      }
    }
    
    piVisited[iTmp]=1;
    for(ii=0; ii<iNumOfRes; ii++)
    {
      if(ppfIntMatrix[iTmp][ii] != -1.0)
      {
        fTmpSum = floorf((pfPathLength[iTmp] + ppfIntMatrix[iTmp][ii]) * 1000 +  0.5) / 1000;
        //if(pfPathLength[iTmp]+ppfIntMatrix[iTmp][ii] < pfPathLength[ii])
        if(fTmpSum < pfPathLength[ii])
        {
          //pfPathLength[ii]   = pfPathLength[iTmp] + ppfIntMatrix[iTmp][ii];
          pfPathLength[ii]   = fTmpSum;
          ppiPrevNode[ii][0] = 1;
          ppiPrevNode[ii][1] = iTmp;
        }
        //else if(pfPathLength[iTmp] + ppfIntMatrix[iTmp][ii] == pfPathLength[ii])
        else if(fTmpSum == pfPathLength[ii])
        {
          ppiPrevNode[ii][0]++;
          ppiPrevNode[ii][ppiPrevNode[ii][0]] = iTmp;
        }
      }
    }
  }
}

void GenPath(int iDest, int iDepth)
{
  // Generates paths
  int     ii;                                                           // Some iterators
  
  iCalcStep = 10;


  if(iLastDepth >= MAXNUMOFPATHS || iDepth >= MAXNUMOFPATHS)
    PSNPathError(-1);

  if(iDest == iRes1)
    iNumOfGenPaths++;

  if(ppiPathMatrix[iLastDepth][iDepth] == -1)
    ppiPathMatrix[iLastDepth][iDepth] = iDest;
  
  else
  {
    iLastDepth++;
    for(ii=0; ii<iDepth; ii++)
      ppiPathMatrix[iLastDepth][ii] = ppiPathMatrix[iLastDepth-1][ii];
      
    ppiPathMatrix[iLastDepth][iDepth] = iDest;
  }

  for(ii=1; ii<=ppiPrevNode[iDest][0]; ii++)
    GenPath(ppiPrevNode[iDest][ii], iDepth+1);

}

void SavePath()
{
  // Save last path(s)
  int     ii, jj, kk;                                                   // Some iterators
  int     iFoundFlag=0;
  int     iPathNum;
  int     iLastRes=-1;
  
  float   fIntSum = 0.0;
  float   fTmpCorrMin1=FLT_MAX, fTmpCorrMax1=FLT_MIN, fTmpCorrAvg1=0.0; // Min, Max and Avg Res 1 Corr Values
  float   fTmpCorrMin2=FLT_MAX, fTmpCorrMax2=FLT_MIN, fTmpCorrAvg2=0.0; // Min, Max and Avg Res 2 Corr Values
  float   fScoreMin=FLT_MAX, fScoreMax=FLT_MIN, fAvgScore=0.0;          // Used by -plt option
  
  iCalcStep = 11;
  
  for(ii=0; ii<MAXNUMOFPATHS; ii++)
  {
    pcTmpPaths[ii][0] = '\0';
    piTmpPathsLen[ii] = 0;
    pfTmpScore[ii]    = 0.0;
    pfTmpWeights[ii]  = 0.0;
  }
  
  iPathNum = -1;

  if(iNumOfGenPaths>0)
  {
    for(ii=0; ii<MAXNUMOFPATHS; ii++)
    {
      if(ppiPathMatrix[ii][0]!=-1)
      {
        iPathNum++;
        strcpy(pcTmpPaths[iPathNum], "");
        piTmpPathsLen[iPathNum]=0;
        pfTmpScore[iPathNum]=0.0;
        
        fTmpCorrMin1 = FLT_MAX;
        fTmpCorrMax1 = FLT_MIN;
        fTmpCorrAvg1 = 0.0;
        fTmpCorrMin2 = FLT_MAX;
        fTmpCorrMax2 = FLT_MIN;
        fTmpCorrAvg2 = 0.0;
        iLastRes     = -1;
        fIntSum      = 0.0;
        
        for(jj=0; jj<iNumOfRes; jj++)
        {
          if(ppiPathMatrix[ii][jj]!=-1)
          {
            
            if(iLastRes != -1)
              fIntSum  = fIntSum + ppfRealIntMatrix[iLastRes][ppiPathMatrix[ii][jj]];
            
            iLastRes = ppiPathMatrix[ii][jj];
            
            piTmpPathsLen[iPathNum]++;
            pfTmpScore[iPathNum]+=piCorrVector[ppiPathMatrix[ii][jj]];

            strcat(pcTmpPaths[iPathNum], pcSequence[ppiPathMatrix[ii][jj]]);
            strcat(pcTmpPaths[iPathNum], "=>");
            
            if(ppfCorMatrix[iRes1][ppiPathMatrix[ii][jj]] >= fCorrCutOff && ppiPathMatrix[ii][jj] != iRes1 && ppiPathMatrix[ii][jj] != iRes2)
            {
              if(fTmpCorrMin1 > ppfCorMatrix[iRes1][ppiPathMatrix[ii][jj]])
                fTmpCorrMin1 = ppfCorMatrix[iRes1][ppiPathMatrix[ii][jj]];
              
              if(fTmpCorrMax1 < ppfCorMatrix[iRes1][ppiPathMatrix[ii][jj]])
                fTmpCorrMax1 = ppfCorMatrix[iRes1][ppiPathMatrix[ii][jj]];
              
              fTmpCorrAvg1 = fTmpCorrAvg1 + ppfCorMatrix[iRes1][ppiPathMatrix[ii][jj]];
            }
            
            if(ppfCorMatrix[iRes2][ppiPathMatrix[ii][jj]] >= fCorrCutOff && ppiPathMatrix[ii][jj] != iRes1 && ppiPathMatrix[ii][jj] != iRes2)
            {
              if(fTmpCorrMin2 > ppfCorMatrix[iRes2][ppiPathMatrix[ii][jj]])
                fTmpCorrMin2 = ppfCorMatrix[iRes2][ppiPathMatrix[ii][jj]];
                
              if(fTmpCorrMax2 < ppfCorMatrix[iRes2][ppiPathMatrix[ii][jj]])
                fTmpCorrMax2 = ppfCorMatrix[iRes2][ppiPathMatrix[ii][jj]];
                
              fTmpCorrAvg2 = fTmpCorrAvg2 + ppfCorMatrix[iRes2][ppiPathMatrix[ii][jj]];
            }
          }
          else
            break;
        }
        
        pcTmpPaths[ii][strlen(pcTmpPaths[iPathNum])-2]='\0';
        pfTmpWeights[ii] = fIntSum;
        
        if(fTmpCorrMin1 == FLT_MAX)
          fTmpCorrMin1 = 0.0;
        if(fTmpCorrMax1 == FLT_MIN)
          fTmpCorrMax1 = 0.0;
          
        if(fTmpCorrMin2 == FLT_MAX)
          fTmpCorrMin2 = 0.0;
        if(fTmpCorrMax2 == FLT_MIN)
          fTmpCorrMax2 = 0.0;
          
        ppfTmpCorrStat[iPathNum][0] = fTmpCorrMin1;
        ppfTmpCorrStat[iPathNum][1] = fTmpCorrMax1;
        ppfTmpCorrStat[iPathNum][2] = (fTmpCorrAvg1 / (pfTmpScore[iPathNum] - 2));
        ppfTmpCorrStat[iPathNum][3] = fTmpCorrMin2;
        ppfTmpCorrStat[iPathNum][4] = fTmpCorrMax2;
        ppfTmpCorrStat[iPathNum][5] = (fTmpCorrAvg2 / (pfTmpScore[iPathNum] - 2));
      }
      else
        break;
    }
  }
  
  else
  {
    iPathNum=0;
    fIntSum      = 0.0;
    strcpy(pcTmpPaths[0], "[NULL_PATH]");
    piTmpPathsLen[0]=0;
    pfTmpScore[0]=0.0;
  }
    
  iPathNum++;
  
  for(kk=0; kk<iPathNum; kk++)
  {
    if(pfTmpScore[kk] == 2.0 )
      continue;

    if(piTmpPathsLen[kk] == 1)
    {
      // Null path :(
      strcpy(pcTmpPaths[kk], "[NULL_PATH]");
      piTmpPathsLen[kk]=0;
      pfTmpScore[kk]=0.0;
    }
    
    // Discards paths smaller than user-defined minimum length
    if(piTmpPathsLen[kk] != 0 && (piTmpPathsLen[kk] - 2) < iMinLen)
    {
      continue;
    }
    
    if(iStatFlag == 1)
    {
      fAvgScore=fAvgScore+pfTmpScore[kk];
      
      if(pfTmpScore[kk]<fScoreMin)
        fScoreMin = pfTmpScore[kk];
        
      if(pfTmpScore[kk]>fScoreMax)
        fScoreMax = pfTmpScore[kk];
    }
    
    if(iFrameFlag == 1)
      fprintf(FFrameFile, "%9d   %s\n", iFrameCont, pcTmpPaths[kk]);
    
    if(iNumOfPaths==0)
    {
      iNumOfPaths=1;
      strcpy(pcPathList[0], "");
      strcpy(pcPathList[0], pcTmpPaths[kk]);
      piPathFreq[0]++;
      piNodesInPath[0]      = piTmpPathsLen[kk];
      pfPathScores[0]       = pfTmpScore[kk];
      ppfCorrStatList[0][0] = ppfTmpCorrStat[kk][0];
      ppfCorrStatList[0][1] = ppfTmpCorrStat[kk][1];
      ppfCorrStatList[0][2] = ppfTmpCorrStat[kk][2];
      ppfCorrStatList[0][3] = ppfTmpCorrStat[kk][3];
      ppfCorrStatList[0][4] = ppfTmpCorrStat[kk][4];
      ppfCorrStatList[0][5] = ppfTmpCorrStat[kk][5];
      pfPathWeights[0]      = pfTmpWeights[kk];
    }

    else
    {
      iFoundFlag=0;
      for(ii=0; ii<iNumOfPaths; ii++)
      {
        if(strlen(pcTmpPaths[kk])==strlen(pcPathList[ii]))
        {
          if(strcmp(pcTmpPaths[kk], pcPathList[ii])==0)
          {
            iFoundFlag=1;
            piPathFreq[ii]++;
            pfPathWeights[ii] = pfPathWeights[ii] + pfTmpWeights[kk];
            break;
          }
        }
      }
      
      if(iFoundFlag==0)
      {
        iNumOfPaths++;
        
        if(iNumOfPaths >= GAIN*MAXNUMOFPATHS)
        {
          fprintf( stderr, "PSNPATH module: Not enough memory\nTry to increase --MEMPATH and/or--MEMGAIN values\n");
          fprintf( stderr, "Actual values :\n--MEMPATH : %d\n--MEMGAIN : %d\n", MAXNUMOFPATHS, GAIN);
          fprintf( stderr, "Actual pair   :\nPair : %d\nRes1 : %d\nRes2 : %d\n", iPairCont, iRes1, iRes2);
          return;
        }
        
        pcPathList[iNumOfPaths-1][0]='\0';
        strcpy(pcPathList[iNumOfPaths-1], pcTmpPaths[kk]);
        piPathFreq[iNumOfPaths-1]++;
        piNodesInPath[iNumOfPaths-1]      = piTmpPathsLen[kk];
        pfPathScores[iNumOfPaths-1]       = pfTmpScore[kk];
        ppfCorrStatList[iNumOfPaths-1][0] = ppfTmpCorrStat[kk][0];
        ppfCorrStatList[iNumOfPaths-1][1] = ppfTmpCorrStat[kk][1];
        ppfCorrStatList[iNumOfPaths-1][2] = ppfTmpCorrStat[kk][2];
        ppfCorrStatList[iNumOfPaths-1][3] = ppfTmpCorrStat[kk][3];
        ppfCorrStatList[iNumOfPaths-1][4] = ppfTmpCorrStat[kk][4];
        ppfCorrStatList[iNumOfPaths-1][5] = ppfTmpCorrStat[kk][5];
        pfPathWeights[iNumOfPaths-1]      = pfTmpWeights[kk];
      }
    }
  }
  
  if(iStatFlag == 1)
    fprintf(FStatFile, "%9d   %9d   %9d   %8.2f   %8.2f   %8.2f   %8.2f\n",
            iFrameCont, iPathNum, piTmpPathsLen[kk-1]-2, (fAvgScore/(float)iPathNum), fScoreMin, fScoreMax, (fScoreMax-fScoreMin));
}

void SavePathLite()
{
  // Save last path(s)
  int     ii, jj;                                                       // Some iterators
  int     iCorrFlag;
  int     iPathLen;
  int     iLastRes;
  
  float   fIntSum = 0.0;
  
  iCalcStep = 16;
  if(iNumOfGenPaths>0)
  {
    for(ii=0; ii<MAXNUMOFPATHS; ii++)
    {
      if(ppiPathMatrix[ii][0]!=-1)
      {
        strcpy(pcPath, "");
        iPathLen  = 0;
        iCorrFlag = 0;
        iLastRes  = -1;
        fIntSum   = 0.0;
        
        for(jj=0; jj<iNumOfRes; jj++)
        {
          if(ppiPathMatrix[ii][jj]!=-1)
          {
            iPathLen++;
            
            if(iLastRes != -1)
              fIntSum  = fIntSum + ppfRealIntMatrix[iLastRes][ppiPathMatrix[ii][jj]];
            
            iLastRes = ppiPathMatrix[ii][jj];

            strcat(pcPath, pcSequence[ppiPathMatrix[ii][jj]]);
            strcat(pcPath, "-");
            
            if(ppfCorMatrix[iRes1][ppiPathMatrix[ii][jj]] >= fCorrCutOff && iRes1 != ppiPathMatrix[ii][jj] && iRes2 != ppiPathMatrix[ii][jj])
              iCorrFlag = 1;
            
            if(ppfCorMatrix[iRes2][ppiPathMatrix[ii][jj]] >= fCorrCutOff && iRes1 != ppiPathMatrix[ii][jj] && iRes2 != ppiPathMatrix[ii][jj])
              iCorrFlag = 1;
          }
          
          else
            break;
        }
        
        pcPath[strlen(pcPath)-1]='\0';
        
        if(iPathLen < 3)
          continue;
        
        if(iCorrFlag == 0)
          continue;
        
        if(iPathLen - 2 < iMinLen)
          continue;
          
        iNumOfPaths++;
        
        fprintf(FFrameFile, "%s,%.2f\n", pcPath, fIntSum);
      }
      
      else
        break;
    }
  }
  
  iCalcStep = 0;
}

void SaveData()
{
  // Save Paths on file
  int     ii;                                                           // Some iterators
  int     iShortPathNum=0, iLongPathNum=0;
  int     iHighFreqNum=0,  iLowFreqNum=0;
  int     iHighScoreNum=0, iLowScoreNum=0;
  int     iShortPath=INT_MAX, iLongPath=0;
  int     iNullPathFlag=0, iNullPathIndex=0;
  int     iNumOfR1CorrRes, iNumOfR2CorrRes;
  int     iPathNum=0, iFinalNumOfPaths=0;
  int     iRealHighFreqNum=0;
  
  float   fHighFreq=FLT_MIN, fLowFreq=FLT_MAX;
  float   fHighScore=FLT_MIN, fLowScore=FLT_MAX;
  float   fHighScoreW=FLT_MIN, fLowScoreW=FLT_MAX;
  float   fNullPathFreq=0.0;
  
  iCalcStep = 12;
  //iFrameCont++;
  
  for(ii=0; ii<iNumOfPaths; ii++)
  {
    // Discards paths with frequencies under user-defined limit
    if(((piPathFreq[ii]*100)/(float)(iFrameCont)) < fMinFreq)
      continue;

    iFinalNumOfPaths++;
    
    if(piNodesInPath[ii]==0)
    {
      iNullPathFlag=1;
      iNullPathIndex=ii+1;
      fNullPathFreq=(float)piPathFreq[ii];
    }

    if(piNodesInPath[ii]<iShortPath && piNodesInPath[ii]!=0)
    {
      iShortPath=piNodesInPath[ii];
      iShortPathNum=iFinalNumOfPaths;
    }
    
    if(piNodesInPath[ii]>iLongPath && piNodesInPath[ii]!=0)
    {
      iLongPath=piNodesInPath[ii];
      iLongPathNum=iFinalNumOfPaths;
    }
    
    if(piPathFreq[ii]<fLowFreq && piNodesInPath[ii]!=0)
    {
      fLowFreq=(float)piPathFreq[ii];
      iLowFreqNum=iFinalNumOfPaths;
    }
    
    if(piPathFreq[ii]>fHighFreq && piNodesInPath[ii]!=0)
    {
      fHighFreq=(float)piPathFreq[ii];
      iHighFreqNum=iFinalNumOfPaths;
      iRealHighFreqNum=ii;
    }
    
    if((pfPathScores[ii]/(float)piNodesInPath[ii])<fLowScoreW && piNodesInPath[ii]!=0)
    {
      fLowScore=pfPathScores[ii]-2;
      fLowScoreW=((pfPathScores[ii]-2)/(float)(piNodesInPath[ii]-2));
      iLowScoreNum=iFinalNumOfPaths;
    }
    
    if((pfPathScores[ii]/(float)piNodesInPath[ii])>fHighScoreW && piNodesInPath[ii]!=0)
    {
      fHighScore=pfPathScores[ii]-2;
      fHighScoreW=((pfPathScores[ii]-2)/(float)(piNodesInPath[ii]-2));
      iHighScoreNum=iFinalNumOfPaths;
    }
  }
  
  if(fHighFreq==FLT_MIN)
    fHighFreq=0.0;
    
  if(fLowFreq==FLT_MAX)
    fLowFreq=0.0;
    
  if(fHighScore==FLT_MIN)
    fHighScore=0.0;
    
  if(fLowScore==FLT_MAX)
    fLowScore=0.0;
    
  if(fHighScoreW==FLT_MIN)
    fHighScoreW=0.0;
    
  if(fLowScoreW==FLT_MAX)
    fLowScoreW=0.0;
    
  if(iShortPath==INT_MAX)
    iShortPath=0;

  time(&time_Finish);

  if(iLogFlag == 1)
  {
    if(((fHighFreq*100)/(float)(iFrameCont)) != 0.0)
    {
      fprintf(FLogFile, "%7d   %7s   %7s      %6.2f         %6.2f    %7d        %6.2f        %6d        %6d   %s\n",
        iPairCont, pcSequence[iRes1], pcSequence[iRes2],
        (float)((iNumOfBadFrames*100)/(float)(iFrameCont)),
        ((fNullPathFreq*100)/(float)(iFrameCont)),
        iFinalNumOfPaths, ((fHighFreq*100)/(float)(iFrameCont)),
        piNodesInPath[iRealHighFreqNum]-2,
        (int)pfPathScores[iRealHighFreqNum]-2,
        pcPathList[iRealHighFreqNum]);
    }
    
    else
    {
      fprintf(FLogFile, "%7d   %7s   %7s      %6.2f         %6.2f    %7d        %6.2f        %6d        %6d   %s\n",
        iPairCont, pcSequence[iRes1], pcSequence[iRes2],
        (float)((iNumOfBadFrames*100)/(float)(iFrameCont)),
        ((fNullPathFreq*100)/(float)(iFrameCont)),
        iFinalNumOfPaths, ((fHighFreq*100)/(float)(iFrameCont)),
        0, 0, pcPathList[iRealHighFreqNum]);
    }
  }
  
  fprintf(FOutFile, "# =========================================\n");
  fprintf(FOutFile, "#       *** WORDOM PSNPATH MODULE ***      \n");
  fprintf(FOutFile, "# =========================================\n");
  fprintf(FOutFile, "#\n");
  fprintf(FOutFile, "# Version        : %s\n",       PSNPATHVERSTRING);
  fprintf(FOutFile, "# License        : GPL 3\n");
  fprintf(FOutFile, "# Copyright      : Fanelli, Felline\n");
  fprintf(FOutFile, "#                  University of Modena\n");
  fprintf(FOutFile, "#                  Modena - Italy\n");
  fprintf(FOutFile, "#\n");
  fprintf(FOutFile, "# Date           : %s",         asctime(localtime(&time_Finish)));
  fprintf(FOutFile, "#\n");
  fprintf(FOutFile, "# Calculation T. : ~ %d sec\n", (int)(difftime(time_Finish, time_Start)));
  fprintf(FOutFile, "#\n");
  
  if(iPSNTypeFlag == 0)
    fprintf(FOutFile, "# Type           : RAW\n");
  else if(iPSNTypeFlag == 1)
    fprintf(FOutFile, "# Type           : AVG\n");
  else if(iPSNTypeFlag == 2)
    fprintf(FOutFile, "# Type           : MIX\n");
  
  if(iPSNTypeFlag == 0)
    fprintf(FOutFile, "# RAW PSN File   : %s\n",       cRawPSNFileName);
  else
    fprintf(FOutFile, "# RAW PSN File   : none\n");
  
  if(iPSNTypeFlag != 0)
    fprintf(FOutFile, "# AVG PSN File   : %s\n",       cAvgPSNFileName);
  else
    fprintf(FOutFile, "# AVG PSN File   : none\n");
  
  fprintf(FOutFile, "# First Res      : %s\n",       pcSequence[iRes1]);
  fprintf(FOutFile, "# Last  Res      : %s\n",       pcSequence[iRes2]);
  fprintf(FOutFile, "# IntMin CutOff  : %.3f\n",     fIntMinCutOff);
  
  if(iPSNTypeFlag != 0)
    fprintf(FOutFile, "# Int Min Freq   : %.3f\n",     fIntMinFreq);
  else
    fprintf(FOutFile, "# Int Min Freq   : none\n");
  
  fprintf(FOutFile, "# Corr File      : %s\n",       cCorrFileName);
  fprintf(FOutFile, "# Corr CutOff    : %.2f\n",     fCorrCutOff);
  fprintf(FOutFile, "# Min Len        : %d\n",       iMinLen);
  fprintf(FOutFile, "# Min Freq       : %.2f\n",     fMinFreq);
  fprintf(FOutFile, "# Seg OffSet     : %s %d\n",    pcSegNameList[0], piSegOffSet[0]);
  for(ii=1; ii<iNumOfSeg-1; ii++)
    fprintf(FOutFile, "#                : %s %d\n",  pcSegNameList[ii], piSegOffSet[ii]);
  
  if(iStatFlag == 0)
    fprintf(FOutFile, "# Stat file      : none\n");
  else
    fprintf(FOutFile, "# Stat file      : %s\n",     cStatFileName);
  
  if(iFrameFlag == 0)
    fprintf(FOutFile, "# Frame file     : none\n");
  else
    fprintf(FOutFile, "# Frame file     : %s\n",     cFrameFileName);
  
  fprintf(FOutFile, "#\n");
  fprintf(FOutFile, "# Total Res      : %d\n",       iNumOfRes);
  fprintf(FOutFile, "# Total Frames   : %d\n",       iNumOfFrames);
  fprintf(FOutFile, "#\n");
  
  fprintf(FOutFile, "# Tot Bad Frames : %.2f%% (%d/%d)\n", (float)((iNumOfBadFrames*100)/(float)(iFrameCont)), iNumOfBadFrames, iFrameCont);

  iNumOfR1CorrRes = 0;
  iNumOfR2CorrRes = 0;
    
  fprintf(FOutFile, "# Res 1 Corr Res : ");
  for(ii=0; ii<iNumOfRes; ii++)
  {
    if(ii == iRes1)
      continue;
    
    if(ppfCorMatrix[iRes1][ii] >= fCorrCutOff)
    {
      iNumOfR1CorrRes++;
      fprintf(FOutFile, "%s %.2f ", pcSequence[ii], ppfCorMatrix[iRes1][ii]);
    }
  }
  fprintf(FOutFile, "\n");
  
  fprintf(FOutFile, "# Res 2 Corr Res : ");
  for(ii=0; ii<iNumOfRes; ii++)
  {
    if(ii == iRes2)
      continue;
    
    if(ppfCorMatrix[iRes2][ii] >= fCorrCutOff)
    {
      iNumOfR2CorrRes++;
      fprintf(FOutFile, "%s %.2f ", pcSequence[ii], ppfCorMatrix[iRes2][ii]);
    }
  }
  fprintf(FOutFile, "\n");
  
  fprintf(FOutFile, "# Total Corr Res : %d (Res 1 = %d; Res 2 = %d)\n", iNumOfCorrRes-2, iNumOfR1CorrRes, iNumOfR2CorrRes);
  
  
  fprintf(FOutFile, "# Total Paths    : %d\n", iFinalNumOfPaths);
  fprintf(FOutFile, "#\n");
  fprintf(FOutFile, "# Shortest  Path : %d (path %d)\n", iShortPath-2, iShortPathNum);
  fprintf(FOutFile, "# Longest   Path : %d (path %d)\n", iLongPath-2, iLongPathNum);
  fprintf(FOutFile, "#\n");
  fprintf(FOutFile, "# Highest   Freq : %.2f %% (path %d)\n", ((fHighFreq*100)/(float)(iFrameCont)),iHighFreqNum);
  fprintf(FOutFile, "# Lowest    Freq : %.2f %% (path %d)\n", ((fLowFreq*100)/(float)(iFrameCont)), iLowFreqNum);
  fprintf(FOutFile, "#\n");
  fprintf(FOutFile, "# Highest  Score : %.2f (%.2f) (path %d)\n", fHighScore, fHighScoreW, iHighScoreNum);
  fprintf(FOutFile, "# Lowest   Score : %.2f (%.2f) (path %d)\n", fLowScore, fLowScoreW, iLowScoreNum);
  
  if(iNullPathFlag == 1)
    fprintf(FOutFile, "#\n# Null Path Freq : %.2f %% (path %d)\n", ((fNullPathFreq*100)/(float)(iFrameCont)), iNullPathIndex);
  
  fprintf(FOutFile, "# =========================================\n\n");
  
  for(ii=0; ii<iNumOfPaths; ii++)
  {
    // Discards paths with frequencies under user-defined limit
    if(((piPathFreq[ii]*100)/(float)(iFrameCont)) < fMinFreq)
      continue;
    
    iPathNum++;
    fprintf(FOutFile, "path # : %d\n", iPathNum);
    fprintf(FOutFile, "links  : %s\n", pcPathList[ii]);
    
    if(piNodesInPath[ii] != 0)
      fprintf(FOutFile, "length : %d\n", piNodesInPath[ii]-2);
    
    else
      fprintf(FOutFile, "length : 0\n");
    
    fprintf(FOutFile, "freq   : %5.2f%% (%d/%d)\n", ((piPathFreq[ii]*100)/(float)(iFrameCont)), piPathFreq[ii], iFrameCont);
    
    if(pfPathScores[ii] != 0.0)
    {
      fprintf(FOutFile, "corr   : Num %d; N/L %.2f; R1Min %.2f; R1Max %.2f; R1Avg %.2f; R2Min %.2f; R2Max %.2f; R2Avg %.2f\n",
                                  (int) pfPathScores[ii]-2, ((pfPathScores[ii]-2)/(float)(piNodesInPath[ii]-2)),
                                  ppfCorrStatList[ii][0], ppfCorrStatList[ii][1], ppfCorrStatList[ii][2],
                                  ppfCorrStatList[ii][3], ppfCorrStatList[ii][4], ppfCorrStatList[ii][5]);
    }

    else
    {
      fprintf(FOutFile, "corr   : Num 0; N/L 0.0; R1Min 0.0; R1Max 0.0; R1Avg 0.0; R2Min 0.0; R2Max 0.0; R2Avg 0.0\n");
    }
    
    fprintf(FOutFile, "weight : %.2f\n\n", (pfPathWeights[ii] / (float) piPathFreq[ii]));
  }
  
  iCalcStep = 0;
}

// === PSN Param Section ===============================================
int InitPSNParam(char **ppcInput, int iInputLineNum)
{
  int     ii;                                                           // Some iterators
  int     iTmpCont1, iTmpCont2;
  int     iOriginalInputLineNum;
  int     iOptionFlag;                                                  // Used to catch invalid options
  
  char    pcInputLine[1024];                                            // Input lines
  char    cFileName[1024], cSeleString[1024];                           // Temporary strings
  char   *cTemp1;
  
  float   fParam;                                                       // Normalization factor
  
  struct inp_psnparam inp_psnparam;
  
  // === Default values ==================
  inp_psnparam.iAvgMode          = 0;
  inp_psnparam.iNumOfIgnore        = 0;
  inp_psnparam.iNumOfMol           = 0;
  inp_psnparam.iNumOfTargetAtoms   = 0;
  inp_psnparam.iNumOfWarning       = 0;
  inp_psnparam.iVerboseFlag        = 0;
  inp_psnparam.fDistCutoff         = 4.5;
  inp_psnparam.fAvgMaxInteractions = 0.0;
  inp_psnparam.cTarget[0]          = '\0';
  inp_psnparam.cVerboseFileName[0] = '\0';
  strcpy(inp_psnparam.cTitle, "NoTitle");
  // =====================================

  memset(pcInputLine, '\0', sizeof(pcInputLine));
  iOriginalInputLineNum = iInputLineNum;
  
  // Process input file directives
  while(strncmp(pcInputLine, "END", 3) != 0)
  {
    iOptionFlag = 0;
    sprintf(pcInputLine, "%s", ppcInput[iInputLineNum]);
    
    if(!strncmp(pcInputLine, "BEGIN", 5) || !strncmp(pcInputLine, "END", 3) || pcInputLine[0] == '#')
      iOptionFlag = 1;
    
    else if(strncmp(pcInputLine, "--TITLE", 7) == 0)
    {
      sscanf(pcInputLine, "--TITLE %s", inp_psnparam.cTitle);
      iOptionFlag = 1;
    }

    else if(strncmp(pcInputLine, "--TARGET", 8) == 0)
    {
      sscanf(pcInputLine, "--TARGET %s", inp_psnparam.cTarget);
      iOptionFlag = 1;
    }

    else if(strncmp(pcInputLine, "--MOL", 5) == 0)
    {
      inp_psnparam.iNumOfMol++;
      iOptionFlag = 1;
    }
    
    else if(strncmp(pcInputLine, "--IGNORE", 8) == 0)
    {
      inp_psnparam.iNumOfIgnore++;
      iOptionFlag = 1;
    }

    else if (strncmp(pcInputLine, "--VERBOSE", 9) == 0)
    {
      sscanf(pcInputLine, "--VERBOSE %d", &inp_psnparam.iVerboseFlag);
      if(inp_psnparam.iVerboseFlag != 0 && inp_psnparam.iVerboseFlag != 1)
      {
        fprintf(stderr, "PSNParam module: Invalid --VERBOSE value: %d, valid values are: 0 and 1\n", inp_psnparam.iVerboseFlag);
        return 1;
      }
      iOptionFlag = 1;
    }

    else if (strncmp(pcInputLine, "--AVGMODE", 9) == 0)
    {
      sscanf(pcInputLine, "--AVGMODE %d", &inp_psnparam.iAvgMode);
      if(inp_psnparam.iAvgMode != 0 && inp_psnparam.iAvgMode != 1)
      {
        fprintf(stderr, "PSNParam module: Invalid --AVGMODE value: %d, valid values are: 0 and 1\n", inp_psnparam.iVerboseFlag);
        return 1;
      }
      iOptionFlag = 1;
    }

    if(iOptionFlag == 0)
    {
      fprintf(stderr, "PSNParam module: Could NOT understand option: %s\n", pcInputLine);
      return 1;
    }
    
    iInputLineNum++;
  }

  // some sanity checks
  if(strlen(inp_psnparam.cTarget) == 0)
  {
    fprintf(stderr, "PSNParam module: Need a --TARGET option\n");
    return 1;
  }
  
  if(inp_psnparam.iNumOfMol == 0)
  {
    fprintf(stderr, "PSNParam module: Need at least one --MOL option\n");
    return 1;
  }
  
  // allocates some vectors
  inp_psnparam.ccMolFileVect = (char **) calloc(inp_psnparam.iNumOfMol, sizeof(char *));
  inp_psnparam.ccMolSeleVect = (char **) calloc(inp_psnparam.iNumOfMol, sizeof(char *));
  for(ii=0; ii<inp_psnparam.iNumOfMol; ++ii)
  {
    inp_psnparam.ccMolFileVect[ii] = (char *) calloc(1024, sizeof(char));
    inp_psnparam.ccMolSeleVect[ii] = (char *) calloc(1024, sizeof(char));
  }
  
  inp_psnparam.ccIgnoreVect = (char **) calloc(inp_psnparam.iNumOfIgnore, sizeof(char *));
  for(ii=0; ii<inp_psnparam.iNumOfIgnore; ++ii)
    inp_psnparam.ccIgnoreVect[ii] = (char *) calloc(50, sizeof(char));
  
  memset(pcInputLine, '\0', sizeof(pcInputLine));
  iInputLineNum = iOriginalInputLineNum;
  
  iTmpCont1 = -1;
  iTmpCont2 = -1;
  // Process input file directives
  while(strncmp(pcInputLine, "END", 3) != 0)
  {
    sprintf(pcInputLine, "%s", ppcInput[iInputLineNum]);
    if(strncmp(pcInputLine, "--MOL ", 6) == 0)
    {
      iTmpCont1++;
      
      strcpy(cSeleString, pcInputLine + 6);                             // get option values
      cSeleString[strlen(cSeleString) - 1] = '\0';                      // delete newline char
      
      cTemp1 = strchr(cSeleString, ' ');                                // find the first space char
      
      if(cTemp1)
      {
        strncpy(cFileName, cSeleString, cTemp1 - cSeleString);          // extract file name
        cFileName[cTemp1 - cSeleString] = '\0';                         // add null
        
        strcpy(cSeleString, cSeleString + strlen(cFileName) + 1);       // extract sele string
        
        //printf("*-* |%s| |%s|\n", cFileName, cSeleString);
        
        
        strcpy(inp_psnparam.ccMolFileVect[iTmpCont1], cFileName);       // update file name vect
        strcpy(inp_psnparam.ccMolSeleVect[iTmpCont1], cSeleString);     // updare sele string vect
      }
      
      else
      {
        // something goes wrong
        fprintf(stderr, "PSNParam module: Invalid --MOL line: %s\n", pcInputLine);
        exit(1);
      }
    }
    
    if(strncmp(pcInputLine, "--IGNORE", 8) == 0)
    {
      iTmpCont2++;
      sscanf(pcInputLine, "--IGNORE %s", cSeleString);
      strcpy(inp_psnparam.ccIgnoreVect[iTmpCont2], cSeleString);
    }
    
    iInputLineNum++;
  }

  if(inp_psnparam.iVerboseFlag == 1)
  {
    sprintf(inp_psnparam.cVerboseFileName, "PSNParam-%s.log", inp_psnparam.cTitle);
    inp_psnparam.FVerboseFile = fopen(inp_psnparam.cVerboseFileName, "w");
    
    time(&inp_psnparam.time_tToday);
    fprintf(inp_psnparam.FVerboseFile, "# ==========================================\n");
    fprintf(inp_psnparam.FVerboseFile, "#       *** WORDOM PSNParam MODULE ***      \n");
    fprintf(inp_psnparam.FVerboseFile, "# ==========================================\n");
    fprintf(inp_psnparam.FVerboseFile, "#\n");
    fprintf(inp_psnparam.FVerboseFile, "# Version         : 0.1a\n");
    fprintf(inp_psnparam.FVerboseFile, "# License         : GPL 3\n");
    fprintf(inp_psnparam.FVerboseFile, "# Copyright       : Fanelli, Felline\n");
    fprintf(inp_psnparam.FVerboseFile, "#                   University of Modena\n");
    fprintf(inp_psnparam.FVerboseFile, "#                   Modena - Italy\n");
    fprintf(inp_psnparam.FVerboseFile, "#\n");
    fprintf(inp_psnparam.FVerboseFile, "# Date            : %s", asctime(localtime(&inp_psnparam.time_tToday)));
    fprintf(inp_psnparam.FVerboseFile, "#\n");
    fprintf(inp_psnparam.FVerboseFile, "# Title           : %s\n", inp_psnparam.cTitle);
    fprintf(inp_psnparam.FVerboseFile, "# Target          : %s\n", inp_psnparam.cTarget);
    fprintf(inp_psnparam.FVerboseFile, "# Average Mode    : %d\n", inp_psnparam.iAvgMode);
    fprintf(inp_psnparam.FVerboseFile, "# Distance CutOff : %f\n", inp_psnparam.fDistCutoff);
    fprintf(inp_psnparam.FVerboseFile, "#\n");
    fprintf(inp_psnparam.FVerboseFile, "# Num Of Mol      : %d\n", inp_psnparam.iNumOfMol);
    for(ii=0; ii<inp_psnparam.iNumOfMol; ++ii)
      fprintf(inp_psnparam.FVerboseFile, "# Mol %d, File %s, Sele %s\n", ii+1, inp_psnparam.ccMolFileVect[ii], inp_psnparam.ccMolSeleVect[ii]);
    fprintf(inp_psnparam.FVerboseFile, "#\n");

    fprintf(inp_psnparam.FVerboseFile, "# Num Of Ignore   : %d\n", inp_psnparam.iNumOfIgnore);
    for(ii=0; ii<inp_psnparam.iNumOfIgnore; ++ii)
      fprintf(inp_psnparam.FVerboseFile, "# Ignore %d, Res %s\n", ii+1, inp_psnparam.ccIgnoreVect[ii]);
    fprintf(inp_psnparam.FVerboseFile, "#\n");
    
    fprintf(inp_psnparam.FVerboseFile, "# ==========================================\n\n");
  }
  
  for(ii=0; ii<inp_psnparam.iNumOfMol; ++ii)
    CalcPSNParam(&inp_psnparam, ii);
  
  fParam = inp_psnparam.fAvgMaxInteractions / (float) inp_psnparam.iNumOfMol;
  
  if(inp_psnparam.iNumOfWarning != 0)
  {
    if(inp_psnparam.iNumOfWarning == 1)
      fprintf(stderr, "There was a warning, run again this analysis with --VERBOSE option setted to 1 and read the log file\n");
    else
      fprintf(stderr, "There were %d warnings, run again this analysis with --VERBOSE option setted to 1 and read the log file\n", inp_psnparam.iNumOfWarning);
  }

  fParam = inp_psnparam.fAvgMaxInteractions / (float) inp_psnparam.iNumOfMol;

  fprintf(stdout, "The Normalization Factor for mol %s is %f\n", inp_psnparam.cTarget, fParam);
  
  if(inp_psnparam.iVerboseFlag == 1)
  {
    fprintf(inp_psnparam.FVerboseFile, "\n\nThe Normalization Factor for mol %s is %f\n", inp_psnparam.cTarget, fParam);
    fprintf(inp_psnparam.FVerboseFile, "\nTo use this normalization factor add the following line in your PSN inputs:\n--PARAM %s %f\n", inp_psnparam.cTarget, fParam);
    fclose(inp_psnparam.FVerboseFile);
  }

  return 0;
}

void CalcPSNParam(struct inp_psnparam *inp_psnparam, int iMolIndex)
{
  int     ii, jj;                                                       // some iterators
  int     iNumOfTarget;                                                 // the number of target residue in this molecule
  int     iNumOfSeleRes;                                                // the number of selected residues
  int     iSetTargetAtmNumFlag;                                         // 1 if the number of target atoms has to be setted
  int     iTargetId;                                                    // progressive target id number
  int     iNumOfInteractions;                                           // the number of target interactions
  int     iNumOfThisTargetAtoms;                                        // the number of this target atoms
  int     iNumOfIgnoredInt;                                             // the number of ingored interactions
  
  char    cLastRes[15], cThisRes[15], cTestRes[15];
  char    cBestTarget[100], cThisTarget[100];
  
  float   fTargetCoordX, fTargetCoordY, fTargetCoordZ;                  // target x, y, and z coordinates
  float   fTestCoordX, fTestCoordY, fTestCoordZ;                        // other residues x, y and z coordinates
  float   fDist;                                                        // Atom-Atom distance
  float   fMaxInteractions;                                             // highest number of target itneractions
  
  // Set the number of target atoms from the first target residue in the first passed molecule
  if(iMolIndex == 0)
    iSetTargetAtmNumFlag = 1;
  else
    iSetTargetAtmNumFlag = 0;

  // read the iMolIndex-th molecule
  inp_psnparam->molecule = ReadMolecule(inp_psnparam->ccMolFileVect[iMolIndex], wrd_whichmoltype(inp_psnparam->ccMolFileVect[iMolIndex]));
  
  // PBC are ignored - structures are supposed to be crystallographic data
  if( inp_psnparam->molecule->coor.pbc_flag != 0 )
  {
    inp_psnparam->molecule->coor.pbc_flag = 0;
    free(inp_psnparam->molecule->coor.pbc);
    inp_psnparam->molecule->coor.pbc = NULL;
    fprintf( inp_psnparam->FVerboseFile, "# Original PBC    : 1\n"); 
  }
  
  // set and get selection
  strcpy(inp_psnparam->sele.selestring, inp_psnparam->ccMolSeleVect[iMolIndex]);
  GetSele(inp_psnparam->sele.selestring, &inp_psnparam->sele, inp_psnparam->molecule);
  
  // check for target residue(s) and set the reference number of target atoms
  iNumOfSeleRes = 0;
  iNumOfTarget  = 0;
  for(ii=0; ii<inp_psnparam->sele.nselatm; ++ii)
  {
    sprintf(cThisRes, "%s:%s%d",
            inp_psnparam->molecule->rawmol.segId[inp_psnparam->sele.selatm[ii]-1],
            inp_psnparam->molecule->rawmol.restype[inp_psnparam->sele.selatm[ii]-1],
            inp_psnparam->molecule->rawmol.resn[inp_psnparam->sele.selatm[ii]-1]);

    if(strcmp(cThisRes, cLastRes) != 0)
    {
      strcpy(cLastRes, cThisRes);
      iNumOfSeleRes++;
    
      if(strcmp(inp_psnparam->molecule->rawmol.restype[inp_psnparam->sele.selatm[ii]-1], inp_psnparam->cTarget) == 0)
        iNumOfTarget++;
    }
    
    if(strcmp(inp_psnparam->molecule->rawmol.restype[inp_psnparam->sele.selatm[ii]-1], inp_psnparam->cTarget) == 0)
      if(iNumOfTarget == 1 && iSetTargetAtmNumFlag == 1)
        inp_psnparam->iNumOfTargetAtoms++;
  }

  if(inp_psnparam->iVerboseFlag == 1 && iSetTargetAtmNumFlag == 1)
    fprintf(inp_psnparam->FVerboseFile, "Target Ref Atm Num                    : %d\n\n", inp_psnparam->iNumOfTargetAtoms);
  
  // some fancy info about this molecule                                       
  if(inp_psnparam->iVerboseFlag == 1)                                          
  {                                                                            
    fprintf(inp_psnparam->FVerboseFile, "Mol  Num                              : %d\n", iMolIndex + 1);
    fprintf(inp_psnparam->FVerboseFile, "File Name                             : %s\n", inp_psnparam->ccMolFileVect[iMolIndex]);
    fprintf(inp_psnparam->FVerboseFile, "Seg Num                               : %d\n", inp_psnparam->molecule->nSeg);
    fprintf(inp_psnparam->FVerboseFile, "PBC Flag                              : %d\n", inp_psnparam->molecule->coor.pbc_flag);
    fprintf(inp_psnparam->FVerboseFile, "Selection                             : %s\n", inp_psnparam->ccMolSeleVect[iMolIndex]);
    
    if(inp_psnparam->sele.nselatm != 0)
      fprintf(inp_psnparam->FVerboseFile, "Num of Sele Atoms                     : %d\n", inp_psnparam->sele.nselatm);
    else
    {
      inp_psnparam->iNumOfWarning++;
      fprintf(inp_psnparam->FVerboseFile, "Num of Sele Atoms                     : WARNING: 0 atom selected!\n");
    }
    
    if(iNumOfSeleRes != 0)
      fprintf(inp_psnparam->FVerboseFile, "Num of Sele Residues                  : %d\n", iNumOfSeleRes);
    else
    {
      inp_psnparam->iNumOfWarning++;
      fprintf(inp_psnparam->FVerboseFile, "Num of Sele Residues                  : WARNING: 0 residue selected\n");
    }
    
    if(iNumOfTarget != 0)
      fprintf(inp_psnparam->FVerboseFile, "Num of Target Residues                : %d\n", iNumOfTarget);
    else
    {
      inp_psnparam->iNumOfWarning++;
      fprintf(inp_psnparam->FVerboseFile, "Num of Target Residues                : WARNING: 0 target in selection\n");
    }
  }
  
  // calculate the number target interactions
  cLastRes[0]           = '\0';
  cThisRes[0]           = '\0';
  cTestRes[0]           = '\0';
  cBestTarget[0]        = '\0';
  iNumOfInteractions    = 0;
  fMaxInteractions      = 0;
  iTargetId             = 0;
  iNumOfThisTargetAtoms = 0;
  iNumOfIgnoredInt      = 0;
  for(ii=0; ii<inp_psnparam->sele.nselatm; ++ii)
  {
    if(strcmp(inp_psnparam->molecule->rawmol.restype[inp_psnparam->sele.selatm[ii]-1], inp_psnparam->cTarget) == 0)
    {
      sprintf(cThisRes, "%s:%s%d",
              inp_psnparam->molecule->rawmol.segId[inp_psnparam->sele.selatm[ii]-1],
              inp_psnparam->molecule->rawmol.restype[inp_psnparam->sele.selatm[ii]-1],
              inp_psnparam->molecule->rawmol.resn[inp_psnparam->sele.selatm[ii]-1]);
      
      if(strcmp(cThisRes, cLastRes) != 0)
      {
        sprintf(cThisTarget, "#%d/%s/%d", iTargetId, cLastRes, iNumOfThisTargetAtoms);
        if(inp_psnparam->iVerboseFlag == 1 && iTargetId != 0) {
          if(iNumOfThisTargetAtoms == inp_psnparam->iNumOfTargetAtoms)
            fprintf(inp_psnparam->FVerboseFile, "Target %-30s : %d (%d ignored)\n", cThisTarget, iNumOfInteractions, iNumOfIgnoredInt);
          else
          {
            inp_psnparam->iNumOfWarning++;
            fprintf(inp_psnparam->FVerboseFile, "Target %-30s : %d (%d ignored) WARNING: target of different size %d atoms instead of %d\n", cThisTarget, iNumOfInteractions, iNumOfThisTargetAtoms, inp_psnparam->iNumOfTargetAtoms, iNumOfIgnoredInt);
          }
        }

        if(inp_psnparam->iAvgMode == 0)
        {
          if(fMaxInteractions < iNumOfInteractions)
          {
            fMaxInteractions = iNumOfInteractions;
            sprintf(cBestTarget, "#%d/%s", iTargetId, cLastRes);
          }
        }
        
        else
          fMaxInteractions = fMaxInteractions + iNumOfInteractions;
        
        strcpy(cLastRes, cThisRes);
        iNumOfInteractions    = 0;
        iNumOfThisTargetAtoms = 0;
        iNumOfIgnoredInt      = 0;
        iTargetId++;
        iNumOfThisTargetAtoms++;
      }
      
      else
        iNumOfThisTargetAtoms++;
      
      for(jj=0; jj<inp_psnparam->sele.nselatm; ++jj)
      {

        sprintf(cTestRes, "%s:%s%d",
              inp_psnparam->molecule->rawmol.segId[inp_psnparam->sele.selatm[jj]-1],
              inp_psnparam->molecule->rawmol.restype[inp_psnparam->sele.selatm[jj]-1],
              inp_psnparam->molecule->rawmol.resn[inp_psnparam->sele.selatm[jj]-1]);
        
        if(strcmp(cTestRes, cLastRes) != 0)
        {
          if(CheckIgnore(inp_psnparam, inp_psnparam->molecule->rawmol.restype[inp_psnparam->sele.selatm[jj]-1]) == 0)
          {
          
            fTargetCoordX = inp_psnparam->molecule->coor.xcoor[inp_psnparam->sele.selatm[ii]-1];
            fTargetCoordY = inp_psnparam->molecule->coor.ycoor[inp_psnparam->sele.selatm[ii]-1];
            fTargetCoordZ = inp_psnparam->molecule->coor.zcoor[inp_psnparam->sele.selatm[ii]-1];
            
            fTestCoordX = inp_psnparam->molecule->coor.xcoor[inp_psnparam->sele.selatm[jj]-1];
            fTestCoordY = inp_psnparam->molecule->coor.ycoor[inp_psnparam->sele.selatm[jj]-1];
            fTestCoordZ = inp_psnparam->molecule->coor.zcoor[inp_psnparam->sele.selatm[jj]-1];
            
            fDist = DistanceCoor(fTargetCoordX, fTargetCoordY, fTargetCoordZ,
                                 fTestCoordX,   fTestCoordY,   fTestCoordZ,
                                 inp_psnparam->molecule->coor.pbc);
            
            if(fDist <= inp_psnparam->fDistCutoff)
              iNumOfInteractions++;
          }
          
          else
            iNumOfIgnoredInt++;
        }
      }
    }
  }
  
  // process the last target
  if(inp_psnparam->iAvgMode == 0)
  {
    if(fMaxInteractions < iNumOfInteractions)
    {
      fMaxInteractions = iNumOfInteractions;
      sprintf(cBestTarget, "#%d/%s", iTargetId, cLastRes);
    }
  }
  
  else
    fMaxInteractions = fMaxInteractions + iNumOfInteractions;
  
  sprintf(cThisTarget, "#%d/%s/%d", iTargetId, cLastRes, iNumOfThisTargetAtoms);
  if(inp_psnparam->iVerboseFlag == 1 && iTargetId != 0)
  {
    if(iNumOfThisTargetAtoms == inp_psnparam->iNumOfTargetAtoms)
      fprintf(inp_psnparam->FVerboseFile, "Target %-30s : %d (%d ignored)\n", cThisTarget, iNumOfInteractions, iNumOfIgnoredInt);
    else
    {
      inp_psnparam->iNumOfWarning++;
      fprintf(inp_psnparam->FVerboseFile, "Target %-30s : %d (%d ignored) WARNING: target of different size %d atoms instead of %d\n", cThisTarget, iNumOfInteractions, iNumOfThisTargetAtoms, inp_psnparam->iNumOfTargetAtoms, iNumOfIgnoredInt);
    }
  }
  
  if(inp_psnparam->iAvgMode == 1)
  {
    strcpy(cBestTarget, "***averaged***");
    fMaxInteractions = fMaxInteractions / (float) iNumOfTarget;
  }

  if(inp_psnparam->iVerboseFlag == 1) {
    if(inp_psnparam->iAvgMode == 0)
      fprintf(inp_psnparam->FVerboseFile, "Max Num Of Interactions               : %d by target %s\n\n", (int) fMaxInteractions, cBestTarget);
    else
      fprintf(inp_psnparam->FVerboseFile, "Max Num Of Interactions               : %f by target %s\n\n", fMaxInteractions, cBestTarget);
  }
  // de-allocate molecule
  DelMolecule(inp_psnparam->molecule);

  inp_psnparam->fAvgMaxInteractions = inp_psnparam->fAvgMaxInteractions + fMaxInteractions;
}

int CheckIgnore(struct inp_psnparam *inp_psnparam, char *cResToCheck)
{
  int     ii;                                                           // some iterators
  
  for(ii=0; ii<inp_psnparam->iNumOfIgnore; ++ii)
    if(strcmp(cResToCheck, inp_psnparam->ccIgnoreVect[ii]) == 0)
      return 1;
  
  return 0;
}
// =====================================================================
