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
 * GEPOL algorithm is inspired with permission from the work of 
 * JL Pascual-Ahuir, E. Silla, and I. Tunon; Journal of Computational Chemistry, 1994, 15(10)
 * ARVO algorithm is inspired with permission from the work of 
 * J. Busa, J. Dzurina, E. Hayryan, S. Hayryan, C.K. Hu, J. Plavka, I. Pokorny, 
 * J. Skrivrnek, and M.C. Wu; Computer Physics Communications, 2005, 165(1), 59–96.
 * */

#include "wordom.h"
#include "tools.h"
#include "fileio.h"
#include "datahandler.h"
#include "surf.h"
#include "time.h"
#include "float.h"

// === ARVO Functions =========================================================================================
void MakeNeighbors(struct inp_surf *inp_surf)
{
  // Determination of neighbors for all atoms
  
  int     ii, jj;

  inp_surf->piNeighborsStartIndex[0]=0;
  for(ii=0; ii<inp_surf->iNumOfAtoms; ii++)
  {
    inp_surf->piNeighborsNumber[ii]=Neighbors(ii, inp_surf);
    if(inp_surf->piNeighborsNumber[ii]<=0)
    {
      // sphere is subset or there are no neighbors
      inp_surf->piNeighborsStartIndex[ii+1]=inp_surf->piNeighborsStartIndex[ii];
    }
    else
    {
      // there are neighbors
      inp_surf->piNeighborsStartIndex[ii+1]=inp_surf->piNeighborsStartIndex[ii]+inp_surf->piNeighborsNumber[ii];

      for(jj=0; jj<inp_surf->piNeighborsNumber[ii]; jj++)
        inp_surf->piNeighborsIndices[inp_surf->piNeighborsStartIndex[ii]+jj]=inp_surf->piIndicesVector[jj+1];

    }
  }
}

int Neighbors(int iIndex, struct inp_surf *inp_surf)
{
  /*
     If ith sphere is a subset of other sphere, index_number(i)=-1 and we change its 
     radius in matrix spheres to -radius !!!
     If some other sphere is subset of ith sphere, than we change its radius to -radius !!!
  */

  int     ii;
  double  dDist;
  int     iNumOfNeighbors;
  
  iNumOfNeighbors=0;
  for(ii=0; ii<inp_surf->iNumOfAtoms; ii++)
  {
    if(ii!=iIndex)
    {
      if(fabs(inp_surf->ppfWithinAtomsCoord[iIndex][0]-inp_surf->ppfWithinAtomsCoord[ii][0])<
        (inp_surf->ppfWithinAtomsCoord[iIndex][3]+inp_surf->ppfWithinAtomsCoord[ii][3]))
      {
        dDist=sqrt(
                pow((inp_surf->ppfWithinAtomsCoord[iIndex][0]-inp_surf->ppfWithinAtomsCoord[ii][0]), 2)+
                pow((inp_surf->ppfWithinAtomsCoord[iIndex][1]-inp_surf->ppfWithinAtomsCoord[ii][1]), 2)+
                pow((inp_surf->ppfWithinAtomsCoord[iIndex][2]-inp_surf->ppfWithinAtomsCoord[ii][2]), 2));
		
        if(dDist<(inp_surf->ppfWithinAtomsCoord[iIndex][3]+inp_surf->ppfWithinAtomsCoord[ii][3]))
        {
          if((dDist+inp_surf->ppfWithinAtomsCoord[iIndex][3])<=inp_surf->ppfWithinAtomsCoord[ii][3])
          {
            //iIndex-th sphere is inside of other sphere
            iNumOfNeighbors=-1;
            break;
          }
          else if((dDist+inp_surf->ppfWithinAtomsCoord[ii][3])>inp_surf->ppfWithinAtomsCoord[iIndex][3])
          {
            //ii-th sphere is neighbor
            iNumOfNeighbors++;
            inp_surf->piIndicesVector[iNumOfNeighbors]=ii;
          }
        }
      }
    }
  }
  return iNumOfNeighbors;
}

int NorthPoleTest(struct inp_surf *inp_surf)
{
  /*
    Here we check, that North Pole of no sphere lies on other neighbor sphere
    
  */
  
  int     ii, jj;
  int     iNPTtest, iTmp;
  double  dMin, dDist;
  double  dNPTCritical;
  
  /*
    dNPTCritical -> Critical value for North Pole test
    if the smallest over all atoms distance of the North Poles
    to the surface of other atoms is smaller than this value
    the molecule is rotated by the random angle
  */
  dNPTCritical=0.000001;
  // dMin -> Square of minimal distance of the North Pole to neighbor sphere surface
  dMin=10000.0;

  for(ii=0; ii<inp_surf->iNumOfAtoms; ii++)
  {

    for(jj=1; jj<inp_surf->piNeighborsNumber[ii]; jj++)
    {
      // jj-th neighbor index
      iTmp=inp_surf->piNeighborsIndices[inp_surf->piNeighborsStartIndex[ii]+jj];
      dDist=fabs(sqrt(
                    pow(inp_surf->ppfWithinAtomsCoord[ii][0]-inp_surf->ppfWithinAtomsCoord[iTmp][0], 2)+
                    pow(inp_surf->ppfWithinAtomsCoord[ii][1]-inp_surf->ppfWithinAtomsCoord[iTmp][1], 2)+
                    pow(inp_surf->ppfWithinAtomsCoord[ii][2]+inp_surf->ppfWithinAtomsCoord[ii][3]-
                    inp_surf->ppfWithinAtomsCoord[iTmp][2], 2))-inp_surf->ppfWithinAtomsCoord[iTmp][3]);
      
      if(dDist<dMin)
        dMin=dDist;
    }
  }

  if(dMin<dNPTCritical)
    iNPTtest=0;

  else
    iNPTtest=1;

  return iNPTtest;
}

void MolRotation(struct inp_surf *inp_surf)
{
  /*
    Random rotation of molecule about the y-axis after bad North Pole test
    Some North Pole is near other spheres surface
  */

  int     ii;
  float   fRotation=0.0, fTmp;
  
  while(fRotation==0.0)
    fRotation=rand() / (RAND_MAX / (360) +1.); // 0..360

  fTmp=sqrt(1.0-fRotation*fRotation);
  for(ii=0;ii<inp_surf->iNumOfAtoms;ii++)
  {
    inp_surf->ppfWithinAtomsCoord[ii][0]=fTmp*inp_surf->ppfWithinAtomsCoord[ii][0]-fRotation*inp_surf->ppfWithinAtomsCoord[ii][2];
    inp_surf->ppfWithinAtomsCoord[ii][2]=fRotation*inp_surf->ppfWithinAtomsCoord[ii][0]+fTmp*inp_surf->ppfWithinAtomsCoord[ii][2];
  }
}

double GetSurf(int iAtomNumber, struct inp_surf *inp_surf)
{
  int      ii, jj;
  double   PI = 3.14159265358979323846264;
  double   dAtomSurf, dAtomCoordZ, dAtomRadius;
  int      iNumOfLocalAtoms, iNumOfArcs, iNumOfPos;
  
  dAtomSurf=dAtomCoordZ=dAtomRadius=0.0;
  iNumOfLocalAtoms=iNumOfArcs=iNumOfPos=0;
  
  // Determination of ii-th atom's neighbors
  if(inp_surf->piNeighborsNumber[iAtomNumber]<0)
  {
    // ii-th atom is subset of other atoms
    // radius will be done negative
    dAtomSurf=0.0;
  }
  else if(inp_surf->piNeighborsNumber[iAtomNumber]==0)
  {
    // there are no neighbors
    dAtomSurf=(4.0 * PI * pow(inp_surf->ppfWithinAtomsCoord[iAtomNumber][3], 2));
  }
  else
  {
    // there are neighbors
    iNumOfLocalAtoms=inp_surf->piNeighborsNumber[iAtomNumber]+1;
    inp_surf->piIndicesVector2[0]=iAtomNumber;
    for(ii=0; ii<iNumOfLocalAtoms-1; ii++)
      inp_surf->piIndicesVector2[ii+1]=inp_surf->piNeighborsIndices[inp_surf->piNeighborsStartIndex[iAtomNumber]+ii];
    
    // fill ppdAtomNeighbors
    for(ii=0; ii<iNumOfLocalAtoms; ii++)
    {
      inp_surf->ppdAtomNeighbors[ii][0]=inp_surf->ppfWithinAtomsCoord[inp_surf->piIndicesVector2[ii]][0];
      inp_surf->ppdAtomNeighbors[ii][1]=inp_surf->ppfWithinAtomsCoord[inp_surf->piIndicesVector2[ii]][1];
      inp_surf->ppdAtomNeighbors[ii][2]=inp_surf->ppfWithinAtomsCoord[inp_surf->piIndicesVector2[ii]][2];
      inp_surf->ppdAtomNeighbors[ii][3]=inp_surf->ppfWithinAtomsCoord[inp_surf->piIndicesVector2[ii]][3];
    }
        
    MakeCircles(inp_surf, iNumOfLocalAtoms);
    iNumOfArcs=CircleToArcs(inp_surf, iNumOfLocalAtoms);
    iNumOfPos=0;
    for(jj=0; jj<iNumOfLocalAtoms-1; jj++)
    {
      if(inp_surf->ppdCircles[jj][3]>0.0)
      {
        iNumOfPos++;
      }
    }
    
    dAtomCoordZ=inp_surf->ppdAtomNeighbors[0][2];
    dAtomRadius=inp_surf->ppdAtomNeighbors[0][3];
    if(iNumOfPos>0)
    {
      // there exists positive oriented circle
      dAtomSurf=dAtomSurf+AvgIntegral(iNumOfArcs, dAtomCoordZ, dAtomRadius, inp_surf);
    }
    else
    {
      // All circles are negative oriented - we compute complement
      dAtomSurf=dAtomSurf+AvgIntegral(iNumOfArcs, dAtomCoordZ, dAtomRadius, inp_surf)+4.0*PI*pow(inp_surf->ppdAtomNeighbors[0][3], 2);
    }
  }

  return dAtomSurf;
}

void MakeCircles(struct inp_surf *inp_surf, int iNumOfLocalAtoms)
{
  /*
     Preparing circles structure for 1st sphere in array circles
     according to the paper Hayrjan, Dzurina, Plavka, Busa
     
     inp_surf->ppdCircles[ii][0]=ti 
     inp_surf->ppdCircles[ii][1]=si    - ith circle's center coordinates
     inp_surf->ppdCircles[ii][2]=ri    - ith circle's radius 
     inp_surf->ppdCircles[ii][3]=+1/-1 - circle orientation 
  */

  int     ii;
  double  dRad, dDX, dDY;
  double  dA, dB, dC, dD;
  
  
  dRad=inp_surf->ppdAtomNeighbors[0][3];
  for(ii=0; ii<iNumOfLocalAtoms-1; ii++)
  {
    dDX=inp_surf->ppdAtomNeighbors[0][0]-inp_surf->ppdAtomNeighbors[ii+1][0];
    dDY=inp_surf->ppdAtomNeighbors[0][1]-inp_surf->ppdAtomNeighbors[ii+1][1];
    dA=((dDX*dDX)+(dDY*dDY)+pow((inp_surf->ppdAtomNeighbors[0][2]+dRad-inp_surf->ppdAtomNeighbors[ii+1][2]), 2)-pow(inp_surf->ppdAtomNeighbors[ii+1][3], 2));
    dB=(double)8.0*dRad*dRad*dDX;
    dC=(double)8.0*dRad*dRad*dDY;
    dD=((double)4.0*dRad*dRad*(dDX*dDX+dDY*dDY+pow((inp_surf->ppdAtomNeighbors[0][2]-dRad-inp_surf->ppdAtomNeighbors[ii+1][2]), 2)-pow(inp_surf->ppdAtomNeighbors[ii+1][3], 2)));
    inp_surf->ppdCircles[ii][0]=-dB/(2.0*dA);
    inp_surf->ppdCircles[ii][1]=-dC/(2.0*dA);
    inp_surf->ppdCircles[ii][2]=sqrt((dB*dB+dC*dC-(double)4.0*dA*dD)/((double)4.0*dA*dA));
    if(dA>0.0)
      inp_surf->ppdCircles[ii][3]=-1;

    else
      inp_surf->ppdCircles[ii][3]=1;
  }
}

int CircleToArcs(struct inp_surf *inp_surf, int iNumOfLocalAtoms)
{
  /*
     Computing integration arcs
     
     inp_surf->ppdArcs[i][1]=ci        - corresponding circle index
     inp_surf->ppdArcs[i][2]=sigma     - starting arc angle 
     inp_surf->ppdArcs[i][3]=delta     - oriented arc angle

     Arcs (with their orientation) are parts of circles, which
     bounds are circles intersection points. If the center of
     arc lies inside all other positive and outside all other
     negative circles, then we will put it inside arcs structure
  */

  int     ii, jj, kk;
  int     iNumOfArcs=0, iNumOfNewArcs=0;
  double  PI = 3.14159265358979323846264;
  
  if(iNumOfLocalAtoms==2)
  {
    // There is only 1 circle    
    iNumOfArcs=1;
    inp_surf->ppdArcs[0][0]=0.0;
    inp_surf->ppdArcs[0][1]=0.0;
    inp_surf->ppdArcs[0][2]=2.0*PI*inp_surf->ppdCircles[0][3];
  }
  else
  {
    for(ii=0; ii<iNumOfLocalAtoms-1; ii++)
    {
      iNumOfNewArcs=MakeNewArcs(ii, iNumOfLocalAtoms, inp_surf);
      if(iNumOfNewArcs>0)
      {
        for(jj=0; jj<iNumOfNewArcs; jj++)
        {
          for(kk=0;kk<3; kk++)
            inp_surf->ppdArcs[iNumOfArcs+jj][kk]=inp_surf->ppdNewArcs[jj][kk];
        }
        iNumOfArcs=iNumOfArcs+iNumOfNewArcs;
      }
    }
  }
  return iNumOfArcs;
}

int MakeNewArcs(int iIndex, int iNumOfLocalAtoms, struct inp_surf *inp_surf)
{
  /*
    Function prepares arcs, which are part of i-th circle
    in circle structure circles.
    Interesting are these arcs, which are inside other positive
    circles or outside other negative circles
    
    Matrix arcsnew in each row has elements

        inp_surf->ppdNewArcs[i][1]=ic    - ic is the index of arc-circle in circle 
        inp_surf->ppdNewArcs[i][2]=sigma - sigma is the starting angle of arc
        inp_surf->ppdNewArcs[i][3]=delta - delta is oriented arc angle

  */
  
  int     ii, jj;
  double  PI = 3.14159265358979323846264;
  int     iNumOfArcs, iNumberCond;
  int     iNumOfAngles, iTmpNumOfAngles;
  double  dRefT, dRefS, dRefR;
  double  dCompT, dCompS, dCompR;
  double  dDist;
  
  
  iNumOfArcs=iNumberCond=iNumOfAngles=iTmpNumOfAngles=0;
  dRefT=dRefS=dRefR=dCompT=dCompS=dCompR=dDist=0.0;
  
  dRefT=inp_surf->ppdCircles[iIndex][0];
  dRefS=inp_surf->ppdCircles[iIndex][1];
  dRefR=inp_surf->ppdCircles[iIndex][2];
  for(ii=0; ii<iNumOfLocalAtoms-1; ii++)
  {
    // Composition of angles vector, consisting of intersection points
    if(ii!=iIndex)
    {
      dCompT=inp_surf->ppdCircles[ii][0];
      dCompS=inp_surf->ppdCircles[ii][1];
      dCompR=inp_surf->ppdCircles[ii][2];
      dDist=sqrt(pow((dRefT-dCompT), 2)+pow((dRefS-dCompS), 2));
      if((dDist < (dCompR+dRefR)) && (fabs(dCompR-dRefR)<dDist))
      {
        // 2 intersection points exist
        CirclesIntersection(iIndex, ii, inp_surf);
        inp_surf->pdAngles[iNumOfAngles]=inp_surf->dIntPointAngleA1;
        inp_surf->pdAngles[iNumOfAngles+1]=inp_surf->dIntPointAngleA2;
        iNumOfAngles=iNumOfAngles+2;
      }
    }
  }
  
  if(iNumOfAngles==0)
  {
    // there are no double intersections of iIndex-th circles with others
    // if iIndex-th circle is inside of all other positive and outside of 
    // all other negative circles, it will be new arc
    iNumberCond=0;
    for(ii=0; ii<iNumOfLocalAtoms-1; ii++)
    {
      if(ii!=iIndex)
        iNumberCond=iNumberCond+CircleInCircle(iIndex, ii, inp_surf);
    }
    if(iNumberCond==(iNumOfLocalAtoms-2))
    {
      // all conditions hold
      iNumOfArcs=1;
      inp_surf->ppdNewArcs[0][0]=iIndex;
      inp_surf->ppdNewArcs[0][1]=0.0;
      inp_surf->ppdNewArcs[0][2]=(double)2.0*PI*inp_surf->ppdCircles[iIndex][3];
    }
  }
  else
  {
    // there are double intersection points
    if(inp_surf->ppdCircles[iIndex][3] > 0.0)
    {
      DirectSort(iNumOfAngles, inp_surf);
    }
    else
    {
      InvertedSort(iNumOfAngles, inp_surf);
    }
    iTmpNumOfAngles=DeleteEqual(iNumOfAngles, inp_surf);
    iNumOfAngles=iTmpNumOfAngles;    
    for(ii=0; ii<(iTmpNumOfAngles); ii++)
    {
      iNumberCond=0;
      for(jj=0; jj<(iNumOfLocalAtoms-1); jj++)
      {
        if(jj!=iIndex)
        {
          dCompT=dRefT+dRefR*cos((inp_surf->pdAngles[ii]+inp_surf->pdAngles[ii+1])/2.0);
          dCompS=dRefS+dRefR*sin((inp_surf->pdAngles[ii]+inp_surf->pdAngles[ii+1])/2.0);
          iNumberCond=iNumberCond+PointInCircle(dCompT, dCompS, jj, inp_surf);
        }
      }
      if(iNumberCond==(iNumOfLocalAtoms-2))
      {
        // all conditions hold
        iNumOfArcs++;
        inp_surf->ppdNewArcs[iNumOfArcs-1][0]=iIndex;
        inp_surf->ppdNewArcs[iNumOfArcs-1][1]=inp_surf->pdAngles[ii];
        inp_surf->ppdNewArcs[iNumOfArcs-1][2]=inp_surf->pdAngles[ii+1]-inp_surf->pdAngles[ii];
      }
    }
    iNumberCond=0;
    for(ii=0; ii<(iNumOfLocalAtoms-1); ii++)
    {
      if(ii!=iIndex)
      {
        dCompT=dRefT+dRefR*cos(((inp_surf->pdAngles[0]+(double)2.0*PI+inp_surf->pdAngles[iTmpNumOfAngles])/(double)2.0));
        dCompS=dRefS+dRefR*sin(((inp_surf->pdAngles[0]+(double)2.0*PI+inp_surf->pdAngles[iTmpNumOfAngles])/(double)2.0));
        iNumberCond=iNumberCond+PointInCircle(dCompT, dCompS, ii, inp_surf);
      }
    }
    if(iNumberCond==(iNumOfLocalAtoms-2))
    {
      // all conditions hold
      iNumOfArcs++;
      inp_surf->ppdNewArcs[iNumOfArcs-1][0]=iIndex;
      inp_surf->ppdNewArcs[iNumOfArcs-1][1]=inp_surf->pdAngles[iTmpNumOfAngles];
      inp_surf->ppdNewArcs[iNumOfArcs-1][2]=inp_surf->pdAngles[0]+inp_surf->ppdCircles[iIndex][3]*(double)2.0*PI-inp_surf->pdAngles[iTmpNumOfAngles];
    }
  }
  
  return iNumOfArcs;
}

void CirclesIntersection(int iIndex1, int iIndex2, struct inp_surf *inp_surf)
{

  /*
    Function returns angles of two intersection points
    of circles with indices ic1 and ic2 in circles structure circles
    (we will use it ONLY IN CASE, WHEN 2 INTERSECTION POINTS EXIST!!!)

    dIntPointAngleA1 and dIntPointAngleA2 are corresponding angles with respect
    to the center of 1st circle and dIntPointAngleB1 and dIntPointAngleB2 are
    corresponding angles with respect to the center of 2nd circle
  */
  
  double      PI = 3.14159265358979323846264;
  double      dIdx1Cent1, dIdx1Cent2, dIdx1Rad;   // 1st Circle center and radius
  double      dIdx2Cent1, dIdx2Cent2, dIdx2Rad;   // 2nd Circle center and radius
  double      dTmpA, dTmpB, dTmpC, dTmpD;
  long double dCirclesCompCritical=9E-13;         // The critical value for comparison of t_1 and t_2 coordinates of
                                                  // two circles in the plane when the intersection points are calculated
  dIdx1Cent1=dIdx1Cent2=dIdx1Rad=0.0;
  dIdx2Cent1=dIdx2Cent2=dIdx2Rad=0.0;
  dTmpA=dTmpB=dTmpC=dTmpD=0.0;
  
  dIdx1Cent1=inp_surf->ppdCircles[iIndex1][0];
  dIdx1Cent2=inp_surf->ppdCircles[iIndex1][1];
  dIdx1Rad  =inp_surf->ppdCircles[iIndex1][2];
  dIdx2Cent1=inp_surf->ppdCircles[iIndex2][0];
  dIdx2Cent2=inp_surf->ppdCircles[iIndex2][1];
  dIdx2Rad  =inp_surf->ppdCircles[iIndex2][2];
						 
  if(fabs(dIdx2Cent1-dIdx1Cent1) < dCirclesCompCritical)
  {
    // dIdx2Cent1 == dIdx1Cent1
    dTmpB=((dIdx1Rad*dIdx1Rad-dIdx2Rad*dIdx2Rad)/(dIdx2Cent2-dIdx1Cent2)-(dIdx2Cent2-dIdx1Cent2))/(double)2.0;
    dTmpA=sqrt(dIdx2Rad*dIdx2Rad-dTmpB*dTmpB);
    if(dTmpB==0)
    {
      inp_surf->dIntPointAngleB1=(double)0.0;
      inp_surf->dIntPointAngleB2=PI;
    }
    else if(dTmpB>0)
    {
      inp_surf->dIntPointAngleB1=atan(fabs(dTmpB/dTmpA));
      inp_surf->dIntPointAngleB2=PI-inp_surf->dIntPointAngleB1;
    }
    else
    {
      inp_surf->dIntPointAngleB1=PI+atan(fabs(dTmpB/dTmpA));
      inp_surf->dIntPointAngleB2=(double)3.0*PI-inp_surf->dIntPointAngleB1;
    }
    
    dTmpB=dTmpB+dIdx2Cent2-dIdx1Cent2;
    if(dTmpB==0)
    {
      inp_surf->dIntPointAngleA1=(double)0.0;
      inp_surf->dIntPointAngleA1=PI;
    }
    else if(dTmpB>0)
    {
      inp_surf->dIntPointAngleA1=atan(fabs(dTmpB/dTmpA));
      inp_surf->dIntPointAngleA2=PI-inp_surf->dIntPointAngleA1;
    }
    else
    {
      inp_surf->dIntPointAngleA1=PI+atan(fabs(dTmpB/dTmpA));
      inp_surf->dIntPointAngleA2=(double)3.0*PI-inp_surf->dIntPointAngleA1;
    }  
  }
  else
  {
    // dIdx2Cent1 != dIdx1Cent1
    dTmpC=((dIdx1Rad*dIdx1Rad-dIdx2Rad*dIdx2Rad-pow((dIdx2Cent2-dIdx1Cent2),2))/(dIdx2Cent1-dIdx1Cent1)-(dIdx2Cent1-dIdx1Cent1))/2.0;
    dTmpD=(dIdx1Cent2-dIdx2Cent2)/(dIdx2Cent1-dIdx1Cent1);
    dTmpB=(-dTmpC*dTmpD+sqrt((dTmpD*dTmpD+1.0)*dIdx2Rad*dIdx2Rad-dTmpC*dTmpC))/(dTmpD*dTmpD+1.0);
    dTmpA=dTmpC+dTmpD*dTmpB;
    
    if(dTmpA==0.0)
    {
      if(dTmpB>0.0)
        inp_surf->dIntPointAngleB1=PI/(double)2.0;
      else
        inp_surf->dIntPointAngleB1=-PI/(double)2.0;
    }
    
    else if(dTmpA>0)
      inp_surf->dIntPointAngleB1=atan(dTmpB/dTmpA);

    else
      inp_surf->dIntPointAngleB1=PI+atan(dTmpB/dTmpA);

    dTmpB=dTmpB+dIdx2Cent2-dIdx1Cent2;
    dTmpA=dTmpA+dIdx2Cent1-dIdx1Cent1;
    if(dTmpA==0.0)
    {
      if(dTmpB>0.0)
        inp_surf->dIntPointAngleA1=PI/2.0;

      else
        inp_surf->dIntPointAngleA1=-PI/2.0;
    }

    else if(dTmpA>0.0)
      inp_surf->dIntPointAngleA1=atan(dTmpB/dTmpA);

    else
      inp_surf->dIntPointAngleA1=PI+atan(dTmpB/dTmpA);
      
    dTmpB=((-dTmpC*dTmpD-sqrt((dTmpD*dTmpD+1.0)*dIdx2Rad*dIdx2Rad-dTmpC*dTmpC))/(dTmpD*dTmpD+1.0));
    dTmpA=dTmpC+dTmpD*dTmpB;
    if(dTmpA==0.0)
    {
      if(dTmpB>0.0)
	      inp_surf->dIntPointAngleB2=PI/2.0;

      else
      	inp_surf->dIntPointAngleB2=-PI/2.0;
    }

    else if(dTmpA>0.0)
      inp_surf->dIntPointAngleB2=atan(dTmpB/dTmpA);

    else
      inp_surf->dIntPointAngleB2=PI+atan(dTmpB/dTmpA);

    dTmpB=dTmpB+dIdx2Cent2-dIdx1Cent2;
    dTmpA=dTmpA+dIdx2Cent1-dIdx1Cent1;
    
    if(dTmpA==0.0)
    {
      if(dTmpB>0.0)
        inp_surf->dIntPointAngleA2=PI/2.0;

      else
        inp_surf->dIntPointAngleA2=-PI/2.0;
    }
    else if(dTmpA>0.0)
      inp_surf->dIntPointAngleA2=atan(dTmpB/dTmpA);

    else
      inp_surf->dIntPointAngleA2=PI+atan(dTmpB/dTmpA);
  }
  
  if (inp_surf->dIntPointAngleA1<0.0)
    inp_surf->dIntPointAngleA1=inp_surf->dIntPointAngleA1+2.0*PI;
    
  if (inp_surf->dIntPointAngleA2<0.0)
    inp_surf->dIntPointAngleA2=inp_surf->dIntPointAngleA2+2.0*PI;
    
  if (inp_surf->dIntPointAngleB1<0.0)
    inp_surf->dIntPointAngleB1=inp_surf->dIntPointAngleB1+2.0*PI;
    
  if (inp_surf->dIntPointAngleB2<0.0)
    inp_surf->dIntPointAngleB2=inp_surf->dIntPointAngleB2+2.0*PI;
  
}


int CircleInCircle(int iIndex1, int iIndex2, struct inp_surf *inp_surf)
{
  /*
     1  - if i-th circle is inside k-th positive circle or outside k-th negative circle
     0  - otherwise

     WE KNOW, THAT CIRCLES HAVE LESS THAN 2 INTERSECTION POINTS !!!
  */
 
  double      dDist=0.0;
  int         iCircleInCircle=0;
   
  dDist=sqrt(pow((inp_surf->ppdCircles[iIndex1][0]+inp_surf->ppdCircles[iIndex1][2]-inp_surf->ppdCircles[iIndex2][0]), 2)+
	pow((inp_surf->ppdCircles[iIndex1][1]-inp_surf->ppdCircles[iIndex2][1]), 2));
  
  if(dDist<inp_surf->ppdCircles[iIndex2][2])
  {
    if(inp_surf->ppdCircles[iIndex2][3]>0.0)
    {
      iCircleInCircle=1;
    }
    else
    {
      iCircleInCircle=0;
    }    
  }
  else if(dDist>inp_surf->ppdCircles[iIndex2][3])
  {
    if(inp_surf->ppdCircles[iIndex2][3]>0.0)
    {
      iCircleInCircle=0;
    }
    else
    {
      iCircleInCircle=1;
    }
  }
  else
  {
    // right point on iIndex2-th circle - touching of circles
    dDist=sqrt(pow((inp_surf->ppdCircles[iIndex1][0]-inp_surf->ppdCircles[iIndex2][0]),2)+
          pow((inp_surf->ppdCircles[iIndex1][1]-inp_surf->ppdCircles[iIndex2][1]), 2));
    if(dDist<inp_surf->ppdCircles[iIndex2][2])
    {
      if(inp_surf->ppdCircles[iIndex2][3]>0.0)
        iCircleInCircle=1;

      else
        iCircleInCircle=0;

    }
    else
    {
      if(inp_surf->ppdCircles[iIndex2][3]>0.0)
        iCircleInCircle=0;

      else
        iCircleInCircle=1;
    }
  }

  return iCircleInCircle;
}

void DirectSort(int iNumOfAngles, struct inp_surf *inp_surf)
{
  /*
    Sorting array angles in increasing order
    iNumOfAngles is the angles array length
  */

  int         ii, jj, iIndex;
  double      dAngleMax;
  
  for(ii=0; ii<iNumOfAngles-1; ii++)
  {
    iIndex=ii;
    dAngleMax=inp_surf->pdAngles[ii];
    
    for(jj=ii+1; jj<iNumOfAngles; jj++)
    {
      if(dAngleMax>inp_surf->pdAngles[jj])
      {
        iIndex=jj;
      	dAngleMax=inp_surf->pdAngles[jj];
      }
    }
    if(iIndex!=ii)
    {
      inp_surf->pdAngles[iIndex]=inp_surf->pdAngles[ii];
      inp_surf->pdAngles[ii]=dAngleMax;
    }
  }
}


void InvertedSort(int iNumOfAngles, struct inp_surf *inp_surf)
{
  /*
    Sorting array angles in increasing order
    iNumOfAngles is the angles array length
  */
  int         ii, jj, iIndex;
  double      dAngleMin;
  
  for(ii=0; ii<iNumOfAngles-1; ii++)
  {
    iIndex=ii;
    dAngleMin=inp_surf->pdAngles[ii];
    for(jj=ii+1; jj<iNumOfAngles; jj++)
    {
      if(dAngleMin<inp_surf->pdAngles[jj])
      {
        iIndex=jj;
      	dAngleMin=inp_surf->pdAngles[jj];
      }
    }
    if(iIndex!=ii)
    {
      inp_surf->pdAngles[iIndex]=inp_surf->pdAngles[ii];
      inp_surf->pdAngles[ii]=dAngleMin;
    }
  }
}

int DeleteEqual(int iNumOfAngles, struct inp_surf *inp_surf)
{

  /*
    Deletion of "equal" (to some precision dAngleCompCritical)
    angles in sorted vector inp_surf->pdAngles
  */

  int         ii, iIndex;
  double      dAngle;
  long double dAngleCompCritical=1E-12;          // The critical value for comparison of two angles in delete_equal function
                                                 // if two points on the circle are close to each other, they are declared equal and only
                                                 // one point is left.

  iIndex=0;
  dAngle=inp_surf->pdAngles[0];
  inp_surf->pdNewAngles[0]=dAngle;
  for(ii=1; ii<iNumOfAngles; ii++)
  {
    if(fabs(inp_surf->pdAngles[ii]-dAngle)>dAngleCompCritical)
    {
      dAngle=inp_surf->pdAngles[ii];
      iIndex++;
      inp_surf->pdNewAngles[iIndex]=dAngle;
    }
  }
  
  for(ii=0; ii<iIndex; ii++)
  {
    inp_surf->pdAngles[ii]=inp_surf->pdNewAngles[ii];
  }
  
  return iIndex;
}

int PointInCircle(double dCompT, double dCompS, int iJJIndex, struct inp_surf *inp_surf)
{

  /*
    1  - if point (t,s) is inside k-th positive circle or outside k-th negative circle
    0  - otherwise
    WE KNOW, THAT POINT IS NOT ON THE CIRCLE !!!
  */

  int         iNumOfPointInCircle;
  double      dDist;
  
  dDist=sqrt(pow((dCompT-inp_surf->ppdCircles[iJJIndex][0]), 2)+pow((dCompS-inp_surf->ppdCircles[iJJIndex][1]), 2));
  if(dDist<inp_surf->ppdCircles[iJJIndex][2])
  {
    if(inp_surf->ppdCircles[iJJIndex][3]>0)
      iNumOfPointInCircle=1;

    else
      iNumOfPointInCircle=0;
  }
  else
  {
    if(inp_surf->ppdCircles[iJJIndex][3]>0)
      iNumOfPointInCircle=0;

    else
      iNumOfPointInCircle=1;
  }
  
  return iNumOfPointInCircle;
}

double AvgIntegral(int iNumOfArcs, double dAtomCoordZ, double dAtomRadius, struct inp_surf *inp_surf)
{
  /*
    Computing integrals over arcs given in arc structure
    according to paper Hayrian, Dzurina, Plavka, Busa
  */
  
  int     kk;
  double  dTmpA, dTmpB, dTmpC, dTmpS, dTmpT, dTmpR, dTmpS2, dTmpRR;
  double  dVI1, dVJ1, dDeltaA;
  double  dAL, dBE;
  double  PI = 3.14159265358979323846264, dAtomSurf;
  double  dTwoPI = 1E-12;
  
  dAtomSurf=0.0;
  
  
  for(kk=0; kk<iNumOfArcs; kk++)
  {
    dTmpT=inp_surf->ppdCircles[(int)inp_surf->ppdArcs[kk][0]][0];
    dTmpS=inp_surf->ppdCircles[(int)inp_surf->ppdArcs[kk][0]][1];
    dTmpR=inp_surf->ppdCircles[(int)inp_surf->ppdArcs[kk][0]][2];

    dTmpA=(4.0*dAtomRadius*dAtomRadius+dTmpT*dTmpT+dTmpS*dTmpS+dTmpR*dTmpR)/2.0;
    dTmpB=dTmpT*dTmpR;
    dTmpC=dTmpS*dTmpR;
    dTmpS2=sqrt(dTmpA*dTmpA-dTmpB*dTmpB-dTmpC*dTmpC);
    dTmpRR=dTmpR*dTmpR-dTmpA;
    
    if(fabs(fabs(inp_surf->ppdArcs[kk][2])-2.0*PI)<dTwoPI)
    {
      dVI1=2.0*PI/dTmpS2;
      //dVI2=2.0*PI*dTmpA/pow(dTmpS2, 3);
      //dVI3=PI*(2.0*dTmpA*dTmpA+dTmpB*dTmpB+dTmpC*dTmpC)/pow(dTmpS2, 5);
      dVJ1=PI+dTmpRR/2.0*dVI1;
      //dVJ2=(dVI1+dTmpRR*dVI2)/4.0;
      //dVJ3=(dVI2+dTmpRR*dVI3)/8.0;
      dDeltaA=2.0*dVJ1*pow(dAtomRadius, 2);
      
      if(inp_surf->ppdArcs[kk][2]<0)
	      dDeltaA=-dDeltaA;
        
      dAtomSurf=dAtomSurf+dDeltaA;
    }
    else
    {
      // integration over arcs
      if(inp_surf->ppdArcs[kk][2]<0)
      {
	      dAL=inp_surf->ppdArcs[kk][1]+inp_surf->ppdArcs[kk][2];
	      dBE=inp_surf->ppdArcs[kk][1];
      }
      else
      {
	      dBE=inp_surf->ppdArcs[kk][1]+inp_surf->ppdArcs[kk][2];
        dAL=inp_surf->ppdArcs[kk][1];
      }
      dVI1=2.0*(PI/2.0-atan((dTmpA*cos((dBE-dAL)/2.0)+dTmpB*cos((dAL+dBE)/2.0)+
           dTmpC*sin((dAL+dBE)/2.0))/(dTmpS2*sin((dBE-dAL)/2.0))))/dTmpS2;
      
      //dSB=sin(dBE);
      //dCB=cos(dBE);
      //dSA=sin(dAL);
      //dCA=cos(dAL);
      
      //dVI2=(Fract(dTmpA, dTmpB, dTmpC, dSB, dCB, 1)-Fract(dTmpA, dTmpB, dTmpC, dSA, dCA, 1)+dTmpA*dVI1)/(dTmpS2*dTmpS2);
      //dVI3=(Fract(dTmpA, dTmpB, dTmpC, dSB, dCB, 2)-Fract(dTmpA, dTmpB, dTmpC, dSA, dCA, 2)+
      //     (Fract(dTmpA, dTmpB, dTmpC, dSB, dCB, 1)-Fract(dTmpA, dTmpB, dTmpC, dSA, dCA, 1))/
      //     dTmpA+(2.0+dTmpA*dTmpA+dTmpB*dTmpB+dTmpC*dTmpC)*dVI2/dTmpA)/(2.0*dTmpS2*dTmpS2);

      dVJ1=((dBE-dAL)+dTmpRR*dVI1)/2.0;
      //dVJ2=(dVI1+dTmpRR*dVI2)/4.0;
      //dVJ3=(dVI2+dTmpRR*dVI3)/8.0;
      dDeltaA=2.0*dVJ1*pow(dAtomRadius, 2);
      
      if(inp_surf->ppdArcs[kk][2]<0)
        dDeltaA=-dDeltaA;

      dAtomSurf=dAtomSurf+dDeltaA;
    }
  }
  
  return dAtomSurf;
}

double Fract(double dA, double dB, double dC, double dSinPhi, double dCosPhi, int iExp)
{
  /*
    Fraction evaluation for integral
  */

  return (-dB*dSinPhi+dC*dCosPhi)/pow((dA+dB*dCosPhi+dC*dSinPhi), iExp);
}

double ComputeARVOSurf(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring) //FILE *OutputFile)
{
  int     ii, jj, iNumOfWithinAtoms, iRotCont;
  float   fDist;
  double  dSurf=0.0, dTmpSurf=0.0;
  
  for(ii=0; ii<inp_surf->iNumOfRealAtoms; ii++)
  {
    inp_surf->piWithinAtoms[ii]=0;
  }

  iNumOfWithinAtoms=-1;
  for(ii=0; ii<inp_surf->sele1.nselatm; ii++)
  {
    iNumOfWithinAtoms++;
    inp_surf->piWithinAtoms[inp_surf->sele1.selatm[ii]-1]=1;
    inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][0]=trj_crd->xcoor[inp_surf->sele1.selatm[ii]-1];
    inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][1]=trj_crd->ycoor[inp_surf->sele1.selatm[ii]-1];
    inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][2]=trj_crd->zcoor[inp_surf->sele1.selatm[ii]-1];
    inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][3]=inp_surf->pdAtomsRadii[inp_surf->sele1.selatm[ii]-1];   
  }
  for(ii=0; ii<inp_surf->sele1.nselatm; ii++)
  {
    for(jj=0; jj<inp_surf->iNumOfRealAtoms; jj++)
    {
      
      if(inp_surf->piWithinAtoms[jj]==0)
      {
        fDist=sqrt(
              pow((trj_crd->xcoor[inp_surf->sele1.selatm[ii]-1]-trj_crd->xcoor[jj]), 2)+
              pow((trj_crd->ycoor[inp_surf->sele1.selatm[ii]-1]-trj_crd->ycoor[jj]), 2)+
              pow((trj_crd->zcoor[inp_surf->sele1.selatm[ii]-1]-trj_crd->zcoor[jj]), 2));
        if(fDist<=50)
        {
          iNumOfWithinAtoms++;
          inp_surf->piWithinAtoms[jj]=1;
          inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][0]=trj_crd->xcoor[jj];
          inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][1]=trj_crd->ycoor[jj];
          inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][2]=trj_crd->zcoor[jj];
          inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][3]=inp_surf->pdAtomsRadii[jj];
        }
      }
    }
  }

  inp_surf->iNumOfAtoms=iNumOfWithinAtoms+1;

  MakeNeighbors(inp_surf);
  
  iRotCont=0;
  while(NorthPoleTest(inp_surf)==0)
  {
    iRotCont++;
    MolRotation(inp_surf);    
  }

  for(ii=0; ii<inp_surf->sele1.nselatm; ii++)
  {
    dTmpSurf=GetSurf(ii, inp_surf);
    if(isnan(dTmpSurf)==1)
    {
      dTmpSurf=0.0;
    }
    dSurf=dSurf+dTmpSurf;
  }
  
  return dSurf;
}
// ============================================================================================================

// === GEPOL Functions ========================================================================================
double ComputeGEPOLSurf(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring) //FILE *OutputFile)
{
  int     ii;
  double  dSurf=0.0;
    
  // Set inp_surf->iNumOfAtoms to inp_surf->iNumOfRealAtoms
  // and copy frame coords to inp_surf->ppfAtomCoord
  inp_surf->iNumOfAtoms=inp_surf->iNumOfRealAtoms;
  for(ii=0; ii<inp_surf->iNumOfRealAtoms; ii++)
  {
    inp_surf->ppfAtomCoord[ii][0]=trj_crd->xcoor[ii];
    inp_surf->ppfAtomCoord[ii][1]=trj_crd->ycoor[ii];
    inp_surf->ppfAtomCoord[ii][2]=trj_crd->zcoor[ii];
  }

  // Set all atoms as ghost
  for(ii=0; ii<inp_surf->iMaxSize; ii++)
    inp_surf->piAtomType[ii]=3;
    
  // Set selected atoms as 1 or 6 in accordance with radii
  for(ii=0; ii<inp_surf->sele1.nselatm; ii++)
  {
    if(inp_surf->pdAtomsRadii[inp_surf->sele1.selatm[ii]-1]==0)
      inp_surf->piAtomType[inp_surf->sele1.selatm[ii]-1]=1;
    
    else
      inp_surf->piAtomType[inp_surf->sele1.selatm[ii]-1]=6;
  }

  if(inp_surf->sele1.nselatm<inp_surf->iNumOfRealAtoms)
    inp_surf->iGhostFlag=1;
  
  // === Surface Calculation ==========================
  if(inp_surf->iSurfType==2)
  {
    if(inp_surf->iGhostFlag==1)
      Shell(inp_surf);
      
    Bulk(inp_surf);
    Clean5(inp_surf);
    GenAtoms(inp_surf);
    Clean5(inp_surf);
  }
  
  if(inp_surf->iSurfType==1 && inp_surf->iSolvRadAddedFalg==0)
  {
    inp_surf->iSolvRadAddedFalg=1;
    AddSolvRad(inp_surf);
  }
  
  GeoCav(inp_surf);
  dSurf=AreaSum(inp_surf);
    
//printf("%8.3f", dSurf);
  return dSurf;
}

void FillJVTArrays(struct inp_surf *inp_surf)
{
  /*
   * This function fills ppiJVT1 & ppiJVT2 arrays
   */
  
  // === JVT1 =========== //
  inp_surf->ppiJVT1[0][0]=0;
  inp_surf->ppiJVT1[0][1]=0;
  inp_surf->ppiJVT1[0][2]=0;
  inp_surf->ppiJVT1[0][3]=0;
  inp_surf->ppiJVT1[0][4]=0;
  inp_surf->ppiJVT1[0][5]=6;
  inp_surf->ppiJVT1[0][6]=7;
  inp_surf->ppiJVT1[0][7]=8;
  inp_surf->ppiJVT1[0][8]=9;
  inp_surf->ppiJVT1[0][9]=10;
  inp_surf->ppiJVT1[0][10]=7;
  inp_surf->ppiJVT1[0][11]=8;
  inp_surf->ppiJVT1[0][12]=9;
  inp_surf->ppiJVT1[0][13]=10;
  inp_surf->ppiJVT1[0][14]=6;
  inp_surf->ppiJVT1[0][15]=6;
  inp_surf->ppiJVT1[0][16]=7;
  inp_surf->ppiJVT1[0][17]=8;
  inp_surf->ppiJVT1[0][18]=9;
  inp_surf->ppiJVT1[0][19]=10;
  inp_surf->ppiJVT1[0][20]=7;
  inp_surf->ppiJVT1[0][21]=8;
  inp_surf->ppiJVT1[0][22]=9;
  inp_surf->ppiJVT1[0][23]=10;
  inp_surf->ppiJVT1[0][24]=6;
  inp_surf->ppiJVT1[0][25]=6;
  inp_surf->ppiJVT1[0][26]=7;
  inp_surf->ppiJVT1[0][27]=8;
  inp_surf->ppiJVT1[0][28]=9;
  inp_surf->ppiJVT1[0][29]=10;
  inp_surf->ppiJVT1[0][30]=21;
  inp_surf->ppiJVT1[0][31]=22;
  inp_surf->ppiJVT1[0][32]=23;
  inp_surf->ppiJVT1[0][33]=24;
  inp_surf->ppiJVT1[0][34]=25;
  inp_surf->ppiJVT1[0][35]=21;
  inp_surf->ppiJVT1[0][36]=22;
  inp_surf->ppiJVT1[0][37]=23;
  inp_surf->ppiJVT1[0][38]=24;
  inp_surf->ppiJVT1[0][39]=25;
  inp_surf->ppiJVT1[0][40]=21;
  inp_surf->ppiJVT1[0][41]=22;
  inp_surf->ppiJVT1[0][42]=23;
  inp_surf->ppiJVT1[0][43]=24;
  inp_surf->ppiJVT1[0][44]=25;
  inp_surf->ppiJVT1[0][45]=21;
  inp_surf->ppiJVT1[0][46]=22;
  inp_surf->ppiJVT1[0][47]=23;
  inp_surf->ppiJVT1[0][48]=24;
  inp_surf->ppiJVT1[0][49]=25;
  inp_surf->ppiJVT1[0][50]=21;
  inp_surf->ppiJVT1[0][51]=22;
  inp_surf->ppiJVT1[0][52]=23;
  inp_surf->ppiJVT1[0][53]=24;
  inp_surf->ppiJVT1[0][54]=25;
  inp_surf->ppiJVT1[0][55]=31;
  inp_surf->ppiJVT1[0][56]=31;
  inp_surf->ppiJVT1[0][57]=31;
  inp_surf->ppiJVT1[0][58]=31;
  inp_surf->ppiJVT1[0][59]=31;
  inp_surf->ppiJVT1[1][0]=5;
  inp_surf->ppiJVT1[1][1]=1;
  inp_surf->ppiJVT1[1][2]=2;
  inp_surf->ppiJVT1[1][3]=3;
  inp_surf->ppiJVT1[1][4]=4;
  inp_surf->ppiJVT1[1][5]=1;
  inp_surf->ppiJVT1[1][6]=2;
  inp_surf->ppiJVT1[1][7]=3;
  inp_surf->ppiJVT1[1][8]=4;
  inp_surf->ppiJVT1[1][9]=5;
  inp_surf->ppiJVT1[1][10]=1;
  inp_surf->ppiJVT1[1][11]=2;
  inp_surf->ppiJVT1[1][12]=3;
  inp_surf->ppiJVT1[1][13]=4;
  inp_surf->ppiJVT1[1][14]=5;
  inp_surf->ppiJVT1[1][15]=11;
  inp_surf->ppiJVT1[1][16]=12;
  inp_surf->ppiJVT1[1][17]=13;
  inp_surf->ppiJVT1[1][18]=14;
  inp_surf->ppiJVT1[1][19]=15;
  inp_surf->ppiJVT1[1][20]=11;
  inp_surf->ppiJVT1[1][21]=12;
  inp_surf->ppiJVT1[1][22]=13;
  inp_surf->ppiJVT1[1][23]=14;
  inp_surf->ppiJVT1[1][24]=15;
  inp_surf->ppiJVT1[1][25]=16;
  inp_surf->ppiJVT1[1][26]=17;
  inp_surf->ppiJVT1[1][27]=18;
  inp_surf->ppiJVT1[1][28]=19;
  inp_surf->ppiJVT1[1][29]=20;
  inp_surf->ppiJVT1[1][30]=11;
  inp_surf->ppiJVT1[1][31]=12;
  inp_surf->ppiJVT1[1][32]=13;
  inp_surf->ppiJVT1[1][33]=14;
  inp_surf->ppiJVT1[1][34]=15;
  inp_surf->ppiJVT1[1][35]=17;
  inp_surf->ppiJVT1[1][36]=18;
  inp_surf->ppiJVT1[1][37]=19;
  inp_surf->ppiJVT1[1][38]=20;
  inp_surf->ppiJVT1[1][39]=16;
  inp_surf->ppiJVT1[1][40]=16;
  inp_surf->ppiJVT1[1][41]=17;
  inp_surf->ppiJVT1[1][42]=18;
  inp_surf->ppiJVT1[1][43]=19;
  inp_surf->ppiJVT1[1][44]=20;
  inp_surf->ppiJVT1[1][45]=27;
  inp_surf->ppiJVT1[1][46]=28;
  inp_surf->ppiJVT1[1][47]=29;
  inp_surf->ppiJVT1[1][48]=30;
  inp_surf->ppiJVT1[1][49]=26;
  inp_surf->ppiJVT1[1][50]=26;
  inp_surf->ppiJVT1[1][51]=27;
  inp_surf->ppiJVT1[1][52]=28;
  inp_surf->ppiJVT1[1][53]=29;
  inp_surf->ppiJVT1[1][54]=30;
  inp_surf->ppiJVT1[1][55]=27;
  inp_surf->ppiJVT1[1][56]=28;
  inp_surf->ppiJVT1[1][57]=29;
  inp_surf->ppiJVT1[1][58]=30;
  inp_surf->ppiJVT1[1][59]=26;
  inp_surf->ppiJVT1[2][0]=1;
  inp_surf->ppiJVT1[2][1]=2;
  inp_surf->ppiJVT1[2][2]=3;
  inp_surf->ppiJVT1[2][3]=4;
  inp_surf->ppiJVT1[2][4]=5;
  inp_surf->ppiJVT1[2][5]=5;
  inp_surf->ppiJVT1[2][6]=1;
  inp_surf->ppiJVT1[2][7]=2;
  inp_surf->ppiJVT1[2][8]=3;
  inp_surf->ppiJVT1[2][9]=4;
  inp_surf->ppiJVT1[2][10]=11;
  inp_surf->ppiJVT1[2][11]=12;
  inp_surf->ppiJVT1[2][12]=13;
  inp_surf->ppiJVT1[2][13]=14;
  inp_surf->ppiJVT1[2][14]=15;
  inp_surf->ppiJVT1[2][15]=1;
  inp_surf->ppiJVT1[2][16]=2;
  inp_surf->ppiJVT1[2][17]=3;
  inp_surf->ppiJVT1[2][18]=4;
  inp_surf->ppiJVT1[2][19]=5;
  inp_surf->ppiJVT1[2][20]=17;
  inp_surf->ppiJVT1[2][21]=18;
  inp_surf->ppiJVT1[2][22]=19;
  inp_surf->ppiJVT1[2][23]=20;
  inp_surf->ppiJVT1[2][24]=16;
  inp_surf->ppiJVT1[2][25]=11;
  inp_surf->ppiJVT1[2][26]=12;
  inp_surf->ppiJVT1[2][27]=13;
  inp_surf->ppiJVT1[2][28]=14;
  inp_surf->ppiJVT1[2][29]=15;
  inp_surf->ppiJVT1[2][30]=16;
  inp_surf->ppiJVT1[2][31]=17;
  inp_surf->ppiJVT1[2][32]=18;
  inp_surf->ppiJVT1[2][33]=19;
  inp_surf->ppiJVT1[2][34]=20;
  inp_surf->ppiJVT1[2][35]=11;
  inp_surf->ppiJVT1[2][36]=12;
  inp_surf->ppiJVT1[2][37]=13;
  inp_surf->ppiJVT1[2][38]=14;
  inp_surf->ppiJVT1[2][39]=15;
  inp_surf->ppiJVT1[2][40]=26;
  inp_surf->ppiJVT1[2][41]=27;
  inp_surf->ppiJVT1[2][42]=28;
  inp_surf->ppiJVT1[2][43]=29;
  inp_surf->ppiJVT1[2][44]=30;
  inp_surf->ppiJVT1[2][45]=17;
  inp_surf->ppiJVT1[2][46]=18;
  inp_surf->ppiJVT1[2][47]=19;
  inp_surf->ppiJVT1[2][48]=20;
  inp_surf->ppiJVT1[2][49]=16;
  inp_surf->ppiJVT1[2][50]=27;
  inp_surf->ppiJVT1[2][51]=28;
  inp_surf->ppiJVT1[2][52]=29;
  inp_surf->ppiJVT1[2][53]=30;
  inp_surf->ppiJVT1[2][54]=26;
  inp_surf->ppiJVT1[2][55]=26;
  inp_surf->ppiJVT1[2][56]=27;
  inp_surf->ppiJVT1[2][57]=28;
  inp_surf->ppiJVT1[2][58]=29;
  inp_surf->ppiJVT1[2][59]=30;

  // === JVT2 =========== //
  inp_surf->ppiJVT2[0][0]=0;
  inp_surf->ppiJVT2[0][1]=4;
  inp_surf->ppiJVT2[0][2]=3;
  inp_surf->ppiJVT2[0][3]=5;
  inp_surf->ppiJVT2[1][0]=4;
  inp_surf->ppiJVT2[1][1]=1;
  inp_surf->ppiJVT2[1][2]=5;
  inp_surf->ppiJVT2[1][3]=3;
  inp_surf->ppiJVT2[2][0]=3;
  inp_surf->ppiJVT2[2][1]=5;
  inp_surf->ppiJVT2[2][2]=2;
  inp_surf->ppiJVT2[2][3]=4;
}

void CompTriaVertexCoord(struct inp_surf *inp_surf)
{
  /*
   * This computes the triangle vertex coordinates for a sphere of radius
   * one, projecting the pentakisdodecahedro onto it.
   * 
   * original name: TES
   */
  
  int     ii, jj;
  int     iCont=0;
  double  dTmpFiv=0.0;
                                         
  double  pdThev[6]={0.65235813978436820,1.10714871779409050,\
                     1.38208579601133450,1.75950685757845870,\
                     2.03444393579570270,2.48923451380542510};
                    
                    
  double  pdFiv[6]={0.62831853071795860,0.0,\
                    0.62831853071795860,0.0,\
                    0.62831853071795860,0.0};
  
  double  dFir=1.25663706143591730;

  inp_surf->ppdCV[0][0] =0.0;
  inp_surf->ppdCV[0][1] =0.0;
  inp_surf->ppdCV[0][2] =1.0;
  inp_surf->ppdCV[31][0]=0.0;
  inp_surf->ppdCV[31][1]=0.0;
  inp_surf->ppdCV[31][2]=-1.0;
  
  iCont=0;
  for(ii=0; ii<6; ii++)
  {
    dTmpFiv=pdFiv[ii];
    for(jj=0; jj<5; jj++)
    {
      dTmpFiv=dTmpFiv+dFir;
      if(jj==0)
      {
        dTmpFiv=pdFiv[ii];
      }
      iCont++;
      inp_surf->ppdCV[iCont][0]=sin(pdThev[ii])*cos(dTmpFiv);
      inp_surf->ppdCV[iCont][1]=sin(pdThev[ii])*sin(dTmpFiv);
      inp_surf->ppdCV[iCont][2]=cos(pdThev[ii]);
    }
  }
}

void Divide(struct inp_surf *inp_surf)
{
  /*
   * This divides the initial 60 spherical
   * triangles to the level indicated by NDIV
   * 
   * original name: DIVIDE
   */
  
  int     ii, jj, kk, mm, nn;
  int     iCont;
  int     iTmpJVT1_1=0, iTmpJVT1_2=0, iTmpJVT1_3=0;

  int     iTmpJVT2_21=0, iTmpJVT2_22=0, iTmpJVT2_23=0;
  int     iTmpJVT2_31=0, iTmpJVT2_32=0, iTmpJVT2_33=0;
  int     iTmpJVT2_41=0, iTmpJVT2_42=0, iTmpJVT2_43=0;
  int     iTmpJVT2_51=0, iTmpJVT2_52=0, iTmpJVT2_53=0;
  
  double  ppdTmpCV2[6][3];
  double  ppdTmpCV3[6][3];
  double  ppdTmpCV4[6][3];
  double  ppdTmpCV5[6][3];
  
  double  dCoord1X, dCoord1Y, dCoord1Z;
  double  dCoord2X, dCoord2Y, dCoord2Z;
  double  dCoord3X, dCoord3Y, dCoord3Z;

  iCont=-1;
  for(ii=0; ii<60; ii++)
  {
    iTmpJVT1_1=inp_surf->ppiJVT1[0][ii];
    iTmpJVT1_2=inp_surf->ppiJVT1[1][ii];
    iTmpJVT1_3=inp_surf->ppiJVT1[2][ii];
    
    dCoord1X=inp_surf->ppdCV[iTmpJVT1_1][0];
    dCoord1Y=inp_surf->ppdCV[iTmpJVT1_1][1];
    dCoord1Z=inp_surf->ppdCV[iTmpJVT1_1][2];

    dCoord2X=inp_surf->ppdCV[iTmpJVT1_2][0];
    dCoord2Y=inp_surf->ppdCV[iTmpJVT1_2][1];
    dCoord2Z=inp_surf->ppdCV[iTmpJVT1_2][2];
    
    dCoord3X=inp_surf->ppdCV[iTmpJVT1_3][0];
    dCoord3Y=inp_surf->ppdCV[iTmpJVT1_3][1];
    dCoord3Z=inp_surf->ppdCV[iTmpJVT1_3][2];
        
    if(inp_surf->iNDIV==1)
    {
      CalcTriaCenter(inp_surf, dCoord1X, dCoord1Y, dCoord1Z, dCoord2X, dCoord2Y, dCoord2Z, dCoord3X, dCoord3Y, dCoord3Z);
      iCont++;
      inp_surf->pfXCoord1[iCont]=inp_surf->pdTriaCenter[0];
      inp_surf->pfYCoord1[iCont]=inp_surf->pdTriaCenter[1];
      inp_surf->pfZCoord1[iCont]=inp_surf->pdTriaCenter[2];
      continue;
    }
    ppdTmpCV2[0][0]=dCoord1X;
    ppdTmpCV2[0][1]=dCoord1Y;
    ppdTmpCV2[0][2]=dCoord1Z;
    
    ppdTmpCV2[1][0]=dCoord2X;
    ppdTmpCV2[1][1]=dCoord2Y;
    ppdTmpCV2[1][2]=dCoord2Z;
    
    ppdTmpCV2[2][0]=dCoord3X;
    ppdTmpCV2[2][1]=dCoord3Y;
    ppdTmpCV2[2][2]=dCoord3Z;
    
    TriaToFour(ppdTmpCV2);
    
    for(jj=0; jj<4; jj++)
    {
      iTmpJVT2_21=inp_surf->ppiJVT2[0][jj];
      dCoord1X=ppdTmpCV2[iTmpJVT2_21][0];
      dCoord1Y=ppdTmpCV2[iTmpJVT2_21][1];
      dCoord1Z=ppdTmpCV2[iTmpJVT2_21][2];
      
      iTmpJVT2_22=inp_surf->ppiJVT2[1][jj];
      dCoord2X=ppdTmpCV2[iTmpJVT2_22][0];
      dCoord2Y=ppdTmpCV2[iTmpJVT2_22][1];
      dCoord2Z=ppdTmpCV2[iTmpJVT2_22][2];
      
      iTmpJVT2_23=inp_surf->ppiJVT2[2][jj];
      dCoord3X=ppdTmpCV2[iTmpJVT2_23][0];
      dCoord3Y=ppdTmpCV2[iTmpJVT2_23][1];
      dCoord3Z=ppdTmpCV2[iTmpJVT2_23][2];
      
      if(inp_surf->iNDIV==2)
      {
        CalcTriaCenter(inp_surf, dCoord1X, dCoord1Y, dCoord1Z, dCoord2X, dCoord2Y, dCoord2Z, dCoord3X, dCoord3Y, dCoord3Z);
        iCont++;
        inp_surf->pfXCoord1[iCont]=inp_surf->pdTriaCenter[0];
        inp_surf->pfYCoord1[iCont]=inp_surf->pdTriaCenter[1];
        inp_surf->pfZCoord1[iCont]=inp_surf->pdTriaCenter[2];
        continue;
      }

      ppdTmpCV3[0][0]=dCoord1X;
      ppdTmpCV3[0][1]=dCoord1Y;
      ppdTmpCV3[0][2]=dCoord1Z;
      
      ppdTmpCV3[1][0]=dCoord2X;
      ppdTmpCV3[1][1]=dCoord2Y;
      ppdTmpCV3[1][2]=dCoord2Z;
      
      ppdTmpCV3[2][0]=dCoord3X;
      ppdTmpCV3[2][1]=dCoord3Y;
      ppdTmpCV3[2][2]=dCoord3Z;
      
      TriaToFour(ppdTmpCV3);
      
      for(kk=0; kk<4; kk++)
      {
        iTmpJVT2_31=inp_surf->ppiJVT2[0][kk];
        dCoord1X=ppdTmpCV3[iTmpJVT2_31][0];
        dCoord1Y=ppdTmpCV3[iTmpJVT2_31][1];
        dCoord1Z=ppdTmpCV3[iTmpJVT2_31][2];
        
        iTmpJVT2_32=inp_surf->ppiJVT2[1][kk];
        dCoord2X=ppdTmpCV3[iTmpJVT2_32][0];
        dCoord2Y=ppdTmpCV3[iTmpJVT2_32][1];
        dCoord2Z=ppdTmpCV3[iTmpJVT2_32][2];
        
        iTmpJVT2_33=inp_surf->ppiJVT2[2][kk];
        dCoord3X=ppdTmpCV3[iTmpJVT2_33][0];
        dCoord3Y=ppdTmpCV3[iTmpJVT2_33][1];
        dCoord3Z=ppdTmpCV3[iTmpJVT2_33][2];
        
        if(inp_surf->iNDIV==3)
        {
          CalcTriaCenter(inp_surf, dCoord1X, dCoord1Y, dCoord1Z, dCoord2X, dCoord2Y, dCoord2Z, dCoord3X, dCoord3Y, dCoord3Z);
          iCont++;
          inp_surf->pfXCoord1[iCont]=inp_surf->pdTriaCenter[0];
          inp_surf->pfYCoord1[iCont]=inp_surf->pdTriaCenter[1];
          inp_surf->pfZCoord1[iCont]=inp_surf->pdTriaCenter[2];
          continue;
        }

        ppdTmpCV4[0][0]=dCoord1X;
        ppdTmpCV4[0][1]=dCoord1Y;
        ppdTmpCV4[0][2]=dCoord1Z;
        
        ppdTmpCV4[1][0]=dCoord2X;
        ppdTmpCV4[1][1]=dCoord2Y;
        ppdTmpCV4[1][2]=dCoord2Z;
        
        ppdTmpCV4[2][0]=dCoord3X;
        ppdTmpCV4[2][1]=dCoord3Y;
        ppdTmpCV4[2][2]=dCoord3Z;

        TriaToFour(ppdTmpCV4);
        
        for(mm=0; mm<4; mm++)
        {
          iTmpJVT2_41=inp_surf->ppiJVT2[0][mm];
          dCoord1X=ppdTmpCV4[iTmpJVT2_41][0];
          dCoord1Y=ppdTmpCV4[iTmpJVT2_41][1];
          dCoord1Z=ppdTmpCV4[iTmpJVT2_41][2];
          
          iTmpJVT2_42=inp_surf->ppiJVT2[1][mm];
          dCoord2X=ppdTmpCV4[iTmpJVT2_42][0];
          dCoord2Y=ppdTmpCV4[iTmpJVT2_42][1];
          dCoord2Z=ppdTmpCV4[iTmpJVT2_42][2];
          
          iTmpJVT2_43=inp_surf->ppiJVT2[2][mm];
          dCoord3X=ppdTmpCV4[iTmpJVT2_43][0];
          dCoord3Y=ppdTmpCV4[iTmpJVT2_43][1];
          dCoord3Z=ppdTmpCV4[iTmpJVT2_43][2];
          
          if(inp_surf->iNDIV==4)
          {
            CalcTriaCenter(inp_surf, dCoord1X, dCoord1Y, dCoord1Z, dCoord2X, dCoord2Y, dCoord2Z, dCoord3X, dCoord3Y, dCoord3Z);
            iCont++;
            inp_surf->pfXCoord1[iCont]=inp_surf->pdTriaCenter[0];
            inp_surf->pfYCoord1[iCont]=inp_surf->pdTriaCenter[1];
            inp_surf->pfZCoord1[iCont]=inp_surf->pdTriaCenter[2];
            continue;
          }

          ppdTmpCV5[0][0]=dCoord1X;
          ppdTmpCV5[0][1]=dCoord1Y;
          ppdTmpCV5[0][2]=dCoord1Z;
          
          ppdTmpCV5[1][0]=dCoord2X;
          ppdTmpCV5[1][1]=dCoord2Y;
          ppdTmpCV5[1][2]=dCoord2Z;
          
          ppdTmpCV5[2][0]=dCoord3X;
          ppdTmpCV5[2][1]=dCoord3Y;
          ppdTmpCV5[2][2]=dCoord3Z;

          TriaToFour(ppdTmpCV5);
          
          for(nn=0; nn<4; nn++)
          {
            iTmpJVT2_51=inp_surf->ppiJVT2[0][nn];
            dCoord1X=ppdTmpCV5[iTmpJVT2_51][0];
            dCoord1Y=ppdTmpCV5[iTmpJVT2_51][1];
            dCoord1Z=ppdTmpCV5[iTmpJVT2_51][2];
            
            iTmpJVT2_52=inp_surf->ppiJVT2[1][nn];
            dCoord2X=ppdTmpCV5[iTmpJVT2_52][0];
            dCoord2Y=ppdTmpCV5[iTmpJVT2_52][1];
            dCoord2Z=ppdTmpCV5[iTmpJVT2_52][2];
            
            iTmpJVT2_53=inp_surf->ppiJVT2[2][nn];
            dCoord3X=ppdTmpCV5[iTmpJVT2_53][0];
            dCoord3Y=ppdTmpCV5[iTmpJVT2_53][1];
            dCoord3Z=ppdTmpCV5[iTmpJVT2_53][2];
            
            CalcTriaCenter(inp_surf, dCoord1X, dCoord1Y, dCoord1Z, dCoord2X, dCoord2Y, dCoord2Z, dCoord3X, dCoord3Y, dCoord3Z);
            iCont++;
            inp_surf->pfXCoord1[iCont]=inp_surf->pdTriaCenter[0];
            inp_surf->pfYCoord1[iCont]=inp_surf->pdTriaCenter[1];
            inp_surf->pfZCoord1[iCont]=inp_surf->pdTriaCenter[2];
          }
        }
      }
    }
  }
}

void CalcTriaCenter(struct inp_surf *inp_surf,
                    double dCoord1X, double dCoord1Y, double dCoord1Z,
                    double dCoord2X, double dCoord2Y, double dCoord2Z,
                    double dCoord3X, double dCoord3Y, double dCoord3Z)
{
  /*
   * This computes the center of a spherical triangle
   * 
   * original name: CALCEN
   */
   
   double dXCoord, dYCoord, dZCoord, dRadius;
   
   dXCoord=(dCoord1X+dCoord2X+dCoord3X)/3.0;
   dYCoord=(dCoord1Y+dCoord2Y+dCoord3Y)/3.0;
   dZCoord=(dCoord1Z+dCoord2Z+dCoord3Z)/3.0;
   dRadius=1.0/sqrt(dXCoord*dXCoord+dYCoord*dYCoord+dZCoord*dZCoord);
      
   inp_surf->pdTriaCenter[0]=dXCoord*dRadius;
   inp_surf->pdTriaCenter[1]=dYCoord*dRadius;
   inp_surf->pdTriaCenter[2]=dZCoord*dRadius;
}

void TriaToFour(double ppdTmpCV[6][3])
{
  /*
   * This divides one triangle into four
   * 
   * original name: CALVER
   */
  
  int     ii;
  int     iCont1=0, iCont2=0;
  
  double  dTmpX, dTmpY, dTmpZ, dTmpR;
  
  for(ii=0; ii<3; ii++)
  {
    iCont1=ii-1;
    iCont2=ii+3;
    if(ii==0)
    {
      iCont1=2;
    }
    dTmpX=((ppdTmpCV[ii][0]+ppdTmpCV[iCont1][0])/2.0);
    dTmpY=((ppdTmpCV[ii][1]+ppdTmpCV[iCont1][1])/2.0);
    dTmpZ=((ppdTmpCV[ii][2]+ppdTmpCV[iCont1][2])/2.0);
    dTmpR=1.0/sqrt(dTmpX*dTmpX+dTmpY*dTmpY+dTmpZ*dTmpZ);
    ppdTmpCV[iCont2][0]=dTmpX*dTmpR;
    ppdTmpCV[iCont2][1]=dTmpY*dTmpR;
    ppdTmpCV[iCont2][2]=dTmpZ*dTmpR;
  }
}

void Shell(struct inp_surf *inp_surf)
{
  /*
   * This determines which ghost spheres are around the real ones.
   * 
   * original name: SHELL
   */
  
    int       ii, jj;
    float     fDist=0.0, fTest=0.0;
    
    for(ii=0; ii<inp_surf->iNumOfAtoms-1; ii++)
    {      
      for(jj=ii+1; jj<inp_surf->iNumOfAtoms; jj++)
      {
        if((inp_surf->piAtomType[ii]==3) && (inp_surf->piAtomType[jj]==6))
        {
          fDist=(pow((inp_surf->ppfAtomCoord[ii][0]-inp_surf->ppfAtomCoord[jj][0]), 2)+
                 pow((inp_surf->ppfAtomCoord[ii][1]-inp_surf->ppfAtomCoord[jj][1]), 2)+
                 pow((inp_surf->ppfAtomCoord[ii][2]-inp_surf->ppfAtomCoord[jj][2]), 2));
          
          fTest=pow((inp_surf->pdAtomsRadii[ii]+inp_surf->pdAtomsRadii[jj]+2*inp_surf->dSolventRadius), 2);
          if(fDist<fTest)
          {
            inp_surf->piAtomType[ii]=5;
          }
        }
        
        else if((inp_surf->piAtomType[ii]==6) && (inp_surf->piAtomType[jj]==3))
        {
          fDist=(pow((inp_surf->ppfAtomCoord[ii][0]-inp_surf->ppfAtomCoord[jj][0]), 2)+
                 pow((inp_surf->ppfAtomCoord[ii][1]-inp_surf->ppfAtomCoord[jj][1]), 2)+
                 pow((inp_surf->ppfAtomCoord[ii][2]-inp_surf->ppfAtomCoord[jj][2]), 2));
          
          fTest=pow((inp_surf->pdAtomsRadii[ii]+inp_surf->pdAtomsRadii[jj]+2*inp_surf->dSolventRadius), 2);
          if(fDist<fTest)
          {
            inp_surf->piAtomType[jj]=5;
          }
        }
      }
    }
}

void Bulk(struct inp_surf *inp_surf)
{
  /*
   * This subroutine creates new spheres.
   * 
   * original name: BULK
   */
  
  int     ii=0, jj=0, kk=0, mm=0, nn=0, oo=0, ww=0;
  
  int     iNumOfGen=0, iNL1=0;
  int     iCont=0,  iCont2=0, iCont3=0; 
  int     iCont4=0, iCont5=0, iCont6=0;
  int     iNumOfStartingAtoms=0;
  
  int     iExitFlag=0;
  
  float   fPairDist1=0.0, fDistTest1=0.0;
  float   fPairDist2=0.0;
  float   fPairDist3=0.0, fDistTest3=0.0;
  float   fDistTest4=0.0;
  float   fPairDist5=0.0, fDistTest5=0.0;
  float   fPairDist6=0.0, fDistTest6=0.0;

  float   fNewCenter=0.0, fTmpOFAC=0.0;
  float   fNewXCoord=0.0, fTmpRadius1=0.0;
  float   fNewYCoord=0.0, fTmpRadius2=0.0;
  float   fNewZCoord=0.0, fTmpRadius3=0.0;
  
  float   fTmpXCoord=0.0, fNewRadius1=0.0;
  float   fTmpYCoord=0.0, fNewRadius2=0.0;
  float   fTmpZCoord=0.0;
  
  fTmpOFAC=1.0-inp_surf->fOFAC*2.0;
  iNL1=1;
  iNumOfStartingAtoms=inp_surf->iNumOfAtoms;
  
  while(iExitFlag!=1)
  {
    //iBreakFlag=0;
    
    // Loop to select the first atom of the pair
    for(ii=iNL1; ii<iNumOfStartingAtoms; ii++)
    {
      if(inp_surf->piAtomType[ii]<=4)
        continue;
      
      fTmpRadius1=inp_surf->pdAtomsRadii[ii];
      
      for(jj=0; jj<inp_surf->iNumOfAtoms; jj++)
      {
        inp_surf->pfDistVect[jj]=sqrtf(((inp_surf->ppfAtomCoord[ii][0]-inp_surf->ppfAtomCoord[jj][0])*
                                        (inp_surf->ppfAtomCoord[ii][0]-inp_surf->ppfAtomCoord[jj][0]))+
                                       ((inp_surf->ppfAtomCoord[ii][1]-inp_surf->ppfAtomCoord[jj][1])*
                                        (inp_surf->ppfAtomCoord[ii][1]-inp_surf->ppfAtomCoord[jj][1]))+
                                       ((inp_surf->ppfAtomCoord[ii][2]-inp_surf->ppfAtomCoord[jj][2])*
                                        (inp_surf->ppfAtomCoord[ii][2]-inp_surf->ppfAtomCoord[jj][2])));
      }
           
      // Loop to select the second atom of the pair
      for(jj=0; jj<ii; jj++)
      {
        if(inp_surf->piAtomType[jj]<=4) continue;

        fTmpRadius2=inp_surf->pdAtomsRadii[jj];
        fPairDist1=inp_surf->pfDistVect[jj];
        
        // If the solvent can pass through the pair this pair is discarded
        if(fPairDist1>=(inp_surf->dSolventRadius+inp_surf->dSolventRadius+fTmpRadius1+fTmpRadius2)) continue;

        // Test of overlapping and assignment of the radius of the new atom
        if(fTmpRadius1>fTmpRadius2)
        {
          fDistTest1=fTmpRadius1+fTmpRadius2*fTmpOFAC;
          if(fPairDist1<fDistTest1) continue;
          fNewRadius1=fTmpRadius1;
        }
        else
        {
          fDistTest1=fTmpRadius2+fTmpRadius1*fTmpOFAC;
          if(fPairDist1<fDistTest1) continue;
          fNewRadius1=fTmpRadius2;
        }
        
        // Computes coordinates of the new atom
        fNewCenter=((fPairDist1-fTmpRadius2+fTmpRadius1)/(fPairDist1+fTmpRadius2-fTmpRadius1));
        fNewXCoord=(inp_surf->ppfAtomCoord[ii][0]+fNewCenter*inp_surf->ppfAtomCoord[jj][0])/(fNewCenter+1.0);
        fNewYCoord=(inp_surf->ppfAtomCoord[ii][1]+fNewCenter*inp_surf->ppfAtomCoord[jj][1])/(fNewCenter+1.0);
        fNewZCoord=(inp_surf->ppfAtomCoord[ii][2]+fNewCenter*inp_surf->ppfAtomCoord[jj][2])/(fNewCenter+1.0);

        fNewRadius2=(fPairDist1-fTmpRadius2+fTmpRadius1)*0.5;
        
        // Test of overlapping for the new atom
        iCont=-1;
        for(kk=0; kk<inp_surf->iNumOfAtoms; kk++)
        {
          //iBreakFlag=0;
          fPairDist2=inp_surf->pfDistVect[kk];
          if(fPairDist2>=(fNewRadius2+fNewRadius1+inp_surf->pdAtomsRadii[kk])) continue;
          
          if(inp_surf->piAtomType[kk]<=3) continue;
          
          fPairDist3=(pow((fNewXCoord-inp_surf->ppfAtomCoord[kk][0]), 2)+
                      pow((fNewYCoord-inp_surf->ppfAtomCoord[kk][1]), 2)+
                      pow((fNewZCoord-inp_surf->ppfAtomCoord[kk][2]), 2));
          
          fTmpRadius3=inp_surf->pdAtomsRadii[kk];
          
          if(fPairDist3>=pow((fTmpRadius3+fNewRadius1), 2)) continue;
          
          fDistTest3=((fNewRadius1-fTmpRadius3)*(fNewRadius1-fTmpRadius3));
          
          if(fPairDist3<=fDistTest3)
          {
            if(fNewRadius1>fTmpRadius3)
            {
              if(fDistTest3<0.0004) goto L603;
              iCont++;
              inp_surf->piEngulfedAtoms[iCont]=kk;
              goto L604;
            }
            else goto L603;
          }
          
          fDistTest4=fTmpRadius3+fNewRadius1*fTmpOFAC;
          if(fDistTest4<0.0) goto L604;
          
          fDistTest4=fDistTest4*fDistTest4;
          if(fPairDist3<=fDistTest4) goto L603;

          L604:;
        } // kk, K, 604
        
        // Find atoms that are overlapped with the new atom
        // Use radius of the atom plus inp_surf->dSolventRadius
        
        iCont2=-1;
        for(mm=0; mm<inp_surf->iNumOfRealAtoms; mm++)      
        {
          if(inp_surf->piAtomType[mm]>=2)
          {
            fPairDist5=(pow((fNewXCoord-inp_surf->ppfAtomCoord[mm][0]), 2)+
                        pow((fNewYCoord-inp_surf->ppfAtomCoord[mm][1]), 2)+
                        pow((fNewZCoord-inp_surf->ppfAtomCoord[mm][2]), 2));
                        
            fDistTest5=(
                        (fNewRadius1+
                         inp_surf->dSolventRadius+
                         inp_surf->dSolventRadius+
                         inp_surf->pdAtomsRadii[mm])*                        
                        (fNewRadius1+
                         inp_surf->dSolventRadius+
                         inp_surf->dSolventRadius+
                         inp_surf->pdAtomsRadii[mm]));
            
            if(fPairDist5<fDistTest5)
            {
              iCont2++;
              inp_surf->piZeroSurf[iCont2]=mm;
            }
          }
        } // mm, K
        
        // Determine if the new atom has accessible surface area == zero
        //iBreakFlag=0;
        iCont3=0;        
        for(nn=0; nn<inp_surf->iNumOfTesserae; nn++)
        {        
          fTmpXCoord=inp_surf->pfXCoord1[nn]*(fNewRadius1+inp_surf->dSolventRadius)+fNewXCoord;
          fTmpYCoord=inp_surf->pfYCoord1[nn]*(fNewRadius1+inp_surf->dSolventRadius)+fNewYCoord;
          fTmpZCoord=inp_surf->pfZCoord1[nn]*(fNewRadius1+inp_surf->dSolventRadius)+fNewZCoord;

          iCont4=iCont3;          
          for(oo=iCont4; oo<iCont2+1; oo++)
          {
            iCont3=oo;
            iCont5=inp_surf->piZeroSurf[oo];
            fPairDist6=(pow((fTmpXCoord-inp_surf->ppfAtomCoord[iCont5][0]), 2)+
                        pow((fTmpYCoord-inp_surf->ppfAtomCoord[iCont5][1]), 2)+
                        pow((fTmpZCoord-inp_surf->ppfAtomCoord[iCont5][2]), 2));
                      
            fDistTest6=pow((inp_surf->pdAtomsRadii[iCont5]+inp_surf->dSolventRadius), 2);
          
            if(fPairDist6<fDistTest6) goto L3;
          }
          
          for(ww=0; ww<iCont4; ww++)
          {
            iCont3=ww;
            iCont6=inp_surf->piZeroSurf[ww];
            fPairDist6=(pow((fTmpXCoord-inp_surf->ppfAtomCoord[iCont6][0]), 2)+
                        pow((fTmpYCoord-inp_surf->ppfAtomCoord[iCont6][1]), 2)+
                        pow((fTmpZCoord-inp_surf->ppfAtomCoord[iCont6][2]), 2));
            
            fDistTest6=pow((inp_surf->pdAtomsRadii[iCont6]+inp_surf->dSolventRadius), 2);
            if(fPairDist6<fDistTest6) goto L3;
          }
          
          goto L603;
          
          L3:;
        } // nn, L, 3
                
        // Marks atoms that are engulfed by the new atom
        for(ww=0; ww<iCont+1; ww++)
        {
          kk=inp_surf->piEngulfedAtoms[ww];
          inp_surf->piAtomType[kk]=2;
        }
        
        // Save information about the new atom
        inp_surf->iNumOfAtoms++;
        if(inp_surf->iNumOfAtoms>=inp_surf->iMaxSize)
        {
          printf("SURF Memory Error\n");
          /*printf("Too many atoms created -> %d\n", inp_surf->iNumOfAtoms);*/
          printf("Increment --MEMSIZE (now %d)\n",  inp_surf->memsize);
          printf("Exiting now\n");
          exit(1);
        }
        
        inp_surf->ppfAtomCoord[inp_surf->iNumOfAtoms-1][0]=fNewXCoord;
        inp_surf->ppfAtomCoord[inp_surf->iNumOfAtoms-1][1]=fNewYCoord;
        inp_surf->ppfAtomCoord[inp_surf->iNumOfAtoms-1][2]=fNewZCoord;
        inp_surf->pdAtomsRadii[inp_surf->iNumOfAtoms-1]=fNewRadius1;
        inp_surf->piAtomType[inp_surf->iNumOfAtoms-1]=6;
        inp_surf->pfDistVect[inp_surf->iNumOfAtoms-1]=fNewRadius2;
        
        L603:;
      } // jj, J, 603
      
      //L602:;
    } // ii, I, 602
    
    // Compute number of generations
    iNumOfGen++;

    
    // Check if there are new atoms
    if(inp_surf->iNumOfAtoms!=iNumOfStartingAtoms)
    {
      // Mark spheres with area zero
      iNumOfStartingAtoms=MZero5(inp_surf, iNumOfStartingAtoms);
      if(inp_surf->iNumOfAtoms!=iNumOfStartingAtoms)
      {
        iNL1=iNumOfStartingAtoms;
        iNumOfStartingAtoms=inp_surf->iNumOfAtoms;
        goto L600;
      }
      else iExitFlag=1;
    }
    else iExitFlag=1;
    
    L600:;
  } // while, 600
  
}

int MZero5(struct inp_surf *inp_surf, int iNumOfStartingAtoms)
{
  /*
   * This discards atoms (piAtomType=2) engulfed by another
   * Marks (piAtomType=4) the atoms with total area zero.
   * 
   * original name: MZERO5
   */
  
  int     ii, jj, kk;
  int     iIndex1=0, iIndex2=0, iIndex3=0;
  int     iIndex4=0, iIndex5=0, iIndex6=0;
  int     iBreakFlag=0;
  
  float   fTmpXCoord=0.0,  fTmpYCoord=0.0,  fTmpZCoord=0.0,  fTmpRadius=0.0;
  float   fTmpXCoord2=0.0, fTmpYCoord2=0.0, fTmpZCoord2=0.0;
  float   fDist=0.0, fDistTest=0.0;

  // Discard atoms that have been engulfed by another (piAtomType=2)
  iIndex1=0;
  iIndex2=0;
  for(ii=inp_surf->iNumOfRealAtoms-1; ii<inp_surf->iNumOfAtoms; ii++)
  {
    iIndex3=ii-iIndex1;
    inp_surf->ppfAtomCoord[iIndex3][0]=inp_surf->ppfAtomCoord[ii][0];
    inp_surf->ppfAtomCoord[iIndex3][1]=inp_surf->ppfAtomCoord[ii][1];
    inp_surf->ppfAtomCoord[iIndex3][2]=inp_surf->ppfAtomCoord[ii][2];
    inp_surf->pdAtomsRadii[iIndex3]=inp_surf->pdAtomsRadii[ii];
    inp_surf->piAtomType[iIndex3]=inp_surf->piAtomType[ii];
    
    if(inp_surf->piAtomType[ii]==2)
    {
      iIndex1++;
      if(ii<=iNumOfStartingAtoms)
      {
        iIndex2++;
      }
    }
  }
  
  //iNew=inp_surf->iNumOfAtoms-iNumOfStartingAtoms;
  inp_surf->iNumOfAtoms=inp_surf->iNumOfAtoms-iIndex1;
  iNumOfStartingAtoms=iNumOfStartingAtoms-iIndex2;

  // Find atoms with zero area
  for(ii=0; ii<iNumOfStartingAtoms; ii++)
  {
    inp_surf->piZeroCheck[ii]=0;
  }

  for(ii=iNumOfStartingAtoms-1; ii<inp_surf->iNumOfAtoms; ii++)
  {
    inp_surf->piZeroCheck[ii]=0;
    if(inp_surf->piAtomType[ii]>4)
    {
      inp_surf->piZeroCheck[ii]=1;
    }
  }
  
  // Start
  for(ii=inp_surf->iNumOfAtoms-1; ii>=0; ii--)
  {
    if(inp_surf->piZeroCheck[ii]==1)
    {
      fTmpXCoord=inp_surf->ppfAtomCoord[ii][0];
      fTmpYCoord=inp_surf->ppfAtomCoord[ii][1];
      fTmpZCoord=inp_surf->ppfAtomCoord[ii][2];
      fTmpRadius=inp_surf->pdAtomsRadii[ii];
      
      // Find atoms that overlap atom ii
      iIndex4=-1;
      for(jj=0; jj<inp_surf->iNumOfAtoms; jj++)
      {
        if(inp_surf->piAtomType[jj]>=3)
        {
          fDist=(pow((fTmpXCoord-inp_surf->ppfAtomCoord[jj][0]), 2)+
                 pow((fTmpYCoord-inp_surf->ppfAtomCoord[jj][1]), 2)+
                 pow((fTmpZCoord-inp_surf->ppfAtomCoord[jj][2]), 2));
          
          fDistTest=pow((fTmpRadius+inp_surf->pdAtomsRadii[jj]), 2);
          if(fDist<fDistTest)
          {
            if(ii!=jj)
            {
              iIndex4++;
              inp_surf->piZeroSurfMZ5[iIndex4]=jj;
            }
          }
        }
      } // jj, J, 7
      
      // Mark atoms that should be checked
      if(ii>iNumOfStartingAtoms-1)
      {
        for(kk=0; kk<iIndex4+1; kk++)
        {
          iIndex1=inp_surf->piZeroSurfMZ5[kk];
          if(inp_surf->piAtomType[iIndex1]>4)
          {
            inp_surf->piZeroCheck[iIndex1]=1;
          }
        }
      }
      
      // Determine if the atom ii has area zero
      iBreakFlag=0;
      iIndex5=0;
      for(jj=0; jj<inp_surf->iNumOfTesserae; jj++)
      {
        fTmpXCoord2=inp_surf->pfXCoord1[jj]*fTmpRadius+fTmpXCoord;
        fTmpYCoord2=inp_surf->pfYCoord1[jj]*fTmpRadius+fTmpYCoord;
        fTmpZCoord2=inp_surf->pfZCoord1[jj]*fTmpRadius+fTmpZCoord;
        
        iIndex6=iIndex5;
        for(kk=iIndex6; kk<iIndex4+1; kk++)
        {
          iIndex5=kk;
          iIndex1=inp_surf->piZeroSurfMZ5[kk];
          fDist=(pow((fTmpXCoord2-inp_surf->ppfAtomCoord[iIndex1][0]), 2)+
                 pow((fTmpYCoord2-inp_surf->ppfAtomCoord[iIndex1][1]), 2)+
                 pow((fTmpZCoord2-inp_surf->ppfAtomCoord[iIndex1][2]), 2));
          
          fDistTest=(inp_surf->pdAtomsRadii[iIndex1]*inp_surf->pdAtomsRadii[iIndex1]);
          
          if(fDist<fDistTest)
          {
            iBreakFlag=1;
            break;
          }
        } // kk, LL
        
        if(iBreakFlag==1)
        {
          iBreakFlag=0;
          continue;
        }
        
        for(kk=0; kk<iIndex6; kk++)
        {
          iIndex5=kk;
          iIndex1=inp_surf->piZeroSurfMZ5[kk];
          fDist=(pow((fTmpXCoord2-inp_surf->ppfAtomCoord[iIndex1][0]), 2)+
                 pow((fTmpYCoord2-inp_surf->ppfAtomCoord[iIndex1][1]), 2)+
                 pow((fTmpZCoord2-inp_surf->ppfAtomCoord[iIndex1][2]), 2));
                 
          fDistTest=(inp_surf->pdAtomsRadii[iIndex1]*inp_surf->pdAtomsRadii[iIndex1]);
          if(fDist<fDistTest)
          {
            iBreakFlag=1;
            break;
          }
        } // kk, LL
        
        if(iBreakFlag==1)
        {
          iBreakFlag=0;
          continue;
        }
        else
        {
          iBreakFlag=1;
          break;
        }
      } // jj, L, 4
            
      if(iBreakFlag==1)
      {
        iBreakFlag=0;
      }
      else
      {
        inp_surf->piAtomType[ii]=4;
      }
    }
  } // ii, I, 5 

  return iNumOfStartingAtoms;
}

void Clean5(struct inp_surf *inp_surf)
{
  /*
   * This discard the atoms of area zero that are not needed for
   * the computation of the surface
   * 
   * original name: CLEAN5
   */
  
  int     ii, jj, kk, mm, oo, ww;
  int     iIndex1=0, iIndex2=0, iIndex3=0, iIndex4=0, iIndex5=0;
  int     iBreakFlag=0;
  
  float   fXCoord1=0.0, fYCoord1=0.0, fZCoord1=0.0, fRadius1=0.0;
  float   fXCoord2=0.0, fYCoord2=0.0, fZCoord2=0.0;
  float   fDist=0.0, fTest=0.0;
    
  for(ii=0; ii<inp_surf->iNumOfRealAtoms; ii++)
  {
    inp_surf->piUsefulAtoms[ii]=1;
  }
  
  for(ii=inp_surf->iNumOfRealAtoms; ii<inp_surf->iNumOfAtoms; ii++)
  {
    if(inp_surf->piAtomType[ii]!=4)
      inp_surf->piUsefulAtoms[ii]=1;
      
    else
      inp_surf->piUsefulAtoms[ii]=0;
  }
    
  for(ii=0; ii<inp_surf->iNumOfAtoms; ii++)
  {
    if(inp_surf->piAtomType[ii]>=4)
    {
      fXCoord1=inp_surf->ppfAtomCoord[ii][0];
      fYCoord1=inp_surf->ppfAtomCoord[ii][1];
      fZCoord1=inp_surf->ppfAtomCoord[ii][2];
      fRadius1=inp_surf->pdAtomsRadii[ii];
      
      // Find atoms that overlap atom ii-th
      iIndex1=-1;
      for(jj=0; jj<inp_surf->iNumOfAtoms; jj++)
      {
        if(inp_surf->piAtomType[jj]>=3)
        {
          fDist=(((fXCoord1-inp_surf->ppfAtomCoord[jj][0])*(fXCoord1-inp_surf->ppfAtomCoord[jj][0]))+
                 ((fYCoord1-inp_surf->ppfAtomCoord[jj][1])*(fYCoord1-inp_surf->ppfAtomCoord[jj][1]))+
                 ((fZCoord1-inp_surf->ppfAtomCoord[jj][2])*(fZCoord1-inp_surf->ppfAtomCoord[jj][2])));
                 
          fTest=((fRadius1+inp_surf->pdAtomsRadii[jj])*(fRadius1+inp_surf->pdAtomsRadii[jj]));
          if(fDist<fTest)
          {
            if(ii!=jj)
            {
              iIndex1++;
              inp_surf->piZeroSurfC5[iIndex1]=jj;
            }
          }
        }
      }
      
      // Find atoms that are needed to discard the triangle kk-th
      iIndex2=0;
      for(kk=0; kk<inp_surf->iNumOfTesserae; kk++)
      {
        iBreakFlag=0;
        fXCoord2=inp_surf->pfXCoord1[kk]*fRadius1+fXCoord1;
        fYCoord2=inp_surf->pfYCoord1[kk]*fRadius1+fYCoord1;        
        fZCoord2=inp_surf->pfZCoord1[kk]*fRadius1+fZCoord1;
        
        // Among the useful atoms
        iIndex3=iIndex2;
        for(mm=iIndex3; mm<iIndex1+1; mm++)
        {
          iIndex2=mm;
          iIndex4=inp_surf->piZeroSurfC5[mm];
          if(inp_surf->piUsefulAtoms[iIndex4]==1)
          {
            fDist=(((fXCoord2-inp_surf->ppfAtomCoord[iIndex4][0])*(fXCoord2-inp_surf->ppfAtomCoord[iIndex4][0]))+
                   ((fYCoord2-inp_surf->ppfAtomCoord[iIndex4][1])*(fYCoord2-inp_surf->ppfAtomCoord[iIndex4][1]))+
                   ((fZCoord2-inp_surf->ppfAtomCoord[iIndex4][2])*(fZCoord2-inp_surf->ppfAtomCoord[iIndex4][2])));
                   
            fTest=(inp_surf->pdAtomsRadii[iIndex4]*inp_surf->pdAtomsRadii[iIndex4]);
            
            if(fDist<fTest)
            {
              iBreakFlag=1;
              break;
            }
          }
        }
        
        if(iBreakFlag==1)
        {
          iBreakFlag=0;
          continue;
        }
        
        for(oo=0; oo<iIndex3+1; oo++)
        {
          iIndex2=oo;
          iIndex4=inp_surf->piZeroSurfC5[oo];
          if(inp_surf->piUsefulAtoms[iIndex4]==1)
          {
            fDist=(((fXCoord2-inp_surf->ppfAtomCoord[iIndex4][0])*(fXCoord2-inp_surf->ppfAtomCoord[iIndex4][0]))+
                   ((fYCoord2-inp_surf->ppfAtomCoord[iIndex4][1])*(fYCoord2-inp_surf->ppfAtomCoord[iIndex4][1]))+
                   ((fZCoord2-inp_surf->ppfAtomCoord[iIndex4][2])*(fZCoord2-inp_surf->ppfAtomCoord[iIndex4][2])));
                   
            fTest=(inp_surf->pdAtomsRadii[iIndex4]*inp_surf->pdAtomsRadii[iIndex4]);
            
            if(fDist<fTest)
            {
              iBreakFlag=1;
              break;
            }
          }
        }
        
        if(iBreakFlag==1)
        {
          iBreakFlag=0;
          continue;
        }
        
        // Among the not useful
        for(ww=0; ww<iIndex1+1; ww++)
        {
          iIndex2=ww;
          iIndex4=inp_surf->piZeroSurfC5[ww];
          if(inp_surf->piUsefulAtoms[iIndex4]==0)
          {
            fDist=(((fXCoord2-inp_surf->ppfAtomCoord[iIndex4][0])*(fXCoord2-inp_surf->ppfAtomCoord[iIndex4][0]))+
                   ((fYCoord2-inp_surf->ppfAtomCoord[iIndex4][1])*(fYCoord2-inp_surf->ppfAtomCoord[iIndex4][1]))+
                   ((fZCoord2-inp_surf->ppfAtomCoord[iIndex4][2])*(fZCoord2-inp_surf->ppfAtomCoord[iIndex4][2])));
                   
            fTest=(inp_surf->pdAtomsRadii[iIndex4]*inp_surf->pdAtomsRadii[iIndex4]);
            
            if(fDist<fTest)
            {
              inp_surf->piUsefulAtoms[iIndex4]=1;
              break;
            }
          }
        }
      }
    }
  }
  
  // Discard atoms totally inside of others
  iIndex4=0;
  for(ii=inp_surf->iNumOfRealAtoms; ii<inp_surf->iNumOfAtoms; ii++)
  {    
    iIndex5=ii-iIndex4;
    inp_surf->ppfAtomCoord[iIndex5][0]=inp_surf->ppfAtomCoord[ii][0];
    inp_surf->ppfAtomCoord[iIndex5][1]=inp_surf->ppfAtomCoord[ii][1];
    inp_surf->ppfAtomCoord[iIndex5][2]=inp_surf->ppfAtomCoord[ii][2];
    inp_surf->pdAtomsRadii[iIndex5]=inp_surf->pdAtomsRadii[ii];
    inp_surf->piAtomType[iIndex5]=inp_surf->piAtomType[ii];
    
    if(inp_surf->piUsefulAtoms[ii]==0) iIndex4++;
  }
  
  inp_surf->iNumOfAtoms=inp_surf->iNumOfAtoms-iIndex4;  
  
}

void GenAtoms(struct inp_surf *inp_surf)
{
  /*
   * This creates the new atoms
   * 
   * original name: CREA
   */
  
  int     ii, jj, kk;
  int     iTmpLargestAtom=0, iTmpSmallestAtom=0;
  int     iIndex1=0;
  int     iFirstAtom=0, iLastAtom=0;
  
  float   fPairDist1=0.0, fPairDist2=0.0, fPairDist3=0.0;
  float   fDistTest=0.0;
  float   fTmpOFAC=0.0, fRMid=0.0, fFC=0.0;
  float   fTmpLargestRad=0.0, fTmpSmallestRad=0.0;
  float   fTmpLargestRadSolv2=0.0, fTmpSmallestRadSolv2=0.0;
  float   fRgn=0.0, fRin=0.0, fRenD2A=0.0, fRenD2C=0.0;
  float   fXCoord1=0.0, fYCoord1=0.0, fZCoord1=0.0;
  float   fRadius1=0.0, fRadius2=0.0;
  
  fTmpOFAC=1.0-inp_surf->fOFAC*2.0;
  fRMid=((inp_surf->fRMIN+inp_surf->dSolventRadius)*(inp_surf->fRMIN+inp_surf->dSolventRadius));
  
  iFirstAtom=1;
  iLastAtom=inp_surf->iNumOfAtoms;
  
  while(1)
  {
    // Loop to select the first atom of the pair
    for(ii=iFirstAtom; ii<iLastAtom; ii++)
    {
      if(inp_surf->piAtomType[ii]<=4) goto L602;
      
      for(jj=0; jj<inp_surf->iNumOfAtoms; jj++)
      {
        inp_surf->pfDistVect[jj]=sqrtf(((inp_surf->ppfAtomCoord[ii][0]-inp_surf->ppfAtomCoord[jj][0])*
                                        (inp_surf->ppfAtomCoord[ii][0]-inp_surf->ppfAtomCoord[jj][0]))+
                                       ((inp_surf->ppfAtomCoord[ii][1]-inp_surf->ppfAtomCoord[jj][1])*
                                        (inp_surf->ppfAtomCoord[ii][1]-inp_surf->ppfAtomCoord[jj][1]))+
                                       ((inp_surf->ppfAtomCoord[ii][2]-inp_surf->ppfAtomCoord[jj][2])*
                                        (inp_surf->ppfAtomCoord[ii][2]-inp_surf->ppfAtomCoord[jj][2])));
      }
      
      // Loop to select the second atom of the pair
      for(jj=0; jj<ii; jj++)
      {
        if(inp_surf->piAtomType[jj]<=4)
          goto L603;
          
        fPairDist1=inp_surf->pfDistVect[jj];
        
        // If the solvent can pass through the pair this pair is discarded
        fDistTest=(inp_surf->dSolventRadius+inp_surf->dSolventRadius+
                   inp_surf->pdAtomsRadii[ii]+inp_surf->pdAtomsRadii[jj]);
                  
        if(fPairDist1>=fDistTest)
          goto L603;
        
        // Determine which is the largest and smallest atom
        if(inp_surf->pdAtomsRadii[ii]>inp_surf->pdAtomsRadii[jj])
        {
          fTmpLargestRad=inp_surf->pdAtomsRadii[ii];
          fTmpSmallestRad=inp_surf->pdAtomsRadii[jj];
          iTmpLargestAtom=ii;
          iTmpSmallestAtom=jj;
        }
        else
        {
          fTmpLargestRad=inp_surf->pdAtomsRadii[jj];
          fTmpSmallestRad=inp_surf->pdAtomsRadii[ii];
          iTmpLargestAtom=jj;
          iTmpSmallestAtom=ii;
        }
        
        fTmpLargestRadSolv2=((fTmpLargestRad+inp_surf->dSolventRadius)*
                             (fTmpLargestRad+inp_surf->dSolventRadius));
        
        // Determine whether the atoms are overlapped
        if(fPairDist1<=(fTmpSmallestRad+fTmpLargestRad))
        {
          // Test of overlapping
          fDistTest=fTmpLargestRad+fTmpSmallestRad*fTmpOFAC;
          if(fPairDist1<fDistTest)
            goto L603;
          
          // Test of the small atom
          fTmpSmallestRadSolv2=((fTmpSmallestRad+inp_surf->dSolventRadius)*
                                (fTmpSmallestRad+inp_surf->dSolventRadius));
                                
          fRgn=(fPairDist1-fTmpSmallestRad+fTmpLargestRad)*0.5;
          fRenD2A=fTmpLargestRadSolv2+fRgn*(fRgn-(fTmpLargestRadSolv2+(fPairDist1*fPairDist1)-fTmpSmallestRadSolv2)/fPairDist1);          
          if(fRenD2A<=fRMid)
            goto L603;
          
          // Atoms A
          fFC=((fPairDist1-fTmpSmallestRad+fTmpLargestRad)/(fPairDist1+fTmpSmallestRad-fTmpLargestRad));
          fXCoord1=((inp_surf->ppfAtomCoord[iTmpLargestAtom][0]+fFC*inp_surf->ppfAtomCoord[iTmpSmallestAtom][0])/(fFC+1.0));
          fYCoord1=((inp_surf->ppfAtomCoord[iTmpLargestAtom][1]+fFC*inp_surf->ppfAtomCoord[iTmpSmallestAtom][1])/(fFC+1.0));
          fZCoord1=((inp_surf->ppfAtomCoord[iTmpLargestAtom][2]+fFC*inp_surf->ppfAtomCoord[iTmpSmallestAtom][2])/(fFC+1.0));
          fRadius1=sqrtf(fRenD2A)-inp_surf->dSolventRadius;
          fRin=(fPairDist1-inp_surf->pdAtomsRadii[jj]+inp_surf->pdAtomsRadii[ii])*0.5;          
        }
        else
        {
          // Test of the small atom
          fTmpSmallestRadSolv2=((fTmpSmallestRad+inp_surf->dSolventRadius)*(fTmpSmallestRad+inp_surf->dSolventRadius));
          fRenD2C=(fTmpLargestRadSolv2+(fTmpLargestRad*fTmpLargestRad)-(fTmpLargestRad/fPairDist1)*
                  (fTmpLargestRadSolv2+(fPairDist1*fPairDist1)-fTmpSmallestRadSolv2));
          
          if(fRenD2C<=fRMid)
            goto L603;
            
          // Calculate radius for atom of kind B and separate B and C
          fRgn=(fPairDist1-fTmpSmallestRad+fTmpLargestRad)*0.5;
          fRenD2A=(fTmpLargestRadSolv2+fRgn*(fRgn-(fTmpLargestRadSolv2+(fPairDist1*fPairDist1)-fTmpSmallestRadSolv2)/fPairDist1));
          
          if(fRenD2A>fRMid)
          {
            // Atoms B
            fFC=((fPairDist1-fTmpSmallestRad+fTmpLargestRad)/
                 (fPairDist1+fTmpSmallestRad-fTmpLargestRad));
            fXCoord1=((inp_surf->ppfAtomCoord[iTmpLargestAtom][0]+fFC*inp_surf->ppfAtomCoord[iTmpSmallestAtom][0])/(fFC+1.0));
            fYCoord1=((inp_surf->ppfAtomCoord[iTmpLargestAtom][1]+fFC*inp_surf->ppfAtomCoord[iTmpSmallestAtom][1])/(fFC+1.0));
            fZCoord1=((inp_surf->ppfAtomCoord[iTmpLargestAtom][2]+fFC*inp_surf->ppfAtomCoord[iTmpSmallestAtom][2])/(fFC+1.0));
            fRadius1=sqrtf(fRenD2A)-inp_surf->dSolventRadius;
            fRin=(fPairDist1-inp_surf->pdAtomsRadii[jj]+inp_surf->pdAtomsRadii[ii])*0.5;            
          }
          else
          {
            // Atoms C
            fFC=fTmpLargestRad/(fPairDist1-fTmpLargestRad);
            fXCoord1=((inp_surf->ppfAtomCoord[iTmpLargestAtom][0]+fFC*inp_surf->ppfAtomCoord[iTmpSmallestAtom][0])/(fFC+1.0));
            fYCoord1=((inp_surf->ppfAtomCoord[iTmpLargestAtom][1]+fFC*inp_surf->ppfAtomCoord[iTmpSmallestAtom][1])/(fFC+1.0));
            fZCoord1=((inp_surf->ppfAtomCoord[iTmpLargestAtom][2]+fFC*inp_surf->ppfAtomCoord[iTmpSmallestAtom][2])/(fFC+1.0));
            fRadius1=sqrtf(fRenD2C)-inp_surf->dSolventRadius;
            
            if(iTmpLargestAtom==ii)
              fRin=fTmpLargestRad;
            
            else
              fRin=(fPairDist1-fTmpLargestRad);
          }
        }
        
        // Test of overlapping for the new atom
        iIndex1=0;
        for(kk=0; kk<inp_surf->iNumOfAtoms; kk++)
        {
          fPairDist2=inp_surf->pfDistVect[kk];
          
          if(fPairDist2>=(fRin+fRadius1+inp_surf->pdAtomsRadii[kk]))
            goto L604;
            
          if(inp_surf->piAtomType[kk]<=3)
            goto L604;
            
          fPairDist3=(((fXCoord1-inp_surf->ppfAtomCoord[kk][0])*(fXCoord1-inp_surf->ppfAtomCoord[kk][0]))+
                      ((fYCoord1-inp_surf->ppfAtomCoord[kk][1])*(fYCoord1-inp_surf->ppfAtomCoord[kk][1]))+
                      ((fZCoord1-inp_surf->ppfAtomCoord[kk][2])*(fZCoord1-inp_surf->ppfAtomCoord[kk][2])));
          
          fRadius2=inp_surf->pdAtomsRadii[kk];
          
          if(fPairDist3>=((fRadius2+fRadius1)*(fRadius2+fRadius1)))
            goto L604;
          
          fDistTest=((fRadius2-fRadius1)*(fRadius2-fRadius1));
          
          if(fPairDist3<=fDistTest)
          {
            if(fRadius1>fRadius2)
            {
              if(fDistTest<0.0004)
                goto L603;

                inp_surf->piEngulfedAtoms[iIndex1]=kk;
                iIndex1++;
                goto L604;
            }
            else
              goto L603;
          }
          
          fDistTest=fRadius2+fRadius1*fTmpOFAC;
          if(fDistTest<0.0)
            goto L604;
          
          fDistTest=fDistTest*fDistTest;
          if(fPairDist3<=fDistTest)
            goto L603;

          L604:;
        } // kk, K, 604

        // Mark atoms engulfed by the new atom
        for(kk=0; kk<iIndex1; kk++)
          inp_surf->piAtomType[inp_surf->piEngulfedAtoms[kk]]=2;
        
        // Creates the new atoms
        if((inp_surf->iNumOfAtoms+1)>=inp_surf->iMaxSize)
        {
          printf("SURF Memory Error\n");
          /*printf("Too many atoms created -> %d\n", inp_surf->iNumOfAtoms);*/
          printf("Increment --MEMSIZE (now %d)\n",  inp_surf->memsize);
          printf("Exiting now\n");
          exit(1);
        }
        
        inp_surf->ppfAtomCoord[inp_surf->iNumOfAtoms][0]=fXCoord1;
        inp_surf->ppfAtomCoord[inp_surf->iNumOfAtoms][1]=fYCoord1;
        inp_surf->ppfAtomCoord[inp_surf->iNumOfAtoms][2]=fZCoord1;
        
        inp_surf->pdAtomsRadii[inp_surf->iNumOfAtoms]=fRadius1;
        inp_surf->pfDistVect[inp_surf->iNumOfAtoms]=fRin;
        inp_surf->piAtomType[inp_surf->iNumOfAtoms]=6;
        
        inp_surf->iNumOfAtoms++;

        L603:;
      } // jj, J, 603

      L602:;
    } // ii, I, 602
    
    if(inp_surf->iNumOfAtoms!=iLastAtom)
    {
      iLastAtom=MZero5(inp_surf, iLastAtom);
      if(inp_surf->iNumOfAtoms!=iLastAtom)
      {
        iFirstAtom=iLastAtom;
        iLastAtom=inp_surf->iNumOfAtoms;
      }
      else
        break;
    }
    else
      break;
    
  } // 600
}

void AddSolvRad(struct inp_surf *inp_surf)
{
  /*
   * Add the solvent radius to atom radii
   * 
   * original name: SUM
   */
  
  int     ii;
  
  for(ii=0; ii<inp_surf->iNumOfAtoms; ii++)
  {
    if(inp_surf->piAtomType[ii]!=1)
    {
      inp_surf->pdAtomsRadii[ii]=inp_surf->pdAtomsRadii[ii]+inp_surf->dSolventRadius;
    }
  }
}

void GeoCav(struct inp_surf *inp_surf)
{
  /*
   * This computes the surface
   * 
   * original name: GEOCAV
   */

  int   ii, jj, kk, ww, zz;
  int   iBreakFlag=0;
  int   iNTRIAN=0;
  int   iIndex1=0, iIndex2=0, iIndex3=0; 
  int   iIndex4=0, iIndex5=0, iIndex6=0;
  int   iIndex7=0, iIndex8=0;
  
  float fTmpXCoord1=0.0, fTmpYCoord1=0.0, fTmpZCoord1=0.0, fTmpRadius1=0.0;
  float fTmpXCoord2=0.0, fTmpYCoord2=0.0, fTmpZCoord2=0.0;
  float fTmpXCoord3=0.0, fTmpYCoord3=0.0, fTmpZCoord3=0.0;
  float fTmpXCoord4=0.0, fTmpYCoord4=0.0, fTmpZCoord4=0.0;
  
  float fFC=0.0, fATS=0.0, fATP=0.0;
  float fTmpDist=0.0, fTmpRadSum1=0.0, fTmpRadSum2=0.0, fTmpRadSum3=0.0;
  float PI=3.14159265358979323846264;
  
  inp_surf->iNumOfPoints=0;
  iNTRIAN=pow(4, (inp_surf->iNDIV-1));

  for(ii=0; ii<inp_surf->iNumOfAtoms; ii++)
  {
    // Select one atom
    if(inp_surf->piAtomType[ii]==6)
    {
      fTmpRadius1=inp_surf->pdAtomsRadii[ii];
      fTmpXCoord1=inp_surf->ppfAtomCoord[ii][0];
      fTmpYCoord1=inp_surf->ppfAtomCoord[ii][1];
      fTmpZCoord1=inp_surf->ppfAtomCoord[ii][2];
      
      fATS=4.0*PI*fTmpRadius1*fTmpRadius1/(float)inp_surf->iNumOfTesserae;
      
      iIndex1=-1;
      // Determines which atoms are linked to atom-iith
      for(jj=0; jj<inp_surf->iNumOfAtoms; jj++)
      {
        if(inp_surf->piAtomType[jj]>=3)
        {
          if(ii!=jj)
          {
            fTmpDist=(((fTmpXCoord1-inp_surf->ppfAtomCoord[jj][0])*(fTmpXCoord1-inp_surf->ppfAtomCoord[jj][0]))+
                      ((fTmpYCoord1-inp_surf->ppfAtomCoord[jj][1])*(fTmpYCoord1-inp_surf->ppfAtomCoord[jj][1]))+
                      ((fTmpZCoord1-inp_surf->ppfAtomCoord[jj][2])*(fTmpZCoord1-inp_surf->ppfAtomCoord[jj][2])));
            
            fTmpRadSum1=((fTmpRadius1+inp_surf->pdAtomsRadii[jj])*(fTmpRadius1+inp_surf->pdAtomsRadii[jj]));
            
            if(fTmpDist<fTmpRadSum1)
            {
              iIndex1++;
              inp_surf->piIndexVect[iIndex1]=jj;
              iIndex2=iIndex1;
            }
          }
        }
      }
      
      // Selects one main triangle
      iIndex3=-1;
      iIndex4=0;
      for(jj=0; jj<60; jj++)
      {
        fTmpXCoord2=0.0;
        fTmpYCoord2=0.0;
        fTmpZCoord2=0.0;
        iIndex5=0;
        iIndex6=iIndex3+1;
        iIndex3=iIndex6+iNTRIAN-1;
        
        // Selects one secondary triangle
        for(kk=iIndex6; kk<iIndex3+1; kk++)
        {
          fTmpXCoord3=inp_surf->pfXCoord1[kk]*fTmpRadius1;
          fTmpYCoord3=inp_surf->pfYCoord1[kk]*fTmpRadius1;
          fTmpZCoord3=inp_surf->pfZCoord1[kk]*fTmpRadius1;

          fTmpXCoord4=fTmpXCoord3+fTmpXCoord1;
          fTmpYCoord4=fTmpYCoord3+fTmpYCoord1;
          fTmpZCoord4=fTmpZCoord3+fTmpZCoord1;
                    
          // Fixes if the secondary triangle is inside or outside
          iIndex7=iIndex4;
          iBreakFlag=0;

          for(ww=iIndex7; ww<iIndex2+1; ww++)
          {
            iIndex4=ww;
            iIndex8=inp_surf->piIndexVect[ww];
            fTmpDist=(((fTmpXCoord4-inp_surf->ppfAtomCoord[iIndex8][0])*(fTmpXCoord4-inp_surf->ppfAtomCoord[iIndex8][0]))+
                      ((fTmpYCoord4-inp_surf->ppfAtomCoord[iIndex8][1])*(fTmpYCoord4-inp_surf->ppfAtomCoord[iIndex8][1]))+
                      ((fTmpZCoord4-inp_surf->ppfAtomCoord[iIndex8][2])*(fTmpZCoord4-inp_surf->ppfAtomCoord[iIndex8][2])));
            
            fTmpRadSum2=(inp_surf->pdAtomsRadii[iIndex8]*inp_surf->pdAtomsRadii[iIndex8]);
            if(fTmpDist<fTmpRadSum2)
            {              
              iBreakFlag=1;
              break;
            }
          }
          
          if(iBreakFlag==1)
          {
            iBreakFlag=0;
            continue;
          }

          iBreakFlag=0;
          for(zz=0; zz<(iIndex7); zz++)
          {
            iIndex4=zz;
            iIndex8=inp_surf->piIndexVect[zz];
            fTmpDist=(((fTmpXCoord4-inp_surf->ppfAtomCoord[iIndex8][0])*(fTmpXCoord4-inp_surf->ppfAtomCoord[iIndex8][0]))+
                      ((fTmpYCoord4-inp_surf->ppfAtomCoord[iIndex8][1])*(fTmpYCoord4-inp_surf->ppfAtomCoord[iIndex8][1]))+
                      ((fTmpZCoord4-inp_surf->ppfAtomCoord[iIndex8][2])*(fTmpZCoord4-inp_surf->ppfAtomCoord[iIndex8][2])));
            
            fTmpRadSum2=(inp_surf->pdAtomsRadii[iIndex8]*inp_surf->pdAtomsRadii[iIndex8]);
            if(fTmpDist<fTmpRadSum2)
            {
              iBreakFlag=1;
              break;
            }
          }
          
          if(iBreakFlag==1)
          {
            iBreakFlag=0;
            continue;
          }         
          
          // Prepares the coordinates for the main triangle
          fTmpXCoord2=fTmpXCoord2+fTmpXCoord3;
          fTmpYCoord2=fTmpYCoord2+fTmpYCoord3;
          fTmpZCoord2=fTmpZCoord2+fTmpZCoord3;
          iIndex5++;
        }
        
        // Reduces the secondary triangles to the main triangle
        if(iIndex5>0)
        {
          fATP=fATS*iIndex5;
          fTmpXCoord2=fTmpXCoord2/(float)iIndex5;
          fTmpYCoord2=fTmpYCoord2/(float)iIndex5;
          fTmpZCoord2=fTmpZCoord2/(float)iIndex5;
          fTmpRadSum3=sqrtf(fTmpXCoord2*fTmpXCoord2+fTmpYCoord2*fTmpYCoord2+fTmpZCoord2*fTmpZCoord2);
          fFC=fTmpRadius1/fTmpRadSum3;
          
          inp_surf->ppfPCoord[inp_surf->iNumOfPoints][0]=fTmpXCoord2*fFC+fTmpXCoord1;
          inp_surf->ppfPCoord[inp_surf->iNumOfPoints][1]=fTmpYCoord2*fFC+fTmpYCoord1;
          inp_surf->ppfPCoord[inp_surf->iNumOfPoints][2]=fTmpZCoord2*fFC+fTmpZCoord1;
          inp_surf->pfAPVect[inp_surf->iNumOfPoints]=fATP;          
          inp_surf->iNumOfPoints++;
          if(inp_surf->iNumOfPoints>=inp_surf->iMaxSize)
          {
          printf("SURF Memory Error\n");
          /*printf("Too many atoms created -> %d\n", inp_surf->iNumOfAtoms);*/
          printf("Increment --MEMSIZE (now %d)\n",  inp_surf->memsize);
          printf("Exiting now\n");
            exit(1);
          }
        }
      }
    }
  }
}

float AreaSum(struct inp_surf *inp_surf)
{
  /* 
   * This calculates total area
   *
   * original name: VOLARE
   */
  
  int     ii;
  float   fArea=0.0;
  
  for(ii=0; ii<inp_surf->iNumOfPoints; ii++)
  {
    fArea=fArea+inp_surf->pfAPVect[ii];
  }

  return fArea;
}

// ============================================================================================================

int Read_SURF(char **input, int inp_index, struct inp_surf *inp_surf, Molecule *molecule, char *outstring) //FILE *OutputFile)
{ 
  int     ii, jj, kk, mm;
  int     iGepolAlgoOptionsCheck=0;
  int     iSizeMul=20;
  int     iWinpOptFlag=0, iTotAtoms=0;

  char    cWordomInpBuffer[1024];
  char    cTmpString[6];
  int     iRadFileFlag=0;
  char    radfilename[128];
  RadList *radlist;
  
  int     unitsize;
  char    membuffer[128];
  float   memvalue;
  
  extern short int     frame_par;
  frame_par = 1;
  
  srand(time(NULL));
  
  // === Default Values ===========
//sprintf(inp_surf->cTitle, "SAS");
  inp_surf->dSolventRadius=1.4;
  inp_surf->iWriteSurfFlag=1;
  inp_surf->iAlgoFlag=0;
  inp_surf->iSurfType=1;
  inp_surf->iNDIV=3;
  inp_surf->fOFAC=-1;
  inp_surf->fRMIN=-1;
  inp_surf->iGhostFlag=0;
  // ==============================
  
  unitsize = molecule->nato * ( sizeof(double) + 2*sizeof(float **) + 8*sizeof(int) + 8*sizeof(float) );
    
  // === example winp file ===
  //BEGIN
  //--TITLE SURFTEST
  //--ALGO GEPOL
  //--CALC ESURF
  //--SOLVRAD 1.4
  //--SELE /*/*/*
  //--MEMSIZE 10
  //END
  // =========================
  
  memset ( cWordomInpBuffer, '\0', sizeof(cWordomInpBuffer));
  while( strncmp (cWordomInpBuffer, "END", 3))
  {
    iWinpOptFlag = 0;
    sprintf( cWordomInpBuffer, "%s", input[inp_index]);
    if( !strncmp(cWordomInpBuffer, "BEGIN", 5) || !strncmp(cWordomInpBuffer, "END", 3) || cWordomInpBuffer[0] == '#')
      iWinpOptFlag = 1;
    else if ( !strncmp(cWordomInpBuffer, "--TITLE", 7))
    {
      sscanf( cWordomInpBuffer, "--TITLE %s", inp_surf->cTitle);
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--SOLVRAD", 9))
    {
      sscanf( cWordomInpBuffer, "--SOLVRAD %f", &inp_surf->dSolventRadius);
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--ALGO", 6))
    {
      sscanf( cWordomInpBuffer, "--ALGO %s", cTmpString);
      iWinpOptFlag = 1;
      if(!strncmp(cTmpString, "ARVO", 4) || !strncmp(cTmpString, "arvo", 4))
      {
        inp_surf->iAlgoFlag=0;
      }
      else if (!strncmp(cTmpString, "GEPOL", 5) || !strncmp(cTmpString, "gepol", 5))
      {
        inp_surf->iAlgoFlag=1;
      }
      else
      {
        fprintf( stderr, "SURF module: Invalid --ALGO option: %s\n", cTmpString);
        exit(5);
      }
    }
    else if ( !strncmp(cWordomInpBuffer, "--CALC", 6))
    {
      sscanf( cWordomInpBuffer, "--CALC %s", cTmpString);
      iWinpOptFlag = 1;
      if(!strncmp(cTmpString, "WSURF", 5) || !strncmp(cTmpString, "wsurf", 5))
      {
        iGepolAlgoOptionsCheck=1;
        inp_surf->iSurfType=0;
      }
      else if (!strncmp(cTmpString, "ASURF", 5) || !strncmp(cTmpString, "asurf", 5))
      {
        iGepolAlgoOptionsCheck=1;
        inp_surf->iSurfType=1;
      }
      else if (!strncmp(cTmpString, "ESURF", 5) || !strncmp(cTmpString, "esurf", 5))
      {
        iGepolAlgoOptionsCheck=1;
        inp_surf->iSurfType=2;
      }
      else
      {
        fprintf( stderr, "SURF module: Invalid --CALC option: %s\n", cTmpString);
        exit(5);
      }
    }
    else if ( !strncmp(cWordomInpBuffer, "--OFAC", 6))
    {
      sscanf( cWordomInpBuffer, "--OFAC %f", &inp_surf->fOFAC);
      if(inp_surf->fOFAC<0)
      {
        fprintf( stderr, "SURF module: --OFAC must be a real number in the range 0.0..1.0\n");
        exit(5);
      }
      iGepolAlgoOptionsCheck=1;
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--RMIN", 6))
    {
      sscanf( cWordomInpBuffer, "--RMIN %f", &inp_surf->fRMIN);
      if(inp_surf->fRMIN<0)
      {
        fprintf( stderr, "SURF module: --RMIN must be a real number larger than 0.0\n");
        exit(5);
      }
      iGepolAlgoOptionsCheck=1;
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--NDIV", 6))
    {
      sscanf( cWordomInpBuffer, "--NDIV %d", &inp_surf->iNDIV);
      if(inp_surf->iNDIV<0 || inp_surf->iNDIV>5)
      {
        fprintf( stderr, "SURF module: --NDIV must be an integer number in the range 0..5\n");
        exit(5);
      }
      iGepolAlgoOptionsCheck=1;
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--MEMSIZE", 9))
    {
      //sscanf( cWordomInpBuffer, "--MEMSIZE %d", &iSizeMul);s
      sscanf( cWordomInpBuffer, "--MEMSIZE %[^\n]", membuffer);
      memvalue = atof(membuffer);
      if( strchr(membuffer, 'k')!=NULL || strchr(membuffer, 'K')!=NULL )
        memvalue *= 1024;
      else if( strchr(membuffer, 'm')!=NULL || strchr(membuffer, 'M')!=NULL )
        memvalue *= 1024*1024;
      else if( strchr(membuffer, 'g')!=NULL || strchr(membuffer, 'G')!=NULL )
        memvalue *= 1024*1024*1024;
      else
        fprintf( stderr, "You did not specify the memory size unit; taken as bytes\n");
      
      iSizeMul = (int) floor(memvalue/unitsize);
      //fprintf( stderr, "DEBUG: memvalue: %f; unitsize: %d; iSizeMul: %d; nato: %d\n", memvalue, unitsize, iSizeMul, molecule->nato);
      if(iSizeMul<1)
      {
        fprintf( stderr, "SURF module: --MEMSIZE must be at _least_ >= %d\n", unitsize);
        exit(5);
      }
      iGepolAlgoOptionsCheck=1;
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--SELE", 6))
    {
      sscanf(cWordomInpBuffer, "--SELE %[^\n]%*c ", inp_surf->sele1.selestring);     
      GetSele ( inp_surf->sele1.selestring, &inp_surf->sele1, molecule);
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--RADFILE", 9))
    {
      sscanf( cWordomInpBuffer, "--RADFILE %s", radfilename);
      iRadFileFlag = 1;
      radlist = ReadRadius(radfilename);
      iWinpOptFlag = 1;
    }
    if( iWinpOptFlag==0 )
    {
      fprintf( stderr, "SURF module: Could NOT understand option: %s\n", cWordomInpBuffer);
      exit(5);
    }
    inp_index++;
  }
  
  inp_surf->memsize = unitsize * iSizeMul;
  
  // check if Gepol Options has been set and Arvo algorithm has been chosen
  if(iGepolAlgoOptionsCheck==1 && inp_surf->iAlgoFlag==0)
  {
    fprintf( stderr, "SURF module:\n");
    fprintf( stderr, "ARVO Algorithm do not use --CALC, --NDIV, --OFAC, --RMIN and --MEMSIZE options\n");
    exit(5);
  }

  // check if OFAC and RMIN has been set in agreement with --CALC ESURF
  if(inp_surf->iSurfType!=2)
  {
    if(inp_surf->fRMIN!=-1 || inp_surf->fOFAC!=-1)
    {
      fprintf( stderr, "SURF module:\n");
      fprintf( stderr, "--OFAC and --RMIN options can be used only if --ALGO is GEPOL and --CALC is ESURF\n");
      exit(5);
    }
  }

  // Calculates the total number of atoms
  iTotAtoms=0;
  for(ii=0; ii<molecule->nSeg; ii++)
  {
    for(jj=0; jj<molecule->segment[ii].nRpS; jj++)
    {
      iTotAtoms=iTotAtoms+molecule->segment[ii].pRes[jj].nApR;
    }
  }
  
  inp_surf->iNumOfRealAtoms=iTotAtoms;                                  // Number of Atoms in PDB file
  inp_surf->iNumOfAtoms=iTotAtoms;
  
  if(inp_surf->iAlgoFlag==0)
  {
    /* ============ *\
    |* === ARVO === *|
    \* ============ */
    
    inp_surf->iNumOfAtoms=iTotAtoms;                                    // Number of Atoms within 10 A from selected atoms
    inp_surf->iCirclesPlaneForAtom=iTotAtoms*80;                        // Maximum number of local circles in the plane for one atom
    inp_surf->iArcsAnglesCircleIntersect=iTotAtoms*80;                  // Maximum number of arcs and angles which arise from the local circles intersections
    inp_surf->iNumOfNeighbors=iTotAtoms*80;                             // Maximum number of neighbors
    
    // This vector stores the radii of all atoms
    inp_surf->pdAtomsRadii=(double *) calloc(iTotAtoms, sizeof(double));

    if( iRadFileFlag)
    {
      RadAssign( radlist, molecule, inp_surf->pdAtomsRadii );
      for( ii=0; ii<molecule->nato; ii++)
        inp_surf->pdAtomsRadii[ii] += inp_surf->dSolventRadius;
    }
    else
    {
      mm=-1;
      for(ii=0; ii<molecule->nSeg; ii++)
      {
        for(jj=0; jj<molecule->segment[ii].nRpS; jj++)
        {
          for(kk=0; kk<molecule->segment[ii].pRes[jj].nApR; kk++)
          {
            mm++;
            inp_surf->pdAtomsRadii[mm]=(double)molecule->segment[ii].pRes[jj].pAto[kk].bFac+inp_surf->dSolventRadius;
          }
        }
      }
    }

    /*
      Allocates piNeighborsNumber
      This vector stores the number of neighbors of all atoms
    */
    inp_surf->piNeighborsNumber=(int*) calloc(iTotAtoms*40, sizeof(int));

    /*
      Allocates piNeighborsStartIndex
      This vector stores 
    */
    inp_surf->piNeighborsStartIndex=(int*) calloc(iTotAtoms*40, sizeof(int));

    /*
      Allocates piNeighborsStartIndex
      This vector stores 
    */
    inp_surf->piNeighborsIndices=(int*) calloc(inp_surf->iNumOfNeighbors, sizeof(int));

    /*
      Allocates piIndicesVector
      This vector stores 
    */
    inp_surf->piIndicesVector=(int*) calloc(iTotAtoms*40, sizeof(int));

    /*
      Allocates piIndicesVector2
      This vector stores 
    */
    inp_surf->piIndicesVector2=(int*) calloc(iTotAtoms*40, sizeof(int));


    /*
      Circles Allocation
    */
    inp_surf->ppdCircles=(double **) calloc(inp_surf->iCirclesPlaneForAtom, sizeof(double *));
    for(ii=0; ii<inp_surf->iCirclesPlaneForAtom; ii++)
    {
      inp_surf->ppdCircles[ii]=(double *) calloc(4, sizeof(double));
    }
    
    /*
      Arcs Allocation
    */
    inp_surf->ppdArcs=(double **) calloc(inp_surf->iArcsAnglesCircleIntersect, sizeof(double *));
    for(ii=0; ii<inp_surf->iArcsAnglesCircleIntersect; ii++)
    {
      inp_surf->ppdArcs[ii]=(double *) calloc(3, sizeof(double));
    }

    /*
      NewArcs Allocation
    */
    inp_surf->ppdNewArcs=(double **) calloc(inp_surf->iArcsAnglesCircleIntersect, sizeof(double *));
    for(ii=0; ii<inp_surf->iArcsAnglesCircleIntersect; ii++)
    {
      inp_surf->ppdNewArcs[ii]=(double *) calloc(3, sizeof(double));
    }
    
    /*
      Allocates pdAngles 
    */
    inp_surf->pdAngles=(double*) calloc(inp_surf->iArcsAnglesCircleIntersect, sizeof(double));

    /*
      Allocates pdNewAngles 
    */
    inp_surf->pdNewAngles=(double*) calloc(inp_surf->iArcsAnglesCircleIntersect, sizeof(double));

    /*
      Allocates AtomNeighbors
    */
    inp_surf->ppdAtomNeighbors=(double **) calloc(iTotAtoms, sizeof(double*));
    for(ii=0; ii<iTotAtoms; ii++)
    {
      // dAtomNeighbors dimesion: 0 -> xcoord, 1 -> ycoord, 2 -> zcoord, 3 -> radius
      inp_surf->ppdAtomNeighbors[ii]=(double *) calloc(4, sizeof(double));
    }

    /*
      Allocates piWithinAtoms
    */
    inp_surf->piWithinAtoms=(int *) calloc(inp_surf->iNumOfAtoms, sizeof(int ));

    inp_surf->ppfWithinAtomsCoord=(float **) calloc(inp_surf->iNumOfAtoms, sizeof(float *));
    for(ii=0; ii<inp_surf->iNumOfAtoms; ii++)
    {
      inp_surf->ppfWithinAtomsCoord[ii]=(float *) calloc(4, sizeof(float));
    }
  }
  
  else if(inp_surf->iAlgoFlag==1)
  {
    /* ============= *\
    |* === GEPOL === *|
    \* ============= */
    
    inp_surf->iSolvRadAddedFalg=0;
    inp_surf->iMaxSize=(inp_surf->iNumOfRealAtoms*iSizeMul);    
    inp_surf->pdAtomsRadii=(double *) calloc(inp_surf->iMaxSize, sizeof(double));
    
    // Gets Atom radii
    if( iRadFileFlag)
    {
      RadAssign( radlist, molecule, inp_surf->pdAtomsRadii );
    }
    else
    {
      mm=-1;
      for(ii=0; ii<molecule->nSeg; ii++)
      {
        for(jj=0; jj<molecule->segment[ii].nRpS; jj++)
        {
          for(kk=0; kk<molecule->segment[ii].pRes[jj].nApR; kk++)
          {
            mm++;
            inp_surf->pdAtomsRadii[mm]=(double)molecule->segment[ii].pRes[jj].pAto[kk].bFac;
          }
        }
      }
    }

    if(inp_surf->iAlgoFlag==1 && inp_surf->iSurfType==2)
    {
      if(inp_surf->fRMIN==-1)
        inp_surf->fRMIN=0.5;
        
      if(inp_surf->fOFAC==-1)
        inp_surf->fOFAC=0.8;
    }
    
    // Allocates ppfAtomCoord
    inp_surf->ppfAtomCoord=(float **) calloc(inp_surf->iMaxSize, sizeof(float *));
    for(ii=0; ii<inp_surf->iMaxSize; ii++)
    {
      inp_surf->ppfAtomCoord[ii]=(float *) calloc(3, sizeof(float));
    }

    // Allocates piAtomType
    inp_surf->piAtomType=(int *) calloc(inp_surf->iMaxSize, sizeof(int));
    
    // Fills JVT1 & JVT2 arrays
    FillJVTArrays(inp_surf);
     
    // Allocates pfXCoord1, pfYCoord1, pfZCoord1
    inp_surf->iNumOfTesserae=60*pow(4, (inp_surf->iNDIV-1));
    inp_surf->pfXCoord1=(float *) calloc(inp_surf->iNumOfTesserae, sizeof(float));
    inp_surf->pfYCoord1=(float *) calloc(inp_surf->iNumOfTesserae, sizeof(float));
    inp_surf->pfZCoord1=(float *) calloc(inp_surf->iNumOfTesserae, sizeof(float));
    
    CompTriaVertexCoord(inp_surf);
    
    Divide(inp_surf);
    
    if(inp_surf->iSurfType==2)
    {
      // Allocates pfDistVect vector, used in Bulk
      inp_surf->pfDistVect=(float *) calloc(inp_surf->iMaxSize, sizeof(float));
      
      // Allocates piEngulfedAtoms
      inp_surf->piEngulfedAtoms=(int *) calloc(inp_surf->iMaxSize, sizeof(int));
      
      // Allocates piZeroSurf, piZeroSurfMZ5 and piZeroSurfC5
      inp_surf->piZeroSurf    =(int *) calloc(inp_surf->iMaxSize, sizeof(int));
      inp_surf->piZeroSurfMZ5 =(int *) calloc(inp_surf->iMaxSize, sizeof(int));
      inp_surf->piZeroSurfC5  =(int *) calloc(inp_surf->iMaxSize, sizeof(int));
      
      // Allocates piUsefulAtoms
      inp_surf->piUsefulAtoms = (int *) calloc(inp_surf->iMaxSize, sizeof(int));

      // Allocates piZeroCheck
      inp_surf->piZeroCheck=(int *) calloc(inp_surf->iMaxSize, sizeof(int));
    }
    
    // Allocates piIndexVect, used in GeoCav
    inp_surf->piIndexVect = (int *) calloc(inp_surf->iMaxSize, sizeof(int));
    
    // Allocated ppfPCoord
    inp_surf->ppfPCoord = (float **) calloc(inp_surf->iMaxSize, sizeof(float *));
    for(ii=0; ii<inp_surf->iMaxSize; ii++)
      inp_surf->ppfPCoord[ii] = (float *) calloc(3, sizeof(float));

    // Allocates pfAPVect
    inp_surf->pfAPVect = (float *) calloc(inp_surf->iMaxSize, sizeof(float));

  }
  if(inp_surf->iWriteSurfFlag==1)
  {
    sprintf( outstring, " %15s ", inp_surf->cTitle);
    return 17;
  }
  
  return 0;
}

// ========================================================================
int Compute_SURF(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring) //FILE *OutputFile)
{
  float  fSurf=0.0;
  
  if(inp_surf->iAlgoFlag==0)
    fSurf=ComputeARVOSurf(inp_surf, molecule, trj_crd, outstring);

  else if(inp_surf->iAlgoFlag==1)
    fSurf=ComputeGEPOLSurf(inp_surf, molecule, trj_crd, outstring);
  
  if(inp_surf->iWriteSurfFlag==1)
    sprintf( outstring, " %15.3f ", fSurf);

  return 17;
}
// ========================================================================
int Read_SURFCLUST(char **input, int inp_index, struct inp_surf *inp_surf, Molecule *molecule, char *outstring, int iNumOfFrames) //FILE *OutputFile, int iNumOfFrames)
{ 
  int     ii, jj, kk, mm;
  int     iGepolAlgoOptionsCheck=0;
  int     iSizeMul=20;
  int     iWinpOptFlag=0, iTotAtoms=0;

  char    cWordomInpBuffer[1024];
  char    cTmpString[6];
  int     iRadFileFlag=0;
  char    radfilename[128];
  RadList *radlist;

  int     unitsize;
  char    membuffer[128];
  float   memvalue;
  
  extern short int     frame_par;
  frame_par = 0;

  srand(time(NULL));
 
  // === Default Values ===========
  sprintf(inp_surf->cTitle, "SAS");
  inp_surf->dSolventRadius=1.4;
  inp_surf->iWriteSurfFlag=1;
  inp_surf->iAlgoFlag=0;
  inp_surf->iSurfType=1;
  inp_surf->iNDIV=3;
  inp_surf->fOFAC=-1;
  inp_surf->fRMIN=-1;
  inp_surf->iGhostFlag=0;
  inp_surf->fClustBin=10.0;
  // ==============================
    
  unitsize = molecule->nato * ( sizeof(double) + 2*sizeof(float **) + 8*sizeof(int) + 8*sizeof(float) );
    
  // === example winp file ===
  //BEGIN SURFCLUST
  //--TITLE SURFTEST
  //--ALGO ARVO
  //--SOLVRAD 1.4
  //--SELE /*/*/*
  //--CLUSTBIN 10.0
  //END
  // =========================
  
  inp_surf->iNumOfFrames=iNumOfFrames;
  inp_surf->iFrameNum=-1;
  inp_surf->fSmallestSurf=FLT_MAX;
  inp_surf->fLargestSurf=FLT_MIN;
  
  while( strncmp (cWordomInpBuffer, "END", 3))
  {
    iWinpOptFlag = 0;
    sprintf( cWordomInpBuffer, "%s", input[inp_index]);
    if( !strncmp(cWordomInpBuffer, "BEGIN", 5) || !strncmp(cWordomInpBuffer, "END", 3) || cWordomInpBuffer[0] == '#')
      iWinpOptFlag = 1;
    else if ( !strncmp(cWordomInpBuffer, "--TITLE", 7))
    {
      sscanf( cWordomInpBuffer, "--TITLE %s", inp_surf->cTitle);
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--SOLVRAD", 9))
    {
      sscanf( cWordomInpBuffer, "--SOLVRAD %f", &inp_surf->dSolventRadius);
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--CLUSTBIN", 10))
    {
      sscanf( cWordomInpBuffer, "--CLUSTBIN %f", &inp_surf->fClustBin);
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--ALGO", 6))
    {
      sscanf( cWordomInpBuffer, "--ALGO %s", cTmpString);
      iWinpOptFlag = 1;
      if(!strncmp(cTmpString, "ARVO", 4) || !strncmp(cTmpString, "arvo", 4))
      {
        inp_surf->iAlgoFlag=0;
      }
      else if (!strncmp(cTmpString, "GEPOL", 5) || !strncmp(cTmpString, "gepol", 5))
      {
        inp_surf->iAlgoFlag=1;
      }
      else
      {
        fprintf( stderr, "SURFCLUST module: Invalid --ALGO option: %s\n", cTmpString);
        exit(5);
      }
    }
    else if ( !strncmp(cWordomInpBuffer, "--CALC", 6))
    {
      sscanf( cWordomInpBuffer, "--CALC %s", cTmpString);
      iWinpOptFlag = 1;
      if(!strncmp(cTmpString, "WSURF", 5) || !strncmp(cTmpString, "wsurf", 5))
      {
        iGepolAlgoOptionsCheck=1;
        inp_surf->iSurfType=0;
      }
      else if (!strncmp(cTmpString, "ASURF", 5) || !strncmp(cTmpString, "asurf", 5))
      {
        iGepolAlgoOptionsCheck=1;
        inp_surf->iSurfType=1;
      }
      else if (!strncmp(cTmpString, "ESURF", 5) || !strncmp(cTmpString, "esurf", 5))
      {
        iGepolAlgoOptionsCheck=1;
        inp_surf->iSurfType=2;
      }
      else
      {
        fprintf( stderr, "SURFCLUST module: Invalid --CALC option: %s\n", cTmpString);
        exit(5);
      }
    }
    else if ( !strncmp(cWordomInpBuffer, "--OFAC", 6))
    {
      sscanf( cWordomInpBuffer, "--OFAC %f", &inp_surf->fOFAC);
      if(inp_surf->fOFAC<0)
      {
        fprintf( stderr, "SURFCLUST module: --OFAC must be a real number in the range 0.0..1.0\n");
        exit(5);
      }
      iGepolAlgoOptionsCheck=1;
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--RMIN", 6))
    {
      sscanf( cWordomInpBuffer, "--RMIN %f", &inp_surf->fRMIN);
      if(inp_surf->fRMIN<0)
      {
        fprintf( stderr, "SURF module: --RMIN must be a real number larger than 0.0\n");
        exit(5);
      }
      iGepolAlgoOptionsCheck=1;
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--NDIV", 6))
    {
      sscanf( cWordomInpBuffer, "--NDIV %d", &inp_surf->iNDIV);
      if(inp_surf->iNDIV<0 || inp_surf->iNDIV>5)
      {
        fprintf( stderr, "SURFCLUST module: --NDIV must be an integer number in the range 0..5\n");
        exit(5);
      }
      iGepolAlgoOptionsCheck=1;
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--MEMSIZE", 9))
    {
    //sscanf( cWordomInpBuffer, "--MEMSIZE %d", &iSizeMul);
      sscanf( cWordomInpBuffer, "--MEMSIZE %[^\n]", membuffer);
      memvalue = atof(membuffer);
      if( strchr(membuffer, 'k')!=NULL || strchr(membuffer, 'K')!=NULL )
        memvalue *= 1024;
      else if( strchr(membuffer, 'm')!=NULL || strchr(membuffer, 'M')!=NULL )
        memvalue *= 1024*1024;
      else if( strchr(membuffer, 'g')!=NULL || strchr(membuffer, 'G')!=NULL )
        memvalue *= 1024*1024*1024;
      else
        fprintf( stderr, "You did not specify the memory size unit; taken as bytes\n");
      
      iSizeMul = (int) floor(memvalue/unitsize);
      if(iSizeMul<1)
      {
        fprintf( stderr, "SURF module: --MEMSIZE must be at _least_ >= %d\n", unitsize);
        exit(5);
      }
      iGepolAlgoOptionsCheck=1;
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--SELE", 6))
    {
      sscanf(cWordomInpBuffer, "--SELE %[^\n]%*c ", inp_surf->sele1.selestring);     
      GetSele ( inp_surf->sele1.selestring, &inp_surf->sele1, molecule);
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--RADFILE", 9))
    {
      sscanf( cWordomInpBuffer, "--RADFILE %s", radfilename);
      iRadFileFlag = 1;
      radlist = ReadRadius(radfilename);
      iWinpOptFlag = 1;
    }
    if( iWinpOptFlag==0 )
    {
      fprintf( stderr, "SURFCLUST module: Could NOT understand option: %s\n", cWordomInpBuffer);
      exit(5);
    }
    inp_index++;
  }

  // check if Gepol Options has been set and Arvo algorithm has been chosen
  if(iGepolAlgoOptionsCheck==1 && inp_surf->iAlgoFlag==0)
  {
    fprintf( stderr, "SURFCLUST module:\n");
    fprintf( stderr, "ARVO Algorithm do not use --CALC, --NDIV, --OFAC and --RMIN options\n");
    exit(5);
  }

  // check if OFAC and RMIN has been set in agreement with --CALC ESURF
  if(inp_surf->iSurfType!=2)
  {
    if(inp_surf->fRMIN!=-1 || inp_surf->fOFAC!=-1)
    {
      fprintf( stderr, "SURFCLUST module:\n");
      fprintf( stderr, "--OFAC and --RMIN options can be used only if --ALGO is GEPOL and --CALC is ESURF\n");
      exit(5);
    }
  }

  // Calculates the total number of atoms
  iTotAtoms=0;
  for(ii=0; ii<molecule->nSeg; ii++)
  {
    for(jj=0; jj<molecule->segment[ii].nRpS; jj++)
    {
      iTotAtoms=iTotAtoms+molecule->segment[ii].pRes[jj].nApR;
    }
  }
  
  inp_surf->iNumOfFrames=iNumOfFrames;
  inp_surf->iNumOfRealAtoms=iTotAtoms;                                 // Number of Atoms in PDB file
  inp_surf->iNumOfAtoms=iTotAtoms;
  inp_surf->pfSele1Surf = (float *) calloc (iNumOfFrames, sizeof(float));
  
  if(inp_surf->iAlgoFlag==0)
  {
    /* ============ *\
    |* === ARVO === *|
    \* ============ */
    
    inp_surf->iNumOfAtoms=iTotAtoms;                                   // Number of Atoms within 10 A from selected atoms
    inp_surf->iCirclesPlaneForAtom=iTotAtoms*80;                       // Maximum number of local circles in the plane for one atom
    inp_surf->iArcsAnglesCircleIntersect=iTotAtoms*80;                 // Maximum number of arcs and angles which arise from the local circles intersections
    inp_surf->iNumOfNeighbors=iTotAtoms*80;                            // Maximum number of neighbors
    
    // This vector stores the radii of all atoms
    inp_surf->pdAtomsRadii=(double *) calloc(iTotAtoms, sizeof(double));

    if( iRadFileFlag)
    {
      RadAssign( radlist, molecule, inp_surf->pdAtomsRadii );
      for( ii=0; ii<molecule->nato; ii++)
        inp_surf->pdAtomsRadii[ii] += inp_surf->dSolventRadius;
    }
    else
    {
      mm=-1;
      for(ii=0; ii<molecule->nSeg; ii++)
      {
        for(jj=0; jj<molecule->segment[ii].nRpS; jj++)
        {
          for(kk=0; kk<molecule->segment[ii].pRes[jj].nApR; kk++)
          {
            mm++;
            inp_surf->pdAtomsRadii[mm]=(double)molecule->segment[ii].pRes[jj].pAto[kk].bFac+inp_surf->dSolventRadius;
          }
        }
      }
    }

    /*
      Allocates piNeighborsNumber
      This vector stores the number of neighbors of all atoms
    */
    inp_surf->piNeighborsNumber=(int*) calloc(iTotAtoms*40, sizeof(int));

    /*
      Allocates piNeighborsStartIndex
      This vector stores 
    */
    inp_surf->piNeighborsStartIndex=(int*) calloc(iTotAtoms*40, sizeof(int));

    /*
      Allocates piNeighborsStartIndex
      This vector stores 
    */
    inp_surf->piNeighborsIndices=(int*) calloc(inp_surf->iNumOfNeighbors, sizeof(int));

    /*
      Allocates piIndicesVector
      This vector stores 
    */
    inp_surf->piIndicesVector=(int*) calloc(iTotAtoms*40, sizeof(int));

    /*
      Allocates piIndicesVector2
      This vector stores 
    */
    inp_surf->piIndicesVector2=(int*) calloc(iTotAtoms*40, sizeof(int));


    /*
      Circles Allocation
    */
    inp_surf->ppdCircles=(double **) calloc(inp_surf->iCirclesPlaneForAtom, sizeof(double *));
    for(ii=0; ii<inp_surf->iCirclesPlaneForAtom; ii++)
    {
      inp_surf->ppdCircles[ii]=(double *) calloc(4, sizeof(double));
    }
    
    /*
      Arcs Allocation
    */
    inp_surf->ppdArcs=(double **) calloc(inp_surf->iArcsAnglesCircleIntersect, sizeof(double *));
    for(ii=0; ii<inp_surf->iArcsAnglesCircleIntersect; ii++)
    {
      inp_surf->ppdArcs[ii]=(double *) calloc(3, sizeof(double));
    }

    /*
      NewArcs Allocation
    */
    inp_surf->ppdNewArcs=(double **) calloc(inp_surf->iArcsAnglesCircleIntersect, sizeof(double *));
    for(ii=0; ii<inp_surf->iArcsAnglesCircleIntersect; ii++)
    {
      inp_surf->ppdNewArcs[ii]=(double *) calloc(3, sizeof(double));
    }
    
    /*
      Allocates pdAngles 
    */
    inp_surf->pdAngles=(double*) calloc(inp_surf->iArcsAnglesCircleIntersect, sizeof(double));

    /*
      Allocates pdNewAngles 
    */
    inp_surf->pdNewAngles=(double*) calloc(inp_surf->iArcsAnglesCircleIntersect, sizeof(double));

    /*
      Allocates AtomNeighbors
    */
    inp_surf->ppdAtomNeighbors=(double **) calloc(iTotAtoms, sizeof(double*));
    for(ii=0; ii<iTotAtoms; ii++)
    {
      // dAtomNeighbors dimesion: 0 -> xcoord, 1 -> ycoord, 2 -> zcoord, 3 -> radius
      inp_surf->ppdAtomNeighbors[ii]=(double *) calloc(4, sizeof(double));
    }

    /*
      Allocates piWithinAtoms
    */
    inp_surf->piWithinAtoms=(int *) calloc(inp_surf->iNumOfAtoms, sizeof(int ));

    inp_surf->ppfWithinAtomsCoord=(float **) calloc(inp_surf->iNumOfAtoms, sizeof(float *));
    for(ii=0; ii<inp_surf->iNumOfAtoms; ii++)
    {
      inp_surf->ppfWithinAtomsCoord[ii]=(float *) calloc(4, sizeof(float));
    }
  }
  
  else if(inp_surf->iAlgoFlag==1)
  {
    /* ============= *\
    |* === GEPOL === *|
    \* ============= */
    
    inp_surf->iSolvRadAddedFalg=0;
    inp_surf->iMaxSize=(inp_surf->iNumOfRealAtoms*iSizeMul);    
    inp_surf->pdAtomsRadii=(double *) calloc(inp_surf->iMaxSize, sizeof(double));
    
    // Gets Atom radii
    if( iRadFileFlag)
    {
      RadAssign( radlist, molecule, inp_surf->pdAtomsRadii );
    }
    else
    {
      mm=-1;
      for(ii=0; ii<molecule->nSeg; ii++)
      {
        for(jj=0; jj<molecule->segment[ii].nRpS; jj++)
        {
          for(kk=0; kk<molecule->segment[ii].pRes[jj].nApR; kk++)
          {
            mm++;
            inp_surf->pdAtomsRadii[mm]=(double)molecule->segment[ii].pRes[jj].pAto[kk].bFac;
          }
        }
      }
    }

    if(inp_surf->iAlgoFlag==1 && inp_surf->iSurfType==2)
    {
      if(inp_surf->fRMIN==-1)
        inp_surf->fRMIN=0.5;
        
      if(inp_surf->fOFAC==-1)
        inp_surf->fOFAC=0.8;
    }
    
    // Allocates ppfAtomCoord
    inp_surf->ppfAtomCoord=(float **) calloc(inp_surf->iMaxSize, sizeof(float *));
    for(ii=0; ii<inp_surf->iMaxSize; ii++)
    {
      inp_surf->ppfAtomCoord[ii]=(float *) calloc(3, sizeof(float));
    }

    // Allocates piAtomType
    inp_surf->piAtomType=(int *) calloc(inp_surf->iMaxSize, sizeof(int));
    
    // Fills JVT1 & JVT2 arrays
    FillJVTArrays(inp_surf);
     
    // Allocates pfXCoord1, pfYCoord1, pfZCoord1
    inp_surf->iNumOfTesserae=60*pow(4, (inp_surf->iNDIV-1));
    inp_surf->pfXCoord1=(float *) calloc(inp_surf->iNumOfTesserae, sizeof(float));
    inp_surf->pfYCoord1=(float *) calloc(inp_surf->iNumOfTesserae, sizeof(float));
    inp_surf->pfZCoord1=(float *) calloc(inp_surf->iNumOfTesserae, sizeof(float));
    
    CompTriaVertexCoord(inp_surf);
    
    Divide(inp_surf);
    
    if(inp_surf->iSurfType==2)
    {
      // Allocates pfDistVect vector, used in Bulk
      inp_surf->pfDistVect=(float *) calloc(inp_surf->iMaxSize, sizeof(float));
      
      // Allocates piEngulfedAtoms
      inp_surf->piEngulfedAtoms=(int *) calloc(inp_surf->iMaxSize, sizeof(int));
      
      // Allocates piZeroSurf, piZeroSurfMZ5 and piZeroSurfC5
      inp_surf->piZeroSurf    =(int *) calloc(inp_surf->iMaxSize, sizeof(int));
      inp_surf->piZeroSurfMZ5 =(int *) calloc(inp_surf->iMaxSize, sizeof(int));
      inp_surf->piZeroSurfC5  =(int *) calloc(inp_surf->iMaxSize, sizeof(int));
      
      // Allocates piUsefulAtoms
      inp_surf->piUsefulAtoms = (int *) calloc(inp_surf->iMaxSize, sizeof(int));

      // Allocates piZeroCheck
      inp_surf->piZeroCheck=(int *) calloc(inp_surf->iMaxSize, sizeof(int));
    }
        
    // Allocates piIndexVect, used in GeoCav
    inp_surf->piIndexVect = (int *) calloc(inp_surf->iMaxSize, sizeof(int));
    
    // Allocated ppfPCoord
    inp_surf->ppfPCoord = (float **) calloc(inp_surf->iMaxSize, sizeof(float *));
    for(ii=0; ii<inp_surf->iMaxSize; ii++)
      inp_surf->ppfPCoord[ii] = (float *) calloc(3, sizeof(float));

    // Allocates pfAPVect
    inp_surf->pfAPVect = (float *) calloc(inp_surf->iMaxSize, sizeof(float));

  }
  
  if(inp_surf->iWriteSurfFlag==1)
  {
    /*fprintf(OutputFile, " %15s ", inp_surf->cTitle);*/
    sprintf( outstring, " %15s ", inp_surf->cTitle);
    return 17;
  }
  return 0;
}

// =====================================================================
int Compute_SURFCLUST(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring) //FILE *OutputFile)
{
  if(inp_surf->iAlgoFlag==0)
    Compute_ARVO_SURFCLUST(inp_surf, molecule, trj_crd, outstring);
    
  else if(inp_surf->iAlgoFlag==1)
    Compute_GEPOL_SURFCLUST(inp_surf, molecule, trj_crd, outstring);
  
  return 17;  
}

void Compute_ARVO_SURFCLUST(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring) //FILE *OutputFile)
{
  int     ii, jj, iNumOfWithinAtoms, iRotCont=0;
  float   fDist;
  double  dSurf1=0.0, dTmpSurf1=0.0;

  for(ii=0; ii<inp_surf->iNumOfRealAtoms; ii++)
  {
    inp_surf->piWithinAtoms[ii]=0;
  }

  iNumOfWithinAtoms=-1;
  for(ii=0; ii<inp_surf->sele1.nselatm; ii++)
  {
    iNumOfWithinAtoms++;
    inp_surf->piWithinAtoms[inp_surf->sele1.selatm[ii]-1]=1;
    inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][0]=trj_crd->xcoor[inp_surf->sele1.selatm[ii]-1];
    inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][1]=trj_crd->ycoor[inp_surf->sele1.selatm[ii]-1];
    inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][2]=trj_crd->zcoor[inp_surf->sele1.selatm[ii]-1];
    inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][3]=inp_surf->pdAtomsRadii[inp_surf->sele1.selatm[ii]-1];   
  }
  for(ii=0; ii<inp_surf->sele1.nselatm; ii++)
  {
    for(jj=0; jj<inp_surf->iNumOfRealAtoms; jj++)
    {
      
      if(inp_surf->piWithinAtoms[jj]==0)
      {
        fDist=sqrt(
              pow((trj_crd->xcoor[inp_surf->sele1.selatm[ii]-1]-trj_crd->xcoor[jj]), 2)+
              pow((trj_crd->ycoor[inp_surf->sele1.selatm[ii]-1]-trj_crd->ycoor[jj]), 2)+
              pow((trj_crd->zcoor[inp_surf->sele1.selatm[ii]-1]-trj_crd->zcoor[jj]), 2));

        if(fDist<=50)
        {
          iNumOfWithinAtoms++;
          inp_surf->piWithinAtoms[jj]=1;
          inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][0]=trj_crd->xcoor[jj];
          inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][1]=trj_crd->ycoor[jj];
          inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][2]=trj_crd->zcoor[jj];
          inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][3]=inp_surf->pdAtomsRadii[jj];
        }
      }
    }
  }

  inp_surf->iNumOfAtoms=iNumOfWithinAtoms+1;
  
  MakeNeighbors(inp_surf); 
  iRotCont=0;
  while(NorthPoleTest(inp_surf)==0)
  {
    iRotCont++;
    MolRotation(inp_surf);
  }
  
  for(ii=0; ii<inp_surf->sele1.nselatm; ii++)
  {
    dTmpSurf1=GetSurf(ii, inp_surf);
    if(isnan(dTmpSurf1)==1)
    {
      dTmpSurf1=0.0;
    }
    dSurf1=dSurf1+dTmpSurf1;
  }
  
  if(dSurf1>inp_surf->fLargestSurf)
    inp_surf->fLargestSurf=dSurf1;
  
  if(dSurf1<inp_surf->fSmallestSurf)
    inp_surf->fSmallestSurf=dSurf1;
  
  inp_surf->iFrameNum++;
  inp_surf->pfSele1Surf[inp_surf->iFrameNum]=dSurf1;
}

void Compute_GEPOL_SURFCLUST(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring) //FILE *OutputFile)
{
  int     ii;
  double  dSurf=0.0;
    
  // Set inp_surf->iNumOfAtoms to inp_surf->iNumOfRealAtoms
  // and copy frame coords to inp_surf->ppfAtomCoord
  inp_surf->iNumOfAtoms=inp_surf->iNumOfRealAtoms;
  for(ii=0; ii<inp_surf->iNumOfRealAtoms; ii++)
  {
    inp_surf->ppfAtomCoord[ii][0]=trj_crd->xcoor[ii];
    inp_surf->ppfAtomCoord[ii][1]=trj_crd->ycoor[ii];
    inp_surf->ppfAtomCoord[ii][2]=trj_crd->zcoor[ii];
  }

  // Set all atoms as ghost
  for(ii=0; ii<inp_surf->iMaxSize; ii++)
    inp_surf->piAtomType[ii]=3;
    
  // Set selected atoms as 1 or 6 in accordance with radii
  for(ii=0; ii<inp_surf->sele1.nselatm; ii++)
  {
    if(inp_surf->pdAtomsRadii[inp_surf->sele1.selatm[ii]-1]==0)
      inp_surf->piAtomType[inp_surf->sele1.selatm[ii]-1]=1;
    
    else
      inp_surf->piAtomType[inp_surf->sele1.selatm[ii]-1]=6;
  }

  if(inp_surf->sele1.nselatm<inp_surf->iNumOfRealAtoms)
    inp_surf->iGhostFlag=1;
  
  // === Surface Calculation ==========================
  if(inp_surf->iSurfType==2)
  {
    if(inp_surf->iGhostFlag==1)
      Shell(inp_surf);
      
    Bulk(inp_surf);
    Clean5(inp_surf);
    GenAtoms(inp_surf);
    Clean5(inp_surf);
  }
  
  if(inp_surf->iSurfType==1 && inp_surf->iSolvRadAddedFalg==0)
  {
    inp_surf->iSolvRadAddedFalg=1;
    AddSolvRad(inp_surf);
  }
  
  GeoCav(inp_surf);
  dSurf=AreaSum(inp_surf);
  
  if(dSurf>inp_surf->fLargestSurf)
    inp_surf->fLargestSurf=dSurf;
  
  if(dSurf<inp_surf->fSmallestSurf)
    inp_surf->fSmallestSurf=dSurf;
    
  inp_surf->iFrameNum++;
  inp_surf->pfSele1Surf[inp_surf->iFrameNum]=dSurf;  
}

int Post_SURFCLUST(struct inp_surf *inp_surf, struct sopt *OPT)
{
  int          ii, jj;
  
  int          iNumOfClust=0;
  int          iClustPop=0;
  int         *piAssignedFrames; 
  float        fSurfLowBound=0.0;
  float        fSurfUpBound=0.0;
  
  char         cOutputFileName[80];
  time_t       time_Today;
  FILE        *OutputFile;
  
  piAssignedFrames=(int *) calloc((inp_surf->iFrameNum+1), sizeof(int));
  
  time(&time_Today);
  sprintf(cOutputFileName, "%s%s", inp_surf->cTitle, ".surfclust");
  OutputFile=O_File(cOutputFileName, "w");

  iNumOfClust=(int)(((inp_surf->fLargestSurf-inp_surf->fSmallestSurf)/inp_surf->fClustBin)+0.5);
  if(iNumOfClust==0)
    iNumOfClust=1;
    
  fprintf(OutputFile, "=========================================\n");
  fprintf(OutputFile, "*** WORDOM SURFace  CLUSTering Module ***\n");
  fprintf(OutputFile, "=========================================\n");
  fprintf(OutputFile, "Copyright      : Francesca Fanelli\n");
  fprintf(OutputFile, "                 Angelo Felline\n");
  fprintf(OutputFile, "                 (2009)\n");
  fprintf(OutputFile, "License        : GPL v 3\n");
  fprintf(OutputFile, "\n");
  fprintf(OutputFile, "Date           : %s", asctime(localtime(&time_Today)));
  fprintf(OutputFile, "\n");
  fprintf(OutputFile, "Files\n-----\n");
  fprintf(OutputFile, "Mol Name       : %s\n", OPT->IMOL_FILE);
  fprintf(OutputFile, "Total Atoms    : %d\n", inp_surf->iNumOfRealAtoms);
  fprintf(OutputFile, "Trj Name       : %s\n", OPT->ITRJ_FILE);
  fprintf(OutputFile, "Total Frames   : %d\n", inp_surf->iNumOfFrames);
  fprintf(OutputFile, "Used Frames    : %d\n", inp_surf->iFrameNum+1);
  fprintf(OutputFile, "\n");
  fprintf(OutputFile, "Options\n-------\n");
  if(inp_surf->iAlgoFlag==0)
    fprintf(OutputFile, "Algorithm      : ARVO\n");
  else if(inp_surf->iAlgoFlag==1)
  {
    fprintf(OutputFile, "Algorithm      : GEPOL\n");
    if(inp_surf->iSurfType==0)
      fprintf(OutputFile, "Surface        : van der Waals Surface Area\n");
      
    else if(inp_surf->iSurfType==1)
      fprintf(OutputFile, "Surface        : Accessible Surface Area\n");
      
    else if(inp_surf->iSurfType==2)
      fprintf(OutputFile, "Surface        : Excluded Surface Area\n");
      
    fprintf(OutputFile, "NDIV           : %d\n", inp_surf->iNDIV);
    if(inp_surf->iSurfType==2)
    {
      fprintf(OutputFile, "OFAC           : %f\n", inp_surf->fOFAC);
      fprintf(OutputFile, "RMIN           : %f\n", inp_surf->fRMIN);
    }
  }
  fprintf(OutputFile, "Solvent Radius : %f\n", inp_surf->dSolventRadius);
  fprintf(OutputFile, "Sele           : %s\n", inp_surf->sele1.selestring);  
  fprintf(OutputFile, "Sele Atoms     : %d\n", inp_surf->sele1.nselatm);
  fprintf(OutputFile, "Cluster Bin    : %f\n", inp_surf->fClustBin);
  fprintf(OutputFile, "\n");
  fprintf(OutputFile, "Summary\n-------\n");
  fprintf(OutputFile, "Smallest Surf  : %f\n", inp_surf->fSmallestSurf);
  fprintf(OutputFile, "Largest Surf   : %f\n", inp_surf->fLargestSurf);
  fprintf(OutputFile, "Num of Clust   : %d\n", iNumOfClust);
  fprintf(OutputFile, "\n");
  fprintf(OutputFile, "Clusters\n--------\n");
  
  for(ii=0; ii<iNumOfClust; ii++)
  {
    iClustPop=0;
    fSurfLowBound=(inp_surf->fSmallestSurf+((ii)*inp_surf->fClustBin));
    fSurfUpBound=(inp_surf->fSmallestSurf+((ii+1)*inp_surf->fClustBin));
    fprintf(OutputFile, "Cluster #%d (from %.3f to %.3f)\n", ii+1, fSurfLowBound, fSurfUpBound);
    
    for(jj=0; jj<(inp_surf->iFrameNum+1); jj++)
    {
      if(piAssignedFrames[jj]==0)
      {
        if(inp_surf->pfSele1Surf[jj]>=fSurfLowBound && inp_surf->pfSele1Surf[jj]<fSurfUpBound)
        {
          fprintf(OutputFile, "%-5d   %15.3f\n", jj+1, inp_surf->pfSele1Surf[jj]);
          piAssignedFrames[jj]=1;
          iClustPop++;
        }
      }
    }
    fprintf(OutputFile, "Total Cluster Population: %d\n\n", iClustPop);
  }
  
  return 0;
}
// =====================================================================
int Read_SURFCORR(char **input, int inp_index, struct inp_surf *inp_surf, Molecule *molecule, char *outstring, int iNumOfFrames) //FILE *OutputFile, int iNumOfFrames)
{
  int     ii, jj, kk, mm;
  int     iGepolAlgoOptionsCheck=0;
  int     iSizeMul=20;
  int     iWinpOptFlag=0, iTotAtoms=0;
  char    cWordomInpBuffer[1024];
  char    cTmpString[6];
  int     iRadFileFlag=0;
  char    radfilename[128];
  RadList *radlist;

  int     unitsize;
  char    membuffer[128];
  float   memvalue;
  
  extern short int     frame_par;
  frame_par = 1;

  srand(time(NULL));
 
  // === Default Values ===========
  sprintf(inp_surf->cTitle, "SURFCORR");
  inp_surf->dSolventRadius=1.4;
  inp_surf->iWriteSurfFlag=1;
  inp_surf->iWriteRawFlag=0;
  inp_surf->iAlgoFlag=0;
  inp_surf->iSurfType=1;
  inp_surf->iNDIV=3;
  inp_surf->fOFAC=-1;
  inp_surf->fRMIN=-1;
  inp_surf->iGhostFlag=0;
  // ==============================
    
  unitsize = molecule->nato * ( sizeof(double) + 2*sizeof(float **) + 8*sizeof(int) + 8*sizeof(float) );
    
  // === example winp file ===
  //BEGIN
  //--TITLE SURFTEST
  //--ALGO GEPOL
  //--CALC ASURF
  //--NDIV 5
  //--SOLVRAD 1.4
  //--SELE1 /*/*/*
  //--SELE2 /*/*/*
  //--MEMSIZE 10
  //END
  // =========================
  
  inp_surf->iNumOfFrames=iNumOfFrames;
  inp_surf->iFrameNum=-1;
  
  while( strncmp (cWordomInpBuffer, "END", 3))
  {
    iWinpOptFlag = 0;
    sprintf( cWordomInpBuffer, "%s", input[inp_index]);
    if( !strncmp(cWordomInpBuffer, "BEGIN", 5) || !strncmp(cWordomInpBuffer, "END", 3) || cWordomInpBuffer[0] == '#')
      iWinpOptFlag = 1;
    else if ( !strncmp(cWordomInpBuffer, "--TITLE", 7))
    {
      sscanf( cWordomInpBuffer, "--TITLE %s", inp_surf->cTitle);
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--SOLVRAD", 9))
    {
      sscanf( cWordomInpBuffer, "--SOLVRAD %f", &inp_surf->dSolventRadius);
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--VERBOSE", 9))
    {
      // sscanf( cWordomInpBuffer, "--VERBOSE %f", &inp_surf->dSolventRadius);
      inp_surf->iWriteRawFlag=1;
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--ALGO", 6))
    {
      sscanf( cWordomInpBuffer, "--ALGO %s", cTmpString);
      iWinpOptFlag = 1;
      if(!strncmp(cTmpString, "ARVO", 4) || !strncmp(cTmpString, "arvo", 4))
      {
        inp_surf->iAlgoFlag=0;
      }
      else if (!strncmp(cTmpString, "GEPOL", 5) || !strncmp(cTmpString, "gepol", 5))
      {
        inp_surf->iAlgoFlag=1;
      }
      else
      {
        fprintf( stderr, "SURFCORR module: Invalid --ALGO option: %s\n", cTmpString);
        exit(5);
      }
    }
    else if ( !strncmp(cWordomInpBuffer, "--CALC", 6))
    {
      sscanf( cWordomInpBuffer, "--CALC %s", cTmpString);
      iWinpOptFlag = 1;
      if(!strncmp(cTmpString, "WSURF", 5) || !strncmp(cTmpString, "wsurf", 5))
      {
        iGepolAlgoOptionsCheck=1;
        inp_surf->iSurfType=0;
      }
      else if (!strncmp(cTmpString, "ASURF", 5) || !strncmp(cTmpString, "asurf", 5))
      {
        iGepolAlgoOptionsCheck=1;
        inp_surf->iSurfType=1;
      }
      else if (!strncmp(cTmpString, "ESURF", 5) || !strncmp(cTmpString, "esurf", 5))
      {
        iGepolAlgoOptionsCheck=1;
        inp_surf->iSurfType=2;
      }
      else
      {
        fprintf( stderr, "SURFCORR module: Invalid --CALC option: %s\n", cTmpString);
        exit(5);
      }
    }
    else if ( !strncmp(cWordomInpBuffer, "--OFAC", 6))
    {
      sscanf( cWordomInpBuffer, "--OFAC %f", &inp_surf->fOFAC);
      if(inp_surf->fOFAC<0)
      {
        fprintf( stderr, "SURF module: --OFAC must be a real number in the range 0.0..1.0\n");
        exit(5);
      }
      iGepolAlgoOptionsCheck=1;
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--RMIN", 6))
    {
      sscanf( cWordomInpBuffer, "--RMIN %f", &inp_surf->fRMIN);
      if(inp_surf->fRMIN<0)
      {
        fprintf( stderr, "SURFCORR module: --RMIN must be a real number larger than 0.0\n");
        exit(5);
      }
      iGepolAlgoOptionsCheck=1;
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--NDIV", 6))
    {
      sscanf( cWordomInpBuffer, "--NDIV %d", &inp_surf->iNDIV);
      if(inp_surf->iNDIV<0 || inp_surf->iNDIV>5)
      {
        fprintf( stderr, "SURF module: --NDIV must be an integer number in the range 0..5\n");
        exit(5);
      }
      iGepolAlgoOptionsCheck=1;
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--MEMSIZE", 9))
    {
    //sscanf( cWordomInpBuffer, "--MEMSIZE %d", &iSizeMul);
      sscanf( cWordomInpBuffer, "--MEMSIZE %[^\n]", membuffer);
      memvalue = atof(membuffer);
      if( strchr(membuffer, 'k')!=NULL || strchr(membuffer, 'K')!=NULL )
        memvalue *= 1024;
      else if( strchr(membuffer, 'm')!=NULL || strchr(membuffer, 'M')!=NULL )
        memvalue *= 1024*1024;
      else if( strchr(membuffer, 'g')!=NULL || strchr(membuffer, 'G')!=NULL )
        memvalue *= 1024*1024*1024;
      else
        fprintf( stderr, "You did not specify the memory size unit; taken as bytes\n");
      
      iSizeMul = (int) floor(memvalue/unitsize);
      if(iSizeMul<1)
      {
        fprintf( stderr, "SURF module: --MEMSIZE must be at _least_ >= %d\n", unitsize);
        exit(5);
      }
      iGepolAlgoOptionsCheck=1;
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--SELE1", 7))
    {
      sscanf(cWordomInpBuffer, "--SELE1 %[^\n]%*c ", inp_surf->sele1.selestring);     
      GetSele ( inp_surf->sele1.selestring, &inp_surf->sele1, molecule);
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--SELE2", 7))
    {
      sscanf(cWordomInpBuffer, "--SELE2 %[^\n]%*c ", inp_surf->sele2.selestring);
      GetSele ( inp_surf->sele2.selestring, &inp_surf->sele2, molecule);      
      iWinpOptFlag = 1;
    }
    else if ( !strncmp(cWordomInpBuffer, "--RADFILE", 9))
    {
      sscanf( cWordomInpBuffer, "--RADFILE %s", radfilename);
      iRadFileFlag = 1;
      radlist = ReadRadius(radfilename);
      iWinpOptFlag = 1;
    }
    if( iWinpOptFlag==0 )
    {
      fprintf( stderr, "SURFCORR module: Could NOT understand option: %s\n", cWordomInpBuffer);
      exit(5);
    }
    inp_index++;
  }

  // check if Gepol Options has been set and Arvo algorithm has been chosen
  if(iGepolAlgoOptionsCheck==1 && inp_surf->iAlgoFlag==0)
  {
    fprintf( stderr, "SURFCORR module:\n");
    fprintf( stderr, "ARVO Algorithm do not use --CALC, --NDIV, --OFAC and --RMIN options\n");
    exit(5);
  }

  // check if OFAC and RMIN has been set in agreement with --CALC ESURF
  if(inp_surf->iSurfType!=2)
  {
    if(inp_surf->fRMIN!=-1 || inp_surf->fOFAC!=-1)
    {
      fprintf( stderr, "SURFCORR module:\n");
      fprintf( stderr, "--OFAC and --RMIN options can be used only if --ALGO is GEPOL and --CALC is ESURF\n");
      exit(5);
    }
  }

  // Calculates the total number of atoms
  iTotAtoms=0;
  for(ii=0; ii<molecule->nSeg; ii++)
  {
    for(jj=0; jj<molecule->segment[ii].nRpS; jj++)
    {
      iTotAtoms=iTotAtoms+molecule->segment[ii].pRes[jj].nApR;
    }
  }
  
  inp_surf->iNumOfRealAtoms=iTotAtoms;                                 // Number of Atoms in PDB file
  // Allocates pfSele1Surf and pfSele2Surf
  inp_surf->pfSele1Surf = (float *) calloc (iNumOfFrames, sizeof(float));
  inp_surf->pfSele2Surf = (float *) calloc (iNumOfFrames, sizeof(float));
  
  if(inp_surf->iAlgoFlag==0)
  {
    /* ============ *\
    |* === ARVO === *|
    \* ============ */
    
    inp_surf->iNumOfAtoms=iTotAtoms;                                   // Number of Atoms within 10 A from selected atoms
    inp_surf->iCirclesPlaneForAtom=iTotAtoms*80;                       // Maximum number of local circles in the plane for one atom
    inp_surf->iArcsAnglesCircleIntersect=iTotAtoms*80;                 // Maximum number of arcs and angles which arise from the local circles intersections
    inp_surf->iNumOfNeighbors=iTotAtoms*80;                            // Maximum number of neighbors
    
    // This vector stores the radii of all atoms
    inp_surf->pdAtomsRadii=(double *) calloc(iTotAtoms, sizeof(double));

    if( iRadFileFlag)
    {
      RadAssign( radlist, molecule, inp_surf->pdAtomsRadii );
      for( ii=0; ii<molecule->nato; ii++)
        inp_surf->pdAtomsRadii[ii] += inp_surf->dSolventRadius;
    }
    else
    {
      mm=-1;
      for(ii=0; ii<molecule->nSeg; ii++)
      {
        for(jj=0; jj<molecule->segment[ii].nRpS; jj++)
        {
          for(kk=0; kk<molecule->segment[ii].pRes[jj].nApR; kk++)
          {
            mm++;
            inp_surf->pdAtomsRadii[mm]=(double)molecule->segment[ii].pRes[jj].pAto[kk].bFac+inp_surf->dSolventRadius;
          }
        }
      }
    }

    /*
      Allocates piNeighborsNumber
      This vector stores the number of neighbors of all atoms
    */
    inp_surf->piNeighborsNumber=(int*) calloc(iTotAtoms*40, sizeof(int));

    /*
      Allocates piNeighborsStartIndex
      This vector stores 
    */
    inp_surf->piNeighborsStartIndex=(int*) calloc(iTotAtoms*40, sizeof(int));

    /*
      Allocates piNeighborsStartIndex
      This vector stores 
    */
    inp_surf->piNeighborsIndices=(int*) calloc(inp_surf->iNumOfNeighbors, sizeof(int));

    /*
      Allocates piIndicesVector
      This vector stores 
    */
    inp_surf->piIndicesVector=(int*) calloc(iTotAtoms*40, sizeof(int));

    /*
      Allocates piIndicesVector2
      This vector stores 
    */
    inp_surf->piIndicesVector2=(int*) calloc(iTotAtoms*40, sizeof(int));


    /*
      Circles Allocation
    */
    inp_surf->ppdCircles=(double **) calloc(inp_surf->iCirclesPlaneForAtom, sizeof(double *));
    for(ii=0; ii<inp_surf->iCirclesPlaneForAtom; ii++)
    {
      inp_surf->ppdCircles[ii]=(double *) calloc(4, sizeof(double));
    }
    
    /*
      Arcs Allocation
    */
    inp_surf->ppdArcs=(double **) calloc(inp_surf->iArcsAnglesCircleIntersect, sizeof(double *));
    for(ii=0; ii<inp_surf->iArcsAnglesCircleIntersect; ii++)
    {
      inp_surf->ppdArcs[ii]=(double *) calloc(3, sizeof(double));
    }

    /*
      NewArcs Allocation
    */
    inp_surf->ppdNewArcs=(double **) calloc(inp_surf->iArcsAnglesCircleIntersect, sizeof(double *));
    for(ii=0; ii<inp_surf->iArcsAnglesCircleIntersect; ii++)
    {
      inp_surf->ppdNewArcs[ii]=(double *) calloc(3, sizeof(double));
    }
    
    /*
      Allocates pdAngles 
    */
    inp_surf->pdAngles=(double*) calloc(inp_surf->iArcsAnglesCircleIntersect, sizeof(double));

    /*
      Allocates pdNewAngles 
    */
    inp_surf->pdNewAngles=(double*) calloc(inp_surf->iArcsAnglesCircleIntersect, sizeof(double));

    /*
      Allocates AtomNeighbors
    */
    inp_surf->ppdAtomNeighbors=(double **) calloc(iTotAtoms, sizeof(double*));
    for(ii=0; ii<iTotAtoms; ii++)
    {
      // dAtomNeighbors dimesion: 0 -> xcoord, 1 -> ycoord, 2 -> zcoord, 3 -> radius
      inp_surf->ppdAtomNeighbors[ii]=(double *) calloc(4, sizeof(double));
    }

    /*
      Allocates piWithinAtoms
    */
    inp_surf->piWithinAtoms=(int *) calloc(inp_surf->iNumOfAtoms, sizeof(int ));

    inp_surf->ppfWithinAtomsCoord=(float **) calloc(inp_surf->iNumOfAtoms, sizeof(float *));
    for(ii=0; ii<inp_surf->iNumOfAtoms; ii++)
    {
      inp_surf->ppfWithinAtomsCoord[ii]=(float *) calloc(4, sizeof(float));
    }
  }
  
  else if(inp_surf->iAlgoFlag==1)
  {
    /* ============= *\
    |* === GEPOL === *|
    \* ============= */
    
    inp_surf->iSolvRadAddedFalg=0;
    inp_surf->iMaxSize=(inp_surf->iNumOfRealAtoms*iSizeMul);    
    inp_surf->pdAtomsRadii=(double *) calloc(inp_surf->iMaxSize, sizeof(double));
    
    // Gets Atom radii
    if( iRadFileFlag)
    {
      RadAssign( radlist, molecule, inp_surf->pdAtomsRadii );
    }
    else
    {
      mm=-1;
      for(ii=0; ii<molecule->nSeg; ii++)
      {
        for(jj=0; jj<molecule->segment[ii].nRpS; jj++)
        {
          for(kk=0; kk<molecule->segment[ii].pRes[jj].nApR; kk++)
          {
            mm++;
            inp_surf->pdAtomsRadii[mm]=(double)molecule->segment[ii].pRes[jj].pAto[kk].bFac;
          }
        }
      }
    }
    
    if(inp_surf->iAlgoFlag==1 && inp_surf->iSurfType==2)
    {
      if(inp_surf->fRMIN==-1)
        inp_surf->fRMIN=0.5;
        
      if(inp_surf->fOFAC==-1)
        inp_surf->fOFAC=0.8;
    }
    
    // Allocates ppfAtomCoord
    inp_surf->ppfAtomCoord=(float **) calloc(inp_surf->iMaxSize, sizeof(float *));
    for(ii=0; ii<inp_surf->iMaxSize; ii++)
    {
      inp_surf->ppfAtomCoord[ii]=(float *) calloc(3, sizeof(float));
    }

    // Allocates piAtomType
    inp_surf->piAtomType=(int *) calloc(inp_surf->iMaxSize, sizeof(int));
    
    // Fills JVT1 & JVT2 arrays
    FillJVTArrays(inp_surf);
     
    // Allocates pfXCoord1, pfYCoord1, pfZCoord1
    inp_surf->iNumOfTesserae=60*pow(4, (inp_surf->iNDIV-1));
    inp_surf->pfXCoord1=(float *) calloc(inp_surf->iNumOfTesserae, sizeof(float));
    inp_surf->pfYCoord1=(float *) calloc(inp_surf->iNumOfTesserae, sizeof(float));
    inp_surf->pfZCoord1=(float *) calloc(inp_surf->iNumOfTesserae, sizeof(float));
    
    CompTriaVertexCoord(inp_surf);
    
    Divide(inp_surf);
    
    if(inp_surf->iSurfType==2)
    {
      // Allocates pfDistVect vector, used in Bulk
      inp_surf->pfDistVect=(float *) calloc(inp_surf->iMaxSize, sizeof(float));
      
      // Allocates piEngulfedAtoms
      inp_surf->piEngulfedAtoms=(int *) calloc(inp_surf->iMaxSize, sizeof(int));
      
      // Allocates piZeroSurf, piZeroSurfMZ5 and piZeroSurfC5
      inp_surf->piZeroSurf    =(int *) calloc(inp_surf->iMaxSize, sizeof(int));
      inp_surf->piZeroSurfMZ5 =(int *) calloc(inp_surf->iMaxSize, sizeof(int));
      inp_surf->piZeroSurfC5  =(int *) calloc(inp_surf->iMaxSize, sizeof(int));
      
      // Allocates piUsefulAtoms
      inp_surf->piUsefulAtoms = (int *) calloc(inp_surf->iMaxSize, sizeof(int));

      // Allocates piZeroCheck
      inp_surf->piZeroCheck=(int *) calloc(inp_surf->iMaxSize, sizeof(int));
    }
        
    // Allocates piIndexVect, used in GeoCav
    inp_surf->piIndexVect = (int *) calloc(inp_surf->iMaxSize, sizeof(int));
    
    // Allocated ppfPCoord
    inp_surf->ppfPCoord = (float **) calloc(inp_surf->iMaxSize, sizeof(float *));
    for(ii=0; ii<inp_surf->iMaxSize; ii++)
      inp_surf->ppfPCoord[ii] = (float *) calloc(3, sizeof(float));

    // Allocates pfAPVect
    inp_surf->pfAPVect = (float *) calloc(inp_surf->iMaxSize, sizeof(float));

  }
  
  /*fprintf(OutputFile, " %15s ", inp_surf->cTitle);*/
  sprintf( outstring, " %15s ", inp_surf->cTitle);
  
  return 17;
}

// =====================================================================
int Compute_SURFCORR(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring) //FILE *OutputFile)
{
  
  if(inp_surf->iAlgoFlag==0)
    Compute_ARVO_SURFCORR(inp_surf, molecule, trj_crd, outstring);

  else if(inp_surf->iAlgoFlag==1)
    Compute_GEPOL_SURFCORR(inp_surf, molecule, trj_crd, outstring);

  return 17;
}

void Compute_ARVO_SURFCORR(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring) //FILE *OutputFile)
{
  int     ii, jj, iNumOfWithinAtoms, iRotCont=0;
  float   fDist;
  double  dSurf1=0.0, dTmpSurf1=0.0;
  double  dSurf2=0.0, dTmpSurf2=0.0;

  for(ii=0; ii<inp_surf->iNumOfRealAtoms; ii++)
  {
    inp_surf->piWithinAtoms[ii]=0;
  }

  iNumOfWithinAtoms=-1;
  for(ii=0; ii<inp_surf->sele1.nselatm; ii++)
  {
    iNumOfWithinAtoms++;
    inp_surf->piWithinAtoms[inp_surf->sele1.selatm[ii]-1]=1;
    inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][0]=trj_crd->xcoor[inp_surf->sele1.selatm[ii]-1];
    inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][1]=trj_crd->ycoor[inp_surf->sele1.selatm[ii]-1];
    inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][2]=trj_crd->zcoor[inp_surf->sele1.selatm[ii]-1];
    inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][3]=inp_surf->pdAtomsRadii[inp_surf->sele1.selatm[ii]-1];   
  }

  for(ii=0; ii<inp_surf->sele1.nselatm; ii++)
  {
    for(jj=0; jj<inp_surf->iNumOfRealAtoms; jj++)
    {
      
      if(inp_surf->piWithinAtoms[jj]==0)
      {
        fDist=sqrt(
              pow((trj_crd->xcoor[inp_surf->sele1.selatm[ii]-1]-trj_crd->xcoor[jj]), 2)+
              pow((trj_crd->ycoor[inp_surf->sele1.selatm[ii]-1]-trj_crd->ycoor[jj]), 2)+
              pow((trj_crd->zcoor[inp_surf->sele1.selatm[ii]-1]-trj_crd->zcoor[jj]), 2));
        if(fDist<=50)
        {
          iNumOfWithinAtoms++;
          inp_surf->piWithinAtoms[jj]=1;
          inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][0]=trj_crd->xcoor[jj];
          inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][1]=trj_crd->ycoor[jj];
          inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][2]=trj_crd->zcoor[jj];
          inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][3]=inp_surf->pdAtomsRadii[jj];
        }
      }
    }
  }

  inp_surf->iNumOfAtoms=iNumOfWithinAtoms+1;
  
  MakeNeighbors(inp_surf);
  
  iRotCont=0;
  while(NorthPoleTest(inp_surf)==0)
  {
    iRotCont++;
    MolRotation(inp_surf);
  }

  for(ii=0; ii<inp_surf->sele1.nselatm; ii++)
  {
    dTmpSurf1=GetSurf(ii, inp_surf);
    if(isnan(dTmpSurf1)==1)
    {
      dTmpSurf1=0.0;
    }
    dSurf1=dSurf1+dTmpSurf1;
  }

  for(ii=0; ii<inp_surf->iNumOfRealAtoms; ii++)
  {
    inp_surf->piWithinAtoms[ii]=0;
  }

  iNumOfWithinAtoms=-1;
  for(ii=0; ii<inp_surf->sele2.nselatm; ii++)
  {
    iNumOfWithinAtoms++;
    inp_surf->piWithinAtoms[inp_surf->sele2.selatm[ii]-1]=1;
    inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][0]=trj_crd->xcoor[inp_surf->sele2.selatm[ii]-1];
    inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][1]=trj_crd->ycoor[inp_surf->sele2.selatm[ii]-1];
    inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][2]=trj_crd->zcoor[inp_surf->sele2.selatm[ii]-1];
    inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][3]=inp_surf->pdAtomsRadii[inp_surf->sele2.selatm[ii]-1];   
  }

  for(ii=0; ii<inp_surf->sele2.nselatm; ii++)
  {
    for(jj=0; jj<inp_surf->iNumOfRealAtoms; jj++)
    {
      
      if(inp_surf->piWithinAtoms[jj]==0)
      {
        fDist=sqrt(
              pow((trj_crd->xcoor[inp_surf->sele2.selatm[ii]-1]-trj_crd->xcoor[jj]), 2)+
              pow((trj_crd->ycoor[inp_surf->sele2.selatm[ii]-1]-trj_crd->ycoor[jj]), 2)+
              pow((trj_crd->zcoor[inp_surf->sele2.selatm[ii]-1]-trj_crd->zcoor[jj]), 2));
        if(fDist<=50)
        {
          iNumOfWithinAtoms++;
          inp_surf->piWithinAtoms[jj]=1;
          inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][0]=trj_crd->xcoor[jj];
          inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][1]=trj_crd->ycoor[jj];
          inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][2]=trj_crd->zcoor[jj];
          inp_surf->ppfWithinAtomsCoord[iNumOfWithinAtoms][3]=inp_surf->pdAtomsRadii[jj];
        }
      }
    }
  }

  inp_surf->iNumOfAtoms=iNumOfWithinAtoms+1;
  
  MakeNeighbors(inp_surf);
  
  iRotCont=0;
  while(NorthPoleTest(inp_surf)==0)
  {
    MolRotation(inp_surf);
    iRotCont++;
  }

  for(ii=0; ii<inp_surf->sele2.nselatm; ii++)
  {
    dTmpSurf2=GetSurf(ii, inp_surf);
    if(isnan(dTmpSurf2)==1)
    {
      dTmpSurf2=0.0;
    }
    dSurf2=dSurf2+dTmpSurf2;
  }

  inp_surf->iFrameNum++;
  inp_surf->pfSele1Surf[inp_surf->iFrameNum]=dSurf1;
  inp_surf->pfSele2Surf[inp_surf->iFrameNum]=dSurf2;
  
  return;
}

void Compute_GEPOL_SURFCORR(struct inp_surf *inp_surf, Molecule *molecule, CoorSet *trj_crd, char *outstring) //FILE *OutputFile)
{
  int     ii;
  double  dSurf1=0.0;
  double  dSurf2=0.0;

  // === 1st Selection =================================================
  
  // Set inp_surf->iNumOfAtoms to inp_surf->iNumOfRealAtoms
  // and copy frame coords to inp_surf->ppfAtomCoord
  inp_surf->iNumOfAtoms=inp_surf->iNumOfRealAtoms;
  for(ii=0; ii<inp_surf->iNumOfRealAtoms; ii++)
  {
    inp_surf->ppfAtomCoord[ii][0]=trj_crd->xcoor[ii];
    inp_surf->ppfAtomCoord[ii][1]=trj_crd->ycoor[ii];
    inp_surf->ppfAtomCoord[ii][2]=trj_crd->zcoor[ii];
  }

  // Set all atoms as ghost
  for(ii=0; ii<inp_surf->iMaxSize; ii++)
    inp_surf->piAtomType[ii]=3;
    
  // Set selected atoms as 1 or 6 in accordance with radii
  for(ii=0; ii<inp_surf->sele1.nselatm; ii++)
  {
    if(inp_surf->pdAtomsRadii[inp_surf->sele1.selatm[ii]-1]==0)
      inp_surf->piAtomType[inp_surf->sele1.selatm[ii]-1]=1;
    
    else
      inp_surf->piAtomType[inp_surf->sele1.selatm[ii]-1]=6;
  }

  if(inp_surf->sele1.nselatm<inp_surf->iNumOfRealAtoms)
    inp_surf->iGhostFlag=1;

  if(inp_surf->iSurfType==2)
  {
    if(inp_surf->iGhostFlag==1)
      Shell(inp_surf);
      
    Bulk(inp_surf);
    Clean5(inp_surf);
    GenAtoms(inp_surf);
    Clean5(inp_surf);
  }
  
  if(inp_surf->iSurfType==1 && inp_surf->iSolvRadAddedFalg==0)
  {
    inp_surf->iSolvRadAddedFalg=1;
    AddSolvRad(inp_surf);
  }
    
  GeoCav(inp_surf);
  dSurf1=AreaSum(inp_surf);

  // === 2nd Selection =================================================
  
  // Set inp_surf->iNumOfAtoms to inp_surf->iNumOfRealAtoms
  // and copy frame coords to inp_surf->ppfAtomCoord
  inp_surf->iNumOfAtoms=inp_surf->iNumOfRealAtoms;
  for(ii=0; ii<inp_surf->iNumOfRealAtoms; ii++)
  {
    inp_surf->ppfAtomCoord[ii][0]=trj_crd->xcoor[ii];
    inp_surf->ppfAtomCoord[ii][1]=trj_crd->ycoor[ii];
    inp_surf->ppfAtomCoord[ii][2]=trj_crd->zcoor[ii];
  }

  // Set all atoms as ghost
  for(ii=0; ii<inp_surf->iMaxSize; ii++)
    inp_surf->piAtomType[ii]=3;
    
  // Set selected atoms as 1 or 6 in accordance with radii
  for(ii=0; ii<inp_surf->sele2.nselatm; ii++)
  {
    if(inp_surf->pdAtomsRadii[inp_surf->sele2.selatm[ii]-1]==0)
      inp_surf->piAtomType[inp_surf->sele2.selatm[ii]-1]=1;
    
    else
      inp_surf->piAtomType[inp_surf->sele2.selatm[ii]-1]=6;
  }

  if(inp_surf->sele2.nselatm<inp_surf->iNumOfRealAtoms)
    inp_surf->iGhostFlag=1;

  if(inp_surf->iSurfType==2)
  {
    if(inp_surf->iGhostFlag==1)
      Shell(inp_surf);
      
    Bulk(inp_surf);
    Clean5(inp_surf);
    GenAtoms(inp_surf);
    Clean5(inp_surf);
  }
  
  if(inp_surf->iSurfType==1 && inp_surf->iSolvRadAddedFalg==0)
  {
    inp_surf->iSolvRadAddedFalg=1;
    AddSolvRad(inp_surf);
  }
    
  GeoCav(inp_surf);
  dSurf2=AreaSum(inp_surf);  
  
  // === Post Calculation ==============================================
  
  inp_surf->iFrameNum++;
  inp_surf->pfSele1Surf[inp_surf->iFrameNum]=dSurf1;
  inp_surf->pfSele2Surf[inp_surf->iFrameNum]=dSurf2;
  
  return;
}

int Post_SURFCORR(struct inp_surf *inp_surf, struct sopt *OPT)
{
  int          ii;
  char         cOutputFileName[80];
  
  float        fAvg1=0.0, fAvg2=0.0;
  float        fVariance1=0.0, fVariance2=0.0;
  float        fLnAvg1=0.0, fLnAvg2=0.0;
  float        fLnVariance1=0.0, fLnVariance2=0.0;
  float        fLow1=FLT_MAX, fLow2=FLT_MAX;
  float        fHigh1=FLT_MIN, fHigh2=FLT_MIN;
  float        fStdDev1=0.0, fStdDev2=0.0;
  float        fCorrR2Linear=0.0;
  float        fCorrR2Logarithmic=0.0;
  float        fCorrR2Exponential=0.0;
  float        fCorrR2Power=0.0;
  float        fCovar=0.0;

  time_t       time_Today;
  FILE        *OutputFile;
    
  // === Averages, Lowest and Higher values ============================
  for(ii=0; ii<(inp_surf->iFrameNum+1); ii++)
  {
    if(inp_surf->pfSele1Surf[ii]<fLow1)
    {
      fLow1=inp_surf->pfSele1Surf[ii];
    }
    if(inp_surf->pfSele1Surf[ii]>fHigh1)
    {
      fHigh1=inp_surf->pfSele1Surf[ii];
    }
    if(inp_surf->pfSele2Surf[ii]<fLow2)
    {
      fLow2=inp_surf->pfSele2Surf[ii];
    }
    if(inp_surf->pfSele2Surf[ii]>fHigh2)
    {
      fHigh2=inp_surf->pfSele2Surf[ii];
    }
    
    fAvg1=fAvg1+inp_surf->pfSele1Surf[ii];
    fAvg2=fAvg2+inp_surf->pfSele2Surf[ii];
    
    fLnAvg1=fLnAvg1+logf(inp_surf->pfSele1Surf[ii]);
    fLnAvg2=fLnAvg2+logf(inp_surf->pfSele2Surf[ii]);
  }
  
  fAvg1=fAvg1/(float)(inp_surf->iFrameNum+1);
  fAvg2=fAvg2/(float)(inp_surf->iFrameNum+1);
  fLnAvg1=fLnAvg1/(float)(inp_surf->iFrameNum+1);
  fLnAvg2=fLnAvg2/(float)(inp_surf->iFrameNum+1);
  // ===================================================================
  
  // === Standard Deviations ===========================================
  for(ii=0; ii<inp_surf->iFrameNum+1; ii++)
  {
    fVariance1=fVariance1+pow((inp_surf->pfSele1Surf[ii]-fAvg1), 2);
    fVariance2=fVariance2+pow((inp_surf->pfSele2Surf[ii]-fAvg2), 2);
    
    fLnVariance1=fLnVariance1+pow((logf(inp_surf->pfSele1Surf[ii])-fLnAvg1), 2);
    fLnVariance2=fLnVariance2+pow((logf(inp_surf->pfSele2Surf[ii])-fLnAvg2), 2);
  }
  
  fStdDev1=sqrt(fVariance1/(float)(inp_surf->iFrameNum+1));
  fStdDev2=sqrt(fVariance2/(float)(inp_surf->iFrameNum+1));
  // ===================================================================
  
  // === Coefficients of determination R2 ==============================
  for(ii=0; ii<(inp_surf->iFrameNum+1); ii++)
  {
    fCorrR2Linear=fCorrR2Linear+((inp_surf->pfSele1Surf[ii]-fAvg1)*
                                 (inp_surf->pfSele2Surf[ii]-fAvg2));
                                 
    fCorrR2Logarithmic=fCorrR2Logarithmic+((logf(inp_surf->pfSele1Surf[ii])-fLnAvg1)*
                                           (inp_surf->pfSele2Surf[ii]-fAvg2));
                                           
    fCorrR2Exponential=fCorrR2Exponential+((inp_surf->pfSele1Surf[ii]-fAvg1)*
                                           (logf(inp_surf->pfSele2Surf[ii])-fLnAvg2));
    
    fCorrR2Power=fCorrR2Power+((logf(inp_surf->pfSele1Surf[ii])-fLnAvg1)*
                               (logf(inp_surf->pfSele2Surf[ii])-fLnAvg2));
  }
  
  fCovar=fCorrR2Linear/(float)(inp_surf->iFrameNum+1);
  fCorrR2Linear=pow((fCorrR2Linear/sqrt((fVariance1*fVariance2))), 2);
  fCorrR2Logarithmic=pow((fCorrR2Logarithmic/sqrt((fLnVariance1*fVariance2))), 2);
  fCorrR2Exponential=pow((fCorrR2Exponential/sqrt((fVariance1*fLnVariance2))), 2);
  fCorrR2Power=pow((fCorrR2Power/sqrt((fLnVariance1*fLnVariance2))), 2);
  // ===================================================================
  
  
  time(&time_Today);
  sprintf(cOutputFileName, "%s%s", inp_surf->cTitle, ".surfcorr");
  OutputFile=O_File(cOutputFileName, "w");
  
  fprintf(OutputFile, "=========================================\n");
  fprintf(OutputFile, "*** WORDOM SURFace CORRelation Module ***\n");
  fprintf(OutputFile, "=========================================\n");
  fprintf(OutputFile, "Copyright      : Francesca Fanelli\n");
  fprintf(OutputFile, "                 Angelo Felline\n");
  fprintf(OutputFile, "                 (2009)\n");
  fprintf(OutputFile, "License        : GPL v 3\n");
  fprintf(OutputFile, "\n");
  fprintf(OutputFile, "Date           : %s", asctime(localtime(&time_Today)));
  fprintf(OutputFile, "\n");
  fprintf(OutputFile, "Files\n-----\n");
  fprintf(OutputFile, "Mol Name       : %s\n", OPT->IMOL_FILE);
  fprintf(OutputFile, "Total Atoms    : %d\n", inp_surf->iNumOfRealAtoms);
  fprintf(OutputFile, "Trj Name       : %s\n", OPT->ITRJ_FILE);
  fprintf(OutputFile, "Total Frames   : %d\n", inp_surf->iNumOfFrames);
  fprintf(OutputFile, "Used Frames    : %d\n", inp_surf->iFrameNum+1);
  fprintf(OutputFile, "\n");
  fprintf(OutputFile, "Options\n-------\n");
  if(inp_surf->iAlgoFlag==0)
    fprintf(OutputFile, "Algorithm      : ARVO\n");
  else if(inp_surf->iAlgoFlag==1)
  {
    fprintf(OutputFile, "Algorithm      : GEPOL\n");
    if(inp_surf->iSurfType==0)
      fprintf(OutputFile, "Surface        : van der Waals Surface Area\n");
      
    else if(inp_surf->iSurfType==1)
      fprintf(OutputFile, "Surface        : Accessible Surface Area\n");
      
    else if(inp_surf->iSurfType==2)
      fprintf(OutputFile, "Surface        : Excluded Surface Area\n");
      
    fprintf(OutputFile, "NDIV           : %d\n", inp_surf->iNDIV);
    if(inp_surf->iSurfType==2)
    {
      fprintf(OutputFile, "OFAC           : %f\n", inp_surf->fOFAC);
      fprintf(OutputFile, "RMIN           : %f\n", inp_surf->fRMIN);
    }
  }
  fprintf(OutputFile, "Solvent Radius : %f\n", inp_surf->dSolventRadius);
  fprintf(OutputFile, "Sele1          : %s\n", inp_surf->sele1.selestring);  
  fprintf(OutputFile, "Sele1 Atoms    : %d\n", inp_surf->sele1.nselatm);
  fprintf(OutputFile, "Sele2          : %s\n", inp_surf->sele2.selestring);
  fprintf(OutputFile, "Sele2 Atoms    : %d\n", inp_surf->sele2.nselatm);
  fprintf(OutputFile, "\n");
  fprintf(OutputFile, "Statistics\n----------\n");
  fprintf(OutputFile, "Sele1 Lowest   : %f\n", fLow1);
  fprintf(OutputFile, "Sele2 Lowest   : %f\n", fLow2);
  fprintf(OutputFile, "\n");
  fprintf(OutputFile, "Sele1 Highest  : %f\n", fHigh1);
  fprintf(OutputFile, "Sele2 Highest  : %f\n", fHigh2);
  fprintf(OutputFile, "\n");
  fprintf(OutputFile, "Sele1 Average  : %f\n", fAvg1);
  fprintf(OutputFile, "Sele2 Average  : %f\n", fAvg2);
  fprintf(OutputFile, "\n");
  fprintf(OutputFile, "Sele1 StdDev   : %f\n", fStdDev1);
  fprintf(OutputFile, "Sele2 StdDev   : %f\n", fStdDev2);
  fprintf(OutputFile, "\n");
  fprintf(OutputFile, "Covariance     : %f\n", fCovar);
  fprintf(OutputFile, "\n");
  fprintf(OutputFile, "Correlations (R2)\n-----------------\n");
  fprintf(OutputFile, "Linear         : %f\n", fCorrR2Linear);
  fprintf(OutputFile, "Logarithmic    : %f\n", fCorrR2Logarithmic);
  fprintf(OutputFile, "Exponential    : %f\n", fCorrR2Exponential);
  fprintf(OutputFile, "Power          : %f\n", fCorrR2Power);
  fprintf(OutputFile, "\n");
  
  if(inp_surf->iWriteRawFlag==1)
  {
    fprintf(OutputFile, "Raw Data\n--------\n");
    fprintf(OutputFile, "%7s   %15s   %15s\n", "#   nFr", "Sele1", "Sele2");
    for(ii=0; ii<(inp_surf->iFrameNum+1); ii++)
    {
      fprintf(OutputFile, "%7d   %15.3f   %15.3f\n", ii+1, inp_surf->pfSele1Surf[ii], inp_surf->pfSele2Surf[ii]);
    }
  }
  
  return 0;
}
