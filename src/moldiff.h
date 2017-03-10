// ------------------------------------------------------------------
// Copyright (C) 2003  University of Zurich
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
/*! \file moldiff.h
 \brief All Molecule similarity (et similar) module
 
 Headers for all molecule similarity (et similar) procedures
*/
#ifndef MOLDIFF
#define MOLDIFF

//=== RMSD PARAMETERS =================================================
//! Rmsd Structure
struct inp_rms
{
  char           title[64];
  Selection       sele;
  Selection       fitsele;
  short          super;
  short          trjwrite;
  short          progressive;
  short          extout;
  short          fit;
  char           trjout[64];
  trjtype        outrj_type;
  FILE          *outtrj;
  Traj          *outtraj;
  CoorSet        crd;
  float        **reference;
  float        **moving;
  float        **fitref;
  float        **fitmov;
  float         *xcmpcoor;
  float         *ycmpcoor;
  float         *zcmpcoor;
};
//----------------------------------------------------------
int Read_iRmsd(char **input, int inp_index, struct inp_rms * , char * , Molecule *, CoorSetData *, sopt * );
int Compute_Rmsd (struct inp_rms *, struct sopt *, CoorSet *, char *  );
//======================================================================
//=== RMSD PARAMETERS =================================================
//! Rmsd Structure
struct inp_rms2
{
  char           title[64];
  Selection       sele;
  Selection       fitsele;
  short          super;
  short          trjwrite;
  short          progressive;
  short          extout;
  short          fit;
  char           trjout[64];
  trjtype        outrj_type;
  FILE          *outtrj;
  Traj          *outtraj;
  CoorSet        crd;
  float        **reference;
  float        **moving;
  float        **fitref;
  float        **fitmov;
  float          centerRef[3];
  float          centerMov[3];
};
//----------------------------------------------------------
int Read_iRmsd2(char **input, int inp_index, struct inp_rms2 * , char * , Molecule *, CoorSetData *, sopt * );
int Compute_Rmsd2 (struct inp_rms2 *, struct sopt *, CoorSet *, char *  );
//======================================================================

//=== DRMS PARAMETERS =================================================
//! Drms Structure
struct inp_drms
{
  int           *selatm;         // Indexes of selected atoms 
  int            nselatm;         // Number of selected atoms
  Selection      sele;
  short          progressive;
  short          cmapd;
  float          cmapd_cutoff;   // Cutoff for the contact map
  float         *distancemat;
};
//----------------------------------------------------------
int Read_idrms ( char **input, int inp_index, struct inp_drms *, char *, Molecule *,struct sopt *);
int Compute_drms ( struct inp_drms *, struct sopt *, CoorSet *, char *  );
int Calc_dist_mtx( struct sopt *, int  , int * , CoorSet * , float * );
//======================================================================

//=== RMSF PARAMETERS =================================================
//! Rmsf Structure
struct inp_rmsf
{
  char           title[64];
  Selection       sele;
  Selection       fitsele;
  short          super;
  short          fit;
  trjtype        outrj_type;
  FILE          *outtrj;
  Traj         *outtraj;
  CoorSet    crd;
  float        **reference;
  float        **moving;
  float        **fitref;
  float        **fitmov;
  float         *xcmpcoor;
  float         *ycmpcoor;
  float         *zcmpcoor;
  float         *outrmsf;
  int            framecounter;
};
// ------------------------------------------------------------------
int Read_iRmsf ( char **input, int inp_index, struct inp_rmsf *inp_rmsf, char *, Molecule *molecule );
int Compute_Rmsf ( struct inp_rmsf *inp_rmsf, struct sopt *OPT, CoorSet *trj_crd, char * );
int Post_Rmsf ( struct inp_rmsf *inp_rmsf, struct sopt *OPT, int nframe, Molecule *molecule );
float fret(float , float );
// ------------------------------------------------------------------

//=== MSDF PARAMETERS =================================================
//! Msdf Structure
struct inp_msdf
{
  char           title[64];
  Selection      sele;
  CoorSet        crd;
  float       ***distmat;
  double       **d_mat;
  double       **sqmat;
  double       **meanmat;
  double       **fluctmat;
  int            framecounter;
  short          pass;
  
  
  // added by me
  char         **pcLabels;                                              // Labels (segname:restype+resnum)
  int            root;
  int            verbose;
  int            iNumOfFrames;
  int            iNumOfRes;
  int            iNumOfAtoms;
  int            iNumOfSelRes;
  int            iNumOfSelAtoms;
  int            iMultiAtomFlag;
  int           *piSelResLen;                                           // Number of atoms per selected residue 
  int           *piProgNum;                                             // Progressive res numbers, used if iMultiAtomFlag == 1
  int            iPBCFlag;
  float        **ppfVirtAtomCoord;                                      // Frame Virtual Atom Coordinates, used if iMultiAtomFlag == 1
  float        **ppfVirtRefAtomCoord;                                   // Reference Virtual Atom Coordinates, used if iMultiAtomFlag == 1
};
// ------------------------------------------------------------------
int Read_iMsdf ( char **input, int inp_index, struct inp_msdf *inp_msdf, char *, Molecule *molecule, int nframe );
int Compute_Msdf ( struct inp_msdf *inp_msdf, struct sopt *OPT, CoorSet *trj_crd, char * );
int Post_Msdf ( struct inp_msdf *inp_msdf, struct sopt *OPT, int nframe, Molecule *molecule );
float fret(float , float );
// ------------------------------------------------------------------


#endif
