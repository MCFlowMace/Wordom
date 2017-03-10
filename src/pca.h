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
/*! \file pca.h
 \brief Principal Component Analysisi (PCA) and projectiosn module
 
 Module for PCA and projections along eigenvectos coming from PCA(/ENM)
*/
#ifndef PCA
#define PCA

// =====================================================================
// ** Principal Component Analysis Structure
struct inp_pca
{
  float           *xavgcoor;
  float           *yavgcoor;
  float           *zavgcoor;
  float          **covmatrix;
  float         ***covmatrices;
  char             title[64];
  int              nprint;
  int              skip;
  int              nsteps;
  int              verbose;
  Selection        sele;
  float          **reference;
  int              missed_frames;
  float            T;
  double         **cov;
  float           *xcoor;
  float           *ycoor;
  float           *zcoor;
  float          **coor;
  float           *cords;
  float           *dx;
  float           *dy;
  float           *dz;
};
//----------------------------------------------------------------------
int Read_ipca (char **input, int inp_index, struct inp_pca * , char * , Molecule * , int );
int Compute_PCA (struct inp_pca *, struct sopt *, CoorSet *, char *, int );
int Post_PCA (struct inp_pca *, struct sopt *, int , Molecule * );
// =====================================================================

// =====================================================================
//! Principal Component Analysis Projection Structure
struct inp_pro
{
  int            nmode;
  float         *modevec;
  int           *modes;
  float        **modesvec;
  int           *sel1sel2;
  float         *xrefcoor;
  float         *yrefcoor;
  float         *zrefcoor;
  float         *xtmpcoor;
  float         *ytmpcoor;
  float         *ztmpcoor;
  char           title[64];
  char           dcdout[64];
  int            dcdflag;
  FILE          *dcdoutfile;
  int            rangeflag;
  float          range;
  int            spacflag;
  float          spacing;
  float          minpro;
  float          maxpro;
  int            projected_frames;
  struct trjh    dcdouthdr;
  CoorSet          coor;
  Selection       sele1, sele2;
};
//----------------------------------------------------------------------
int Read_ipro (char **input, int inp_index, struct inp_pro * , char * , Molecule * );
int Compute_pro (struct inp_pro *, struct sopt *, CoorSet *, char *  );
int Post_pro ( struct inp_pro *inp_pro, struct sopt *OPT, int nframe, Molecule *molecule );
// ---------------------------------------------------------------------

#endif
