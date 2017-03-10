// ------------------------------------------------------------------
// Copyright (C) 2009  University of Modena
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
/*! \file ssa.h
 \brief Secondary Structure Assignment module
 
 Headers for all secondary structure assignment procedures
*/
#ifndef SSA
#define SSA

// ------------------------------------------------------------------
// ------------------------------------------------------------------
// ** Secondary Structure Assignment
struct inp_ss
{
  char           title[64];
  Selection       sele;
  int            npairs;
  int            ntotres;
  int           *olist;
  int           *hlist;
  int           *clist;
  int           *nlist;
  int           *calist;
  float         *xcoor1, *xcoor2, *xcoor3, *xcoor4, *xcoor5;
  float         *ycoor1, *ycoor2, *ycoor3, *ycoor4, *ycoor5;
  float         *zcoor1, *zcoor2, *zcoor3, *zcoor4, *zcoor5;
  float         *distON;
  float         *distCH;
  float         *distOH;
  float         *distCN;
  float        **intE;
  short        **hbonds;
  float          cutoffs[10];
  int            weights[10];
  
  char          *secstruct;
  char          *ssnotes;
  int            nLines;
  
  int            dssp;
  int            cont;
  
  char         **secstructS;
  char         **ssnotesS;
};
// ------------------------------------------------------------------

// ------------------------------------------------------------------
// ------------------------------------------------------------------
int Read_iSSA ( char **input, int inp_index, struct inp_ss *inp_ss, char *printout, Molecule *molecule, int nframe );
//------------------------------------------------------------------------------
int Compute_SSA ( struct inp_ss *inp_ss, struct sopt *OPT, Molecule *molecule, CoorSet *trj_crd, char *output );
//------------------------------------------------------------------------------
int computessDistances( struct inp_ss *, CoorSet *dcd_crd );
//------------------------------------------------------------------------------
int mkintEmatrix( struct inp_ss *inp_ss );
//------------------------------------------------------------------------------
int sStruct( struct inp_ss *inp_ss, char *secstruct );
//------------------------------------------------------------------------------


#endif
