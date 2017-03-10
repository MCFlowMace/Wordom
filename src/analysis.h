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
/*! \file analysis.h
 \brief 
 
 
*/
#ifndef ANALYSIS
#define ANALYSIS



// ------------------------------------------------------------------
typedef struct _coorsetdata      // coordinate sets data
{
   int          nato;
   int          nframe;
   int          begframe;
   int          endframe;
   int          skip;
   int          i_nUsedFrames;
   TrjList      trj_list;
   short        issingletraj;
   short        issinglestr;
   short        islistoftrajs;
   short        islistofstrs;
   Traj         traj;
   int          pbc_flag;
   Molecule     molecule;
   int          frbefthistrj;   // # of frames in trajs before the one opened
   int          currenttrj;
   short        isfrlist;
   struct intlist frame_list;
   char        *datasetname;
} CoorSetData;
// ------------------------------------------------------------------
// ** Mawaf Input Structure
struct inp_A
{
   int            nBEG;
   int           *begIndex;
   int           *A_type;
   union input  **input;
   short         *bonoffmodulescompute;
   short         *bonoffmodulespost;
   char         **output;
   char        ***printout;
   int            outputstringsize;
   int            nframe;
   short          threaded;
   int            nthreads;
   int            thread_rank;
   FILE          *oA_f;
   sopt          *opt;
   Molecule      *molecule;
};
struct  threading
{
   int            thread_rank;
   int            thread_size;
   struct inp_A   iA_data;
   CoorSetData  *coorsetdata;
   char         **printout;
};
struct temp_ll
{
  struct temp_ll  *next;
  char  *title;
  int    size;
} ;
// ------------------------------------------------------------------

#include "geometry.h"
#include "moldiff.h"
#include "cluster.h"
#include "pca.h"
#include "qentropy.h"
#include "kga.h"
#include "ssa.h"
#include "enm.h"
#include "psn.h"
#include "corr.h"
#include "surf.h"
#include "ring.h"
#include "com.h"
#include "samplemod.h"
#include "volumes.h"
#include "tilt.h"
#include "twist.h"
#include "hole.h"
#include "flux.h"

//===================================================================
//------------- Input Structures for Analysis Functions -------------
//===================================================================
// ------------------------------------------------------------------
// ------------------------------------------------------------------
// ** Kinetic Analysis Pre-processing Structure
struct inp_pproc
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
  Selection         sele;
  float        **reference;
  int        missed_frames;
  float                 T;
  double       **cov;
};
// ------------------------------------------------------------------
// ------------------------------------------------------------------
// ** Correlation Matrix Module data Structure
struct inp_corrmat
{
  char           title[64];
  Selection       sele;
  short          someflag;
  short          extout;
  float         *tempdist;      // temporary buffer for distance calculations
  float         *distance;      // distance of each atom from its average
  float        **reference;
  float        **thisframe;
  float        **matnum;        // matrix of numerators
  float        **matden;        // matrix of denominators
};
// ------------------------------------------------------------------
// ------------------------------------------------------------------
// ------------------------------------------------------------------
// ** Input Files Union
union input
{
    struct inp_angle    inp_angle;
    struct inp_Q        inp_Q;
    struct inp_contacts inp_contacts;
    struct inp_dihe     inp_dihe;
    struct inp_P        inp_P;
    struct inp_dist     inp_dist;
    struct inp_drms     inp_drms;
    struct inp_pca      inp_pca;
    struct inp_pro      inp_pro;
    struct inp_rms      inp_rms;
    struct inp_rms2     inp_rms2;
    struct inp_rmsf     inp_rmsf;
    struct inp_msdf     inp_msdf;
    struct inp_rgyr     inp_rgyr;
    struct inp_Cluster  inp_cluster;
    struct inp_CAssign  inp_cassign;
    struct inp_hb       inp_hb;
    struct inp_pproc    inp_pproc;
    struct inp_ss       inp_ss;
    struct inp_corrmat  inp_corrmat;
    struct inp_enm      inp_enm;
    struct inp_psn      inp_psn;
    struct inp_corr     inp_corr;
    struct inp_surf     inp_surf;
    struct inp_sample   inp_sample;
    struct inp_within   inp_within;
    struct inp_com      inp_com;
    struct inp_tilt     inp_tilt;
    struct inp_twist    inp_twist;
    struct inp_hole     inp_hole;
    struct inp_flux     inp_flux;
};
// ------------------------------------------------------------------
//==============================================================================
//   ANALYSIS FUNCTION
//==============================================================================
void superCalc ( struct sopt   *OPT  );
// ---------------------------------------------------------------------
void iA_Calc ( struct sopt 	*    ) ;
// ---------------------------------------------------------------------
void * iA_Calc_thread ( void 	*    ) ;
// ---------------------------------------------------------------------
void Calc2 ( struct sopt * );
//==============================================================================
int GetThisCoorSet( CoorSetData *coorsetdata, int number, CoorSet *crdset ) ;
// ---------------------------------------------------------------------
void run_iA ( int ii, int intex, CoorSet *trj_crd, struct inp_A *iA_data, char *printout, Molecule *molecule );
// ---------------------------------------------------------------------
void * Run_iA ( void * );
// ------------------------------------------------------------------
char * Read_iA( struct sopt *, struct inp_A *, int * , Molecule *, CoorSetData *coorsetdata );
// ------------------------------------------------------------------
void Read_ia ( struct sopt  *OPT, struct inp_A *iA_data, FILE *outfile, Molecule *molecule, CoorSetData *coorsetdata );
// ------------------------------------------------------------------
int whichModule( char *typestring );
// ------------------------------------------------------------------
int whichEModule( char *typestring );
// ------------------------------------------------------------------
int Get_nBEG ( char **, int ) ;
//----------------------------------------------------------
int Get_nSeles ( char **input, int intex ) ;
//----------------------------------------------------------
void spotBeg ( char **input, int nlines, int *begIndex );
//----------------------------------------------------------
void Get_Atype ( char **, int, int * , int ) ;
void Get_Etype ( char **, int, int * , int ) ;
#endif
