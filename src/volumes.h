// ------------------------------------------------------------------
// Copyright (C) 2012  University of Modena and Reggio Emilia and
//                     University of Zurich
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

/*! \file volumes.h
 \brief volumes/superposition module
 
 Headers for the Volumes/Superposition module.
*/
#ifndef VOLUMES
#define VOLUMES
// ------------------------------------------------------------------
// ** volumes inp Structure
struct inp_vol
{
  char            supermol[128];
  char            mollist[128];
  char            boxfile[128];
  int             boxflag;
  int             centerboxflag;
  int             cornerboxflag;
  int             stepflag;
  float           step;
  int             nt;
  int             mthread;
  float           box_xcenter;
  float           box_ycenter;
  float           box_zcenter;
  float           box_xcorner;
  float           box_ycorner;
  float           box_zcorner;
  float           box_xsize;
  float           box_ysize;
  float           box_zsize;
};
typedef struct _moldata
{
  int       natoms;
  float    *radii;
  float   **coords;
} moldata;
typedef struct _boxdim
{
  float       xorg;
  float       yorg;
  float       zorg;
  float       xsize;
  float       ysize;
  float       zsize;
  int         xbins;
  int         ybins;
  int         zbins;
  int         totalbins;
  float       step;
} BoxDim;

int CalcVolumes( char **input, int inp_index, FILE *output );

int linecount( FILE *inpfile );
#endif
