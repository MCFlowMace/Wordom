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

/*! \file datahandler.h
 \brief Data manipulation functions
 
 Headers for all data manipulation functions. 
*/

#ifndef DATAHANDLER
#define DATAHANDLER
// ------------------------------------------------------------------
//! Concatenates trajectories
/*!
 Merges (cancatenates) different trajectories., listed in a text file
 Al input data are contained in a sopt structure (see wordom.h)
*/
void Trj_Merge ( struct sopt * ) ;
// ------------------------------------------------------------------
void Trj_sum   ( struct sopt * ) ;
// ------------------------------------------------------------------
void Trj_append ( struct sopt * ) ;
// ------------------------------------------------------------------
void Trj_combine ( struct sopt * ) ;
// ------------------------------------------------------------------
void TrjExtract ( struct sopt * ) ;
// ------------------------------------------------------------------
void MolExtract ( struct sopt * ) ;
// ------------------------------------------------------------------
void MultiMolExtract ( struct sopt * ) ;
// ------------------------------------------------------------------
void Mol_append ( struct sopt * ) ;
// ------------------------------------------------------------------
void AvgExtract ( struct sopt * ) ;
// ------------------------------------------------------------------
void TrjHeader_Print ( struct sopt * ) ;
// ------------------------------------------------------------------
void MolInfo_Print( char *filename );
// ------------------------------------------------------------------
void Head_Mod        ( struct sopt * ) ;
// ------------------------------------------------------------------
void xyzExtract ( struct sopt * ) ;
// ------------------------------------------------------------------
void MkMono     ( struct sopt * ) ;
// ------------------------------------------------------------------
void MolConv ( struct sopt *OPT ) ;
// ------------------------------------------------------------------
void CheckSele ( struct sopt   *OPT ) ;
// ------------------------------------------------------------------
void Dcd2Xtc ( struct sopt * ) ;
// ------------------------------------------------------------------
void Xtc2Dcd ( struct sopt * ) ;
// ------------------------------------------------------------------
void mkResList ( struct sopt   *OPT );
//==============================================================================
//   ANALYSIS FUNCTION
//==============================================================================
void iA_Calc ( struct sopt 	*    ) ;
//==============================================================================
void ia_Calc ( char **argv    ) ;
//==============================================================================
#endif
