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
/*! \file geometry.h
 \brief Simple geometrical analyses module
 
 Module for simple geometrical analyses, such as distances, contacts,
 angles, dihedrals, Root Mean Square Deviation (RMSD), Distance Root
 Mean Square Deviation (DRMS) etc.
*/
#ifndef GEOM
#define GEOM


//=== DISTANCE =========================================================
//! Distance Structure
struct inp_dist
{
  _distvcts  vcts;
  int        nLines;
  int      **ind0, **ind1;
  short      mass;
  float     *masses;
  
  Selection *seleA;
  Selection *seleB;
};
//--- distance functions -----------------------------------------------
int Read_iDist ( char **input, int inp_index, struct inp_dist *inp_dist , struct temp_ll *temp_ll, Molecule *molecule );
int Compute_Dist ( struct inp_dist *, struct sopt *, Molecule *, CoorSet *, char * );   // options struct, molecule struct, coordinates struct, output file
//======================================================================

//=== CONTACTS =========================================================
//! Contact Structure
struct inp_Q
{

  _distvcts     vcts;
  int           nLines;
  
  int         **seleA, **seleB;
  int          *counter;
  int          *contacts;
  int          *nContacts;
  float        *cutoff;
  short         mass;
  float        *masses;
};
//--- contacts functions -----------------------------------------------
int Read_iQ (char **input, int inp_index, struct inp_Q * , char * , Molecule * );
int Compute_Q ( struct inp_Q *, struct sopt *, Molecule *, CoorSet *, char *outstring );   // options struct, molecule struct, coordinates struct, output file
//======================================================================

//=== CONTACTS =========================================================
//! Contact Structure
struct inp_contacts
{

  _distvcts     vcts;
  int           nLines;
  
  int          *counter;
  int          *contacts;
  int          *nContacts;
  float        *cutoff;
  short         mass;
  float        *masses;
  Selection *seleA;
  Selection *seleB;
};
//--- contacts functions -----------------------------------------------
//int Read_iQ (char **input, int inp_index, struct inp_Q * , char * , Molecule * );
//int Compute_Q ( struct inp_Q *, struct sopt *, Molecule *, CoorSet *, char *outstring );   // options struct, molecule struct, coordinates struct, output file
//int Read_iContacts (char **input, int inp_index, struct inp_contacts * , char * , Molecule * );
int Read_iContacts ( char **input, int inp_index, struct inp_contacts *inp_contacts , struct temp_ll *temp_ll, Molecule *molecule );
int Compute_Contacts ( struct inp_contacts *, struct sopt *, Molecule *, CoorSet *, char *outstring );   // options struct, molecule struct, coordinates struct, output file
//======================================================================

//=== ANGLE ============================================================
//! Angle structure
/*!
  Structure for angle-between-atoms computation.
*/
struct inp_angle
{
   int  atm1, atm2, atm3;              /*<! \brief id of atoms defining the angle*/
   int  someflag;
   int  extout;
};
//--- angle functions --------------------------------------------------
int Read_iAngle ( char **input, int inp_index, struct inp_angle *inp_angle, char *printout, Molecule *molecule );
int Compute_Angle ( struct inp_angle *inp_angle, struct sopt *OPT, CoorSet *trj_crd, char *outstring );
void Post_Angle ( struct inp_angle *inp_angle, struct sopt *OPT, int nframe, Molecule *molecule );
//======================================================================

//=== DIHEDRALS ========================================================
//! Dihedral Structure
struct inp_dihe
{
  int            nLines;
  int            atm1, atm2, atm3, atm4;
};
//--- dihedrals functions ----------------------------------------------
int Read_iDihe (char **input, int inp_index, struct inp_dihe * , char * , Molecule * );
int Compute_Dihe ( struct inp_dihe *, struct sopt *, Molecule *, CoorSet *, char *outstring );   // options struct, molecule struct, coordinates struct, output file
//======================================================================

//=== RADIUS OF GYRATION ===============================================
//! Rgyr Structure
struct inp_rgyr
{
  char           title[64];
  Selection       sele;
  CoorSet    crd;
  float        **moving;
  short          mass;
  float         *masses;
};
//--- radius of gyration functions -------------------------------------
int Read_iRgyr(char **input, int inp_index, struct inp_rgyr * , char * , Molecule *, CoorSetData * );
int Compute_Rgyr (struct inp_rgyr *, struct sopt *, CoorSet *, char * );
//----------------------------------------------------------

//=== ORDER PARAMETERS =================================================
//! Order Parameter Structure
struct inp_P
{
  struct line3
   {
    int tail_sID;       // tail atom segment IDs 
    int tail_rID;       // tail atom residue IDs 
    char tail_aT[5];    // tail atom types

    int head_sID;       // head atom segment IDs 
    int head_rID;       // head atom residue IDs 
    char head_aT[5];    // head atom types
   }           **inps;

  int           nLines;
  int          *atm1;
  int          *atm2;
};
//--- o.p. functions ----------------------------------------------
int Read_iP (char **input, int inp_index, struct inp_P * , char * , Molecule * );
int Compute_P ( struct inp_P *, struct sopt *, Molecule *, CoorSet *, char * );   // options struct, molecule struct, coordinates struct, output file
//======================================================================

//=== HYDROGEN BONDS ===================================================
//! Hydrogen Bonds Structure
struct inp_hb
{
  Selection       sele1, sele2, sele3;
  int           *ind0, *ind1, *ind2;
  int            ntitles;
  int            nLines;
  int           *titlen;
  
  float         *xcoor1, *xcoor2, *xcoor3, *ycoor1, *ycoor2, *ycoor3, *zcoor1, *zcoor2, *zcoor3;
  float         *dist;
  
  float          hbdist;
  float          hbangle;
  float          hbrad;
  
//const float    hbdist =   2.4;        // distance donor-acceptor
//const float    hbangl = 130.0;        // angle donor-hydrogen-acceptor
  
};
// --- h.b. functions --------------------------------------------------
int Read_iHB(char **input, int inp_index, struct inp_hb * , char * , Molecule * );
int Compute_hb (struct inp_hb *, struct sopt *, CoorSet *, char *  );
//======================================================================

//=== WITHIN MODULE ====================================================
//! Within Structure
struct inp_within
{
  int          iVerboseFlag;                                            // Used to check --VERBOSE option
  int          iLevelFlag;                                              // Used to check --LEVEL option, 0 = ATM, 1 = RES
  
  char         cTitle[80];                                              // Output column name and verbose file name
  char         cOutPutFileName[80];                                     // Verbose file name
  
  FILE        *FOutputFile;                                             // Verbose file handler
  
  Selection    sele;                                                    // Selection structure
  
  Molecule    *TmpMolecule;                                             // Used to store frame coords
};
// --- Within functions ------------------------------------------------
int Init_Within(char **input, int input_index, struct inp_within *inp_within, Molecule *molecule, char *outstring);
int Compute_Within(struct inp_within *inp_within, Molecule *molecule, CoorSet *trj_crd, char *outstring, int intex);
//======================================================================

#endif
