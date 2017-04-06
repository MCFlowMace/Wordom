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

/*! \file enm.h
 \brief Elastic Network Models module
 
 Headers for the Elastic Network Models (ENM) analysis module
*/
#ifndef ENM
#define ENM

typedef struct load_dist 
{
 float   drdx;    // difference beween x coor divided by distance
 float   drdy;    // difference beween y coor divided by distance
 float   drdz;    // difference beween z coor divided by distance
 float   dist;    // distance between particle 1 and 2
 int     part1;   // particle1
 int     part2;   // particle2
 float   cost;    // Spring force constant
 int     block1;  // particle 1 's block
 int     block2;  // particle 2 's block
 int     Id;      // Identification number
} _dr;


//RTB METHOD BLOCKS STRUCTURE

typedef struct block
{
  int	       blockid;    /*Block ID number*/
  int	       natom;      /*Number of atoms for each block*/
  int       *atom;	     /*Atoms ID sequential numbers*/
  int       *sindex;	   /*Atoms index in selatm array (for coor collection)*/
  float      x;
  float      y;
  float      z;  
  float     *atom_mass;  /*Atom mass*/
  float      block_mass; /*Block mass*/
  int        resn;
  Selection  sele;
} _blk;

// ------------------------------------------------------------------
//ELASTIC NETWORK MODELS MODULE Structure
struct inp_enm
{
////ENM INPUT FILE
  char         title[128];      
  Selection    sele1;   //Selection for the construction of the starting HEssian
  Selection    vsasele; //Selection for the subsytem
  Selection    env;     //Selection for the environment
  Selection    sele2;   //Selection for a second molecule to be compared
  Selection    selean;  //Selection for analyses module
  Selection    rtblevan;//Selection that specify the level of analysis for the RTB method - refers to selmol molecule
  Selection    rtblevan_org;//Selection that specify the level of analysis for the RTB method - refers to original molecule
  Molecule    *molecule;
  Molecule    *molecule2;
  Molecule    *selmol;
  int          sele1_flag;
  int          sele2_flag;
  int	         vsa_flag;
  char         rtb_method[64];
  int          rtb_flag;
  int          rtblevan_flag;
  float        cutoff;
  char         type[64];
  char         mol2[1024];
  char         vecfile2[1024];
  int	         mol2_flag;
  int	         vecfile2_flag;
  char        *beta;
  int	         beta_flag;
  char	       temp[64];
  char        *perturb;
  int	         perturb_flag;
  char        *corr;
  int	         corr_flag;
  char        *defen;
  int	         defen_flag;  
  int	         perturb_norm;
  char        *nmodes;  
  char         matrix[1024];
  int          matrix_chk;
  int          network_print;
  char        *proj;
  int          proj_flag;
  int          dcd_flag;
  int          mass_flag;
  int          essential;
  int          essential_flag;
  char         enm_level[64];
  int          enm_level_flag;

/////ENM BASIC    
  int	         nint;//Pairwise interactions number
  int          print; 
  int          current_frame;        
  float      **He;//Hessian of the system selected by --SELE           
 _dr	        *dr;
 _eigen       *eigen;
  FILE        *evectrj;
  CoorSet     *evecoor;           
  struct trjh  evechdr;
  float       *distances;
  float      **coord;
  float      **virt_coor;
  int          nselpart;
    
//ENM PERTURBATION
  float      **Hp;//Perturbed Hessian of the system selected by --SELE
//  float    **Hs;//He-Hp
  float       *scalar;
  float     ***pertmat;
  float      **delta_omega;
  float      **delta_omega_plot;
  int          pertmat_flag;

//ENM CORRELATION & FLUCTUATION  
  float      **cov;
  float      **cov2;
  float      **avg;
  float       *fluc;

//ENM COMPARISON BETWEEN NORMAL MODES AND CONFORMATIONAL TRANSITION VECTOR  
  float      **reference;
  float      **moving_sys;
  float      **moving_selean;
  float       *bfact;
  float       *deviation;
  float       *overlap; 
  float      **vecmat2;
  int          matvec_size;
  float      **mat_overlap;	
  float       *norm;

//ENM PROJECTIONS
  FILE        *projdcdfile;
  char         projdcd[64];
  struct trjh  projdcdhdr;
  CoorSet     *coor; 
  
//ENM DEFORMATION ENERGY
//  float       *deformation;
  int          distmat_flag;

//ENM VSA METHOD 
  int        **sys_list;       
  int        **subsys_list;
  int        **env_list;
  float      **Hss;//Subsystem Hessian selected with --SUBSELE
  float      **Hee;//Environment Hessian
  float      **Hes;//Interface Hessian between environment/subsystem;
  float      **invHee;//inverse of the Environment Hessian;
 _eigen       *subeigen;
  _dr         *drvsa;
  int          nintvsa;

//ENM RTB METHOD 
  int          nblock;
 _blk         *blk;
  float     ***RT;
 _eigen       *rtb2full_eigen;
  int         *atm2blk;
  
};

//------------------------------------------------------------------------------

/*ENM input file reading function*/
int Read_iEnm ( char **input, int inp_index, struct inp_enm *inp_enm, char *printout, Molecule *molecule );

/*Basic ENM computation module*/
int Compute_Enm ( struct inp_enm *inp_enm, Molecule *molecule, CoorSet *trj_crd, char *outprint, int intex );

/*ENM related functions*/
int Get_Dr ( struct inp_enm *inp_enm, _dr *dr, CoorSet *trj_crd, Molecule *molecule ); /*Function to compute 2nd derivative of the Harmonic Potential function*/
int Get_DrVSA ( struct inp_enm *inp_enm, _dr *dr, _dr *drvsa );
int Build_Hessian ( struct inp_enm *inp_enm, float **hessian , _dr *dr, int );/*Hessian Building function*/
int Build_Perthessian ( struct inp_enm *inp_enm, float **hessian , _dr *dr, int, int *pert_list, int pert_list_size );

/*ENM analysis functions*/
int Enm_Correlation ( struct inp_enm *inp_enm, _eigen *eig,  Molecule *molecule );
int Enm_Fluc        ( struct inp_enm *inp_enm,  Molecule *molecule, _eigen *eig );
int Enm_Perturb     ( struct inp_enm *inp_enm, Molecule *molecule, _eigen *eig, _dr *dr, CoorSet *trj_crd );
int Enm_Compare     ( struct inp_enm *inp_enm, Molecule *molecule, CoorSet *trj_crd, _eigen *eig );
int Enm_GetVsaMat   ( struct inp_enm *inp_enm, float **input_matrix , float **output_matrix);
int Enm_DiagVsaMat  ( struct inp_enm *inp_enm, _eigen * );
int Enm_Rtb         ( struct inp_enm *inp_enm, Molecule *molecule );
int Enm_Rtb2full    ( struct inp_enm *inp_enm, Molecule *molecule );
int Enm_Getblock    ( struct inp_enm *inp_enm, Molecule *molecule );
int Enm_Defen       ( struct inp_enm *inp_enm,  Molecule *molecule, _eigen *eig, CoorSet *trj_crd );
int Enm_Proj        ( struct inp_enm *inp_enm,  _eigen *eig, Molecule *molecule, CoorSet *trj_crd );
int Enm_Variance    ( struct inp_enm *inp_enm, _eigen *eig);
/*To be moved into tools */
/*Matrix inverse function*/
int GeomCentre      (struct inp_enm *inp_enm, CoorSet *crd, Molecule *molecule);
extern void sgetrf_ ( int *, int *, float *, int *, int *, int * );         /*LU factorization*/
extern void sgetri_ ( int *, float *, int *, int *, float *, int *, int * );/*Inverse of a matrix*/
void InvMatrix      ( int, float **, float **);
/* Graham-Schmidt orthonormalization function*/
void OrtoNorm       (int, int, float **);

#endif
