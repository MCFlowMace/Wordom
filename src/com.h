/******************************************************************************
 *
 *  COM Module for WORDOM
 *  Calculate and print the center of mass of a selection along the trajectory
 *
 *  Copyright (C) 2013 Simone Conti
 *
 *  COM Module for WORDOM is free software: you can redistribute it and/or 
 *  modify it under the terms of the GNU General Public License as published 
 *  by the  *  Free Software Foundation, either version 3 of the License, or 
 *  (at your option) any later version.
 *
 *  COM Module for WORDOM is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/


#ifndef _MOD_COM_
#define _MOD_COM_

struct inp_com {
	int		pbc;
	float		*masses;
	Selection	sele;
};

int Read_iCom ( char **input, int inp_index, struct inp_com *inp_com , char *printout, Molecule *molecule, int pbcflag );
int Compute_Com ( struct inp_com *inp_com, struct sopt *OPT, Molecule *molecule, CoorSet *trj_crd, char *outstring );

// Calculate the position of the center of mass of the selection. 
// It is a wrapper function that select between the other four functions. 
// If *masses==NULL it selects the Calc_?_Nomass functions for weighted averages.
// The weights are read from the beta-factor column in the reference pdb.
// If pbc==1 it selects the Calc_Pbc_? functions for periodic boundary conditions.
// Only orthorhombic unit cell is supported; box dimensions are read from CoorSet->Pbc
void Calc_Com(float c[3], CoorSet *trj, Selection sele, float *masses, int pbc);

// Calculate the center of mass using 1 as weights and not using PBC.
void Calc_Com_NopbcNomass(float c[3], CoorSet *trj, Selection sele);

// Calculate the COM with weight but not PBC
void Calc_Com_NopbcMass(float c[3], CoorSet *trj, Selection sele, float *masses);

// Calculate the COM with PBC and no weights
void Calc_Com_PbcNomass(float c[3], CoorSet *trj, Selection sele);

// Calculate the COM with PBC and weights
void Calc_Com_PbcMass(float c[3], CoorSet *trj, Selection sele, float *masses);

#endif

