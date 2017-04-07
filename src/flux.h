/******************************************************************************
 *
 *  FLUX Module for WORDOM
 *  Calculate and print the flux angle between two selections
 *
 *  Copyright (C) 2015 Nicolas Martin, Simone Conti
 *  Copyright (C) 2015 Unviersité de Strasbourg 
 *
 *  FLUX Module for WORDOM is free software: you can redistribute it and/or 
 *  modify it under the terms of the GNU General Public License as published 
 *  by the Free Software Foundation, either version 3 of the License, or 
 *  (at your option) any later version.
 *
 *  FLUX Module for WORDOM is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/


#ifndef _MOD_FLUX_
#define _MOD_FLUX_


/* Structure of flux, filled by the function Read_iFlux */

struct inp_flux {
  Selection  sele;
  Selection  seleCenter; // Optional selection to define the center of the channel
  double     Rcut;
  double     delta;
  double     upper;
  double     lower;
  int       *status;
  FILE      *transout;
};

int Read_iFlux ( char **input, int inp_index, struct inp_flux *inp_flux , char *title, Molecule *molecule );
int Compute_Flux ( struct inp_flux *inp_flux, CoorSet *trj_crd, char *outstring, int frame );

#endif /* _MOD_FLUX_ */
