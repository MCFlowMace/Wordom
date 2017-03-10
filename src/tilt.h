/******************************************************************************
 *
 *  TILT Module for WORDOM
 *  Calculate and print the tilt angle between two selections
 *
 *  Copyright (C) 2014, 2015 Nicolas Martin, Simone Conti, Florian Blanc
 *
 *  TILT Module for WORDOM is free software: you can redistribute it and/or 
 *  modify it under the terms of the GNU General Public License as published 
 *  by the  *  Free Software Foundation, either version 3 of the License, or 
 *  (at your option) any later version.
 *
 *  TILT Module for WORDOM is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#ifndef _MOD_TILT_
#define _MOD_TILT_

/* Structure of tilt, filled by the function read_iTilt */

struct inp_tilt {
  Selection  seleRef;
  Selection  sele;
  int decompose;
  int noorient;
};

int Read_iTilt ( char **input, int inp_index, struct inp_tilt *inp_tilt , char *title, Molecule *molecule );
int Compute_Tilt ( struct inp_tilt *inp_tilt, CoorSet *trj_crd, char *outstring );

float Compute_Inertia(CoorSet *trj_crd, Selection *sele, float ***mtx);
float compute_angle(float dot);
//float compute_angle_no_orientation(float dot);

#endif /* _MOD_TILT_ */
