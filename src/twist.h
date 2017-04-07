/******************************************************************************
 *
 *  TWIST Module for WORDOM
 *  Calculate and print the twist angle between two selections relative to
 *  a delocated axis
 *
 *  Copyright (C) 2014 Nicolas Martin, Simone Conti
 *
 *  TWIST Module for WORDOM is free software: you can redistribute it and/or 
 *  modify it under the terms of the GNU General Public License as published 
 *  by the  *  Free Software Foundation, either version 3 of the License, or 
 *  (at your option) any later version.
 *
 *  TWIST Module for WORDOM is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/


#ifndef _MOD_TWIST_
#define _MOD_TWIST_



/* Structure of twist, filled by the function read_itwist*/

struct inp_twist 
{
  Selection  seleA;
  Selection  seleB;
  Selection  seleAxis;
};

int Read_iTwist ( char **input, int inp_index, struct inp_twist *inp_twist , char *title, Molecule *molecule );
int Compute_Twist ( struct inp_twist *inp_twist, CoorSet *trj_crd, char *outstring );

#endif /* _MOD_TWIST_ */
