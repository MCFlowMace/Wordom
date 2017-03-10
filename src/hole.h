/******************************************************************************
 *
 *  HOLE Module for WORDOM
 *
 *  Copyright (C) 2015 Simone Conti
 *  Copyright (C) 2015 Unviersité de Strasbourg 
 *
 *  HOLE Module for WORDOM is free software: you can redistribute it and/or 
 *  modify it under the terms of the GNU General Public License as published 
 *  by the Free Software Foundation, either version 3 of the License, or 
 *  (at your option) any later version.
 *
 *  HOLE Module for WORDOM is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#ifndef _MOD_HOLE_
#define _MOD_HOLE_

struct inp_hole {
    Selection sele;
    int seed;
    int Nmax;
    int nframe;
    double Pzero[3];
    double V[3];
    double Dmax;
    double K;
    double scale;
    double d;
    double cutoff;
    double Rmax;
    double **dataX;
    double **dataY;
    double **dataZ;
    double **dataR;
    double *vdw;
    int *dataN;
    int *dataS;
    char title[64];
    char prof[64];
    char traj[64];
    char visu[64];
};

int Read_iHole(char **input, int inp_index, struct inp_hole *inp_hole, char *title, Molecule *molecule );
int Compute_Hole(struct inp_hole *inp_hole, CoorSet *trj_crd, char *outstring);
int Post_Hole(struct inp_hole *inp_hole, struct sopt *OPT, int nframe, Molecule *molecule );
    
#endif /* _MOD_HOLE_ */
