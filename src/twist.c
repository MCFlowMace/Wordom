/******************************************************************************
 *
 *  TWIST Module for WORDOM
 *  Calculate and print the twist angle between two selections and a delocalized
 *  axis.
 * 
 *  Copyright (C) 2014 Nicolas Martin, Simone Conti
 *
 *  TWIST Module for WORDOM is free software: you can redistribute it and/or 
 *  modify it under the terms of the GNU General Public License as published 
 *  by the Free Software Foundation, either version 3 of the License, or 
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


#include "wordom.h"
#include "tools.h"
#include "analysis.h"
#include "twist.h"
#include "tilt.h"


/*Read the command line wordom and parse it for the twist selection */

int 
Read_iTwist ( char **input, int inp_index, struct inp_twist *inp_twist , char *printout, Molecule *molecule )
{

    char          buffer[256];
    char          title[64]="Twist";
    int           gotit;

    memset ( buffer, '\0', sizeof(buffer));
    while( strncmp (buffer, "END", 3)) {
        gotit = 0;
        memset ( buffer, '\0', sizeof(buffer));
        sprintf( buffer, "%s", input[inp_index]);

        if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
            gotit = 1;

        else if ( !strncmp(buffer, "--TITLE", 7)) {
            memset( title, '\0', 64);
            sscanf( buffer, "--TITLE %s", title);
            gotit = 1;
        }

        else if ( !strncmp(buffer, "--SELE1",7 ) ) {
            sscanf( buffer, "--SELE1 %[^:]", inp_twist->seleA.selestring);
            GetSele ( inp_twist->seleA.selestring, &inp_twist->seleA, molecule);
            if( inp_twist->seleA.nselatm < 2 ) {
                fprintf( stderr, "Check selection: %s - %d atoms selected\n", inp_twist->seleA.selestring, inp_twist->seleA.nselatm);
                exit(0);
            }
            gotit = 1;
        }

        else if ( !strncmp(buffer, "--SELE2", 7) ) {
            sscanf( buffer, "--SELE2 %[^:]", inp_twist->seleB.selestring);
            GetSele ( inp_twist->seleB.selestring, &inp_twist->seleB, molecule);
            if( inp_twist->seleB.nselatm > 1 )  
                fprintf( stderr, "Warning: GC doesn't work with PBC yet; centers will be computed w/o PBC, distances with PBC\n");
    
            if( inp_twist->seleB.nselatm < 2 ) {
                fprintf( stderr, "Check selection: %s - %d atoms selected\n", inp_twist->seleB.selestring, inp_twist->seleB.nselatm);
                exit(0);
            }
            gotit = 1;
        }

        else if ( !strncmp(buffer, "--SELEAXIS", 10) ) {
            sscanf( buffer, "--SELEAXIS %[^:]", inp_twist->seleAxis.selestring);
            GetSele ( inp_twist->seleAxis.selestring, &inp_twist->seleAxis, molecule);
            if( inp_twist->seleAxis.nselatm < 2 ) {
                fprintf( stderr, "Check selection: %s - %d atoms selected\n", inp_twist->seleAxis.selestring, inp_twist->seleAxis.nselatm);
                exit(0);
            } 
            gotit = 1;
        }
   
        if( gotit==0 ) {
            fprintf( stderr, "Could not understand option: %s\n", buffer);
            exit(5);
        }
        inp_index++;
    }
  
    sprintf( printout, " %10s ", title);
    return 12;
}  



// calculate le twist angle from the two selections specified relative to the axis defined by seleAxis
    
int 
Compute_Twist ( struct inp_twist *inp_twist, CoorSet *trj_crd, char *outstring)
{

    float twist, scalar;
    float v1[3], v2[3];
    float proj1, proj2;
    float **mtxAxis;
    float v1diff[3], v2diff[3], *vAll;
    float com1[3], com2[3], com3[3]; 

    /*printf("# Coord of sele\n");
    for( jj=0; jj<inp_twist->sele->nselatm; jj++ ) {
        printf ("%f\t",trj_crd->xcoor[inp_twist->sele->selatm[jj]-1]);
        printf ("%f\t",trj_crd->ycoor[inp_twist->sele->selatm[jj]-1]);
        printf ("%f\t",trj_crd->zcoor[inp_twist->sele->selatm[jj]-1]);
        printf("\n");
    }
    */

    // 1) find principal axes of rotation of SELEAXIS to use them as reference
    Compute_Inertia(trj_crd, &inp_twist->seleAxis, &mtxAxis);

    vAll = mtxAxis[2]; // eigenvec with the lowest eigenval
    w_norm(vAll);

    // 2) compute COM for SELE1 and SELE2 and SELEAXIS

    Calc_Com_NopbcNomass(com1, trj_crd, inp_twist->seleA); 
    Calc_Com_NopbcNomass(com2, trj_crd, inp_twist->seleB);
    Calc_Com_NopbcNomass(com3, trj_crd, inp_twist->seleAxis);

    // 3) projection of the COM on Z axis 
  
    // v1diff and v2diff are the vectors that goes from COM of seleAxis to COM of Sele1 and sele2 respectively
    v1diff[0] = com1[0] - com3[0];
    v1diff[1] = com1[1] - com3[1];
    v1diff[2] = com1[2] - com3[2];
    w_norm(v1diff);

    v2diff[0] = com2[0] - com3[0];
    v2diff[1] = com2[1] - com3[1];
    v2diff[2] = com2[2] - com3[2];
    w_norm(v2diff);

    // 4) projection of v1diff and v2 diff on plan orthogonal to vAll

    proj1 = dot_prod(v1diff,vAll);
    proj2 = dot_prod(v2diff,vAll);

    v1[0] = v1diff[0] - proj1 * vAll[0];
    v1[1] = v1diff[1] - proj1 * vAll[1];
    v1[2] = v1diff[2] - proj1 * vAll[2];

    v2[0] = v2diff[0] - proj2 * vAll[0];
    v2[1] = v2diff[1] - proj2 * vAll[1];
    v2[2] = v2diff[2] - proj2 * vAll[2];

    w_norm(v1);
    w_norm(v2);

    // 5) angle calculation

    scalar = dot_prod(v1,v2);
    twist = compute_angle(scalar);  
    sprintf( outstring, " %10.5f ", twist);

    return 12;
}

