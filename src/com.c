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


#include "wordom.h"
#include "tools.h"
#include "analysis.h"
#include "com.h"


int Read_iCom ( char **input, int inp_index, struct inp_com *inp_com , char *printout, Molecule *molecule, int pbcflag ) {

    int ii;
    int gotit;
    char *buffer;
    char title[64]="";

    inp_com->pbc    = 0;
    inp_com->masses = NULL;

    buffer=input[inp_index];

    while (strncmp (buffer, "END", 3)) {

        gotit = 0;

        if (!strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#') {
            gotit = 1;
        }

        else if (!strncmp(buffer, "--TITLE", 7)) {
            sscanf(buffer, "--TITLE %s", title);
            gotit = 1;
        }

        else if (!strncmp(buffer, "--SELE", 6) ) {
            sscanf(buffer, "--SELE %[^\n]%*c ", inp_com->sele.selestring);
            GetSele ( inp_com->sele.selestring, &inp_com->sele, molecule);
            gotit = 1;
        }

        else if (!strncmp(buffer, "--MASS", 6)) {
            inp_com->masses = (float *)malloc(molecule->nato * sizeof(float) );
            for( ii=0; ii<molecule->nato; ii++ ) {
                inp_com->masses[ii] = molecule->rawmol.bFac[ii];
            }
            gotit = 1;
        }

        else if (!strncmp(buffer, "--PBC", 5)) {
            if (pbcflag != 1) {
                fprintf( stderr, "ERROR: --PBC used, but no PBC found in mol/trj file nor command line\n" );
                exit(0);
            }
            inp_com->pbc = 1;
            gotit = 1;
        }

        if (gotit==0) {
            fprintf( stderr, "Could not understand option: %s\n", buffer);
            exit(EXIT_FAILURE);
        }

        inp_index++;
        buffer=input[inp_index];
    }


    if (inp_com->sele.nselatm > 1 && inp_com->pbc==0) {
        fprintf( stderr, "Warning: Your selection contains more than 1 atom but you are not using periodic boundary conditions. It may create artifacts if the selection cross the border of the box. To use periodic boundary conditions, just add in your input the --PBC option\n");
    } else if (inp_com->sele.nselatm == 0) {
        fprintf( stderr, "Check selection: %s - 0 atoms selected\n", inp_com->sele.selestring);
        exit(EXIT_FAILURE);
    }

    if (title[0]=='\0') {
        sprintf(printout, "%s", "     COM_X      COM_Y       COM_Z  ");
    } else {
        sprintf(printout, " %8.8s_X  %8.8s_Y  %8.8s_Z ", title, title, title);
    }
	return 36;
}



int Compute_Com ( struct inp_com *inp_com, struct sopt *OPT, Molecule *molecule, CoorSet *trj_crd, char *outstring ) {

    float c[3];

    Calc_Com(c, trj_crd, inp_com->sele, inp_com->masses, inp_com->pbc);
    sprintf(outstring, " %10.5f  %10.5f  %10.5f ", c[0], c[1], c[2]);

    return 36;
}


void Calc_Com(float c[3], CoorSet *trj, Selection sele, float *masses, int pbc) {

	if (sele.nselatm>1 ) {
		if (pbc==1) {
			if (masses!=NULL) {
				Calc_Com_PbcMass(c,trj,sele,masses);
			} else {
				Calc_Com_PbcNomass(c,trj,sele);
				return;
			}
		} else {
			if (masses!=NULL) {
				Calc_Com_NopbcMass(c,trj,sele,masses);
				return;
			} else {
				Calc_Com_NopbcNomass(c,trj,sele);
				return;
			}
		}
	} else {	
		c[0] = trj->xcoor[sele.selatm[0]-1];
		c[1] = trj->ycoor[sele.selatm[0]-1];
		c[2] = trj->zcoor[sele.selatm[0]-1];
		return;
	}
}


void Calc_Com_NopbcNomass(float c[3], CoorSet *trj, Selection sele) {

	int i;

	c[0]=c[1]=c[2]=0;
	
	for (i=0; i<sele.nselatm; i++) {
		c[0] += trj->xcoor[sele.selatm[i]-1];
		c[1] += trj->ycoor[sele.selatm[i]-1];
		c[2] += trj->zcoor[sele.selatm[i]-1];
	}

	c[0] /= sele.nselatm;
	c[1] /= sele.nselatm;
	c[2] /= sele.nselatm;

	return;	
}

void Calc_Com_NopbcMass(float c[3], CoorSet *trj, Selection sele, float *masses) {

	int i;
	float totmass=0;

	c[0]=c[1]=c[2]=0;

	for (i=0; i<sele.nselatm; i++) {
		c[0] += trj->xcoor[sele.selatm[i]-1] * masses[sele.selatm[i]-1];
		c[1] += trj->ycoor[sele.selatm[i]-1] * masses[sele.selatm[i]-1];
		c[2] += trj->zcoor[sele.selatm[i]-1] * masses[sele.selatm[i]-1];
		totmass += masses[sele.selatm[i]-1];
	}

	if (totmass<1) {
		fprintf(stderr, "Warning: your total mass is less than 1. Is that right?\n");
	}

	c[0] /= totmass;
	c[1] /= totmass;
	c[2] /= totmass;
	
	return;
}

void Calc_Com_PbcNomass(float c[3], CoorSet *trj, Selection sele) {

	/* Use the method described in 
		Calculating Center of Mass in an Unbounded 2D Environment
		Journal of Graphics, GPU, and Game Tools
		Linge Bai & David Breen
		http://dx.doi.org/10.1080/2151237X.2008.10129266
	*/

	int i;
	float pbc_x, pbc_y, pbc_z;
	float theta, xi_x=0, zeta_x=0, xi_y=0, zeta_y=0, xi_z=0, zeta_z=0;
	float twopi=0.5/M_PI; // M_PI is pi in math.h, this is 1/(2*pi)

	// Check if the cell is orthorhombic or not
	//if (trj->pbc->angle1!=90 || trj->pbc->angle2!=90 || trj->pbc->angle3!=90) {
	if (trj->pbc->angle1!=0 || trj->pbc->angle2!=0 || trj->pbc->angle3!=0) {
		fprintf(stderr, "Error: only orthorhombic cells are supported.\n");
		fprintf(stderr,"Angles: %f %f %f\n", trj->pbc->angle1,  trj->pbc->angle2,  trj->pbc->angle3); 
	} //Quit or not??

	pbc_x = trj->pbc->a_size*twopi;
	pbc_y = trj->pbc->b_size*twopi;
	pbc_z = trj->pbc->c_size*twopi;

	for (i=0; i<sele.nselatm; i++) {

		theta   = trj->xcoor[sele.selatm[i]-1] / pbc_x;
		xi_x   += pbc_x * cos(theta);
		zeta_x += pbc_x * sin(theta);

		theta   = trj->ycoor[sele.selatm[i]-1] / pbc_y;
		xi_y   += pbc_y * cos(theta);
		zeta_y += pbc_y * sin(theta);

		theta   = trj->zcoor[sele.selatm[i]-1] / pbc_z;
		xi_z   += pbc_z * cos(theta);
		zeta_z += pbc_z * sin(theta);

	}

	float n=1.0/sele.nselatm;

	theta = atan2(-zeta_x/n,-xi_x/n)+M_PI;
	c[0]=pbc_x*theta;
	theta = atan2(-zeta_y/n,-xi_y/n)+M_PI;
	c[1]=pbc_y*theta;
	theta = atan2(-zeta_z/n,-xi_z/n)+M_PI;
	c[2]=pbc_z*theta;

	return;

}



void Calc_Com_PbcMass(float c[3], CoorSet *trj, Selection sele, float *masses) {

	/* Use the method described in 
		Calculating Center of Mass in an Unbounded 2D Environment
		Journal of Graphics, GPU, and Game Tools
		Linge Bai & David Breen
		http://dx.doi.org/10.1080/2151237X.2008.10129266
	*/

	int i;
	float pbc_x, pbc_y, pbc_z;
	float theta, xi_x=0, zeta_x=0, xi_y=0, zeta_y=0, xi_z=0, zeta_z=0;
	float twopi=0.5/M_PI; // M_PI is pi in math.h, this is 1/(2*pi)
	float totmass=0;

	// Check if the cell is orthorhombic or not
	//if (trj->pbc->angle1!=90 || trj->pbc->angle2!=90 || trj->pbc->angle3!=90) {
	if (trj->pbc->angle1!=0 || trj->pbc->angle2!=0 || trj->pbc->angle3!=0) {
		fprintf(stderr, "Error: only orthorhombic cells are supported.\n");
		fprintf(stderr,"Angles: %f %f %f\n", trj->pbc->angle1,  trj->pbc->angle2,  trj->pbc->angle3); 
	} //Quit or not??

	pbc_x = trj->pbc->a_size*twopi;
	pbc_y = trj->pbc->b_size*twopi;
	pbc_z = trj->pbc->c_size*twopi;

	for (i=0; i<sele.nselatm; i++) {

		theta   = trj->xcoor[sele.selatm[i]-1] / pbc_x;
		xi_x   += pbc_x * cos(theta) * masses[sele.selatm[i]-1];
		zeta_x += pbc_x * sin(theta) * masses[sele.selatm[i]-1];

		theta   = trj->ycoor[sele.selatm[i]-1] / pbc_y;
		xi_y   += pbc_y * cos(theta) * masses[sele.selatm[i]-1];
		zeta_y += pbc_y * sin(theta) * masses[sele.selatm[i]-1];

		theta   = trj->zcoor[sele.selatm[i]-1] / pbc_z;
		xi_z   += pbc_z * cos(theta) * masses[sele.selatm[i]-1];
		zeta_z += pbc_z * sin(theta) * masses[sele.selatm[i]-1];
		
		totmass += masses[sele.selatm[i]-1];
	}

	float n=1.0/totmass;

	theta = atan2(-zeta_x/n,-xi_x/n)+M_PI;
	c[0]=pbc_x*theta;
	theta = atan2(-zeta_y/n,-xi_y/n)+M_PI;
	c[1]=pbc_y*theta;
	theta = atan2(-zeta_z/n,-xi_z/n)+M_PI;
	c[2]=pbc_z*theta;

	return;

}
