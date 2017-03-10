/******************************************************************************
 *
 *  FLUX Module for WORDOM
 *  Calculate and print the number of molecule experiencing a transition 
 *  through a proteic channel.
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

#include "wordom.h"
#include "tools.h"
#include "analysis.h"
#include "flux.h"

int 
Read_iFlux (char **input, int inp_index, struct inp_flux *inp_flux, char *printout, Molecule *molecule)
{
    char *buffer;
    char title[12];
    char transout[64];
    int  gotit, i;

    /* Set default values */
    strncpy(title, "flux", 4);
    strncpy(inp_flux->sele.selestring, "/*/*/OH2", 8);
    inp_flux->upper=INFINITY;
    inp_flux->lower=INFINITY;
    inp_flux->delta=5.0;
    strncpy(inp_flux->seleCenter.selestring, "/*/*/CA", 7);
    inp_flux->Rcut=0;

    /* Parsing options */
    buffer=input[inp_index];
    while (strncmp(buffer, "END", 3)!=0) {
        gotit = 0;

        if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#') {
            gotit = 1;
        } 

        else if ( !strncmp(buffer, "--TITLE", 7)) {
            sscanf(buffer, "--TITLE %10s", title);
            gotit = 1;
        }

        else if ( !strncmp(buffer, "--SELE", 6) ) {
            sscanf( buffer, "--SELE %[^:]", inp_flux->sele.selestring);
            gotit = 1;
        }
 
        else if ( !strncmp(buffer, "--REF", 5) ) {
            sscanf( buffer, "--REF %[^:]", inp_flux->seleCenter.selestring);
            gotit = 1;
        }
          
        else if ( !strncmp(buffer, "--UPPER", 7) ) {
            sscanf( buffer, "--UPPER %lf", &(inp_flux->upper));
            gotit = 1;
        }
      
        else if ( !strncmp(buffer, "--LOWER", 7) ) {
            sscanf( buffer, "--LOWER %lf", &(inp_flux->lower));
            gotit = 1;
        }
        else if ( !strncmp(buffer, "--DELTA", 7) ) {
            sscanf( buffer, "--DELTA %lf", &(inp_flux->delta));
            gotit = 1;
        }
        else if ( !strncmp(buffer, "--RCUT", 6) ) {
            sscanf( buffer, "--RCUT %lf", &(inp_flux->Rcut));
            gotit = 1;
        }
        if( gotit==0 ) {
            fprintf( stderr, "Could not understand option: %s\n", buffer);
            exit(5);
        }
        inp_index++;
        buffer=input[inp_index];
    }


    /* Getting selections */
    if (inp_flux->sele.selestring[strlen(inp_flux->sele.selestring)-1]=='\n') {
        inp_flux->sele.selestring[strlen(inp_flux->sele.selestring)-1]='\0';
    }
    if (inp_flux->seleCenter.selestring[strlen(inp_flux->seleCenter.selestring)-1]=='\n') {
        inp_flux->seleCenter.selestring[strlen(inp_flux->seleCenter.selestring)-1]='\0';
    }
    GetSele(inp_flux->seleCenter.selestring, &inp_flux->seleCenter, molecule);
    GetSele(inp_flux->sele.selestring, &inp_flux->sele, molecule);

    /* Print/check parsed options */
    fprintf(stderr, "FLUX computation:\n");
    fprintf(stderr, "  Upper boundary     : %10.2f\n" , inp_flux->upper);
    if (inp_flux->upper == INFINITY) {
        fprintf(stderr, "\n    ERROR! You MUST set something with --UPPER!\n\n");
        exit(-10);
    }
    fprintf(stderr, "  Lower boundary     : %10.2f\n" , inp_flux->lower);
    if (inp_flux->lower == -INFINITY) {
        fprintf(stderr, "\n    ERROR! You MUST set something with --LOWER!\n\n");
        exit(-10);
    }
    fprintf(stderr, "  Delta distance     : %10.2f\n" , inp_flux->delta);
    fprintf(stderr, "  Tracked molecules  : %10d (selection: %s)\n", inp_flux->sele.nselatm, inp_flux->sele.selestring);
    if (inp_flux->sele.nselatm<=0) {
        fprintf(stderr, "\n ERROR! Your selection is empty! What should I track?\n\n");
        exit(-10);
    }
    if (inp_flux->Rcut<=0) {
        fprintf(stderr, "  Not using any cutoff cylinder.\n");
    } else {
        fprintf(stderr, "  Cutoff distance    : %10.2f\n" , inp_flux->Rcut);
        fprintf(stderr, "  Cutoff selection   : %10d (selection: %s)\n", inp_flux->seleCenter.nselatm, inp_flux->seleCenter.selestring);
        if (inp_flux->seleCenter.nselatm<=0) {
            fprintf(stderr, "\n WARNING! Your selection is empty!\n\n");
        }
    }
    strncpy(transout, title, 64);
    strncat(transout, "-trans.dat", 64);
    inp_flux->transout = fopen(transout, "w");
    fprintf(stderr, "  Writing list of transitions to %s\n", transout);
    fprintf(stderr, "\n");

    /* Write header in transition file */
    fprintf(inp_flux->transout, "#%9s %10s\n", "Frame", "AtomIndex");

    /* Allocate and initialize status array */
    inp_flux->status = malloc(inp_flux->sele.nselatm * sizeof(int));
    for (i=0; i<inp_flux->sele.nselatm; i++) {
        inp_flux->status[i] = 3;
    }
    sprintf( printout, " %10s ", title );
    return 12;
}

int Compute_Flux(struct inp_flux *inp_flux, CoorSet *trj, char *outstring, int frame)
{
    int   i, statusCur, nTransitions=0;
    double upper = inp_flux->upper;
    double lower = inp_flux->lower;
    double delta = inp_flux->delta;
    double Rcut  = inp_flux->Rcut;
    double dx, dy;
    float *xcoor=trj->xcoor;
    float *ycoor=trj->ycoor;
    float *zcoor=trj->zcoor;
    float comCenter[3];
    Selection sele=inp_flux->sele;

    /* If the user specified a selection for the definition of the channel,
         caulcate its center of geometry */
    if (Rcut>0) {
        Calc_Com_NopbcNomass(comCenter, trj, inp_flux->seleCenter);
    }

    /* Loop over all tracked atoms */
    for (i=0; i<sele.nselatm; i++) {

        statusCur=-1;

        /* If the channel is speficied, filter all atoms out of the cylinder of radius Rcut */
        if (Rcut>0) {
            dx = xcoor[sele.selatm[i]-1]-comCenter[0];
            dy = ycoor[sele.selatm[i]-1]-comCenter[1];
            if ( dx*dx+dy*dy > Rcut*Rcut ) {
                statusCur = 3;
            }
        }

        /* If not filetered out, let see its status */
        if (statusCur!=3) {
            if (zcoor[sele.selatm[i]-1] >= upper+delta) {
                statusCur = 3;
            } else if (zcoor[sele.selatm[i]-1] >= upper) {
                statusCur = 1;
            } else if (zcoor[sele.selatm[i]-1] <= lower-delta) {
                statusCur = 3;
            } else if (zcoor[sele.selatm[i]-1] <= lower) {
                statusCur = 0;
            } else {
                statusCur = inp_flux->status[i];
            }
        }

        /* Check if there is any transition */
        if (statusCur!=3) {
            if (inp_flux->status[i]!=3) {
                if (statusCur != inp_flux->status[i]) {
                    /* Transiton spotted! Print frame and atom index to transout file 
                        and increment number of transitions */
                    fprintf(inp_flux->transout, "%10d %10d\n", frame, sele.selatm[i]);
                    nTransitions += 1;
                }
            }
        }

        /* Update status */
        inp_flux->status[i] = statusCur;
    }

	sprintf(outstring, " %10d ", nTransitions);
    return 12;
}

