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


#include "wordom.h"
#include "tools.h"
#include "analysis.h"
#include "hole.h"


int 
Read_iHole(char **input, int inp_index, struct inp_hole *inp_hole, char *printout, Molecule *molecule)
{
    char *buffer;
    int  gotit, i;

    /* Set defaults */

    /* Title */
    inp_hole->title[0]='\0';

    /* Initial point to origin */
    inp_hole->Pzero[0] = 0;
    inp_hole->Pzero[1] = 0;
    inp_hole->Pzero[2] = 0;

    /* Channel vector to z axis */
    inp_hole->V[0] = 0;
    inp_hole->V[1] = 0;
    inp_hole->V[2] = 1;

    /* Vector of vdw radii */
    inp_hole->vdw = NULL;

    /* Seed for random number generator */
    inp_hole->seed=0;

    /* Maximum displacement in Monte Carlo search */
    inp_hole->Dmax = 0.1;

    /* "kBT" in Monte Carlo search */
    inp_hole->K = 0.1;

    /* Maximum number of iterations in Monte Carlo search */
    inp_hole->Nmax = 1000;

    /* Step scale in Monte Carlo search */
    inp_hole->scale = 0.9;

    /* Displacement along channel vector */
    inp_hole->d = 0.1;

    /* Maximum radius which stop the run */
    inp_hole->Rmax = 10;

    /* Cutoff for neighbour list */
    inp_hole->cutoff = 5;

    /* Read input */

    buffer = input[inp_index];
    while(strncmp(buffer, "END", 3)) {

        gotit = 0;

        if(!strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#') {
            gotit = 1;
        } 

        else if (!strncmp(buffer, "--TITLE", 7)) {
            sscanf(buffer, "--TITLE %s", inp_hole->title);
            gotit = 1;
        }

        else if (!strncmp(buffer, "--SELE", 6) ) {
            sscanf(buffer, "--SELE %[^:]", inp_hole->sele.selestring);
            GetSele(inp_hole->sele.selestring, &inp_hole->sele, molecule);
            if( inp_hole->sele.nselatm < 1 ) {
                fprintf(stderr, "\nERROR! Selection %s has less than 1 atom (have %d)\n\n", 
                            inp_hole->sele.selestring, inp_hole->sele.nselatm);
                exit(0);
            }
            if (inp_hole->sele.selestring[strlen(inp_hole->sele.selestring)-1]=='\n') {
                inp_hole->sele.selestring[strlen(inp_hole->sele.selestring)-1]='\0';
            }
            gotit = 1;
            fprintf(stderr, "Number of selected atoms: %d\n" ,  inp_hole->sele.nselatm);
        }

        else if (!strncmp(buffer, "--RMAX", 6)) {
            sscanf(buffer, "--RMAX %lf", &(inp_hole->Rmax));
            if (inp_hole->Rmax<=0) {
                fprintf(stderr, "\nERROR! RMAX must be greater than zero! (is %lf)\n\n", inp_hole->Rmax);
                exit(0);
            }
            gotit = 1;
        }

        else if (!strncmp(buffer, "--DISP", 6)) {
            sscanf(buffer, "--DISP %lf", &(inp_hole->d));
            if (inp_hole->d<=0) {
                fprintf(stderr, "\nERROR! DISP must be greater than zero! (is %lf)\n\n", inp_hole->d);
                exit(0);
            }
            gotit = 1;
        }

        else if (!strncmp(buffer, "--SEED", 6)) {
            sscanf(buffer, "--SEED %d", &(inp_hole->seed));
            if (inp_hole->seed<=0) {
                fprintf(stderr, "\nERROR! SEED must be greater than zero! (is %d)\n\n", inp_hole->seed);
                exit(0);
            }
            gotit = 1;
        }

        else if (!strncmp(buffer, "--NMAX", 6)) {
            sscanf(buffer, "--NMAX %d", &(inp_hole->Nmax));
            if (inp_hole->Nmax<=0) {
                fprintf(stderr, "\nERROR! NMAX must be greater than zero! (is %d)\n\n", inp_hole->Nmax);
                exit(0);
            }
            gotit = 1;
        }

        else if (!strncmp(buffer, "--KBT", 5)) {
            sscanf(buffer, "--KBT %lf", &(inp_hole->K));
            if (inp_hole->K<=0) {
                fprintf(stderr, "\nERROR! KBT must be greater than zero! (is %lf)\n\n", inp_hole->K);
                exit(0);
            }
            gotit = 1;
        }

        else if (!strncmp(buffer, "--DMAX", 6)) {
            sscanf(buffer, "--DMAX %lf", &(inp_hole->Dmax));
            if (inp_hole->Dmax<=0) {
                fprintf(stderr, "\nERROR! DMAX must be greater than zero! (is %lf)\n\n", inp_hole->Dmax);
                exit(0);
            }
            gotit = 1;
        }

        else if (!strncmp(buffer, "--SCALE", 7)) {
            sscanf(buffer, "--SCALE %lf", &(inp_hole->scale));
            if (inp_hole->scale<=0) {
                fprintf(stderr, "\nERROR! SCALE must be greater than zero! (is %lf)\n\n", inp_hole->scale);
                exit(0);
            }
            gotit = 1;
        }

        else if (!strncmp(buffer, "--CUTOFF", 8)) {
            sscanf(buffer, "--CUTOFF %lf", &(inp_hole->cutoff));
            if (inp_hole->cutoff<=0) {
                fprintf(stderr, "\nERROR! CUTOFF must be greater than zero! (is %lf)\n\n", inp_hole->cutoff);
                exit(0);
            }
            gotit = 1;
        }

        else if (!strncmp(buffer, "--POINT", 7)) {
            sscanf(buffer, "--POINT (%lf:%lf:%lf)", &(inp_hole->Pzero[0]), &(inp_hole->Pzero[1]), &(inp_hole->Pzero[2]));
            gotit = 1;
        }

        else if (!strncmp(buffer, "--VECT", 6)) {
            sscanf(buffer, "--VECT (%lf:%lf:%lf)", &(inp_hole->V[0]), &(inp_hole->V[1]), &(inp_hole->V[2]));
            double norm = sqrt( inp_hole->V[0]*inp_hole->V[0] + inp_hole->V[1]*inp_hole->V[1] + inp_hole->V[2]*inp_hole->V[2] );
            inp_hole->V[0] /= norm;
            inp_hole->V[1] /= norm;
            inp_hole->V[2] /= norm;
            gotit = 1;
        }

        
        if( gotit==0 ) {
            fprintf( stderr, "Could not understand option: %s\n", buffer);
            exit(5);
        }
        inp_index++;
        buffer=input[inp_index];
    }


    /* Setup atomic radii -- TODO Improve reading a rad file */
    /* Use Bondi, A. (1964) J. Phys. Chem. 68:441-451 */
    inp_hole->vdw = malloc(inp_hole->sele.nselatm*sizeof(double));
    for (i=0; i<inp_hole->sele.nselatm; i++) {
        switch (molecule->rawmol.atmtype[inp_hole->sele.selatm[i]-1][0]) {
            case 'C': inp_hole->vdw[i] =  1.70; break;
            case 'O': inp_hole->vdw[i] =  1.52; break;
            case 'S': inp_hole->vdw[i] =  1.80; break;
            case 'N': inp_hole->vdw[i] =  1.55; break;
            case 'H': inp_hole->vdw[i] =  1.20; break;
            case 'P': inp_hole->vdw[i] =  1.80; break;
            default: 
                inp_hole->vdw[i]=2; 
                fprintf(stderr, "Warning! Unknown element %s in assign vdw radius. Using 2ang.\n\n", molecule->rawmol.atmtype[i]);
                break;
        }
    }

    /* Initialize random number generator */
    srand(inp_hole->seed);

    /* Check title */
    if (inp_hole->title[0]=='\0') {
        strncpy(inp_hole->title, "pore", 64);
    }
    strncpy(inp_hole->prof, inp_hole->title, 64);
    strncpy(inp_hole->traj, inp_hole->title, 64);
    strncpy(inp_hole->visu, inp_hole->title, 64);
    strncat(inp_hole->prof, "-prof.dat", 64);
    strncat(inp_hole->traj, "-traj.vtf", 64);
    strncat(inp_hole->visu, "-visu.tcl", 64);

    /* Allocate some memory... */
    inp_hole->dataX = malloc(sizeof(double*));
    inp_hole->dataY = malloc(sizeof(double*));
    inp_hole->dataZ = malloc(sizeof(double*));
    inp_hole->dataR = malloc(sizeof(double*));
    inp_hole->dataN = malloc(sizeof(int));
    inp_hole->dataS = malloc(sizeof(int));
    inp_hole->nframe = 1;


    if (inp_hole->seed==0) inp_hole->seed=time(NULL);

    /* Print parsed input */
    fprintf(stderr, "\n");
    fprintf(stderr, "Parsed HOLE input\n");
    fprintf(stderr, "  Selected atoms:     %10d (%s)\n", inp_hole->sele.nselatm, inp_hole->sele.selestring);
    fprintf(stderr, "  Initial point:      %10.4lf %10.4lf %10.4lf\n", inp_hole->Pzero[0], inp_hole->Pzero[1], inp_hole->Pzero[2]);
    fprintf(stderr, "  Channel vector:     %10.4lf %10.4lf %10.4lf\n", inp_hole->V[0], inp_hole->V[1], inp_hole->V[2]);
    fprintf(stderr, "  Step displacement:  %10.4lf\n", inp_hole->d);
    fprintf(stderr, "  Maximum radius:     %10.4lf\n", inp_hole->Rmax);
    fprintf(stderr, "  Neighbour cutoff:   %10.4lf\n", inp_hole->cutoff);
    fprintf(stderr, "  Monte Carlo search setup:\n");
    fprintf(stderr, "    Seed for RNG:         %10d\n", inp_hole->seed);
    fprintf(stderr, "    Maximum displament:   %10.4lf\n", inp_hole->Dmax);
    fprintf(stderr, "    Temperature:          %10.4lf\n", inp_hole->K);
    fprintf(stderr, "    Number of iterations: %10.4d\n", inp_hole->Nmax);
    fprintf(stderr, "    Annealing scale:      %10.4lf\n", inp_hole->scale);
    fprintf(stderr, "  Output files:\n");
    fprintf(stderr, "    Pore profile:      %s\n", inp_hole->prof);
    fprintf(stderr, "    Pore trajectory:   %s\n", inp_hole->traj);
    fprintf(stderr, "    VMD visualization: %s\n", inp_hole->visu);
    fprintf(stderr, "\n");
     

    return 0;
}


int Compute_Hole(struct inp_hole *inp_hole, CoorSet *trj_crd, char *outstring)
{

    int i, j;
    int natm = inp_hole->sele.nselatm;
    int Nmax = inp_hole->Nmax;
    int dispid=0, dispidmax=0, dispidmin=0;
    int nnbond;
    double *xcoor = malloc(natm*sizeof(double));
    double *ycoor = malloc(natm*sizeof(double));
    double *zcoor = malloc(natm*sizeof(double));
    double *vdw = malloc(natm*sizeof(double));
    double Rmax=inp_hole->Rmax, Rnew=Rmax, Rold=-99999, Rtmp;
    double Poldx, Poldy, Poldz, Pnewx, Pnewy, Pnewz;
    double Yrandx, Yrandy, Yrandz;
    double Vx, Vy, Vz;
    double disp=inp_hole->d;
    double Dmax=inp_hole->Dmax;
    double scale=inp_hole->scale;
    double K=inp_hole->K;
    double dot, norm;
    double cutoff;
    double *dataRf=NULL, *dataRb=NULL;
    double *dataXf=NULL, *dataXb=NULL;
    double *dataYf=NULL, *dataYb=NULL;
    double *dataZf=NULL, *dataZb=NULL;

    Vx = inp_hole->V[0];
    Vy = inp_hole->V[1];
    Vz = inp_hole->V[2];

    Poldx = inp_hole->Pzero[0];
    Poldy = inp_hole->Pzero[1];
    Poldz = inp_hole->Pzero[2];

    while (1) {

        /* Create non bonded list */
        cutoff=Rnew+inp_hole->cutoff;
        nnbond = 0;
        for (i=0; i<natm; i++) {
            Rtmp  = (Poldx-trj_crd->xcoor[inp_hole->sele.selatm[i]-1])*(Poldx-trj_crd->xcoor[inp_hole->sele.selatm[i]-1]);
            Rtmp += (Poldy-trj_crd->ycoor[inp_hole->sele.selatm[i]-1])*(Poldy-trj_crd->ycoor[inp_hole->sele.selatm[i]-1]);
            Rtmp += (Poldz-trj_crd->zcoor[inp_hole->sele.selatm[i]-1])*(Poldz-trj_crd->zcoor[inp_hole->sele.selatm[i]-1]);
            if (Rtmp<cutoff*cutoff) {
                xcoor[nnbond] = trj_crd->xcoor[inp_hole->sele.selatm[i]-1];
                ycoor[nnbond] = trj_crd->ycoor[inp_hole->sele.selatm[i]-1];
                zcoor[nnbond] = trj_crd->zcoor[inp_hole->sele.selatm[i]-1];
                vdw[nnbond]   = inp_hole->vdw[i];
                nnbond++;
            }
        }

        /* Do Monte Carlo search on this plane */
        for (j=0; j<Nmax; j++) {

            /* Get random displacement in plane perpendicular to V */
            Yrandx = 2.0*rand()/RAND_MAX-1;
            Yrandy = 2.0*rand()/RAND_MAX-1;
            Yrandz = 2.0*rand()/RAND_MAX-1;
            norm = sqrt(Yrandx*Yrandx + Yrandy*Yrandy + Yrandz*Yrandz);
            Yrandx /= norm;
            Yrandy /= norm;
            Yrandz /= norm;
            dot = Yrandx*Vx + Yrandy*Vy + Yrandz*Vz;
            Yrandx -= dot*Vx;
            Yrandy -= dot*Vy;
            Yrandz -= dot*Vz;

            /* Get new point P */
            dot=1.0*rand()/RAND_MAX;
            Pnewx = Poldx + Yrandx*Dmax*dot;
            Pnewy = Poldy + Yrandy*Dmax*dot;
            Pnewz = Poldz + Yrandz*Dmax*dot;

            /* Calculate R at Pnew */
            Rnew = Rmax+1;
            for (i=0; i<nnbond; i++) {
                Rtmp  = (Pnewx-xcoor[i])*(Pnewx-xcoor[i]);
                Rtmp += (Pnewy-ycoor[i])*(Pnewy-ycoor[i]);
                Rtmp += (Pnewz-zcoor[i])*(Pnewz-zcoor[i]);
                if (Rtmp < (Rnew+vdw[i])*(Rnew+vdw[i])) { 
                    Rnew = sqrt(Rtmp)-vdw[i];
                }
            }

            /* See if to accept the move or not */
            if ( (Rnew>Rold) || (exp((Rnew-Rold)/K)>1.0*rand()/RAND_MAX) ) {
                Rold = Rnew;
                Poldx = Pnewx;
                Poldy = Pnewy;
                Poldz = Pnewz;
            }
    
            K *= scale;
        }

        /* Save position and radius of the hole */
        if (dispid>=0) {
            if (dispid>=dispidmax) {
                dispidmax += 10;
                dataRf = realloc(dataRf, dispidmax*sizeof(double));
                dataXf = realloc(dataXf, dispidmax*sizeof(double));
                dataYf = realloc(dataYf, dispidmax*sizeof(double));
                dataZf = realloc(dataZf, dispidmax*sizeof(double));
            }
            dataRf[dispid] = Rold;
            dataXf[dispid] = Poldx;
            dataYf[dispid] = Poldy;
            dataZf[dispid] = Poldz;
            dispid++;
        } else {
            if (-dispid>=dispidmin) {
                dispidmin += 10;
                dataRb = realloc(dataRb, dispidmin*sizeof(double));
                dataXb = realloc(dataXb, dispidmin*sizeof(double));
                dataYb = realloc(dataYb, dispidmin*sizeof(double));
                dataZb = realloc(dataZb, dispidmin*sizeof(double));
            }
            dataRb[-dispid] = Rold;
            dataXb[-dispid] = Poldx;
            dataYb[-dispid] = Poldy;
            dataZb[-dispid] = Poldz;
            dispid--;
        }

        /* Update for next step? */
        if (Rold>Rmax) {
            if (disp>0) {
                disp *= -1;
                Poldx = inp_hole->Pzero[0]+Vx*disp;
                Poldy = inp_hole->Pzero[1]+Vy*disp;
                Poldz = inp_hole->Pzero[2]+Vz*disp;
                dispidmax = dispid;
                dispid = -1;
            } else {
                dispidmin = -dispid;
                break;
            }
        } else {
            Poldx += Vx*disp;
            Poldy += Vy*disp;
            Poldz += Vz*disp;
        }
        Rold = -999999;
        K = inp_hole->K;
    }

    /* Save profile in order */
    inp_hole->dataN = realloc(inp_hole->dataN, (inp_hole->nframe+1)*sizeof(int));
    inp_hole->dataS = realloc(inp_hole->dataS, (inp_hole->nframe+1)*sizeof(int));
    inp_hole->dataR = realloc(inp_hole->dataR, (inp_hole->nframe+1)*sizeof(double*));
    inp_hole->dataX = realloc(inp_hole->dataX, (inp_hole->nframe+1)*sizeof(double*));
    inp_hole->dataY = realloc(inp_hole->dataY, (inp_hole->nframe+1)*sizeof(double*));
    inp_hole->dataZ = realloc(inp_hole->dataZ, (inp_hole->nframe+1)*sizeof(double*));
    inp_hole->dataN[inp_hole->nframe] = dispidmin+dispidmax-1;
    inp_hole->dataR[inp_hole->nframe] = malloc((dispidmin+dispidmax-1)*sizeof(double));
    inp_hole->dataX[inp_hole->nframe] = malloc((dispidmin+dispidmax-1)*sizeof(double));
    inp_hole->dataY[inp_hole->nframe] = malloc((dispidmin+dispidmax-1)*sizeof(double));
    inp_hole->dataZ[inp_hole->nframe] = malloc((dispidmin+dispidmax-1)*sizeof(double));
    int npt=0;
    for (i=dispidmin-1; i>0; i--) {
        inp_hole->dataR[inp_hole->nframe][npt] = dataRb[i];
        inp_hole->dataX[inp_hole->nframe][npt] = dataXb[i];
        inp_hole->dataY[inp_hole->nframe][npt] = dataYb[i];
        inp_hole->dataZ[inp_hole->nframe][npt] = dataZb[i];
        npt++;
    }
    for (i=0; i<dispidmax; i++) {
        inp_hole->dataR[inp_hole->nframe][npt] = dataRf[i];
        inp_hole->dataX[inp_hole->nframe][npt] = dataXf[i];
        inp_hole->dataY[inp_hole->nframe][npt] = dataYf[i];
        inp_hole->dataZ[inp_hole->nframe][npt] = dataZf[i];
        npt++;
    }
    inp_hole->dataS[inp_hole->nframe] = 1-dispidmin;
    inp_hole->nframe++;

    /* Clean memory and exit */
    free(dataRf);
    free(dataXf);
    free(dataYf);
    free(dataZf);
    free(dataRb);
    free(dataXb);
    free(dataYb);
    free(dataZb);
    free(xcoor);
    free(ycoor);
    free(zcoor);
    free(vdw);

    return 0;
}

int Post_Hole(struct inp_hole *inp_hole, struct sopt *OPT, int nframe, Molecule *molecule )
{

    int i, j, dispidmin, nptmax, diff, diffmax;


    /* Get maximum number of points per frame and minimum displacement */
    nptmax=inp_hole->dataN[1];
    dispidmin=inp_hole->dataS[1];
    for (i=1; i<inp_hole->nframe; i++) {
        dispidmin = inp_hole->dataS[i] < dispidmin ? inp_hole->dataS[i] : dispidmin;
        nptmax    = inp_hole->dataN[i] > nptmax    ? inp_hole->dataN[i] : nptmax;
    }
    diffmax=0;
    for (i=1; i<inp_hole->nframe; i++) {
        diff = inp_hole->dataS[i]-dispidmin;
        diffmax = diff>diffmax ? diff : diffmax;
    }


    /* Print profiles */
    fprintf(stderr, "Printing profiles along channel vector to %s\n", inp_hole->prof);
    FILE *prof=fopen(inp_hole->prof, "w");
    fprintf(prof, "%10s ", "#displ");
    for (i=1; i<inp_hole->nframe; i++) {
        fprintf(prof, "%9df ", i);
    }
    fprintf(prof, "\n");
    for (j=0; j<nptmax+diffmax; j++) {
        fprintf(prof, "%10.2lf ", (j+dispidmin)*inp_hole->d);
        for (i=1; i<inp_hole->nframe; i++) {
            diff = inp_hole->dataS[i]-dispidmin;
            if (j<diff || j-diff>=inp_hole->dataN[i]) {
                fprintf(prof, "%10.2lf ", INFINITY);
            } else {
                fprintf(prof, "%10.2lf ", inp_hole->dataR[i][j-diff]);
            }
        }
        fprintf(prof, "\n");
    }
    fclose(prof);


    /* Print VTF trajectory */
    fprintf(stderr, "Printing hole trajectory to %s\n", inp_hole->traj);
    FILE *traj=fopen(inp_hole->traj, "w");
    fprintf(traj, "atom 0:%d name O\n\n", nptmax+diffmax-1);
    for (i=1; i<inp_hole->nframe; i++) {
        fprintf(traj, "timestep\n");
        for (j=0; j<nptmax+diffmax; j++) {
            diff = inp_hole->dataS[i]-dispidmin;
            if (j<diff || j-diff>=inp_hole->dataN[i]) {
                fprintf(traj, "%10.2lf %10.2lf %10.2lf %10.2lf\n", 0.0, 0.0, 0.0, 0.0);
            } else {
                fprintf(traj, "%10.2lf %10.2lf %10.2lf %10.2lf\n", inp_hole->dataX[i][j-diff], 
                inp_hole->dataY[i][j-diff], inp_hole->dataZ[i][j-diff], inp_hole->dataR[i][j-diff]);
            }
        }
        fprintf(traj, "\n");
    }
    fclose(traj);


    /* Print VMD reading script */
    fprintf(stderr, "Printing VMD script to visualize to %s\n", inp_hole->visu);
    FILE *vmd=fopen(inp_hole->visu, "w");
    fprintf(vmd, " \n\
package require vtftools \n\
proc change_radius_callback { name molid op } { \n\
    global vmd_frame radius_update_sel \n\
    set pattern \"$molid.step$vmd_frame($molid).\" \n\
    set radii {} \n\
    $radius_update_sel update \n\
    foreach pid [ $radius_update_sel list ] { \n\
        lappend radii $::VTFTools::userdata($pattern$pid) \n\
    } \n\
    $radius_update_sel set radius $radii \n\
    $radius_update_sel set beta $radii \n\
} \n\
proc sample_load { filename } { \n\
    global radius_update_sel \n\
    vtf_load $filename \n\
    set molid [ molinfo \"top\" ] \n\
    set radius_update_sel [ uplevel atomselect $molid \"all\" ] \n\
    uplevel trace add variable vmd_frame($molid) write change_radius_callback \n\
    change_radius_callback {} $molid {} \n\
    return $molid \n\
} \n\
set molid [sample_load %s] \n\
mol modstyle 0 $molid vdw \n\
mol modcolor 0 $molid beta \n\
mol new %s \n\
mol addfile %s \n\
animate delete beg 0 end 0 skip 0 1 \n\
mol modstyle 0 1 NewCartoon 0.300000 10.000000 4.100000 0 \n\
mol modcolor 0 1 Chain \n\
    ", inp_hole->traj, OPT->IMOL_FILE, OPT->ITRJ_FILE);
    fclose(vmd);


    /* Clean memory and exit */
    free(inp_hole->vdw); inp_hole->vdw=NULL;
    free(inp_hole->dataN); inp_hole->dataN=NULL;
    free(inp_hole->dataS); inp_hole->dataS=NULL;
    for (i=1; i<inp_hole->nframe; i++) {
        free(inp_hole->dataX[i]); inp_hole->dataX[i]=NULL;
        free(inp_hole->dataY[i]); inp_hole->dataY[i]=NULL;
        free(inp_hole->dataZ[i]); inp_hole->dataZ[i]=NULL;
        free(inp_hole->dataR[i]); inp_hole->dataR[i]=NULL;
    }
    free(inp_hole->dataX); inp_hole->dataX=NULL;
    free(inp_hole->dataY); inp_hole->dataY=NULL;
    free(inp_hole->dataZ); inp_hole->dataZ=NULL;
    free(inp_hole->dataR); inp_hole->dataR=NULL;
    
    return 0;
}
