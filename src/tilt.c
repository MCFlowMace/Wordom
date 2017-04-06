/******************************************************************************
 *
 *  TILT Module for WORDOM
 *  Calculate and print the tilt angle between two selections
 *
 *  Copyright (C) 2014, 2015 Nicolas Martin, Simone Conti, Florian Blanc
 *
 *  TILT Module for WORDOM is free software: you can redistribute it and/or 
 *  modify it under the terms of the GNU General Public License as published 
 *  by the Free Software Foundation, either version 3 of the License, or 
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


#include "wordom.h"
#include "tools.h"
#include "analysis.h"
#include "tilt.h"


int 
Read_iTilt ( char **input, int inp_index, struct inp_tilt *inp_tilt , char *printout, Molecule *molecule )
{
    //char buffer[256];
    char *buffer;
    char title[64]="";
    int  gotit;

    inp_tilt->decompose=0; /* If 1 request polar decomposition */
    inp_tilt->noorient=0;  /* If 1 request the NOORIENT option */

    //memset ( buffer, '\0', sizeof(buffer));
    buffer = input[inp_index];
    while( strncmp (buffer, "END", 3)) {

        gotit = 0;
        //memset ( buffer, '\0', sizeof(buffer));
        //sprintf( buffer, "%s", input[inp_index]);

        if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#') {
            gotit = 1;
        } 

        else if ( !strncmp(buffer, "--TITLE", 7)) {
            memset( title, '\0', 64);
            sscanf( buffer, "--TITLE %s", title);
            gotit = 1;
        }

        else if ( !strncmp(buffer, "--SELEREF", 9) ) {
            sscanf( buffer, "--SELEREF %[^:]", inp_tilt->seleRef.selestring);
            GetSele ( inp_tilt->seleRef.selestring, &inp_tilt->seleRef, molecule);
            if( inp_tilt->seleRef.nselatm < 2 ) {
                fprintf( stderr, "\nERROR! Selection %s must have more than 2 atoms (have %d)\n\n", inp_tilt->seleRef.selestring, inp_tilt->seleRef.nselatm);
                exit(0);
            }
            gotit = 1;
        }

        else if ( !strncmp(buffer, "--SELE", 6) ) {
            sscanf( buffer, "--SELE %[^:]", inp_tilt->sele.selestring);
            GetSele ( inp_tilt->sele.selestring, &inp_tilt->sele, molecule);
            if( inp_tilt->sele.nselatm > 1 ) {
                fprintf( stderr, "Warning: GC doesn't work with PBC yet; centers will be computed w/o PBC, distances with PBC\n");
            }
            if( inp_tilt->sele.nselatm < 2 ) {
                fprintf( stderr, "\nERROR! Selection %s must have more than 2 atoms (have %d)\n\n", inp_tilt->sele.selestring, inp_tilt->sele.nselatm);
                exit(0);
            }
            gotit = 1;
        }
        
        /* Florian 20 May 2015 */
        /* --DEC option (decompose) requires the polar decomposition of the angle */
        else if ( !strncmp(buffer, "--DEC", 5) ) {
            inp_tilt->decompose = 1;
            gotit = 1;
        }
        
        /* Florian 20 May 2015 */
        /* --NOORIENT option computes the tilt angle without caring about the relative orientation of the two selections */
        else if ( !strncmp(buffer, "--NOORIENT", 10) ) {
            inp_tilt->noorient = 1;
            gotit = 1;
        }
                
        if( gotit==0 ) {
            fprintf( stderr, "Could not understand option: %s\n", buffer);
            exit(5);
        }
        inp_index++;
        buffer=input[inp_index];
    }

    /* DEC and NOORIENT options are mutually exclusive! */
    if ( inp_tilt->decompose == 1 && inp_tilt->noorient == 1 ) {
        fprintf( stderr, "\nERROR! NOORIENT is incompatible with DEC.");
        exit(5);
    }

    /* Initialize Zprev to 0 0 1 */
    inp_tilt->Zprev[0]=0;
    inp_tilt->Zprev[1]=0;
    inp_tilt->Zprev[2]=1;

    /* Print header and return */
    if (inp_tilt->decompose == 1) {
        if (title[0]=='\0') {
    	    sprintf( printout, " %10s  %10s ", "theta", "phi");
        } else {
    	    sprintf( printout, " %4.2s-theta  %6.2s-phi ", title, title);
        }
        return 24;
    } else {
        if (title[0]=='\0') {
            sprintf( printout, " %10s ","theta");
        } else {
            sprintf( printout, " %4.2s-theta ", title);
        }
        return 12;
    }
}


/* Compute the angle from the dot product between two vectors */
float compute_angle(float dot) {

    float AC, AB;
    AC = dot;
    if(AC<0.9999999){
        AB = sqrt(1 - AC*AC);
        return atan2(AB, AC)*180.0/M_PI;
    }
    else{
        return 0.00;
    }
}

/*Compute moment of inertia of a selection of atoms */
float Compute_Inertia(CoorSet *trj_crd, Selection *sele, float ***mtx)
{
    int   jj;
    float  xx, yy, zz, xy, yz, xz;
    _eigen I;
    float com[3];  
    float *xcoor = trj_crd->xcoor;
    float *ycoor = trj_crd->ycoor;
    float *zcoor = trj_crd->zcoor;

    /* Compute the center of geometry of the selection */
    Calc_Com_NopbcNomass(com, trj_crd, *sele);

    if( sele->nselatm <= 1 ) {
        fprintf( stderr, "You must have more than one atom to compute the inertia moments"); 
        exit(0);
    }
    xx = yy = zz = xy = yz = xz = 0;
    for( jj=0; jj<sele->nselatm; jj++ ) {  
        xx += (xcoor[sele->selatm[jj]-1] - com[0]) * (xcoor[sele->selatm[jj]-1] - com[0]);
        yy += (ycoor[sele->selatm[jj]-1] - com[1]) * (ycoor[sele->selatm[jj]-1] - com[1]);
        zz += (zcoor[sele->selatm[jj]-1] - com[2]) * (zcoor[sele->selatm[jj]-1] - com[2]);
        xy += (xcoor[sele->selatm[jj]-1] - com[0]) * (ycoor[sele->selatm[jj]-1] - com[1]);
        yz += (ycoor[sele->selatm[jj]-1] - com[1]) * (zcoor[sele->selatm[jj]-1] - com[2]);
        xz += (xcoor[sele->selatm[jj]-1] - com[0]) * (zcoor[sele->selatm[jj]-1] - com[2]);
    }

    /* Set the intertia tensor */
    I.size = 3;
    I.inpmat = wrd_fmatAlloc( 3, 3 );
    
    I.inpmat[0][0] = yy + zz ; 
    I.inpmat[1][1] = xx + zz ;  
    I.inpmat[2][2] = xx + yy ;
    I.inpmat[0][1] = - xy ;
    I.inpmat[1][0] = - xy ;
    I.inpmat[0][2] = - xz ;
    I.inpmat[2][0] = - xz ;
    I.inpmat[1][2] = - yz ;
    I.inpmat[2][1] = - yz ;

    /* Diagonalize the inertia tensor */
    DiagMatrix( &I );

    *mtx =  I.eigvec;
    
    return 0;    
}


// Calculate le tilt angle from the two selections specified
    
int Compute_Tilt ( struct inp_tilt *inp_tilt, CoorSet *trj_crd, char *outstring)
{

    float theta, phi, scalar, angle;
    float *v1;
    float X[3], Y[3], *Z;
    float projv1X, projv1Y, projv1Z;
    float **mtxref, **mtx;
    float V_theta[3], V_phi[3], VX[3], VY[3], VZ[3];
    float com1[3], com2[3], com1com2[3]; 
    float com1com2proj;
    float endpoint1[3];                    // end-point vector for orientation check
    float endpoint1X, endpoint1Y,endpoint1Z; // coordinates of the end-point vector
    float endpoint2[3];                    
    float endpoint2X, endpoint2Y,endpoint2Z; 
    float endpoint_dot_eigenvec1; 
    float endpoint_dot_eigenvec2;         // scalar product of the eigenvector and the end-point vector
    
    
    Selection seleRef = inp_tilt->seleRef;
    Selection sele = inp_tilt->sele;

    
    /* Added by Florian: if the --DEC keyword is requested: */
     
    if ( inp_tilt->decompose == 1 ) {
     
        /* 1) find principal axes of rotation of SELEREF to use them as reference */

        Compute_Inertia(trj_crd, &inp_tilt->seleRef, &mtxref);

        Z = mtxref[2]; /* eigenvec with the highest eigenval */
        w_norm(Z);

        /* Check random flipping of the Z axis. */
        /* If the Z axis rotates more than 178 degree */
        /* Still not 100% accurate. Random flips are still possible */
        angle = compute_angle(dot_prod(Z, inp_tilt->Zprev));
        if (angle<-178 || angle>178) {
            Z[0]*=-1;
            Z[1]*=-1;
            Z[2]*=-1;
        }
        inp_tilt->Zprev[0] = Z[0];
        inp_tilt->Zprev[1] = Z[1];
        inp_tilt->Zprev[2] = Z[2];

        /* X defined as vector between COM of sele and COM of seleRef projected on the plan orthogonal to Z axis */

        Calc_Com_NopbcNomass(com1, trj_crd, inp_tilt->seleRef); 
        Calc_Com_NopbcNomass(com2, trj_crd, inp_tilt->sele);

        com1com2[0] = com2[0] - com1[0];
        com1com2[1] = com2[1] - com1[1];
        com1com2[2] = com2[2] - com1[2];

        com1com2proj = dot_prod(com1com2,Z);

        X[0] = com1com2[0] - com1com2proj*Z[0];
        X[1] = com1com2[1] - com1com2proj*Z[1];
        X[2] = com1com2[2] - com1com2proj*Z[2];

        w_norm(X);

        /* Y is then calculated upon Z and X */

        cross_prod (Z , X , Y);
        w_norm(Y);

        /* 2) find vector that bestfits SELE */

        Compute_Inertia(trj_crd, &inp_tilt->sele, &mtx);

        /* 3) find X Y and Z composante of vector v1 by projection */

        v1 = mtx[2];

        projv1X = dot_prod(v1,X);
        VX[0] = projv1X * X[0];
        VX[1] = projv1X * X[1];
        VX[2] = projv1X * X[2];

        projv1Y = dot_prod(v1,Y);
        VY[0] = projv1Y * Y[0];
        VY[1] = projv1Y * Y[1];
        VY[2] = projv1Y * Y[2];

        projv1Z = dot_prod(v1,Z);
        VZ[0] = projv1Z * Z[0];
        VZ[1] = projv1Z * Z[1];
        VZ[2] = projv1Z * Z[2];
        w_norm(VZ);

        /* theta = V - VY  (project v1 on XZ plane) */
        /* phi   = V - VX  (project v1 on YZ plane) */

        V_theta[0] = v1[0] - VY[0];
        V_theta[1] = v1[1] - VY[1];
        V_theta[2] = v1[2] - VY[2];
        w_norm(V_theta);

        V_phi[0]   = v1[0] - VX[0];
        V_phi[1]   = v1[1] - VX[1];
        V_phi[2]   = v1[2] - VX[2]; 
        w_norm(V_phi);

        /* ..compute theta */
        scalar = dot_prod(V_theta, VZ);
        theta = compute_angle(scalar);

        /* ..compute phi */
        scalar = dot_prod(V_phi, VZ);
        phi = compute_angle(scalar);

        /* Assign sign base on quadrants */
        if (projv1X<0) theta*=-1;
        if (projv1Y<0) phi*=-1;

        sprintf( outstring, " %10.5f  %10.5f ", theta, phi);

        return 24;
        
    } else if ( inp_tilt->decompose == 0) {

        /* if the polar/azimuthal decomposition is not requested,
         * we simply compute the angle between the principal axes of
         * the two selections.
         * In this case SELE and SELEREF are equivalent.
         */
     
        /* Two cases are to be distinguished:
         * 1. NOORIENT is selected: the relative orientation of 
         * the two selections is ignored and the returned angle is
         * comprised between 0 and 90 degrees.
         *
         * 2. Default behaviour: the overall orientation of each selection
         * is first estimated using the end-point vector; if necessary the
         * orientation of the principal eigenvector is modified to match 
         * the orientation of the end-point vector. It is useful when 
         * the groups of atoms have a natural orientation, such as alpha
         * helices or beta-strands. The returned angle is comprised between
         * 0 and 180 degrees.
         */
      
        if ( inp_tilt->noorient == 1 ) {
      
     
            /* 1. Principal axis of SELEREF */
            Compute_Inertia(trj_crd, &inp_tilt->seleRef, &mtxref);

            Z = mtxref[2]; /* eigenvec with the highest eigenval */
            w_norm(Z);     /* normalized to an unitary vector */
            
            /* 2. Principal axis of SELE */
            Compute_Inertia(trj_crd, &inp_tilt->sele, &mtx);
            
            v1 = mtx[2]; /* principal axis of SELE */
            w_norm(v1);  /* normalization */
            
            /* 3. Angle Computation */
            scalar = dot_prod(Z,v1); 
            theta  = compute_angle(scalar);
            if (theta > 90.0) {theta =  180.0 - theta;}

            /* 4. Write result ! */
            sprintf( outstring, " %10.5f ", theta);
            return 12;
        }
        
        else if ( inp_tilt->noorient == 0) {

            /* 1. Principal eigenvector of SELEREF */
            Compute_Inertia(trj_crd, &inp_tilt->seleRef, &mtxref);
            
            Z = mtxref[2]; /* eigenvec with the highest eigenval */
            w_norm(Z);     /* normalized to an unitary vector */
            
            /* Orientation check: */
            /*   end-point vector: */
            endpoint1X = trj_crd->xcoor[seleRef.selatm[seleRef.nselatm-1]] - trj_crd->xcoor[seleRef.selatm[0]];
            endpoint1Y = trj_crd->ycoor[seleRef.selatm[seleRef.nselatm-1]] - trj_crd->ycoor[seleRef.selatm[0]];
            endpoint1Z = trj_crd->zcoor[seleRef.selatm[seleRef.nselatm-1]] - trj_crd->zcoor[seleRef.selatm[0]];
            
            endpoint1[0] = endpoint1X;
            endpoint1[1] = endpoint1Y;
            endpoint1[2] = endpoint1Z;
            
            w_norm(endpoint1);
            
            /* Projection on the eigenvector */
            endpoint_dot_eigenvec1 = dot_prod(endpoint1,Z);
            
            /* If the projection is negative, the eigenvector has to be reversed */
            if ( endpoint_dot_eigenvec1 < 0. ) {
                Z[0] = -Z[0];
                Z[1] = -Z[1];
                Z[2] = -Z[2];
            }
            
            /* And we carry on ! */
            
            /* 2. Principal axis of SELE */
            Compute_Inertia(trj_crd, &inp_tilt->sele, &mtx);
            
            v1 = mtx[2]; /* principal axis of SELE */
            w_norm(v1);  /* normalization */
            
            /* Orientation check: */
            /*   end-point vector: */
            endpoint2X = trj_crd->xcoor[sele.selatm[sele.nselatm-1]] - trj_crd->xcoor[sele.selatm[0]];
            endpoint2Y = trj_crd->ycoor[sele.selatm[sele.nselatm-1]] - trj_crd->ycoor[sele.selatm[0]];
            endpoint2Z = trj_crd->zcoor[sele.selatm[sele.nselatm-1]] - trj_crd->zcoor[sele.selatm[0]];
            
            endpoint2[0] = endpoint2X;
            endpoint2[1] = endpoint2Y;
            endpoint2[2] = endpoint2Z;
            
            w_norm(endpoint2);
            
            /* Projection on the eigenvector */
            endpoint_dot_eigenvec2 = dot_prod(endpoint2,v1);
            
            /* If the projection is negative, the eigenvector has to be reversed */
            if ( endpoint_dot_eigenvec2 < 0. ) {
                v1[0] = -v1[0];
                v1[1] = -v1[1];
                v1[2] = -v1[2];
            }
            
            /* Let's go on. */
            
            /* 3. Angle Computation */
            scalar = dot_prod(Z,v1);
            theta = compute_angle(scalar);
            
            /* 4. Write result ! */
            sprintf( outstring, " %10.5f ", theta);
        
            return 12;
        }    
    }
}

