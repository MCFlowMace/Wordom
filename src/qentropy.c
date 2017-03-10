// ------------------------------------------------------------------
// Copyright (C) 2003  University of Zurich
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ------------------------------------------------------------------
/*! \file geometry.h
 \brief Simple geometrical analyses module
 
 Module for simple geometrical analyses, such as distances, contacts,
 angles, dihedrals, Root Mean Square Deviation (RMSD), Distance Root
 Mean Square Deviation (DRMS) etc.
*/
#include "wordom.h"
#include "tools.h"
#include "analysis.h"
#include "qentropy.h"

// ------------------------------------------------------------------
// Qh Entropy module
// ------------------------------------------------------------------
int Read_iEntropy ( char **input, int inp_intex, struct inp_pca *inp_pca, char *printout, Molecule *molecule )
{
   int   ii, jj, kk, ll, gotit;
   char  buffer[5120];
   char  title[64];
   
   extern short int     no_frame_par;
   no_frame_par = 1;

   inp_pca->T = 300;
   
   memset ( title, '\0', sizeof(title));
   memset ( buffer, '\0', sizeof(buffer));
   while( strncmp (buffer, "END", 3))
   {
     gotit = 0;
//   fgets(buffer, 512, inpfile);
     sprintf( buffer, "%s", input[inp_intex]);
//   if ( (buffer[0] != '$') && (buffer[0] != '#') && (buffer[0] != '\n') )
     if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
       gotit = 1;
     else if ( !strncmp(buffer, "--TITLE", 7))
     {
       sscanf( buffer, "--TITLE %s", title );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--SELE", 6) )
     {
       sscanf(buffer, "--SELE %[^\n]%*c", inp_pca->sele.selestring);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--TEMP", 6))
     {
       sscanf( buffer, "--TEMP %f", &inp_pca->T);
       gotit = 1;
     }
     if( gotit==0 )
     {
       fprintf( stderr, "Could not understand option: %s\n", buffer);
       exit(5);
     }
     inp_intex++;
   }
   
   GetSele ( inp_pca->sele.selestring, &inp_pca->sele, molecule);
   sscanf(title, "%s", inp_pca->title);
   
   inp_pca->xavgcoor = (float *)walloc( inp_pca->sele.nselatm, sizeof(float));
   inp_pca->yavgcoor = (float *)walloc( inp_pca->sele.nselatm, sizeof(float));
   inp_pca->zavgcoor = (float *)walloc( inp_pca->sele.nselatm, sizeof(float));
   inp_pca->reference = (float **)walloc( inp_pca->sele.nselatm, sizeof(float*));
   for (ii=0;ii<inp_pca->sele.nselatm;ii++)
     inp_pca->reference[ii] = (float *)walloc( 3, sizeof(float));
    
   ll=0;
   for(ii=0; ii<molecule->nSeg; ii++)
   {
     for(jj=0; jj<molecule->segment[ii].nRpS; jj++)
     {
       for(kk=0; kk<molecule->segment[ii].pRes[jj].nApR; kk++)
       {
         if(inp_pca->sele.selatm[ll] == molecule->segment[ii].pRes[jj].pAto[kk].atomn)
         {
           inp_pca->xavgcoor[ll] = 0.;
           inp_pca->yavgcoor[ll] = 0.;
           inp_pca->zavgcoor[ll] = 0.;
	   inp_pca->reference[ll][0] = molecule->segment[ii].pRes[jj].pAto[kk].xCoor;
	   inp_pca->reference[ll][1] = molecule->segment[ii].pRes[jj].pAto[kk].yCoor;
	   inp_pca->reference[ll][2] = molecule->segment[ii].pRes[jj].pAto[kk].zCoor;
           ll++;
           if(ll>=inp_pca->sele.nselatm)
            break;
         }
       }
       if(ll>=inp_pca->sele.nselatm)
        break;
     }
     if(ll>=inp_pca->sele.nselatm)
      break;
   }
    
   inp_pca->cov = (double **)walloc( inp_pca->sele.nselatm*3, sizeof(double *));
   for (ii=0; ii<inp_pca->sele.nselatm*3; ii++)
   {
     inp_pca->cov[ii] = (double *)walloc( inp_pca->sele.nselatm*3, sizeof(double));
     for ( jj=0; jj<inp_pca->sele.nselatm*3; jj++)
      inp_pca->cov[ii][jj] = 0.0;
   }

   return 0;
}
// ------------------------------------------------------------------
int Compute_Entropy ( struct inp_pca *inp_pca, struct sopt *OPT, CoorSet *dcd_crd, char *outprint )
{
   int ii, jj, kk, ll;
   static float     **coor, *scoo;
   float          rotmatrix[3][3], transvec[3];
   float          aa, bb, cc;
   
   if(!coor){
     coor = (float **)walloc ( inp_pca->sele.nselatm, sizeof(float *)); 
     for( ii=0; ii<inp_pca->sele.nselatm; ii++)
      {
       coor[ii] = (float *)walloc ( 3, sizeof(float));
       for ( jj=0; jj<3; jj++)
        coor[ii][jj] = 0.0;
      }
    }
   
   for( ii=0; ii<inp_pca->sele.nselatm; ii++)
    {
     coor[ii][0] = dcd_crd->xcoor[ (inp_pca->sele.selatm[ii])-1 ];
     coor[ii][1] = dcd_crd->ycoor[ (inp_pca->sele.selatm[ii])-1 ];
     coor[ii][2] = dcd_crd->zcoor[ (inp_pca->sele.selatm[ii])-1 ];
    }
   //Superimpose the structure to the reference 
   superimpose( inp_pca->sele.nselatm , inp_pca->reference, coor , rotmatrix , transvec );
   
   //Check that the supreimposition worked otherwise jump to next frame
   for (ii=0;ii<3;ii++){
     if(isnan(transvec[ii])) {inp_pca->missed_frames++; return 0;}
     for (jj=0;jj<3;jj++) if(isnan(rotmatrix[ii][jj])) {inp_pca->missed_frames++; return 0;}
   }
     
   for ( ii=0; ii<inp_pca->sele.nselatm; ii++)
     {
       aa = coor[ii][0]*rotmatrix[0][0] + coor[ii][1]*rotmatrix[0][1] + coor[ii][2]*rotmatrix[0][2] ;
       bb = coor[ii][0]*rotmatrix[1][0] + coor[ii][1]*rotmatrix[1][1] + coor[ii][2]*rotmatrix[1][2] ;
       cc = coor[ii][0]*rotmatrix[2][0] + coor[ii][1]*rotmatrix[2][1] + coor[ii][2]*rotmatrix[2][2] ;
       coor[ii][0] = aa + transvec[0]; //+ xcenter
       coor[ii][1] = bb + transvec[1]; //+ ycenter
       coor[ii][2] = cc + transvec[2]; //+ zcenter
     }
   //Add the rotated structure to the average
    
   for ( ii=0; ii<inp_pca->sele.nselatm; ii++)    inp_pca->xavgcoor[ii] += ( coor[ii][0] - inp_pca->reference[ii][0] );
   for ( ii=0; ii<inp_pca->sele.nselatm; ii++)    inp_pca->yavgcoor[ii] += ( coor[ii][1] - inp_pca->reference[ii][1] );
   for ( ii=0; ii<inp_pca->sele.nselatm; ii++)    inp_pca->zavgcoor[ii] += ( coor[ii][2] - inp_pca->reference[ii][2] );
   
   if(!scoo) scoo = (float *)walloc( inp_pca->sele.nselatm*3, sizeof(float));
   
   for ( ii=0; ii<inp_pca->sele.nselatm; ii++){
       scoo[ii*3+0] = ( coor[ii][0] - inp_pca->reference[ii][0] ) ;
       scoo[ii*3+1] = ( coor[ii][1] - inp_pca->reference[ii][1] ) ;
       scoo[ii*3+2] = ( coor[ii][2] - inp_pca->reference[ii][2] ) ;
   }
    

   for ( kk=0; kk<inp_pca->sele.nselatm*3; kk++ )
     for ( ll=0; ll<inp_pca->sele.nselatm*3; ll++ )
       inp_pca->cov[kk][ll] += ((double) scoo[kk]) * ((double) scoo[ll]) ;

   return 0;
}
// ------------------------------------------------------------------
int Post_Entropy ( struct inp_pca *inp_pca, int nframe , Molecule * molecule )
{
   FILE  *output;
   int    ii, jj;
   char   filename[128];
   float *coor , * mass, amass;  

    // SSYEV variables
    char JOBZ = 'V';
    char UPLO = 'U';

    int No;
    int LDA;
    int INFO;

    float *W;
    float *A;
    float *WORK;
    float LWORK;
   
   float *omega,*Bw;
   double k,kT,R,NA,hbar,Sho,Eho,Cvho;
   float T;
   
   
   nframe-=inp_pca->missed_frames;
   
   coor = (float *)walloc( inp_pca->sele.nselatm*3, sizeof(float)) ;
   omega = (float *)walloc( inp_pca->sele.nselatm*3, sizeof(float)) ;
   Bw = (float *)walloc( inp_pca->sele.nselatm*3, sizeof(float)) ;
   
   //normalize the average of the superimposed structures
   for ( ii=0; ii<inp_pca->sele.nselatm; ii++)    inp_pca->xavgcoor[ii] /= nframe ;
   for ( ii=0; ii<inp_pca->sele.nselatm; ii++)    inp_pca->yavgcoor[ii] /= nframe ;
   for ( ii=0; ii<inp_pca->sele.nselatm; ii++)    inp_pca->zavgcoor[ii] /= nframe ;

   for ( ii=0; ii<inp_pca->sele.nselatm; ii++){
       coor[ii*3+0] = inp_pca->xavgcoor[ii]  ;
       coor[ii*3+1] = inp_pca->yavgcoor[ii]  ;
       coor[ii*3+2] = inp_pca->zavgcoor[ii]  ;
   }
   
   for ( ii=0; ii<inp_pca->sele.nselatm; ii++)
    {
     molecule->coor.xcoor[inp_pca->sele.selatm[ii]-1] = inp_pca->xavgcoor[ii] + inp_pca->reference[ii][0];
     molecule->coor.ycoor[inp_pca->sele.selatm[ii]-1] = inp_pca->yavgcoor[ii] + inp_pca->reference[ii][1];
     molecule->coor.zcoor[inp_pca->sele.selatm[ii]-1] = inp_pca->zavgcoor[ii] + inp_pca->reference[ii][2];
   }
   
   mass = (float *)walloc ( inp_pca->sele.nselatm*3, sizeof(float));
   for ( ii=0; ii<inp_pca->sele.nselatm; ii++)    {
     amass=molecule->rawmol.bFac[inp_pca->sele.selatm[ii]-1];

     amass=sqrt(amass);
     mass[ii*3+0] = amass ;
     mass[ii*3+1] = amass ;
     mass[ii*3+2] = amass ;
   }
   

   sprintf( filename, "%s-avg.pdb", inp_pca->title);
   WritePdb_unstr ( filename , molecule ) ;

   
   sprintf( filename, "%s-matrix.txt", inp_pca->title);
   output = O_File ( filename, "w");
   
   for( ii=0; ii<inp_pca->sele.nselatm*3; ii++)
    {
     for( jj=0; jj<inp_pca->sele.nselatm*3; jj++)
      {
       inp_pca->cov[ii][jj] /= nframe;
       inp_pca->cov[ii][jj] -= ((double) coor[ii]) * ((double) coor[jj]) ;
       inp_pca->cov[ii][jj] *= mass[ii] * mass[jj] ;
       fprintf(output, "%7.3f ", inp_pca->cov[ii][jj]);
      }
     fprintf( output, "\n");
    }
   
   fclose(output);
   memset ( filename, '\0', sizeof(filename));
   
    No = (inp_pca->sele.nselatm)*3;
    LDA = (inp_pca->sele.nselatm)*3;
    LWORK = LDA*No;
    WORK=calloc(LWORK,sizeof(float));
    W=calloc(No,sizeof(float));
    INFO=0;

    A = (float *)walloc ( (inp_pca->sele.nselatm*3)*(inp_pca->sele.nselatm*3),  sizeof(float));
    for( ii=0; ii<inp_pca->sele.nselatm*3; ii++)
      for ( jj=0; jj<inp_pca->sele.nselatm*3; jj++)
        A[ii*(inp_pca->sele.nselatm*3)+jj] = inp_pca->cov[jj][ii];
    
    #ifdef LAPACK
     ssyev_(&JOBZ, &UPLO, &No, A, &LDA, W, WORK, &LWORK, &INFO);
    #else
     fprintf( stderr, "Sorry, all modules requiring matrix diagonalization with lapack libraries\n have been disabled in this executable\n");
     exit(0);
    #endif
   
   sprintf( filename, "%s-eigval.txt", inp_pca->title);
   output = O_File ( filename, "w");
   hbar=1.05457148e-34;                    // Planck's hbar constant
   k=1.3806504e-23;                        // Boltzmann's constant   J / K
   NA=6.0221418e23;                        // Avogadro's Number      1 / mol
   R=k*NA;                                 // Gas constant           J / K / mol
   T=inp_pca->T;                           //Absolute temperature (provided as input in kelvin)
   kT=k*T;
   for( ii=6; ii<inp_pca->sele.nselatm*3; ii++) omega[ii]=sqrt(R*T*1e23/W[ii]);  // W is in Angstrom**2 * amu=Angstrom**2 * g / mol while R is in J /K /mol
   for( ii=6; ii<inp_pca->sele.nselatm*3; ii++) Bw[ii]=hbar/kT*omega[ii];
   for( Sho=0.,ii=6; ii<inp_pca->sele.nselatm*3; ii++) Sho+=R*(Bw[ii]/(exp(Bw[ii])-1)-log(1-exp(-Bw[ii])));
   for( Eho=0.,ii=6; ii<inp_pca->sele.nselatm*3; ii++) Eho+=hbar*omega[ii]*(0.5+1./(exp(Bw[ii]-1)));
   for( Cvho=0.,ii=6; ii<inp_pca->sele.nselatm*3; ii++) Cvho+=R*Bw[ii]*Bw[ii]*exp(Bw[ii])/(exp(Bw[ii])-1)/(exp(Bw[ii])-1);
   fprintf(output, "Sho = %e J/K/mol, Eho = %e J/mol, Cvho = %e J/K/mol\n",Sho,Eho*NA,Cvho);
   
   for( ii=0; ii<inp_pca->sele.nselatm*3; ii++)
    fprintf(output, "%d\t%7.6f\n", ii+1, W[ii]);
   
   fclose(output);
   memset ( filename, '\0', sizeof(filename));
   
   sprintf( filename, "%s-eigvec.txt", inp_pca->title);
   output = O_File ( filename, "w");
   for( ii=0; ii<inp_pca->sele.nselatm*3; ii++)
    {
     for( jj=0; jj<inp_pca->sele.nselatm*3; jj++)
      fprintf(output, "%9.5f", A[ii+jj*(inp_pca->sele.nselatm*3)]);
     fprintf(output, "\n");
    }
   fclose(output);
   
   return 0;
}
// ------------------------------------------------------------------
