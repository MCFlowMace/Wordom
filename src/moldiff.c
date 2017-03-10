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
/*! \file moldiff.c
 \brief configuration similarity and variability measures
 
 Several methods to compute similarity between configurations and
 conformational variability along trajectories.
*/
#include "wordom.h"
#include "fileio.h"
#include "tools.h"
#include "datahandler.h"
#include "analysis.h"
#include "moldiff.h"
#include "qcprot.h"

// ------------------------------------------------------------------
// Rmsd module
// ------------------------------------------------------------------
int Read_iRmsd ( char **input, int inp_intex, struct inp_rms *inp_rms, char *printout, Molecule *molecule, CoorSetData *coorsetdata, sopt *OPT )
{
   int         ii;
   char        buffer[1024];
   char        title[64];
   int         gotit;
   extern short int     frame_par;
   frame_par = 1;
   
   memset ( title, '\0', sizeof(title));
   memset ( buffer, '\0', sizeof(buffer));

   inp_rms->trjwrite    = 0;
   inp_rms->progressive = 0;
   inp_rms->super       = 1;
   inp_rms->extout      = 0;
   sprintf( inp_rms->sele.selestring, "/*/*/*");
   while( strncmp (buffer, "END", 3))
   {
    gotit = 0;
    sprintf( buffer, "%s", input[inp_intex]);
    if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
      gotit = 1;
    else if ( !strncmp(buffer, "--SELE", 6))
    {
      memset( inp_rms->sele.selestring, '\0', 16);
      sscanf(buffer, "--SELE %[^\n]%*c ", inp_rms->sele.selestring);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--TITLE", 7))
    {
     sscanf( buffer, "--TITLE %s", title);
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--NOSUPER", 9))
    {
     inp_rms->super=0;
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--PROGRESSIVE", 13))
    {
     inp_rms->progressive=1;
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--EXTOUT", 8))
    {
     inp_rms->extout=1;
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--FIT", 5))
    {
     sscanf(buffer, "--FIT %[^\n]%*c ", inp_rms->fitsele.selestring);
     inp_rms->fit=1;
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--TRJOUT", 8))
    {
     inp_rms->trjwrite=1;
     inp_rms->outtraj = malloc(sizeof(Traj));
     sscanf( buffer, "--TRJOUT %s", inp_rms->trjout);
     OpenTrj ( inp_rms->trjout, inp_rms->outtraj, "w" );
     CopyTrjHeader( coorsetdata->traj.hdr, inp_rms->outtraj->hdr );
     inp_rms->outtraj->hdr->sixtyfour = 0;
     if( OPT->NOPBC_FLAG == 1 )
       inp_rms->outtraj->hdr->varpbc = 0;
     WriteTrjHeader ( inp_rms->outtraj, "all" );
     
     CalloCoor( &inp_rms->crd, molecule->nato);
     
     // if pbc are present provide 
     if( inp_rms->outtraj->hdr->varpbc )
     {
       inp_rms->crd.pbc_flag = 1;
       inp_rms->crd.pbc = calloc( 1, sizeof(Pbc));
     }
     
     gotit = 1;
    }
    
    if( gotit==0 )
    {
     fprintf( stderr, "Sorry, could not understand option: %s\n", buffer);
     exit(5);
    }
    inp_intex++;
   }

   GetSele ( inp_rms->sele.selestring, &inp_rms->sele, molecule);
   if( inp_rms->sele.nselatm < 3 )
   {
     fprintf( stderr, "Rmsd needs at the very least 3 atoms to work predictably\n");
     exit(0);
   }
   
   inp_rms->reference    = calloc (   inp_rms->sele.nselatm, sizeof (float * ));
   inp_rms->reference[0] = calloc ( 3*inp_rms->sele.nselatm, sizeof ( float  ));
   for ( ii =0; ii<inp_rms->sele.nselatm; ii++)
    inp_rms->reference[ii] = inp_rms->reference[0] + 3*ii;

   for ( ii=0; ii<inp_rms->sele.nselatm; ii++)
    {
     inp_rms->reference[ii][0] = molecule->coor.xcoor[inp_rms->sele.selatm[ii]-1];
     inp_rms->reference[ii][1] = molecule->coor.ycoor[inp_rms->sele.selatm[ii]-1];
     inp_rms->reference[ii][2] = molecule->coor.zcoor[inp_rms->sele.selatm[ii]-1];
    }

   inp_rms->moving    = calloc (   inp_rms->sele.nselatm, sizeof (float * ));
   inp_rms->moving[0] = calloc ( 3*inp_rms->sele.nselatm, sizeof ( float  ));
   for ( ii =0; ii<inp_rms->sele.nselatm; ii++)
    inp_rms->moving[ii] = inp_rms->moving[0] + 3*ii;

   // if FIT is called, prapare also space for alt fitting
   
   if( inp_rms->fit )
   {
    GetSele ( inp_rms->fitsele.selestring, &inp_rms->fitsele, molecule);
   
    inp_rms->fitref    = calloc (   inp_rms->fitsele.nselatm, sizeof (float * ));
    inp_rms->fitref[0] = calloc ( 3*inp_rms->fitsele.nselatm, sizeof ( float  ));
    for ( ii =0; ii<inp_rms->fitsele.nselatm; ii++)
     inp_rms->fitref[ii] = inp_rms->fitref[0] + 3*ii;

    for ( ii=0; ii<inp_rms->fitsele.nselatm; ii++)
    {
     inp_rms->fitref[ii][0] = molecule->coor.xcoor[inp_rms->fitsele.selatm[ii]-1];
     inp_rms->fitref[ii][1] = molecule->coor.ycoor[inp_rms->fitsele.selatm[ii]-1];
     inp_rms->fitref[ii][2] = molecule->coor.zcoor[inp_rms->fitsele.selatm[ii]-1];
    }

    inp_rms->fitmov    = calloc (   inp_rms->fitsele.nselatm, sizeof (float * ));
    inp_rms->fitmov[0] = calloc ( 3*inp_rms->fitsele.nselatm, sizeof ( float  ));
    for ( ii =0; ii<inp_rms->fitsele.nselatm; ii++)
     inp_rms->fitmov[ii] = inp_rms->fitmov[0] + 3*ii;
   }
   
   sprintf( printout, " %10s ", title);
   
   return 12;
}
// ------------------------------------------------------------------
int Compute_Rmsd ( struct inp_rms *inp_rms, struct sopt *OPT, CoorSet *trj_crd, char *output )
{
   int            ii;
   float          rotmatrix[3][3], transvec[3];
   float          aa, bb, cc;
   float          dist=0.0;


   for ( ii=0; ii<inp_rms->sele.nselatm; ii++)
   {
    inp_rms->moving[ii][0] = trj_crd->xcoor[inp_rms->sele.selatm[ii]-1];
    inp_rms->moving[ii][1] = trj_crd->ycoor[inp_rms->sele.selatm[ii]-1];
    inp_rms->moving[ii][2] = trj_crd->zcoor[inp_rms->sele.selatm[ii]-1];
   }
   if ( inp_rms->super )
   {
    if( inp_rms->fit )
    {
     for ( ii=0; ii<inp_rms->fitsele.nselatm; ii++)
     {
      inp_rms->fitmov[ii][0] = trj_crd->xcoor[inp_rms->fitsele.selatm[ii]-1];
      inp_rms->fitmov[ii][1] = trj_crd->ycoor[inp_rms->fitsele.selatm[ii]-1];
      inp_rms->fitmov[ii][2] = trj_crd->zcoor[inp_rms->fitsele.selatm[ii]-1];
     }
     superimpose ( inp_rms->fitsele.nselatm, inp_rms->fitref, inp_rms->fitmov, rotmatrix, transvec );
    }
    else
     superimpose ( inp_rms->sele.nselatm, inp_rms->reference, inp_rms->moving, rotmatrix, transvec );

    // apply rot_transl trovate
    for ( ii=0; ii<inp_rms->sele.nselatm; ii++)
     {
      aa = inp_rms->moving[ii][0]*rotmatrix[0][0] + inp_rms->moving[ii][1]*rotmatrix[0][1] + inp_rms->moving[ii][2]*rotmatrix[0][2] ;
      bb = inp_rms->moving[ii][0]*rotmatrix[1][0] + inp_rms->moving[ii][1]*rotmatrix[1][1] + inp_rms->moving[ii][2]*rotmatrix[1][2] ;
      cc = inp_rms->moving[ii][0]*rotmatrix[2][0] + inp_rms->moving[ii][1]*rotmatrix[2][1] + inp_rms->moving[ii][2]*rotmatrix[2][2] ;
      inp_rms->moving[ii][0] = aa + transvec[0]; //+ xcenter
      inp_rms->moving[ii][1] = bb + transvec[1]; //+ ycenter
      inp_rms->moving[ii][2] = cc + transvec[2]; //+ zcenter
     }
   }
   //printf("==============\n");
   //printf("%8.3f %8.3f %8.3f %8.3f\n", rotmatrix[0][0], rotmatrix[1][0], rotmatrix[2][0], transvec[0]);
   //printf("%8.3f %8.3f %8.3f %8.3f\n", rotmatrix[0][1], rotmatrix[1][1], rotmatrix[2][1], transvec[1]);
   //printf("%8.3f %8.3f %8.3f %8.3f\n", rotmatrix[0][2], rotmatrix[1][2], rotmatrix[2][2], transvec[2]);
   //printf("--------------\n");
   
   for ( ii=0; ii<inp_rms->sele.nselatm; ii++ )
    dist += ( (inp_rms->reference[ii][0]-inp_rms->moving[ii][0])*(inp_rms->reference[ii][0]-inp_rms->moving[ii][0]) + (inp_rms->reference[ii][1]-inp_rms->moving[ii][1])*(inp_rms->reference[ii][1]-inp_rms->moving[ii][1]) + (inp_rms->reference[ii][2]-inp_rms->moving[ii][2])*(inp_rms->reference[ii][2]-inp_rms->moving[ii][2]) );
   
   dist /= inp_rms->sele.nselatm;
   dist = sqrt ( dist );
   
   if ( inp_rms->trjwrite )
    {
     inp_rms->outtraj->hdr->nframe++;
     for ( ii=0; ii<inp_rms->outtraj->hdr->nato; ii++ )
      {
       inp_rms->crd.xcoor[ii] = trj_crd->xcoor[ii]*rotmatrix[0][0] + trj_crd->ycoor[ii]*rotmatrix[0][1] + trj_crd->zcoor[ii]*rotmatrix[0][2] + transvec[0] ;
       inp_rms->crd.ycoor[ii] = trj_crd->xcoor[ii]*rotmatrix[1][0] + trj_crd->ycoor[ii]*rotmatrix[1][1] + trj_crd->zcoor[ii]*rotmatrix[1][2] + transvec[1] ;
       inp_rms->crd.zcoor[ii] = trj_crd->xcoor[ii]*rotmatrix[2][0] + trj_crd->ycoor[ii]*rotmatrix[2][1] + trj_crd->zcoor[ii]*rotmatrix[2][2] + transvec[2] ;
      }

     if( trj_crd->pbc_flag == 1 && inp_rms->outtraj->hdr->varpbc == 1)
     {
       inp_rms->crd.pbc->a_size = trj_crd->pbc->a_size;
       inp_rms->crd.pbc->b_size = trj_crd->pbc->b_size;
       inp_rms->crd.pbc->c_size = trj_crd->pbc->c_size;
       inp_rms->crd.pbc->angle1 = trj_crd->pbc->angle1;
       inp_rms->crd.pbc->angle2 = trj_crd->pbc->angle2;
       inp_rms->crd.pbc->angle3 = trj_crd->pbc->angle3;
     }
     // pbc conditions are checked _after_ read_irms has been called so this must be done here
     else if( trj_crd->pbc_flag == -1 )// inefficient (once would be enough) but should work
     {
       inp_rms->crd.pbc_flag = -1;          
       inp_rms->outtraj->hdr->varpbc = 0;   
     }
     WriteTrjCoor ( &inp_rms->crd, inp_rms->outtraj );
    }
   
   if ( inp_rms->progressive )
    {
     for ( ii=0; ii<inp_rms->sele.nselatm; ii++)
      {
       inp_rms->reference[ii][0] = trj_crd->xcoor[inp_rms->sele.selatm[ii]-1];
       inp_rms->reference[ii][1] = trj_crd->ycoor[inp_rms->sele.selatm[ii]-1];
       inp_rms->reference[ii][2] = trj_crd->zcoor[inp_rms->sele.selatm[ii]-1];
      }
     
    }
   
   sprintf( output, " %10.6f ", dist ) ;
   
   return 12;
}
// ------------------------------------------------------------------
// ------------------------------------------------------------------
// Drms module
// ------------------------------------------------------------------
int Read_idrms ( char **input, int inp_intex, struct inp_drms *inp_drms, char *printout , Molecule *molecule , struct sopt * OPT )
{
   int         ii, jj, kk, ll;
   char        buffer[1024];
   char        title[64];
   CoorSet    *crd ;
   int         gotit;
   extern short int     frame_par;
   frame_par = 1;
   
   memset ( title, '\0', sizeof(title));
   memset ( buffer, '\0', sizeof(buffer));
   memset ( inp_drms->sele.selestring, '\0', sizeof(inp_drms->sele.selestring));

   inp_drms->progressive = 0;
   inp_drms->cmapd = 0;
   inp_drms->cmapd_cutoff = 0;
   inp_drms->selatm = NULL;
   inp_drms->nselatm = 0;
   while( strncmp (buffer, "END", 3))
   {
     gotit = 0;
     sprintf( buffer, "%s", input[inp_intex]);
     if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
       gotit = 1;
     else if ( !strncmp(buffer, "--TITLE", 7))
     {
       sscanf( buffer, "--TITLE %s", title);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--SELE", 6) )
     {
       sscanf(buffer, "--SELE %[^\n]%*c ", inp_drms->sele.selestring);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--PROGRESSIVE", 13))
     {
       inp_drms->progressive = 1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--CMAPD", 7))
     {
       inp_drms->cmapd = 1;
       sscanf(buffer, "--CMAPD %f", &(inp_drms->cmapd_cutoff));
       gotit = 1;
     }

     if( gotit==0 )
     {
       fprintf( stderr, "Could not understand option: %s\n", buffer);
       exit(5);
     }
     inp_intex++;
   }
   
   if( inp_drms->sele.selestring[0] == '\0' )
   {
     fprintf( stderr, "Missing selectin in drms module: please provide one\n");
     exit(0);
   }
   
   GetSele ( inp_drms->sele.selestring, &inp_drms->sele, molecule);
   inp_drms->nselatm = inp_drms->sele.nselatm;
   inp_drms->selatm  = calloc ( inp_drms->nselatm , sizeof ( int ) );
   for(ii=0; ii<inp_drms->sele.nselatm; ii++)
    inp_drms->selatm[ii] = inp_drms->sele.selatm[ii]-1;
   
   crd = InitCoor( molecule->nato );
   
   for(ii=0; ii<molecule->nSeg; ii++)
    for(jj=0; jj<molecule->segment[ii].nRpS; jj++)
     for(kk=0; kk<molecule->segment[ii].pRes[jj].nApR; kk++)
     {
      crd->xcoor[molecule->segment[ii].pRes[jj].pAto[kk].atomn-1] = molecule->segment[ii].pRes[jj].pAto[kk].xCoor;
      crd->ycoor[molecule->segment[ii].pRes[jj].pAto[kk].atomn-1] = molecule->segment[ii].pRes[jj].pAto[kk].yCoor;
      crd->zcoor[molecule->segment[ii].pRes[jj].pAto[kk].atomn-1] = molecule->segment[ii].pRes[jj].pAto[kk].zCoor;
     }
   
   ll = inp_drms->nselatm ; 
   kk=(int)(ll*(ll-1)/2);
   inp_drms->distancemat = calloc ( kk , sizeof (float));
   Calc_dist_mtx( OPT, inp_drms->nselatm , inp_drms->selatm, crd , inp_drms->distancemat );

   DelCoor(crd);

   sprintf( printout, " %10s ", title);
   return 12;
}
// ------------------------------------------------------------------
int Compute_drms ( struct inp_drms *inp_drms, struct sopt *OPT, CoorSet *dcd_crd, char *output )
{
    float                *dist_mtx=NULL;
    float                *cudi_mtx=NULL;
    float                 drms=0 ;
    int                   n ;
    int                   k;
    int                   ii;

    drms = 0.0;
    n = inp_drms->nselatm;
    k = (int)(n*(n-1)/2);
    dist_mtx = calloc ( k , sizeof (float));
    cudi_mtx = calloc ( k , sizeof (float));

    for( ii=0; ii<k; ii++ )
      cudi_mtx[ii] = inp_drms->distancemat[ii];

    Calc_dist_mtx ( OPT, n , inp_drms->selatm , dcd_crd , dist_mtx ) ;

    if( inp_drms->progressive )
      for( ii=0; ii<k; ii++ )
        inp_drms->distancemat[ii] = dist_mtx[ii];

    if( inp_drms->cmapd )
    {
      for( ii=0; ii<k; ii++ )
        dist_mtx[ii] = fret( dist_mtx[ii], inp_drms->cmapd_cutoff );
      for( ii=0; ii<k; ii++ )
        cudi_mtx[ii] = fret( cudi_mtx[ii], inp_drms->cmapd_cutoff );
	  }
    
    Drms ( k , cudi_mtx , dist_mtx , & drms ) ;

    free( dist_mtx );
    free( cudi_mtx );
    
    sprintf( output, " %10.6f ", drms ) ;
    
    return 12;
}
//-------------------------------------------------------------------
//Compute the fret efficiency (sigmoidal function) at distance d with
//r0 as forrester radius
float fret(float d, float r0)
{
  if( d==0.0 && r0==0.0 )
  {
  	fprintf(stderr,"fret: both d and r0 are null\n");
  	exit(EXIT_FAILURE);
  }
  d *= d*d;
  d *= d;
  r0 *= r0*r0;
  r0 *= r0;
  return (r0/(r0+d));
}
// ------------------------------------------------------------------
// Compute the distance matrix
//void Calc_dist_mtx( struct sopt *OPT, int  n , int * selatm , CoorSet *dcd_crd , float * dist_mtx )
int Calc_dist_mtx( struct sopt *OPT, int  n , int * selatm , CoorSet *dcd_crd , float * dist_mtx )
{
  int j,i,k;
  float **coor;

  coor = calloc (6, sizeof(float *));
  for (i=0;i<6;i++)
    coor[i] = calloc ( (int)(n*(n-1)/2) , sizeof(float));
  
  for(k=0,i=0;i<n;i++){
    for (j=0;j<i;j++){
     coor[0][k] = dcd_crd->xcoor[selatm[i]];   // Assign the coordinates of the relevant pairs of
     coor[1][k] = dcd_crd->ycoor[selatm[i]];   //  atoms to the vectors for distance calculation
     coor[2][k] = dcd_crd->zcoor[selatm[i]];
     coor[3][k] = dcd_crd->xcoor[selatm[j]];
     coor[4][k] = dcd_crd->ycoor[selatm[j]];
     coor[5][k] = dcd_crd->zcoor[selatm[j]];
      k++;
    }
  }
  
  // Compute the distance matrix as a vector having k = n * (n-1) / 2 components
//  if(OPT->PBC_FLAG)
//    vDistancePBC( k, dist_mtx, coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], OPT->PBCBOX); 
  if( OPT->PBC_FLAG )
    vDistancePBC( k, dist_mtx, coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], OPT->PBCBOX ); 
  else
    vDistance( k, dist_mtx, coor[0], coor[1], coor[2], coor[3], coor[4], coor[5]); 
  
  for (i=0;i<6;i++)
    free(coor[i]);
  free(coor);
  
  return k;
}
// ------------------------------------------------------------------

// ------------------------------------------------------------------
// Rmsf module
// ------------------------------------------------------------------
int Read_iRmsf ( char **input, int inp_index, struct inp_rmsf *inp_rmsf, char *printout, Molecule *molecule )
{
   int         ii;
   char        buffer[1024];
   char        title[64];
   int         gotit;
   
   extern short int     frame_par;
   frame_par = 0;
   
   memset ( title, '\0', sizeof(title));
   memset ( buffer, '\0', sizeof(buffer));

   inp_rmsf->super       = 1;
   while( strncmp (buffer, "END", 3))
   {
    gotit = 0;
    sprintf( buffer, "%s", input[inp_index]);
    if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
      gotit = 1;
    else if ( !strncmp(buffer, "--TITLE", 7))
    {
     sscanf( buffer, "--TITLE %s", inp_rmsf->title);
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--SELE", 6))
    {
      sscanf(buffer, "--SELE %[^\n]%*c ", inp_rmsf->sele.selestring);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--NOSUPER", 9))
    {
     inp_rmsf->super=0;
     gotit = 1;
     fprintf( stderr, "--NOSUPER is useless: trajectory is never aligned inside the RMSF module\n");
    }
    else if ( !strncmp(buffer, "--FIT", 5))
    {
     sscanf(buffer, "--FIT %[^\n]%*c ", inp_rmsf->fitsele.selestring);
     inp_rmsf->fit=1;
     fprintf( stderr, "Error: --FIT is useless since trajectory is never aligned inside the RMSF module\n");
     gotit = 1;
    }
    
    if( gotit==0 )
    {
     fprintf( stderr, "Sorry, could not understand option: %s\n", buffer);
     exit(5);
    }
    inp_index++;
   }
   
   if( inp_rmsf->sele.selestring == NULL )
   {
     fprintf( stderr, "You should really specify a selection\n");
     exit(5);
   }

   GetSele ( inp_rmsf->sele.selestring, &inp_rmsf->sele, molecule);
   
   inp_rmsf->xcmpcoor = (float *)walloc( inp_rmsf->sele.nselatm, sizeof(float));
   inp_rmsf->ycmpcoor = (float *)walloc( inp_rmsf->sele.nselatm, sizeof(float));
   inp_rmsf->zcmpcoor = (float *)walloc( inp_rmsf->sele.nselatm, sizeof(float));
   
   for ( ii=0; ii<inp_rmsf->sele.nselatm; ii++)
    {
     inp_rmsf->xcmpcoor[ii] = molecule->coor.xcoor[inp_rmsf->sele.selatm[ii]-1];
     inp_rmsf->ycmpcoor[ii] = molecule->coor.ycoor[inp_rmsf->sele.selatm[ii]-1];
     inp_rmsf->zcmpcoor[ii] = molecule->coor.zcoor[inp_rmsf->sele.selatm[ii]-1];
    }
   
   inp_rmsf->reference    = calloc (   inp_rmsf->sele.nselatm, sizeof (float * ));
   inp_rmsf->reference[0] = calloc ( 3*inp_rmsf->sele.nselatm, sizeof ( float  ));
   for ( ii =0; ii<inp_rmsf->sele.nselatm; ii++)
    inp_rmsf->reference[ii] = inp_rmsf->reference[0] + 3*ii;

   for ( ii=0; ii<inp_rmsf->sele.nselatm; ii++)
    {
     inp_rmsf->reference[ii][0] = inp_rmsf->xcmpcoor[ii];
     inp_rmsf->reference[ii][1] = inp_rmsf->ycmpcoor[ii];
     inp_rmsf->reference[ii][2] = inp_rmsf->zcmpcoor[ii];
    }

   inp_rmsf->outrmsf = calloc( inp_rmsf->sele.nselatm, sizeof (float ));
   inp_rmsf->framecounter = 0;
   
   return 0;
}
// ------------------------------------------------------------------
int Compute_Rmsf ( struct inp_rmsf *inp_rmsf, struct sopt *OPT, CoorSet *trj_crd, char *output )
{
  int          ii;
  
  for( ii=0; ii<inp_rmsf->sele.nselatm; ii++)
  {
    inp_rmsf->outrmsf[ii] += (inp_rmsf->xcmpcoor[ii] - trj_crd->xcoor[inp_rmsf->sele.selatm[ii]-1]) * (inp_rmsf->xcmpcoor[ii] - trj_crd->xcoor[inp_rmsf->sele.selatm[ii]-1]) +\
                             (inp_rmsf->ycmpcoor[ii] - trj_crd->ycoor[inp_rmsf->sele.selatm[ii]-1]) * (inp_rmsf->ycmpcoor[ii] - trj_crd->ycoor[inp_rmsf->sele.selatm[ii]-1]) +\
                             (inp_rmsf->zcmpcoor[ii] - trj_crd->zcoor[inp_rmsf->sele.selatm[ii]-1]) * (inp_rmsf->zcmpcoor[ii] - trj_crd->zcoor[inp_rmsf->sele.selatm[ii]-1]);
  }
  
  inp_rmsf->framecounter ++;
  
  return 0;
}
// ------------------------------------------------------------------
int Post_Rmsf ( struct inp_rmsf *inp_rmsf, struct sopt *OPT, int nframe, Molecule *molecule )
{
  int          ii;
  FILE        *rmsfout;
  
  rmsfout = O_File( inp_rmsf->title, "w");
  
  for( ii=0; ii<inp_rmsf->sele.nselatm; ii++)
    fprintf( rmsfout, "%7d %10.5f\n", inp_rmsf->sele.selatm[ii], sqrt(inp_rmsf->outrmsf[ii]/inp_rmsf->framecounter));
  
  return 0;
}
// ------------------------------------------------------------------
// Msfd module (old)
// ------------------------------------------------------------------
/*int Read_iMsdf_old ( char **input, int inp_index, struct inp_msdf *inp_msdf, char *printout, Molecule *molecule, int nframe )
{
  int         ii;
  char        buffer[1024];
  char        title[64];
  int         gotit;
  
  int         iResNum, iResCont;
  char        cTmpString[10];
  char        cLastRes[50];
  char        cResCode1[3];
  
  extern short int     frame_par;
  frame_par = 0;
  
  memset ( title, '\0', sizeof(title));
  memset ( buffer, '\0', sizeof(buffer));
  
  inp_msdf->root = 0;
  inp_msdf->verbose = 0;
  inp_msdf->iMultiAtomFlag = 0;
  
  while( strncmp (buffer, "END", 3))
  {
   gotit = 0;
   sprintf( buffer, "%s", input[inp_index]);
   if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
     gotit = 1;
   else if ( !strncmp(buffer, "--TITLE", 7))
   {
    sscanf( buffer, "--TITLE %s", inp_msdf->title);
    gotit = 1;
   }
   else if ( !strncmp(buffer, "--VERBOSE", 9))
   {
    inp_msdf->verbose = 1;
    gotit = 1;
   }
   else if ( !strncmp(buffer, "--SELE", 6))
   {
     sscanf(buffer, "--SELE %[^\n]%*c ", inp_msdf->sele.selestring);
     gotit = 1;
   }
   else if ( !strncmp(buffer, "--ROOT", 6))
   {
     inp_msdf->root = 1;
     gotit = 1;
   }
   
   if( gotit==0 )
   {
    fprintf( stderr, "Sorry, could not understand option: %s\n", buffer);
    exit(5);
   }
   inp_index++;
  }
  
  if( inp_msdf->sele.selestring == NULL )
  {
    fprintf( stderr, "You should really specify a selection\n");
    exit(5);
  }

  GetSele(inp_msdf->sele.selestring, &inp_msdf->sele, molecule);
  
  inp_msdf->distmat = (float ***)calloc( nframe, sizeof( float ** ));
  inp_msdf->framecounter = 0;
  
  // === General Settings ==============================================
  if(inp_msdf->verbose == 1)
  {
    inp_msdf->iNumOfSelAtoms = inp_msdf->sele.nselatm;
    inp_msdf->pcLabels    = (char **) calloc(inp_msdf->iNumOfSelAtoms, sizeof(char *));
    for(ii=0; ii<inp_msdf->iNumOfSelAtoms; ii++)
      inp_msdf->pcLabels[ii] = (char *) calloc(20, sizeof(char));
  }
  
  inp_msdf->iNumOfFrames=nframe;
  
  inp_msdf->iNumOfRes=molecule->nRes;
  inp_msdf->iNumOfAtoms=molecule->nato;

  inp_msdf->iNumOfSelRes=0;
  inp_msdf->iNumOfSelAtoms=inp_msdf->sele.nselatm;
  inp_msdf->iMultiAtomFlag=0;
  
  cTmpString[0]='\0';
  cLastRes[0]='\0';
  
  for(ii=0; ii<inp_msdf->sele.nselatm; ii++)
  {
    sprintf(cTmpString, "%s%d", molecule->rawmol.segId[inp_msdf->sele.selatm[ii]-1], 
                                molecule->rawmol.resn[inp_msdf->sele.selatm[ii]-1]);
    if(strcmp(cTmpString, cLastRes)!=0)
    {
      inp_msdf->iNumOfSelRes++;
      strcpy(cLastRes, cTmpString);
    }
    else
      inp_msdf->iMultiAtomFlag=1;                                      // More than 1 atom per residue
  }
  
  if(inp_msdf->iMultiAtomFlag == 1)                                    // Calculates the number of atoms per selected residue
  {
    inp_msdf->piSelResLen = (int   *) calloc(inp_msdf->iNumOfSelRes, sizeof(int   ));
    inp_msdf->piProgNum   = (int   *) calloc(inp_msdf->iNumOfSelRes, sizeof(int   ));
    inp_msdf->pcLabels    = (char **) calloc(inp_msdf->iNumOfSelRes, sizeof(char *));
    for(ii=0; ii<inp_msdf->iNumOfSelRes; ii++)
      inp_msdf->pcLabels[ii] = (char *) calloc(20, sizeof(char));
    
    cTmpString[0]='\0';
    cLastRes[0]='\0';
    sprintf(cLastRes, "%s%d", molecule->rawmol.segId[inp_msdf->sele.selatm[0]-1], 
                              molecule->rawmol.resn[inp_msdf->sele.selatm[0]-1]);
    
    iResNum=-1;
    iResCont=-1;

    Res3ToRes1(molecule->rawmol.restype[inp_msdf->sele.selatm[0]-1], cResCode1);
    sprintf(inp_msdf->pcLabels[0], "%s:%s%d", molecule->rawmol.segId[inp_msdf->sele.selatm[0]-1], cResCode1,
                                              molecule->rawmol.resn[inp_msdf->sele.selatm[0]-1]);
                                                        
    inp_msdf->piProgNum[0]=molecule->rawmol.presn[inp_msdf->sele.selatm[0]-1];
    
    for(ii=0; ii<inp_msdf->sele.nselatm; ii++)
    {
      iResCont++;
      sprintf(cTmpString, "%s%d", molecule->rawmol.segId[inp_msdf->sele.selatm[ii]-1], 
                                  molecule->rawmol.resn[inp_msdf->sele.selatm[ii]-1]);
      
      if(strcmp(cTmpString, cLastRes)!=0)
      {
        iResNum++;
        inp_msdf->piSelResLen[iResNum]=iResCont;
        strcpy(cLastRes, cTmpString);
        Res3ToRes1(molecule->rawmol.restype[inp_msdf->sele.selatm[ii]-1], cResCode1);
        sprintf(inp_msdf->pcLabels[iResNum+1], "%s:%s%d", molecule->rawmol.segId[inp_msdf->sele.selatm[ii]-1], cResCode1,
                                                          molecule->rawmol.resn[inp_msdf->sele.selatm[ii]-1]);
                                                        
        inp_msdf->piProgNum[iResNum+1]=molecule->rawmol.presn[inp_msdf->sele.selatm[ii]-1];
        iResCont=0;
      }
    }
    
    iResNum++;
    inp_msdf->piSelResLen[iResNum]=iResCont+1;
    inp_msdf->ppfVirtAtomCoord    = (float **) calloc(inp_msdf->iNumOfSelRes, sizeof(float *));
    for(ii=0; ii<inp_msdf->iNumOfSelRes; ii++)                         // Stores Geometrical centre coordinates
    {
      inp_msdf->ppfVirtAtomCoord[ii]    = (float *) calloc(3, sizeof(float));
    }
  }

  // =================================================================

  return 0;
}*/
// ------------------------------------------------------------------
/*int Compute_Msdf_old ( struct inp_msdf *inp_msdf, struct sopt *OPT, CoorSet *trj_crd, char *output )
{
  int          ii, jj, thisframe, nato, counter=0;
  
  // ===
  int          iResNum         = 0;
  int          iResCont        = 0;
  int          iVirtAtom       = -1;
  
  float        fVirtAtomXCoord = 0.0;
  float        fVirtAtomYCoord = 0.0;
  float        fVirtAtomZCoord = 0.0;
  
  inp_msdf->iPBCFlag = trj_crd->pbc_flag;
  // ===
  
  nato = inp_msdf->sele.nselatm;
  thisframe = inp_msdf->framecounter;
  
  inp_msdf->distmat[thisframe] = (float **)calloc( nato , sizeof(float *) );
  inp_msdf->distmat[thisframe][0] = (float *)calloc( (nato*nato-nato)/2 , sizeof(float ) );
  for( ii=1; ii<nato; ii++ )
  {
    inp_msdf->distmat[thisframe][ii] = inp_msdf->distmat[thisframe][0] + counter + ii - 1;
    counter += ii-1;
  }
  
  if(inp_msdf->iMultiAtomFlag == 0)
  {
    for( ii=0; ii<nato; ii++)
      for( jj=0; jj<ii; jj++ )
        inp_msdf->distmat[thisframe][ii][jj] = sqrt( (trj_crd->xcoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->xcoor[inp_msdf->sele.selatm[jj]-1] ) * (trj_crd->xcoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->xcoor[inp_msdf->sele.selatm[jj]-1] ) +
                                                      (trj_crd->ycoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->ycoor[inp_msdf->sele.selatm[jj]-1] ) * (trj_crd->ycoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->ycoor[inp_msdf->sele.selatm[jj]-1] ) +
                                                      (trj_crd->zcoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->zcoor[inp_msdf->sele.selatm[jj]-1] ) * (trj_crd->zcoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->zcoor[inp_msdf->sele.selatm[jj]-1] ) );
  }
  
  else
  {
    // More than one atom per residue
    
    // === Calculates geo centres of current frame res =================
    iResNum=0;
    iResCont=0;
    iVirtAtom=-1;
    
    fVirtAtomXCoord = 0.0;
    fVirtAtomYCoord = 0.0;
    fVirtAtomZCoord = 0.0;
    
    for(ii=0; ii<inp_msdf->sele.nselatm; ii++)
    {
      if(iResCont == inp_msdf->piSelResLen[iResNum])
      {
        fVirtAtomXCoord = (fVirtAtomXCoord/(float)inp_msdf->piSelResLen[iResNum]);
        fVirtAtomYCoord = (fVirtAtomYCoord/(float)inp_msdf->piSelResLen[iResNum]);
        fVirtAtomZCoord = (fVirtAtomZCoord/(float)inp_msdf->piSelResLen[iResNum]);

        iResNum++;
        iResCont=0;
        iVirtAtom++;
        
        inp_msdf->ppfVirtAtomCoord[iVirtAtom][0]=fVirtAtomXCoord;
        inp_msdf->ppfVirtAtomCoord[iVirtAtom][1]=fVirtAtomYCoord;
        inp_msdf->ppfVirtAtomCoord[iVirtAtom][2]=fVirtAtomZCoord;
        
        fVirtAtomXCoord = 0.0;
        fVirtAtomYCoord = 0.0;
        fVirtAtomZCoord = 0.0;
      }
      fVirtAtomXCoord += trj_crd->xcoor[inp_msdf->sele.selatm[ii]-1];
      fVirtAtomYCoord += trj_crd->ycoor[inp_msdf->sele.selatm[ii]-1];
      fVirtAtomZCoord += trj_crd->zcoor[inp_msdf->sele.selatm[ii]-1];
      
      iResCont++;
    }
    
    // Process the last residue
    fVirtAtomXCoord = (fVirtAtomXCoord/(float)inp_msdf->piSelResLen[iResNum]);
    fVirtAtomYCoord = (fVirtAtomYCoord/(float)inp_msdf->piSelResLen[iResNum]);
    fVirtAtomZCoord = (fVirtAtomZCoord/(float)inp_msdf->piSelResLen[iResNum]);

    iVirtAtom++;
    inp_msdf->ppfVirtAtomCoord[iVirtAtom][0]=fVirtAtomXCoord;
    inp_msdf->ppfVirtAtomCoord[iVirtAtom][1]=fVirtAtomYCoord;
    inp_msdf->ppfVirtAtomCoord[iVirtAtom][2]=fVirtAtomZCoord;

    for(ii=0; ii<inp_msdf->iNumOfSelRes; ii++)
    {
      for( jj=0; jj<ii; jj++ )
      {
        inp_msdf->distmat[thisframe][ii][jj] = sqrt( (inp_msdf->ppfVirtAtomCoord[ii][0] - inp_msdf->ppfVirtAtomCoord[jj][0] ) * (inp_msdf->ppfVirtAtomCoord[ii][0] - inp_msdf->ppfVirtAtomCoord[jj][0] ) +
                                                      (inp_msdf->ppfVirtAtomCoord[ii][1] - inp_msdf->ppfVirtAtomCoord[jj][1] ) * (inp_msdf->ppfVirtAtomCoord[ii][1] - inp_msdf->ppfVirtAtomCoord[jj][1] ) +
                                                      (inp_msdf->ppfVirtAtomCoord[ii][2] - inp_msdf->ppfVirtAtomCoord[jj][2] ) * (inp_msdf->ppfVirtAtomCoord[ii][2] - inp_msdf->ppfVirtAtomCoord[jj][2] ) );
      }
    }
  }

  
  inp_msdf->framecounter ++;
  
  return 0;
}*/
// ------------------------------------------------------------------
/*int Post_Msdf_old ( struct inp_msdf *inp_msdf, struct sopt *OPT, int nframe, Molecule *molecule )
{
  int          ii, jj;
  int          nato, ndist, counter=0;
  double     **meanmat;
  double     **fluctmat;
  FILE        *rmsfdout;
  
  char         cLabelA[20], cLabelB[20], cResCode1[3];
  int          iProgNumA, iProgNumB;
  double       dRMSFd;
  time_t       time_Today;
  
  nato = inp_msdf->iNumOfSelRes;
  
  meanmat = (double **)calloc( nato , sizeof(double *) );
  meanmat[0] = (double *)calloc( (nato*nato-nato)/2 , sizeof(double ) );
  for( ii=1; ii<nato; ii++ )
  {
    meanmat[ii] = meanmat[0] + counter + ii - 1;
    counter += ii-1;
  }
  counter = 0;
  fluctmat = (double **)calloc( nato , sizeof(double *) );
  fluctmat[0] = (double *)calloc( (nato*nato-nato)/2 , sizeof(double ) );
  for( ii=1; ii<nato; ii++ )
  {
    fluctmat[ii] = fluctmat[0] + counter + ii - 1;
    counter += ii-1;
  }
  
  ndist = ((nato*nato)-nato)/2;
  for( ii=0; ii<nframe; ii++)
   for( jj=0; jj<ndist; jj++ )
     meanmat[0][jj] += (double)inp_msdf->distmat[ii][0][jj]/nframe;
  
  for( ii=0; ii<nframe; ii++)
   for( jj=0; jj<ndist; jj++ )
     fluctmat[0][jj] += ((double)inp_msdf->distmat[ii][0][jj]-meanmat[0][jj])*((double)inp_msdf->distmat[ii][0][jj]-meanmat[0][jj])/nframe;
  
  if( inp_msdf->root )
    for( ii=0; ii<ndist; ii++)
      fluctmat[0][ii] = sqrt( fluctmat[0][ii] );
  
  rmsfdout = O_File( inp_msdf->title, "w");
  
  if(inp_msdf->verbose == 0)
  {
    fprintf( rmsfdout, "     ");
    for( ii=0; ii<nato; ii++)
      fprintf( rmsfdout, "%8d ", inp_msdf->sele.selatm[ii]);
    fprintf( rmsfdout, "\n");
    for( ii=0; ii<nato; ii++)
    {
      fprintf( rmsfdout, "%4d ", inp_msdf->sele.selatm[ii]);
      for( jj=0; jj<ii; jj++ )
        fprintf( rmsfdout, "%8.3f ", fluctmat[ii][jj] );
      fprintf( rmsfdout, "\n");
    }
  }
  
  else
  {
    time(&time_Today);
    fprintf(rmsfdout, "# =================================================================\n");
    fprintf(rmsfdout, "# ***                     WORDOM MSDF MODULE                    ***\n");
    fprintf(rmsfdout, "# =================================================================\n");
    fprintf(rmsfdout, "#\n");
    fprintf(rmsfdout, "# Version   : 0.1a\n");
    fprintf(rmsfdout, "# License   : GPL 3\n");
    fprintf(rmsfdout, "# Copyright : Michele Seeber, Angelo Felline, Francesca Fanelli\n");
    fprintf(rmsfdout, "#             University of Modena and Reggio Emilia\n");
    fprintf(rmsfdout, "#             Modena - Italy\n");
    fprintf(rmsfdout, "#\n");
    fprintf(rmsfdout, "# Date      : %s", asctime(localtime(&time_Today)));
    fprintf(rmsfdout, "#\n");
    fprintf(rmsfdout, "# Mol File  : %s\n", OPT->IMOL_FILE);
    fprintf(rmsfdout, "# Res Num   : %d\n", inp_msdf->iNumOfRes);
    fprintf(rmsfdout, "# Traj File : %s\n", OPT->ITRJ_FILE);
    fprintf(rmsfdout, "# Frame Num : %d\n", inp_msdf->iNumOfFrames);
    
    if(inp_msdf->iPBCFlag == 0)
      fprintf(rmsfdout, "# PBC       : No\n");
    else
      fprintf(rmsfdout, "# PBC       : Yes\n");
    
    fprintf(rmsfdout, "#\n");
    fprintf(rmsfdout, "# Title     : %s\n", inp_msdf->title);
    fprintf(rmsfdout, "# Sele      : %s\n", inp_msdf->sele.selestring);
    
    fprintf(rmsfdout, "# Sele Res  : %d\n", inp_msdf->iNumOfSelRes);
    
    if(inp_msdf->iMultiAtomFlag == 1)
      fprintf(rmsfdout, "# MultiAtom : True\n");
    else
      fprintf(rmsfdout, "# MultiAtom : False\n");
    
    fprintf(rmsfdout, "#\n");
    fprintf(rmsfdout, "# Warning!  : PBC not taken into account\n");
    fprintf(rmsfdout, "# =================================================================\n");
    
    fprintf(rmsfdout, "#%8s   %8s   %15s   %15s   %8s\n", "ProgNumA", "ProgNumB", "LabelA", "LabelB", "RMSFd");
    
    for(ii=0; ii<inp_msdf->iNumOfSelRes; ii++)
    {
      for(jj=0; jj<inp_msdf->iNumOfSelRes; jj++)
      {
        if(inp_msdf->iMultiAtomFlag==0)
        {
          Res3ToRes1(molecule->rawmol.restype[inp_msdf->sele.selatm[ii]-1], cResCode1);
          sprintf(cLabelA, "%s:%s%d", molecule->rawmol.segId[inp_msdf->sele.selatm[ii]-1], cResCode1,
                                      molecule->rawmol.resn[inp_msdf->sele.selatm[ii]-1]);
        
          Res3ToRes1(molecule->rawmol.restype[inp_msdf->sele.selatm[jj]-1], cResCode1);
          sprintf(cLabelB, "%s:%s%d", molecule->rawmol.segId[inp_msdf->sele.selatm[jj]-1], cResCode1,
                                      molecule->rawmol.resn[inp_msdf->sele.selatm[jj]-1]);
          
          iProgNumA=molecule->rawmol.presn[inp_msdf->sele.selatm[ii]-1];
          iProgNumB=molecule->rawmol.presn[inp_msdf->sele.selatm[jj]-1];
        }
        
        else
        {
          strcpy(cLabelA, inp_msdf->pcLabels[ii]);
          strcpy(cLabelB, inp_msdf->pcLabels[jj]);
          iProgNumA=inp_msdf->piProgNum[ii];
          iProgNumB=inp_msdf->piProgNum[jj];
        }
        
        if(ii == jj)
          dRMSFd = 0.0;
        
        else if(ii < jj)
          dRMSFd = fluctmat[jj][ii];
        
        else if (ii > jj)
          dRMSFd = fluctmat[ii][jj];

        fprintf(rmsfdout, " %8d   %8d   %15s   %15s   %8.3f\n", iProgNumA, iProgNumB, cLabelA, cLabelB, dRMSFd);
      }
    }
  }
  
  return 0;
}*/
// ------------------------------------------------------------------
// Msfd module
// ------------------------------------------------------------------
int Read_iMsdf ( char **input, int inp_index, struct inp_msdf *inp_msdf, char *printout, Molecule *molecule, int nframe )
{
  int         ii;
  char        buffer[1024];
  char        title[64];
  int         gotit;
  
  int         iResNum, iResCont, nato, counter;
  char        cTmpString[10];
  char        cLastRes[50];
  char        cResCode1[3];
  
  extern short int     frame_par;
  frame_par = 0;
  
  memset ( title, '\0', sizeof(title));
  memset ( buffer, '\0', sizeof(buffer));
  
  inp_msdf->pass = 1;
  inp_msdf->root = 0;
  inp_msdf->verbose = 0;
  inp_msdf->iMultiAtomFlag = 0;
  
  while( strncmp (buffer, "END", 3))
  {
   gotit = 0;
   sprintf( buffer, "%s", input[inp_index]);
   if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
     gotit = 1;
   else if ( !strncmp(buffer, "--TITLE", 7))
   {
    sscanf( buffer, "--TITLE %s", inp_msdf->title);
    gotit = 1;
   }
   else if ( !strncmp(buffer, "--VERBOSE", 9))
   {
    inp_msdf->verbose = 1;
    gotit = 1;
   }
   else if ( !strncmp(buffer, "--SELE", 6))
   {
     sscanf(buffer, "--SELE %[^\n]%*c ", inp_msdf->sele.selestring);
     gotit = 1;
   }
   else if ( !strncmp(buffer, "--ROOT", 6))
   {
     inp_msdf->root = 1;
     gotit = 1;
   }
   
   if( gotit==0 )
   {
    fprintf( stderr, "Sorry, could not understand option: %s\n", buffer);
    exit(5);
   }
   inp_index++;
  }
  
  if( inp_msdf->sele.selestring == NULL )
  {
    fprintf( stderr, "You should really specify a selection\n");
    exit(5);
  }

  GetSele(inp_msdf->sele.selestring, &inp_msdf->sele, molecule);
  
  inp_msdf->distmat = (float ***)calloc( nframe, sizeof( float ** ));
  inp_msdf->framecounter = 0;
  
  // === General Settings ==============================================
  if(inp_msdf->verbose == 1)
  {
    inp_msdf->iNumOfSelAtoms = inp_msdf->sele.nselatm;
    inp_msdf->pcLabels    = (char **) calloc(inp_msdf->iNumOfSelAtoms, sizeof(char *));
    for(ii=0; ii<inp_msdf->iNumOfSelAtoms; ii++)
      inp_msdf->pcLabels[ii] = (char *) calloc(20, sizeof(char));
  }
  
  inp_msdf->iNumOfFrames=nframe;
  
  inp_msdf->iNumOfRes=molecule->nRes;
  inp_msdf->iNumOfAtoms=molecule->nato;

  inp_msdf->iNumOfSelRes=0;
  inp_msdf->iNumOfSelAtoms=inp_msdf->sele.nselatm;
  inp_msdf->iMultiAtomFlag=0;
  
  cTmpString[0]='\0';
  cLastRes[0]='\0';
  
  for(ii=0; ii<inp_msdf->sele.nselatm; ii++)
  {
    sprintf(cTmpString, "%s%d", molecule->rawmol.segId[inp_msdf->sele.selatm[ii]-1], 
                                molecule->rawmol.resn[inp_msdf->sele.selatm[ii]-1]);
    if(strcmp(cTmpString, cLastRes)!=0)
    {
      inp_msdf->iNumOfSelRes++;
      strcpy(cLastRes, cTmpString);
    }
    else
      inp_msdf->iMultiAtomFlag=1;                                      // More than 1 atom per residue
  }
  
  if(inp_msdf->iMultiAtomFlag == 1)                                    // Calculates the number of atoms per selected residue
  {
    inp_msdf->piSelResLen = (int   *) calloc(inp_msdf->iNumOfSelRes, sizeof(int   ));
    inp_msdf->piProgNum   = (int   *) calloc(inp_msdf->iNumOfSelRes, sizeof(int   ));
    inp_msdf->pcLabels    = (char **) calloc(inp_msdf->iNumOfSelRes, sizeof(char *));
    for(ii=0; ii<inp_msdf->iNumOfSelRes; ii++)
      inp_msdf->pcLabels[ii] = (char *) calloc(20, sizeof(char));
    
    cTmpString[0]='\0';
    cLastRes[0]='\0';
    sprintf(cLastRes, "%s%d", molecule->rawmol.segId[inp_msdf->sele.selatm[0]-1], 
                              molecule->rawmol.resn[inp_msdf->sele.selatm[0]-1]);
    
    iResNum=-1;
    iResCont=-1;

    Res3ToRes1(molecule->rawmol.restype[inp_msdf->sele.selatm[0]-1], cResCode1);
    sprintf(inp_msdf->pcLabels[0], "%s:%s%d", molecule->rawmol.segId[inp_msdf->sele.selatm[0]-1], cResCode1,
                                              molecule->rawmol.resn[inp_msdf->sele.selatm[0]-1]);
                                                        
    inp_msdf->piProgNum[0]=molecule->rawmol.presn[inp_msdf->sele.selatm[0]-1];
    
    for(ii=0; ii<inp_msdf->sele.nselatm; ii++)
    {
      iResCont++;
      sprintf(cTmpString, "%s%d", molecule->rawmol.segId[inp_msdf->sele.selatm[ii]-1], 
                                  molecule->rawmol.resn[inp_msdf->sele.selatm[ii]-1]);
      
      if(strcmp(cTmpString, cLastRes)!=0)
      {
        iResNum++;
        inp_msdf->piSelResLen[iResNum]=iResCont;
        strcpy(cLastRes, cTmpString);
        Res3ToRes1(molecule->rawmol.restype[inp_msdf->sele.selatm[ii]-1], cResCode1);
        sprintf(inp_msdf->pcLabels[iResNum+1], "%s:%s%d", molecule->rawmol.segId[inp_msdf->sele.selatm[ii]-1], cResCode1,
                                                          molecule->rawmol.resn[inp_msdf->sele.selatm[ii]-1]);
                                                        
        inp_msdf->piProgNum[iResNum+1]=molecule->rawmol.presn[inp_msdf->sele.selatm[ii]-1];
        iResCont=0;
      }
    }
    
    iResNum++;
    inp_msdf->piSelResLen[iResNum]=iResCont+1;
    inp_msdf->ppfVirtAtomCoord    = (float **) calloc(inp_msdf->iNumOfSelRes, sizeof(float *));
    for(ii=0; ii<inp_msdf->iNumOfSelRes; ii++)                         // Stores Geometrical centre coordinates
    {
      inp_msdf->ppfVirtAtomCoord[ii]    = (float *) calloc(3, sizeof(float));
    }
  }

  // =================================================================
  nato = inp_msdf->sele.nselatm;
  counter = 0;
  inp_msdf->meanmat = (double **)calloc( nato , sizeof(double *) );
  inp_msdf->meanmat[0] = (double *)calloc( (nato*nato-nato)/2 , sizeof(double ) );
  for( ii=1; ii<nato; ii++ )
  {
    inp_msdf->meanmat[ii] = inp_msdf->meanmat[0] + counter + ii - 1;
    counter += ii-1;
  }
  inp_msdf->d_mat = (double **)calloc( nato , sizeof(double *) );
  inp_msdf->d_mat[0] = (double *)calloc( (nato*nato-nato)/2 , sizeof(double ) );
  for( ii=1; ii<nato; ii++ )
  {
    inp_msdf->d_mat[ii] = inp_msdf->d_mat[0] + counter + ii - 1;
    counter += ii-1;
  }
  inp_msdf->sqmat = (double **)calloc( nato , sizeof(double *) );
  inp_msdf->sqmat[0] = (double *)calloc( (nato*nato-nato)/2 , sizeof(double ) );
  for( ii=1; ii<nato; ii++ )
  {
    inp_msdf->sqmat[ii] = inp_msdf->sqmat[0] + counter + ii - 1;
    counter += ii-1;
  }

  return 0;
}
// ------------------------------------------------------------------
int Compute_Msdf ( struct inp_msdf *inp_msdf, struct sopt *OPT, CoorSet *trj_crd, char *output )
{
  int          ii, jj, nato;
  
  // ===
  int          iResNum         = 0;
  int          iResCont        = 0;
  int          iVirtAtom       = -1;
  
  float        fVirtAtomXCoord = 0.0;
  float        fVirtAtomYCoord = 0.0;
  float        fVirtAtomZCoord = 0.0;
  double       thisfluct;
  
  inp_msdf->iPBCFlag = trj_crd->pbc_flag;
  // ===
  
  if( inp_msdf->pass == 1)
  {
    nato = inp_msdf->sele.nselatm;
    /*
    thisframe = inp_msdf->framecounter;
    
    inp_msdf->distmat[thisframe] = (float **)calloc( nato , sizeof(float *) );
    inp_msdf->distmat[thisframe][0] = (float *)calloc( (nato*nato-nato)/2 , sizeof(float ) );
    for( ii=1; ii<nato; ii++ )
    {
      inp_msdf->distmat[thisframe][ii] = inp_msdf->distmat[thisframe][0] + counter + ii - 1;
      counter += ii-1;
    }
    */
    if(inp_msdf->iMultiAtomFlag == 0)
    {
      for( ii=0; ii<nato; ii++)
        for( jj=0; jj<ii; jj++ )
          //inp_msdf->distmat[thisframe][ii][jj] = sqrt( (trj_crd->xcoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->xcoor[inp_msdf->sele.selatm[jj]-1] ) * (trj_crd->xcoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->xcoor[inp_msdf->sele.selatm[jj]-1] ) +
          //                                              (trj_crd->ycoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->ycoor[inp_msdf->sele.selatm[jj]-1] ) * (trj_crd->ycoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->ycoor[inp_msdf->sele.selatm[jj]-1] ) +
          //                                              (trj_crd->zcoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->zcoor[inp_msdf->sele.selatm[jj]-1] ) * (trj_crd->zcoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->zcoor[inp_msdf->sele.selatm[jj]-1] ) );
          inp_msdf->meanmat[ii][jj] += sqrt( (trj_crd->xcoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->xcoor[inp_msdf->sele.selatm[jj]-1] ) * (trj_crd->xcoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->xcoor[inp_msdf->sele.selatm[jj]-1] ) +
                                                        (trj_crd->ycoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->ycoor[inp_msdf->sele.selatm[jj]-1] ) * (trj_crd->ycoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->ycoor[inp_msdf->sele.selatm[jj]-1] ) +
                                                        (trj_crd->zcoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->zcoor[inp_msdf->sele.selatm[jj]-1] ) * (trj_crd->zcoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->zcoor[inp_msdf->sele.selatm[jj]-1] ) );
    }
    else
    {
      // More than one atom per residue
      
      // === Calculates geo centres of current frame res =================
      iResNum=0;
      iResCont=0;
      iVirtAtom=-1;
      
      fVirtAtomXCoord = 0.0;
      fVirtAtomYCoord = 0.0;
      fVirtAtomZCoord = 0.0;
      
      for(ii=0; ii<inp_msdf->sele.nselatm; ii++)
      {
        if(iResCont == inp_msdf->piSelResLen[iResNum])
        {
          fVirtAtomXCoord = (fVirtAtomXCoord/(float)inp_msdf->piSelResLen[iResNum]);
          fVirtAtomYCoord = (fVirtAtomYCoord/(float)inp_msdf->piSelResLen[iResNum]);
          fVirtAtomZCoord = (fVirtAtomZCoord/(float)inp_msdf->piSelResLen[iResNum]);

          iResNum++;
          iResCont=0;
          iVirtAtom++;
          
          inp_msdf->ppfVirtAtomCoord[iVirtAtom][0]=fVirtAtomXCoord;
          inp_msdf->ppfVirtAtomCoord[iVirtAtom][1]=fVirtAtomYCoord;
          inp_msdf->ppfVirtAtomCoord[iVirtAtom][2]=fVirtAtomZCoord;
          
          fVirtAtomXCoord = 0.0;
          fVirtAtomYCoord = 0.0;
          fVirtAtomZCoord = 0.0;
        }
        fVirtAtomXCoord += trj_crd->xcoor[inp_msdf->sele.selatm[ii]-1];
        fVirtAtomYCoord += trj_crd->ycoor[inp_msdf->sele.selatm[ii]-1];
        fVirtAtomZCoord += trj_crd->zcoor[inp_msdf->sele.selatm[ii]-1];
        
        iResCont++;
      }
      
      // Process the last residue
      fVirtAtomXCoord = (fVirtAtomXCoord/(float)inp_msdf->piSelResLen[iResNum]);
      fVirtAtomYCoord = (fVirtAtomYCoord/(float)inp_msdf->piSelResLen[iResNum]);
      fVirtAtomZCoord = (fVirtAtomZCoord/(float)inp_msdf->piSelResLen[iResNum]);

      iVirtAtom++;
      inp_msdf->ppfVirtAtomCoord[iVirtAtom][0]=fVirtAtomXCoord;
      inp_msdf->ppfVirtAtomCoord[iVirtAtom][1]=fVirtAtomYCoord;
      inp_msdf->ppfVirtAtomCoord[iVirtAtom][2]=fVirtAtomZCoord;

      for(ii=0; ii<inp_msdf->iNumOfSelRes; ii++)
      {
        for( jj=0; jj<ii; jj++ )
        {
          //inp_msdf->distmat[thisframe][ii][jj] = sqrt( (inp_msdf->ppfVirtAtomCoord[ii][0] - inp_msdf->ppfVirtAtomCoord[jj][0] ) * (inp_msdf->ppfVirtAtomCoord[ii][0] - inp_msdf->ppfVirtAtomCoord[jj][0] ) +
          //                                              (inp_msdf->ppfVirtAtomCoord[ii][1] - inp_msdf->ppfVirtAtomCoord[jj][1] ) * (inp_msdf->ppfVirtAtomCoord[ii][1] - inp_msdf->ppfVirtAtomCoord[jj][1] ) +
          //                                              (inp_msdf->ppfVirtAtomCoord[ii][2] - inp_msdf->ppfVirtAtomCoord[jj][2] ) * (inp_msdf->ppfVirtAtomCoord[ii][2] - inp_msdf->ppfVirtAtomCoord[jj][2] ) );
          inp_msdf->meanmat[ii][jj] += sqrt( (inp_msdf->ppfVirtAtomCoord[ii][0] - inp_msdf->ppfVirtAtomCoord[jj][0] ) * (inp_msdf->ppfVirtAtomCoord[ii][0] - inp_msdf->ppfVirtAtomCoord[jj][0] ) +
                                                        (inp_msdf->ppfVirtAtomCoord[ii][1] - inp_msdf->ppfVirtAtomCoord[jj][1] ) * (inp_msdf->ppfVirtAtomCoord[ii][1] - inp_msdf->ppfVirtAtomCoord[jj][1] ) +
                                                        (inp_msdf->ppfVirtAtomCoord[ii][2] - inp_msdf->ppfVirtAtomCoord[jj][2] ) * (inp_msdf->ppfVirtAtomCoord[ii][2] - inp_msdf->ppfVirtAtomCoord[jj][2] ) );
        }
      }
    }

    
    inp_msdf->framecounter ++;
  }
  else if( inp_msdf->pass == 2 )
  {
    //fprintf( stdout, "2 pass is on\n"); fflush(stdout);
    nato = inp_msdf->sele.nselatm;
    if(inp_msdf->iMultiAtomFlag == 0)
    {
      for( ii=1; ii<nato; ii++)
        for( jj=0; jj<ii; jj++ )
        {
          thisfluct = inp_msdf->meanmat[ii][jj] - sqrt( (trj_crd->xcoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->xcoor[inp_msdf->sele.selatm[jj]-1] ) * (trj_crd->xcoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->xcoor[inp_msdf->sele.selatm[jj]-1] ) +
                                                        (trj_crd->ycoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->ycoor[inp_msdf->sele.selatm[jj]-1] ) * (trj_crd->ycoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->ycoor[inp_msdf->sele.selatm[jj]-1] ) +
                                                        (trj_crd->zcoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->zcoor[inp_msdf->sele.selatm[jj]-1] ) * (trj_crd->zcoor[inp_msdf->sele.selatm[ii]-1] - trj_crd->zcoor[inp_msdf->sele.selatm[jj]-1] ) ) ;
          inp_msdf->fluctmat[ii][jj] += thisfluct*thisfluct;
        }
    }
    else
    {
      // More than one atom per residue
      
      // === Calculates geo centres of current frame res =================
      iResNum=0;
      iResCont=0;
      iVirtAtom=-1;
      
      fVirtAtomXCoord = 0.0;
      fVirtAtomYCoord = 0.0;
      fVirtAtomZCoord = 0.0;
      
      for(ii=0; ii<inp_msdf->sele.nselatm; ii++)
      {
        if(iResCont == inp_msdf->piSelResLen[iResNum])
        {
          fVirtAtomXCoord = (fVirtAtomXCoord/(float)inp_msdf->piSelResLen[iResNum]);
          fVirtAtomYCoord = (fVirtAtomYCoord/(float)inp_msdf->piSelResLen[iResNum]);
          fVirtAtomZCoord = (fVirtAtomZCoord/(float)inp_msdf->piSelResLen[iResNum]);

          iResNum++;
          iResCont=0;
          iVirtAtom++;
          
          inp_msdf->ppfVirtAtomCoord[iVirtAtom][0]=fVirtAtomXCoord;
          inp_msdf->ppfVirtAtomCoord[iVirtAtom][1]=fVirtAtomYCoord;
          inp_msdf->ppfVirtAtomCoord[iVirtAtom][2]=fVirtAtomZCoord;
          
          fVirtAtomXCoord = 0.0;
          fVirtAtomYCoord = 0.0;
          fVirtAtomZCoord = 0.0;
        }
        fVirtAtomXCoord += trj_crd->xcoor[inp_msdf->sele.selatm[ii]-1];
        fVirtAtomYCoord += trj_crd->ycoor[inp_msdf->sele.selatm[ii]-1];
        fVirtAtomZCoord += trj_crd->zcoor[inp_msdf->sele.selatm[ii]-1];
        
        iResCont++;
      }
      
      // Process the last residue
      fVirtAtomXCoord = (fVirtAtomXCoord/(float)inp_msdf->piSelResLen[iResNum]);
      fVirtAtomYCoord = (fVirtAtomYCoord/(float)inp_msdf->piSelResLen[iResNum]);
      fVirtAtomZCoord = (fVirtAtomZCoord/(float)inp_msdf->piSelResLen[iResNum]);

      iVirtAtom++;
      inp_msdf->ppfVirtAtomCoord[iVirtAtom][0]=fVirtAtomXCoord;
      inp_msdf->ppfVirtAtomCoord[iVirtAtom][1]=fVirtAtomYCoord;
      inp_msdf->ppfVirtAtomCoord[iVirtAtom][2]=fVirtAtomZCoord;

      for(ii=0; ii<inp_msdf->iNumOfSelRes; ii++)
      {
        for( jj=0; jj<ii; jj++ )
        {
          thisfluct = inp_msdf->meanmat[ii][jj] - sqrt( (inp_msdf->ppfVirtAtomCoord[ii][0] - inp_msdf->ppfVirtAtomCoord[jj][0] ) * (inp_msdf->ppfVirtAtomCoord[ii][0] - inp_msdf->ppfVirtAtomCoord[jj][0] ) +
                                                        (inp_msdf->ppfVirtAtomCoord[ii][1] - inp_msdf->ppfVirtAtomCoord[jj][1] ) * (inp_msdf->ppfVirtAtomCoord[ii][1] - inp_msdf->ppfVirtAtomCoord[jj][1] ) +
                                                        (inp_msdf->ppfVirtAtomCoord[ii][2] - inp_msdf->ppfVirtAtomCoord[jj][2] ) * (inp_msdf->ppfVirtAtomCoord[ii][2] - inp_msdf->ppfVirtAtomCoord[jj][2] ) );
          inp_msdf->fluctmat[ii][jj] += thisfluct*thisfluct;
        }
      }
    }
  }
  else
  {
    fprintf( stderr, "What happened in MSDF?\n");
    fflush(stderr);
    exit(1);
  }
  
  return 0;
}
// ------------------------------------------------------------------
int Post_Msdf ( struct inp_msdf *inp_msdf, struct sopt *OPT, int nframe, Molecule *molecule )
{
  int          ii, jj;
  int          nato, ndist, counter=0;
  FILE        *rmsfdout;
  
  char         cLabelA[20], cLabelB[20], cResCode1[3];
  int          iProgNumA, iProgNumB;
  double       dRMSFd;
  time_t       time_Today;
  
  nato = inp_msdf->iNumOfSelRes;
  ndist = ((nato*nato)-nato)/2;
  
  /*
  meanmat = (double **)calloc( nato , sizeof(double *) );
  meanmat[0] = (double *)calloc( (nato*nato-nato)/2 , sizeof(double ) );
  for( ii=1; ii<nato; ii++ )
  {
    meanmat[ii] = meanmat[0] + counter + ii - 1;
    counter += ii-1;
  }
  counter = 0;
  fluctmat = (double **)calloc( nato , sizeof(double *) );
  fluctmat[0] = (double *)calloc( (nato*nato-nato)/2 , sizeof(double ) );
  for( ii=1; ii<nato; ii++ )
  {
    fluctmat[ii] = fluctmat[0] + counter + ii - 1;
    counter += ii-1;
  }
  
  ndist = ((nato*nato)-nato)/2;
  for( ii=0; ii<nframe; ii++)
   for( jj=0; jj<ndist; jj++ )
     meanmat[0][jj] += (double)inp_msdf->distmat[ii][0][jj];
     //meanmat[0][jj] += (double)inp_msdf->distmat[ii][0][jj]/nframe;
  //for( ii=0; ii<nframe; ii++)
  for( jj=0; jj<ndist; jj++ )
     meanmat[0][jj] /= (double) nframe;
  
  for( ii=0; ii<nframe; ii++)
   for( jj=0; jj<ndist; jj++ )
     fluctmat[0][jj] += ((double)inp_msdf->distmat[ii][0][jj]-meanmat[0][jj])*((double)inp_msdf->distmat[ii][0][jj]-meanmat[0][jj])/nframe;
  
  */ 
  if( inp_msdf->pass == 1 )
  {
    for( jj=0; jj<ndist; jj++ )
       inp_msdf->meanmat[0][jj] /= (double) nframe;
    
    counter = 0;
    inp_msdf->fluctmat = (double **)calloc( nato , sizeof(double *) );
    inp_msdf->fluctmat[0] = (double *)calloc( (nato*nato-nato)/2 , sizeof(double ) );
    for( ii=1; ii<nato; ii++ )
    {
      inp_msdf->fluctmat[ii] = inp_msdf->fluctmat[0] + counter + ii - 1;
      counter += ii-1;
    }
    
    inp_msdf->pass = 2;
    return 1;
  }
  
  for( jj=0; jj<ndist; jj++ )
    inp_msdf->fluctmat[0][jj] /= nframe;
  
  if( inp_msdf->root )
    for( ii=0; ii<ndist; ii++)
      inp_msdf->fluctmat[0][ii] = sqrt( inp_msdf->fluctmat[0][ii] );
  
  rmsfdout = O_File( inp_msdf->title, "w");
  
  if(inp_msdf->verbose == 0)
  {
    fprintf( rmsfdout, "     ");
    for( ii=0; ii<nato; ii++)
      fprintf( rmsfdout, "%8d ", inp_msdf->sele.selatm[ii]);
    fprintf( rmsfdout, "\n");
    for( ii=0; ii<nato; ii++)
    {
      fprintf( rmsfdout, "%4d ", inp_msdf->sele.selatm[ii]);
      for( jj=0; jj<ii; jj++ )
        fprintf( rmsfdout, "%8.3f ", inp_msdf->fluctmat[ii][jj] );
      fprintf( rmsfdout, "\n");
    }
  }
  else
  {
    time(&time_Today);
    fprintf(rmsfdout, "# =================================================================\n");
    fprintf(rmsfdout, "# ***                     WORDOM MSDF MODULE                    ***\n");
    fprintf(rmsfdout, "# =================================================================\n");
    fprintf(rmsfdout, "#\n");
    fprintf(rmsfdout, "# Version   : 0.1a\n");
    fprintf(rmsfdout, "# License   : GPL 3\n");
    fprintf(rmsfdout, "# Copyright : Michele Seeber, Angelo Felline, Francesca Fanelli\n");
    fprintf(rmsfdout, "#             University of Modena and Reggio Emilia\n");
    fprintf(rmsfdout, "#             Modena - Italy\n");
    fprintf(rmsfdout, "#\n");
    fprintf(rmsfdout, "# Date      : %s", asctime(localtime(&time_Today)));
    fprintf(rmsfdout, "#\n");
    fprintf(rmsfdout, "# Mol File  : %s\n", OPT->IMOL_FILE);
    fprintf(rmsfdout, "# Res Num   : %d\n", inp_msdf->iNumOfRes);
    fprintf(rmsfdout, "# Traj File : %s\n", OPT->ITRJ_FILE);
    fprintf(rmsfdout, "# Frame Num : %d\n", inp_msdf->iNumOfFrames);
    
    if(inp_msdf->iPBCFlag == 0)
      fprintf(rmsfdout, "# PBC       : No\n");
    else
      fprintf(rmsfdout, "# PBC       : Yes\n");
    
    fprintf(rmsfdout, "#\n");
    fprintf(rmsfdout, "# Title     : %s\n", inp_msdf->title);
    fprintf(rmsfdout, "# Sele      : %s\n", inp_msdf->sele.selestring);
    
    fprintf(rmsfdout, "# Sele Res  : %d\n", inp_msdf->iNumOfSelRes);
    
    if(inp_msdf->iMultiAtomFlag == 1)
      fprintf(rmsfdout, "# MultiAtom : True\n");
    else
      fprintf(rmsfdout, "# MultiAtom : False\n");
    
    fprintf(rmsfdout, "#\n");
    fprintf(rmsfdout, "# Warning!  : PBC not taken into account\n");
    fprintf(rmsfdout, "# =================================================================\n");
    
    fprintf(rmsfdout, "#%8s   %8s   %15s   %15s   %8s\n", "ProgNumA", "ProgNumB", "LabelA", "LabelB", "RMSFd");
    
    for(ii=0; ii<inp_msdf->iNumOfSelRes; ii++)
    {
      for(jj=0; jj<inp_msdf->iNumOfSelRes; jj++)
      {
        if(inp_msdf->iMultiAtomFlag==0)
        {
          Res3ToRes1(molecule->rawmol.restype[inp_msdf->sele.selatm[ii]-1], cResCode1);
          sprintf(cLabelA, "%s:%s%d", molecule->rawmol.segId[inp_msdf->sele.selatm[ii]-1], cResCode1,
                                      molecule->rawmol.resn[inp_msdf->sele.selatm[ii]-1]);
        
          Res3ToRes1(molecule->rawmol.restype[inp_msdf->sele.selatm[jj]-1], cResCode1);
          sprintf(cLabelB, "%s:%s%d", molecule->rawmol.segId[inp_msdf->sele.selatm[jj]-1], cResCode1,
                                      molecule->rawmol.resn[inp_msdf->sele.selatm[jj]-1]);
          
          iProgNumA=molecule->rawmol.presn[inp_msdf->sele.selatm[ii]-1];
          iProgNumB=molecule->rawmol.presn[inp_msdf->sele.selatm[jj]-1];
        }
        
        else
        {
          strcpy(cLabelA, inp_msdf->pcLabels[ii]);
          strcpy(cLabelB, inp_msdf->pcLabels[jj]);
          iProgNumA=inp_msdf->piProgNum[ii];
          iProgNumB=inp_msdf->piProgNum[jj];
        }
        
        if(ii == jj)
          dRMSFd = 0.0;
        
        else if(ii < jj)
          dRMSFd = inp_msdf->fluctmat[jj][ii];
        
        else if (ii > jj)
          dRMSFd = inp_msdf->fluctmat[ii][jj];

        fprintf(rmsfdout, " %8d   %8d   %15s   %15s   %8.3f\n", iProgNumA, iProgNumB, cLabelA, cLabelB, dRMSFd);
      }
    }
  }
  
  return 0;
}
