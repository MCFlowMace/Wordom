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
/*! \file samplemod.c
 \brief sample analysis module source file
 
 Sample analysis module source file. When adding a new analysis, use
 this as a template, rename it and include it and call it from
 analysis.c
*/
#include "wordom.h"
#include "tools.h"
#include "pca.h"

// ------------------------------------------------------------------
// PCA module
// ------------------------------------------------------------------
int Read_ipca ( char **input, int inp_index, struct inp_pca *inp_pca, char *printout, Molecule *molecule , int totnframe)
{
   int   ii, jj, kk, ll;
   char  buffer[256];
   char  title[64];
   int   gotit;
   extern short int  no_frame_par;
   no_frame_par = 1;
    
    memset ( title, '\0', sizeof(title));
    memset ( buffer, '\0', sizeof(buffer));
    inp_pca->skip    = 0;
    inp_pca->verbose = 0;
    inp_pca->nprint  = 0;
    
    while( strncmp (buffer, "END", 3))
    {
     gotit = 0;
//   fgets(buffer, 512, inpfile);
     sprintf( buffer, "%s", input[inp_index]);
     if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
      gotit = 1;
     if ( (buffer[0] != '-') && (buffer[0] != '#') && (buffer[0] != '\n') && (buffer[0] != 'B') && (buffer[0] != 'E') )
     {
      fprintf( stderr, "Missing parameter flag in: \n%s maybe --SELE?\n\n", buffer);
      exit(0);
     }
     else if ( !strncmp(buffer, "--TITLE", 7))
     {
      sscanf( buffer, "--TITLE %s", title);
      gotit = 1;
     }
     else if ( !strncmp(buffer, "--SELE", 6))
     {
      sscanf( buffer, "--SELE %[^\n]%*c", inp_pca->sele.selestring );
      gotit = 1;
     }
     else if ( !strncmp(buffer, "--NPRINT", 8))
     {
      sscanf( buffer, "--NPRINT %d", &inp_pca->nprint );
      gotit = 1;
     }
     else if ( !strncmp(buffer, "--VERBOSE", 9))
     {
      inp_pca->verbose = 1;
      gotit = 1;
     }
     else if ( !strncmp(buffer, "--PROGRESSIVE", 13))
     {
      sscanf( buffer, "--PROGRESSIVE %d", &inp_pca->skip );
      gotit = 1;
      if( inp_pca->skip == 0)
      {
       fprintf( stderr, "Wrong value for skip: check --PROGRESSIVE flag\n");
       exit(0);
      }
     }
     if( gotit==0 )
     {
      fprintf( stderr, "Could not understand option: %s\n", buffer);
      exit(5);
     }
     inp_index++;
    }
   
    sscanf(title, "%s", inp_pca->title);
    
    GetSele ( inp_pca->sele.selestring, &inp_pca->sele, molecule);
    if( inp_pca->nprint > (inp_pca->sele.nselatm*3) )
    {
      fprintf( stderr, "You asked for %d printed eigenvectors but only selected\n \
               %d atom (n_vectors = 3*n_atoms))", inp_pca->nprint, inp_pca->sele.nselatm);
      exit(0);
    }
    
    inp_pca->xavgcoor = (float *)walloc( inp_pca->sele.nselatm, sizeof(float));
    inp_pca->yavgcoor = (float *)walloc( inp_pca->sele.nselatm, sizeof(float));
    inp_pca->zavgcoor = (float *)walloc( inp_pca->sele.nselatm, sizeof(float));
    for( ii=0; ii<inp_pca->sele.nselatm; ii++ )
     inp_pca->xavgcoor[ii] = inp_pca->yavgcoor[ii] = inp_pca->zavgcoor[ii] = 0;
   
   ll=0;
   for(ii=0; ii<molecule->nSeg; ii++)
   {
    for(jj=0; jj<molecule->segment[ii].nRpS; jj++)
    {
     for(kk=0; kk<molecule->segment[ii].pRes[jj].nApR; kk++)
     {
      if(inp_pca->sele.selatm[ll] == molecule->segment[ii].pRes[jj].pAto[kk].atomn)
      {
       inp_pca->xavgcoor[ll] = molecule->segment[ii].pRes[jj].pAto[kk].xCoor;
       inp_pca->yavgcoor[ll] = molecule->segment[ii].pRes[jj].pAto[kk].yCoor;
       inp_pca->zavgcoor[ll] = molecule->segment[ii].pRes[jj].pAto[kk].zCoor;
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
    
    inp_pca->covmatrix = (float **)walloc( inp_pca->sele.nselatm*3, sizeof(float *));
    for (ii=0; ii<inp_pca->sele.nselatm*3; ii++)
    {
     inp_pca->covmatrix[ii] = (float *)walloc( inp_pca->sele.nselatm*3, sizeof(float));
     for ( jj=0; jj<inp_pca->sele.nselatm*3; jj++)
      inp_pca->covmatrix[ii][jj] = 0;
    }
    
    if( inp_pca->skip )
    {
     inp_pca->nsteps = totnframe/inp_pca->skip;
     inp_pca->covmatrices = calloc( inp_pca->nsteps+1, sizeof(float **));
     for( ii=0; ii<inp_pca->nsteps+1; ii++ )
     {
      inp_pca->covmatrices[ii] = calloc( inp_pca->sele.nselatm*3, sizeof(float *));
      for( jj=0; jj<inp_pca->sele.nselatm*3; jj++ )
      {
       inp_pca->covmatrices[ii][jj] = calloc( inp_pca->sele.nselatm*3, sizeof(float *));
       for( kk=0; kk<inp_pca->sele.nselatm*3; kk++ )
        inp_pca->covmatrices[ii][jj][kk] = 0;
      }
     }
    }
    inp_pca->cords = calloc( 3*inp_pca->sele.nselatm, sizeof(float));
    inp_pca->dx = calloc( inp_pca->sele.nselatm, sizeof(float));
    inp_pca->dy = calloc( inp_pca->sele.nselatm, sizeof(float));
    inp_pca->dz = calloc( inp_pca->sele.nselatm, sizeof(float));
    
    return 0;
}
// ------------------------------------------------------------------
int Compute_PCA ( struct inp_pca *inp_pca, struct sopt *OPT, CoorSet *trj_crd, char *outprint , int nframe)
{
   int            ii, jj, kk, ll, index;
   float         *xcoor, *ycoor, *zcoor;
   float         *dx, *dy, *dz, *distances;
   
   xcoor = &inp_pca->cords[0];
   ycoor = &inp_pca->cords[inp_pca->sele.nselatm];
   zcoor = &inp_pca->cords[2*inp_pca->sele.nselatm];
   for( ii=0; ii<inp_pca->sele.nselatm; ii++)
   {
     xcoor[ii] = trj_crd->xcoor[ (inp_pca->sele.selatm[ii])-1 ];
     ycoor[ii] = trj_crd->ycoor[ (inp_pca->sele.selatm[ii])-1 ];
     zcoor[ii] = trj_crd->zcoor[ (inp_pca->sele.selatm[ii])-1 ];
   }
   dx = inp_pca->dx;
   dy = inp_pca->dy;
   dz = inp_pca->dz;
   
   if( trj_crd->pbc == NULL )
   {
     for( kk=0; kk<inp_pca->sele.nselatm; kk++ )
     {
       dx[kk] = xcoor[kk] - inp_pca->xavgcoor[kk];
       dy[kk] = ycoor[kk] - inp_pca->yavgcoor[kk];
       dz[kk] = zcoor[kk] - inp_pca->zavgcoor[kk];
     }
   }
   else
   {
     distances = calloc( 3, sizeof(float));
     for( kk=0; kk<inp_pca->sele.nselatm; kk++ )
     {
       DistanceAxis( distances, xcoor[kk], ycoor[kk], zcoor[kk], inp_pca->xavgcoor[kk], inp_pca->yavgcoor[kk], inp_pca->zavgcoor[kk], trj_crd->pbc);
       dx[kk] = distances[0];
       dy[kk] = distances[1];
       dz[kk] = distances[2];
     }
   }
   
   for ( kk=0; kk<inp_pca->sele.nselatm*3; kk+=3 )
   {
     index = kk/3;
     for ( ll=0; ll<inp_pca->sele.nselatm*3; ll+=3 )
       inp_pca->covmatrix[kk][ll] += dx[index]*dx[ll/3];
     for ( ll=1; ll<inp_pca->sele.nselatm*3; ll+=3 )
       inp_pca->covmatrix[kk][ll] += dx[index]*dy[(ll-1)/3];
     for ( ll=2; ll<inp_pca->sele.nselatm*3; ll+=3 )
       inp_pca->covmatrix[kk][ll] += dx[index]*dz[(ll-2)/3];
   }
   for ( kk=1; kk<inp_pca->sele.nselatm*3; kk+=3 )
   {  
     index = (kk-1)/3;
     for ( ll=0; ll<inp_pca->sele.nselatm*3; ll+=3 )
       inp_pca->covmatrix[kk][ll] += dy[index]*dx[ll/3];
     for ( ll=1; ll<inp_pca->sele.nselatm*3; ll+=3 )
       inp_pca->covmatrix[kk][ll] += dy[index]*dy[(ll-1)/3];
     for ( ll=2; ll<inp_pca->sele.nselatm*3; ll+=3 )
       inp_pca->covmatrix[kk][ll] += dy[index]*dz[(ll-2)/3];
   }
   for ( kk=2; kk<inp_pca->sele.nselatm*3; kk+=3 )
   {  
     index = (kk-2)/3;
     for ( ll=0; ll<inp_pca->sele.nselatm*3; ll+=3 )
       inp_pca->covmatrix[kk][ll] += dz[index]*dx[ll/3];
     for ( ll=1; ll<inp_pca->sele.nselatm*3; ll+=3 )
       inp_pca->covmatrix[kk][ll] += dz[index]*dy[(ll-1)/3];
     for ( ll=2; ll<inp_pca->sele.nselatm*3; ll+=3 )
       inp_pca->covmatrix[kk][ll] += dz[index]*dz[(ll-2)/3];
   }

   if( inp_pca->skip )
   {
    if( !fmod(nframe, inp_pca->skip) )
    {
     index = (nframe/inp_pca->skip);
     for( ii=0; ii<inp_pca->sele.nselatm*3; ii++ )
      for( jj=0; jj<inp_pca->sele.nselatm*3; jj++ )
       inp_pca->covmatrices[index][ii][jj] = inp_pca->covmatrix[ii][jj];
    }
   }
   
   return 0;
}
// ------------------------------------------------------------------
int Post_PCA ( struct inp_pca *inp_pca, struct sopt *OPT, int nframe, Molecule *molecule )
{
   FILE  *output;
   int    ii, jj, kk, ll;
   char   filename[128];

   _eigen  eigen;
   float bFac_jj;
   
   _eigen  *eigens;
   
   float **projection;
   
   float fTotVar, fThisVar, fProgVar; // used to compute the fraction of variance explained by eigenvectors
   
   if( inp_pca->skip )
   {
    
    eigens = (_eigen *) calloc( inp_pca->nsteps+1, sizeof(_eigen));
    
    projection = (float **) calloc( inp_pca->nprint*inp_pca->nprint, sizeof( float *));
    for( ii=0; ii<inp_pca->nprint; ii++ )
    {
     projection[ii] = (float *) calloc( inp_pca->nprint, sizeof(float));
     for( jj=0; jj<inp_pca->nprint; jj++ )
      projection[ii][jj] = 0.0;
    }
    
    for( ll=1; ll<inp_pca->nsteps+1; ll++ )
    {
     for( ii=0;  ii<inp_pca->sele.nselatm*3; ii++)
      for( jj=0; jj<inp_pca->sele.nselatm*3; jj++)
       inp_pca->covmatrices[ll][ii][jj] /= inp_pca->skip*ll;
     
     eigens[ll].size = inp_pca->sele.nselatm*3;
     eigens[ll].inpmat = inp_pca->covmatrices[ll];
     DiagMatrix ( &eigens[ll] );
     
     if( inp_pca->verbose && (inp_pca->skip*ll)>0 );
     {
      sprintf( filename, "%s-eigval_%d.txt", inp_pca->title, inp_pca->skip*ll);
      output = O_File ( filename, "w");
      for( ii=0; ii<inp_pca->sele.nselatm*3; ii++)
       fprintf(output, "%d\t%7.6f\n", ii+1, eigens[ll].eigval[ii]);
   
      fclose(output);
      memset ( filename, '\0', sizeof(filename));
   
      // write out eigenvectors; aligned in columns
      sprintf( filename, "%s-eigvec_%d.txt", inp_pca->title, inp_pca->skip*ll);
      output = O_File ( filename, "w");
      for( ii=0; ii<inp_pca->sele.nselatm*3; ii++)
      {
       for( jj=0; jj<inp_pca->sele.nselatm*3; jj++)
        //  aligned in columns
        fprintf(output, "%9.5f", eigens[ll].eigvec[jj][ii]);
       fprintf(output, "\n");
      }
      fclose(output);
     }
    }
    
    for( ll=2; ll<inp_pca->nsteps+1; ll++ )
    {
     printf("comparing PCA at %d with PCA at %d\n", inp_pca->skip*ll, inp_pca->skip*(ll-1));
     for( ii=0; ii<inp_pca->nprint; ii++ )
     {
      for( jj=0; jj<inp_pca->nprint; jj++ )
      {
       for( kk=0; kk<inp_pca->sele.nselatm*3; kk++ )
        projection[ii][jj] += eigens[ll].eigvec[ii][kk]*eigens[ll-1].eigvec[jj][kk];
       printf("%10.5f ", projection[ii][jj]);
      }
      printf("\n");
     }
     
     for( ii=0; ii<inp_pca->nprint; ii++ )
      for( jj=0; jj<inp_pca->nprint; jj++ )
       projection[ii][jj] = 0.0;
     
    }
    
   }
   
   sprintf( filename, "%s-matrix.txt", inp_pca->title);
   output = O_File ( filename, "w");
   
   for( ii=0; ii<inp_pca->sele.nselatm*3; ii++)
    {
     for( jj=0; jj<inp_pca->sele.nselatm*3; jj++)
      {
       inp_pca->covmatrix[ii][jj] /= nframe;
       fprintf(output, "%10.5f ", inp_pca->covmatrix[ii][jj]);
      }
     fprintf( output, "\n");
    }
   
   fclose(output);
   memset ( filename, '\0', sizeof(filename));
   
   if( !inp_pca->skip)
   {
    eigen.size   = inp_pca->sele.nselatm*3;
    eigen.inpmat = inp_pca->covmatrix;
    DiagMatrix( &eigen);
   }
   else if (inp_pca->skip)
   {
    eigen.eigval = eigens[inp_pca->nsteps].eigval;
    eigen.eigvec = eigens[inp_pca->nsteps].eigvec;
   }
   
   sprintf( filename, "%s-eigval.txt", inp_pca->title);
   output = O_File ( filename, "w");
   for( ii=0; ii<inp_pca->sele.nselatm*3; ii++)
    fprintf(output, "%d\t%7.6f\n", ii+1, eigen.eigval[ii]);
   
   fclose(output);
   memset ( filename, '\0', sizeof(filename));
   
   // === *-* ANF *-* ==================================================
   sprintf( filename, "%s-fractvar.txt", inp_pca->title);
   output = O_File ( filename, "w");
   fprintf(output, "%5s   %8s   %8s\n", "EVect", "Var%", "ProgVar%");
   
   fTotVar  = 0.0;
   fThisVar = 0.0;
   fProgVar = 0.0;
   for( ii=0; ii<inp_pca->sele.nselatm*3; ii++)
     fTotVar += fabsf(eigen.eigval[ii]);
     
   for( ii=0; ii<inp_pca->sele.nselatm*3; ii++)
   {
     fThisVar = ((eigen.eigval[ii]/fTotVar)*100.0);
     fProgVar += fThisVar;
     fprintf(output, "%5d   %8.3f   %8.3f\n", ii+1, fThisVar, fProgVar);
   }
   fclose(output);
   
   memset ( filename, '\0', sizeof(filename));
   // === *-* ANF *-* ==================================================
   
   // write out eigenvctors
   sprintf( filename, "%s-eigvec.txt", inp_pca->title);
   output = O_File ( filename, "w");
   for( ii=0; ii<inp_pca->sele.nselatm*3; ii++)
   {
    for( jj=0; jj<inp_pca->sele.nselatm*3; jj++)
     fprintf(output, "%9.5f", eigen.eigvec[jj][ii]);
    fprintf(output, "\n");
   }
   fclose(output);
   
   for( ii=0; ii<inp_pca->nprint; ii++ )
   {
    for( jj=0; jj< molecule->nSeg; jj++)
     for( kk=0; kk< molecule->segment[jj].nRpS; kk++)
      for( ll=0; ll< molecule->segment[jj].pRes[kk].nApR; ll++)
       molecule->segment[jj].pRes[kk].pAto[ll].bFac = 0.0;
    
    for( jj=0; jj< molecule->nato; jj++)
     molecule->rawmol.bFac[jj] = 0.0;
     
    for( jj=0; jj<inp_pca->sele.nselatm; jj++)
    {
     // compute the module of the projection of the ii-th PC on jj-th atom coordinate system
     bFac_jj = sqrt( (eigen.eigvec[ii][jj*3+0]) *(eigen.eigvec[ii][jj*3+0])  + 
                     (eigen.eigvec[ii][jj*3+1]) *(eigen.eigvec[ii][jj*3+1])  +
                     (eigen.eigvec[ii][jj*3+2]) *(eigen.eigvec[ii][jj*3+2])  );
     // write the eigenvector components as beta-factors
     molecule->rawmol.bFac[inp_pca->sele.selatm[jj]-1] = 100 * bFac_jj ;//fabs(100*(A[ii+jj*(inp_pca->sele.nselatm*3)]));
    }
    
    sprintf( filename, "%s-eig%d.pdb", inp_pca->title, ii+1);
    WritePdb_unstr( filename, molecule );
   }
   
   free(eigen.eigvec);
   free(eigen.eigval);
   
   return 0;
}
// ------------------------------------------------------------------
// Projection module (accessory to PCA)
// ------------------------------------------------------------------
int Read_ipro ( char **input, int inp_index, struct inp_pro *inp_pro, char *printout, Molecule *molecule )
{
   FILE  *matrix;
   int    ii, jj, index;
   char   buffer[256];
   char   title[64];
   char   filename[64];
   float *matbuffer;
   // some data struct for dcd out writing
   int         gotit;
   extern short int  no_frame_par;
   no_frame_par = 1;
   
   memset ( title, '\0', sizeof(title));
   memset ( buffer, '\0', sizeof(buffer));
   inp_pro->dcdflag = 0;
   inp_pro->rangeflag = 0;
   inp_pro->spacflag = 0;
   while( strncmp (buffer, "END", 3))
   {
    gotit = 0;
//  fgets(buffer, 512, inpfile);
    sprintf( buffer, "%s", input[inp_index]);
    if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
     gotit = 1;
    if ( !strncmp(buffer, "--TITLE", 7))
    {
     sscanf( buffer, "--TITLE %s", inp_pro->title);
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--FILE", 6) )
    {
     sscanf(buffer, "--FILE %[^\n]%*c", filename);
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--VECTOR", 8) )
    {
     sscanf(buffer, "--VECTOR %d%*[^\n]%*c", &inp_pro->nmode);
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--SELE", 6))
    {
     sscanf(buffer, "--SELE %[^:]:%[^\n]%*c ", inp_pro->sele1.selestring, inp_pro->sele2.selestring);
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--DCDOUT", 8))
    {
     inp_pro->dcdflag = 1;
     sscanf( buffer, "--DCDOUT %[^\n]%*c ", inp_pro->dcdout);
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--RANGEFILL", 11))
    {
     inp_pro->spacflag = 1;
     sscanf( buffer, "--RANGEFILL %f", &inp_pro->spacing);
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--RANGE", 7))
    {
     inp_pro->rangeflag = 1;
     sscanf( buffer, "--RANGE %f", &inp_pro->range);
     gotit = 1;
    }
    if( gotit==0 )
    {
     fprintf( stderr, "Could not understand option: %s\n", buffer);
     exit(5);
    }
    inp_index++;
   }
   
   if( inp_pro->title[0] == '\0' || filename[0] == '\0' || inp_pro->nmode == 0 || inp_pro->sele1.selestring[0] == '\0' || inp_pro->sele2.selestring[0] == '\0' )
   {
    printf("projection module missing some parameters: check \n");
    exit(1);
   }

   GetSele ( inp_pro->sele1.selestring, &inp_pro->sele1, molecule);
   GetSele ( inp_pro->sele2.selestring, &inp_pro->sele2, molecule);
   
   inp_pro->sel1sel2 = (int *)walloc ( inp_pro->sele2.nselatm, sizeof(int));
   for( ii=0; ii<inp_pro->sele2.nselatm; ii++ )
    inp_pro->sel1sel2[ii] = 0;
   index = 0;
   for( ii=0; ii<inp_pro->sele2.nselatm; ii++)
   {
    for( jj=index; jj<inp_pro->sele1.nselatm; jj++)
    {
     if (inp_pro->sele2.selatm[ii] == inp_pro->sele1.selatm[jj])
     {
      inp_pro->sel1sel2[ii] = jj;
      index++;
      break;
     }
    }
   }
   
   inp_pro->xtmpcoor = calloc( molecule->nato, sizeof(float));
   inp_pro->ytmpcoor = calloc( molecule->nato, sizeof(float));
   inp_pro->ztmpcoor = calloc( molecule->nato, sizeof(float));
   inp_pro->xrefcoor = calloc( molecule->nato, sizeof(float));
   inp_pro->yrefcoor = calloc( molecule->nato, sizeof(float));
   inp_pro->zrefcoor = calloc( molecule->nato, sizeof(float));
   
   for( ii=0; ii<molecule->nato; ii++)
   {
    inp_pro->xrefcoor[ii] = molecule->coor.xcoor[ii];
    inp_pro->yrefcoor[ii] = molecule->coor.ycoor[ii];
    inp_pro->zrefcoor[ii] = molecule->coor.zcoor[ii];
   }
   
   if( inp_pro->dcdflag )
   {
    inp_pro->projected_frames = 0;
    if( inp_pro->sele1.nselatm != inp_pro->sele2.nselatm )
    {
     fprintf( stderr, "Dcd out implemented only with sele1 = sele2 up to now");
     exit(0);
    }
    inp_pro->dcdoutfile = O_File( inp_pro->dcdout, "w" );
    FillDcdHeader( molecule, &inp_pro->dcdouthdr );
    WriteDcdHeader( inp_pro->dcdoutfile, &inp_pro->dcdouthdr, "all");
    CalloCoor( &inp_pro->coor, molecule->nato );
   }
   
   matbuffer = calloc ( inp_pro->sele1.nselatm*3, sizeof(float) );
   inp_pro->modevec = (float *)walloc ( inp_pro->sele1.nselatm*3, sizeof(float));
   for( ii=0; ii<inp_pro->sele1.nselatm*3; ii++ )
     inp_pro->modevec[ii] = 0.0;
   
   matrix = O_File ( filename, "r");

   for ( ii=0; ii<inp_pro->sele1.nselatm*3; ii++)
   {
    for ( jj=0; jj<inp_pro->nmode+1; jj++)
    {
     fscanf(matrix, "%f", &matbuffer[jj]);
    }
    inp_pro->modevec[ii] = matbuffer[inp_pro->nmode-1];
    fscanf( matrix, "%*[^\n]*c");
   }

   sprintf( printout, " %10s ", inp_pro->title);
   return 12;
}   
// ------------------------------------------------------------------
int Compute_pro ( struct inp_pro *inp_pro, struct sopt *OPT, CoorSet *dcd_crd, char *outprint )
{
   int   ii;
   float projection = 0;
   
   for( ii=0; ii<inp_pro->sele2.nselatm; ii++)
    {
     inp_pro->xtmpcoor[ii] = dcd_crd->xcoor[ (inp_pro->sele2.selatm[ii])-1 ];
     inp_pro->ytmpcoor[ii] = dcd_crd->ycoor[ (inp_pro->sele2.selatm[ii])-1 ];
     inp_pro->ztmpcoor[ii] = dcd_crd->zcoor[ (inp_pro->sele2.selatm[ii])-1 ];
    }
   for( ii=0; ii<inp_pro->sele2.nselatm; ii++)
    {
     inp_pro->xtmpcoor[ii] -= inp_pro->xrefcoor[(inp_pro->sele2.selatm[ii])-1];
     inp_pro->xtmpcoor[ii] *= inp_pro->modevec[(inp_pro->sel1sel2[ii]*3)+0];
     inp_pro->ytmpcoor[ii] -= inp_pro->yrefcoor[(inp_pro->sele2.selatm[ii])-1];
     inp_pro->ytmpcoor[ii] *= inp_pro->modevec[(inp_pro->sel1sel2[ii]*3)+1];
     inp_pro->ztmpcoor[ii] -= inp_pro->zrefcoor[(inp_pro->sele2.selatm[ii])-1];
     inp_pro->ztmpcoor[ii] *= inp_pro->modevec[(inp_pro->sel1sel2[ii]*3)+2];
     inp_pro->xtmpcoor[ii] += inp_pro->ytmpcoor[ii];
     inp_pro->xtmpcoor[ii] += inp_pro->ztmpcoor[ii];
    }
   
   for( ii=0; ii<inp_pro->sele2.nselatm; ii++)
    projection += inp_pro->xtmpcoor[ii];
   
   if( projection<inp_pro->minpro )
    inp_pro->minpro = projection;
   if( projection>inp_pro->maxpro )
    inp_pro->maxpro = projection;
   
   sprintf( outprint, " %10.5f ", projection ) ;
   
   // now use xcoor, ycoor and zcoor to store coordinates to be written in the outdcd
   if( inp_pro->dcdflag )
   {
    for( ii=0; ii<inp_pro->dcdouthdr.nato; ii++)
    {
     inp_pro->coor.xcoor[ii] = inp_pro->xrefcoor[ii];
     inp_pro->coor.ycoor[ii] = inp_pro->yrefcoor[ii];
     inp_pro->coor.zcoor[ii] = inp_pro->zrefcoor[ii];
    }
    
    for( ii=0; ii<inp_pro->sele2.nselatm; ii++ )
    {
     inp_pro->coor.xcoor[(inp_pro->sele2.selatm[ii])-1] = projection * inp_pro->modevec[(inp_pro->sel1sel2[ii]*3)+0];
     inp_pro->coor.xcoor[(inp_pro->sele2.selatm[ii])-1] += inp_pro->xrefcoor[(inp_pro->sele2.selatm[ii])-1];
     inp_pro->coor.ycoor[(inp_pro->sele2.selatm[ii])-1] = projection * inp_pro->modevec[(inp_pro->sel1sel2[ii]*3)+1];
     inp_pro->coor.ycoor[(inp_pro->sele2.selatm[ii])-1] += inp_pro->yrefcoor[(inp_pro->sele2.selatm[ii])-1];
     inp_pro->coor.zcoor[(inp_pro->sele2.selatm[ii])-1] = projection * inp_pro->modevec[(inp_pro->sel1sel2[ii]*3)+2];
     inp_pro->coor.zcoor[(inp_pro->sele2.selatm[ii])-1] += inp_pro->zrefcoor[(inp_pro->sele2.selatm[ii])-1];
    }
    WriteDcdCoor( inp_pro->dcdoutfile, &inp_pro->coor, &inp_pro->dcdouthdr);
    inp_pro->projected_frames ++;
   }

   return 12;
}
// ------------------------------------------------------------------
int Post_pro ( struct inp_pro *inp_pro, struct sopt *OPT, int nframe, Molecule *molecule )
{
   int          ii, jj;
   float        projection = 0, step = 0;
   char         title[256];
   
   // if called, update the nframe in the projections trajectory
   if( inp_pro->dcdflag )
   {
    inp_pro->dcdouthdr.nframe = inp_pro->projected_frames;
    WriteDcdHeader( inp_pro->dcdoutfile, &inp_pro->dcdouthdr, "nframe" );
    fclose(inp_pro->dcdoutfile);
   }
   
   // if called, prepare traj with range min-max projection for nice paper pictures
   if( inp_pro->spacflag )
   {
    sprintf(title, "%s_rangefill.dcd", inp_pro->title);
    inp_pro->dcdoutfile = O_File( title, "w" );
    FillDcdHeader( molecule, &inp_pro->dcdouthdr );
    WriteDcdHeader( inp_pro->dcdoutfile, &inp_pro->dcdouthdr, "all");
    
    if( !inp_pro->dcdflag )
    {
     inp_pro->coor.xcoor = calloc( molecule->nato, sizeof(float));
     inp_pro->coor.ycoor = calloc( molecule->nato, sizeof(float));
     inp_pro->coor.zcoor = calloc( molecule->nato, sizeof(float));
    }
    
    jj=0;
    
    if( inp_pro->rangeflag )
    {
      inp_pro->minpro = -1*fabs(inp_pro->range);
      inp_pro->maxpro =    fabs(inp_pro->range);
    }
    if( inp_pro->spacing<0 )
      step = -(inp_pro->maxpro-inp_pro->minpro)/inp_pro->spacing;
    else
      step = inp_pro->spacing;
    for( projection = inp_pro->minpro; projection<=inp_pro->maxpro; projection+=step)
    {
     for( ii=0; ii<inp_pro->dcdouthdr.nato; ii++)
     {
      inp_pro->coor.xcoor[ii] = inp_pro->xrefcoor[ii];
      inp_pro->coor.ycoor[ii] = inp_pro->yrefcoor[ii];
      inp_pro->coor.zcoor[ii] = inp_pro->zrefcoor[ii];
     }
     
     for( ii=0; ii<inp_pro->sele2.nselatm; ii++ )
     {
      inp_pro->coor.xcoor[(inp_pro->sele2.selatm[ii])-1] = projection * inp_pro->modevec[(inp_pro->sel1sel2[ii]*3)+0];
      inp_pro->coor.xcoor[(inp_pro->sele2.selatm[ii])-1] += inp_pro->xrefcoor[(inp_pro->sele2.selatm[ii])-1];
      inp_pro->coor.ycoor[(inp_pro->sele2.selatm[ii])-1] = projection * inp_pro->modevec[(inp_pro->sel1sel2[ii]*3)+1];
      inp_pro->coor.ycoor[(inp_pro->sele2.selatm[ii])-1] += inp_pro->yrefcoor[(inp_pro->sele2.selatm[ii])-1];
      inp_pro->coor.zcoor[(inp_pro->sele2.selatm[ii])-1] = projection * inp_pro->modevec[(inp_pro->sel1sel2[ii]*3)+2];
      inp_pro->coor.zcoor[(inp_pro->sele2.selatm[ii])-1] += inp_pro->zrefcoor[(inp_pro->sele2.selatm[ii])-1];
     }
     WriteDcdCoor( inp_pro->dcdoutfile, &inp_pro->coor, &inp_pro->dcdouthdr);
     jj++;
    }
    
    inp_pro->dcdouthdr.nframe = jj;
    WriteDcdHeader( inp_pro->dcdoutfile, &inp_pro->dcdouthdr, "nframe" );
    fclose(inp_pro->dcdoutfile);
   }
   

   return 0;
}
