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
#include "fileio.h"
#include "analysis.h"
#include "geometry.h"

// ------------------------------------------------------------------
// Distance module
// ------------------------------------------------------------------
int Read_iDist ( char **input, int inp_index, struct inp_dist *inp_dist , struct temp_ll *temp_ll, Molecule *molecule )
{
  // read distance module input
  int           ii;
  char          buffer[256];
  char          title[64];
  int           gotit;
  int           title_offset = 0;
  char         *printout;

  extern short int     frame_par;
  frame_par = 1;
  
  inp_dist->nLines = Get_nSeles ( input, inp_index);
  
  temp_ll->title = calloc( 12*(inp_dist->nLines+1), sizeof(char));
  memset( temp_ll->title, '\0',  12*(inp_dist->nLines+1));
  printout = temp_ll->title;

  inp_dist->seleA = calloc ( inp_dist->nLines, sizeof ( Selection ));
  inp_dist->seleB = calloc ( inp_dist->nLines, sizeof ( Selection ));
  
  inp_dist->mass = 0;
  
  memset ( buffer, '\0', sizeof(buffer));
  ii=0;
  while( strncmp (buffer, "END", 3))
  {
   gotit = 0;
   memset ( buffer, '\0', sizeof(buffer));
   sprintf( buffer, "%s", input[inp_index]);
   if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
    gotit = 1;
   else if ( !strncmp(buffer, "--TITLE", 7))
   {
     memset( title, '\0', 64);
     sscanf( buffer, "--TITLE %s", title);
     sprintf( &printout[title_offset], " %10s ", title);
     title_offset +=12;
     gotit = 1;
   }
   else if ( !strncmp(buffer, "--SELE", 6) )
   {
    sscanf( buffer, "--SELE %[^:]:%[^\n]", inp_dist->seleA[ii].selestring, inp_dist->seleB[ii].selestring);
    GetSele ( inp_dist->seleA[ii].selestring, &inp_dist->seleA[ii], molecule);
    GetSele ( inp_dist->seleB[ii].selestring, &inp_dist->seleB[ii], molecule);
    if( inp_dist->seleA[ii].nselatm > 1 || inp_dist->seleB[ii].nselatm > 1 )
      fprintf( stderr, "Warning: GC doesn't work with PBC yet; centers will be computed w/o PBC, distances with PBC\n");
    
    if( inp_dist->seleA[ii].nselatm == 0 )
    {
     fprintf( stderr, "Check selection: %s - 0 atoms selected\n", inp_dist->seleA[ii].selestring);
     exit(0);
    }
    if( inp_dist->seleB[ii].nselatm == 0 )
    {
     fprintf( stderr, "Check selection: %s - 0 atoms selected\n", inp_dist->seleB[ii].selestring);
     exit(0);
    }
    
    ii++;
    gotit = 1;
   }
   else if ( !strncmp(buffer, "--MASS", 6))
   {
    inp_dist->mass = 1;
    gotit = 1;
   }
   
   if( gotit==0 )
   {
    fprintf( stderr, "Could not understand option: %s\n", buffer);
    exit(5);
   }
   inp_index++;
  }
  
  if( inp_dist->mass )
  {
    inp_dist->masses = (float *) walloc( molecule->nato, sizeof( float) );
    for( ii=0; ii<molecule->nato; ii++ )
      inp_dist->masses[ii] = molecule->rawmol.bFac[ii];
  }
  
  if( inp_dist->nLines != title_offset/12 )
  {
    fprintf( stderr, "Number of --SELE fields differs from # of --TITLE fields!");
    exit(0);
  }
  
  return title_offset;
}
// ------------------------------------------------------------------
int Compute_Dist ( struct inp_dist *inp_dist, struct sopt *OPT, Molecule *molecule, CoorSet *trj_crd, char *outstring )
{
    int         ii, jj;
    float       thismass=0;
    float       xcoor1, ycoor1, zcoor1, xcoor2, ycoor2, zcoor2;
    
    float     **distances = NULL;
    int         counter=0;
    distances = calloc( inp_dist->nLines, sizeof( float *));
    // computing distances
    for ( ii=0; ii<inp_dist->nLines; ii++ )
    {
      if( inp_dist->seleA[ii].nselatm == 1 && inp_dist->seleB[ii].nselatm == 1 )
      {
        distances[ii] = DistanceSelCoor( trj_crd, &inp_dist->seleA[ii], trj_crd, &inp_dist->seleB[ii]);
        continue;
      }
      distances[ii] = calloc( 1, sizeof(float));
      if( inp_dist->seleA[ii].nselatm>1 )
      {
        counter = 0;
        xcoor1 = ycoor1 = zcoor1 = 0;
        if( inp_dist->mass == 0 )
          for( jj=0; jj<inp_dist->seleA[ii].nselatm; jj++ )
          {
            xcoor1 += trj_crd->xcoor[inp_dist->seleA[ii].selatm[jj]-1];
            ycoor1 += trj_crd->ycoor[inp_dist->seleA[ii].selatm[jj]-1];
            zcoor1 += trj_crd->zcoor[inp_dist->seleA[ii].selatm[jj]-1];
            counter ++;
          }
        else if ( inp_dist->mass == 1 )
          for( jj=0; jj<inp_dist->seleA[ii].nselatm; jj++ )
          {
            thismass = inp_dist->masses[inp_dist->seleA[ii].selatm[jj]-1];
            xcoor1 += trj_crd->xcoor[inp_dist->seleA[ii].selatm[jj]-1] * thismass;
            ycoor1 += trj_crd->ycoor[inp_dist->seleA[ii].selatm[jj]-1] * thismass;
            zcoor1 += trj_crd->zcoor[inp_dist->seleA[ii].selatm[jj]-1] * thismass;
            counter += thismass;
          }
        xcoor1 /= counter;
        ycoor1 /= counter;
        zcoor1 /= counter;
      }
      else
      {
        xcoor1 = trj_crd->xcoor[inp_dist->seleA[ii].selatm[0]];
        ycoor1 = trj_crd->ycoor[inp_dist->seleA[ii].selatm[0]];
        zcoor1 = trj_crd->zcoor[inp_dist->seleA[ii].selatm[0]];
      }
      if( inp_dist->seleB[ii].nselatm>1 )
      {
        counter = 0;
        xcoor2 = ycoor2 = zcoor2 = 0;
        if( inp_dist->mass == 0 )
          for( jj=0; jj<inp_dist->seleB[ii].nselatm; jj++ )
          {
            xcoor2 += trj_crd->xcoor[inp_dist->seleB[ii].selatm[jj]-1];
            ycoor2 += trj_crd->ycoor[inp_dist->seleB[ii].selatm[jj]-1];
            zcoor2 += trj_crd->zcoor[inp_dist->seleB[ii].selatm[jj]-1];
            counter ++;
           }
        else if ( inp_dist->mass == 1 )
          for( jj=0; jj<inp_dist->seleB[ii].nselatm; jj++ )
          {
            thismass = inp_dist->masses[inp_dist->seleB[ii].selatm[jj]-1];
            xcoor2 += trj_crd->xcoor[inp_dist->seleB[ii].selatm[jj]-1] * thismass;
            ycoor2 += trj_crd->ycoor[inp_dist->seleB[ii].selatm[jj]-1] * thismass;
            zcoor2 += trj_crd->zcoor[inp_dist->seleB[ii].selatm[jj]-1] * thismass;
            counter += thismass;
          }
        xcoor2 /= counter;
        ycoor2 /= counter;
        zcoor2 /= counter;
      }
      else
      {
        xcoor2 = trj_crd->xcoor[inp_dist->seleB[ii].selatm[0]];
        ycoor2 = trj_crd->ycoor[inp_dist->seleB[ii].selatm[0]];
        zcoor2 = trj_crd->zcoor[inp_dist->seleB[ii].selatm[0]];
      }
      distances[ii][0] = DistanceCoor( xcoor1, ycoor1, zcoor1, xcoor2, ycoor2, zcoor2, trj_crd->pbc);
    }
    for( ii=0; ii<inp_dist->nLines; ii++ )
    {
      sprintf( &outstring[12*ii], " %10.5f ", distances[ii][0]);
      //free(distances[ii]); // DEBUG : this is ridiculous
    }

    //free( distances );  // this causes a segfault on 32-bits architecture. 
    //distances = NULL;   // thus commented - memleak is anyway very small.
    
    return 12*inp_dist->nLines;
    
}
// ------------------------------------------------------------------
// Contact Module
// ------------------------------------------------------------------
int Read_iQ ( char **input, int inp_index, struct inp_Q *inp_Q , char *printout, Molecule *molecule )
{
    // read contacts module input
    int   ii, jj;
    char  buffer[512];
    char  title[64];
    Selection  sele1, sele2;
    int   gotit = 0;
   extern short int     frame_par;
   frame_par = 1;
    
    inp_Q->mass = 0;
    
//  inp_Q->nLines = Get_nSeles ( inpfile);
    inp_Q->nLines = Get_nSeles ( input, inp_index);
    inp_Q->seleA = calloc ( inp_Q->nLines, sizeof ( int *));
    inp_Q->seleB = calloc ( inp_Q->nLines, sizeof ( int *));
    inp_Q->cutoff   = calloc ( inp_Q->nLines, sizeof (float));
    inp_Q->contacts = calloc ( inp_Q->nLines, sizeof ( int ));
    inp_Q->counter  = calloc ( inp_Q->nLines+1, sizeof ( int ));
    for ( ii=0; ii< inp_Q->nLines; ii++)
     inp_Q->counter[ii]=0;
    
    inp_Q->vcts.ndata = inp_Q->nLines;
    FillDistvcts ( &inp_Q->vcts );
    
    memset ( buffer, '\0', sizeof(buffer));
    ii=0;
    while( strncmp (buffer, "END", 3))
    {
      memset ( buffer, '\0', sizeof(buffer));
      sprintf( buffer, "%s", input[inp_index]);
      gotit = 0;
      if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
        gotit = 1;
      else if ( !strncmp (buffer, "--TITLE", 7))
      {
        sscanf( buffer, "--TITLE %s", title);
        sprintf( printout, " %8s ", title);
        inp_Q->counter[0]++;
        gotit = 1;
      }
      else if ( !strncmp(buffer, "--SELE", 6) )
      {
        sscanf( buffer, "--SELE %[^:]:%[^:]:%f", sele1.selestring, sele2.selestring, &inp_Q->cutoff[ii]);
        GetSele ( sele1.selestring, &sele1, molecule);
        GetSele ( sele2.selestring, &sele2, molecule);
        if( sele1.nselatm == 0 || sele2.nselatm == 0 )
        {
          fprintf( stderr, "Error! no atoms selected (%d, %d)(%s - %s)\n", sele1.nselatm, sele2.nselatm, sele1.selestring, sele2.selestring);
          exit(0);
        }
        inp_Q->seleA[ii] = calloc ( sele1.nselatm+1, sizeof (int) );
        inp_Q->seleB[ii] = calloc ( sele2.nselatm+1, sizeof (int) );
        inp_Q->seleA[ii][0] = sele1.nselatm;
        inp_Q->seleB[ii][0] = sele2.nselatm;
        for ( jj=0; jj<sele1.nselatm; jj++)
          inp_Q->seleA[ii][jj+1] = sele1.selatm[jj];
        for ( jj=0; jj<sele2.nselatm; jj++)
          inp_Q->seleB[ii][jj+1] = sele2.selatm[jj];
        inp_Q->counter[inp_Q->counter[0]]++;
        ii++;
        gotit = 1;
      }
      else if ( !strncmp(buffer, "--MASS", 6) )
      {
        inp_Q->mass = 1;
        gotit = 1;
      }
      if( gotit==0 )
      {
        fprintf( stderr, "Could not understand option: %s\n", buffer);
        exit(5);
      }
      inp_index++;
    }
    
    if( inp_Q->mass )
    {
    inp_Q->masses = (float *) walloc( molecule->nato, sizeof( float) );
    for( ii=0; ii<molecule->nato; ii++ )
      inp_Q->masses[ii] = molecule->rawmol.bFac[ii];
    }
    
    return 10*inp_Q->counter[0];
}
// ------------------------------------------------------------------
// ------------------------------------------------------------------
int Compute_Q ( struct inp_Q *inp_Q, struct sopt *OPT, Molecule *molecule, CoorSet *dcd_crd, char *outstring )
{
    int         ii, jj, kk;
    int         nContact=0;
    float       totweight;
    
     // first the coordinates assignements
    if( !inp_Q->mass )
      for( ii=0; ii<inp_Q->nLines; ii++)
      {
        inp_Q->vcts.x0[ii] = 0;
        inp_Q->vcts.y0[ii] = 0;
        inp_Q->vcts.z0[ii] = 0;
        inp_Q->vcts.x1[ii] = 0;
        inp_Q->vcts.y1[ii] = 0;
        inp_Q->vcts.z1[ii] = 0;
        
        for( jj=0; jj<inp_Q->seleA[ii][0]; jj++)
         inp_Q->vcts.x0[ii] += dcd_crd->xcoor[ inp_Q->seleA[ii][jj+1]-1 ];
        inp_Q->vcts.x0[ii] /= inp_Q->seleA[ii][0];
        for( jj=0; jj<inp_Q->seleA[ii][0]; jj++)
         inp_Q->vcts.y0[ii] += dcd_crd->ycoor[ inp_Q->seleA[ii][jj+1]-1 ];
        inp_Q->vcts.y0[ii] /= inp_Q->seleA[ii][0];
        for( jj=0; jj<inp_Q->seleA[ii][0]; jj++)
         inp_Q->vcts.z0[ii] += dcd_crd->zcoor[ inp_Q->seleA[ii][jj+1]-1 ];
        inp_Q->vcts.z0[ii] /= inp_Q->seleA[ii][0];

        for( jj=0; jj<inp_Q->seleB[ii][0]; jj++)
         inp_Q->vcts.x1[ii] += dcd_crd->xcoor[ inp_Q->seleB[ii][jj+1]-1 ];
        inp_Q->vcts.x1[ii] /= inp_Q->seleB[ii][0];
        for( jj=0; jj<inp_Q->seleB[ii][0]; jj++)
         inp_Q->vcts.y1[ii] += dcd_crd->ycoor[ inp_Q->seleB[ii][jj+1]-1 ];
        inp_Q->vcts.y1[ii] /= inp_Q->seleB[ii][0];
        for( jj=0; jj<inp_Q->seleB[ii][0]; jj++)
         inp_Q->vcts.z1[ii] += dcd_crd->zcoor[ inp_Q->seleB[ii][jj+1]-1 ];
        inp_Q->vcts.z1[ii] /= inp_Q->seleB[ii][0];
      }
    else if( inp_Q->mass )
      for( ii=0; ii<inp_Q->nLines; ii++)
      {
        inp_Q->vcts.x0[ii] = 0;
        inp_Q->vcts.y0[ii] = 0;
        inp_Q->vcts.z0[ii] = 0;
        inp_Q->vcts.x1[ii] = 0;
        inp_Q->vcts.y1[ii] = 0;
        inp_Q->vcts.z1[ii] = 0;
        
        totweight = 0;
        for( jj=0; jj<inp_Q->seleA[ii][0]; jj++)
          totweight += inp_Q->masses[inp_Q->seleA[ii][jj+1]-1];
        for( jj=0; jj<inp_Q->seleA[ii][0]; jj++)
          inp_Q->vcts.x0[ii] += dcd_crd->xcoor[ inp_Q->seleA[ii][jj+1]-1 ] * inp_Q->masses[jj];
        inp_Q->vcts.x0[ii] /= totweight;
        for( jj=0; jj<inp_Q->seleA[ii][0]; jj++)
          inp_Q->vcts.y0[ii] += dcd_crd->ycoor[ inp_Q->seleA[ii][jj+1]-1 ] * inp_Q->masses[jj];
        inp_Q->vcts.y0[ii] /= totweight;
        for( jj=0; jj<inp_Q->seleA[ii][0]; jj++)
          inp_Q->vcts.z0[ii] += dcd_crd->zcoor[ inp_Q->seleA[ii][jj+1]-1 ] * inp_Q->masses[jj];
        inp_Q->vcts.z0[ii] /= totweight;

        totweight = 0;
        for( jj=0; jj<inp_Q->seleB[ii][0]; jj++)
          totweight += inp_Q->masses[inp_Q->seleB[ii][jj+1]-1];
        for( jj=0; jj<inp_Q->seleB[ii][0]; jj++)
         inp_Q->vcts.x1[ii] += dcd_crd->xcoor[ inp_Q->seleB[ii][jj+1]-1 ] * inp_Q->masses[jj];
        inp_Q->vcts.x1[ii] /= totweight;
        for( jj=0; jj<inp_Q->seleB[ii][0]; jj++)
         inp_Q->vcts.y1[ii] += dcd_crd->ycoor[ inp_Q->seleB[ii][jj+1]-1 ] * inp_Q->masses[jj];
        inp_Q->vcts.y1[ii] /= totweight;
        for( jj=0; jj<inp_Q->seleB[ii][0]; jj++)
         inp_Q->vcts.z1[ii] += dcd_crd->zcoor[ inp_Q->seleB[ii][jj+1]-1 ] * inp_Q->masses[jj];
        inp_Q->vcts.z1[ii] /= totweight;
      }
    

  // now the real distance calculations - vectorial
    if ( OPT->PBC_FLAG )
     vdistancePBC ( OPT->PBCBOX, &inp_Q->vcts );
    else
     vDistance ( inp_Q->vcts.ndata, inp_Q->vcts.dist, inp_Q->vcts.x0, inp_Q->vcts.y0, inp_Q->vcts.z0, inp_Q->vcts.x1, inp_Q->vcts.y1, inp_Q->vcts.z1);
    
    for ( ii=0; ii<inp_Q->nLines; ii++ )
    {
     if (inp_Q->vcts.dist[ii] < inp_Q->cutoff[ii])
      inp_Q->contacts[ii] = 1;
     else if (inp_Q->vcts.dist[ii] >= inp_Q->cutoff[ii])
      inp_Q->contacts[ii] = 0;
    }
    
    kk=0;
    for ( ii=1; ii<=inp_Q->counter[0]; ii++)
    {
     for ( jj=0; jj<inp_Q->counter[ii]; jj++)
     {
      nContact += inp_Q->contacts[kk];
      kk++;
     }
     sprintf( &outstring[ii*10], " %8d ", nContact);
     nContact = 0;
    }
    if (kk != inp_Q->nLines)
    {
     fprintf( stderr, "Error #96 in Distance module\n\n");
     exit(96);
    }
    
    return 10*inp_Q->counter[0];
}
// ------------------------------------------------------------------
// Contact Module (2) - using DistanceCoor functions
// ------------------------------------------------------------------
int Read_iContacts ( char **input, int inp_index, struct inp_contacts *inp_contacts , struct temp_ll *temp_ll, Molecule *molecule )
{
  // read contacts module input
   int   ii;
   char  buffer[512];
   char  title[64];
   int   gotit = 0;
   int   title_offset=0;
   int   counter = 0;
   char *printout;
   extern short int     frame_par;
   frame_par = 1;
    
   inp_contacts->mass = 0;
    
   inp_contacts->nLines = Get_nSeles ( input, inp_index);
   
   temp_ll->title = calloc( 12*(inp_contacts->nLines+1), sizeof(char));
   memset( temp_ll->title, '\0',  12*(inp_contacts->nLines+1));
   printout = temp_ll->title;

   inp_contacts->seleA = calloc ( inp_contacts->nLines, sizeof ( Selection ));
   inp_contacts->seleB = calloc ( inp_contacts->nLines, sizeof ( Selection ));
   inp_contacts->cutoff   = calloc ( inp_contacts->nLines, sizeof (float));
   inp_contacts->contacts = calloc ( inp_contacts->nLines, sizeof ( int ));
   inp_contacts->counter  = calloc ( inp_contacts->nLines+1, sizeof ( int ));

   memset ( buffer, '\0', sizeof(buffer));
   ii=0;
   while( strncmp (buffer, "END", 3))
   {
     gotit = 0;
     memset ( buffer, '\0', sizeof(buffer));
     sprintf( buffer, "%s", input[inp_index]);
     if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
       gotit = 1;
     else if ( !strncmp(buffer, "--TITLE", 7))
     {
       sscanf( buffer, "--TITLE %s", title);
       sprintf( &printout[title_offset], " %10s ", title);
       title_offset +=12;
       counter ++;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--SELE", 6) )
     {
       sscanf( buffer, "--SELE %[^:]:%[^:]:%f", inp_contacts->seleA[ii].selestring, inp_contacts->seleB[ii].selestring, &inp_contacts->cutoff[ii]);
       GetSele ( inp_contacts->seleA[ii].selestring, &inp_contacts->seleA[ii], molecule);
       GetSele ( inp_contacts->seleB[ii].selestring, &inp_contacts->seleB[ii], molecule);
       if( inp_contacts->seleA[ii].nselatm > 1 || inp_contacts->seleB[ii].nselatm > 1 )
         fprintf( stderr, "Warning: GC doesn't work fully with PBC yet; centers will be computed w/o PBC, distances with PBC\n");
       
       if( inp_contacts->seleA[ii].nselatm == 0 )
       {
        fprintf( stderr, "Check selection: %s - 0 atoms selected\n", inp_contacts->seleA[ii].selestring);
        exit(0);
       }
       if( inp_contacts->seleB[ii].nselatm == 0 )
       {
        fprintf( stderr, "Check selection: %s - 0 atoms selected\n", inp_contacts->seleB[ii].selestring);
        exit(0);
       }
     
       inp_contacts->counter[counter]++;
       ii++;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--MASS", 6))
     {
       inp_contacts->mass = 1;
       gotit = 1;
     }
     
     if( gotit==0 )
     {
       fprintf( stderr, "Could not understand option: %s\n", buffer);
       exit(5);
     }
     inp_index++;
   }
   
   if( inp_contacts->mass )
   {
     inp_contacts->masses = (float *) walloc( molecule->nato, sizeof( float) );
     for( ii=0; ii<molecule->nato; ii++ )
       inp_contacts->masses[ii] = molecule->rawmol.bFac[ii];
   }
   
   inp_contacts->counter[0] = counter;
   
   return title_offset;
}
// ------------------------------------------------------------------
// ------------------------------------------------------------------
int Compute_Contacts ( struct inp_contacts *inp_contacts, struct sopt *OPT, Molecule *molecule, CoorSet *trj_crd, char *outstring )
{
    int         ii, jj, kk;
    int         nContacts=0;
    float       xcoor1, ycoor1, zcoor1, xcoor2, ycoor2, zcoor2;
    float       counter, thismass;
    float     **distances2 = NULL;
    
    distances2 = calloc( inp_contacts->nLines, sizeof( float *));
    // computing distances2
    for ( ii=0; ii<inp_contacts->nLines; ii++ )
    {
      if( inp_contacts->seleA[ii].nselatm == 1 && inp_contacts->seleB[ii].nselatm == 1 )
      {
        distances2[ii] = DistanceSelCoor( trj_crd, &inp_contacts->seleA[ii], trj_crd, &inp_contacts->seleB[ii]);
        continue;
      }
      distances2[ii] = calloc( 1, sizeof(float));
      if( inp_contacts->seleA[ii].nselatm>1 )
      {
        counter = 0;
        xcoor1 = ycoor1 = zcoor1 = 0;
        if( inp_contacts->mass == 0 )
          for( jj=0; jj<inp_contacts->seleA[ii].nselatm; jj++ )
          {
            xcoor1 += trj_crd->xcoor[inp_contacts->seleA[ii].selatm[jj]-1];
            ycoor1 += trj_crd->ycoor[inp_contacts->seleA[ii].selatm[jj]-1];
            zcoor1 += trj_crd->zcoor[inp_contacts->seleA[ii].selatm[jj]-1];
            counter ++;
          }
        else if ( inp_contacts->mass == 1 )
          for( jj=0; jj<inp_contacts->seleA[ii].nselatm; jj++ )
          {
            thismass = inp_contacts->masses[inp_contacts->seleA[ii].selatm[jj-1]];
            xcoor1 += trj_crd->xcoor[inp_contacts->seleA[ii].selatm[jj]-1] * thismass;
            ycoor1 += trj_crd->ycoor[inp_contacts->seleA[ii].selatm[jj]-1] * thismass;
            zcoor1 += trj_crd->zcoor[inp_contacts->seleA[ii].selatm[jj]-1] * thismass;
            counter += thismass;
          }
        xcoor1 /= counter;
        ycoor1 /= counter;
        zcoor1 /= counter;
      }
      else
      {
        xcoor1 = trj_crd->xcoor[inp_contacts->seleA[ii].selatm[0]];
        ycoor1 = trj_crd->ycoor[inp_contacts->seleA[ii].selatm[0]];
        zcoor1 = trj_crd->zcoor[inp_contacts->seleA[ii].selatm[0]];
      }

      if( inp_contacts->seleB[ii].nselatm>1 )
      {
        counter = 0;
        xcoor2 = ycoor2 = zcoor2 = 0;
        if( inp_contacts->mass == 0 )
          for( jj=0; jj<inp_contacts->seleB[ii].nselatm; jj++ )
          {
            xcoor2 += trj_crd->xcoor[inp_contacts->seleB[ii].selatm[jj]-1];
            ycoor2 += trj_crd->ycoor[inp_contacts->seleB[ii].selatm[jj]-1];
            zcoor2 += trj_crd->zcoor[inp_contacts->seleB[ii].selatm[jj]-1];
            counter ++;
          }
        else if ( inp_contacts->mass == 1 )
          for( jj=0; jj<inp_contacts->seleB[ii].nselatm; jj++ )
          {
            thismass = inp_contacts->masses[inp_contacts->seleB[ii].selatm[jj-1]];
            xcoor2 += trj_crd->xcoor[inp_contacts->seleB[ii].selatm[jj]-1] * thismass;
            ycoor2 += trj_crd->ycoor[inp_contacts->seleB[ii].selatm[jj]-1] * thismass;
            zcoor2 += trj_crd->zcoor[inp_contacts->seleB[ii].selatm[jj]-1] * thismass;
            counter += thismass;
          }
        xcoor2 /= counter;
        ycoor2 /= counter;
        zcoor2 /= counter;
      }
      else
      {
        xcoor2 = trj_crd->xcoor[inp_contacts->seleB[ii].selatm[0]];
        ycoor2 = trj_crd->ycoor[inp_contacts->seleB[ii].selatm[0]];
        zcoor2 = trj_crd->zcoor[inp_contacts->seleB[ii].selatm[0]];
      }
      distances2[ii][0] = DistanceCoor( xcoor1, ycoor1, zcoor1, xcoor2, ycoor2, zcoor2, trj_crd->pbc);
    }
    
    for ( ii=0; ii<inp_contacts->nLines; ii++ )
    {
      if ( distances2[ii][0] < inp_contacts->cutoff[ii])
        inp_contacts->contacts[ii] = 1;
      else // if ( distances2[ii][0] >= inp_contacts->cutoff[ii])
        inp_contacts->contacts[ii] = 0;
      free(distances2[ii]);
    }
    kk = 0;
    for( ii=0; ii<inp_contacts->counter[0]; ii++ )
    {
      nContacts = 0;
      for( jj=0; jj<inp_contacts->counter[ii+1]; jj++)
      {
        nContacts += inp_contacts->contacts[kk];
        kk++;
      }
      sprintf( &outstring[12*ii], " %10d ", nContacts);
    }
    free( distances2 );
    distances2 = NULL;

    return 12*inp_contacts->counter[0];
    
}
// ------------------------------------------------------------------
// Angle module
// ------------------------------------------------------------------
int Read_iAngle ( char **input, int inp_index, struct inp_angle *inp_angle, char *printout, Molecule *molecule )
{
   char         buffer[1024];
   char         title[64];
   int          gotit;
   Selection     sele1, sele2, sele3;
   
   extern short int     frame_par;
   frame_par = 1;
   
   memset ( title, '\0', sizeof(title));
   memset ( buffer, '\0', sizeof(buffer));

   inp_angle->someflag    = 0;
   inp_angle->extout      = 0;
   while( strncmp (buffer, "END", 3))
   {
     gotit = 0;
     sprintf( buffer, "%s", input[inp_index]);
     if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
       gotit = 1;
     else if ( !strncmp(buffer, "--SELE", 6))
     {
       sscanf(buffer, "--SELE %[^:]:%[^:]:%[^\n]%*c", sele1.selestring, sele2.selestring, sele3.selestring);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--TITLE", 7))
     {
       sscanf( buffer, "--TITLE %s", title);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--SOMEFLAG", 10))
     {
       inp_angle->someflag=1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--EXTOUT", 8))
     {
       inp_angle->extout=1;
       gotit = 1;
     }
     if( gotit==0 )
     {
       fprintf( stderr, "Could not understand option: %s\n", buffer);
       exit(5);
     }
     inp_index++;
   }

   GetSele ( sele1.selestring, &sele1, molecule);
   GetSele ( sele2.selestring, &sele2, molecule);
   GetSele ( sele3.selestring, &sele3, molecule);
   
   if( sele1.nselatm != 1 || sele2.nselatm != 1 || sele3.nselatm != 1 )
   {
     fprintf( stderr, "Ambiguous selection (%d, %d, %d)\n", sele1.nselatm, sele2.nselatm, sele3.nselatm);
     exit(29);
   }
   
   inp_angle->atm1 = sele1.selatm[0];
   inp_angle->atm2 = sele2.selatm[0];
   inp_angle->atm3 = sele3.selatm[0];

   // write header to normal output file (extout means extended output)
   sprintf( printout, " %10s ", title);
   
   return 12;
}
// ------------------------------------------------------------------
int Compute_Angle ( struct inp_angle *inp_angle, struct sopt *OPT, CoorSet *trj_crd, char *outstring )
{

    double a[3], b[3], c[3], angle;
    
    a[0] = trj_crd->xcoor[inp_angle->atm1 -1];
    a[1] = trj_crd->ycoor[inp_angle->atm1 -1];
    a[2] = trj_crd->zcoor[inp_angle->atm1 -1];
    b[0] = trj_crd->xcoor[inp_angle->atm2 -1];
    b[1] = trj_crd->ycoor[inp_angle->atm2 -1];
    b[2] = trj_crd->zcoor[inp_angle->atm2 -1];
    c[0] = trj_crd->xcoor[inp_angle->atm3 -1];
    c[1] = trj_crd->ycoor[inp_angle->atm3 -1];
    c[2] = trj_crd->zcoor[inp_angle->atm3 -1];

    angle = AngleCalc(a,b,c);

    sprintf( outstring, " %10.5f ", angle);

    return 12;
}
// ------------------------------------------------------------------
// Dihedral module
// ------------------------------------------------------------------
int Read_iDihe ( char **input, int inp_index, struct inp_dihe *inp_dihe, char *printout, Molecule *molecule )
{
  // read dihedral module input
  char          buffer[512];
  char          title[64];
  Selection      sele1, sele2, sele3, sele4;
  int         gotit;
 
   extern short int     frame_par;
   frame_par = 1;
  
  inp_dihe->nLines = Get_nSeles ( input, inp_index);

  memset ( buffer, '\0', sizeof(buffer));
  while( strncmp (buffer, "END", 3))
  {
   gotit = 0;
   memset ( buffer, '\0', sizeof(buffer));
   sprintf( buffer, "%s", input[inp_index]);

   if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
    gotit = 1;
   else if ( !strncmp(buffer, "--TITLE", 7))
   {
    sscanf( buffer, "--TITLE %s", title);
    sprintf( printout, " %10s ", title);
    gotit = 1;
   }
   else if ( !strncmp(buffer, "--SELE", 6) )
   {
    sscanf( buffer, "--SELE %[^:]:%[^:]:%[^:]:%[^\n]", sele1.selestring, sele2.selestring, sele3.selestring, sele4.selestring);
    GetSele ( sele1.selestring, &sele1, molecule);
    GetSele ( sele2.selestring, &sele2, molecule);
    GetSele ( sele3.selestring, &sele3, molecule);
    GetSele ( sele4.selestring, &sele4, molecule);
    if ( sele1.nselatm !=1 || sele2.nselatm !=1 || sele3.nselatm !=1 || sele4.nselatm !=1 )
    {
     if (sele1.nselatm !=1)
      printf("Error in selection: %s\n",sele1.selestring );
     else if (sele2.nselatm !=1)
      printf("Error in selection: %s\n",sele2.selestring );
     else if (sele3.nselatm !=1)
      printf("Error in selection: %s\n",sele3.selestring );
     else if (sele4.nselatm !=1)
      printf("Error in selection: %s\n",sele4.selestring );
     
     exit(85);
    }
    inp_dihe->atm1 = sele1.selatm[0];
    inp_dihe->atm2 = sele2.selatm[0];
    inp_dihe->atm3 = sele3.selatm[0];
    inp_dihe->atm4 = sele4.selatm[0];
    
    gotit = 1;
   }
   if( gotit==0 )
   {
    fprintf( stderr, "Could not understand option: %s\n", buffer);
    exit(5);
   }
   inp_index++;
  }
  
  return 12;
}
// ------------------------------------------------------------------
int Compute_Dihe ( struct inp_dihe *inp_dihe, struct sopt *OPT, Molecule *molecule, CoorSet *dcd_crd, char *outstring )
{
    double a[3];
    double b[3];
    double c[3];
    double d[3];

    double dihe;

    a[0] = (double)dcd_crd->xcoor[ inp_dihe->atm1 -1 ];
    a[1] = (double)dcd_crd->ycoor[ inp_dihe->atm1 -1 ];
    a[2] = (double)dcd_crd->zcoor[ inp_dihe->atm1 -1 ];
    b[0] = (double)dcd_crd->xcoor[ inp_dihe->atm2 -1 ];
    b[1] = (double)dcd_crd->ycoor[ inp_dihe->atm2 -1 ];
    b[2] = (double)dcd_crd->zcoor[ inp_dihe->atm2 -1 ];
    c[0] = (double)dcd_crd->xcoor[ inp_dihe->atm3 -1 ];
    c[1] = (double)dcd_crd->ycoor[ inp_dihe->atm3 -1 ];
    c[2] = (double)dcd_crd->zcoor[ inp_dihe->atm3 -1 ];
    d[0] = (double)dcd_crd->xcoor[ inp_dihe->atm4 -1 ];
    d[1] = (double)dcd_crd->ycoor[ inp_dihe->atm4 -1 ];
    d[2] = (double)dcd_crd->zcoor[ inp_dihe->atm4 -1 ];

    dihe = DiheCalc(a, b, c, d);
    sprintf( outstring, " %10.5f ", dihe);

    return 12;
}

//------------------------------------------------------------------------------
// Radius of Gyration module
//------------------------------------------------------------------------------
int Read_iRgyr ( char **input, int inp_intex, struct inp_rgyr *inp_rgyr, char *printout, Molecule *molecule, CoorSetData *coorsetdata )
{
   int         ii;
   char        buffer[1024];
   char        title[64];
   int         gotit;
   
   extern short int     frame_par;
   frame_par = 1;
   
   memset ( title, '\0', sizeof(title));
   memset ( buffer, '\0', sizeof(buffer));

   inp_rgyr->mass = 0;
   
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
    else if ( !strncmp(buffer, "--SELE", 6))
    {
      sscanf(buffer, "--SELE %[^\n]%*c ", inp_rgyr->sele.selestring);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--MASS", 6))
    {
      inp_rgyr->mass = 1;
      gotit = 1;
    }
    if( gotit==0 )
    {
     fprintf( stderr, "Sorry, could not understand option: %s\n", buffer);
     exit(5);
    }
    inp_intex++;
   }

   GetSele ( inp_rgyr->sele.selestring, &inp_rgyr->sele, molecule);
   
   inp_rgyr->moving    = calloc (   inp_rgyr->sele.nselatm, sizeof (float * ));
   inp_rgyr->moving[0] = calloc ( 3*inp_rgyr->sele.nselatm, sizeof ( float  ));
   for ( ii =0; ii<inp_rgyr->sele.nselatm; ii++)
     inp_rgyr->moving[ii] = inp_rgyr->moving[0] + 3*ii;
   
   if( inp_rgyr->mass )
   {
     inp_rgyr->masses = (float *) walloc( inp_rgyr->sele.nselatm, sizeof(float));
     for( ii=0; ii<inp_rgyr->sele.nselatm; ii++)
       inp_rgyr->masses[ii] = molecule->rawmol.bFac[inp_rgyr->sele.selatm[ii]-1];
   }
   
   sprintf( printout, " %10s ", title);
   
   return 12;
}
// ------------------------------------------------------------------
int Compute_Rgyr ( struct inp_rgyr *inp_rgyr, struct sopt *OPT, CoorSet *dcd_crd, char *outstring )
{
   int            ii;
   float          xcenter=0.0, ycenter=0.0,zcenter=0.0;
   float          xxcenter=0.0, yycenter=0.0,zzcenter=0.0;
   float          rgyr=0.0;
   float          totweight=0.0;
   float          divider;

   if( !inp_rgyr->mass )
   {
     for ( ii=0; ii<inp_rgyr->sele.nselatm; ii++)              //Separate loops privide more cache hits than a single loop and room for vectorization
       xcenter += dcd_crd->xcoor[inp_rgyr->sele.selatm[ii]-1];     
     for ( ii=0; ii<inp_rgyr->sele.nselatm; ii++)
       xxcenter += dcd_crd->xcoor[inp_rgyr->sele.selatm[ii]-1]*dcd_crd->xcoor[inp_rgyr->sele.selatm[ii]-1];     
     for ( ii=0; ii<inp_rgyr->sele.nselatm; ii++)
       ycenter += dcd_crd->ycoor[inp_rgyr->sele.selatm[ii]-1];     
     for ( ii=0; ii<inp_rgyr->sele.nselatm; ii++)
       yycenter += dcd_crd->ycoor[inp_rgyr->sele.selatm[ii]-1]*dcd_crd->ycoor[inp_rgyr->sele.selatm[ii]-1];     
     for ( ii=0; ii<inp_rgyr->sele.nselatm; ii++)
       zcenter += dcd_crd->zcoor[inp_rgyr->sele.selatm[ii]-1];     
     for ( ii=0; ii<inp_rgyr->sele.nselatm; ii++)
       zzcenter += dcd_crd->zcoor[inp_rgyr->sele.selatm[ii]-1]*dcd_crd->zcoor[inp_rgyr->sele.selatm[ii]-1];     
     
     divider = inp_rgyr->sele.nselatm;
     xcenter /= divider;
     xxcenter /= divider;
     ycenter /= divider;
     yycenter /= divider;
     zcenter /= divider;
     zzcenter /= divider;
     rgyr = sqrt( xxcenter - xcenter*xcenter + yycenter - ycenter*ycenter + zzcenter - zcenter*zcenter );
   }
   else if( inp_rgyr->mass )
   {
     for ( ii=0; ii<inp_rgyr->sele.nselatm; ii++)              //Separate loops provide more cache hits than a single loop and room for vectorization
       xcenter += dcd_crd->xcoor[inp_rgyr->sele.selatm[ii]-1] * inp_rgyr->masses[ii];     
     for ( ii=0; ii<inp_rgyr->sele.nselatm; ii++)
       xxcenter += dcd_crd->xcoor[inp_rgyr->sele.selatm[ii]-1] * inp_rgyr->masses[ii] * dcd_crd->xcoor[inp_rgyr->sele.selatm[ii]-1] * inp_rgyr->masses[ii];     
     for ( ii=0; ii<inp_rgyr->sele.nselatm; ii++)
       ycenter += dcd_crd->ycoor[inp_rgyr->sele.selatm[ii]-1] * inp_rgyr->masses[ii];     
     for ( ii=0; ii<inp_rgyr->sele.nselatm; ii++)
       yycenter += dcd_crd->ycoor[inp_rgyr->sele.selatm[ii]-1] * inp_rgyr->masses[ii] * dcd_crd->ycoor[inp_rgyr->sele.selatm[ii]-1] * inp_rgyr->masses[ii];     
     for ( ii=0; ii<inp_rgyr->sele.nselatm; ii++)
       zcenter += dcd_crd->zcoor[inp_rgyr->sele.selatm[ii]-1] * inp_rgyr->masses[ii];     
     for ( ii=0; ii<inp_rgyr->sele.nselatm; ii++)
       zzcenter += dcd_crd->zcoor[inp_rgyr->sele.selatm[ii]-1] * inp_rgyr->masses[ii] * dcd_crd->zcoor[inp_rgyr->sele.selatm[ii]-1] * inp_rgyr->masses[ii];     
     for ( ii=0; ii<inp_rgyr->sele.nselatm; ii++)
       totweight += inp_rgyr->masses[ii];

     xcenter /= totweight;
     ycenter /= totweight;
     zcenter /= totweight;
     
     rgyr = 0.0;
     for( ii=0; ii<inp_rgyr->sele.nselatm; ii++)
     {
       rgyr += ((dcd_crd->xcoor[inp_rgyr->sele.selatm[ii]-1] - xcenter)*(dcd_crd->xcoor[inp_rgyr->sele.selatm[ii]-1] - xcenter) +
                (dcd_crd->ycoor[inp_rgyr->sele.selatm[ii]-1] - ycenter)*(dcd_crd->ycoor[inp_rgyr->sele.selatm[ii]-1] - ycenter) +
                (dcd_crd->zcoor[inp_rgyr->sele.selatm[ii]-1] - zcenter)*(dcd_crd->zcoor[inp_rgyr->sele.selatm[ii]-1] - zcenter)) * inp_rgyr->masses[ii];
     }
     rgyr = sqrt( rgyr/totweight );
   }
   
   sprintf( outstring, " %10.5f ", rgyr);
   
   return 12;
}
// ------------------------------------------------------------------

// ------------------------------------------------------------------
// Order Parameters module
// ------------------------------------------------------------------
int Read_iP ( char **input, int inp_intex, struct inp_P *inp_P, char *printout, Molecule *molecule )
{
    int       ii=0;
    char      buffer[256];
    char      title[64];
    Selection  sele1, sele2;
   int         gotit;
   extern short int     frame_par;
   frame_par = 1;
    
    inp_P->nLines = Get_nSeles ( input, inp_intex);

    inp_P->atm1 = calloc ( inp_P->nLines, sizeof ( int ) );
    inp_P->atm2 = calloc ( inp_P->nLines, sizeof ( int ) );

    memset ( title, '\0', sizeof(title));
    memset ( buffer, '\0', sizeof(buffer));
    ii =0;
    while( strncmp (buffer, "END", 3))
    {
     gotit = 0;
     sprintf( buffer, "%s", input[inp_intex]);
     if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
      gotit = 1;
     if ( !strncmp(buffer, "--SELE", 6) )
     {
      sscanf(buffer, "--SELE %[^:]:%[^\n]", sele1.selestring, sele2.selestring );
      GetSele ( sele1.selestring, &sele1, molecule);
      GetSele ( sele2.selestring, &sele2, molecule);
      if ( sele1.nselatm !=1 || sele2.nselatm != 1 )
      {
       printf("\n\nError: selection not correct in P module call %s - (%s-%s =>%d,%d)\n\n", title, sele1.selestring, sele2.selestring, sele1.nselatm, sele2.nselatm);
       exit(94);
      }
      inp_P->atm1[ii] = sele1.selatm[0];
      inp_P->atm2[ii] = sele2.selatm[0];
      ii++; 
      gotit = 1;
     }
     else if ( !strncmp(buffer, "--TITLE", 7))
     {
      sscanf( buffer, "--TITLE %s", title);
      sprintf( printout, " %8s-1  %8s-2 ", title, title);
      gotit = 1;
     }
     if( gotit==0 )
     {
      fprintf( stderr, "Could not understand option: %s\n", buffer);
      exit(5);
     }
     inp_intex++;
    }

    return 23;
}
// ------------------------------------------------------------------
int Compute_P ( struct inp_P *inp_P, struct sopt *OPT, Molecule *molecule, CoorSet *dcd_crd, char *outstring )
{
    int nVec;
    int ii,jj,kk;

    float xH, yH, zH;
    float xT, yT, zT;

//  float q[3][3];
    float **q;

    float p1 = 0.0;
    float p2 = 0.0;
    float  n[3];
    float  sc1, sc2, scalar;

    float ** vec;
//  float       *result;
    // SSYEV variables
    char JOBZ = 'V';
    char UPLO = 'U';

    int No = 3;
    int LDA = 3;
    int INFO;

    float * W;
    float * A;
    float * WORK;
    float LWORK;


    q = calloc( 3, sizeof(float *));
    for(ii=0; ii<3; ii++)
      q[ii] = calloc( 3, sizeof(float));
//  for(ii=0; ii<3; ii++)
//   for(jj=0; jj<3; jj++)
//    q[ii][jj] = 0 ;
    
    nVec = inp_P->nLines;

    // Orientational vectors memory allocation
    vec = calloc (nVec, sizeof(float *) );
    for(jj=0; jj<nVec; jj++)	
	vec[jj] = calloc (3, sizeof(float) );

    A = calloc (9, sizeof(float *) );

    for ( ii=0; ii<nVec; ii++)
    {
     xH = dcd_crd->xcoor[ inp_P->atm1[ii] -1 ];
     yH = dcd_crd->ycoor[ inp_P->atm1[ii] -1 ];
     zH = dcd_crd->zcoor[ inp_P->atm1[ii] -1 ];
     xT = dcd_crd->xcoor[ inp_P->atm2[ii] -1 ];
     yT = dcd_crd->ycoor[ inp_P->atm2[ii] -1 ];
     zT = dcd_crd->zcoor[ inp_P->atm2[ii] -1 ];
     
     vec[ii][0] = xH - xT;
     vec[ii][1] = yH - yT;
     vec[ii][2] = zH - zT;

     w_norm(vec[ii]);
    }

    for(ii=0; ii<nVec; ii++)
     for(jj=0; jj<3; jj++)     
      for(kk=0; kk<3; kk++)
      {
       q[jj][kk] += 3*vec[ii][jj]*vec[ii][kk];
       if ( jj == kk )
        q[jj][kk] -= 1.0;
      }

    for( ii=0; ii<3; ii++)
     for( jj=0; jj<3; jj++)
      A[ii+jj*3] = q[ii][jj]/( 2.0*nVec);

    LWORK=3*No+1;
    WORK=calloc(LWORK,sizeof(float));
    W=calloc(No,sizeof(float));
    INFO=0;

    #ifdef LAPACK
     ssyev_(&JOBZ, &UPLO, &No, A, &LDA, W, WORK, &LWORK, &INFO);
    #else
     fprintf( stderr, "Sorry, all modules requiring matrix diagonalization with lapack libraries\n have been disabled in this executable\n");
     exit(0);
    #endif
//  ssyev_(&JOBZ, &UPLO, &No, A, &LDA, W, WORK, &LWORK, &INFO);

    n[0] = A[6];
    n[1] = A[7];
    n[2] = A[8];

    sc1=0.0;
    sc2=0.0;
    // Compute Order Param.
    for(ii=0; ii<nVec; ii++)
    {
     scalar = dot_prod(n, vec[ii]);
     sc1 += scalar;
     sc2 += (scalar * scalar);
    }

    p1 = sc1/nVec;
    p2 = 1.5*(sc2/nVec) - 0.5;

    // Free Allocated Memory
    free(A);
    free(W);
    free(WORK);

    for (jj=0; jj<nVec; jj++)
	free( vec[jj] );
    free(vec);

    for(ii=0; ii<3; ii++)
      free(q[ii]);
    free(q);
    
    sprintf( outstring, " %10.5f %10.5f ", p1, p2);
    
    return 23;
}
// ------------------------------------------------------------------
// Hydrogen Bonds module
// ------------------------------------------------------------------
int Read_iHB ( char **input, int inp_intex, struct inp_hb *inp_hb , char *printout, Molecule *molecule )
{
    int       ii;
    int       gotit;
    char      buffer[2048];
    char      title[64];
    Selection  sele1, sele2, sele3;
   
    extern short int     frame_par;
    frame_par = 1;
 
    inp_hb->hbdist = 3.6;
    inp_hb->hbangle = 130.0;
    
    inp_hb->nLines = Get_nSeles ( input, inp_intex);
    
    inp_hb->xcoor1 = (float *)walloc ( inp_hb->nLines, sizeof(float) );
    inp_hb->xcoor2 = (float *)walloc ( inp_hb->nLines, sizeof(float) );
    inp_hb->xcoor3 = (float *)walloc ( inp_hb->nLines, sizeof(float) );
    inp_hb->ycoor1 = (float *)walloc ( inp_hb->nLines, sizeof(float) );
    inp_hb->ycoor2 = (float *)walloc ( inp_hb->nLines, sizeof(float) );
    inp_hb->ycoor3 = (float *)walloc ( inp_hb->nLines, sizeof(float) );
    inp_hb->zcoor1 = (float *)walloc ( inp_hb->nLines, sizeof(float) );
    inp_hb->zcoor2 = (float *)walloc ( inp_hb->nLines, sizeof(float) );
    inp_hb->zcoor3 = (float *)walloc ( inp_hb->nLines, sizeof(float) );
    
    inp_hb->dist   = (float *)walloc ( inp_hb->nLines, sizeof(float) );

    inp_hb->ind0 = calloc ( inp_hb->nLines, sizeof ( int ));
    inp_hb->ind1 = calloc ( inp_hb->nLines, sizeof ( int ));
    inp_hb->ind2 = calloc ( inp_hb->nLines, sizeof ( int ));
    // size of titl (title length) is nlines rather than ntitles because ntitles is not yet known
    inp_hb->titlen = calloc ( inp_hb->nLines, sizeof ( int ));
    for ( ii=0; ii<inp_hb->nLines; ii++)
     inp_hb->titlen[ii] = 0;
    
    memset ( buffer, '\0', sizeof(buffer));
    ii=0;
    while( strncmp (buffer, "END", 3))
    {
      gotit = 0;
      memset ( buffer, '\0', sizeof(buffer));
      sprintf( buffer, "%s", input[inp_intex]);
      if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
        gotit = 1;
      else if ( !strncmp(buffer, "--SELE", 6) )
      {
        sscanf( buffer, "--SELE %[^:]:%[^:]:%[^\n]", sele1.selestring, sele2.selestring, sele3.selestring);

        GetSele ( sele1.selestring, &sele1, molecule);
        GetSele ( sele2.selestring, &sele2, molecule);
        GetSele ( sele3.selestring, &sele3, molecule);
        if ( sele1.nselatm!=1 || sele2.nselatm!=1 || sele3.nselatm!=1 )
        {
          fprintf(stderr, "Ambiguous selection(s)! other than 1 atom selected per field!\n");
          exit(99);
        }
        inp_hb->ind0[ii] = sele1.selatm[0];
        inp_hb->ind1[ii] = sele2.selatm[0];
        inp_hb->ind2[ii] = sele3.selatm[0];

        inp_hb->titlen[inp_hb->ntitles-1] ++;
        ii++;
        gotit = 1;
      }
      else if ( !strncmp(buffer, "--TITLE", 7))
      {
        sscanf( buffer, "--TITLE %s", title);
        sprintf( printout, " %8s ", title);
        inp_hb->ntitles++;
        gotit = 1;
      }
      else if ( !strncmp(buffer, "--ANGLE", 7))
      {
        sscanf( buffer, "--ANGLE %f", &inp_hb->hbangle);
        gotit = 1;
      }
      else if ( !strncmp(buffer, "--DIST", 6))
      {
        sscanf( buffer, "--DIST %f", &inp_hb->hbdist);
        gotit = 1;
      }
      if( gotit==0 )
      {
        fprintf( stderr, "Could not understand option: %s\n", buffer);
        exit(5);
      }
     inp_intex++;
    }
    
    inp_hb->hbrad = (inp_hb->hbangle/180)*3.14159;
    
    return 10*inp_hb->ntitles;
}
// ------------------------------------------------------------------
int Compute_hb ( struct inp_hb *inp_hb, struct sopt *OPT, CoorSet *dcd_crd, char *outstring )
{
   int    ii, jj, kk;
   int    output;
   float  R21[3], R23[3];
   float  angle;
   
   for ( ii=0; ii<inp_hb->nLines; ii++ )
    {
     inp_hb->xcoor1[ii] = dcd_crd->xcoor[inp_hb->ind0[ii]-1];
     inp_hb->ycoor1[ii] = dcd_crd->ycoor[inp_hb->ind0[ii]-1];
     inp_hb->zcoor1[ii] = dcd_crd->zcoor[inp_hb->ind0[ii]-1];
     inp_hb->xcoor2[ii] = dcd_crd->xcoor[inp_hb->ind1[ii]-1];
     inp_hb->ycoor2[ii] = dcd_crd->ycoor[inp_hb->ind1[ii]-1];
     inp_hb->zcoor2[ii] = dcd_crd->zcoor[inp_hb->ind1[ii]-1];
     inp_hb->xcoor3[ii] = dcd_crd->xcoor[inp_hb->ind2[ii]-1];
     inp_hb->ycoor3[ii] = dcd_crd->ycoor[inp_hb->ind2[ii]-1];
     inp_hb->zcoor3[ii] = dcd_crd->zcoor[inp_hb->ind2[ii]-1];
    }
   
   if ( OPT->PBC_FLAG )
    vDistancePBC ( inp_hb->nLines, inp_hb->dist, inp_hb->xcoor1, inp_hb->ycoor1, inp_hb->zcoor1, inp_hb->xcoor3, inp_hb->ycoor3, inp_hb->zcoor3, OPT->PBCBOX );
   else
    vDistance ( inp_hb->nLines, inp_hb->dist, inp_hb->xcoor1, inp_hb->ycoor1, inp_hb->zcoor1, inp_hb->xcoor3, inp_hb->ycoor3, inp_hb->zcoor3);
   
   kk=0;
   for ( ii=0; ii<inp_hb->ntitles; ii++)
   {
     output=0;
     for ( jj=0; jj<inp_hb->titlen[ii]; jj++ )
     {
       if ( inp_hb->dist[kk] < inp_hb->hbdist )
       {
         R21[0] = inp_hb->xcoor2[kk] - inp_hb->xcoor1[kk];
         R21[1] = inp_hb->ycoor2[kk] - inp_hb->ycoor1[kk];
         R21[2] = inp_hb->zcoor2[kk] - inp_hb->zcoor1[kk];
         R23[0] = inp_hb->xcoor2[kk] - inp_hb->xcoor3[kk];
         R23[1] = inp_hb->ycoor2[kk] - inp_hb->ycoor3[kk];
         R23[2] = inp_hb->zcoor2[kk] - inp_hb->zcoor3[kk];
         
         w_norm(R21);
         w_norm(R23);
         
         angle = acos( dot_prod( R21, R23) );
         if ( angle>inp_hb->hbrad )
          output++;
        }
        kk++;
     }

     sprintf( &outstring[10*ii], " %8d ", output );
   }
   
   return 10*inp_hb->ntitles;
}
// ------------------------------------------------------------------
int Init_Within(char **input, int input_index, struct inp_within *inp_within, Molecule *molecule, char *outstring)
{
  int     iWinpOptFlag=0;
  
  char    cWordomInpBuffer[1024], cTmpString[6];

  // === example winp file ===
  //BEGIN
  //--TITLE WITHIN
  //--SELE /*/*/* [n]
  //--LEVEL RES
  //--VERBOSE 1
  //END
  // =========================
  
  // === Default Values ======
  strcpy(inp_within->cTitle, "WITHIN");
  strcpy(inp_within->sele.selestring, "/*/*/*");
  inp_within->iLevelFlag=0;
  inp_within->iVerboseFlag=0;
  // =========================
  
  memset( cWordomInpBuffer, '\0', 1024 );
  while( strncmp (cWordomInpBuffer, "END", 3))
  {
    iWinpOptFlag = 0;
    sprintf( cWordomInpBuffer, "%s", input[input_index]);

    if( !strncmp(cWordomInpBuffer, "BEGIN", 5) || !strncmp(cWordomInpBuffer, "END", 3) || cWordomInpBuffer[0] == '#')
      iWinpOptFlag = 1;

    else if ( !strncmp(cWordomInpBuffer, "--TITLE", 7))
    {
      sscanf( cWordomInpBuffer, "--TITLE %s", inp_within->cTitle);
      iWinpOptFlag = 1;
    }

    else if ( !strncmp(cWordomInpBuffer, "--LEVEL", 7))
    {
      sscanf( cWordomInpBuffer, "--LEVEL %s", cTmpString);
      if(strcmp(cTmpString, "RES")==0)
        inp_within->iLevelFlag=1;
      else if(strcmp(cTmpString, "ATM")==0)
        inp_within->iLevelFlag=0;
      else
      {
        fprintf( stderr, "WITHIN module: Could NOT understand --LEVEL option: %s\n", cTmpString);
        exit(5);
      }
      iWinpOptFlag = 1;
    }

    else if ( !strncmp(cWordomInpBuffer, "--SELE", 6))
    {
      sscanf(cWordomInpBuffer, "--SELE %[^\n]%*c ", inp_within->sele.selestring);
      iWinpOptFlag = 1;
    }
    
    else if ( !strncmp(cWordomInpBuffer, "--VERBOSE", 9))
    {
      sscanf(cWordomInpBuffer, "--VERBOSE %d ", &inp_within->iVerboseFlag);
      
      if(inp_within->iVerboseFlag != 0 && inp_within->iVerboseFlag != 1)
      {
        fprintf( stderr, "WITHIN module: Invalid --VERBOSE value: %s\n", cWordomInpBuffer);
        exit(5);
      }
      
      iWinpOptFlag = 1;
    }
    
    if( iWinpOptFlag==0 )
    {
      fprintf( stderr, "WITHIN module: Could NOT understand option: %s\n", cWordomInpBuffer);
      exit(5);
    }
    input_index++;
  }
  
  if(inp_within->iVerboseFlag == 1)
  {
    sprintf(inp_within->cOutPutFileName, "%s%s", inp_within->cTitle, ".within");
    inp_within->FOutputFile = fopen(inp_within->cOutPutFileName, "w");
    if(inp_within->iLevelFlag == 0)
      fprintf(inp_within->FOutputFile, "#   nFr   NumOfAtm   AtomList\n");
    else
      fprintf(inp_within->FOutputFile, "#   nFr   NumOfRes   ResList\n");
  }
  
  inp_within->TmpMolecule = CopyMol(molecule);
  sprintf( outstring, " %10s ", inp_within->cTitle);
  return 12;
  
}

int Compute_Within(struct inp_within *inp_within, Molecule *molecule, CoorSet *trj_crd, char *outstring, int intex)
{
  int     ii;
  int     iNumOfRes=0;
  int     iResNum=-1, iLastResNum=-1;
  
  char    cListItem[20], cResCode[3];
  
  for(ii=0; ii<molecule->nato; ii++)
  {
    inp_within->TmpMolecule->coor.xcoor[ii] = trj_crd->xcoor[ii];
    inp_within->TmpMolecule->coor.ycoor[ii] = trj_crd->ycoor[ii];
    inp_within->TmpMolecule->coor.zcoor[ii] = trj_crd->zcoor[ii];
  }
  
  GetSele(inp_within->sele.selestring, &inp_within->sele, inp_within->TmpMolecule);
  
  if(inp_within->iLevelFlag == 0)
  {
    // --LEVEL ATM
    sprintf( outstring, " %10d ", inp_within->sele.nselatm);
    
    if(inp_within->iVerboseFlag == 1)
    {
      fprintf(inp_within->FOutputFile, " %6d   %8d   ", intex, inp_within->sele.nselatm);
      
      for(ii=0; ii<inp_within->sele.nselatm; ii++)
      {
        strcpy(cListItem, "");
        Res3ToRes1(molecule->rawmol.restype[inp_within->sele.selatm[ii]-1], cResCode);
        
        sprintf(cListItem, "%s:%s%d:%s", molecule->rawmol.segId[inp_within->sele.selatm[ii]-1],
                                         cResCode,
                                         molecule->rawmol.resn[inp_within->sele.selatm[ii]-1],
                                         molecule->rawmol.atmtype[inp_within->sele.selatm[ii]-1]);
        
        fprintf(inp_within->FOutputFile, "%s ", cListItem);
      }
      
      fprintf(inp_within->FOutputFile, "\n");
    }
  }
  
  else
  {
    // --LEVEL RES
    for(ii=0; ii<inp_within->sele.nselatm; ii++)
    {
      iResNum = molecule->rawmol.resn[inp_within->sele.selatm[ii]-1];
      if(iResNum != iLastResNum)
      {
        iLastResNum = iResNum;
        iNumOfRes++;
      }
    }
    sprintf( outstring, " %10d ", iNumOfRes);
    
    if(inp_within->iVerboseFlag == 1)
    {
      fprintf(inp_within->FOutputFile, " %6d   %8d   ", intex, iNumOfRes);
      
      iLastResNum = -1;
      for(ii=0; ii<inp_within->sele.nselatm; ii++)
      {
        iResNum = molecule->rawmol.resn[inp_within->sele.selatm[ii]-1];
        
        if(iResNum != iLastResNum)
        {
          iLastResNum = iResNum;
          strcpy(cListItem, "");
          Res3ToRes1(molecule->rawmol.restype[inp_within->sele.selatm[ii]-1], cResCode);
          sprintf(cListItem, "%s:%s%d", molecule->rawmol.segId[inp_within->sele.selatm[ii]-1],
                                        cResCode,
                                        molecule->rawmol.resn[inp_within->sele.selatm[ii]-1]);
          
          fprintf(inp_within->FOutputFile, "%s ", cListItem);
        }
      }
      fprintf(inp_within->FOutputFile, "\n");
    }
  }
  
  return 12;
}
