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
#include "fileio.h"
#include "tools.h"
#include "datahandler.h"
#include "analysis.h"

int Read_iSample ( char **input, int inp_index, struct inp_sample *inp_sample, char *printout, Molecule *molecule )
{
   int         ii, jj, kk;
   char        buffer[1024];
   char        title[64];
   int         gotit;
   extern short int     frame_par;
   frame_par = 1;
   
   memset ( title, '\0', sizeof(title));
   memset ( buffer, '\0', sizeof(buffer));

   inp_sample->someflag    = 0;
   inp_sample->extout      = 0;
   while( strncmp (buffer, "END", 3))
   {
     gotit = 0;
     sprintf( buffer, "%s", input[inp_index]);
     if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
       gotit = 1;
     if ( !strncmp(buffer, "--SELE", 6))
     {
       sscanf(buffer, "--SELE %[^\n]%*c ", inp_sample->sele.selestring);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--TITLE", 7))
     {
       sscanf( buffer, "--TITLE %s", title);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--SOMEFLAG", 10))
     {
       inp_sample->someflag=1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--EXTOUT", 8))
     {
       inp_sample->extout=1;
       gotit = 1;
     }
     if( gotit==0 )
     {
       fprintf( stderr, "Could not understand option: %s\n", buffer);
       exit(5);
     }
     inp_index++;
   }

   GetSele ( inp_sample->sele.selestring, &inp_sample->sele, molecule);
   
   // get some coordinates from the reference structure if needed
   inp_sample->reference    = (float **)walloc (   inp_sample->sele.nselatm, sizeof (float * ));
   inp_sample->reference[0] = (float *)walloc ( 3*inp_sample->sele.nselatm, sizeof ( float  ));
   for ( ii =0; ii<inp_sample->sele.nselatm; ii++)
    inp_sample->reference[ii] = inp_sample->reference[0] + 3*ii;

   for ( ii=0; ii<inp_sample->sele.nselatm; ii++)
   {
     inp_sample->reference[ii][0] = molecule->coor.xcoor[inp_sample->sele.selatm[ii]-1];
     inp_sample->reference[ii][1] = molecule->coor.ycoor[inp_sample->sele.selatm[ii]-1];
     inp_sample->reference[ii][2] = molecule->coor.zcoor[inp_sample->sele.selatm[ii]-1];
   }

   // prepare storage for frame data, if needed
   inp_sample->moving    = calloc (   inp_sample->sele.nselatm, sizeof (float * ));
   inp_sample->moving[0] = calloc ( 3*inp_sample->sele.nselatm, sizeof ( float  ));
   for ( ii =0; ii<inp_sample->sele.nselatm; ii++)
     inp_sample->moving[ii] = inp_sample->moving[0] + 3*ii;

   // allocate space for anything that is used in every cycle
   // turn on/off flags
   // open additional output files
   // etc
   
   // write header to normal output file
   sprintf( printout, " %12s ", title);
   
   return 12;
}
// ------------------------------------------------------------------
int Compute_Sample ( struct inp_sample *inp_sample, struct sopt *OPT, CoorSet *trj_crd, char *outprint )
{
   int            ii, jj;
   int            avariable;
   float          anothervariable;
   float          morevariables[9];
   float          value=0.0;

   // copy data from frame to local data structure for use
   for ( ii=0; ii<inp_sample->sele.nselatm; ii++)
   {
    inp_sample->moving[ii][0] = trj_crd->xcoor[inp_sample->sele.selatm[ii]-1];
    inp_sample->moving[ii][1] = trj_crd->ycoor[inp_sample->sele.selatm[ii]-1];
    inp_sample->moving[ii][2] = trj_crd->zcoor[inp_sample->sele.selatm[ii]-1];
   }
   
   // run your very own algorithm
   // here, a sum of all coordinates is computed
   for ( ii=0; ii<inp_sample->sele.nselatm; ii++)
   {
    value += inp_sample->moving[ii][0];
    value += inp_sample->moving[ii][1];
    value += inp_sample->moving[ii][2];
   }
   
   sprintf( outprint, " %10.5f ", value);

   return 12;
}
// ------------------------------------------------------------------
int Post_Sample ( struct inp_sample *inp_sample, struct sopt *OPT, int nframe, Molecule *molecule )
{
   // if you need further processing
   // once data has been collected from all frames,
   // this is the place to do it.
   return 0;
}
