// ------------------------------------------------------------------
// Copyright (C) 2009  University of Modena
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
/*! \file ssa.c
 \brief secondary structure assignment source file
 
 Secondary Structure Assignment module source code
*/
#include "wordom.h"
#include "fileio.h"
#include "tools.h"
#include "datahandler.h"
#include "analysis.h"
#include "ssa.h"

// ------------------------------------------------------------------
//------------------------------------------------------------------------------
// Secondary Structure Assignment module
//------------------------------------------------------------------------------
int Read_iSSA ( char **input, int inp_intex, struct inp_ss *inp_ss, char *printout, Molecule *molecule, int nframe )
{
   int           ii;
   char          buffer[1024];
   int           ntotres=0;
   int           gotit;
   char         *lastsel;
   int           p_counter;
   
   extern short int     frame_par;
   frame_par = 1;
   
   memset ( buffer, '\0', sizeof(buffer));

   inp_ss->dssp = inp_ss->cont = 0;
   sprintf( inp_ss->sele.selestring, "/*/*/*");
   while( strncmp (buffer, "END", 3))
   {
     gotit = 0;
     memset ( buffer, '\0', sizeof(buffer));
     sprintf( buffer, "%s", input[inp_intex]);
     if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
       gotit = 1;
     else if ( !strncmp(buffer, "--TITLE", 7))
     {
      sscanf( buffer, "--TITLE %s", inp_ss->title);
      sprintf( printout, " %8s ", inp_ss->title);
      gotit = 1;
     }
     else if ( !strncmp(buffer, "--DCLIKE", 8))
     {
      inp_ss->cont = 1;
      gotit = 1;
     }
     else if ( !strncmp(buffer, "--DLIKE", 7))
     {
      inp_ss->dssp = 1;
      gotit = 1;
     }
     else if ( !strncmp(buffer, "--SELE", 6))
     {
      sscanf(buffer, "--SELE %[^\n]%*c ", inp_ss->sele.selestring);
      gotit = 1;
      fprintf( stderr, "selection not implemented yet for ssa\n");
      exit(0);
     }
     if( gotit==0 )
     {
      fprintf( stderr, "Could not understand option : \"%s\"\n", buffer);
      exit(5);
     }
     inp_intex++;
   }
   
   if( (inp_ss->dssp + inp_ss->cont) != 1 )
   {
    fprintf( stderr, "Error! you must choose one (and only one) method for ss prediction\n");
    exit(47);
   }
   
   lastsel = strrchr( inp_ss->sele.selestring, '/');
   if( lastsel[1] != '*' )
   {
     fprintf(stderr, "you can't really select %s for ssa\n", lastsel);
     exit(0);
   }
   
   GetSele( inp_ss->sele.selestring, &inp_ss->sele, molecule);
   
   ntotres = 0;
   
   if( !strcmp( inp_ss->sele.selestring, "/*/*/*") )
     for( ii=0; ii<molecule->nSeg; ii++ )
       ntotres += molecule->segment[ii].nRpS;
   else
     for( ii=0; ii<molecule->nSeg; ii++ )
       ntotres += molecule->segment[ii].nRpS;
   
   inp_ss->ntotres = ntotres;
   
   inp_ss->clist = calloc( ntotres, sizeof(int));
   inp_ss->olist = calloc( ntotres, sizeof(int));
   inp_ss->nlist = calloc( ntotres, sizeof(int));
   inp_ss->hlist = calloc( ntotres, sizeof(int));
   inp_ss->calist = calloc( ntotres, sizeof(int));
   inp_ss->intE  = calloc( ntotres, sizeof(float *));
   inp_ss->intE[0] = calloc( ntotres*ntotres, sizeof( float));
   for( ii=0; ii<ntotres; ii++ )
    inp_ss->intE[ii] = inp_ss->intE[0] + ii*ntotres;
   
   inp_ss->hbonds  = calloc( ntotres, sizeof(short *));
   inp_ss->hbonds[0] = calloc( ntotres*ntotres, sizeof( short));
   for( ii=0; ii<ntotres; ii++ )
    inp_ss->hbonds[ii] = inp_ss->hbonds[0] + ii*ntotres;
   
   inp_ss->xcoor1 = (float *) walloc( ntotres*ntotres, sizeof( float));
   inp_ss->xcoor2 = (float *) walloc( ntotres*ntotres, sizeof( float));
   inp_ss->xcoor3 = (float *) walloc( ntotres*ntotres, sizeof( float));
   inp_ss->xcoor4 = (float *) walloc( ntotres*ntotres, sizeof( float));
   inp_ss->xcoor5 = (float *) walloc( ntotres*ntotres, sizeof( float));
   inp_ss->ycoor1 = (float *) walloc( ntotres*ntotres, sizeof( float));
   inp_ss->ycoor2 = (float *) walloc( ntotres*ntotres, sizeof( float));
   inp_ss->ycoor3 = (float *) walloc( ntotres*ntotres, sizeof( float));
   inp_ss->ycoor4 = (float *) walloc( ntotres*ntotres, sizeof( float));
   inp_ss->ycoor5 = (float *) walloc( ntotres*ntotres, sizeof( float));
   inp_ss->zcoor1 = (float *) walloc( ntotres*ntotres, sizeof( float));
   inp_ss->zcoor2 = (float *) walloc( ntotres*ntotres, sizeof( float));
   inp_ss->zcoor3 = (float *) walloc( ntotres*ntotres, sizeof( float));
   inp_ss->zcoor4 = (float *) walloc( ntotres*ntotres, sizeof( float));
   inp_ss->zcoor5 = (float *) walloc( ntotres*ntotres, sizeof( float));
   
   inp_ss->distON = (float *) walloc( ntotres*ntotres, sizeof( float));
   inp_ss->distCH = (float *) walloc( ntotres*ntotres, sizeof( float));
   inp_ss->distOH = (float *) walloc( ntotres*ntotres, sizeof( float));
   inp_ss->distCN = (float *) walloc( ntotres*ntotres, sizeof( float));
   
   inp_ss->secstruct = (char *) walloc( ntotres, sizeof(char));
   inp_ss->ssnotes   = (char *) walloc( ntotres, sizeof(char));
   
   p_counter = 0;
   for( ii=0; ii<molecule->nato; ii++ )
   {
     if( !strcmp(molecule->rawmol.atmtype[ii], "N") )
       inp_ss->nlist[p_counter] = molecule->rawmol.atomn[ii];//-1;
     else if( !strcmp(molecule->rawmol.atmtype[ii], "HN") ||  !strcmp(molecule->rawmol.atmtype[ii], "H"))
       inp_ss->hlist[p_counter] = molecule->rawmol.atomn[ii];//-1;
     else if( !strcmp(molecule->rawmol.atmtype[ii], "C") )
       inp_ss->clist[p_counter] = molecule->rawmol.atomn[ii];//-1;
     else if( !strcmp(molecule->rawmol.atmtype[ii], "O") || !strcmp(molecule->rawmol.atmtype[ii], "OT1"))
       inp_ss->olist[p_counter] = molecule->rawmol.atomn[ii];//-1;
     else if( !strcmp(molecule->rawmol.atmtype[ii], "CA") || !strcmp(molecule->rawmol.atmtype[ii], "OT1"))
       inp_ss->calist[p_counter] = molecule->rawmol.atomn[ii];//-1;
     if( (molecule->rawmol.resn[ii] != molecule->rawmol.resn[ii+1]) || strncmp(molecule->rawmol.restype[ii], molecule->rawmol.restype[ii+1], 3) || strncmp(molecule->rawmol.segId[ii], molecule->rawmol.segId[ii+1], 4))
       p_counter++;
   }
   
   if( p_counter!=ntotres )
   {
    fprintf( stderr, "Uh-oh! error in residues number\n");
    exit(57);
   }
   
   if( inp_ss->cont )
   {
    inp_ss->cutoffs[0] = -0.2  ;
    inp_ss->cutoffs[1] = -0.25 ;
    inp_ss->cutoffs[2] = -0.3  ;
    inp_ss->cutoffs[3] = -0.4  ;
    inp_ss->cutoffs[4] = -0.5  ;
    inp_ss->cutoffs[5] = -0.6  ;
    inp_ss->cutoffs[6] = -0.7  ;
    inp_ss->cutoffs[7] = -0.8  ;
    inp_ss->cutoffs[8] = -0.9  ;
    inp_ss->cutoffs[9] = -1    ;
   
    inp_ss->weights[0] = 4  ;
    inp_ss->weights[1] = 7  ;
    inp_ss->weights[2] = 23 ;
    inp_ss->weights[3] = 23 ;
    inp_ss->weights[4] = 22 ;
    inp_ss->weights[5] = 14 ;
    inp_ss->weights[6] = 7  ;
    inp_ss->weights[7] = 4  ;
    inp_ss->weights[8] = 2  ;
    inp_ss->weights[9] = 1  ;
    
    inp_ss->secstructS = (char **) walloc( 10, sizeof(char *));
    inp_ss->secstructS[0] = (char *) walloc( 10*ntotres, sizeof(char));
    for( ii=0; ii<10; ii++ )
     inp_ss->secstructS[ii] = inp_ss->secstructS[0] + ii*ntotres;
    inp_ss->ssnotesS   = (char **) walloc( ntotres, sizeof(char *));
    inp_ss->ssnotesS[0] = (char *) walloc( 10*ntotres, sizeof(char));
    for( ii=0; ii<10; ii++ )
     inp_ss->ssnotesS[ii] = inp_ss->ssnotesS[0] + ii*ntotres;
   }
   
   return (inp_ss->ntotres + 2);
}
//------------------------------------------------------------------------------
int Compute_SSA ( struct inp_ss *inp_ss, struct sopt *OPT, Molecule *molecule, CoorSet *trj_crd, char * output )
{
   int          ii, jj, ci;
   float        cutoff;
   
   int          coeff[8];       // coefficient for each sstype (used for cont)
   int          max, maxss;
   
   // this computes the secondary structure following the approach described in:
   // Kabsch W. and Sander C. ; Biopolymers, 1983, 22, 2577-2637
   
   // first compute distances CO_ii-HN_jj for each ii-jj pair
   computessDistances( inp_ss, trj_crd );
   
   // computation of interaction energy (intE) for each putative hbond
   mkintEmatrix( inp_ss );
   
   for( ii=0; ii<inp_ss->ntotres; ii++ )
    inp_ss->secstruct[ii] = 'X';
    
   if( inp_ss->dssp == 1 )
   {
    cutoff = -0.5;
    
     // get hbonds matrix from interaction matrix by comparing 
     // with cutoff ( if intE[ii]<cutoff ->hbond[ii] = 1 );
    for( ii=0; ii < (inp_ss->ntotres*inp_ss->ntotres); ii++ )
     inp_ss->hbonds[0][ii] = ( inp_ss->intE[0][ii] < cutoff ) ? 1 : 0 ;
    
   // we are NOT interested in H-bonds with adiacent residues - let's set them to 0
    for( ii=1; ii<inp_ss->ntotres-1; ii++ )
     inp_ss->hbonds[ii][ii-1] = inp_ss->hbonds[ii][ii+1] = 0;
    inp_ss->hbonds[0][1] = inp_ss->hbonds[inp_ss->ntotres-1][inp_ss->ntotres-2] = 0;
    
    sStruct( inp_ss, inp_ss->secstruct );
   }
   else if( inp_ss->cont == 1 )
   {
    for( ii=0; ii<10*inp_ss->ntotres; ii++ )
     inp_ss->secstructS[0][ii] = 'X';
    // SS PRED WITH CONT-like PROCEDURE
    
    for( ci=0; ci<10; ci++ )
    {
     cutoff = inp_ss->cutoffs[ci];
     
     for( ii=0; ii < (inp_ss->ntotres*inp_ss->ntotres); ii++ )
      inp_ss->hbonds[0][ii] = ( inp_ss->intE[0][ii] < cutoff ) ? 1 : 0 ;
    
     for( ii=1; ii<inp_ss->ntotres-1; ii++ )
      inp_ss->hbonds[ii][ii-1] = inp_ss->hbonds[ii][ii+1] = 0;
     inp_ss->hbonds[0][1] = inp_ss->hbonds[inp_ss->ntotres-1][inp_ss->ntotres-2] = 0;
    
     sStruct( inp_ss, inp_ss->secstructS[ci] );
    }
    
    for( ii=0; ii<inp_ss->ntotres; ii++ )
    {
     for( jj=0; jj<8; jj++ )
      coeff[jj] = 0;
     for( ci=0; ci<10; ci++ )
     {
      switch(inp_ss->secstructS[ci][ii])
      {
       case 'X':
        coeff[0] += inp_ss->weights[ci];
        break;
       case 'H':
        coeff[1] += inp_ss->weights[ci];
        break;
       case 'B':
        coeff[2] += inp_ss->weights[ci];
        break;
       case 'E':
        coeff[3] += inp_ss->weights[ci];
        coeff[2] += inp_ss->weights[ci];        // B and E together for consensus purposes
       break;
       case 'G':
        coeff[4] += inp_ss->weights[ci];
        break;
       case 'I':
        coeff[5] += inp_ss->weights[ci];
        break;
       case 'T':
        coeff[6] += inp_ss->weights[ci];
        break;
       case 'S':
        coeff[7] += inp_ss->weights[ci];
        break;
       case 'L':
        coeff[0] += inp_ss->weights[ci];        // X and L together for consensus purposes
        break;
       default:
        abort ();
      }
     }

     // find the bigger coeff (L and X together, B and E together)
     max = 0;
     maxss = 0;
     for( jj=0; jj<8; jj++ )
     {
      if( coeff[jj] > max )
      {
       max = coeff[jj];
       maxss = jj;
      }
     }
     switch (maxss)
     {
      case 0:
       inp_ss->secstruct[ii] = 'L';
       break;
      case 1:
       inp_ss->secstruct[ii] = 'H';
       break;
      case 2:
       inp_ss->secstruct[ii] = 'B';
       break;
      case 3:
       inp_ss->secstruct[ii] = 'E';
       break;
      case 4:
       inp_ss->secstruct[ii] = 'G';
       break;
      case 5:
       inp_ss->secstruct[ii] = 'I';
       break;
      case 6:
       inp_ss->secstruct[ii] = 'T';
       break;
      case 7:
       inp_ss->secstruct[ii] = 'S';
       break;
     }
     
     if( ii>0 && (inp_ss->secstruct[ii]=='B' || inp_ss->secstruct[ii]=='E') && (inp_ss->secstruct[ii-1]=='B' || inp_ss->secstruct[ii-1]=='E'))
      inp_ss->secstruct[ii-1] = inp_ss->secstruct[ii] = 'E';
    }
   }

   sprintf( output, " %s ", inp_ss->secstruct);
   
   return (inp_ss->ntotres + 2);
}
// ------------------------------------------------------------------
int computessDistances( struct inp_ss *inp_ss, CoorSet *dcd_crd)
{
   int  ii, jj, kk;
   int  cc, oo, nn, hh, ca;
   int  npairs;
   
   // compute distances CO_ii-HN_jj for each ii-jj pair
   kk = 0;
   for( ii=0; ii<inp_ss->ntotres; ii++ )
    for( jj=0; jj<inp_ss->ntotres; jj++ )
    {
     cc = inp_ss->clist[ii]-1;
     oo = inp_ss->olist[ii]-1;
     nn = inp_ss->nlist[jj]-1;
     hh = inp_ss->hlist[jj]-1;
     ca = inp_ss->calist[jj]-1;
     
     if( cc == -1 )
     {
      inp_ss->xcoor1[kk] = 9999.0;
      inp_ss->ycoor1[kk] = 9999.0;
      inp_ss->zcoor1[kk] = 9999.0;
     }
     else
     {
      inp_ss->xcoor1[kk] = dcd_crd->xcoor[cc];
      inp_ss->ycoor1[kk] = dcd_crd->ycoor[cc];
      inp_ss->zcoor1[kk] = dcd_crd->zcoor[cc];
     }
     
     if( oo == -1 )
     {
      inp_ss->xcoor2[kk] = 9999.0;
      inp_ss->ycoor2[kk] = 9999.0;
      inp_ss->zcoor2[kk] = 9999.0;
     }
     else
     {
      inp_ss->xcoor2[kk] = dcd_crd->xcoor[oo];
      inp_ss->ycoor2[kk] = dcd_crd->ycoor[oo];
      inp_ss->zcoor2[kk] = dcd_crd->zcoor[oo];
     }
     
     if( nn == -1 )
     {
      inp_ss->xcoor3[kk] = 9999.0;
      inp_ss->ycoor3[kk] = 9999.0;
      inp_ss->zcoor3[kk] = 9999.0;
     }
     else
     {
      inp_ss->xcoor3[kk] = dcd_crd->xcoor[nn];
      inp_ss->ycoor3[kk] = dcd_crd->ycoor[nn];
      inp_ss->zcoor3[kk] = dcd_crd->zcoor[nn];
     }
     
     if( hh == -1 )
     {
      inp_ss->xcoor4[kk] = 9999.0;
      inp_ss->ycoor4[kk] = 9999.0;
      inp_ss->zcoor4[kk] = 9999.0;
     }
     else
     {
      inp_ss->xcoor4[kk] = dcd_crd->xcoor[hh];
      inp_ss->ycoor4[kk] = dcd_crd->ycoor[hh];
      inp_ss->zcoor4[kk] = dcd_crd->zcoor[hh];
     }
     
     if( ca == -1 )
     {
      inp_ss->xcoor5[kk] = 9999.0;
      inp_ss->ycoor5[kk] = 9999.0;
      inp_ss->zcoor5[kk] = 9999.0;
     }
     else
     {
      inp_ss->xcoor5[kk] = dcd_crd->xcoor[ca];
      inp_ss->ycoor5[kk] = dcd_crd->ycoor[ca];
      inp_ss->zcoor5[kk] = dcd_crd->zcoor[ca];
     }
     kk ++;
    }
   npairs = kk;
   vDistance ( npairs, inp_ss->distON, inp_ss->xcoor2, inp_ss->ycoor2, inp_ss->zcoor2, inp_ss->xcoor3, inp_ss->ycoor3, inp_ss->zcoor3);
   vDistance ( npairs, inp_ss->distCH, inp_ss->xcoor1, inp_ss->ycoor1, inp_ss->zcoor1, inp_ss->xcoor4, inp_ss->ycoor4, inp_ss->zcoor4);
   vDistance ( npairs, inp_ss->distOH, inp_ss->xcoor2, inp_ss->ycoor2, inp_ss->zcoor2, inp_ss->xcoor4, inp_ss->ycoor4, inp_ss->zcoor4);
   vDistance ( npairs, inp_ss->distCN, inp_ss->xcoor1, inp_ss->ycoor1, inp_ss->zcoor1, inp_ss->xcoor3, inp_ss->ycoor3, inp_ss->zcoor3);
   
   return npairs;
}
// ------------------------------------------------------------------
int mkintEmatrix( struct inp_ss *inp_ss )
{
   // computation of interaction energy (intE) for each putative hbond
   // stems from E = q_1*q_2(1/r(ON)+1/r(CH)-1/r(OH)-1/r(CN))*f
   // with q_1 = 0.42e , q_2 = 0.20e and f = 332
   // as reported in Kabsch and Sander
   int  ii, jj, kk;
   kk = 0;
   for( ii=0; ii<inp_ss->ntotres; ii++ )
    for( jj=0; jj<inp_ss->ntotres; jj++ )
    {
     if ( ii!=jj )
      inp_ss->intE[ii][jj] = 27.888 * ( (1.0/inp_ss->distON[kk]) + (1.0/inp_ss->distCH[kk]) - (1.0/inp_ss->distOH[kk]) - (1.0/inp_ss->distCN[kk]) );
     kk ++;
    }
   
   return 0;
}
// ------------------------------------------------------------------
int sStruct( struct inp_ss *inp_ss, char *secstruct )
{
   int  ii, jj;
   float Rkj[3], Rkl[3];
   float angle;
   int   nhbonds;
  
    //  *** START OF THE SECONDARY STRUCTURE PREDICTION ***  //
    
   for( ii=0; ii<inp_ss->ntotres; ii++ )
    {
     if( secstruct[ii] == 'H')
      continue;
    // >>> check for H-HELICES <<< ============================================
//   if( (ii>1 && ii<inp_ss->ntotres-3) && inp_ss->hbonds[ii-1][ii+2] && inp_ss->hbonds[ii][ii+3] ) // formula according to paper
     if( (ii>1 && ii<inp_ss->ntotres-4) && inp_ss->hbonds[ii-1][ii+3] && inp_ss->hbonds[ii][ii+4] ) // formula according to sense
     {
      secstruct[ii] = 'H';
      
      secstruct[ii+1] = 'H';
      secstruct[ii+2] = 'H';
      secstruct[ii+3] = 'H';
      
      continue;
     }
    // aa still helix if 1 bond is missing and prev 2 and following 2 aa are H
     if( (ii>2 && ii<inp_ss->ntotres-6) &&                                      // if we are inside the limits
          secstruct[ii-1]=='H' && secstruct[ii-2]=='H' &&       // previous 2 aa are H
         (inp_ss->hbonds[ii-1][ii+3] || inp_ss->hbonds[ii][ii+4] ) &&           // and only 1 hbond is missing
         (inp_ss->hbonds[ii+1][ii+5] || inp_ss->hbonds[ii+2][ii+6] ) )          // and there's another H then
      {
      secstruct[ii] = 'H'; //hgbhtfytdfytgrd and next is proper helix!!!!
      continue;
     }
     if( (ii>2 && ii<inp_ss->ntotres-3) &&                                      // if we are inside a H-helix couple
          secstruct[ii-1]=='H' && secstruct[ii-2]=='H' &&       // >>44<<
         ((inp_ss->hbonds[ii-1][ii+3] && inp_ss->hbonds[ii-2][ii+2] ) ||        // 
          (inp_ss->hbonds[ii-3][ii+1] && inp_ss->hbonds[ii-2][ii+2] )) )        // 
      {
      secstruct[ii] = 'H'; //hgbhtfytdfytgrd and next is proper helix!!!!
      continue;
     }
    // aaand aa is helix also if it is on the receiving end of the bond
     if( (ii>4 && ii<inp_ss->ntotres-1) && inp_ss->hbonds[ii-4][ii] && (inp_ss->hbonds[ii-3][ii+1]))// || inp_ss->hbonds[ii-5][ii-1]))
//   if( (ii>4 && ii<inp_ss->ntotres-1) && inp_ss->hbonds[ii-4][ii] && (inp_ss->hbonds[ii-3][ii+1] || inp_ss->hbonds[ii-5][ii-1]))
     {
      secstruct[ii] = 'H';
      continue;
     }
     
    // >>> check for BRIDGES (B or E) <<< =============================================
     for( jj=1; jj<inp_ss->ntotres-1; jj++ )
     {
      if( (ii>0 && ii<inp_ss->ntotres-1) && (jj < ii-1 || jj > ii+1) &&         // avoid out-of-memory checks AND stretches must not overlap to define a bridge
          ( (inp_ss->hbonds[ii-1][jj] && inp_ss->hbonds[jj][ii+1]) ||
            (inp_ss->hbonds[jj-1][ii] && inp_ss->hbonds[ii][jj+1]) ) )
      {
       secstruct[ii] = 'B';
       inp_ss->ssnotes[ii] = 'P';       // parallel flag
       break;
      }
      if( (ii>0 && ii<inp_ss->ntotres-1) && (jj < ii-1 || jj > ii+1) &&         // avoid out-of-memory checks AND stretches must not overlap to define a bridge
          ( (inp_ss->hbonds[ii][jj] && inp_ss->hbonds[jj][ii])     ||
            (inp_ss->hbonds[ii-1][jj+1] && inp_ss->hbonds[jj-1][ii+1]) ) )
      {
       secstruct[ii] = 'B';
       inp_ss->ssnotes[ii] = 'A';       // antiparallel flag
       break;
      }
     } 
     if( secstruct[ii] == 'B' &&
         (secstruct[ii-1] == 'B' || secstruct[ii-1] == 'E') &&
         (secstruct[ii-2] == 'B' || secstruct[ii-2] == 'E')    )
     {
      secstruct[ii] = 'E';
      secstruct[ii-1] = 'E';
      secstruct[ii-2] = 'E';
     }
     // alternate option: 2 BBs are enough to set all to EE
     if( secstruct[ii] == 'B' &&
         (secstruct[ii-1] == 'B' || secstruct[ii-1] == 'E') )
     {
      secstruct[ii] = 'E';
      secstruct[ii-1] = 'E';
     }
     
     if( secstruct[ii] == 'B' && 
        (secstruct[ii-2] == 'B' || secstruct[ii-2] == 'E') &&
         inp_ss->ssnotes[ii] == inp_ss->ssnotes[ii-2] )
     {
      secstruct[ii] = 'E';
      secstruct[ii-1] = 'E';
      secstruct[ii-2] = 'E';
     }
     
     if( secstruct[ii] == 'B' || secstruct[ii] == 'E')
      continue;
     
    // >>> check for ALT.HELICES <<< ===========================================
    // G (3_10 helix)
//   if( (ii>1 && ii<inp_ss->ntotres-2) &&  inp_ss->hbonds[ii-1][ii+1] && inp_ss->hbonds[ii][ii+2] ) // formula according to paper
     if( (ii>1 && ii<inp_ss->ntotres-3) &&  inp_ss->hbonds[ii-1][ii+2] && inp_ss->hbonds[ii][ii+3] ) // formula according to sense
     {
      secstruct[ii] = 'G';
      continue;
     }
     if( (ii>2 && ii<inp_ss->ntotres-2) &&  inp_ss->hbonds[ii-2][ii+1] && inp_ss->hbonds[ii-1][ii+2] ) // formula according to sense
     {
      secstruct[ii] = 'G';
      continue;
     }
     if( (ii>3 && ii<inp_ss->ntotres-1) &&  inp_ss->hbonds[ii-3][ii] && inp_ss->hbonds[ii-2][ii+1] ) // formula according to sense
     {
      secstruct[ii] = 'G';
      continue;
     }
     
    // I (pi helix)
     if( (ii>1 && ii<inp_ss->ntotres-5) &&  inp_ss->hbonds[ii-1][ii+4] && inp_ss->hbonds[ii][ii+5] )
     {
      secstruct[ii] = 'I';
      continue;
     }
     if( ( ii>4 && ii<inp_ss->ntotres-1 && inp_ss->hbonds[ii-4][ii+1] ) ||
         ( ii>3 && ii<inp_ss->ntotres-1 && inp_ss->hbonds[ii-3][ii+1] ) || 
         ( ii>3 && ii<inp_ss->ntotres-2 && inp_ss->hbonds[ii-3][ii+2] ) || 
         ( ii>2 && ii<inp_ss->ntotres-1 && inp_ss->hbonds[ii-2][ii+1] ) || 
         ( ii>2 && ii<inp_ss->ntotres-2 && inp_ss->hbonds[ii-2][ii+2] ) || 
         ( ii>2 && ii<inp_ss->ntotres-3 && inp_ss->hbonds[ii-2][ii+3] ) || 
         ( ii>1 && ii<inp_ss->ntotres-2 && inp_ss->hbonds[ii-1][ii+2] ) || 
         ( ii>1 && ii<inp_ss->ntotres-3 && inp_ss->hbonds[ii-1][ii+3] ) || 
         ( ii>1 && ii<inp_ss->ntotres-4 && inp_ss->hbonds[ii-1][ii+4] )    )
     {
      secstruct[ii] = 'T';
      continue;
     }
     
    // check for BENDS  (S)
     if( ii>1 && ii<inp_ss->ntotres-2)
     {
      Rkl[0] = inp_ss->xcoor5[ii-2] - inp_ss->xcoor5[ii];
      Rkl[1] = inp_ss->ycoor5[ii-2] - inp_ss->ycoor5[ii];
      Rkl[2] = inp_ss->zcoor5[ii-2] - inp_ss->zcoor5[ii];
      Rkj[0] = inp_ss->xcoor5[ii]   - inp_ss->xcoor5[ii+2];
      Rkj[1] = inp_ss->ycoor5[ii]   - inp_ss->ycoor5[ii+2];
      Rkj[2] = inp_ss->zcoor5[ii]   - inp_ss->zcoor5[ii+2];
      w_norm(Rkl);
      w_norm(Rkj);
      angle = acos( dot_prod(Rkl, Rkj) );
      angle*= 180/3.1415927;
      if( angle > 70 )
      {
       secstruct[ii] = 'S';
       continue;
      }
     }

     nhbonds = 0;
     for( jj=0; jj<inp_ss->ntotres; jj++ )
      if( ii!=jj )
        nhbonds += inp_ss->hbonds[ii][jj] + inp_ss->hbonds[jj][ii];
     if( nhbonds == 0 )
      secstruct[ii] = 'L';
    }

   return 1;
}
// ------------------------------------------------------------------
