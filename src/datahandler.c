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
/*! \file datahandler.c
 \brief Data manipulation functions
 
 source code for all data manipulation functions. 
*/
#include "wordom.h"
#include "fileio.h"
#include "tools.h"
#include "datahandler.h"

extern short int verbose;

// ------------------------------------------------------------------
void Trj_Merge ( struct sopt   *OPT    )
{
   int          ii = 0, jj = 0;
   int          step = 1;
   Traj        *traj1, *traj2, *traj3;
   TrjList      trj_list;
   CoorSet     *trj_crd;
   
   if( !wrd_isdattxt(OPT->ITRJ_FILE) )
   {
     fprintf( stderr, "Error: attempting trajectory merge but -itrj argument is not list (.txt/.dat)\n");
     exit(0);
   }
   ReadTrjList    ( OPT->ITRJ_FILE, &trj_list);
   printf("Found %d trj files to merge\n", trj_list.ntrj);
           
   traj1 = InitTrj ( trj_list.trjname[0], "r" );
   traj3 = InitTrj ( OPT->OTRJ_FILE, "wb+" );
   
   CopyTrjHeader ( traj1->hdr, traj3->hdr);
   traj3->hdr->sixtyfour = 0;
   
   if ( OPT->SKIP_FLAG )
   {
    if ( fmod( OPT->SKIP_STEP, traj1->hdr->skpstp ) )
    {
     printf ( "Wrong skip step! check with wordom -head\n" );
     exit(97);
    }
    step = ( OPT->SKIP_STEP / traj1->hdr->skpstp ) ;
   }

   WriteTrjHeader ( traj3, "all" );

   trj_crd = InitCoor(traj1->hdr->nato);
   
   for ( ii = 0 ; ii <traj1->hdr->nframe ; ii+=step )
   {
     ReadTrjCoor  ( traj1, trj_crd, ii+1 );
     WriteTrjCoor ( trj_crd, traj3 );
   }
   
   if ( trj_list.ntrj > 1)
   {
     for ( ii=1 ; ii < trj_list.ntrj ; ii++)
     {
       traj2 = InitTrj ( trj_list.trjname[ii], "r" );
       
       if (traj1->hdr->nato   != traj2->hdr->nato   ||
           traj1->hdr->skpstp != traj2->hdr->skpstp ||
           traj1->hdr->tmstp  != traj2->hdr->tmstp    )
       {
         fprintf ( stderr, "Different Trajectories! better check... and don't use -skip!\n" ) ;
         exit ( 98 ) ;
       }
       
       traj3->hdr->nframe += traj2->hdr->nframe;
       traj3->hdr->numstp += traj2->hdr->numstp;
       
       for ( jj = 0 ; jj < traj2->hdr->nframe ; jj+=step )
       {
        ReadTrjCoor  ( traj2, trj_crd, jj+1 );
        WriteTrjCoor ( trj_crd, traj3 );
       }
            
       CloseTrj ( traj2 );
     }
   }

   if ( OPT->SKIP_FLAG )
   {
    traj3->hdr->nframe /= step ;
    traj3->hdr->skpstp  = OPT->SKIP_STEP;
   }
   
   printf("Correcting nframe (new nframe = %d)\n", traj3->hdr->nframe);
   WriteTrjHeader ( traj3, "nframe" );
   WriteTrjHeader ( traj3, "skpstp" );
   WriteTrjHeader ( traj3, "numstp" );
   
   DelCoor( trj_crd );
   
   CloseTrj( traj1 );
   CloseTrj( traj3 );
   
   return;
}
// ------------------------------------------------------------------
void Trj_append ( struct sopt   *OPT    )
{
   int          ii = 0;
   int          step = 1;
   Traj        *traj1, *traj2;
   CoorSet     *trj_crd;
   
   traj1 = InitTrj( OPT->ATRJ_FILE, "r" );
   traj2 = InitTrj( OPT->OTRJ_FILE, "rb+" );
   
   if( traj1->xtc || traj2->xtc )
   {
     fprintf( stderr, "appending with xtc files not perfected yet\n");
     exit( 17 );
   }
   
   if ( OPT->SKIP_FLAG )
   {
     if ( fmod( OPT->SKIP_STEP, traj1->hdr->skpstp ) )
     {
       fprintf ( stderr, "Wrong skip step! check with wordom -head\n" );
       exit(97);
     }
     step = ( OPT->SKIP_STEP / traj1->hdr->skpstp ) ;
   }

   if( (traj1->hdr->nato!= traj2->hdr->nato) )
    fprintf( stderr, "Number of atoms differs in trajectories! Watch out!\n" );
   
   if( traj1->hdr->fixatm != 0 || traj2->hdr->fixatm != 0 )
   {
     fprintf( stderr, "Fixed atoms present! sorry, Wordom cannot handle this...\n" );
     fprintf( stderr, "You might convert (with wordom) to no-fixed-atoms TRJs and try again\n" );
     exit(0);
   }

   trj_crd = InitCoor( traj1->hdr->nato );
   
   for ( ii = 0 ; ii < traj1->hdr->nframe ; ii+=step )
   {
    ReadTrjCoor  ( traj1, trj_crd, ii+1 );
    WriteTrjCoor ( trj_crd, traj2 );
   }

   traj2->hdr->nframe += traj1->hdr->nframe/step;
   traj2->hdr->numstp += traj1->hdr->numstp/step;
   WriteTrjHeader ( traj2, "nframe" );
   WriteTrjHeader ( traj2, "skpstp" );
   WriteTrjHeader ( traj2, "numstp" );
   
   DelCoor( trj_crd );
   
   CloseTrj ( traj1 );
   CloseTrj ( traj2 );
   
   return;
}
// ------------------------------------------------------------------
void Trj_sum ( struct sopt *OPT )
{
   // each frame is the sum of the corrisponding frames from the listed dcds
   // the coordinates of the specified structure are substracted from
   // each dcd coordinates and added to the final sum
   int          ii, jj, kk;
   TrjList      trj_list;
   Molecule    *molecule;
   Traj       **trajs;
   Traj        *outtraj;
   CoorSet     *trj_crd1, *trj_crd2;
   Selection    sele1;
   trjtype      type;
   
   molecule = ReadMolecule ( OPT->IMOL_FILE, OPT->IMOL_TYPE );
   GetSele ( OPT->SELE_STRING, &sele1, molecule);
   
   ReadTrjList( OPT->TRJS_FILE, &trj_list);
   
   type = wrd_whichtrjtype( trj_list.trjname[0] );
   if( type == xtc )
   {
     fprintf( stderr, "Sorry, summing xtc files still gives some problems. Try converting to dcd first\n");
     exit(0);
   }

   trajs = calloc( trj_list.ntrj, sizeof( Traj *) );
   
   for( ii=0; ii<trj_list.ntrj; ii++)
     trajs[ii] = InitTrj ( trj_list.trjname[ii], "r");
   
   for( ii=1; ii<trj_list.ntrj; ii++)
    if(trajs[ii]->hdr->nato != trajs[ii-1]->hdr->nato)
    {
     fprintf(stderr, "Number of atoms differ in trajectories! Check!\n" );
     exit(4);
    }
   
   trj_crd1 = InitCoor(trajs[0]->hdr->nato);
   trj_crd2 = InitCoor(trajs[0]->hdr->nato);
   
   outtraj = InitTrj ( OPT->OTRJ_FILE, "wb");
   CopyTrjHeader ( trajs[0]->hdr, outtraj->hdr);
   outtraj->hdr->sixtyfour = 0;
   WriteTrjHeader ( outtraj, "all" );
   
   for( kk=0; kk<trajs[0]->hdr->nframe; kk++)
   {
   
    for(jj=0; jj<3*trajs[0]->hdr->nato; jj++)
     trj_crd2->cords[jj] = molecule->coor.cords[jj];
    
    for( ii=0; ii<trj_list.ntrj; ii++)
    {
     ReadTrjCoor ( trajs[ii], trj_crd1, kk+1 );
     for(jj=0; jj<sele1.nselatm; jj++)
     {
      trj_crd2->xcoor[sele1.selatm[jj]-1] += (trj_crd1->xcoor[sele1.selatm[jj]-1] - molecule->coor.xcoor[sele1.selatm[jj]-1]);
      trj_crd2->ycoor[sele1.selatm[jj]-1] += (trj_crd1->ycoor[sele1.selatm[jj]-1] - molecule->coor.ycoor[sele1.selatm[jj]-1]);
      trj_crd2->zcoor[sele1.selatm[jj]-1] += (trj_crd1->zcoor[sele1.selatm[jj]-1] - molecule->coor.zcoor[sele1.selatm[jj]-1]);
     }
    }
   
    WriteTrjCoor ( trj_crd2, outtraj );
   }
   
   for( ii=0; ii<trj_list.ntrj; ii++)
     CloseTrj ( trajs[ii] );
   CloseTrj ( outtraj );

   DelCoor( trj_crd1 );
   DelCoor( trj_crd2 );
   
   return;
}
// ------------------------------------------------------------------
void Trj_combine ( struct sopt *OPT )
{
  // each frame has the atoms of the corrisponding frames from the listed dcds
  // final number of atoms for each frame is the sum of the number of atoms
  // for each listed dcd - now using max of 2 trjs
  
  int          ii, jj, kk, ll;
  TrjList      trj_list;
  Traj       **trajs;
  Traj        *outtraj;
  CoorSet    **trj_crds;
  CoorSet     *out_crd;
  
  ReadTrjList( OPT->TRJC_FILE, &trj_list);
  
  trajs = calloc( trj_list.ntrj, sizeof( Traj *) );
  trj_crds = calloc( trj_list.ntrj, sizeof(CoorSet *));
  
  for( ii=0; ii<trj_list.ntrj; ii++)
    trajs[ii] = InitTrj ( trj_list.trjname[ii], "r");
  
  for( ii=1; ii<trj_list.ntrj; ii++)
   if(trajs[ii]->hdr->nframe != trajs[ii-1]->hdr->nframe)
   {
    fprintf(stderr, "Number of frames differ in trajectories! Check!\n" );
    exit(4);
   }
  
  for( ii=0; ii<trj_list.ntrj; ii++)
    trj_crds[ii] = InitCoor( trajs[ii]->hdr->nato );
  
  outtraj = InitTrj ( OPT->OTRJ_FILE, "wb");
  CopyTrjHeader ( trajs[0]->hdr, outtraj->hdr);
  for( ii=1; ii<trj_list.ntrj; ii++)
    outtraj->hdr->nato += trajs[ii]->hdr->nato;
  fprintf( stdout, "# Output trajectory will have %d atoms\n", outtraj->hdr->nato);
  outtraj->hdr->sixtyfour = 0;
  WriteTrjHeader ( outtraj, "all" );
  out_crd = InitCoor( outtraj->hdr->nato );
  out_crd->pbc = calloc( 1, sizeof( Pbc ));
  out_crd->pbc->cpt = 0;
  out_crd->pbc->a_size = 0;
  out_crd->pbc->b_size = 0;
  out_crd->pbc->c_size = 0;
  out_crd->pbc->angle1 = 0;
  out_crd->pbc->angle2 = 0;
  out_crd->pbc->angle3 = 0;
  out_crd->pbc->zval = 0;
  
  for( ii=0; ii<trajs[0]->hdr->nframe; ii++)
  {
    ll = 0;
    for( jj=0; jj<trj_list.ntrj; jj++ )
    {
      ReadTrjCoor ( trajs[jj], trj_crds[jj], ii+1 );
      for(kk=0; kk<trajs[jj]->hdr->nato; kk++)
      {
        //out_crd->cords[ll] = trj_crds[jj]->cords[kk];
        out_crd->xcoor[ll] = trj_crds[jj]->xcoor[kk];
        out_crd->ycoor[ll] = trj_crds[jj]->ycoor[kk];
        out_crd->zcoor[ll] = trj_crds[jj]->zcoor[kk];
        ll += 1;
        //fprintf( stdout, "%5d", ll );
      }
    }
    WriteTrjCoor ( out_crd, outtraj );
  }
   
  for( ii=0; ii<trj_list.ntrj; ii++)
  {
    CloseTrj ( trajs[ii] );
    DelCoor( trj_crds[ii] );
  }
  CloseTrj ( outtraj );
  free( trj_crds );
  DelCoor( out_crd );
  
  return;
}
// ------------------------------------------------------------------
void Mol_append ( struct sopt   *OPT    )
{
   int          step = 1;
   Molecule    *molecule;
   Traj        *traj;
   
   molecule = ReadMolecule ( OPT->AMOL_FILE, OPT->AMOL_TYPE );

   traj = InitTrj ( OPT->OTRJ_FILE, "rb+" );
   if( (traj->hdr->nato!= molecule->nato) )
   {
     fprintf( stderr, "Number of atoms differ in trajectory(%d) and pdb file(%d)! Check!\n", traj->hdr->nato, molecule->nato);
     exit(0);
   }

   WriteTrjCoor ( &molecule->coor, traj );

   if ( OPT->SKIP_FLAG )
   {
    traj->hdr->nframe /= step ;
    traj->hdr->skpstp =  OPT->SKIP_STEP;
   }
   traj->hdr->nframe ++;
   traj->hdr->numstp += traj->hdr->tmstp;
   WriteTrjHeader ( traj, "nframe" );
   WriteTrjHeader ( traj, "numstp" );
   
   CloseTrj( traj );
   
   return;
}
// ------------------------------------------------------------------
void Head_Mod ( struct sopt   *OPT    )
{
   Traj   *traj;
   char    field  [16];
   int     value;
   float   valuef;
   
   traj = InitTrj ( OPT->ITRJ_FILE, "rb+" );
   
   if ( !strstr (OPT->MOD_STRING, "=") )
    {
     printf("Argument of -mod is not in the format field=value\n");
     exit (0);
    }
   sscanf (OPT->MOD_STRING, "%[^=]%*c", field);
   
   if ( strcmp (field, "begstp") && strcmp (field, "numstp") && strcmp (field, "nframe") && strcmp (field, "tmstp") && strcmp (field, "skpstp") && strcmp (field, "varpbc") && strcmp (field, "64made") )
    {
      printf("Unknown field in -mod option: %s\n", field);
      exit(0);
    }
   printf("Field is: %s\n", field);
   if ( (!strcmp (field, "begstp")) || !strcmp (field, "numstp") || !strcmp (field, "nframe") || !strcmp (field, "skpstp") || !strcmp (field, "varpbc") )
    {
     sscanf( OPT->MOD_STRING, "%*[^=]%*c%d", &value );
     traj->hdr->nframe = traj->hdr->begstp = traj->hdr->numstp = traj->hdr->skpstp = traj->hdr->varpbc = value;
     printf ("I read %s: %d\n", field, value);
    }
   else if ( !strcmp (field, "tmstp"))
    {
     sscanf( OPT->MOD_STRING, "%*[^=]%*c%f", &valuef );
     traj->hdr->tmstp = valuef;
     printf ("I read  tmstp: %f\n", valuef);
     printf ("I write tmstp: %f\n", traj->hdr->tmstp);
    }
    else if ( !strcmp (field, "64made"))
    {
     sscanf( OPT->MOD_STRING, "%*[^=]%*c%d", &value );
     traj->hdr->sixtyfour = value;
     printf ("I read 64made: %d\n", value);
    }
  
   WriteTrjHeader ( traj, field );
   CloseTrj ( traj );
   
   return;
}
// ------------------------------------------------------------------
void MkMono ( struct sopt *OPT )
{
   Traj        *traj;
   Molecule    *molecule;
   
   molecule = ReadMolecule ( OPT->IMOL_FILE, OPT->IMOL_TYPE );
   traj = InitTrj ( OPT->OTRJ_FILE, "w" );

   sprintf ( traj->hdr->type, "%s", "CORD" );
   traj->hdr->nframe =   1 ;
   traj->hdr->begstp = 100 ;
   traj->hdr->skpstp = 100 ;
   traj->hdr->numstp = 100 ;
   traj->hdr->tmstp  =   0.04091 ;
   traj->hdr->nato   = molecule->nato ;
   traj->hdr->sixtyfour    = 0;
   
   WriteTrjHeader ( traj, "all" );
   WriteTrjCoor   ( &molecule->coor, traj  );

   return;
}
// ------------------------------------------------------------------
void MolConv ( struct sopt *OPT )
{
   Molecule    *molecule, *selmol;
   Selection    selestr;
   int          ii;
   char        *outname, suffix[4];
   
   molecule = ReadMolecule( OPT->IMOL_FILE, OPT->IMOL_TYPE );
   
   outname = calloc( strlen(OPT->OMOL_FILE)+3, sizeof(char));
   if( OPT->OMOL_TYPE == 1 )
     sprintf( suffix, "pdb");
   else if( OPT->OMOL_TYPE == 2 )
     sprintf( suffix, "crd");
   
   if(OPT->SELE_FLAG)
   {
     GetSele( OPT->SELE_STRING, &selestr, molecule);
     if (selestr.nselatm == 0)
     {
       fprintf( stderr, "Error: no atom selected\n");
       exit(95);
     }
     selmol = MakeSeleMol ( &selestr , molecule); //, &selmol );
     
     for ( ii=0; ii<selmol->nato; ii++)
     {
       selmol->coor.xcoor[ii] = molecule->coor.xcoor[selestr.selatm[ii]-1];
       selmol->coor.ycoor[ii] = molecule->coor.ycoor[selestr.selatm[ii]-1];
       selmol->coor.zcoor[ii] = molecule->coor.zcoor[selestr.selatm[ii]-1];
     }
     
     sprintf( outname, "%s", OPT->OMOL_FILE );
     if( molecule->altLocB )
       sprintf( &outname[strlen(OPT->OMOL_FILE)-4], "_A.%s", suffix);
     WriteMolecule_unstr ( selmol, outname, OPT->OMOL_TYPE );
     
     if( molecule->altLocB )
     {
       sprintf( &outname[strlen(OPT->OMOL_FILE)-4], "_B.%s", suffix);
       for( ii=0; ii<selmol->nato; ii++ )
       {
         selmol->coor.xcoor[ii] = molecule->coor2.xcoor[selestr.selatm[ii]-1];
         selmol->coor.ycoor[ii] = molecule->coor2.ycoor[selestr.selatm[ii]-1];
         selmol->coor.zcoor[ii] = molecule->coor2.zcoor[selestr.selatm[ii]-1];
       }
       WriteMolecule_unstr ( selmol, outname, OPT->OMOL_TYPE );
     }
     
     if( molecule->altLocC )
     {
       sprintf( &outname[strlen(OPT->OMOL_FILE)-4], "_C.%s", suffix);
       for( ii=0; ii<selmol->nato; ii++ )
       {
         selmol->coor.xcoor[ii] = molecule->coor3.xcoor[selestr.selatm[ii]-1];
         selmol->coor.ycoor[ii] = molecule->coor3.ycoor[selestr.selatm[ii]-1];
         selmol->coor.zcoor[ii] = molecule->coor3.zcoor[selestr.selatm[ii]-1];
       }
       WriteMolecule_unstr ( selmol, outname, OPT->OMOL_TYPE );
     }
   }
   else
   {
     sprintf( outname, "%s", OPT->OMOL_FILE );
     if( molecule->altLocB )
       sprintf( &outname[strlen(OPT->OMOL_FILE)-4], "_A.%s", suffix);
     WriteMolecule_unstr ( molecule, outname, OPT->OMOL_TYPE );
     
     if( molecule->altLocB )
     {
       sprintf( &outname[strlen(OPT->OMOL_FILE)-4], "_B.%s", suffix);
       for( ii=0; ii<3*molecule->nato; ii++ )
         molecule->coor.cords[ii] = molecule->coor2.cords[ii];
       WriteMolecule_unstr ( molecule, outname, OPT->OMOL_TYPE );
     }
     if( molecule->altLocC )
     {
       sprintf( &outname[strlen(OPT->OMOL_FILE)-4], "_C.%s", suffix);
       for( ii=0; ii<3*molecule->nato; ii++ )
         molecule->coor.cords[ii] = molecule->coor3.cords[ii];
       WriteMolecule_unstr ( molecule, outname, OPT->OMOL_TYPE );
     }
   }
   DelMolecule(molecule);
   
   return;
}
// ------------------------------------------------------------------
void Dcd2Xtc ( struct sopt   *OPT )
{
   char        *func="Dcd2Xtc";
   FILE        *idcd_f;
   Trjh         dcd_hdr;   // header
   CoorSet     *dcd_crd;   // coordinates
   int          aa;        // counter
   Traj        *traj;

   Molecule    *molecule;
   int          ii;
   
   double       tmstp;
   int          skpstp;

   // The following belongs to the xtc file
   int           natom;
   int           step;
   double        time;
   matrix        box;
   float         fa,fb,fc;       // box elements in nm
   float         alpha,beta,gamma; //
   rvec          *x;
   float         precision=10000000;//precision=10000;
   int           result;

   // The following does not belong to the xtc file
   double        cosa,cosb,cosg; 
   double        sing; 
     
   molecule = ReadMolecule ( OPT->IMOL_FILE, OPT->IMOL_TYPE);
   idcd_f = O_File ( OPT->ITRJ_FILE, "r" );
   ReadDcdHeader  ( idcd_f, &dcd_hdr );
   
   if ( dcd_hdr.nato != molecule->nato)
   {
     fprintf( stderr, "ERROR:\n Number of atoms in dcd differs from that in pdb!\n");
     exit(98);
   } 

   step  = dcd_hdr.begstp;
   tmstp = floor(dcd_hdr.tmstp*AKMA2FS*10+0.5)/10000;
   skpstp = dcd_hdr.skpstp;
   if (step==0)
       time=0;
   else 
       time  = step * tmstp;
   natom = molecule->nato; // change this later according to the selection
   dcd_crd = InitCoor( natom );
   x=calloc(natom,sizeof(*x));

// clear_mat(box);
   box[XX][XX]=box[XX][YY]=box[XX][ZZ]=0.0;
   box[YY][XX]=box[YY][YY]=box[YY][ZZ]=0.0;
   box[ZZ][XX]=box[ZZ][YY]=box[ZZ][ZZ]=0.0;

   if (molecule->rawmol.b_cryst) {
       // setting up the box
       alpha=molecule->rawmol.alpha;
       beta=molecule->rawmol.beta;
       gamma=molecule->rawmol.gamma;
       fa = (molecule->rawmol.a) * 0.1;
       fb = (molecule->rawmol.b) * 0.1;
       fc = (molecule->rawmol.c) * 0.1; 
       box[XX][XX] = fa;
       if ( (alpha != 90.0) || 
	    (beta  != 90.0) || 
	    (gamma != 90.0) ) 
       {
         if (alpha != 90.0) {
	   cosa = cos(alpha*DEG2RAD);
	 } else {
	   cosa = 0;
	 }
	 if (beta != 90.0) {
	   cosb = cos(beta*DEG2RAD);
	 } else {
	   cosb = 0;
	 }
	 if (gamma != 90.0) {
	   cosg = cos(gamma*DEG2RAD);
	   sing = sin(gamma*DEG2RAD);
	 } else {
	   cosg = 0;
	   sing = 1;
	 }
	 box[YY][XX] = fb*cosg;
	 box[YY][YY] = fb*sing;
	 box[ZZ][XX] = fc*cosb;
	 box[ZZ][YY] = fc*(cosa - cosb*cosg)/sing;
	 box[ZZ][ZZ] = sqrt(fc*fc
			     -box[ZZ][XX]*box[ZZ][XX]-box[ZZ][YY]*box[ZZ][YY]);
	} else {
	  box[YY][YY] = fb;
	  box[ZZ][ZZ] = fc;
	}
   }

   if ( OPT->BEG_FLAG == 0 )
   {
      if ( OPT->SKIP_FLAG == 0 )
      OPT->TRJN_BEG = 1;
      if ( OPT->SKIP_FLAG == 1 )
      OPT->TRJN_BEG = OPT->SKIP_STEP;
   }

   if ( OPT->END_FLAG == 0 )
     OPT->TRJN_END = dcd_hdr.nframe;   
   
   traj = InitTrj( OPT->OTRJ_FILE, "w" );

   for ( aa=OPT->TRJN_BEG; aa<=OPT->TRJN_END; aa=aa+OPT->SKIP_STEP )
   {
     ReadDcdCoor ( idcd_f, dcd_crd, &dcd_hdr, aa );
     for ( ii=0; ii< natom; ii++) {
	 x[ii][XX]=dcd_crd->xcoor[ii]*0.1;
	 x[ii][YY]=dcd_crd->ycoor[ii]*0.1;
	 x[ii][ZZ]=dcd_crd->zcoor[ii]*0.1;
     }
     /* added to cope with var-pbc in dcd */
     if( dcd_crd->pbc->cpt == 1 )
     {
       fa = (dcd_crd->pbc->a_size) * 0.1;
       fb = (dcd_crd->pbc->b_size) * 0.1;
       fc = (dcd_crd->pbc->c_size) * 0.1; 
       box[XX][XX] = fa;
       
     }
     
     result=write_xtc(traj->xtc, natom, step, time, box, x, precision);
     if (0 != result) {
	 printf ("Error in %s. Could not write frame to file %s. Expect trouble\n",func,OPT->OTRJ_FILE);
	 break;
     }
     step += skpstp;
     time += skpstp * tmstp;
   }
   CloseTrj(traj);
   fclose        ( idcd_f );
   rfree(x);
   DelCoor(dcd_crd);

   return;
}
 // ------------------------------------------------------------------
void Xtc2Dcd ( struct sopt   *OPT )
{
    char        *func="Xtc2Dcd";

     FILE       *odcd_f;
     Traj       *traj; 
     Trjh        dcd_hdr;
     
     int        nframes=0;
     float      tmstp;
     int        begstp=0;
     CoorSet    this_coor;
     int        ii;

     // The following belongs to the xtc file
     int           natom;
     int           step;
     float         time;
     matrix        box;
     rvec          *x;
     float         precision;
     int           result;

     traj = InitTrj( OPT->ITRJ_FILE, "r" );
     odcd_f = O_File ( OPT->OTRJ_FILE, "w" );
     xdrfile_close(traj->xtc);
     read_xtc_natoms(traj->filename, &natom);
     x = calloc(natom,sizeof(*x));
     traj->xtc = xdrfile_open(traj->filename,"r");
     result=read_xtc(traj->xtc,natom,&step,&time,box,x,&precision);
     if (exdrOK != result) {
       printf ("Wordom: %s\n\tError! could not read the first XTC coordinate\n",func);
       rfree(x);
       exit(99);
     }
     begstp = step;
     // initialising DCD header
     sprintf ( dcd_hdr.type, "%s", "CORD" );
     dcd_hdr.nframe = nframes ;
     dcd_hdr.varpbc = 0 ; // all info about PBC is lost in conversion
     dcd_hdr.begstp = begstp;     
     dcd_hdr.nato = natom;
     dcd_hdr.sixtyfour    = 0;
     WriteDcdHeader ( odcd_f, &dcd_hdr, "all" );
     while (exdrOK==result) {
       // memalloc for the coordinates
       this_coor.xcoor =  calloc ( natom, sizeof(float) );
       this_coor.ycoor =  calloc ( natom, sizeof(float) );
       this_coor.zcoor =  calloc ( natom, sizeof(float) );
       // end of memory allocation
       for ( ii=0; ii< natom; ii++) {
	 this_coor.xcoor[ii]=x[ii][XX]/0.1;
	 this_coor.ycoor[ii]=x[ii][YY]/0.1;
	 this_coor.zcoor[ii]=x[ii][ZZ]/0.1;
       }
       WriteDcdCoor ( odcd_f, &this_coor, &dcd_hdr );
       free(this_coor.xcoor);
       free(this_coor.ycoor);
       free(this_coor.zcoor);
       nframes++;
       result=read_xtc(traj->xtc, natom, &step, &time, box, x, &precision);
     }
     xdrfile_close(traj->xtc);

     if (nframes > 1)
	 dcd_hdr.skpstp = (int) ( (step-begstp)/ (nframes -1) );
     else 
	 dcd_hdr.skpstp = 1;
     if (time==0 || step ==0) {
	 tmstp = 0;
     }
     else {
	 tmstp = time / (step * AKMA2PS);
     }
     dcd_hdr.tmstp = tmstp;
     dcd_hdr.nframe = nframes;
     dcd_hdr.numstp = step-begstp;

     // correcing dcd header
     WriteDcdHeader( odcd_f, &dcd_hdr, "nframe" );
     WriteDcdHeader( odcd_f, &dcd_hdr, "skpstp" );
     WriteDcdHeader( odcd_f, &dcd_hdr, "numstp" ); 
     WriteDcdHeader( odcd_f, &dcd_hdr, "tmstp" ); 

     fclose( odcd_f );
}
// ------------------------------------------------------------------
void TrjExtract ( struct sopt   *OPT    )
{
   Traj        *traj1, *traj2, *tmptraj;
   CoorSet     *trj_crd, *selcrd;
   IntList      frame_list;
   TrjList      trj_list;
   Selection    selestr;
   int          ii, jj, kk, ll, totpframes, totframesintrjs=0, sortedlist;
   
   Molecule    *molecule;
   
   if( wrd_isdattxt(OPT->ITRJ_FILE) )
   {
     ReadTrjList    ( OPT->ITRJ_FILE, &trj_list);
     traj1 = InitTrj ( trj_list.trjname[0], "r" );
     totframesintrjs = traj1->hdr->nframe;
     for( ii=1; ii<trj_list.ntrj; ii++)
     {  
       tmptraj = InitTrj( trj_list.trjname[ii], "r" );
       totframesintrjs += tmptraj->hdr->nframe;
       CloseTrj( tmptraj );
     }
   }
   else
   {
     traj1 = InitTrj ( OPT->ITRJ_FILE, "r" );
     totframesintrjs = traj1->hdr->nframe;
   }
   
   traj2 = InitTrj ( OPT->OTRJ_FILE, "wb+" );
   CopyTrjHeader ( traj1->hdr, traj2->hdr);
   traj2->hdr->sixtyfour = 0;
   
   if( !strncmp( OPT->FRL_FILE, "all", 3 ))
   {
     if( totframesintrjs == 0 )
     {
       fprintf( stderr, "Frame number = 0 found.\nNote: the \"all\" keyword won't work with xtc files\nUse \"range\", a frame list or first convert to a dcd with the -conv option\n");
       exit(0);
     }
     frame_list.nframes = totframesintrjs;
     frame_list.nframe = malloc( frame_list.nframes * sizeof(int));
     for( ii=0; ii<frame_list.nframes; ii++ )
       frame_list.nframe[ii] = ii+1;
   }
   else if( !strcmp( OPT->FRL_FILE, "range" ))
   {
     if( !OPT->BEG_FLAG || !OPT->END_FLAG )
     {
       fprintf( stderr, "-beg and -end are needed with range flag\n");
       exit(0);
     }
     frame_list.nframes = (OPT->TRJN_END - OPT->TRJN_BEG) +1;
     frame_list.nframe = malloc( frame_list.nframes * sizeof(int));
     jj = 0;
     for( ii=OPT->TRJN_BEG; ii<=OPT->TRJN_END; ii++ )
     {
       frame_list.nframe[jj] = ii;
       jj++;
     }
   }
   else
   {
     ReadFramesList ( OPT->FRL_FILE, &frame_list);
     // check whether the frame list is ordered or not
     sortedlist = checkSorted( &frame_list );
   }
   
   if(OPT->SELE_FLAG)
    {
     molecule = ReadMolecule ( OPT->IMOL_FILE, OPT->IMOL_TYPE );
     if( molecule->nato != traj1->hdr->nato )
     {
       fprintf( stderr, "# of atoms in molecule (%d) different from # in traj (%d)\n", molecule->nato, traj1->hdr->nato);
       exit(1);
     }
     GetSele( OPT->SELE_STRING, &selestr, molecule);
     traj2->hdr->nato = selestr.nselatm;
     if (selestr.nselatm == 0)
     {
      printf("Error: no atom selected\n");
      exit(95);
     }
     selcrd = InitCoor(selestr.nselatm);
     if( traj2->hdr->varpbc )
     {
       selcrd->pbc_flag = 1;
       selcrd->pbc = calloc(1, sizeof(Pbc));
     }
    }
   
   WriteTrjHeader ( traj2, "all" );
   
   trj_crd = InitCoor(traj1->hdr->nato);
   totpframes=0;
   kk=0;
   
   for( ii=0 ; ii < frame_list.nframes ; ii++)
   {
     ll=frame_list.nframe[ii];
     if( wrd_isdattxt(OPT->ITRJ_FILE) )
     {
       if( ! sortedlist )
       {
         totpframes=0;
         kk=0;
         CloseTrj( traj1 );
         traj1 = InitTrj( trj_list.trjname[0], "r" );
       }
       ll=frame_list.nframe[ii]-totpframes;
       while(ll > traj1->hdr->nframe)
       {
         totpframes += traj1->hdr->nframe;
         ll = frame_list.nframe[ii] - totpframes;
         kk++;
         CloseTrj( traj1 );
         if( kk == trj_list.ntrj )
         {
           fprintf( stderr, "[EE] total number of frames is %d but you asked for %d\n", totpframes, OPT->FRN_N);
           exit(99);
         }
         traj1 = InitTrj( trj_list.trjname[kk], "r" );
       }
     }
        
     ReadTrjCoor ( traj1, trj_crd, ll );
     
     if(OPT->SELE_FLAG)
     {
       for (jj=0; jj<selestr.nselatm; jj++)
       {
         selcrd->xcoor[jj] = trj_crd->xcoor[selestr.selatm[jj]-1];
         selcrd->ycoor[jj] = trj_crd->ycoor[selestr.selatm[jj]-1];
         selcrd->zcoor[jj] = trj_crd->zcoor[selestr.selatm[jj]-1];
       }
       if( traj2->hdr->varpbc )
       {
         selcrd->pbc->a_size = trj_crd->pbc->a_size;
         selcrd->pbc->angle1 = trj_crd->pbc->angle1;
         selcrd->pbc->b_size = trj_crd->pbc->b_size;
         selcrd->pbc->angle2 = trj_crd->pbc->angle2;
         selcrd->pbc->angle3 = trj_crd->pbc->angle3;
         selcrd->pbc->c_size = trj_crd->pbc->c_size;
       }
       WriteTrjCoor ( selcrd, traj2 );
     }
     else
       WriteTrjCoor ( trj_crd, traj2 );
     
   }
   
   DelCoor( trj_crd);
   if(OPT->SELE_FLAG)
     DelCoor(selcrd);
   
   traj2->hdr->nframe = frame_list.nframes;
   traj2->hdr->numstp = frame_list.nframes * traj2->hdr->skpstp;
   WriteTrjHeader ( traj2, "nframe" );
   WriteTrjHeader ( traj2, "numstp" );
   
   if(OPT->SELE_FLAG)
     WriteTrjHeader ( traj2, "nato" );
   
   CloseTrj( traj1 );
   CloseTrj( traj2 );

   return;
}
// ------------------------------------------------------------------
void MolExtract ( struct sopt   *OPT    )
{
   Traj        *traj;
   CoorSet     *trj_crd;
   TrjList      trj_list;
   Molecule    *molecule, *selmol;
   Selection    selestr;
   
   int          ii, jj, kk, ll, totpframes;
   

   molecule = ReadMolecule ( OPT->IMOL_FILE, OPT->IMOL_TYPE );
   
/*   traj.hdetail = traj.filetype==xtc ? LITE : FULL; // in case trj is xtc, we only need "lite" header*/
   
   if(OPT->SELE_FLAG)
   {
     GetSele( OPT->SELE_STRING, &selestr, molecule);
     if (selestr.nselatm == 0)
     {
       fprintf( stderr, "Error: no atom selected\n");
       exit(95);
     }
     selmol = MakeSeleMol ( &selestr , molecule);
   }

   if( wrd_isdattxt(OPT->ITRJ_FILE) )
   {
     ReadTrjList( OPT->ITRJ_FILE, &trj_list);
     traj = InitTrj( trj_list.trjname[0], "r" );
   }
   else
     traj = InitTrj( OPT->ITRJ_FILE, "r" );

   if ( traj->hdr->nato != molecule->nato)
   {
     fprintf( stderr, "ERROR:\n Number of atoms in trj differs from that in mol (%d vs %d)\n", traj->hdr->nato, molecule->nato);
     exit(98);
   }
   
   totpframes=0;
   kk=0;

   ll=OPT->FRN_N-totpframes;
   
   trj_crd = InitCoor( traj->hdr->nato);
   
   if( wrd_isdattxt(OPT->ITRJ_FILE) )
   {
     while(ll > traj->hdr->nframe)
     {
       totpframes += traj->hdr->nframe;
       ll = OPT->FRN_N - totpframes;
       kk ++;
       CloseTrj ( traj );
       if( kk == trj_list.ntrj )
       {
         fprintf( stderr, " [EE] total number of frames is %d but you asked for %d\n", totpframes, OPT->FRN_N);
         exit(99);
       }
       traj = InitTrj( trj_list.trjname[kk], "r" );
     }
   }
   
   ReadTrjCoor ( traj, trj_crd, ll );
   
   for ( ii=0; ii< 3*molecule->nato; ii++)
     molecule->coor.cords[ii] = trj_crd->cords[ii];

   if(OPT->SELE_FLAG)
   {
     for ( jj=0; jj<selmol->nato; jj++)
     {
       selmol->coor.xcoor[jj] = molecule->coor.xcoor[selestr.selatm[jj]-1];
       selmol->coor.ycoor[jj] = molecule->coor.ycoor[selestr.selatm[jj]-1];
       selmol->coor.zcoor[jj] = molecule->coor.zcoor[selestr.selatm[jj]-1];
     }
   }

   if(OPT->SELE_FLAG)
     WriteMolecule_unstr ( selmol, OPT->OMOL_FILE, OPT->OMOL_TYPE );
   else
     WriteMolecule_unstr ( molecule, OPT->OMOL_FILE, OPT->OMOL_TYPE );

   DelCoor( trj_crd );

   CloseTrj ( traj );
   
   return;
}
// ------------------------------------------------------------------
void MultiMolExtract ( struct sopt   *OPT    )
{
   Traj        *traj1, *tmptraj;
   CoorSet     *trj_crd;
   TrjList      trj_list;
   Molecule    *molecule, *selmol;
   IntList      frame_list;
   Selection    selestr;
   
   int          ii, jj, kk, ll, totpframes, totframesintrjs;
   char         buffer[128];
   char         basename[128];
   char        *loc;
   int          baselength;
   int          sortedlist=0;
   
   molecule = ReadMolecule ( OPT->IMOL_FILE, OPT->IMOL_TYPE );

/*   traj.hdetail = traj.filetype==xtc ? LITE : FULL; // in case trj is xtc, we only need "lite" header */
   
   if( wrd_isdattxt(OPT->ITRJ_FILE) )
   {
     ReadTrjList    ( OPT->ITRJ_FILE, &trj_list);
     traj1 = InitTrj ( trj_list.trjname[0], "r" );
     totframesintrjs = traj1->hdr->nframe;
     for( ii=1; ii<trj_list.ntrj; ii++)
     {  
       tmptraj = InitTrj( trj_list.trjname[ii], "r" );
       totframesintrjs += tmptraj->hdr->nframe;
       CloseTrj( tmptraj );
     }
   }
   else
   {
     traj1 = InitTrj ( OPT->ITRJ_FILE, "r" );
     totframesintrjs = traj1->hdr->nframe;
   }

   if( !strncmp( OPT->FRL_FILE, "all", 3 ))
   {
     frame_list.nframes = totframesintrjs;
     frame_list.nframe = malloc( frame_list.nframes * sizeof(int));
     for( ii=0; ii<frame_list.nframes; ii++ )
       frame_list.nframe[ii] = ii+1;
   }
   else if( !strcmp( OPT->FRL_FILE, "range" ))
   {
     if( !OPT->BEG_FLAG || !OPT->END_FLAG )
     {
       fprintf( stderr, "-beg and -end are needed with range flag\n");
       exit(0);
     }
     frame_list.nframes = (OPT->TRJN_END - OPT->TRJN_BEG) +1;
     frame_list.nframe = malloc( frame_list.nframes * sizeof(int));
     jj = 0;
     for( ii=OPT->TRJN_BEG; ii<=OPT->TRJN_END; ii++ )
     {
       frame_list.nframe[jj] = ii;
       jj++;
     }
   }
   else
   {
     ReadFramesList ( OPT->FRL_FILE, &frame_list);
     // check whether the frame list is ordered or not
     sortedlist = checkSorted( &frame_list );
   }
   
   if(OPT->SELE_FLAG)
   {
     GetSele( OPT->SELE_STRING, &selestr, molecule);
     if (selestr.nselatm == 0)
     {
       fprintf( stderr, "Error: no atom selected\n");
       exit(95);
     }
     selmol = MakeSeleMol ( &selestr , molecule);
   }

   if ( traj1->hdr->nato != molecule->nato)
   {
     fprintf( stderr, "ERROR:\n Number of atoms in trajectory(%d) differs from that in molecule(%d)!\n", traj1->hdr->nato, molecule->nato);
     exit(98);
   }
   
   memset( basename, '\0', 128 );
   loc = strrchr( OPT->OMOL_FILE, '.' );
   baselength = strlen(OPT->OMOL_FILE ) - strlen( loc );
   for( jj=0; jj<baselength; jj++)
     basename[jj] = OPT->OMOL_FILE[jj];
   
   trj_crd = InitCoor( traj1->hdr->nato );
   totpframes=0;
   kk=0;
   
   for( ii=0; ii<frame_list.nframes; ii++ )
   {
     ll=frame_list.nframe[ii];
     if( wrd_isdattxt(OPT->ITRJ_FILE) )
     {
       if( ! sortedlist )
       {
         totpframes=0;
         kk=0;
         CloseTrj( traj1 );
         traj1 = InitTrj( trj_list.trjname[0], "r" );
       }
       ll=frame_list.nframe[ii]-totpframes;
       while( ll > traj1->hdr->nframe )
       {
         totpframes += traj1->hdr->nframe;
         ll = frame_list.nframe[ii] - totpframes;
         kk++;
         CloseTrj ( traj1 );
         if( kk == trj_list.ntrj )
         {
           fprintf( stderr, "[EE] total number of frames is %d but you asked for %d\n", totpframes, OPT->FRN_N);
           exit(99);
         }
         traj1 = InitTrj( trj_list.trjname[kk], "r" );
       }
     }
   
     ReadTrjCoor ( traj1, trj_crd, ll );
   
     for ( jj=0; jj<3*molecule->nato; jj++)
       molecule->coor.cords[jj] = trj_crd->cords[jj];
     
     if(OPT->SELE_FLAG)
     {
       for ( jj=0; jj<selmol->nato; jj++)
       {
         selmol->coor.xcoor[jj] = molecule->coor.xcoor[selestr.selatm[jj]-1];
         selmol->coor.ycoor[jj] = molecule->coor.ycoor[selestr.selatm[jj]-1];
         selmol->coor.zcoor[jj] = molecule->coor.zcoor[selestr.selatm[jj]-1];
       }
     }

     memset( buffer, '\0', 128 );
     if(OPT->OMOL_TYPE == PDB)
       sprintf( buffer, "%s_%d.pdb", basename, frame_list.nframe[ii]);
     else if(OPT->OMOL_TYPE == CRD)
       sprintf( buffer, "%s_%d.crd", basename, frame_list.nframe[ii]);
     
     if(OPT->SELE_FLAG)
       WriteMolecule_unstr ( selmol, buffer, OPT->OMOL_TYPE );
     else
       WriteMolecule_unstr ( molecule, buffer, OPT->OMOL_TYPE );
     
     if( wrd_isdattxt(OPT->ITRJ_FILE) )
     {
       CloseTrj ( traj1 );
       traj1 = InitTrj( trj_list.trjname[0], "r" );
     }
   }

   DelCoor( trj_crd );

   CloseTrj ( traj1 );
   
   return;
}
// ------------------------------------------------------------------
void AvgExtract ( struct sopt   *OPT    )
{
   Molecule    *molecule;
   Traj        *traj, *tmptraj;
   CoorSet     *trj_crd;
   TrjList      trj_list;
   int          ii, jj, kk, hh;
   int          begframe, endframe, skip=1, nframe, begcount, usedframes;

   molecule = ReadMolecule ( OPT->IMOL_FILE, OPT->IMOL_TYPE );
   
   for ( ii=0; ii<3*molecule->nato; ii++ )
     molecule->coor.cords[ii] = 0;
   
   begframe = 0;
   
   trj_list.ntrj=1;
   if( wrd_isdattxt(OPT->ITRJ_FILE) )
   {
     ReadTrjList    ( OPT->ITRJ_FILE, &trj_list);
     traj = InitTrj( trj_list.trjname[0], "r" );
     endframe = traj->hdr->nframe;
     for( ii=1; ii<trj_list.ntrj; ii++)
     {  
       tmptraj = InitTrj( trj_list.trjname[ii], "r" );
       endframe += tmptraj->hdr->nframe;
       CloseTrj( tmptraj );
     }
   }
   else
   {
     traj = InitTrj( OPT->ITRJ_FILE, "r" );
     endframe = traj->hdr->nframe;
   }

   if ( traj->hdr->nato != molecule->nato)
   {
     fprintf( stderr, "ERROR:\n Number of atoms in trj differs from that in mol!\n");
     exit(98);
   }
   
   if(OPT->BEG_FLAG)
   {
    if(OPT->TRJN_BEG >= 0)
      begframe = OPT->TRJN_BEG;
    else
    {
      fprintf( stderr, "Starting frame number must be >= 0\n");
      exit(0);
    }
   }

   if(OPT->END_FLAG)
   {
    if(OPT->TRJN_END <= endframe)
      endframe = OPT->TRJN_END;
    else
    {
      fprintf( stderr, "Last frame must be <= %d (total number of frames)\n", traj->hdr->nframe);
      exit(0);
    }
   }
   
   if(OPT->SKIP_FLAG)
     skip = OPT->SKIP_STEP;
    
   hh=1;
   nframe=0;
   usedframes = 0;
   if(hh<trj_list.ntrj)
     hh = trj_list.ntrj;
   
   trj_crd = InitCoor( traj->hdr->nato );
   for ( kk=0; kk<hh; kk++ )
   {
    if( kk>0 && wrd_isdattxt(OPT->ITRJ_FILE) )
    {
     CloseTrj ( traj );
     traj = InitTrj( trj_list.trjname[kk], "r" );
    }
    if( kk == 0 )
    {
      begcount = begframe;
      nframe   = begframe;
    }
    else
      begcount=0;
    for( ii=begcount; ii < traj->hdr->nframe && nframe < endframe; ii += skip)
    {
      ReadTrjCoor ( traj, trj_crd, ii+1 );
     
      for ( jj=0; jj<3*traj->hdr->nato; jj++)
        molecule->coor.cords[jj] += trj_crd->cords[jj];
     
      usedframes += 1;
      nframe+=skip;                
    }
   }
   DelCoor( trj_crd );
   
   for ( jj=0; jj<3*traj->hdr->nato; jj++)
     molecule->coor.cords[jj] /= usedframes;

   WriteMolecule_unstr ( molecule, OPT->OMOL_FILE, OPT->OMOL_TYPE );
   
   CloseTrj ( traj );
   
   return;
}
// ------------------------------------------------------------------
void TrjHeader_Print ( struct sopt *OPT ) 
{
   Traj        *traj;
   
   traj = InitTrj( OPT->ITRJ_FILE, "r" );
   traj->hdetail = FULL;
   ShowTrjHeader ( traj );
   CloseTrj ( traj );
   
   return;
}        
// ------------------------------------------------------------------
void MolInfo_Print( char *filename )
{
  Molecule *mol1;
  
  mol1 = ReadMolecule( filename, wrd_whichmoltype( filename ));
  ShowMolInfo( mol1 );
  
}
// ------------------------------------------------------------------
void xyzExtract ( struct sopt   *OPT )
{
   FILE         *oxyz_f;   // xyz output file
   Traj         *traj;
   CoorSet      *trj_crd;   // coordinates
   Selection     sele;      // selection structure
   int           aa;        // counter
   int          *numatms;   // numbers of atoms to extract
   int           nnumb=0;
   char          buffer[1024];
   
   Molecule     *molecule;
   int           ii;
   
   FILE         *list_f;
   
   molecule = ReadMolecule ( OPT->IMOL_FILE, OPT->IMOL_TYPE );
   
/*   traj.hdetail = FULL;*/
   
   if( !OPT->SELE_FLAG)
   {
    OPT->SELE_STRING = malloc ( 10*sizeof( char) );
    //memset( OPT->SELE_STRING, '\0', sizeof (OPT->SELE_STRING) );
    memset( OPT->SELE_STRING, '\0', 10 );
    sprintf( OPT->SELE_STRING, "/*/*/*");
   }
   GetSele( OPT->SELE_STRING, &sele, molecule);

   traj = InitTrj( OPT->ITRJ_FILE, "r" );
   if( traj->filetype == xtc )
     ReadTrjHeader ( traj );
   
   if ( traj->hdr->nato != molecule->nato)
   {
     fprintf( stderr, "ERROR:\n Number of atoms in trj differs from that in mol!\n");
     exit(98);
   }
   
   if( OPT->SKIP_FLAG == 0 ) 
     OPT->SKIP_STEP = 1;
   
   if ( OPT->BEG_FLAG == 0 )
   {
     if ( OPT->SKIP_FLAG == 0 )
       OPT->TRJN_BEG = 1;
     if ( OPT->SKIP_FLAG == 1 )
       OPT->TRJN_BEG = OPT->SKIP_STEP;
   }

   if ( OPT->END_FLAG == 0 )
     OPT->TRJN_END = traj->hdr->nframe;
   
   if ( OPT->ATMNUM_FLAG )
   {
     list_f = O_File ( OPT->ATMNUM_FILE, "r" );

     while ( !feof(list_f))
     {
       ++nnumb;
       if(fgets(buffer, 90, list_f)==NULL && ( !feof(list_f) || ferror(list_f) ))
       {
         fprintf(stderr, "Warning! Premature end of file reached!\n");
       }
     }
     rewind(list_f);

     numatms = calloc ( nnumb, sizeof (int));
       
     for ( ii = 0; ii < nnumb; ii++ )
       fscanf ( list_f, "%d", &numatms[ii]);
   }
   
   trj_crd = InitCoor( traj->hdr->nato );
   if( traj->hdr->varpbc )
   {
     trj_crd->pbc_flag = 1;
     trj_crd->pbc = calloc( 1, sizeof(Pbc));
   }
   
   oxyz_f = O_File ( OPT->OMOL_FILE, "w" );
   for ( aa=OPT->TRJN_BEG; aa<=OPT->TRJN_END; aa=aa+OPT->SKIP_STEP )
   {
     ReadTrjCoor ( traj, trj_crd, aa );
     fprintf ( oxyz_f, "XYZ %d\n", aa );
     if( trj_crd->pbc == NULL )
       fprintf( oxyz_f, "pbc: none\n" );
     else
       fprintf( oxyz_f, "(pbc: %f %f %f)\n", trj_crd->pbc->angle1, trj_crd->pbc->angle2, trj_crd->pbc->angle3 );
     for ( ii=0; ii< sele.nselatm; ii++)
       fprintf ( oxyz_f, "%8.3f %8.3f %8.3f\n", trj_crd->xcoor[sele.selatm[ii]-1], trj_crd->ycoor[sele.selatm[ii]-1], trj_crd->zcoor[sele.selatm[ii]-1] );
   }
   DelCoor( trj_crd );

   CloseTrj ( traj );
   fclose   ( oxyz_f );

   return;
}
// ------------------------------------------------------------------
void CheckSele ( struct sopt   *OPT )
{
   Selection     sele;      // selection structure
   Molecule      molecule;
   int           ii;
   
   GetMolecule ( OPT->IMOL_FILE, &molecule, OPT->IMOL_TYPE );
   if( OPT->NOPBC_FLAG == 1 )
     molecule.coor.pbc_flag = 0;

   GetSele( OPT->SELE_STRING, &sele, &molecule);
   
   fprintf( stdout, "#N of selected atoms (%s): %d\n", sele.selestring, sele.nselatm);
   for( ii=0; ii<sele.nselatm; ii++ )
    fprintf( stdout, "%d\n", sele.selatm[ii]);
   
   return;
}
// ------------------------------------------------------------------
void mkResList ( struct sopt   *OPT )
{
  Molecule      molecule;
  int           ii, jj, ll, mm;
  FILE          *resfile1;
  FILE          *resfile2;
  
  GetMolecule ( OPT->IMOL_FILE, &molecule, OPT->IMOL_TYPE );
  resfile1 = O_File( "reslist1.txt", "w" );
  resfile2 = O_File( "reslist2.txt", "w" );
  
  for( ii=0; ii< molecule.nSeg; ii++ )
  {
    for( jj=0; jj<molecule.segment[ii].nRpS; jj++ )
    {
      for( ll=0; ll< molecule.nSeg; ll++ )
      {
        for( mm=0; mm<molecule.segment[ll].nRpS; mm++ )
        {
          if( ii!=ll || jj!=mm )
          {
            fprintf( resfile1, "%s:%d\n", molecule.segment[ii].segName, molecule.segment[ii].pRes[jj].resn );
            fprintf( resfile2, "%s:%d\n", molecule.segment[ll].segName, molecule.segment[ll].pRes[mm].resn );
          }
        }
      }
    }
  }
  
  fclose( resfile1 );
  fclose( resfile2 );
  
  return;
}
// ===================================================================
