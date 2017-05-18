// ------------------------------------------------------------------
// Copyright (C) 2009  University of Modena and Reggio Emilia and
//                     University of Zurich
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
/*! \file analysis.c
 \brief 
 
 
*/
#include <ctype.h>
#include "wordom.h"
#include "fileio.h"
#include "tools.h"
#include "datahandler.h"
#include "analysis.h"

extern short int verbose;

// ------------------------------------------------------------------
int CoorSetInit(sopt *OPT, CoorSetData *coorsetdata )
{
   int          ii;
   int          nsets, tmp_nato = 0;
   Traj         traj;
   Molecule     molecule;
   
   coorsetdata->nframe = 0;
   coorsetdata->begframe = 0;
   coorsetdata->endframe = 0;
   coorsetdata->skip = 1;
   coorsetdata->issinglestr = 0;
   coorsetdata->issingletraj = 0;
   coorsetdata->islistoftrajs = 0;
   coorsetdata->islistofstrs = 0;
   coorsetdata->isfrlist = 0;
   
   coorsetdata->datasetname = calloc( strlen(OPT->ITRJ_FILE)+1, sizeof(char));
   sprintf( coorsetdata->datasetname, "%s", OPT->ITRJ_FILE);
   
   molecule.ext_flag = 0;
   molecule.coor.pbc = NULL;
   molecule.coor.pbc_flag = 0;
   molecule.filled = calloc( 4, sizeof(char));
   memset( molecule.filled, '\0', 4);
   
   if(OPT->ITRJ_TYPE == txt)
   {
     ReadTrjList    ( OPT->ITRJ_FILE, &coorsetdata->trj_list);
     
     for( ii=1; ii<coorsetdata->trj_list.ntrj; ii++)
       if( wrd_whichtrjtype(coorsetdata->trj_list.trjname[ii-1]) != wrd_whichtrjtype(coorsetdata->trj_list.trjname[ii]))
       {
         fprintf( stderr, "itrj list must be homogeneous\n");
         exit(98);
       }
     
     if( wrd_whichtrjtype(coorsetdata->trj_list.trjname[0]) == dcd ||
         wrd_whichtrjtype(coorsetdata->trj_list.trjname[0]) == xtc    ) 
     {
       for( ii=0; ii<coorsetdata->trj_list.ntrj; ii++ )
       {
         OpenTrj( coorsetdata->trj_list.trjname[ii], &traj, "r" );
         ReadTrjHeader ( &traj );
         coorsetdata->nframe += traj.hdr->nframe;
         if( tmp_nato != 0 && tmp_nato != traj.hdr->nato )
         {
           fprintf( stderr, "Error in trj list: trajectories do not have the same number of atoms\n");
           exit(0);
         }
         tmp_nato = traj.hdr->nato;
         CloseTrj ( &traj );
       }
       OpenTrj( coorsetdata->trj_list.trjname[0], &coorsetdata->traj, "r" );
       ReadTrjHeader ( &coorsetdata->traj );
       coorsetdata->nato = coorsetdata->traj.hdr->nato;
       coorsetdata->islistoftrajs = 1;
       coorsetdata->frbefthistrj = 0;
       coorsetdata->currenttrj = 0;
       nsets = coorsetdata->nframe;
     }
     else if( wrd_whichtrjtype(coorsetdata->trj_list.trjname[0]) == pdb ||
              wrd_whichtrjtype(coorsetdata->trj_list.trjname[0]) == crd ||
              wrd_whichtrjtype(coorsetdata->trj_list.trjname[0]) == cor    )
     {
       GetMolecule( coorsetdata->trj_list.trjname[0], &molecule,  wrd_whichmoltype(coorsetdata->trj_list.trjname[0]) );
       coorsetdata->traj.hdr = malloc(sizeof(struct trjh));
       FillDcdHeader( &coorsetdata->molecule, coorsetdata->traj.hdr);
       coorsetdata->nframe = coorsetdata->trj_list.ntrj;
       coorsetdata->nato = molecule.nato;
       coorsetdata->islistofstrs = 1;
       nsets = coorsetdata->trj_list.ntrj;
     }
     else
     {
       fprintf( stderr, "Did not recognize format of files in itrj list\n");
       exit(0);
     }
   }
   else
   { 
     if( wrd_whichtrjtype( OPT->ITRJ_FILE ) == dcd ||
         wrd_whichtrjtype( OPT->ITRJ_FILE ) == xtc   ) 
     {
       OpenTrj( OPT->ITRJ_FILE, &coorsetdata->traj, "r" );
       ReadTrjHeader ( &coorsetdata->traj );
       coorsetdata->nframe = coorsetdata->traj.hdr->nframe;
       coorsetdata->nato = coorsetdata->traj.hdr->nato;
       coorsetdata->issingletraj = 1;
       nsets = coorsetdata->nframe;
     }
     if( wrd_whichtrjtype( OPT->ITRJ_FILE ) == pdb ||
         wrd_whichtrjtype( OPT->ITRJ_FILE ) == crd || 
         wrd_whichtrjtype( OPT->ITRJ_FILE ) == cor   ) 
     {
       GetMolecule( OPT->ITRJ_FILE, &coorsetdata->molecule,  wrd_whichmoltype(OPT->ITRJ_FILE) );
       coorsetdata->traj.hdr = malloc(sizeof(struct trjh));
       FillDcdHeader( &coorsetdata->molecule, coorsetdata->traj.hdr);
       coorsetdata->nframe = 1;
       coorsetdata->nato = coorsetdata->molecule.nato;
       coorsetdata->issinglestr = 1;
       nsets = 1;
     }
   }
   
   if( OPT->FRL_FLAG )
   {
     coorsetdata->isfrlist = 1;
     ReadFramesList( OPT->FRL_FILE, &coorsetdata->frame_list );
     if( coorsetdata->nframe < coorsetdata->frame_list.nframe[coorsetdata->frame_list.nframes-1] )
     {
       fprintf( stderr, "Frame list exceeds available frames in supplied trj(s)\n");
       exit(0);
     }
     coorsetdata->nframe = coorsetdata->frame_list.nframes;
     nsets = coorsetdata->nframe;
   }
   
   return nsets;
}
// ---------------------------------------------------------------------
int CheckRanges(sopt *OPT, Molecule *molecule, struct inp_A *iA_data, CoorSetData *coorsetdata)
{
   
   if ( coorsetdata->nato != molecule->nato)
   {
     fprintf(stderr, "ERROR:\n Number of atoms in trj (%d) differs from that in mol(%d)!\n", coorsetdata->nato, molecule->nato);
     exit(98);
   }
   
   coorsetdata->begframe = 1;
   coorsetdata->endframe = coorsetdata->nframe;
   
   if(OPT->BEG_FLAG)
   {
     if(OPT->TRJN_BEG > 0)
      coorsetdata->begframe = OPT->TRJN_BEG;
     else
     {
      fprintf(stderr, "Starting frame number must be > 0\n");
      exit(0);
     }
   }

   if(OPT->END_FLAG)
   {
    if( OPT->TRJN_END <= coorsetdata->nframe )
     coorsetdata->endframe = OPT->TRJN_END;
    else
    {
     fprintf(stderr, "Last frame must be <= %d (total number of frames)\n", coorsetdata->nframe);
     exit(0);
    }
   }
   
   if(OPT->SKIP_FLAG)
    coorsetdata->skip = OPT->SKIP_STEP;
   
   coorsetdata->i_nUsedFrames = ( coorsetdata->endframe - coorsetdata->begframe + 1)/ coorsetdata->skip;
   
   if( coorsetdata->isfrlist && ( OPT->SKIP_FLAG || OPT->TRJN_END || OPT->TRJN_BEG ))
   {
     fprintf( stderr, "Com\'on, don't use frame list with beg/end/skip. Edit your file\n");
     exit(0);
   }
   
   if( coorsetdata->issinglestr )
     coorsetdata->endframe = 1;
   
   return 1;
}
// ---------------------------------------------------------------------
int GetThisCoorSet( CoorSetData *coorsetdata, int number, CoorSet *crdset )
{
   int          ii;
   int          intex=0, found=0;
   
   if( coorsetdata->isfrlist )
     number = coorsetdata->frame_list.nframe[number];
   
   if( coorsetdata->issingletraj )
   {
     intex = number + (coorsetdata->begframe-1);
     ReadTrjCoor( &coorsetdata->traj, crdset, intex+1 );
     return intex +1;
   }
   else if( coorsetdata->islistoftrajs )
   {
     intex = number + (coorsetdata->begframe-1);
     if( intex < (coorsetdata->frbefthistrj + coorsetdata->traj.hdr->nframe) )
     {
       ReadTrjCoor( &coorsetdata->traj, crdset, (intex - coorsetdata->frbefthistrj)+1 );
     }
     else
     {
       found = 0;
       while( !found )
       {
         coorsetdata->frbefthistrj += coorsetdata->traj.hdr->nframe;
         CloseTrj( &coorsetdata->traj );
         coorsetdata->currenttrj ++;
         if( coorsetdata->currenttrj > coorsetdata->trj_list.ntrj )
         {
           fprintf( stderr, "Max number of available frames exceeded while looking for frame #%d\n", intex);
           exit(0);
         }
         OpenTrj( coorsetdata->trj_list.trjname[coorsetdata->currenttrj], &coorsetdata->traj, "r" );
         ReadTrjHeader ( &coorsetdata->traj );
         if( intex < (coorsetdata->frbefthistrj + coorsetdata->traj.hdr->nframe) )
         {
           ReadTrjCoor( &coorsetdata->traj, crdset, (intex - coorsetdata->frbefthistrj)+1 );
           found = 1;
         }
       }
     }
     return intex+1;
   }
   else if( coorsetdata->issinglestr )
   {
     for( ii=0; ii<coorsetdata->nato; ii++ )
     {
       crdset->xcoor[ii] = coorsetdata->molecule.coor.xcoor[ii];
       crdset->ycoor[ii] = coorsetdata->molecule.coor.ycoor[ii];
       crdset->zcoor[ii] = coorsetdata->molecule.coor.zcoor[ii];
     }
     intex = 0;
     return intex+1;
   }
   else if( coorsetdata->islistofstrs)
   {
     intex = number + (coorsetdata->begframe-1);
     GetMolecule( coorsetdata->trj_list.trjname[intex], &coorsetdata->molecule, wrd_whichmoltype(coorsetdata->trj_list.trjname[intex]));
     for( ii=0; ii<coorsetdata->nato; ii++ )
     {
       crdset->xcoor[ii] = coorsetdata->molecule.coor.xcoor[ii];
       crdset->ycoor[ii] = coorsetdata->molecule.coor.ycoor[ii];
       crdset->zcoor[ii] = coorsetdata->molecule.coor.zcoor[ii];
     }
     CleanMolecule( &coorsetdata->molecule );
     return intex+1;
   }
   
   fprintf( stderr, "warning: mess in GetThisCoorSet function\n");
   return intex+1;
}
// ---------------------------------------------------------------------
int CoorSetRewind( CoorSetData *coorsetdata )
{
  int       ii;
  char      tmptrjname[1024];
  
  if( coorsetdata->islistoftrajs )
  {
    coorsetdata->frbefthistrj = 0;
    coorsetdata->currenttrj = 0;
    CloseTrj( &coorsetdata->traj );
    OpenTrj( coorsetdata->trj_list.trjname[coorsetdata->currenttrj], &coorsetdata->traj, "r" );
    ReadTrjHeader ( &coorsetdata->traj );
  }
  else if( coorsetdata->traj.filetype == 2) {
    for( ii=0; ii<1024; ii++ )
      tmptrjname[ii] = coorsetdata->traj.filename[ii];
    CloseTrj( &coorsetdata->traj );
    OpenTrj( tmptrjname, &coorsetdata->traj, "r" );
  }
  
  return 0;
}
// ---------------------------------------------------------------------
int GetThisIndex( CoorSetData *coorsetdata, int number )
{
   int          intex;
   
   if( coorsetdata->isfrlist )
     number = coorsetdata->frame_list.nframe[number];
   
   if( coorsetdata->issingletraj )
     intex = number + (coorsetdata->begframe-1);
   else if( coorsetdata->islistoftrajs )
     intex = number + (coorsetdata->begframe-1);
   else if( coorsetdata->issinglestr )
     intex = 0;
   else if( coorsetdata->islistofstrs)
     intex = number + (coorsetdata->begframe-1);

   return intex+1;
}
// ---------------------------------------------------------------------
void superCalc ( struct sopt   *OPT  )
{
   int                  *piExitStatus;
   struct threading     *T_data;
   Molecule              molecule;
   CoorSetData          coorsetdata;
   CoorSet              trj_crd;
   int                   ii, jj, TT, papersize;
   FILE                 *oA_f;
   char                 *header;
   
   extern short int     frame_par;
   extern short int     module_par;
   extern short int     no_frame_par;
   
  #ifdef THREADED
   pthread_t            *pThreads;
   molecule.ext_flag = 0;
   molecule.coor.pbc = NULL;
   molecule.coor.pbc_flag = 0;
   molecule.filled = calloc( 4, sizeof(char));
   memset( molecule.filled, '\0', 4);
   if( OPT->THR_FLAG )
   {
     // Open Output File
     if ( OPT->OA_FLAG )
       oA_f = O_File ( OPT->OA_FILE, "w" );
     else
       oA_f = stdout;
     
     // get molecule for everybody
     GetMolecule ( OPT->IMOL_FILE, &molecule, OPT->IMOL_TYPE );

     // prepare threading structures
     pThreads = (pthread_t *) malloc ( OPT->THR_NUM * sizeof(pthread_t));
     piExitStatus = (int *) malloc ( OPT->THR_NUM * sizeof(int));
     T_data = malloc( OPT->THR_NUM*sizeof(struct threading));
     
     // and data structure for each thread
     T_data[0].iA_data.molecule = &molecule;
     
     // coorsetdata and input read for first thread only (the others will do by themselves)
     T_data[0].coorsetdata = calloc( 1, sizeof(CoorSetData));
     T_data[0].coorsetdata->traj.hdetail = FULL;
     CoorSetInit( OPT, T_data[0].coorsetdata );
     CheckRanges( OPT, &molecule, &T_data[0].iA_data, T_data[0].coorsetdata) ;
     
     T_data[0].iA_data.nframe = T_data[0].coorsetdata->i_nUsedFrames;
     
     fprintf( oA_f, "#Fr     ");
     header = Read_iA ( OPT, &T_data[0].iA_data, &papersize, &molecule, &coorsetdata);
     fprintf(oA_f,"%s\n", header);
     free(header);
     
     // sanity check: can threading be really used/useful ?
     if( no_frame_par == 1 && module_par == 1 )
     {
       fprintf( stderr, "Framewise threading not possible with one of the modules you selected\n");
       fprintf( stderr, "Use the --NC option in the input to use module-wise threading\n");
       exit(0);
     }
     if( no_frame_par == 1 && module_par == 0 )
     {
       fprintf( stderr, "Threading not possible with one of the modules you selected\n");
       fprintf( stderr, "\n");
       exit(0);
     }
     else if( frame_par == 0 )
     {
       fprintf( stderr, "Framewise threading not available with any of the modules you selected\n");
       fprintf( stderr, "Wordom would not run faster then in single-thread mode anyway\n");
       exit(0);
     }
     
     T_data[0].printout = calloc( T_data[0].coorsetdata->nframe+1, sizeof(char *));
     for( ii=0; ii<T_data[0].coorsetdata->nframe+1; ii++)
       T_data[0].printout[ii] = calloc( papersize, sizeof(char));
     
     for( TT=0; TT<OPT->THR_NUM; TT++)
     {
       T_data[TT].iA_data.threaded = 1;
       T_data[TT].iA_data.thread_rank = TT;
       T_data[TT].iA_data.nthreads = OPT->THR_NUM;
       T_data[TT].iA_data.opt = OPT;
       T_data[TT].thread_size = OPT->THR_NUM;
       T_data[TT].thread_rank = TT;
       T_data[TT].printout = T_data[0].printout;
       piExitStatus[TT] = pthread_create( &pThreads[TT], NULL, iA_Calc_thread, (void *) &T_data[TT] );
     }
     
     for(TT=0; TT<OPT->THR_NUM; TT++)
       pthread_join(pThreads[TT], NULL);
     
     // when compute itself is over, finally print all results to file
     for( ii=0; ii<T_data[0].iA_data.nframe; ii++)
       fprintf( oA_f, "%7d %s\n", ii+1, T_data[0].printout[ii]);

   /* Now some post-processing for those who need it (hard to do for threaded stuff - take care) */
     for (jj=0; jj<T_data[0].iA_data.nBEG; jj++)
     {
      switch (T_data[0].iA_data.A_type[jj])
      {
       case   9 :
        if (T_data[0].iA_data.input[jj]->inp_rms.trjwrite)
          WriteTrjHeader ( T_data[0].iA_data.input[jj]->inp_rms.outtraj, "nframe");
        break;
       case  10 :
        Post_Rmsf ( &T_data[0].iA_data.input[jj]->inp_rmsf, OPT, T_data[0].iA_data.nframe, &molecule );
        break;
       case  11 :
        Post_PCA ( &T_data[0].iA_data.input[jj]->inp_pca, OPT, T_data[0].iA_data.nframe, &molecule );
        break;
       case  12 :
        Post_pro ( &T_data[0].iA_data.input[jj]->inp_pro, OPT, T_data[0].iA_data.nframe, &molecule );
        break;
       case  13 :
        Post_Cluster ( &T_data[0].iA_data.input[jj]->inp_cluster );
        break;
        #ifdef CUDA
       case  84 :
        Post_Gcluster ( &T_data[0].iA_data.input[jj]->inp_cluster );
        break;
        #endif
       case  15 :
        Post_Entropy ( &T_data[0].iA_data.input[jj]->inp_pca, T_data[0].iA_data.nframe, &molecule );
        break;
       case  18 :
        Post_PSG (&T_data[0].iA_data.input[jj]->inp_psn, T_data[0].iA_data.nframe, &molecule, OPT );
        break;
       case  19 : // WARNING: trj_crd never used in this scope !!! not working !!!
        Post_CORR ( &T_data[0].iA_data.input[jj]->inp_corr, &molecule, OPT, &trj_crd);
        break;
       case  24 :
        Post_Msdf ( &T_data[0].iA_data.input[jj]->inp_msdf, OPT, T_data[0].iA_data.nframe, &molecule );
        break;
       //case 99 :
       // Post_Sample ( &T_data[0].iA_data.input[jj]->inp_sample, OPT, T_data[0].iA_data.nframe, &molecule );
       // break;
       default :
        break;
      }
     }


     return;
   }
   
   // WARNING : ADD POST_STUFF HERE
  #endif

   fflush(stdout);
   iA_Calc( OPT );
   return;
}
// ===================================================================
void run_iA ( int ii, int intex, CoorSet *trj_crd, struct inp_A *iA_data, char *printout, Molecule *molecule)
{
   int                  jj, offset;
   
   ii = ii+1;
   offset = 0;
   for(jj=0; jj<iA_data->nBEG; jj++)
   {
     if( iA_data->bonoffmodulespost[jj] == 0 )
       continue;
     switch (iA_data->A_type[jj])
     {
       case  1:
        offset += Compute_Dist ( &iA_data->input[jj]->inp_dist, iA_data->opt, molecule, trj_crd, &printout[offset] );
        break;
       case  2:
        offset += Compute_Contacts ( &iA_data->input[jj]->inp_contacts, iA_data->opt, molecule, trj_crd, &printout[offset] );
        break;
       case  3:
        offset += Compute_Angle ( &iA_data->input[jj]->inp_angle, iA_data->opt, trj_crd, &printout[offset] );
        break;
       case  4:
        offset += Compute_Dihe ( &iA_data->input[jj]->inp_dihe, iA_data->opt, molecule, trj_crd, &printout[offset] );
        break;
       case  5:
        offset += Compute_Rgyr ( &iA_data->input[jj]->inp_rgyr, iA_data->opt, trj_crd, &printout[offset] );
        break;
       case  6:
        offset += Compute_hb ( &iA_data->input[jj]->inp_hb, iA_data->opt, trj_crd, &printout[offset] );
        break;
       case  7:
        offset += Compute_P ( &iA_data->input[jj]->inp_P, iA_data->opt, molecule, trj_crd, &printout[offset] );
        break;
       case  8:
        offset += Compute_drms ( &iA_data->input[jj]->inp_drms, iA_data->opt, trj_crd, &printout[offset] );
        break;
       case  9:
        offset += Compute_Rmsd ( &iA_data->input[jj]->inp_rms, iA_data->opt, trj_crd, &printout[offset] );
        break;
       case 10:
        offset += Compute_Rmsf ( &iA_data->input[jj]->inp_rmsf, iA_data->opt, trj_crd, &printout[offset] );
        break;
       case 11:
        offset += Compute_PCA ( &iA_data->input[jj]->inp_pca, iA_data->opt, trj_crd, &printout[offset], intex );
        break;
       case 12:
        offset += Compute_pro ( &iA_data->input[jj]->inp_pro, iA_data->opt, trj_crd, &printout[offset] );
        break;
       case 13:
        offset += Compute_Cluster ( &iA_data->input[jj]->inp_cluster, iA_data->opt, trj_crd, &printout[offset], intex );
        break;
        #ifdef CUDA
       case 84:
        offset += Compute_Gcluster ( &iA_data->input[jj]->inp_cluster, iA_data->opt, trj_crd, &printout[offset], intex );
        break;
        #endif
       case 14:
        offset += Compute_CAssign ( &iA_data->input[jj]->inp_cassign, iA_data->opt, trj_crd, &printout[offset], intex );
        break;
       case 15:
        offset += Compute_Entropy ( &iA_data->input[jj]->inp_pca, iA_data->opt, trj_crd,  &printout[offset]);
        break;
       case 16:
        offset += Compute_SSA ( &iA_data->input[jj]->inp_ss, iA_data->opt, molecule, trj_crd, &printout[offset] );
        break;
       case 17:
        offset += Compute_Enm( &iA_data->input[jj]->inp_enm, molecule, trj_crd, &printout[offset], intex);
        break;
       case 18:
        offset += Compute_PSG ( &iA_data->input[jj]->inp_psn, iA_data->opt, trj_crd, &printout[offset], molecule );
        break;
       case 19:
        offset += Compute_CORR ( &iA_data->input[jj]->inp_corr, molecule, trj_crd, &printout[offset] );
        break;
       case 20 :
        offset += Compute_SURFCLUST ( &iA_data->input[jj]->inp_surf, molecule, trj_crd, &printout[offset] );
        break;
       case 21 :
        offset += Compute_SURFCORR ( &iA_data->input[jj]->inp_surf, molecule, trj_crd, &printout[offset] );
        break;
       case 22 :
        offset += Compute_SURF ( &iA_data->input[jj]->inp_surf, molecule, trj_crd, &printout[offset] );
        break;
       case 23 :
        offset += Compute_Within( &iA_data->input[jj]->inp_within, molecule, trj_crd, &printout[offset], intex);
        break;
       case 24:
        offset += Compute_Msdf ( &iA_data->input[jj]->inp_msdf, iA_data->opt, trj_crd, &printout[offset] );
        break;
       case 71:
         offset += Compute_Com ( &iA_data->input[jj]->inp_com, iA_data->opt, molecule, trj_crd, &printout[offset] );
         break;
       case 80:
         offset += Compute_Tilt ( &iA_data->input[jj]->inp_tilt, trj_crd, &printout[offset]);
         break;
       case 81:
         offset += Compute_Twist ( &iA_data->input[jj]->inp_twist, trj_crd, &printout[offset]);
         break;
       case 82:
         offset += Compute_Hole ( &iA_data->input[jj]->inp_hole, trj_crd, &printout[offset]);
         break;
       case 83:
         offset += Compute_Flux ( &iA_data->input[jj]->inp_flux, trj_crd, &printout[offset], ii);
         break;
       //case 99:
       // offset += Compute_Sample ( &iA_data->input[jj]->inp_sample, iA_data->opt, trj_crd, &printout[offset] );
       // break;
       default :
        break;
     } 
   }
   
   return;
}
// ---------------------------------------------------------------------
int CoreCalc( struct sopt *OPT, Molecule *molecule, CoorSetData * coorsetdata, struct inp_A  * iA_data, FILE * oA_f )
{
   int                  ii, jj;
   int                  intex;
   CoorSet              trj_crd;
   clock_t              start, end;
   double               elapsed;
   char                *tmp_print;
   int                  tmp_print_size;
   
   tmp_print_size = iA_data->outputstringsize;
   tmp_print = (char *)calloc( tmp_print_size+1, sizeof( char));
   CalloCoor( &trj_crd, coorsetdata->nato );
   
   if( OPT->PBC_FLAG )
   {
     trj_crd.pbc_flag = 1;
     trj_crd.pbc = calloc( 1, sizeof(Pbc));
     trj_crd.pbc->a_size = OPT->PBCBOX[0];
     trj_crd.pbc->b_size = OPT->PBCBOX[1];
     trj_crd.pbc->c_size = OPT->PBCBOX[2];
     trj_crd.pbc->angle1 = 0.0;
     trj_crd.pbc->angle2 = 0.0;
     trj_crd.pbc->angle3 = 0.0;
   }
   else if( coorsetdata->traj.hdr->varpbc )
   {
     trj_crd.pbc_flag = 1;
     trj_crd.pbc = calloc( 1, sizeof(Pbc));
   }
   else if( molecule->coor.pbc_flag == 1 )
   {
     trj_crd.pbc_flag = 1;
     trj_crd.pbc = calloc( 1, sizeof(Pbc));
     trj_crd.pbc->a_size = molecule->coor.pbc->a_size;
     trj_crd.pbc->b_size = molecule->coor.pbc->b_size;
     trj_crd.pbc->c_size = molecule->coor.pbc->c_size;
     trj_crd.pbc->angle1 = molecule->coor.pbc->angle1;
     trj_crd.pbc->angle2 = molecule->coor.pbc->angle2;
     trj_crd.pbc->angle3 = molecule->coor.pbc->angle3;
   }
   
   // unless user specifically tells us to neglect PBC ...
   if( trj_crd.pbc_flag == 1 && OPT->NOPBC_FLAG )
   {
     trj_crd.pbc_flag = -1;
     free( trj_crd.pbc );
     trj_crd.pbc = NULL;
   }
   

   /*  run over all coordinates set
    *  goes from 0 to the final number because begframe is defaulted to 1 */
   if (OPT->VERBOSE_FLAG) {
    start = clock();
   }
   for( ii=0; ii<=(coorsetdata->endframe-coorsetdata->begframe); ii+=coorsetdata->skip )
   {
	 intex = GetThisCoorSet( coorsetdata, ii, &trj_crd );
     fprintf(oA_f,"%7d ", intex); fflush(oA_f);
     memset( tmp_print, '\0', tmp_print_size);
     run_iA( ii, intex, &trj_crd, iA_data, tmp_print, molecule );
     fprintf( oA_f, "%s\n", tmp_print);
   }
   
   if(OPT->VERBOSE_FLAG)
   {
     end = clock();
     elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
     printf("Compute phase completed in \t%8.2f\"\n", elapsed);
     start = end;
   }
   fflush( oA_f );
   
   /* turn all modules OFF once computing is over */
   for( ii=0; ii<iA_data->nBEG; ii++ )
     iA_data->bonoffmodulescompute[ii] = 0;
   
   /* Now some post-processing for those who need it */
   for (jj=0; jj<iA_data->nBEG; jj++)
   {
     if( iA_data->bonoffmodulespost[jj] == 0 )
       continue;
     iA_data->bonoffmodulespost[jj] = 0;
     switch (iA_data->A_type[jj])
     {
       case  9 :
        if (iA_data->input[jj]->inp_rms.trjwrite)
          WriteTrjHeader ( iA_data->input[jj]->inp_rms.outtraj, "nframe");
        break;
       case  10 :
        iA_data->bonoffmodulespost[jj] = Post_Rmsf ( &iA_data->input[jj]->inp_rmsf, OPT, iA_data->nframe, molecule );
        break;
       case  11 :
        iA_data->bonoffmodulespost[jj] = Post_PCA ( &iA_data->input[jj]->inp_pca, OPT, iA_data->nframe, molecule );
        break;
       case  12 :
        iA_data->bonoffmodulespost[jj] = Post_pro ( &iA_data->input[jj]->inp_pro, OPT, iA_data->nframe, molecule );
        break;
       case  13 :
        iA_data->bonoffmodulespost[jj] = Post_Cluster ( &iA_data->input[jj]->inp_cluster );
        break;
        #ifdef CUDA
       case  84 :
        iA_data->bonoffmodulespost[jj] = Post_Gcluster ( &iA_data->input[jj]->inp_cluster, oA_f );
        break;
        #endif
       case  15 :
        iA_data->bonoffmodulespost[jj] = Post_Entropy ( &iA_data->input[jj]->inp_pca, iA_data->nframe, molecule );
        break;
       case  18 :
        iA_data->bonoffmodulespost[jj] = Post_PSG (&iA_data->input[jj]->inp_psn, iA_data->nframe, molecule, OPT );
       break;
       case  19 :
        iA_data->bonoffmodulespost[jj] = Post_CORR (&iA_data->input[jj]->inp_corr, molecule, OPT, &trj_crd);
        break;
       case  20 :
        iA_data->bonoffmodulespost[jj] = Post_SURFCLUST( &iA_data->input[jj]->inp_surf, OPT);
        break;
       case  21 :
        iA_data->bonoffmodulespost[jj] = Post_SURFCORR( &iA_data->input[jj]->inp_surf, OPT);
        break;
       case  24 :
        iA_data->bonoffmodulespost[jj] = Post_Msdf ( &iA_data->input[jj]->inp_msdf, OPT, iA_data->nframe, molecule );
        break;
       case  82 :
        iA_data->bonoffmodulespost[jj] = Post_Hole ( &iA_data->input[jj]->inp_hole, OPT, iA_data->nframe, &molecule );
        break;
       //case  99 :
       // iA_data->bonoffmodulespost[jj] = Post_Sample ( &iA_data->input[jj]->inp_sample, OPT, iA_data->nframe, molecule );
       // break;
       default :
        break;
     }
   }
   
   /* compute modules' flags are set according to what happened 
    * to the post modules' flags during "post" operations */
   for( ii=0; ii<iA_data->nBEG; ii++ )
     iA_data->bonoffmodulescompute[ii] = iA_data->bonoffmodulespost[ii];
   
   return 0;
}
// ---------------------------------------------------------------------
int CoreCalc_thread( struct sopt *OPT, Molecule *molecule, CoorSetData * coorsetdata, struct inp_A  * iA_data, struct threading *T_data, int begin, int end, int skipstep )
{
   int                  ii;
   int                  intex;
   CoorSet              trj_crd;
   
   //fprintf( stdout, "DEBUG: thread %d - c - inside corecalc\n", T_data->thread_rank);
   CalloCoor( &trj_crd, coorsetdata->nato );
   if( OPT->PBC_FLAG )
   {
     trj_crd.pbc_flag = 1;
     trj_crd.pbc = calloc( 1, sizeof(Pbc));
     trj_crd.pbc->a_size = OPT->PBCBOX[0];
     trj_crd.pbc->b_size = OPT->PBCBOX[1];
     trj_crd.pbc->c_size = OPT->PBCBOX[2];
     trj_crd.pbc->angle1 = 0.0;
     trj_crd.pbc->angle2 = 0.0;
     trj_crd.pbc->angle3 = 0.0;
   }
   else if( coorsetdata->traj.hdr->varpbc )
   {
     trj_crd.pbc_flag = 1;
     trj_crd.pbc = calloc( 1, sizeof(Pbc));
   }
   else if( molecule->coor.pbc_flag == 1 )
   {
     trj_crd.pbc_flag = 1;
     trj_crd.pbc = calloc( 1, sizeof(Pbc));
     trj_crd.pbc->a_size = molecule->coor.pbc->a_size;
     trj_crd.pbc->b_size = molecule->coor.pbc->b_size;
     trj_crd.pbc->c_size = molecule->coor.pbc->c_size;
     trj_crd.pbc->angle1 = molecule->coor.pbc->angle1;
     trj_crd.pbc->angle2 = molecule->coor.pbc->angle2;
     trj_crd.pbc->angle3 = molecule->coor.pbc->angle3;
   }
   
   // unless user specifically tells us to neglect PBC ...
   if( trj_crd.pbc_flag == 1 && OPT->NOPBC_FLAG )
   {
     trj_crd.pbc_flag = -1;
     free( trj_crd.pbc );
     trj_crd.pbc = NULL;
   }
   

   /* run over all coordinates set
    * goes from 0 to the final number because begframe is defaulted to 1 */
   for( ii = begin; ii < end; ii += skipstep )
   {
     intex = GetThisCoorSet( coorsetdata, ii, &trj_crd );
     run_iA( ii, intex, &trj_crd, iA_data, T_data->printout[ii], molecule );
     //fprintf( stdout, "DEBUG : %d %d (cct)\n", T_data->thread_rank, intex );
   }

   
   return 0;
}
// ---------------------------------------------------------------------
int isthereONmodule( short *modulelist, int nmodules )
{
  int   ii;
  
  for( ii=0; ii<nmodules; ii++ )
    if( modulelist[ii] > 0 )
      return 1;
  
  return 0;
}
// ---------------------------------------------------------------------
#ifdef THREADED
void * iA_Calc_thread ( void * pT_data ) 
{
   struct threading    *T_data;
   struct inp_A        *iA_data;       // input file structures
   Molecule            *molecule=NULL;
   CoorSetData         *coorsetdata;
   sopt                *opt;
   
   int                  ii;
   CoorSet              trj_crd;
   int                  firstframe, lastframe, skipstep, begin, end, papersize;
   char                *header;
   
   T_data       = (struct threading *) pT_data;
   
   iA_data      = &T_data->iA_data;
   opt          =  T_data->iA_data.opt;
   
   //fprintf( stdout, "DEBUG: entering thread %d\n", T_data->thread_rank);
   
   if( T_data->thread_rank == 0 )
   {
     molecule    =  T_data->iA_data.molecule;
     coorsetdata =  T_data->coorsetdata;
   }
   else
   {
     coorsetdata = calloc( 1, sizeof( CoorSetData));

     GetMolecule( opt->IMOL_FILE, molecule, opt->IMOL_TYPE );
     iA_data->molecule = molecule;
     
     coorsetdata->traj.hdetail = FULL;
     CoorSetInit( opt, coorsetdata );
     CheckRanges( opt, molecule, iA_data, coorsetdata);
   
     iA_data->nframe = coorsetdata->i_nUsedFrames;

     header = Read_iA ( opt, iA_data, &papersize, molecule, coorsetdata);
     free(header);
   }
   
   CalloCoor( &trj_crd, coorsetdata->nato );

   //fprintf( stdout, "DEBUG: thread %d - a\n", T_data->thread_rank);
   if( opt->PBC_FLAG )
   {
     trj_crd.pbc_flag = 1;
     trj_crd.pbc = calloc( 1, sizeof(Pbc));
     trj_crd.pbc->a_size = opt->PBCBOX[0];
     trj_crd.pbc->b_size = opt->PBCBOX[1];
     trj_crd.pbc->c_size = opt->PBCBOX[2];
     trj_crd.pbc->angle1 = 0.0;
     trj_crd.pbc->angle2 = 0.0;
     trj_crd.pbc->angle3 = 0.0;
   }
   else if( coorsetdata->traj.hdr->varpbc )
   {
     trj_crd.pbc_flag = 1;
     trj_crd.pbc = calloc( 1, sizeof(Pbc));
   }
   else if( molecule->coor.pbc_flag == 1 )
   {
     trj_crd.pbc_flag = 1;
     trj_crd.pbc = calloc( 1, sizeof(Pbc));
     trj_crd.pbc->a_size = molecule->coor.pbc->a_size;
     trj_crd.pbc->b_size = molecule->coor.pbc->b_size;
     trj_crd.pbc->c_size = molecule->coor.pbc->c_size;
     trj_crd.pbc->angle1 = molecule->coor.pbc->angle1;
     trj_crd.pbc->angle2 = molecule->coor.pbc->angle2;
     trj_crd.pbc->angle3 = molecule->coor.pbc->angle3;
   }
   
   // unless user specifically tells us to neglect PBC ...
   if( trj_crd.pbc_flag == 1 && opt->NOPBC_FLAG )
   {
     trj_crd.pbc_flag = -1;
     free( trj_crd.pbc );
     trj_crd.pbc = NULL;
   }
   
   begin       = T_data->thread_rank;
   firstframe  = coorsetdata->begframe; // + T_data->thread_rank;
   lastframe   = coorsetdata->endframe;
   end         = lastframe - firstframe + 1;
   skipstep    = coorsetdata->skip * T_data->thread_size;
   
   iA_data->threaded = 1;
   
   //fprintf( stdout, "DEBUG: thread %d - b\n", T_data->thread_rank);
   /* turn all found computing modules ON */
   iA_data->bonoffmodulescompute = (short *)calloc( iA_data->nBEG, sizeof( short) );
   for( ii=0; ii<iA_data->nBEG; ii++ )
     iA_data->bonoffmodulescompute[ii] = 1;
   iA_data->bonoffmodulespost = (short *)calloc( iA_data->nBEG, sizeof( short) );
   for( ii=0; ii<iA_data->nBEG; ii++ )
     iA_data->bonoffmodulespost[ii] = 1;
   
   if( isthereONmodule( iA_data->bonoffmodulescompute, iA_data->nBEG ) )
   {
     CoreCalc_thread( opt, molecule, coorsetdata, iA_data, T_data, begin, end, skipstep );
     CoorSetRewind( coorsetdata );
   }
   
   //fprintf( stdout, "DEBUG: thread %d - c\n", T_data->thread_rank);
   
   /* run over all coordinates set
    * goes from 0 to the final number because begframe is defaulted to 1 */
   //for( ii = begin; ii < end; ii += skipstep )
   //{
     //intex = GetThisCoorSet( coorsetdata, ii, &trj_crd );
     //run_iA( ii, intex, &trj_crd, iA_data, T_data->printout[ii], molecule );
     //fprintf( stdout, "DEBUG : %d %d\n", T_data->thread_rank, intex );
   //}

   return pT_data;
}
#endif
// ---------------------------------------------------------------------
void iA_Calc ( struct sopt   *OPT  ) 
{
   FILE                *oA_f;
   
   Molecule             molecule;
   CoorSet              trj_crd;

   struct inp_A         iA_data;        // input file structures
   //struct inp_A        *iA_datas;       // input file structures
   int                  ii;
   int                  papersize=0;
   CoorSetData          coorsetdata;
   //char                *tmp_print;
   char                *header;
   
   extern short int     frame_par;
   extern short int     module_par;
   
   molecule.ext_flag = 0;
   molecule.coor.pbc = NULL;
   molecule.coor.pbc_flag = 0;
   molecule.filled = calloc( 4, sizeof(char));
   memset( molecule.filled, '\0', 4);
   GetMolecule ( OPT->IMOL_FILE, &molecule, 0 );
   iA_data.molecule = &molecule;
   
   coorsetdata.traj.hdetail = FULL;

   CoorSetInit( OPT, &coorsetdata );
   
   CheckRanges( OPT, &molecule, &iA_data, &coorsetdata);
   CalloCoor( &trj_crd, coorsetdata.nato );
   
   if( OPT->PBC_FLAG )
   {
     trj_crd.pbc_flag = 1;
     trj_crd.pbc = calloc( 1, sizeof(Pbc));
     trj_crd.pbc->a_size = OPT->PBCBOX[0];
     trj_crd.pbc->b_size = OPT->PBCBOX[1];
     trj_crd.pbc->c_size = OPT->PBCBOX[2];
     trj_crd.pbc->angle1 = 0.0;
     trj_crd.pbc->angle2 = 0.0;
     trj_crd.pbc->angle3 = 0.0;
   }
   else if( coorsetdata.traj.hdr->varpbc )
   {
     trj_crd.pbc_flag = 1;
     trj_crd.pbc = calloc( 1, sizeof(Pbc));
   }
   else if( molecule.coor.pbc_flag == 1 )
   {
     trj_crd.pbc_flag = 1;
     trj_crd.pbc = calloc( 1, sizeof(Pbc));
     trj_crd.pbc->a_size = molecule.coor.pbc->a_size;
     trj_crd.pbc->b_size = molecule.coor.pbc->b_size;
     trj_crd.pbc->c_size = molecule.coor.pbc->c_size;
     trj_crd.pbc->angle1 = molecule.coor.pbc->angle1;
     trj_crd.pbc->angle2 = molecule.coor.pbc->angle2;
     trj_crd.pbc->angle3 = molecule.coor.pbc->angle3;
   }
   
   // unless user specifically tells us to neglect PBC ...
   if( trj_crd.pbc_flag == 1 && OPT->NOPBC_FLAG )
   {
     trj_crd.pbc_flag = -1;
     free( trj_crd.pbc );
     trj_crd.pbc = NULL;
   }
   coorsetdata.pbc_flag = trj_crd.pbc_flag;
   
   // Open Output File, Read Input File
   if ( OPT->OA_FLAG )
     oA_f = O_File ( OPT->OA_FILE, "w" );
   else
     oA_f = stdout;
    
   
   fprintf(oA_f, "#Fr     ");

  // reading input file
   iA_data.nframe = coorsetdata.i_nUsedFrames;
   header = Read_iA ( OPT, &iA_data, &papersize, &molecule, &coorsetdata);
   fprintf(oA_f,"%s\n", header); fflush( oA_f);
   free(header);
   //tmp_print = calloc( papersize, sizeof(char));
   iA_data.outputstringsize = papersize;
   
   iA_data.oA_f = oA_f;
   iA_data.opt  = OPT;
   
   /* turn all found computing modules ON */
   iA_data.bonoffmodulescompute = (short *)calloc( iA_data.nBEG, sizeof( short) );
   for( ii=0; ii<iA_data.nBEG; ii++ )
     iA_data.bonoffmodulescompute[ii] = 1;
   iA_data.bonoffmodulespost = (short *)calloc( iA_data.nBEG, sizeof( short) );
   for( ii=0; ii<iA_data.nBEG; ii++ )
     iA_data.bonoffmodulespost[ii] = 1;
   
   // Write Output File Headers
   //fprintf(oA_f, "#   nFr ");

   while( isthereONmodule( iA_data.bonoffmodulescompute, iA_data.nBEG ) )
   {
     //fprintf( stdout, "DEBUG: *** running CoreCalc !\n"); fflush(stdout);
     CoreCalc( OPT, &molecule, &coorsetdata, &iA_data, oA_f );
     CoorSetRewind( &coorsetdata );
   }
   
   ///*  run over all coordinates set
    //*  goes from 0 to the final number because begframe is defaulted to 1 */
   //for( ii=0; ii<=(coorsetdata.endframe-coorsetdata.begframe); ii+=coorsetdata.skip )
   //{
     //intex = GetThisCoorSet( &coorsetdata, ii, &trj_crd );
     //fprintf(oA_f,"%7d ", intex); fflush(oA_f);
     //memset( tmp_print, '\0', papersize);
     //run_iA( ii, intex, &trj_crd, &iA_data, tmp_print, &molecule );
     //fprintf( oA_f, "%s\n", tmp_print);
   //}
   
   //if(OPT->VERBOSE_FLAG)
   //{
     //end = clock();
     //elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
     //printf("Compute phase completed in \t%8.2f\"\n", elapsed);
     //start = end;
   //}
   //fflush( oA_f );
   
   ///* turn all modules OFF once computing is over */
   //for( ii=0; ii<iA_data.nBEG; ii++ )
     //iA_data.bonoffmodulescompute[ii] = 0;
   
   ///* Now some post-processing for those who need it */
   //for (jj=0; jj<iA_data.nBEG; jj++)
   //{
     //if( iA_data.bonoffmodulespost[jj] == 0 )
       //continue;
     //switch (iA_data.A_type[jj])
     //{
       //case  9 :
        //if (iA_data.input[jj]->inp_rms.trjwrite)
          //WriteTrjHeader ( iA_data.input[jj]->inp_rms.outtraj, "nframe");
        //break;
       //case  10 :
        //Post_Rmsf ( &iA_data.input[jj]->inp_rmsf, OPT, iA_data.nframe, &molecule );
        //break;
       //case  11 :
        //Post_PCA ( &iA_data.input[jj]->inp_pca, OPT, iA_data.nframe, &molecule );
        //break;
       //case  12 :
        //Post_pro ( &iA_data.input[jj]->inp_pro, OPT, iA_data.nframe, &molecule );
        //break;
       //case  13 :
        //Post_Cluster ( &iA_data.input[jj]->inp_cluster );
        //break;
       //case  15 :
        //Post_Entropy ( &iA_data.input[jj]->inp_pca, iA_data.nframe, &molecule );
        //break;
       //case  18 :
        //Post_PSG (&iA_data.input[jj]->inp_psn, iA_data.nframe, &molecule, OPT );
       //break;
       //case  19 :
        //Post_CORR ( &iA_data.input[jj]->inp_corr, &molecule, OPT, &trj_crd);
        //break;
     //case  24 :
      //Post_Rmsfd ( &iA_data.input[jj]->inp_rmsfd, OPT, iA_data.nframe, &molecule );
      //break;
       //case  99 :
        //Post_Sample ( &iA_data.input[jj]->inp_sample, OPT, iA_data.nframe, &molecule );
        //break;
       //default :
        //break;
     //}
   //}
   
   //if(OPT->VERBOSE_FLAG)
   //{
     //end = clock();
     //elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
     //printf("Post phase completed in \t%8.2f\"\n", elapsed);
     //start = end;
   //}
   
   fclose (oA_f);
   
   return;
}                
// ------------------------------------------------------------------
char * Read_iA( struct sopt  *OPT, struct inp_A *iA_data, int *papersize, Molecule *molecule, CoorSetData *coorsetdata )
{
    FILE       *iA_f;
    char      **input_text;
    int         n_inputlines, offset;
    int         ii, jj;
    char        type[16], textline[20480];
    _iaOpt     *current;
    char       *printout;
    
    struct temp_ll *temp_ll, *begin, *prev;
    
    if( OPT->IA_FLAG )
    { 
      iA_f = O_File ( OPT->IA_FILE, "r");
      n_inputlines = 0;
      while ( !feof(iA_f))
      {
        if ( textline[0] != '#' && textline[0] != '\n' )
          ++n_inputlines;
        if(fgets(textline, 1024, iA_f)==NULL && ( !feof(iA_f) || ferror(iA_f) ))
        {
          fprintf(stderr, "Warning! Premature end of file reached!\n");
        }
      }
      n_inputlines--;
      rewind(iA_f);
      input_text = calloc( n_inputlines, sizeof(char * ));
      for( ii=0; ii< n_inputlines;  )
      {
        memset( textline, '\0', 20480);
        if(fgets(textline, 20480, iA_f)==NULL && ( !feof(iA_f) || ferror(iA_f) ))
        {
          fprintf(stderr, "Warning! Premature end of file reached!\n");
        }
        if ( textline[0] != '#' && textline[0] != '\n' )
        {
          input_text[ii] = calloc( 20480, sizeof(char));
          sprintf( input_text[ii], "%s", textline);
          ii++;
        }
      }
      fclose( iA_f );
      iA_data->nBEG = Get_nBEG ( input_text, n_inputlines );
      iA_data->begIndex = calloc( iA_data->nBEG, sizeof( int));
      spotBeg( input_text, n_inputlines, iA_data->begIndex);
      iA_data->A_type = calloc ( iA_data->nBEG, sizeof ( int ));
      Get_Atype ( input_text, n_inputlines, iA_data->A_type, iA_data->nBEG );
    }
    else if( OPT->Ia_FLAG )
    {	
      n_inputlines = OPT->iaopt->noptions+2;
      input_text = calloc( n_inputlines, sizeof(char * ));
      for( ii=0; ii< n_inputlines; ii++ )
        input_text[ii] = calloc( 20480, sizeof(char));
      sprintf( input_text[0], "BEGIN %s\n", OPT->Ia_STRING);
      current = OPT->iaopt;
      for( ii=0; ii< OPT->iaopt->noptions; ii++ )
      {
        sprintf( input_text[ii+1], "%s %s\n", current->option, current->value );
        current = current->next;
      }
      sprintf( input_text[ii+1], "END\n");
      
      iA_data->nBEG = 1;
      iA_data->begIndex = calloc( 1, sizeof(int));
      iA_data->begIndex[0] = 0;
      iA_data->A_type = calloc ( 1, sizeof ( int ));
      sscanf(OPT->Ia_STRING, "%15s", type);
      for( jj=0; jj<strlen(type); jj++ )
        type[jj] = tolower(type[jj]);
      iA_data->A_type[0] = whichModule( type );
    }
    else
    {
       fprintf( stderr, "Hey, no Ix_FLAG is set!\n");
    }
    
    iA_data->input  = calloc ( iA_data->nBEG, sizeof(union input * ));
    iA_data->printout = calloc ( iA_data->nBEG, sizeof(char ** ));
    for( ii=0; ii<iA_data->nBEG; ii++)
      iA_data->printout[ii] = calloc ( iA_data->nframe, sizeof(char * ));
    
    offset = 0;
    papersize[0] = 0;

    temp_ll = calloc( 1, sizeof(struct temp_ll));
    temp_ll->next = NULL;
    temp_ll->title = calloc( 2048, sizeof(char));
    begin = temp_ll;
    
    for ( ii=0; ii<iA_data->nBEG; ii++)
    {
      iA_data->input[ii] = calloc ( 1, sizeof(union input));
      switch (iA_data->A_type[ii])
      {
        case '0':
         printf("Error in input file: could not recognize option #%d\n", ii+1);
         exit(0);
         break;
        case  1:
         free(temp_ll->title);
         temp_ll->size = Read_iDist( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_dist , temp_ll , molecule );
         break;
        case  2:
         free(temp_ll->title);
         temp_ll->size = Read_iContacts( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_contacts , temp_ll , molecule );
         break;
        case  3:
         temp_ll->size = Read_iAngle( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_angle , temp_ll->title , molecule );
         break;
        case  4:
         temp_ll->size = Read_iDihe( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_dihe , temp_ll->title , molecule );
         break;
        case  5:
         temp_ll->size = Read_iRgyr( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_rgyr , temp_ll->title , molecule, coorsetdata );
         break;
        case  6:
         temp_ll->size = Read_iHB( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_hb , temp_ll->title , molecule );
         break;
        case  7:
         temp_ll->size = Read_iP( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_P , temp_ll->title , molecule );
         break;
        case  8:
         temp_ll->size = Read_idrms( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_drms , temp_ll->title , molecule , OPT);
         break;
        case  9:
         temp_ll->size = Read_iRmsd( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_rms , temp_ll->title , molecule, coorsetdata, OPT );
         break;
        case 10:
         temp_ll->size = Read_iRmsf( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_rmsf , temp_ll->title , molecule );
         break;
        case 11:
         temp_ll->size = Read_ipca( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_pca , temp_ll->title , molecule , iA_data->nframe);
         break;
        case 12:
         temp_ll->size = Read_ipro( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_pro , temp_ll->title , molecule );
         break;
        case 13:
         temp_ll->size = Read_iCluster( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_cluster , temp_ll->title , molecule, coorsetdata, iA_data->nframe );
         iA_data->input[ii]->inp_cluster.pbcflag = OPT->PBC_FLAG;
         iA_data->input[ii]->inp_cluster.pbcbox  = OPT->PBCBOX;
         if( OPT->SKIP_STEP != 1 )
           iA_data->input[ii]->inp_cluster.step  = OPT->SKIP_STEP;
         break;
         #ifdef CUDA
        case 84:
         temp_ll->size = Read_iGcluster( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_cluster , temp_ll->title , molecule, coorsetdata, iA_data->nframe );
         iA_data->input[ii]->inp_cluster.pbcflag = OPT->PBC_FLAG;
         iA_data->input[ii]->inp_cluster.pbcbox  = OPT->PBCBOX;
         if( OPT->SKIP_STEP != 1 )
           iA_data->input[ii]->inp_cluster.step  = OPT->SKIP_STEP;
         break;
         #endif
        case 14:
         temp_ll->size = Read_iCAssign( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_cassign , temp_ll->title , molecule, coorsetdata, iA_data->nframe );
         iA_data->input[ii]->inp_cassign.pbcflag = OPT->PBC_FLAG;
         iA_data->input[ii]->inp_cassign.pbcbox  = OPT->PBCBOX;
         iA_data->input[ii]->inp_cassign.skip  = OPT->SKIP_STEP;
         break;
        case 15:
         temp_ll->size = Read_iEntropy( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_pca , temp_ll->title , molecule );
         break;
        case 16:
         temp_ll->size = Read_iSSA( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_ss , temp_ll->title , molecule , iA_data->nframe);
         break;
        case 17:
         temp_ll->size = Read_iEnm( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_enm, temp_ll->title, molecule  );
         break;
        case 18:
         temp_ll->size = Read_iPSG( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_psn , temp_ll->title , molecule, iA_data->nframe, OPT );
         break;
        case 19:
         temp_ll->size = Read_CORR( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_corr , molecule , temp_ll->title, iA_data->nframe );
         break;
        case 20 :
         temp_ll->size = Read_SURFCLUST( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_surf , molecule , temp_ll->title, iA_data->nframe );
         break;
        case 21 :
         temp_ll->size = Read_SURFCORR( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_surf , molecule , temp_ll->title, iA_data->nframe );
         break;
        case 22 :
         temp_ll->size = Read_SURF( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_surf , molecule , temp_ll->title );
         break;
        case 23 :
         temp_ll->size = Init_Within(input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_within , molecule , temp_ll->title);
         break;
        case 24:
         temp_ll->size = Read_iMsdf( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_msdf , temp_ll->title , molecule, iA_data->nframe );
         break;
        case  71:
         temp_ll->size = Read_iCom( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_com , temp_ll->title , molecule, coorsetdata->pbc_flag );
          break;
        case 80: 
          temp_ll->size = Read_iTilt ( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_tilt, temp_ll->title, molecule );
          break;
        case 81: 
          temp_ll->size = Read_iTwist ( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_twist, temp_ll->title, molecule );
          break;
        case 82: 
          temp_ll->size = Read_iHole ( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_hole, temp_ll->title, molecule );
          break;
        case 83: 
          temp_ll->size = Read_iFlux ( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_flux, temp_ll->title, molecule );
          break;
        //case 99 :
        // temp_ll->size = Read_iSample( input_text, iA_data->begIndex[ii], &iA_data->input[ii]->inp_sample , temp_ll->title , molecule );
        // break;
        default :
         printf ("Didn't find use for option: %c\n", iA_data->A_type[ii]);
         exit(0);
         break;
      }
      papersize[0] += temp_ll->size;
      temp_ll->next = calloc( 1, sizeof(struct temp_ll));
      temp_ll = temp_ll->next;
      temp_ll->next = NULL;
      temp_ll->title = calloc( 2048, sizeof(char));
    }
    /* copy titles in string to return */
    printout = calloc( papersize[0]+2, sizeof(char));
    temp_ll = begin;
    offset = 0;
    for( ii=0; ii<iA_data->nBEG; ii++)
    {
      sprintf( &printout[offset], "%s", temp_ll->title);
      offset += temp_ll->size;
      /* skip to next title and free previous area */
      prev = temp_ll;
      temp_ll = temp_ll->next;
      free(prev);
    }
    
    if( OPT->PBC_FLAG )
    {
     for ( ii=0; ii<iA_data->nBEG; ii++)
     {
      if( iA_data->A_type[ii] == 10 )
      {
       printf("Error: pbc not implemented yet with rmsd\n");
       exit(0);
      }
      if( iA_data->A_type[ii] == 12 &&  iA_data->input[ii]->inp_cluster.distance == 1)
      {
       printf("Error: pbc not implemented yet with rmsd-based clustering\n");
       exit(0);
      }
     }
    }

    return printout;
}
// ------------------------------------------------------------------
// ------------------------------------------------------------------
int Get_nBEG ( char **input, int nlines ) 
{
    int nBEG = 0;
    int nEND = 0;
    int ii;
    

    for( ii=0; ii<nlines; ii++ )
    {
      if( !strncmp( input[ii], "BEGIN", 5 ))
      {
        nBEG++;
      }
      if( !strncmp( input[ii], "END", 3 ))
      {
        nEND++;
      }
    }
    
    if( nBEG != nEND )
    {
     fprintf( stderr, "Syntax Error! number of BEGIN (%d) differs from number of END (%d)\n", nBEG, nEND);
     exit(0);
    }
    
    return nBEG;
}	
// ------------------------------------------------------------------
void spotBeg ( char **input, int nlines, int *begIndex ) 
{
    int ii, jj;

    jj = 0;
    for( ii=0; ii<nlines-1; ii++ )
    {
      if( !strncmp( input[ii], "BEGIN", 5 ))
      {
        begIndex[jj] = ii;
        jj++;
      }
    }
    
    return;
}	
// ------------------------------------------------------------------
int whichModule( char *typestring )
{
   if(!strncmp( typestring, "distance", 8))
    return  1;
   if(!strncmp( typestring, "contacts", 8))
    return  2;
   if(!strncmp( typestring, "angle", 5))
    return  3;
   if(!strncmp( typestring, "dihedral", 8))
    return  4;
   if(!strncmp( typestring, "rgyr", 4))
    return  5;
   if(!strncmp( typestring, "hbond", 5))
    return  6;
   if(!strncmp( typestring, "orienta", 7))
    return  7;
   if(!strncmp( typestring, "drms", 4))
    return  8;
   if(!strncmp( typestring, "rmsd", 4))
    return  9;
   if(!strncmp( typestring, "msdf", 4))
    return 24;
   if(!strncmp( typestring, "rmsf", 4))
    return 10;
   if(!strncmp( typestring, "pca", 3))
    return 11;
   if(!strncmp( typestring, "project", 7))
    return 12;
   if(!strncmp( typestring, "cluster", 7))
    return 13;
   if(!strncmp( typestring, "cassign", 7))
    return 14;
   if(!strncmp( typestring, "entropy", 7))
    return 15;
   if(!strncmp( typestring, "ssa", 7))
    return 16;
   if(!strncmp( typestring, "enm", 3))
    return 17;
   if(!strncmp( typestring, "psn", 3))
    return 18;
   if(!strncmp( typestring, "corr", 4))
    return 19 ;
   if(!strncmp( typestring, "surfcluster", 11))
    return 20 ;
   if(!strncmp( typestring, "surfcorr", 8))
    return 21 ;
   if(!strncmp( typestring, "surf", 4))
    return 22 ;
   if(!strncmp( typestring, "within", 6))
    return 23 ;
   if(!strncmp( typestring, "com", 3))
    return 71 ;
   if(!strncmp( typestring, "tilt", 4))
    return 80 ;
   if(!strncmp( typestring, "twist", 3))
    return 81 ;
   if(!strncmp( typestring, "hole", 3))
    return 82 ;
   if(!strncmp( typestring, "flux", 3))
    return 83 ;
   if(!strncmp( typestring, "gcluster", 3)) {
	   #ifdef CUDA
		return 84 ;
		#endif
		fprintf(stderr, "gcluster requires compilation on a system with CUDA!\n");
	}
   if(!strncmp( typestring, "sample", 6))
    return 99 ;
   
   fprintf( stderr, "Unknown analysis flag %s\n", typestring);
   exit(0);
}
// ------------------------------------------------------------------
int whichEModule( char *typestring )
{
   if(!strncmp( typestring, "logbin", 6))
    return  1;
   if(!strncmp( typestring, "ka", 2))
    return  2;
   if(!strncmp( typestring, "basin", 5))
    return  3;
   if(!strncmp( typestring, "mfptnet", 7))
    return  4;
   if(!strncmp( typestring, "mfpt", 4))
    return  5;
   if(!strncmp( typestring, "pfoldfnet", 9))
    return  6;
   if(!strncmp( typestring, "pfoldf", 6))
    return  7;
   if(!strncmp( typestring, "equil", 5))
    return  8;
   if(!strncmp( typestring, "psnpath", 10))
    return  9;
   if(!strncmp( typestring, "elec", 4))
    return 10;
   if(!strncmp( typestring, "cluster", 7))
    return 13;
   if(!strncmp( typestring, "volume", 6))
    return 14;
   if(!strncmp( typestring, "psnparam", 8))
    return  15;
   if(!strncmp( typestring, "ring", 4))
    return 70;
   if(!strncmp( typestring, "sample", 6))
    return 99 ;
   
   fprintf( stderr, "Unknown analysis flag %s\n", typestring);
   exit(0);
}
// ------------------------------------------------------------------
void Get_Atype ( char **input_text, int nlines, int *A_type, int nBEG ) 
{
    int   ii = 0, jj, kk;
    char  type[16];
    
    memset( type, '\0', 16);
    jj = 0;
    for( ii=0; ii<nlines; ii++ )
    {
      if( !strncmp(input_text[ii],"BEGIN", 5) )
      {
        A_type[jj] = 0;
        sscanf(input_text[ii], "BEGIN %s", type);
        for( kk=0; kk<16; kk++ )
          type[kk] = tolower(type[kk]);
        A_type[jj] = whichModule( type );
        jj++;
      }
    }

    if ( jj != nBEG )
    {
      fprintf( stderr, "Error in Input files (%d-%d): number of BEG unclear?!?\n", jj, nBEG);
      exit(0);
    }

    return;
}
// ------------------------------------------------------------------
void Get_Etype ( char **input_text, int nlines, int *E_type, int nBEG ) 
{
    int   ii = 0, jj, kk;
    char  type[16];
    
    memset( type, '\0', 16);
    jj = 0;
    for( ii=0; ii<nlines; ii++ )
    {
      if( !strncmp(input_text[ii],"BEGIN", 5) )
      {
        E_type[jj] = 0;
        sscanf(input_text[ii], "BEGIN %s", type);
        for( kk=0; kk<16; kk++ )
          type[kk] = tolower(type[kk]);
        E_type[jj] = whichEModule( type );
        jj++;
      }
    }

    if ( jj != nBEG )
    {
      fprintf( stderr, "Error in Input files (%d-%d): number of BEG unclear?!?\n", jj, nBEG);
      exit(0);
    }

    return;
}
// ------------------------------------------------------------------
int Get_nSeles ( char **input, int intex ) 
{
    int   nSeles = 0;
    
    while( strncmp ( input[intex], "END", 3) )
    {
      if( !strncmp ( input[intex], "--SELE", 6) )
        nSeles ++;
      intex++;
    }
    
    return nSeles;
}
// ------------------------------------------------------------------
/* what if module does not use molecule or trajectory data: */
void Calc2 ( struct sopt *OPT )
{
   int          ii, jj;
   FILE        *iA_f, *outfile;
   char       **input_text;
   int          n_inputlines;
   char         type[16], textline[1024];
   _iaOpt      *current;
   struct inp_A         iA_data;        // input file structures
   
   outfile = NULL;
   
   if( OPT->IE_FLAG )
   {
     iA_f = O_File ( OPT->IA_FILE, "r");
     n_inputlines = 0;
     memset( textline, '\0', 1024);
     if(fgets(textline, 1024, iA_f)==NULL && ( !feof(iA_f) || ferror(iA_f) ))
     {
       fprintf(stderr, "Warning! Premature end of file reached!\n");
     }
     while ( !feof(iA_f))
     {
       if ( textline[0] != '#' && textline[0] != '\n' )
         ++n_inputlines;
       if(fgets(textline, 1024, iA_f)==NULL && ( !feof(iA_f) || ferror(iA_f) ))
       {
         fprintf(stderr, "Warning! Premature end of file reached!\n");
       }
     }
     
     rewind(iA_f);
     input_text = calloc( n_inputlines, sizeof(char * ));
     for( ii=0; ii< n_inputlines;  )
     {
       memset( textline, '\0', 1024);
       if(fgets(textline, 1024, iA_f)==NULL && ( !feof(iA_f) || ferror(iA_f) ))
       {
         fprintf(stderr, "Warning! Premature end of file reached!\n");
       }
       if ( textline[0] != '#' && textline[0] != '\n' )
       {
         input_text[ii] = calloc( 1024, sizeof(char));
         sprintf( input_text[ii], "%s", textline);
         ii++;
       }
     }
     fclose( iA_f );
     iA_data.nBEG = Get_nBEG ( input_text, n_inputlines );
     iA_data.begIndex = calloc( iA_data.nBEG, sizeof( int));
     spotBeg( input_text, n_inputlines, iA_data.begIndex);
     iA_data.A_type = calloc ( iA_data.nBEG, sizeof ( int ));
     for( ii=0; ii<iA_data.nBEG; ii++ )
       iA_data.A_type[ii] = 0;
     Get_Etype ( input_text, n_inputlines, iA_data.A_type, iA_data.nBEG );
   }
   else if( OPT->Ie_FLAG )
   {
     n_inputlines = OPT->iaopt->noptions+3;
     input_text = calloc( n_inputlines, sizeof(char * ));
     for( ii=0; ii< n_inputlines; ii++ )
       input_text[ii] = calloc( 1024, sizeof(char));
     sprintf( input_text[0], "BEGIN %s\n", OPT->Ia_STRING);
     current = OPT->iaopt;
     for( ii=0; ii< OPT->iaopt->noptions; ii++ )
     {
       sprintf( input_text[ii+1], "%s %s\n", current->option, current->value );
       current = current->next;
     }
     sprintf( input_text[ii+1], "END\n");
     
     iA_data.nBEG = 1;
     iA_data.begIndex = calloc( 1, sizeof(int));
     iA_data.begIndex[0] = 0;
     iA_data.A_type = calloc ( 1, sizeof ( int ));
     iA_data.A_type[0] = 0;
     sscanf(OPT->Ia_STRING, "%s", type);
     for( jj=0; jj<sizeof(type); jj++ )
       type[jj] = tolower(type[jj]);
     iA_data.A_type[0] = whichEModule( type );
   }
   else
   {
     fprintf( stderr, "Calc2 error: not iE nor ie\n");
     exit(0);
   }
   
   for (ii=0; ii<iA_data.nBEG; ii++)
   {
     switch(iA_data.A_type[ii])
     {
       case  1:
         LogBin( input_text, iA_data.begIndex[ii], outfile );
         break;
       case  2:
         Ka( input_text, iA_data.begIndex[ii], outfile );
         break;
       case  3:
         Basin( input_text, iA_data.begIndex[ii], outfile );
         break;
       case  4:
         MFPTnet( input_text, iA_data.begIndex[ii], outfile );
         break;
       case  5:
         MFPT( input_text, iA_data.begIndex[ii], outfile );
         break;
       case  6:
         PFoldfnet( input_text, iA_data.begIndex[ii], outfile );
         break;
       case  7:
         PFoldf( input_text, iA_data.begIndex[ii], outfile );
         break;
       case  8:
         Equil( input_text, iA_data.begIndex[ii], outfile );
         break;
       case  9 :
         PSNPATH( input_text, iA_data.begIndex[ii] );
         break;
       case 10 :
         ELEC( input_text, iA_data.begIndex[ii] );
         break;
       case 11 :
         break;
       case 12 :
         break;
       case 13 :
         FileClustering( input_text, iA_data.begIndex[ii] );
         break;
       case 14 :
         CalcVolumes( input_text, iA_data.begIndex[ii], outfile );
         break;
       case 15 :
         InitPSNParam( input_text, iA_data.begIndex[ii] );
         break;
       case 70 :
         Compute_Ring( input_text, iA_data.begIndex[ii] );
         break;
     }
   }
   return;
}
// ------------------------------------------------------------------
