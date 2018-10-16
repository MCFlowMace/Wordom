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
/*! \file main.c
 \brief Main function of wordom
 
 The main function in main.c collects command-line options, fills the
 relevant structure and runs the appropriate procedure. Hearth and
 interface of the program
*/
#include "wordom.h"
#include "messages.h"
#include "fileio.h"
#include "tools.h"
#include "datahandler.h"
#include "analysis.h"

short int is_par       = 0;
short int frame_par    = 0;
short int module_par   = 0;
short int no_frame_par = 0;
extern short int verbose;
// ------------------------------------------------------------------

int main (int argc, char **argv)
{
	
  int           c;
  int           ii;
  int           iaoptions;
  
  sopt OPT;
  
  setOPT (&OPT);
//  int opterr;
  // =========================================
  // === first, let's sort the options out ===
  // =========================================
  
  // if analysis is called from the command line
  // extraneous options have to be collected
  // to be passed to the read_ia functions
  OPT.Ia_FLAG = 0;
  iaoptions = 0;
	
  for( ii=0; ii<argc; ii++ )
  {
	
    if( !strncmp( argv[ii], "-ia", 3) )
    {
      OPT.Ia_FLAG = 1; // found! take notice!
      OPT.Ia_STRING = argv[ii+1];
      OPT.iaopt = malloc( sizeof( _iaOpt));
      OPT.iaopt->noptions = 0;
      opterr = 0;
      break;
    }
    else if( !strncmp( argv[ii], "-ie", 3) )
    {
      OPT.Ie_FLAG = 1; // found! take notice!
      OPT.Ia_STRING = argv[ii+1];
      OPT.iaopt = malloc( sizeof( _iaOpt));
      OPT.iaopt->noptions = 0;
      opterr = 0;
      break;
    }
  }

  while (1)
  {
    static struct option long_options[] =
    {
     {"itrj", required_argument, 0, 'T'},
     {"imol", required_argument, 0, 'M'},
     {"otrj", required_argument, 0, 't'},
     {"omol", required_argument, 0, 'm'},
     {"atrj", required_argument, 0, 'r'},
     {"amol", required_argument, 0, 'l'},
     {"skip", required_argument, 0, 'k'},
     {"beg" , required_argument, 0, 'b'},
     {"end" , required_argument, 0, 'e'},
     {"head", required_argument, 0, 'c'},
     {"help", optional_argument, 0, 'h'},
     {"mod" , required_argument, 0, 'd'},
     {"nc"  , required_argument, 0, 'p'},
     {"mono",       no_argument, 0, 'g'},
     {"conv",       no_argument, 0, 'n'},
     {"combine", required_argument, 0, '4'},
     {"avg" ,       no_argument, 0, 'v'},
     {"verbose",    no_argument, 0, 'V'},
     {"info",       no_argument, 0, 'i'},
//   analysis options
     { "iA" , required_argument, 0, 'A'},
     { "iE" , required_argument, 0, 'E'},
     {"otxt", required_argument, 0, 'o'},
     {"xpbc", required_argument, 0, '0'},
     {"ypbc", required_argument, 0, '1'},
     {"zpbc", required_argument, 0, '2'},
     {"nopbc",      no_argument, 0, '3'},
     {"sele", required_argument, 0, 's'},
     {"checksele" , required_argument, 0, 'j'},
     {"mkreslist" , required_argument, 0, '9'},
     {"debug" ,     no_argument, 0, 'X'},
     {0, 0, 0, 0}
    };

    int option_index = 0;
    c = getopt_long_only (argc, argv, "hivVF:f:M:m:S:d:p:s:e:b:t:0:1:2:34:r:E:j:k:X",
                     long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
    {
     case 0:
       fprintf ( stderr, "Curious, there're some problems with getopt_long_only\n");
       exit(0);

     case 'V':
       OPT.VERBOSE_FLAG = 1;
       verbose = 1;
       break;

     case 'h':
       Print_Help    ( &OPT );
       exit          (99);
       break;

     case 'i':
       OPT.INFO_FLAG = 1;
       break;

     case 'F':
       OPT.FRL_FLAG = 1;
       OPT.FRL_FILE = optarg ;
       break;

     case 'f':
       fprintf( stderr, "Warning : the -f flag is being phased out.\nPlease use the -F option with the frame number you wish to extract\n");
       OPT.FRN_FLAG = 1;
       OPT.FRN_N = atoi( optarg );
       break;

     case 'E':
       OPT.IE_FLAG = 1;
       OPT.IA_FILE = optarg ;
       break;

     case 'b':
       OPT.BEG_FLAG = 1;
       OPT.TRJN_BEG  = atoi( optarg ) ;
       break;

     case 'e':
       OPT.END_FLAG = 1;
       OPT.TRJN_END  = atoi( optarg ) ;
       break;

     case 'S':
       OPT.TRJS_FLAG = 1;
       OPT.TRJS_FILE = optarg ;
       break;

     case '4':
       OPT.TRJC_FLAG = 1;
       OPT.TRJC_FILE = optarg ;
       break;

     case 'T':
       OPT.ITRJ_FLAG = 1;
       OPT.ITRJ_FILE = optarg ;
       OPT.ITRJ_TYPE = wrd_whichtrjtype(optarg);
       break;

     case 't':
       OPT.OTRJ_FLAG = 1;
       OPT.OTRJ_FILE = optarg ;
       OPT.OTRJ_TYPE = wrd_whichtrjtype(optarg);
       break;

     case 'M':                  // input structure file
       OPT.IMOL_FLAG = 1;
       OPT.IMOL_FILE = optarg ;
       OPT.IMOL_TYPE = wrd_whichmoltype(optarg);
       break;

     case 'm':                  // output structure file
       OPT.OMOL_FLAG = 1;
       OPT.OMOL_FILE = optarg ;
       OPT.OMOL_TYPE = wrd_whichmoltype(optarg);
       break;

     case 'r':
       OPT.ATRJ_FLAG = 1;
       OPT.ATRJ_FILE = optarg ;
       OPT.ATRJ_TYPE = wrd_whichtrjtype(optarg);
       break;

     case 'l':                  // append a structure file
       OPT.AMOL_FLAG = 1;
       OPT.AMOL_FILE = optarg ;
       OPT.AMOL_TYPE = wrd_whichmoltype(optarg);
       break;

     case 'c':
       OPT.HEAD_FLAG = 1;
       OPT.ITRJ_FILE = optarg ;
       OPT.ITRJ_TYPE = wrd_whichtrjtype(optarg);
       break;

     case 'g':
       OPT.MONO_FLAG = 1;
       break;

     case 'p':
       OPT.THR_FLAG = 1;
       OPT.THR_NUM = atoi( optarg );
       break;

     case 'n':
       OPT.CONV_FLAG = 1;
       break;

     case  'v':
       OPT.AVG_FLAG  = 1;
       break;

     case 'k':
       OPT.SKIP_FLAG = 1;
       OPT.SKIP_STEP = atoi( optarg );
       break;

     case 'A':
       OPT.IA_FLAG = 1;
       OPT.IA_FILE = optarg ;
       break;

     case 'o':
       OPT.OA_FLAG = 1;
       OPT.OA_FILE = optarg ;
       break;

     case 'd':
       OPT.MOD_FLAG = 1;
       OPT.MOD_STRING = optarg ;
       break;

     case 's':
       OPT.SELE_FLAG = 1;
       OPT.SELE_STRING = optarg ;
       break;

      case 'j':
       OPT.CSELE_FLAG = 1;
       OPT.SELE_STRING = optarg ;
       break;

     case '0':
       OPT.PBC_FLAG = 1;
       OPT.PBCBOX[0] =  atof( optarg );
       if (! OPT.PBCBOX[1])
        OPT.PBCBOX[1] =  atof( optarg );
       if (! OPT.PBCBOX[2])
        OPT.PBCBOX[2] =  atof( optarg );
       break;

     case '1':
       OPT.PBC_FLAG = 1;
       OPT.PBCBOX[1] =  atof( optarg );
       break;

     case '2':
       OPT.PBC_FLAG = 1;
       OPT.PBCBOX[2] =  atof( optarg );
       break;

     case '3':
       OPT.NOPBC_FLAG = 1;
       break;

     case 'X':
       OPT.DEBUG_FLAG = 1;
       break;

     case '9':                  // input structure file
       OPT.IMOL_FLAG = 1;
       OPT.IMOL_FILE = optarg ;
       OPT.IMOL_TYPE = wrd_whichmoltype(optarg);
       mkResList( &OPT );
       exit(0);
     
     case 'z':
       printf ("option -z with value `%s'\n", optarg);
       break;


     case '?':
       /* getopt_long would have already printed an error message. */
       if( OPT.Ia_FLAG || OPT.Ie_FLAG )
       {
         if( !strcmp( argv[optind-1], "-ia") || !strcmp( argv[optind-1], "-ie"))
         {
           OPT.Ia_STRING = argv[optind];
           break;
         }
         else if( iaoptions == 0 )    //OPT.current = NULL )
         {
           OPT.current = OPT.iaopt;
           OPT.iaopt->noptions = 1;
           iaoptions = 1;
         }
         else
         {
           OPT.current->next = malloc( sizeof( _iaOpt ));
           OPT.current = OPT.current->next;
           OPT.iaopt->noptions ++;
         }
         OPT.current->option = argv[optind-1];
         OPT.current->value = argv[optind];
       }
       break;

     default:
       fprintf( stderr, "Whoa, went right through the options\n");
       abort ();
    }
  }



  // ================================
  // === now let the action begin ===
  // ================================


  OPT.TCONV_FLAG = 
      OPT.ITRJ_FLAG && OPT.OTRJ_FLAG && OPT.CONV_FLAG;

  if( (OPT.IA_FLAG + OPT.Ia_FLAG + OPT.IE_FLAG + OPT.Ie_FLAG) > 1 )
  {
    fprintf( stderr, "IA, Ia, IE and Ie are mutually exclusive\n");
    exit(0);
  }


  if ( OPT.Ia_FLAG )
  {
    if ( OPT.IA_FLAG )
    {
      fprintf( stderr, "Cannot use both -iA and -ia !\n");
      exit(0);
    }
    if ( (!OPT.ITRJ_FLAG ) || (!OPT.IMOL_FLAG) )
    {
      fprintf ( stderr, "Something's missing: -itrj (or -E) and -imol are required!\a\n" );
      exit ( 99 );
    }

    superCalc( &OPT );
    exit(0);
  }

  
  if ( OPT.IA_FLAG == 1 )
  {
    if ( (!OPT.ITRJ_FLAG ) || (!OPT.IMOL_FLAG) )
    {
      fprintf ( stderr, "Something's missing: -itrj and -imol are required!\a\n" );
      exit ( 99 );
    }
//  iA_Calc( &OPT );
    superCalc( &OPT );
    exit (0);
  }

  if ( OPT.IE_FLAG || OPT.Ie_FLAG )
  {
    Calc2( &OPT );
    exit (0);
  }

  if ( OPT.FRN_FLAG == 1)
  {
    if ( (!OPT.ITRJ_FLAG) || (!OPT.IMOL_FLAG)  || (!OPT.OMOL_FLAG) )
    {
      printf ( "Something's missing: an input trajectory and structure and an output structure are required!\a\n" );
      exit ( 99 );
    }

    MolExtract( &OPT );
    exit (0);
  }

  if ( OPT.FRL_FLAG == 1)
  {
    if ( OPT.IA_FLAG != 1 )     // we're not dealing with an analysis run
    {
      if ( !OPT.ITRJ_FLAG )
      {
        printf ( "Something's missing: -itrj required!\a\n" );
        exit ( 99 );
      }
      if ( OPT.OTRJ_FLAG )
      {
        TrjExtract( &OPT );
        exit (0);
      }
      else if ( OPT.OMOL_FLAG )
      {
        if ( !OPT.IMOL_FLAG )
        {
          printf ( "Something's missing: -imol required!\a\n" );
          exit ( 99 );
        }
        if(wrd_isnumber(OPT.FRL_FILE))
        {
          OPT.FRN_FLAG = 1;
          OPT.FRN_N = atoi(OPT.FRL_FILE);
          MolExtract( &OPT );
        }
        else
          MultiMolExtract ( &OPT );
        exit (0);
      }
      else
      {
        fprintf( stderr, "Something's missing: -omol or -otrj required!");
      }
    }
  }

  if ( OPT.ITRJ_FLAG && OPT.OTRJ_FLAG && !(OPT.IMOL_FLAG || OPT.OMOL_FLAG || OPT.TCONV_FLAG))
   {
     fprintf( stdout, "Merging trajectories\n");
     Trj_Merge( &OPT );
     exit (0);
   }

  if ( OPT.TRJS_FLAG == 1)
   {
     if ( !OPT.OTRJ_FLAG || !OPT.SELE_FLAG )
      {
         printf ( "Something's missing: -otrj and -sele required!\a\n" );
         exit ( 99 );
      }
     Trj_sum( &OPT );
     exit (0);
   }

  if ( OPT.TRJC_FLAG == 1)
   {
     if ( !OPT.OTRJ_FLAG )
      {
         printf ( "Something's missing: -otrj required!\a\n" );
         exit ( 99 );
      }
     Trj_combine( &OPT );
     exit (0);
   }

  if ( OPT.ATRJ_FLAG == 1)
   {
     if ( !OPT.OTRJ_FLAG )
      {
         printf ( "Something's missing: -otrj required!\a\n" );
         exit ( 99 );
      }
     Trj_append( &OPT );
     exit (0);
   }

  if ( OPT.AMOL_FLAG == 1 )
   {
     if ( !OPT.OTRJ_FLAG )
      {
         printf ( "Something's missing: -otrj required!\a\n" );
         exit ( 99 );
      }
     Mol_append( &OPT );
     exit (0);
   }

 #ifdef GMX
  if ( OPT.TCONV_FLAG == 1 )
  {
   if( !OPT.ITRJ_FLAG || !OPT.OTRJ_FLAG )
   {
    printf("Something's missing: -itrj) and -otrj are required!\a\n" );
    exit(99);
   }
   
   if ( OPT.ITRJ_TYPE == xtc)
   {
    if ( OPT.OTRJ_TYPE == xtc ) 
    {
     printf ( "Something's missing: -otrj should not be a xtc (again)!\a\n" );
     exit ( 99 );
    }
    Xtc2Dcd( &OPT );
    exit (0);
   }
   
   if ( OPT.ITRJ_TYPE == dcd)
   {
    if ( ! OPT.IMOL_FLAG || OPT.OTRJ_TYPE != xtc  )
    {
      printf ( "Something's missing: -imol is required and -otrj should not be a dcd\a\n" );
      exit ( 99 );
    }
    Dcd2Xtc( &OPT );
    exit (0);
   }
   
  }
 #endif
  
  if ( OPT.CONV_FLAG == 1 && OPT.ITRJ_FLAG && OPT.OMOL_TYPE == 4)
   {
     if ( !OPT.IMOL_FLAG || !OPT.ITRJ_FLAG )
     {
       printf ( "Something's missing: input trajectory and input reference structure file required!\a\n" );
       exit ( 99 );
     }
     fprintf( stderr, "Calling xyzextract\n");
     xyzExtract( &OPT );
     exit (0);
   }

  if ( OPT.CONV_FLAG == 1)
   {
     if (!OPT.IMOL_FLAG || !OPT.OMOL_FLAG)
      {
         printf ( "Something's missing: -imol and -omol required!\a\n" );
         exit ( 99 );
      }
     MolConv( &OPT );
     exit (0);
   }
  
  if ( OPT.HEAD_FLAG == 1)
   {
     TrjHeader_Print ( &OPT );
     exit(0);
   }
  
  if ( OPT.CSELE_FLAG == 1)
   {
     if ( !OPT.IMOL_FLAG )
      {
         printf ( "Something's missing: a structure is required!\a\n" );
         exit ( 99 );
      }
     CheckSele( &OPT );
     exit (0);
   }

  if ( OPT.MOD_FLAG == 1)
   {
     if ( !OPT.ITRJ_FLAG )
      {
         printf ( "Something's missing: -itrj!\a\n" );
         exit ( 99 );
      }
     Head_Mod( &OPT );
     exit (0);
   }

  if ( OPT.MONO_FLAG == 1)
   {
     if ( (!OPT.IMOL_FLAG) || !OPT.OTRJ_FLAG )
      {
         printf ( "Something's missing: -imol and -otrj required!\a\n" );
         exit ( 99 );
      }
     MkMono( &OPT );
     exit (0);
   }
  
  if ( OPT.AVG_FLAG == 1)
   {
     if ( (!OPT.IMOL_FLAG) || (!OPT.ITRJ_FLAG) || (!OPT.OMOL_FLAG) )
      {
         printf ( "Something's missing: -imol, -itrj and -omol required!\a\n" );
         exit ( 99 );
      }
     AvgExtract( &OPT );
     exit (0);
   }
  
  if ( OPT.Ia_FLAG == 1)
  {
    if ( (!OPT.ITRJ_FLAG) || (!OPT.IMOL_FLAG) )
    {
       printf ( "Something's missing: -itrj (or -E) and -imol are required!\a\n" );
       exit ( 99 );
    }
    iA_Calc( &OPT );
    exit (0);
  }
  
  if( OPT.DEBUG_FLAG == 1 || OPT.IMOL_FLAG == 1 )
  {
    fprintf( stdout, "Debug function called\n");
    if ( !OPT.IMOL_FLAG )
    {
      fprintf( stderr, "Input mol file (-imol) needed!\n" );
      exit(0);
    }
    MolInfo_Print( OPT.IMOL_FILE );
    exit(0);
  }
  
  Print_About ();
  exit (0);

}
void setOPT(struct sopt *OPT)
{
       OPT->DEBUG_FLAG =  0;

         OPT->FRN_FLAG =  0;
            OPT->FRN_N =  0;

        OPT->HEAD_FLAG =  0;

        OPT->ITRJ_FLAG =  0;
        OPT->OTRJ_FLAG =  0;
        OPT->ATRJ_FLAG =  0;
        OPT->IMOL_FLAG =  0;
        OPT->OMOL_FLAG =  0;
        OPT->AMOL_FLAG =  0;
         OPT->FRL_FLAG =  0;
        OPT->T2Ms_FLAG =  0;
        OPT->MONO_FLAG =  0;
       OPT->TCONV_FLAG =  0;
         OPT->BEG_FLAG =  0;
         OPT->END_FLAG =  0;
        OPT->TRJL_FLAG =  0;
        OPT->TRJS_FLAG =  0;
        OPT->TRJC_FLAG =  0;
        OPT->SKIP_FLAG =  0;
          OPT->IA_FLAG =  0;
          OPT->OA_FLAG =  0;
          OPT->Ia_FLAG =  0;
          OPT->IE_FLAG =  0;
          OPT->Ie_FLAG =  0;
        OPT->CONV_FLAG =  0;
         OPT->PBC_FLAG =  0;
       OPT->NOPBC_FLAG =  0;
      OPT->ATMNUM_FLAG =  0;
        OPT->HELP_FLAG =  0;
     OPT->VERBOSE_FLAG =  0;
         OPT->MOD_FLAG =  0;
         OPT->THR_FLAG =  0;
         OPT->AVG_FLAG =  0;
     OPT->MERGEAN_FLAG =  0;
        OPT->SELE_FLAG =  0;
       OPT->CSELE_FLAG =  0;
        OPT->INFO_FLAG =  0;

   	   OPT->ITRJ_FILE =  NULL;
   	   OPT->OTRJ_FILE =  NULL;
       OPT->ATRJ_FILE =  NULL;
       OPT->IMOL_FILE =  NULL;
       OPT->OMOL_FILE =  NULL;
       OPT->AMOL_FILE =  NULL;
   	   OPT->TRJL_FILE =  NULL;
   	   OPT->TRJS_FILE =  NULL;
   	   OPT->TRJC_FILE =  NULL;
   	    OPT->FRL_FILE =  NULL;
   	   OPT->T2Ms_FILE =  NULL;
   	     OPT->IA_FILE =  NULL;
   	   OPT->Ia_STRING =  NULL;
   	     OPT->OA_FILE =  NULL;
     OPT->ATMNUM_FILE =  NULL;
   	    OPT->HELP_ARG =  NULL;
      OPT->MOD_STRING =  NULL;
     OPT->SELE_STRING =  NULL;
    OPT->MERGEAN_FILE =  NULL;
   	   OPT->DCDN_BASE =  NULL;
  OPT->SUBSELE_STRING =  NULL;
   
       OPT->IMOL_TYPE =  99;
       OPT->OMOL_TYPE =  99;
       OPT->AMOL_TYPE =  99;
       OPT->ITRJ_TYPE =  99;
       OPT->OTRJ_TYPE =  99;
       OPT->ATRJ_TYPE =  99;
  
        OPT->TRJN_BEG =  0;
        OPT->TRJN_END =  0;

        OPT->TYPE_NUM =  0; 
        OPT->TYPE_ATM =  NULL;

       OPT->SKIP_STEP =  1;
        OPT->BEG_STEP =  0;
        OPT->END_STEP =  0;
         OPT->THR_NUM =  0;
       OPT->PBCBOX[0] =  0.0;
       OPT->PBCBOX[1] =  0.0;
       OPT->PBCBOX[2] =  0.0;
   
           OPT->iaopt =  NULL;
         OPT->current = NULL;
   
   return;
}
// ---------------------------------------------------------------------
void copyOPT( struct sopt *OPTsource, struct sopt *OPTdest )
{
       OPTdest->DEBUG_FLAG =     OPTsource->DEBUG_FLAG;
                          
         OPTdest->FRN_FLAG =       OPTsource->FRN_FLAG;
            OPTdest->FRN_N =          OPTsource->FRN_N;
                          
        OPTdest->HEAD_FLAG =      OPTsource->HEAD_FLAG;
                          
        OPTdest->ITRJ_FLAG =      OPTsource->ITRJ_FLAG;
        OPTdest->OTRJ_FLAG =      OPTsource->OTRJ_FLAG;
        OPTdest->ATRJ_FLAG =      OPTsource->ATRJ_FLAG;
        OPTdest->IMOL_FLAG =      OPTsource->IMOL_FLAG;
        OPTdest->OMOL_FLAG =      OPTsource->OMOL_FLAG;
        OPTdest->AMOL_FLAG =      OPTsource->AMOL_FLAG;
         OPTdest->FRL_FLAG =       OPTsource->FRL_FLAG;
        OPTdest->T2Ms_FLAG =      OPTsource->T2Ms_FLAG;
        OPTdest->MONO_FLAG =      OPTsource->MONO_FLAG;
       OPTdest->TCONV_FLAG =     OPTsource->TCONV_FLAG;
         OPTdest->BEG_FLAG =       OPTsource->BEG_FLAG;
         OPTdest->END_FLAG =       OPTsource->END_FLAG;
        OPTdest->TRJL_FLAG =      OPTsource->TRJL_FLAG;
        OPTdest->TRJS_FLAG =      OPTsource->TRJS_FLAG;
        OPTdest->TRJC_FLAG =      OPTsource->TRJC_FLAG;
        OPTdest->SKIP_FLAG =      OPTsource->SKIP_FLAG;
          OPTdest->IA_FLAG =        OPTsource->IA_FLAG;
          OPTdest->OA_FLAG =        OPTsource->OA_FLAG;
          OPTdest->Ia_FLAG =        OPTsource->Ia_FLAG;
          OPTdest->IE_FLAG =        OPTsource->IE_FLAG;
          OPTdest->Ie_FLAG =        OPTsource->Ie_FLAG;
        OPTdest->CONV_FLAG =      OPTsource->CONV_FLAG;
         OPTdest->PBC_FLAG =       OPTsource->PBC_FLAG;
       OPTdest->NOPBC_FLAG =     OPTsource->NOPBC_FLAG;
      OPTdest->ATMNUM_FLAG =    OPTsource->ATMNUM_FLAG;
        OPTdest->HELP_FLAG =      OPTsource->HELP_FLAG;
     OPTdest->VERBOSE_FLAG =   OPTsource->VERBOSE_FLAG;
         OPTdest->MOD_FLAG =       OPTsource->MOD_FLAG;
         OPTdest->THR_FLAG =       OPTsource->THR_FLAG;
         OPTdest->AVG_FLAG =       OPTsource->AVG_FLAG;
     OPTdest->MERGEAN_FLAG =   OPTsource->MERGEAN_FLAG;
        OPTdest->SELE_FLAG =      OPTsource->SELE_FLAG;
       OPTdest->CSELE_FLAG =     OPTsource->CSELE_FLAG;
        OPTdest->INFO_FLAG =      OPTsource->INFO_FLAG;

   	   OPTdest->ITRJ_FILE =   	   OPTsource->ITRJ_FILE;
   	   OPTdest->OTRJ_FILE =   	   OPTsource->OTRJ_FILE;
       OPTdest->ATRJ_FILE =       OPTsource->ATRJ_FILE;
       OPTdest->IMOL_FILE =       OPTsource->IMOL_FILE;
       OPTdest->OMOL_FILE =       OPTsource->OMOL_FILE;
       OPTdest->AMOL_FILE =       OPTsource->AMOL_FILE;
   	   OPTdest->TRJL_FILE =   	   OPTsource->TRJL_FILE;
   	   OPTdest->TRJS_FILE =   	   OPTsource->TRJS_FILE;
   	   OPTdest->TRJC_FILE =   	   OPTsource->TRJC_FILE;
   	    OPTdest->FRL_FILE =   	    OPTsource->FRL_FILE;
   	   OPTdest->T2Ms_FILE =   	   OPTsource->T2Ms_FILE;
   	     OPTdest->IA_FILE =   	     OPTsource->IA_FILE;
   	   OPTdest->Ia_STRING =   	   OPTsource->Ia_STRING;
   	     OPTdest->OA_FILE =   	     OPTsource->OA_FILE;
     OPTdest->ATMNUM_FILE =     OPTsource->ATMNUM_FILE;
   	    OPTdest->HELP_ARG =   	    OPTsource->HELP_ARG;
      OPTdest->MOD_STRING =      OPTsource->MOD_STRING;
     OPTdest->SELE_STRING =     OPTsource->SELE_STRING;
    OPTdest->MERGEAN_FILE =    OPTsource->MERGEAN_FILE;
   	   OPTdest->DCDN_BASE =   	   OPTsource->DCDN_BASE;
  OPTdest->SUBSELE_STRING =  OPTsource->SUBSELE_STRING;
   
       OPTdest->IMOL_TYPE =  OPTsource->IMOL_TYPE;
       OPTdest->OMOL_TYPE =  OPTsource->OMOL_TYPE;
       OPTdest->AMOL_TYPE =  OPTsource->AMOL_TYPE;
       OPTdest->ITRJ_TYPE =  OPTsource->ITRJ_TYPE;
       OPTdest->OTRJ_TYPE =  OPTsource->OTRJ_TYPE;
       OPTdest->ATRJ_TYPE =  OPTsource->ATRJ_TYPE;
                         
        OPTdest->TRJN_BEG =   OPTsource->TRJN_BEG;
        OPTdest->TRJN_END =   OPTsource->TRJN_END;
                         
        OPTdest->TYPE_NUM =   OPTsource->TYPE_NUM; 
        OPTdest->TYPE_ATM =   OPTsource->TYPE_ATM;
                         
       OPTdest->SKIP_STEP =  OPTsource->SKIP_STEP;
        OPTdest->BEG_STEP =   OPTsource->BEG_STEP;
        OPTdest->END_STEP =   OPTsource->END_STEP;
         OPTdest->THR_NUM =    OPTsource->THR_NUM;
       OPTdest->PBCBOX[0] =  OPTsource->PBCBOX[0];
       OPTdest->PBCBOX[1] =  OPTsource->PBCBOX[1];
       OPTdest->PBCBOX[2] =  OPTsource->PBCBOX[2];
                         
           OPTdest->iaopt =   OPTsource->iaopt;
         OPTdest->current =   OPTsource->current;
   
   return;
}
