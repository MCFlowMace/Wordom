// ------------------------------------------------------------------
// Copyright (C) 2003-2009  Univ. of Zurich and Univ. of Modena & Reggio Emilia
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

/*! \file fileio.c
 \brief the input/output functions and structures
 
 fileio provides the low-level input/output functions
 to access data in structure (crd/pdb) and trajectory(dcd/xtc) files
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <unistd.h>
#include <fnmatch.h>
#include <ctype.h>
/*#ifdef MACOSX
 #include "/usr/include/malloc/malloc.h"
#else
 #include <malloc.h>
#endif*/

#include "fileio.h"
#include "tools.h"

/*! \def _FILE_OFFSET_BITS
 these define are used to give big file (>4GB) support
*/
#define _FILE_OFFSET_BITS 64    /*!< \brief used to give big file (>4GB) support */
#define MAX_LENGTH_SUFFIX 4
#define O_GMX_NOT_COMPILED fprintf(stderr, "Sorry, this version was not compiled with XTC files support\n");

#define PIE 3.14159265358979323846

// defines for fnm_match
/* Bits set in the FLAGS argument to `fnmatch'.  */
#define	FNM_PATHNAME	(1 << 0) /* No wildcard can ever match `/'.  */
#define	FNM_NOESCAPE	(1 << 1) /* Backslashes don't quote special chars.  */
#define	FNM_PERIOD	(1 << 2) /* Leading `.' is matched only explicitly.  */

# define FNM_FILE_NAME	 FNM_PATHNAME	/* Preferred GNU name.  */
# define FNM_LEADING_DIR (1 << 3)	/* Ignore `/...' after a match.  */
# define FNM_CASEFOLD	 (1 << 4)	/* Compare without regard to case.  */
# define FNM_EXTMATCH	 (1 << 5)	/* Use ksh-like extended matching. */

/* Value returned by `fnmatch' if STRING does not match PATTERN.  */
#define	FNM_NOMATCH	1

short int verbose = 0;

/* ------------------------------------------------------------------ */
int isincharlist(char *list, int nitems, char item )
{
  int ii;
  
  for( ii=0; ii<nitems; ii++ )
  {
    if( list[ii] == item )
      return ii;
    else if( list[ii] == '\0' )
      return -1;
  }
  return -1;
}
/* ------------------------------------------------------------------ */
int isinstrlist(char **list, int nitems, char *item )
{
  int ii;
  
  for( ii=0; ii<nitems; ii++ )
  {
    if( !strcmp( list[ii], item) )
      return ii;
    else if( list[ii][0] == '\0' )
      return -1;
  }
  return -1;
}
/* ------------------------------------------------------------------ */
FILE *O_File ( char *file, char *mode )
{
   FILE *stream;
   char *func="O_File";
   if (( stream = fopen ( file, mode ) ) == NULL)
    {
      fprintf ( stderr,"%s: (EE) %s: \"%s\"\n", func, strerror( errno ), file );

      exit( 99 );
    }
   return ( stream );
}
/* ------------------------------------------------------------------ */
char *swap(int n,char *byte)
{
   int i,m;
   char c;

   switch (n)
    {
     case 2:
       c=byte[0]; 
       byte[0]=byte[1]; 
       byte[1]=c;
       break;
     case 4:
       c=byte[0]; 
       byte[0]=byte[3]; 
       byte[3]=c;
       c=byte[1]; 
       byte[1]=byte[2]; 
       byte[2]=c;
       break;
     case 8:
       c=byte[0]; 
       byte[0]=byte[7]; 
       byte[7]=c;
       c=byte[1]; 
       byte[1]=byte[6]; 
       byte[6]=c;
       c=byte[2]; 
       byte[2]=byte[5]; 
       byte[5]=c;
       c=byte[3]; 
       byte[3]=byte[4]; 
       byte[4]=c;
       break;
     default:
       m=n/2;
       for (i=0;i<m;i++)
        {
          c=byte[i];
          byte[i]=byte[n-i-1];
          byte[n-i-1]=c;
        }
       break;
    }
   return(byte);
}
/* ------------------------------------------------------------------ */

/* Get the suffix (or extension) of a file name */
int 
wrd_getsuffix(char *file_name, char **suffix)
{
    char *loc;
    size_t len;
   
    /* Get last position of char '.' */
    loc = strrchr(file_name, '.');

    /* If we don't find it, print an error and exit */
    if (loc==NULL) {
        fprintf( stderr, "Troubles to get the suffix of %s\n", file_name);
        exit(0);
    } 

    /* Otherwise, allocate memory in suffix and copy the suffix in it */
    else {
        loc++;
        len = strlen(loc);
        *suffix = malloc((len+1)*sizeof(char));
        if (*suffix!=NULL) {
            strncpy(*suffix, loc, len);
            (*suffix)[len]='\0';
        } else {
            fprintf(stderr, "Error allocating memory\n");
            exit(0);
        }
    }

    return 0;
}


/* Check if the supplied filename is a DAT of TXT file */
int 
wrd_isdattxt(char *filename)
{
    int     ii;
    char   *suffix=NULL;
    int     islist=0;
   
    /* Get the suffix */
    wrd_getsuffix(filename, &suffix);
    if (suffix == NULL) {
        fprintf(stderr, "File %s has no suffix: could not determine filetype\n", filename);
        exit(0);
    }
   
    /* Move to uppercase */
    for (ii=0; ii<strlen(suffix); ii++ ) {
        suffix[ii] = toupper(suffix[ii]);
    }
   
    /* Check if it a TXT or DAT file */
    if (strncmp(suffix, "TXT", 3)==0) {
        islist = 1;
    } else if (strncmp(suffix, "DAT", 3)==0) {
        islist = 1;
    }
   
    /* Free memory */
    free(suffix);
   
    return islist;
}


/* ------------------------------------------------------------------ */

/* Guess the file format based on the extension of the filename for molecule files */
moltype 
wrd_whichmoltype(char *filename)
{
    int          ii;
    char         *suffix=NULL;
    moltype      type=0;
   
    /* Get the suffix */
    wrd_getsuffix(filename, &suffix);
    if (suffix == NULL) {
        fprintf(stderr, "File %s has no suffix: could not determine filetype\n", filename);
        exit(0);
    }

    /* Move to uppercase */
    for (ii=0; ii<strlen(suffix); ii++) {
        suffix[ii] = toupper(suffix[ii]);
    }

    /* Compare with supported filetypes */
    if (strncmp(suffix, "PDB", 3)==0) {
        type = PDB;
    } else if (strncmp(suffix, "CRD", 3)==0) {
        type = CRD;
    } else if (strncmp(suffix, "COR", 3)==0) {
        type = COR;
    } else if (strncmp(suffix, "XYZ", 3)==0) {
        type = XYZ;
    } else if (strncmp(suffix, "GRO", 3)==0) {
        type = GRO;
    } else {
        fprintf( stderr, "Unknown mol format %s \n", suffix);
        exit(0);
    }

    /* Free memory */
    free(suffix);
   
    return type;
}


/* ------------------------------------------------------------------ */
/* Guess the file format based on the extension of the filename for trajectory files */
trjtype 
wrd_whichtrjtype(char *filename)
{
    int          ii;
    char         *suffix=NULL;
    trjtype      type=0;
    
    /* Get the suffix */
    wrd_getsuffix(filename, &suffix);
    if (suffix == NULL) {
        fprintf(stderr, "File %s has no suffix: could not determine filetype\n", filename);
        exit(0);
    }

    /* Move to uppercase */
    for (ii=0; ii<strlen(suffix); ii++) {
        suffix[ii] = toupper(suffix[ii]);
    }

    /* Compare with supported filetypes */
    if (strncmp(suffix, "DCD", 3)==0 || strncmp(suffix, "TRJ", 3)==0) {
        type = dcd;
    } else if (strncmp(suffix, "XTC", 3)==0) {
        type = xtc;
    } else if (strncmp(suffix, "TXT", 3)==0) {
        type = txt;
    } else if (strncmp(suffix, "PDB", 3)==0) {
        type = pdb;
    } else if (strncmp(suffix, "CRD", 3)==0) {
        type = crd;
    } else if (strncmp(suffix, "COR", 3)==0) {
        type = cor;
    } else {
        fprintf( stderr, "Unknown trj format %s \n", suffix);
        exit(0);
    }
   
    /* Free memory */
    free(suffix);
   
   return type;
}


/* ------------------------------------------------------------------ */
int OpenTrj ( char *filename, Traj *traj , char *mode ) 
{
   trjtype      type;
   
   traj->xtc = 0;
   type = wrd_whichtrjtype( filename );
   if ( type == dcd )
   {
    traj->filetype = dcd;
    traj->file =  O_File (filename, mode);
   }
   else if ( type == xtc )
   {
    traj->filetype = xtc;
    #ifdef GMX
     traj->xtc = xdrfile_open(filename, mode);
     if (traj->xtc == NULL) {
      fprintf( stderr, "Cannot open file %s \n", filename);
      exit(99);	 
     }
     traj->previousframe = 0;
    #else
     O_GMX_NOT_COMPILED;
    #endif
   }
   else
   {
    fprintf( stderr, "Unknown file extension \"%s\"", filename);
    exit(99);
   }
   
   strcpy( traj->filename, filename );
   
   if( !strncmp( mode, "w", 1) )
     traj->hdr = malloc(sizeof(struct trjh));
   
   traj->open = 1;
   
   return 0;
}
/* ------------------------------------------------------------------ */
Traj * InitTrj ( char *filename, char *mode ) 
{
   trjtype      type;
   Traj        *traj=NULL;
   
   traj = calloc( 1, sizeof(Traj));
   if( traj == NULL )
    return NULL;
   
   traj->xtc = 0;
   type = wrd_whichtrjtype( filename );
   if ( type == dcd )
   {
    traj->filetype = dcd;
    traj->file =  O_File (filename, mode);
   }
   else if ( type == xtc )
   {
     traj->filetype = xtc;
     traj->xtc = xdrfile_open(filename, mode);
     if (traj->xtc == NULL)
     {
       fprintf( stderr, "Cannot open file %s \n", filename);
       exit(99);	 
     }
     
     traj->previousframe = 0;
     traj->hdetail = LITE;
   }
   else
   {
    fprintf( stderr, "Unknown file extension \"%s\"", filename);
    exit(99);
   }
   
   strcpy( traj->filename, filename );
   
   if( mode[0] == 'w' )
     traj->hdr = malloc(sizeof(struct trjh));
   else if( mode[0] == 'r' )
   {
     ReadTrjHeader( traj );
     traj->crdset = InitCoor( traj->hdr->nato );
   }
   
   traj->open = 1;
   
   return traj;
}
// ------------------------------------------------------------------
void CloseTrj ( Traj *traj ) 
{
   if( traj->filetype == dcd )
    fclose(traj->file);
   else if( traj->filetype == xtc )
    xdrfile_close(traj->xtc);
   memset( traj->filename, '\0', sizeof(traj->filename) );
   traj->filetype = 99;
   // free(traj->hdr); DEBUG: it shouldn't be so, but this something generates a double free error
   traj->open = 0;
   
   return;
}
// ------------------------------------------------------------------
void ReadDcdHeader ( FILE *dcd, Trjh *hdr )      // dcd file, dcd header
{
   int       ii;
   //char      buffer[10240];
   char     *buffer;
   int       buf_size = 10240;
   char      check1[4], check2[4], check3[4];
   int       nato, comnum;
   long long filesize;
   float     real_nframe;
   int       header_size;
   int       frame_size, frame1_size;
   int       sex = 0;                   // 0 for HOMO, 1 for HETERO (needs swap)
   int       garbage=0;
   
   int       field1, field2;
   
   CoorSet  *crd;
   long long offset;
   char         recl[8];
   
   buffer = calloc( buf_size, sizeof(char));
   
   hdr->sixtyfour=0;

   fseeko ( dcd, 0 , SEEK_END ); 
   filesize = ftello ( dcd );
   fseeko ( dcd, 0 , SEEK_SET );
   //printf("DEBUG: filesize: %d\n", filesize );
   if ( filesize == 0 )
    {
     printf( "*** Error: empty file! ***\n\n" );
     exit(0);
    }

   memset ( check1, '\0', sizeof(check1));
   memset ( check2, '\0', sizeof(check2));
   memset ( check3, '\0', sizeof(check3));
   
   fread ( check1, 4, 1, dcd);
   fread ( check2, 4, 1, dcd);
   fread ( check3, 4, 1, dcd);
   
   field1 = *(int *)(check1+0);
   field2 = *(int *)(check2+0);
   
   if ( field1==84 && (!strncmp( check2, "CORD", 4) || !strncmp( check2, "VELD", 4)) )
   {
    sex=0;
    hdr->sixtyfour=0;
   }
   else if ( field1==84 && field2==0 && (!strncmp( check3, "CORD", 4) || !strncmp( check3, "VELD", 4)) )
   {
    sex=0;
    hdr->sixtyfour=1;
   }
   else if ( field1==1409286144 && (!strncmp( check2, "CORD", 4) || !strncmp( check2, "VELD", 4)) )
   {
    sex=1;
    hdr->sixtyfour=0;
   }
   else if ( field1==0 && field2==1409286144 && (!strncmp( check3, "CORD", 4) || !strncmp( check3, "VELD", 4)) )
   {
    sex=1;
    hdr->sixtyfour=1;
   }
   else
   {
    fprintf( stderr, "*** This is not a proper dcd file ***\n\n");
    exit(0);
   }

   memset ( buffer, '\0', buf_size);

   fseeko ( dcd, 4*(1+hdr->sixtyfour) , SEEK_SET );

   fread(hdr->type, 1, 4, dcd);
   hdr->type[4]='\0';

   fread(buffer, 1, 80, dcd);
        
//  swap if different sex (byteorder)
   if ( sex )
    {
     swap ( 4, buffer+0  );
     swap ( 4, buffer+4  );
     swap ( 4, buffer+8  );
     swap ( 4, buffer+12 );
     swap ( 4, buffer+32 );
     swap ( 4, buffer+36 );
     swap ( 4, buffer+40 );
    }

 
   hdr->c_nframe =  *(int *) (buffer+ 0);
   hdr->begstp   =  *(int *) (buffer+ 4);
   hdr->skpstp   =  *(int *) (buffer+ 8);
   hdr->c_numstp =  *(int *) (buffer+12);
   hdr->fixatm   =  *(int *) (buffer+32);
   hdr->tmstp    = *(float *)(buffer+36);
   hdr->varpbc   =  *(int *) (buffer+40);

   if( hdr->varpbc != 0 && hdr->varpbc != 1)
   {
     fprintf( stderr, "varpbc field is wrong (%d)\n", hdr->varpbc );
     exit(0);
   }
   //printf("DEBUG: varpbc: %d\n", hdr->varpbc );
   fread(buffer, 1, 2*4*(1+hdr->sixtyfour), dcd);

   fread(buffer, 1, 4, dcd);
   if ( sex)
    swap ( 4, buffer);
   
   comnum = *(int *)(buffer);
   
   if ( comnum > 12 )
    {
     printf (" Whoa, verbose comments(%d lines)! I'm only retaining the first 12 lines!\n", comnum);
     garbage = comnum-12;
     comnum=12;
    }
   for(ii=0; ii<comnum; ii++)
    {
     fread(hdr->comm[ii], 1, 80, dcd);
     hdr->comm[ii][80]='\0';
    }
   for(ii=0; ii<garbage; ii++)
    {
     fread(buffer, 1, 80, dcd);
    }

   fread(buffer, 1, 2*4*(1+hdr->sixtyfour), dcd);
   
   fread(buffer, 1, 4, dcd);
   
   if ( sex )
    swap ( 4, buffer);

   nato = *(int *)(buffer);

   hdr->nato   =   nato;
   hdr->comnum = comnum;

 /* End of headers */
 
   hdr->sex = sex;
 
 /* Check number of frames */
 
   fseeko ( dcd, 0 , SEEK_END );
 
   if ( !hdr->fixatm)
    {
      header_size = ( 6*4*(1+hdr->sixtyfour) + 23*4 + comnum*80 );
      frame_size  = ( 6*4*(1+hdr->sixtyfour) + 3*4*nato + 56*(hdr->varpbc) );
      real_nframe = ( (filesize - header_size ) / frame_size );
    }

   if ( hdr->fixatm )
    {
      header_size = ( 6*4*(1+hdr->sixtyfour) + 23*4 +  comnum*80 + 
                 (2*4*(1+hdr->sixtyfour) + 4*(nato - hdr->fixatm)) );
      frame1_size = ( 6*4*(1+hdr->sixtyfour) + 3*4*nato + 56*(hdr->varpbc) );
      frame_size  = ( 6*4*(1+hdr->sixtyfour) + 3*4*(nato - hdr->fixatm) + 56*(hdr->varpbc) );
      real_nframe = ( (filesize - header_size - frame1_size) / frame_size ) + 1;
    }
 
   hdr->nframe = real_nframe;
   hdr->numstp = real_nframe * hdr->skpstp;

   hdr->size = header_size;
   
   /* if fixatm, store coords of fixed atoms */
   if( hdr->fixatm )
   {
     hdr->fixatm_list = calloc( hdr->fixatm, sizeof(int) );
     hdr->fixcoor     = calloc( hdr->fixatm, sizeof(float));
     // first, get the atom numbers for the moving atoms
     offset = (long long) 6*4*(1+hdr->sixtyfour) + 23*4 + hdr->comnum*80 + 4*(1+hdr->sixtyfour);
     fseeko(dcd, offset, SEEK_SET);
     if ( !hdr->sex )
       fread ( hdr->fixatm_list, 4, hdr->fixatm, dcd);
     if ( hdr->sex )
     {
       fread ( buffer, 1, hdr->fixatm*4, dcd);
       for ( ii=0; ii<hdr->fixatm; ii++)
       {
         swap ( 4, buffer + (ii*4));
         hdr->fixatm_list[ii] = *(int *) (buffer + (ii*4));
       }
     }
    
    // now, get the coordinates of ALL atoms (first frame)
     crd = InitCoor( hdr->nato );
     offset = (long long) 6*4*(1+hdr->sixtyfour) + 23*4 + hdr->comnum*80 + 
             2*4*(1+hdr->sixtyfour) + hdr->fixatm*4;
     fseeko(dcd, offset, SEEK_SET);
     if ( !hdr->sex )
     {
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
       fread(crd->xcoor, 4, hdr->nato, dcd);
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
       fread(crd->ycoor, 4, hdr->nato, dcd);
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
       fread(crd->zcoor, 4, hdr->nato, dcd);
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
     }
    
     if ( hdr->sex )
     {
       free(buffer);
       buffer = calloc ( hdr->nato * 4, sizeof(char) );

       // buffer the x-coor, swap the byte order and get them
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
       fread( buffer, 1, hdr->nato*4, dcd);
       for ( ii = 0; ii < hdr->nato; ii++)
       {
         swap ( 4, buffer + ( ii*4 ));
         crd->xcoor[ii] = * (float *) (buffer + ii*4);
       }
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);

       // same for the y-coor block
       fread( buffer, 1, hdr->nato*4, dcd);
       for ( ii = 0; ii < hdr->nato; ii++)
       {
         swap ( 4, buffer + ( ii*4 ));
         crd->ycoor[ii] = * (float *) (buffer + ii*4);
       }
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);

       // last but not least, the z-coordinates block
       fread( buffer, 1, hdr->nato*4, dcd);
       for ( ii = 0; ii < hdr->nato; ii++)
       {
         swap ( 4, buffer + ( ii*4 ));
         crd->zcoor[ii] = * (float *) (buffer + ii*4);
       }
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
       
       free ( buffer );
     }
   }
   
   return;
}
// ------------------------------------------------------------------
void ReadTrjHeader ( Traj *traj)
{
   char *func="ReadTrjHeader";
   
   traj->hdr = calloc( 1, sizeof( struct trjh));
   switch ( traj->filetype )
   {
    case 0:                     // read a dcd
     fprintf( stderr, ".ooo trj type reading not implemented yet\n");
     exit(0);
    case 1:                     // read a dcd
     ReadDcdHeader ( traj->file, traj->hdr );
     break;
    case 2:                     // read a xtc
     if (traj->hdetail == LITE) 
       ReadXtcHeaderLITE ( traj );
     else if (traj->hdetail == FULL)
       ReadXtcHeaderFULL ( traj );
     else
     {
      traj->hdetail = FULL;
      ReadXtcHeaderFULL ( traj);
     }
     break;
    default:                     // read a what?
     fprintf( stderr, "Error while calling %s! check filetype\n", func);
     fprintf( stderr, "this trj type reading not implemented yet\n");
     exit(78);
     break;
   }
   
   return;

}
// ------------------------------------------------------------------
void ShowTrjHeader (const Traj *traj)
{
   char *func="ReadTrjHeader";
   
   switch ( traj->filetype )
   {
    case 0:                     // read a dcd
     fprintf( stderr, ".ooo trj type reading not implemented yet\n");
     exit(0);
     break;
    case 1:                     // read a dcd
     ShowDcdHeader (traj->hdr);
     break;
    case 2:                     // read a xtc
     ShowXtcHeader (traj->hdr);
     break;
    case 3:                     // read a what?
     fprintf( stderr, "this trj type reading not implemented yet\n");
     exit(0);
     break;
    default:
     fprintf( stderr, "Error while calling %s! check filetype\n", func);
     exit(78);
     break;
   }   
   return;
}
// ------------------------------------------------------------------
void ShowDcdHeader( Trjh *hdr )         // dcd header
{
   float tmstp;
   int       i;
   // convert timestep from CHARMM units in fs
   tmstp = floor ( (hdr -> tmstp * 488.8780249) + 0.5)/10;
   
   printf ( "Type of DCD\t\t%s\n",      hdr->type      );
   printf ( "Number frames\t%-20d\t",   hdr->c_nframe  );
   printf ( "Real Number frames\t%d\n", hdr->nframe    );
   printf ( "Begin  step\t%-20d\t",     hdr->begstp    );
   printf ( "Skip   step\t\t%d\n",      hdr->skpstp    );
   printf ( "Number steps \t%-20u\t",   hdr->numstp    );
   printf ( "Time   step (fs)\t%2.1f\n",     tmstp     );
   printf ( "Number atoms\t%-20d\t",    hdr->nato      );
   printf ( "Fixed atoms\t\t%d\n",      hdr->fixatm    );
   printf ( "Variable PBC\t%d\t\t\t",   hdr->varpbc    );
   printf ( "Made on 64-bit\t\t%d\n",   hdr->sixtyfour );
   printf ( "Comments (%d lines):\n",   hdr->comnum    );
   for(i=0;i<hdr->comnum;i++)
    printf ( "%80s"        ,      hdr->comm[i]   );
   printf ( "END\n\n" );
   
   return;
}
// ------------------------------------------------------------------
void FillDcdHeader ( Molecule   *molecule,    // pdb structure
                     Trjh       *hdr      )   // dcd header
{
   sprintf(hdr->type, "CORD");
   hdr->nframe  = 1;
   hdr->begstp  = 1000;
   hdr->skpstp  = 1000;
   hdr->numstp  = 1;
   hdr->tmstp   = 0.002;
   hdr->varpbc  = 0;
   hdr->nato    = molecule->nato;
   hdr->fixatm  = 0;
   hdr->fixatm_list = NULL;
   hdr->fixcoor     = NULL;
   hdr->sex     = 0;
   hdr->size    = 0;  /* to be computed later */
   hdr->sixtyfour = 0;
   hdr->c_nframe  = 1;
   hdr->c_numstp  = 1;
   
   return;
}        
// ------------------------------------------------------------------
void CopyTrjHeader ( Trjh *source,         // source header
                     Trjh *target  )       // target header
{
   sprintf(target->type, "CORD");
   target->nframe     = source->nframe;  
   target->begstp     = source->begstp; 
   target->skpstp     = source->skpstp;
   target->numstp     = source->numstp;
   target->tmstp      = source->tmstp;
   target->nato       = source->nato;
   target->sixtyfour  = source->sixtyfour;
   target->varpbc     = source->varpbc;
   
   return;
}        
// ------------------------------------------------------------------
void WriteDcdHeader( FILE        *dcd,         // dcd file
                     Trjh        *hdr,         // dcd header
                     char        *field )      // field to modify
{
   int offset;
   int dum84 = 84;
   int dum0 = 0;
   int dum164 = 164;
   int dum2 = 2;
   int dum4 = 4;
   
   
   
   if (strstr (field, "all"))
   {
     
     fwrite(       &dum84, 4, 1, dcd);
     fwrite(   &hdr->type, 4, 1, dcd);
     fwrite( &hdr->nframe, 4, 1, dcd);
     fwrite( &hdr->begstp, 4, 1, dcd);
     fwrite( &hdr->skpstp, 4, 1, dcd);
     fwrite( &hdr->numstp, 4, 1, dcd);
     fwrite(        &dum0, 4, 1, dcd);
     fwrite(        &dum0, 4, 1, dcd);
     fwrite(        &dum0, 4, 1, dcd);
     fwrite(        &dum0, 4, 1, dcd);
     fwrite(        &dum0, 4, 1, dcd);
     fwrite(  &hdr->tmstp, 4, 1, dcd);
     fwrite( &hdr->varpbc, 4, 1, dcd);
     fwrite(        &dum0, 4, 1, dcd);
     fwrite(        &dum0, 4, 1, dcd);
     fwrite(        &dum0, 4, 1, dcd);
     fwrite(        &dum0, 4, 1, dcd);
     fwrite(        &dum0, 4, 1, dcd);
     fwrite(        &dum0, 4, 1, dcd);
     fwrite(        &dum0, 4, 1, dcd);
     fwrite(        &dum0, 4, 1, dcd);
     fwrite(        &dum2, 4, 1, dcd);
     fwrite(       &dum84, 4, 1, dcd);
     fwrite(      &dum164, 4, 1, dcd);
     fwrite(        &dum2, 4, 1, dcd);
     fwrite(        "*DCD File Generated by Wordom, the friendly trajectory manipulation tool        ", 80, 1, dcd);
     fwrite(        "*Have fun and may the source be with you...                                     ", 80, 1, dcd);
     fwrite(      &dum164, 4, 1, dcd);
     fwrite(        &dum4, 4, 1, dcd);
     fwrite(   &hdr->nato, 4, 1, dcd);
     fwrite(        &dum4, 4, 1, dcd);
     
   }
   
   if ( strstr(field, "nframe"))
   {
     offset = 8;
     if ( hdr->sixtyfour == 1 )
      offset = 12;
     fseeko( dcd, offset, SEEK_SET);
     fwrite( &hdr->nframe, 4, 1, dcd);
   }

   if(strstr(field, "nato"))
   {
     offset = 268;
     if ( hdr->sixtyfour == 1 )
      offset = 288;
     fseeko( dcd, offset, SEEK_SET);
     fwrite( &hdr->nato, 4, 1, dcd);
   }

   if(strstr(field, "skpstp"))
   {
     offset = 16;
     if ( hdr->sixtyfour == 1 )
      offset = 20;
     fseeko( dcd, offset, SEEK_SET);
     fwrite( &hdr->skpstp, 4, 1, dcd);
   }

   if(strstr(field, "numstp"))
   {
     offset = 20;
     if ( hdr->sixtyfour == 1 )
      offset = 24;
     fseeko( dcd, offset, SEEK_SET);
     fwrite( &hdr->numstp, 4, 1, dcd);
   }

   if(strstr(field, "begstp"))
   {
     offset = 12;
     if ( hdr->sixtyfour == 1 )
      offset = 16;
     fseeko( dcd, offset, SEEK_SET);
     fwrite( &hdr->begstp, 4, 1, dcd);
   }

   if(strstr(field, "tmstp"))
   {
     offset = 44;
     if ( hdr->sixtyfour == 1 )
      offset = 48;
     fseeko( dcd, offset, SEEK_SET);
     fwrite( &hdr->tmstp, 4, 1, dcd);
   }

   if(strstr(field, "varpbc"))
   {
     offset = 48;
     if ( hdr->sixtyfour == 1 )
      offset = 52;
     fseeko( dcd, offset, SEEK_SET);
     fwrite( &hdr->varpbc, 4, 1, dcd);
   }

   return;
}        
// ------------------------------------------------------------------
void WriteTrjHeader( Traj     *traj,
                     char     *field )
{
   char *func="WriteTrjHeader";
   switch ( traj->filetype )
   {
    case 0:                     // write to a dcd
     fprintf( stderr, ".ooo trj type writing not implemented yet\n");
     exit(0);
    case 1:                     // write to a dcd
     WriteDcdHeader ( traj->file, traj->hdr, field );
     break;
    case 2:                     // write to a xtc
      fprintf(stderr, "%s: No separate entity to write a header to an XTC file\n",func);
     break;
    case 3:                     // write a what?
     fprintf( stderr, "this trj type writing not implemented yet\n");
     exit(0);
    default:
     fprintf( stderr, "Error while calling WriteTRJ_Header! check filetype\n");
     exit(78);
     break;
   }
   
   return;
}
// ------------------------------------------------------------------
void ReadDcdCoor ( FILE         *dcd,   // dcd file
                   CoorSet      *crd,   // coordinates
                   Trjh         *hdr,   // dcd header
                   int           frn )  // frame number
{
   int          ii;
   long long    offset;
   char         recl[8];
   char        *buffer;
   int          numovatm;
   int         *movatm;
   float       *movxcoor, *movycoor, *movzcoor;
   
   if(frn == 0 || frn > hdr->nframe)
   {
     fprintf ( stderr," Error: frame number must be between 1 and %d, you asked for # %d\n", hdr->nframe, frn );
     exit(0);
   }
   
   if( hdr->varpbc )
     if( crd->pbc_flag == 0 )
     {
       crd->pbc = calloc( 1, sizeof( Pbc ));
       crd->pbc_flag = 1;
     }
   if (!hdr->fixatm)
   {
    offset  = (long long) (6*4*(1+hdr->sixtyfour) + 23*4  + hdr->comnum*80) ; // size of headers
    offset += (long long) (frn-1) * ( 6*4*(1+hdr->sixtyfour) + 3*4*(hdr->nato) );
    if ( hdr->varpbc)
      offset += (long long) ((frn-1) * (48+8*(1+hdr->sixtyfour)));

    fseeko(dcd, offset, SEEK_SET);

    if ( !hdr->sex )
    {
      if( hdr->varpbc )
      {
        if ( crd->pbc_flag != -1 )
        {
          fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
          fread(&crd->pbc->a_size, 8, 1, dcd);
          fread(&crd->pbc->angle1, 8, 1, dcd);
          fread(&crd->pbc->b_size, 8, 1, dcd);
          fread(&crd->pbc->angle2, 8, 1, dcd);
          fread(&crd->pbc->angle3, 8, 1, dcd);
          fread(&crd->pbc->c_size, 8, 1, dcd);
          fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
        }
        else
        {
          buffer = calloc( (6*8+2*4*(1+hdr->sixtyfour) + 10), sizeof( char) );
          fread( buffer, 1, (6*8+2*4*(1+hdr->sixtyfour)), dcd );
          free(buffer);
        }
      }
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
      fread(crd->xcoor, 4, hdr->nato, dcd);
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
      fread(crd->ycoor, 4, hdr->nato, dcd);
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
      fread(crd->zcoor, 4, hdr->nato, dcd);
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
    }
  
    if ( hdr->sex )
    {
      buffer = calloc ( hdr->nato * 4, sizeof(char) );
    
      if( hdr->varpbc )
      {
        fread( recl  , 1,4*(1+hdr->sixtyfour), dcd);
        fread( buffer, 1, 6*8, dcd);
        if ( crd->pbc_flag != -1 )
        {
          swap ( 8, buffer+0  );
          swap ( 8, buffer+8  );
          swap ( 8, buffer+16 );
          swap ( 8, buffer+24 );
          swap ( 8, buffer+32 );
          swap ( 8, buffer+40 );
          crd->pbc->a_size = (float ) (buffer[0] );
          crd->pbc->angle1 = (float ) (buffer[8] );
          crd->pbc->b_size = (float ) (buffer[16]);
          crd->pbc->angle2 = (float ) (buffer[24]);
          crd->pbc->angle3 = (float ) (buffer[32]);
          crd->pbc->c_size = (float ) (buffer[40]);
        }
        fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
      }
      // buffer the x-coor, swap the byte order and get them
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
      fread( buffer, 1, hdr->nato*4, dcd);
      for ( ii = 0; ii < hdr->nato; ii++)
       {
         swap ( 4, buffer + ( ii*4 ));
         crd->xcoor[ii] = * (float *) (buffer + ii*4);
       }
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);

      // same for the y-coor block
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
      fread( buffer, 1, hdr->nato*4, dcd);
      for ( ii = 0; ii < hdr->nato; ii++)
       {
         swap ( 4, buffer + ( ii*4 ));
         crd->ycoor[ii] = * (float *) (buffer + ii*4);
       }
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);

      // last but not least, the z-coordinates block
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
      fread( buffer, 1, hdr->nato*4, dcd);
      for ( ii = 0; ii < hdr->nato; ii++)
       {
         swap ( 4, buffer + ( ii*4 ));
         crd->zcoor[ii] = * (float *) (buffer + ii*4);
       }
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
      
      free ( buffer );
    }
  }
  
  if ( hdr->fixatm )
   {
    numovatm = hdr->nato - hdr->fixatm;   //number of moving atoms
    
    // alloc the memory for atom numbers and coordinates
    movatm   = calloc ( numovatm, sizeof (int)  );
    movxcoor = calloc ( numovatm, sizeof(float) );
    movycoor = calloc ( numovatm, sizeof(float) );
    movzcoor = calloc ( numovatm, sizeof(float) );
    
    // first, get the atom numbers for the moving atoms
    offset = (long long) 6*4*(1+hdr->sixtyfour) + 23*4 + hdr->comnum*80 + 4*(1+hdr->sixtyfour);
    fseeko(dcd, offset, SEEK_SET);
    if ( !hdr->sex )
      fread ( movatm, 4, numovatm, dcd);
    if ( hdr->sex )
    {
      fread ( buffer, 1, numovatm*4, dcd);
      for ( ii=0; ii<numovatm; ii++)
      {
        swap ( 4, buffer + (ii*4));
        movatm[ii] = *(int *) (buffer + (ii*4));
      }
    }
    
    // now, get the coordinates of ALL atoms (first frame)
    offset = (long long) 6*4*(1+hdr->sixtyfour) + 23*4 + hdr->comnum*80 + 
             2*4*(1+hdr->sixtyfour) + numovatm*4;
    if ( hdr->varpbc)
    {
      offset += (long long) ((frn-1) * (48+8*(1+hdr->sixtyfour)));
    }
    fseeko(dcd, offset, SEEK_SET);
    if ( !hdr->sex )
     {
       if( hdr->varpbc )
       {
         if ( crd->pbc_flag != -1 )
         {
           fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
           fread(&crd->pbc->a_size, 8, 1, dcd);
           fread(&crd->pbc->angle1, 8, 1, dcd);
           fread(&crd->pbc->b_size, 8, 1, dcd);
           fread(&crd->pbc->angle2, 8, 1, dcd);
           fread(&crd->pbc->angle3, 8, 1, dcd);
           fread(&crd->pbc->c_size, 8, 1, dcd);
           fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
         }
         else
         {
           buffer = calloc( (6*8+2*4*(1+hdr->sixtyfour) + 10), sizeof( char) );
           fread( buffer, 1, (6*8+2*4*(1+hdr->sixtyfour)), dcd );
           free(buffer);
         }
       }
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
       fread(crd->xcoor, 4, hdr->nato, dcd);
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
       fread(crd->ycoor, 4, hdr->nato, dcd);
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
       fread(crd->zcoor, 4, hdr->nato, dcd);
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
     }
    
    if ( hdr->sex )
     {
       buffer = calloc ( hdr->nato * 4, sizeof(char) );
       
       if( hdr->varpbc )
       {
         fread( recl  , 1,4*(1+hdr->sixtyfour), dcd);
         fread( buffer, 1, 6*8, dcd);
         if ( crd->pbc_flag != -1 )
         {
           swap ( 8, buffer+0  );
           swap ( 8, buffer+8  );
           swap ( 8, buffer+16 );
           swap ( 8, buffer+24 );
           swap ( 8, buffer+32 );
           swap ( 8, buffer+40 );
           crd->pbc->a_size = (float ) (buffer[0] );
           crd->pbc->angle1 = (float ) (buffer[8] );
           crd->pbc->b_size = (float ) (buffer[16]);
           crd->pbc->angle2 = (float ) (buffer[24]);
           crd->pbc->angle3 = (float ) (buffer[32]);
           crd->pbc->c_size = (float ) (buffer[40]);
         }
         fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
       }
       // buffer the x-coor, swap the byte order and get them
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
       fread( buffer, 1, hdr->nato*4, dcd);
       for ( ii = 0; ii < hdr->nato; ii++)
        {
          swap ( 4, buffer + ( ii*4 ));
          crd->xcoor[ii] = * (float *) (buffer + ii*4);
        }
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);

       // same for the y-coor block
       fread( buffer, 1, hdr->nato*4, dcd);
       for ( ii = 0; ii < hdr->nato; ii++)
        {
          swap ( 4, buffer + ( ii*4 ));
          crd->ycoor[ii] = * (float *) (buffer + ii*4);
        }
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);

       // last but not least, the z-coordinates block
       fread( buffer, 1, hdr->nato*4, dcd);
       for ( ii = 0; ii < hdr->nato; ii++)
        {
          swap ( 4, buffer + ( ii*4 ));
          crd->zcoor[ii] = * (float *) (buffer + ii*4);
        }
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
       
       free ( buffer );

     }

    // if it's the first frame we've been asked, that's it!
    if ( frn == 1 )
       return;

    // otherwise, get the coordinates of the moving atoms and update them
    offset  = (long long) 6*4*(1+hdr->sixtyfour) + 23*4 + hdr->comnum*80 + 
              ( 2*4*(1+hdr->sixtyfour) + 4*numovatm); // size of headers
    offset += (long long) (frn-2) * ( (6*4*(1+hdr->sixtyfour)) + 3*4*(numovatm) ) + 
              ((6*4*(1+hdr->sixtyfour)) + 3*4*hdr->nato);
    if ( hdr->varpbc)
      offset += (long long) ((frn-1) * (48+8*(1+hdr->sixtyfour)));

    fseeko(dcd, offset, SEEK_SET);
     
    if ( !hdr->sex )
     {
       if( hdr->varpbc )
       {
         if ( crd->pbc_flag != -1 )
         {
           fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
           fread(&crd->pbc->a_size, 8, 1, dcd);
           fread(&crd->pbc->angle1, 8, 1, dcd);
           fread(&crd->pbc->b_size, 8, 1, dcd);
           fread(&crd->pbc->angle2, 8, 1, dcd);
           fread(&crd->pbc->angle3, 8, 1, dcd);
           fread(&crd->pbc->c_size, 8, 1, dcd);
           fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
         }
         else
         {
           buffer = calloc( (6*8+2*4*(1+hdr->sixtyfour) + 10), sizeof( char) );
           fread( buffer, 1, (6*8+2*4*(1+hdr->sixtyfour)), dcd );
           free(buffer);
         }
       }
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
      fread(movxcoor, 4, numovatm, dcd);
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
      fread(movycoor, 4, numovatm, dcd);
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
      fread(movzcoor, 4, numovatm, dcd);
      fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
     }
           
    if ( hdr->sex )
     {
       buffer = calloc ( numovatm * 4, sizeof(char) );
    
       if( hdr->varpbc )
       {
         fread( recl  , 1,4*(1+hdr->sixtyfour), dcd);
         fread( buffer, 1, 6*8, dcd);
         if ( crd->pbc_flag != -1 )
         {
           swap ( 8, buffer+0  );
           swap ( 8, buffer+8  );
           swap ( 8, buffer+16 );
           swap ( 8, buffer+24 );
           swap ( 8, buffer+32 );
           swap ( 8, buffer+40 );
           crd->pbc->a_size = (float ) (buffer[0] );
           crd->pbc->angle1 = (float ) (buffer[8] );
           crd->pbc->b_size = (float ) (buffer[16]);
           crd->pbc->angle2 = (float ) (buffer[24]);
           crd->pbc->angle3 = (float ) (buffer[32]);
           crd->pbc->c_size = (float ) (buffer[40]);
         }
         fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
       }
       // buffer the x-coor, swap the byte order and get them
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);
       fread( buffer, 1, numovatm*4, dcd);
       for ( ii = 0; ii < numovatm; ii++)
        {
          swap ( 4, buffer + ( ii*4 ));
          movxcoor[ii] = * (float *) (buffer + ii*4);
        }
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);

       // same for the y-coor block
       fread( buffer, 1, numovatm*4, dcd);
       for ( ii = 0; ii < numovatm; ii++)
        {
          swap ( 4, buffer + ( ii*4 ));
          movycoor[ii] = * (float *) (buffer + ii*4);
        }
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);

       // last but not least, the z-coordinates block
       fread( buffer, 1, numovatm*4, dcd);
       for ( ii = 0; ii < numovatm; ii++)
        {
          swap ( 4, buffer + ( ii*4 ));
          movzcoor[ii] = * (float *) (buffer + ii*4);
        }
       fread(recl     , 1,4*(1+hdr->sixtyfour), dcd);

       free (buffer);
     }
       
    // here is the 'update' part
    for ( ii=0; ii < numovatm; ii++)
     {
      crd->xcoor[movatm[ii]-1] = movxcoor[ii];
      crd->ycoor[movatm[ii]-1] = movycoor[ii];
      crd->zcoor[movatm[ii]-1] = movzcoor[ii];
     }
    free(movatm);
    free(movxcoor);
    free(movycoor);
    free(movzcoor);
   }
   
   return;
}
// ------------------------------------------------------------------
void ReadTrjCoor ( Traj         *traj,
                   CoorSet      *crd,           // coordinates (to be filled)
                   int           frn)           // frame to get
{
   switch ( traj->filetype )
   {
     case 0:
       fprintf( stderr, ".ooo trj type reading not implemented yet\n");
       exit(0);
     case 1:                     // read a dcd
       ReadDcdCoor ( traj->file, crd, traj->hdr, frn );
       break;
     case 2:                     // read a xtc
       ReadXtcCoor ( traj, crd, frn);
       break;
     default:                     // read a what?
       fprintf( stderr, "Error while calling ReadTRJ_Coor! check filetype\n");
       fprintf( stderr, "this trj type reading not implemented yet\n");
       exit(78);
       break;
   }
   
   return;
}
// ------------------------------------------------------------------
void WriteDcdCoor ( FILE        *dcd,         // dcd file
                    CoorSet     *crd,         // coordinates
                    Trjh        *hdr  )       // dcd header
{
   int          temp;
   int          t_varpbc;

   temp = 4*hdr->nato;
   fseeko( dcd, 0, SEEK_END);

   if( !hdr->sixtyfour )
   {
     if( hdr->varpbc && crd->pbc_flag != -1 )
     {
       t_varpbc = 48;
       fwrite( &t_varpbc,    4, 1, dcd);
       fwrite( &crd->pbc->a_size, 8, 1, dcd);
       fwrite( &crd->pbc->angle1, 8, 1, dcd);
       fwrite( &crd->pbc->b_size, 8, 1, dcd);
       fwrite( &crd->pbc->angle2, 8, 1, dcd);
       fwrite( &crd->pbc->angle3, 8, 1, dcd);
       fwrite( &crd->pbc->c_size, 8, 1, dcd);
       fwrite( &t_varpbc,    4, 1, dcd);
     }
     
     fwrite( &temp, 4, 1, dcd);
     fwrite( crd->xcoor, hdr->nato*4, 1, dcd);
     fwrite( &temp, 4, 1, dcd);

     fwrite( &temp, 4, 1, dcd);
     fwrite( crd->ycoor, hdr->nato*4, 1, dcd);
     fwrite( &temp, 4, 1, dcd);

     fwrite( &temp, 4, 1, dcd);
     fwrite( crd->zcoor, hdr->nato*4, 1, dcd);
     fwrite( &temp, 4, 1, dcd);
   }
   else if ( hdr->sixtyfour )
   {
     if( hdr->varpbc && crd->pbc_flag != -1 )
     {
       t_varpbc = 48;
       fwrite( &t_varpbc,    8, 1, dcd);
       fwrite( &crd->pbc->a_size, 8, 1, dcd);
       fwrite( &crd->pbc->angle1, 8, 1, dcd);
       fwrite( &crd->pbc->b_size, 8, 1, dcd);
       fwrite( &crd->pbc->angle2, 8, 1, dcd);
       fwrite( &crd->pbc->angle3, 8, 1, dcd);
       fwrite( &crd->pbc->c_size, 8, 1, dcd);
       fwrite( &t_varpbc,    8, 1, dcd);
     }
     
     fwrite( &temp, 8, 1, dcd);
     fwrite( crd->xcoor, hdr->nato*4, 1, dcd);
     fwrite( &temp, 8, 1, dcd);

     fwrite( &temp, 8, 1, dcd);
     fwrite( crd->ycoor, hdr->nato*4, 1, dcd);
     fwrite( &temp, 8, 1, dcd);

     fwrite( &temp, 8, 1, dcd);
     fwrite( crd->zcoor, hdr->nato*4, 1, dcd);
     fwrite( &temp, 8, 1, dcd);
   }
   
   return;
}
// ------------------------------------------------------------------
void WriteTrjCoor ( CoorSet     *crd,   // coordinates to be written
                    Traj        *traj)
{
   switch ( traj->filetype )
   {
    case 0:                     // write to a dcd
     fprintf( stderr, ".ooo trj type writing not implemented yet\n");
     exit(0);
    case 1:                     // write to a dcd
     WriteDcdCoor ( traj->file, crd, traj->hdr );
     break;
    case 2:                     // write to a xtc
     WriteXtcCoor ( traj->xtc, crd, traj->hdr );
     break;
    default:                     // write a what?
     fprintf( stderr, "Error while calling WriteTrjCoor! check filetype\n");
     fprintf( stderr, "this trj type writing not implemented yet\n");
     exit(78);
     break;
   }
   
   return;
}
// ------------------------------------------------------------------
void ReadFramesList ( char         *filename,     // list file name
                      struct intlist  *str   )    
{
   int          ii = 0;
   char         c;
   FILE        *list;

   list = O_File ( filename, "r" );
   
   while ( !feof(list))
    {
     if (c == '\n')
       ++ii;
     c = fgetc (list);
    }
   rewind(list);

   str->nframes = ii;
   str->nframe = calloc ( (str->nframes), sizeof(int) );
   
   for(ii = 0; ii < str->nframes; ii++)
    fscanf( list, "%d\n", &str->nframe[ii]);
   
   fclose(list);
   
   return;
}
// ------------------------------------------------------------------
void ReadTrjList ( char         *filename,     // list file name
                   TrjList      *str   )       // list structure
{
   int          ii = 0;
   char         c;
   FILE         *list;

   list = O_File ( filename, "r" );
  
   while ( !feof(list))
    {
     if (c == '\n')
      ++ii;
     c = fgetc (list);
    }
   rewind(list);

   str->ntrj = ii;

   str->trjname  = (char **)calloc ( str->ntrj, sizeof(char *)  );

   for ( ii=0; ii < str->ntrj ; ii++)
    str->trjname[ii] = (char *)calloc ( 256, sizeof(char));

   for(ii = 0; ii < str->ntrj; ii++)
    fscanf( list, "%s\n", str->trjname[ii]);

   fclose(list);
   
   return;
}
// ------------------------------------------------------------------
void StringCp( char *source, char *dest, int size)
{
   int n;
   for( n=0 ; n<size ; n++)
    dest[n]=source[n];
   return;
}
// ------------------------------------------------------------------
void CleanRawmol( int nato, Molecule *molecule)
{
   int          ii;
   
   free(molecule->coor.cords);
   free(molecule->rawmol.atomn);
   free(molecule->rawmol.chainId);
   free(molecule->rawmol.altLoc);
   free(molecule->rawmol.iCode);
   free(molecule->rawmol.resn);
   free(molecule->rawmol.occup);
   free(molecule->rawmol.bFac);
   free(molecule->rawmol.presn);
//   free(molecule->rawmol.ipresn);
//   free(molecule->rawmol.cpresn);
   free(molecule->rawmol.weight);
   for( ii=0; ii<molecule->nato; ii++)
     free(molecule->rawmol.atmtype[ii]);
   free(molecule->rawmol.atmtype);
   for( ii=0; ii<molecule->nato; ii++)
     free(molecule->rawmol.atmId[ii]);
   free(molecule->rawmol.atmId);
   for( ii=0; ii<molecule->nato; ii++)
     free(molecule->rawmol.restype[ii]);
   free(molecule->rawmol.restype);
   for( ii=0; ii<molecule->nato; ii++)
     free(molecule->rawmol.segId[ii]);
   free(molecule->rawmol.segId);
   for( ii=0; ii<molecule->nato; ii++)
     free(molecule->rawmol.element[ii]);
   free(molecule->rawmol.element);
   for( ii=0; ii<molecule->nato; ii++)
     free(molecule->rawmol.charge[ii]);
   free(molecule->rawmol.charge);
   
   return;
}
// ------------------------------------------------------------------
void Raw2Struct(Molecule *molecule)
{
  int           ii, jj, kk;
  char         *tmpchainlist=NULL;
  char        **tmpsegid;
  char        **totseglist;
  char          iditem[6]="\0\0\0\0\0\0";
  int          *tmpchainpop=NULL;
  int           tmp_flag, previouschainssegs=0;
  int           iProgResNum, iIndex;
  int           segindex, next_segindex, *resindex=NULL, *atmindex=NULL;
  int           thisindex;
  int           atmcounter;
  int           rescounter;
  int           seg_rescounter;
  int           begres, endres;
  
  /*
  printf( "=== DEBUG: list right _inside_ Raw2Struct (%s) ===\n", molecule->rawmol.name );
  for( ii=1 ; ii<molecule->nato; ii++ )
  {
    printf( "DEBUG: atom %d has presn %d\n", ii+1, molecule->rawmol.presn[ii] );
  }
  printf( "==============  end inside list (%s)  ===========\n", molecule->rawmol.name );
  */

// to fill up segment structures, gather data on the molecule(s)   
// Count the number of segments per molecule (nSeg) - done before 
// Chains in order to have them ready to link
  molecule->nSeg=0;
  totseglist = (char **) calloc( molecule->nato, sizeof( char *));
  for( jj=0; jj<molecule->nato; jj++ )
    totseglist[jj] = NULL;
  //totseglist[0] = (char *) calloc( 5, sizeof(char));

  for ( ii=0; ii<molecule->nato; ii++)
  {
    sprintf( iditem, "%c%4s", molecule->rawmol.chainId[ii], molecule->rawmol.segId[ii] );
    if( isinstrlist( totseglist, molecule->nSeg, iditem ) == -1 )
    {
      totseglist[molecule->nSeg] = (char *)calloc( 5+1, sizeof(char));
      sprintf( totseglist[molecule->nSeg], "%5.5s", iditem );
      molecule->nSeg ++;
    }
  }
  //printf( "DEBUG tot number of segments: %d\n", molecule->nSeg);
  molecule->segment = (Segm *)calloc( molecule->nSeg, sizeof( Segm ));
  
//  Count the number of chains
  molecule->nChains = 0;
  //printf("Pointer %x\n", tmpchainlist);
  tmpchainlist = (char *)calloc( molecule->nato, sizeof(char));
  tmpchainpop = (int *)calloc( molecule->nato, sizeof(int));
  //printf("Pointer %x\n", tmpchainlist);
  for ( ii=0; ii<molecule->nato; ii++)
  {
    tmp_flag = isincharlist( tmpchainlist, molecule->nato, molecule->rawmol.chainId[ii] );
    if ( tmp_flag == -1 )
    {
      tmpchainlist[molecule->nChains] = molecule->rawmol.chainId[ii];
      tmpchainpop[molecule->nChains] += 1;
      molecule->nChains ++;
      //printf( "len(%s): %d (chains:%d) (max:%d)\n", tmpchainlist, strlen(tmpchainlist), molecule->nChains, molecule->nato); fflush(stdout);
    }
    else
      tmpchainpop[tmp_flag] += 1;
  }
  molecule->chains = calloc( molecule->nChains, sizeof( wrd_Chain ));
  molecule->nApC = calloc( molecule->nChains, sizeof( int )) ;
  molecule->nRpC = calloc( molecule->nChains, sizeof( int )) ;
  molecule->nSpC = calloc( molecule->nChains, sizeof( int )) ;
  previouschainssegs = 0;
  for( ii=0; ii<molecule->nChains; ii++ )
  {
    molecule->chains[ii].chainId = tmpchainlist[ii];
    molecule->chains[ii].nApC = tmpchainpop[ii];
    molecule->chains[ii].nRpC=1;
    molecule->chains[ii].nSeg=0;
    tmpsegid = (char **)calloc( molecule->chains[ii].nApC, sizeof(char *));
    for( jj=0; jj<molecule->chains[ii].nApC; jj++ )
      tmpsegid[jj] = NULL;
    //tmpsegid[0] = (char *)calloc( 4, sizeof(char));
    for ( jj=0; jj<molecule->nato; jj++)
      if( molecule->rawmol.chainId[jj] == molecule->chains[ii].chainId )
       //if ( strncmp( molecule->rawmol.segId[jj], molecule->rawmol.segId[jj+1], 4) )
        if ( isinstrlist( tmpsegid, molecule->chains[ii].nSeg, molecule->rawmol.segId[jj]) == -1 )
        {
          tmpsegid[molecule->chains[ii].nSeg] = (char *)calloc( 4+1, sizeof(char));
          sprintf( tmpsegid[molecule->chains[ii].nSeg], "%.4s", molecule->rawmol.segId[jj] );
          molecule->nSpC[ii] ++;
          molecule->chains[ii].nSeg ++;
          // molecule->nSeg ++;
        }
    //molecule->chains[ii].pSegm = (Segm *) calloc( molecule->chains[ii].nSeg, sizeof( Segm ));
    molecule->chains[ii].pSegm = &molecule->segment[previouschainssegs];
    previouschainssegs += molecule->chains[ii].nSeg;
    for( jj=0; jj<molecule->chains[ii].nSeg; jj++ )
    {
      sprintf( molecule->chains[ii].pSegm[jj].segName, "%4.4s", tmpsegid[jj]);
      //printf("DEBUG (chain - segname) : >%c< - >%s<\n", molecule->chains[ii].chainId, molecule->chains[ii].pSegm[jj].segName);
    }
  }

  molecule->nApS    = calloc ( molecule->nSeg, sizeof (  int  ));
  molecule->nRpS    = calloc ( molecule->nSeg, sizeof (  int  ));
  for ( ii=0; ii<molecule->nSeg; ii++)
  {
    molecule->segment[ii].nApS = 0;
    molecule->segment[ii].nRpS = 0;
    molecule->nApS[ii] = 0;
    molecule->nRpS[ii] = 0;
  }
  
//  Count the residues and atoms in each segment; get segment names
  sprintf( iditem, "%c%4s", molecule->rawmol.chainId[0], molecule->rawmol.segId[0] );
  next_segindex = isinstrlist( totseglist, molecule->nSeg, iditem );
  for ( ii=0; ii<molecule->nato-1; ii++)
  {
    segindex = next_segindex;
    sprintf( iditem, "%c%4s", molecule->rawmol.chainId[ii+1], molecule->rawmol.segId[ii+1] );
    next_segindex = isinstrlist( totseglist, molecule->nSeg, iditem );
    StringCp ( molecule->rawmol.segId[ii], molecule->segment[segindex].segName, 4);
    molecule->nApS[segindex]++;
    molecule->segment[segindex].nApS++;
    if( molecule->rawmol.resn[ii] != molecule->rawmol.resn[ii+1] || segindex != next_segindex )
    {
      molecule->nRpS[segindex]++;
      molecule->segment[segindex].nRpS++;
      //printf( "DEBUG: (%d) atom %d - resn %d (%f)\n", segindex, ii, molecule->nRpS[segindex], molecule->coor.xcoor[ii] );
    }
  }
  //pre_segindex = segindex;
  /* take care of last atom... */
  sprintf( iditem, "%c%4s", molecule->rawmol.chainId[ii], molecule->rawmol.segId[ii] );
  segindex = isinstrlist( totseglist, molecule->nSeg, iditem );
  StringCp ( molecule->rawmol.segId[ii], molecule->segment[segindex].segName, 4);
  molecule->nApS[segindex]++;
  molecule->segment[segindex].nApS++;
  /*
  if( molecule->rawmol.resn[ii] != molecule->rawmol.resn[ii-1] || segindex != pre_segindex )
  {
    molecule->nRpS[segindex]++;
    molecule->segment[segindex].nRpS++;
    printf( "DEBUG: atom %d - resn %d (%f)\n", ii, molecule->nRpS[segindex], molecule->coor.xcoor[ii] );
  }
  * always add a residue at the end ... the last residue is over
  */
  molecule->nRpS[segindex]++;
  molecule->segment[segindex].nRpS++;
  //printf( "DEBUG: atom %d - resn %d (%f)\n", ii, molecule->nRpS[segindex], molecule->coor.xcoor[ii] );
  
  // count overall residue number and put it in molecule->nRes
  molecule->nRes=0;
  for( ii=0; ii<molecule->nSeg; ii++)
    molecule->nRes += molecule->nRpS[ii];
  
  molecule->seq1 = (char * )calloc(molecule->nRes, sizeof(char  ));
  molecule->seq3 = (char **)calloc(molecule->nRes, sizeof(char *));
  for( ii=0; ii<molecule->nRes; ii++ )
    molecule->seq3[ii] = (char *)calloc( 3, sizeof(char));
  
  molecule->nApR = calloc (molecule->nRes, sizeof(int));
  atmcounter = 1;
  rescounter = 0;
  begres = 1;
  endres = 0;
  for( ii=1; ii<molecule->nato; ii++ )
  {
    if( strcmp( molecule->rawmol.restype[ii], molecule->rawmol.restype[ii-1] ) ||
        molecule->rawmol.resn[ii] != molecule->rawmol.resn[ii-1] )
      begres = 0;
    if( molecule->rawmol.chainId[ii] != molecule->rawmol.chainId[ii-1] ||
        strcmp( molecule->rawmol.segId[ii], molecule->rawmol.segId[ii-1] ))
    {
      begres = 1;
      endres = 1;
    }
    if( molecule->rawmol.chainId[ii] != molecule->rawmol.chainId[ii-1] ||
        strcmp( molecule->rawmol.segId[ii], molecule->rawmol.segId[ii-1] ) ||
        strcmp( molecule->rawmol.restype[ii], molecule->rawmol.restype[ii-1] ) ||
        molecule->rawmol.resn[ii] != molecule->rawmol.resn[ii-1] )
      {
        molecule->nApR[rescounter] = atmcounter;
        //fprintf( stdout, "atm %4d, %s, res%3d : %3d\n", ii, molecule->rawmol.atmtype[ii], rescounter, molecule->nApR[rescounter]); fflush(stdout);
        Res3ToRes1( molecule->rawmol.restype[ii-1], &molecule->seq1[rescounter] );
        molecule->seq3[rescounter][0] = molecule->rawmol.restype[ii-1][0];
        molecule->seq3[rescounter][1] = molecule->rawmol.restype[ii-1][1];
        molecule->seq3[rescounter][2] = molecule->rawmol.restype[ii-1][2];
        rescounter ++;
        if(rescounter == molecule->nRes)
        {
          fprintf( stderr, "Error counting residues\n" );
          exit(0);
        }
        if( strcmp( molecule->rawmol.segId[ii], molecule->rawmol.segId[ii-1] ) )
        {
          molecule->rawmol.segend[ii-1] = 1;
          molecule->rawmol.segbeg[ii] = 1;
        }
        atmcounter = 1;
      }
    else
      atmcounter ++;
    if( begres )
    {
      molecule->rawmol.segbeg[ii] = 1;
    }
    if( endres )
    {
      if(rescounter > 0)
        for( jj=0; jj<molecule->nApR[rescounter-1]; jj++ )
          molecule->rawmol.segend[ii-jj-1] = 1;
      endres = 0;
    }
  }
  molecule->rawmol.segbeg[0] = 1;
  molecule->nApR[rescounter] = atmcounter;
  for( jj=0; jj<molecule->nApR[rescounter]; jj++ )
    molecule->rawmol.segend[ii-jj-1] = 1;
  //for( ii=0; ii<molecule->nato; ii++ )
  //{
  //  fprintf( stdout, "ii: %4d seg: %4s res: %3d atm: %4d segbeg: %d segend: %d\n", ii, molecule->rawmol.segId[ii], molecule->rawmol.resn[ii], molecule->rawmol.atomn[ii], molecule->rawmol.segbeg[ii], molecule->rawmol.segend[ii] );
  //}

  // now, prepare also pointer array for all residues
  molecule->pRes = (Res **)calloc( molecule->nRes, sizeof(Res *));
  rescounter = 0;
  
  molecule->nApR[molecule->nRes-1] = atmcounter;
  //fprintf( stdout, "nRes : %3d\n", molecule->nRes);
  //for( ii=0; ii<molecule->nRes; ii++)
  //  fprintf( stdout, "res%3d : %3d\n", ii+1, molecule->nApR[ii]);
  //fflush(stdout);
  
//  Allocate space for all residues and set in each the number of atoms to 1
  for ( ii=0; ii< molecule->nSeg; ii++)
  {
    //printf( "DEBUG: nRpS:%d (%d)\n", molecule->segment[ii].nRpS,  molecule->nRpS[ii]);
    molecule->segment[ii].pRes = calloc ( molecule->segment[ii].nRpS, sizeof( Res ) );
    seg_rescounter = 1;
    for ( jj=0; jj<molecule->segment[ii].nRpS; jj++)
    {
      molecule->segment[ii].pRes[jj].nApR = 0;
      molecule->segment[ii].pRes[jj].pAto = NULL;
      molecule->segment[ii].pRes[jj].presn = rescounter+1;
      //molecule->pRes[rescounter] = molecule->segment[ii].pRes[jj];
      //fprintf( stdout, "DEBUG %s %d\n", molecule->segment[ii].segName, molecule->segment[ii].pRes[jj].presn ); fflush(stdout);
      seg_rescounter ++;
      rescounter ++;
    }
  }
  
  rescounter = 0;
//  RE-Count the residues and atoms in each segment; assign them 
  resindex = (int *)calloc( molecule->nSeg, sizeof( int ) );
  for( ii=0; ii< molecule->nSeg; ii++ )
    resindex[ii] = 0;
  sprintf( iditem, "%c%4s", molecule->rawmol.chainId[0], molecule->rawmol.segId[0] );
  next_segindex = isinstrlist( totseglist, molecule->nSeg, iditem );
  for ( ii=0; ii<molecule->nato-1; ii++)
  {
    segindex = next_segindex;
    sprintf( iditem, "%c%4s", molecule->rawmol.chainId[ii+1], molecule->rawmol.segId[ii+1] );
    next_segindex = isinstrlist( totseglist, molecule->nSeg, iditem );
    molecule->segment[segindex].pRes[resindex[segindex]].nApR++;
    if( molecule->rawmol.resn[ii] != molecule->rawmol.resn[ii+1] || segindex != next_segindex )
    {
      StringCp( molecule->rawmol.restype[ii], molecule->segment[segindex].pRes[resindex[segindex]].resType, 3);
      molecule->segment[segindex].pRes[resindex[segindex]].resn  = molecule->rawmol.resn[ii];
      // trying to set up molecule-level residue array
      molecule->pRes[rescounter] = &molecule->segment[segindex].pRes[resindex[segindex]];
      //fprintf( stdout, "DEBUG %s %s %d %d\n", molecule->segment[segindex].segName, molecule->pRes[rescounter]->resType, molecule->pRes[rescounter]->resn, molecule->pRes[rescounter]->presn);
      resindex[segindex]++;
      rescounter ++;
      //printf( "DEBUG: [%d] atom %d - resn %d (%f)\n", segindex, ii, molecule->nRpS[segindex], molecule->coor.xcoor[ii] );
    }
  }
  //pre_segindex = segindex;
  /* take care of last atom... */
  sprintf( iditem, "%c%4s", molecule->rawmol.chainId[ii], molecule->rawmol.segId[ii] );
  segindex = isinstrlist( totseglist, molecule->nSeg, iditem );
  StringCp ( molecule->rawmol.segId[ii], molecule->segment[segindex].segName, 4);
  molecule->segment[segindex].pRes[resindex[segindex]].nApR++;
  /*
  if( molecule->rawmol.resn[ii] != molecule->rawmol.resn[ii-1] || segindex != pre_segindex )
  {
    StringCp( molecule->rawmol.restype[ii], molecule->segment[segindex].pRes[resindex[segindex]].resType, 3);
    molecule->segment[segindex].pRes[resindex[segindex]].resn  = molecule->rawmol.resn[ii];
    resindex[segindex]++;
    //printf( "DEBUG: [%d] atom %d - resn %d (%f)\n", segindex, ii, molecule->nRpS[segindex], molecule->coor.xcoor[ii] );
  }
  * as before, always add a residue at the end ...
  */
  StringCp( molecule->rawmol.restype[ii], molecule->segment[segindex].pRes[resindex[segindex]].resType, 3);
  molecule->segment[segindex].pRes[resindex[segindex]].resn  = molecule->rawmol.resn[ii];
  // trying to set up molecule-level residue array
  molecule->pRes[rescounter] = &molecule->segment[segindex].pRes[resindex[segindex]];
  resindex[segindex]++;
  
  // check if final number of res-per-seg fits 
  for( ii=0; ii< molecule->nSeg; ii++ )
  {
    if ( resindex[ii] != molecule->segment[ii].nRpS )
    {
      fprintf( stderr, "Error in nRpS for segment # %d (id %s)\n", ii+1, molecule->segment[ii].segName );
      fprintf( stderr, "resindex[%d] was %d while nRpS[%d] is %d\n", ii, resindex[ii], ii, molecule->segment[ii].nRpS );
      exit(98);
    }
  }
  
  for ( ii=0; ii< molecule->nSeg; ii++ )
    for ( jj=0; jj< molecule->segment[ii].nRpS; jj++)
      molecule->segment[ii].pRes[jj].pAto = calloc ( molecule->segment[ii].pRes[jj].nApR, sizeof( Atom ) );
  
  for( ii=0; ii< molecule->nSeg; ii++ )
    resindex[ii] = 0;
  atmindex = (int *)calloc( molecule->nSeg, sizeof(int));
  for( ii=0; ii< molecule->nSeg; ii++ )
    atmindex[ii] = 0;
  
  sprintf( iditem, "%c%4s", molecule->rawmol.chainId[0], molecule->rawmol.segId[0] );
  next_segindex = isinstrlist( totseglist, molecule->nSeg, iditem );
  for ( ii=0; ii<molecule->nato-1; ii++)
  {
    segindex = next_segindex;
    sprintf( iditem, "%c%4s", molecule->rawmol.chainId[ii+1], molecule->rawmol.segId[ii+1] );
    next_segindex = isinstrlist( totseglist, molecule->nSeg, iditem );
    //atmindex[segindex]++;
    //printf("DEBUG; ii:%d, segindex:%d, resindex:%d, armindex:%d\n", ii, segindex, resindex[segindex], atmindex[segindex]); fflush(stdout);
    molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].hetatm   = molecule->rawmol.hetatm[ii];
    molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].atomn    = molecule->rawmol.atomn[ii];
    StringCp( molecule->rawmol.atmId[ii],   molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].atmId,   5);
    StringCp( molecule->rawmol.atmtype[ii], molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].atmType, 5);
    molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].altLoc  =  molecule->rawmol.altLoc[ii];
    molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].chainId = molecule->rawmol.chainId[ii];
    molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].iCode   =   molecule->rawmol.iCode[ii];
    molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].xCoor   =   molecule->coor.xcoor[ii];
    molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].yCoor   =   molecule->coor.ycoor[ii];
    molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].zCoor   =   molecule->coor.zcoor[ii];
    molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].occup   =   molecule->rawmol.occup[ii];
    molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].bFac    =    molecule->rawmol.bFac[ii];
    molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].presn   =    molecule->rawmol.presn[ii];
    StringCp( molecule->rawmol.segId[ii],   molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].segId, 4);
    StringCp( molecule->rawmol.element[ii], molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].element, 2);
    StringCp( molecule->rawmol.charge[ii],  molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].charge,  2);
    
    //fprintf( stdout, "DEBUG: %d %d %d\n", segindex, resindex[segindex], atmindex[segindex] );
    //fprintf( stdout, "ATOM  %5d  %-4s%-4s%c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s\n",
                    //molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].atomn, 
                    //molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].atmType, 
                    //molecule->segment[segindex].pRes[resindex[segindex]].resType, 
                    //molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].chainId, 
                    //molecule->segment[segindex].pRes[resindex[segindex]].resn, 
                    //molecule->coor.xcoor[ii], 
                    //molecule->coor.ycoor[ii], 
                    //molecule->coor.zcoor[ii], 
                    //molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].occup, 
                    //molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].bFac , 
                    //molecule->segment[segindex].segName , 
                    //molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].element, 
                    //molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].charge);
    
    atmindex[segindex]++;
    if( molecule->rawmol.resn[ii] != molecule->rawmol.resn[ii+1] || segindex != next_segindex )
    {
      resindex[segindex]++;
      atmindex[segindex] = 0;
    }
  }
  //pre_segindex = segindex;
  /* take care of last atom... */
  ii = molecule->nato-1;
  sprintf( iditem, "%c%4s", molecule->rawmol.chainId[ii], molecule->rawmol.segId[ii] );
  segindex = isinstrlist( totseglist, molecule->nSeg, iditem );
  //atmindex[segindex]++;
  molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].hetatm   = molecule->rawmol.hetatm[ii];
  molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].atomn    = molecule->rawmol.atomn[ii];
  StringCp( molecule->rawmol.atmId[ii],   molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].atmId,   5);
  StringCp( molecule->rawmol.atmtype[ii], molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].atmType, 5);
  molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].altLoc  =  molecule->rawmol.altLoc[ii];
  molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].chainId = molecule->rawmol.chainId[ii];
  molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].iCode   =   molecule->rawmol.iCode[ii];
  molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].xCoor   =   molecule->coor.xcoor[ii];
  molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].yCoor   =   molecule->coor.ycoor[ii];
  molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].zCoor   =   molecule->coor.zcoor[ii];
  molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].occup   =   molecule->rawmol.occup[ii];
  molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].bFac    =    molecule->rawmol.bFac[ii];
  molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].presn   =    molecule->rawmol.presn[ii];
  StringCp( molecule->rawmol.segId[ii],   molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].segId, 4);
  StringCp( molecule->rawmol.element[ii], molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].element, 2);
  StringCp( molecule->rawmol.charge[ii],  molecule->segment[segindex].pRes[resindex[segindex]].pAto[atmindex[segindex]].charge,  2);
  
  //ShowMolInfo( molecule );
  // end experimental code 
/* Pdb does not have the ipresn (internal (per segment) progressive residue number) field: let's fill it */
  //molecule->rawmol.ipresn = calloc( molecule->nato, sizeof( int )); // redundant? necessary?
  iProgResNum=0;
  iIndex=-1;
  for(ii=0; ii<molecule->nSeg; ii++ )
  {
    for(jj=0; jj<molecule->segment[ii].nRpS; jj++)
    {
      iProgResNum++;
      //printf ( "DEBUG: %s %s %d %d : %d (%d)\n", molecule->segment[ii].segName, molecule->segment[ii].pRes[jj].resType, molecule->segment[ii].pRes[jj].resn, molecule->segment[ii].pRes[jj].nApR, iProgResNum, molecule->segment[ii].pRes[jj].pAto[0].atomn-1 ); fflush(stdout);
      for(kk=0; kk<molecule->segment[ii].pRes[jj].nApR; kk++)
      {
        iIndex++;
        thisindex = molecule->segment[ii].pRes[jj].pAto[kk].atomn-1;
        //printf( "DEBUG: accessing atom # %d ( presn %d)\n", thisindex, molecule->rawmol.presn[thisindex] );
        //molecule->rawmol.ipresn[iIndex]=iProgResNum;
        //molecule->rawmol.presn[molecule->segment[ii].pRes[jj].pAto[kk].atomn-1] = iProgResNum;
        //fprintf( stdout, "DEBUG - thisindex: %d (nato:%d)\n", thisindex, molecule->nato );
        molecule->rawmol.presn[thisindex] = iProgResNum;
        molecule->segment[ii].pRes[jj].pAto[kk].presn = iProgResNum;
      }
    }
  }
  
  return;
}
// ------------------------------------------------------------------
void old_Raw2Struct(Molecule *molecule)
{
   int        ii, jj, kk;
   int        newpiece, prev, *prevres, *rescount;
   char     **tmpsegnames;

//  Count the number of segments
   molecule->nSeg=1;
   tmpsegnames = (char **)calloc( 1, sizeof(char *));
   for( ii=0; ii<molecule->nato; ii++)
     tmpsegnames[ii] = (char *)calloc( 4, sizeof(char));
   sprintf( tmpsegnames[0], "%s", molecule->rawmol.segId[0]);
   for ( ii=1; ii<molecule->nato-1; ii++)
   {
     newpiece = 1;
     for( jj=0; jj<molecule->nSeg; jj++)
       if( !strcmp( molecule->rawmol.segId[ii], tmpsegnames[jj] ) )
       {
         newpiece = 0;
         break;
       }
     if( newpiece )
     {
       sprintf( tmpsegnames[molecule->nSeg], "%s", molecule->rawmol.segId[ii]);
       molecule->nSeg ++;
     }
   }
   printf("# of segments: %d\n", molecule->nSeg); 
   for( ii=0; ii<molecule->nSeg; ii++)
     printf("%s\n", tmpsegnames[ii]);

// Allocate space for the segments and set in each the number of atoms/residues to 0
   molecule->segment = NULL;
   molecule->segment = calloc ( molecule->nSeg, sizeof ( Segm ));
   if ( molecule->segment == NULL )
    printf(" ?!? Did not actually calloc!\n");
    
   molecule->nApS    = calloc ( molecule->nSeg, sizeof (  int  ));
   molecule->nRpS    = calloc ( molecule->nSeg, sizeof (  int  ));
   for ( ii=0; ii<molecule->nSeg; ii++)
    {
     molecule->segment[ii].nApS = 0;
     molecule->segment[ii].nRpS = 0;
     molecule->nApS[ii] = 0;
     molecule->nRpS[ii] = 0;
     sprintf( molecule->segment[ii].segName, "%s", tmpsegnames[ii] );
    }

//  Count the residues and atoms in each segment; get segment names
   prev = -1;
   for( ii=0; ii<molecule->nSeg; ii++ )
   {
     for( jj=0; jj<molecule->nato; jj++)
     {
       if(!strncmp(molecule->segment[ii].segName, molecule->rawmol.segId[jj], 4))
       {
         molecule->nApS[ii]++;
         molecule->segment[ii].nApS++;
         if( molecule->rawmol.resn[ii] != prev )
         {
           molecule->nRpS[ii]++;
           molecule->segment[ii].nRpS++;
           prev = molecule->rawmol.resn[ii];
         }
       }
     }
   }
   for( ii=0; ii<molecule->nSeg; ii++)
     printf("seg# %d (%s); nres:%d; natom: %d\n", ii+1, molecule->segment[ii].segName, molecule->segment[ii].nRpS, molecule->segment[ii].nApS);
   exit(0);

   // count overall residue number and put it in molecule->nRes
   molecule->nRes=0;
   for( ii=0; ii<molecule->nSeg; ii++)
     molecule->nRes += molecule->nRpS[ii];
   // WARNING: allocated, not filled 
   molecule->nApR = calloc (molecule->nRes, sizeof(int));
   
//  Allocate space for all residues and set in each the number of atoms to 1
   for ( ii=0; ii< molecule->nSeg; ii++)
   {
     molecule->segment[ii].pRes = calloc ( molecule->segment[ii].nRpS, sizeof( Res ) );
     for ( jj=0; jj<molecule->segment[ii].nRpS; jj++)
       molecule->segment[ii].pRes[jj].nApR = 0;
   }
   
//  Count atoms in each residue and set res.res_Type, res.resn and res.resId
   prevres = calloc( molecule->nSeg, sizeof(int));
   rescount = calloc( molecule->nSeg, sizeof(int));
   for ( ii=0; ii<molecule->nato-1; ii++)
   {
     for( jj=0; jj<molecule->nSeg; jj++)
     {
       if(!strcmp(molecule->rawmol.segId[ii], molecule->segment[jj].segName))
       {
         if( molecule->rawmol.resn[ii] == prevres[jj])
           molecule->segment[jj].pRes[rescount[jj]].nApR++;
         else
         {
           rescount[jj] ++;
           molecule->segment[jj].pRes[rescount[jj]].nApR++;
           prevres[jj] = molecule->rawmol.resn[ii];
         }
       }
     }
   }
   jj = 0;
   kk = 0;
   for ( ii=0; ii<molecule->nato-1; ii++)
   {
     if ( molecule->rawmol.resn[ii] == molecule->rawmol.resn[ii+1] && !strncmp(molecule->rawmol.segId[ii], molecule->rawmol.segId[ii+1], 4) )
      molecule->segment[jj].pRes[kk].nApR++;
     if ( !strncmp(molecule->rawmol.segId[ii], molecule->rawmol.segId[ii+1], 4) && (molecule->rawmol.resn[ii] != molecule->rawmol.resn[ii+1]) )
      {
       StringCp( molecule->rawmol.restype[ii], molecule->segment[jj].pRes[kk].resType, 3);
       molecule->segment[jj].pRes[kk].resn  = molecule->rawmol.resn[ii];
       kk++ ;
      }
     if (  strncmp(molecule->rawmol.segId[ii], molecule->rawmol.segId[ii+1], 4) )
      {
       if ( kk+1 != molecule->segment[jj].nRpS )
        {
         fprintf( stderr, "Error in nRpS for segment # %d (%s - %s)- aborting\n", jj+1, molecule->rawmol.segId[ii], molecule->rawmol.segId[ii+1] );
         fprintf( stderr, "kk+1 was %d while nRpS is %d\n", kk+1, molecule->segment[jj].nRpS );
         exit(98);
        }
       StringCp( molecule->rawmol.restype[ii], molecule->segment[jj].pRes[kk].resType, 3);
       molecule->segment[jj].pRes[kk].resn  = molecule->rawmol.resn[ii];
       jj++;
       kk = 0;
      }
   }
   
   return;
}
// ------------------------------------------------------------------
void ReadPdb ( char *name, Molecule *molecule )
{
  FILE        *pdb;
  int          ii;
  char         buffer[80], varbuffer[13];
  int          skewed = 0;
  int         *tmpchainpop=NULL;
  int          tmp_presn = 0;
  
  //fprintf( stdout, "DEBUG: called ReadPdb on %s\n", name ); fflush(stdout);
  pdb = O_File( name, "r");
  molecule->nModels = 0;
  molecule->altLocB = 0;
  molecule->altLocC = 0;
  molecule->rawmol.name = calloc( strlen(name)+1, sizeof(char));
  strcpy( molecule->rawmol.name, name );

// first, get how many atoms are in the structure file
  molecule->nato = 0;
  molecule->rawmol.b_cryst=0;
  memset( buffer, '\0', sizeof(buffer));
  while ( !feof(pdb))
  {
    if ( !strncmp(buffer, "CRYST1", 6) ) 
    {
      memset( varbuffer, '\0', sizeof(varbuffer));
      copyvar( buffer, 6, 15, varbuffer, 12);
      sscanf( varbuffer, "%f", &molecule->rawmol.a);
      memset( varbuffer, '\0', sizeof(varbuffer));
      copyvar( buffer, 15, 24, varbuffer, 12);
      sscanf( varbuffer, "%f", &molecule->rawmol.b);
      memset( varbuffer, '\0', sizeof(varbuffer));
      copyvar( buffer, 24, 33, varbuffer, 12);
      sscanf( varbuffer, "%f", &molecule->rawmol.c);
      memset( varbuffer, '\0', sizeof(varbuffer));
      copyvar( buffer, 33, 40, varbuffer, 12);
      sscanf( varbuffer, "%f", &molecule->rawmol.alpha);
      memset( varbuffer, '\0', sizeof(varbuffer));
      copyvar( buffer, 40, 47, varbuffer, 12);
      sscanf( varbuffer, "%f", &molecule->rawmol.beta);
      memset( varbuffer, '\0', sizeof(varbuffer));
      copyvar( buffer, 47, 54, varbuffer, 12);
      sscanf( varbuffer, "%f", &molecule->rawmol.gamma);
      
      copyvar( buffer, 55, 66, molecule->rawmol.s_group, 11);
      memset( varbuffer, '\0', sizeof(varbuffer));
      copyvar( buffer, 66, 70, varbuffer, 12);
      sscanf( varbuffer, "%d", &molecule->rawmol.z);

      molecule->rawmol.b_cryst=1;
    }
    if ( !strncmp(buffer, "ATOM  ", 6) || !strncmp(buffer, "HETATM", 6) )
    {
      if( buffer[16]==' ' || buffer[16]=='A' )
        ++molecule->nato;
      else if( buffer[16]=='B' )
      {
        //fprintf( stderr, "DEBUG: alternate location (B) data found:\n%s\n", buffer);
        molecule->altLocB = 1;
      }
      else if( buffer[16]=='C' )
      {
        //fprintf( stderr, "DEBUG: alternate location (C) data found\n");
        molecule->altLocC = 1;
      }
      else
        fprintf( stderr, "Warning: found alternate location >C (ignored)\n");
    }
    if( buffer[0] == '*' )
    {
      fprintf( stderr, "Found field starting with * : maybe %s is actually a crd file?\n", name);
      fflush(stdout);
    }
    if( !strncmp( buffer, "ENDMDL", 6 ) || !strncmp( buffer, "ENDMDL", 6 ) )
    {
      molecule->nModels++;
    }
    if(fgets(buffer, 80, pdb)==NULL && ( !feof(pdb) || ferror(pdb) ))
    {
       fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
  }
  rewind(pdb);
  if( molecule->nModels == 0 )
    molecule->nModels = 1;
  else if( molecule->nModels > 1 )
  {
    molecule->cords = (CoorSet *)calloc( molecule->nModels, sizeof( CoorSet ));
  }
  
  molecule->nato /= molecule->nModels;
  molecule->rawmol.nato = molecule->nato;
  molecule->coor.cords = walloc ( 3*molecule->nato, sizeof(float));
  molecule->coor.xcoor = &molecule->coor.cords[0];
  molecule->coor.ycoor = &molecule->coor.cords[molecule->nato];
  molecule->coor.zcoor = &molecule->coor.cords[2*molecule->nato];
  molecule->coor.coors    =  calloc( 3, sizeof(float *));
  molecule->coor.coors[0] =  molecule->coor.cords;
  molecule->coor.coors[1] = &molecule->coor.cords[molecule->nato];
  molecule->coor.coors[2] = &molecule->coor.cords[2*molecule->nato];
  
  if( molecule->altLocB )
    CalloCoor( &molecule->coor2, molecule->nato );
  if( molecule->altLocC )
    CalloCoor( &molecule->coor3, molecule->nato );
  
  molecule->rawmol.hetatm  = calloc ( molecule->nato, sizeof(short) );
  molecule->rawmol.atomn   = calloc ( molecule->nato, sizeof( int ) );
  molecule->rawmol.chainId = calloc ( molecule->nato, sizeof( char) );
  molecule->rawmol.altLoc  = calloc ( molecule->nato, sizeof( char) );
  molecule->rawmol.iCode   = calloc ( molecule->nato, sizeof( char) );
//  molecule->rawmol.ipresn  = calloc ( molecule->nato, sizeof( int ) );
  molecule->rawmol.presn   = calloc ( molecule->nato, sizeof( int ) );
  molecule->rawmol.resn    = calloc ( molecule->nato, sizeof( int ) );
  molecule->rawmol.occup   = calloc ( molecule->nato, sizeof(float) );
  molecule->rawmol.bFac    = calloc ( molecule->nato, sizeof(float) );
  molecule->rawmol.segend  = calloc ( molecule->nato, sizeof(short) );
  molecule->rawmol.segbeg  = calloc ( molecule->nato, sizeof(short) );
  
  // this two just for crd-compatibility
//  molecule->rawmol.cpresn  = calloc ( molecule->nato, sizeof( int ) );
  molecule->rawmol.weight  = calloc ( molecule->nato, sizeof(float) );

  molecule->rawmol.atmtype    = calloc ( molecule->nato, sizeof(char *));
  for ( ii=0; ii < molecule->nato ; ii++)
   molecule->rawmol.atmtype[ii] = calloc( 4+1, sizeof(char));

  molecule->rawmol.atmId      = calloc ( molecule->nato, sizeof(char *));
  for ( ii=0; ii < molecule->nato ; ii++)
   molecule->rawmol.atmId[ii]   = calloc( 5+1, sizeof(char));

  molecule->rawmol.restype    = calloc ( molecule->nato, sizeof(char *));
  for ( ii=0; ii < molecule->nato ; ii++)
   molecule->rawmol.restype[ii] = calloc( 3+1, sizeof(char));

  molecule->rawmol.segId      = calloc ( molecule->nato, sizeof(char *));
  for ( ii=0; ii < molecule->nato ; ii++)
   molecule->rawmol.segId[ii]   = calloc( 4+1, sizeof(char));

  molecule->rawmol.element    = calloc ( molecule->nato, sizeof(char *));
  for ( ii=0; ii < molecule->nato ; ii++)
   molecule->rawmol.element[ii] = calloc( 2+1, sizeof(char));

  molecule->rawmol.charge     = calloc ( molecule->nato, sizeof(char *));
  for ( ii=0; ii < molecule->nato ; ii++)
   molecule->rawmol.charge[ii] = calloc( 2+1, sizeof(char));
  
// fill up the pbds and coor structures (simple serial reading) - only time data is acquired from source (pdb file)
  for ( ii=0; ii<molecule->nato; )
  {
    if(fgets(buffer, 80, pdb)==NULL && ( !feof(pdb) || ferror(pdb) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
    if ( !strncmp(buffer, "ATOM  ", 6) || !strncmp(buffer, "HETATM", 6) )
    {
      if( buffer[16]==' ' || buffer[16]=='A' )
      {
        if( !strncmp(buffer, "ATOM  ", 6) )
          molecule->rawmol.hetatm[ii] = 0;
        else if ( !strncmp(buffer, "HETATM", 6) )
          molecule->rawmol.hetatm[ii] = 1;
        memset( varbuffer, '\0', sizeof(varbuffer));
        copyvar( buffer, 6, 11, varbuffer, 10);
        sscanf( varbuffer, "%5s", molecule->rawmol.atmId[ii]);
        copyvar( buffer, 12, 16, varbuffer, 10);
        sscanf( varbuffer, "%4s", molecule->rawmol.atmtype[ii]);
        molecule->rawmol.altLoc[ii] = buffer[16];
        copyvar( buffer, 17, 20, varbuffer, 10);
        sscanf( varbuffer, "%3s", molecule->rawmol.restype[ii]);
        molecule->rawmol.chainId[ii] = buffer[21];
        copyvar( buffer, 22, 26, varbuffer, 10);
        sscanf( varbuffer, "%4d", &molecule->rawmol.resn[ii]);
        molecule->rawmol.iCode[ii] = buffer[26];
        copyvar( buffer, 30, 38, varbuffer, 10);
        sscanf( varbuffer, "%8f", &molecule->coor.xcoor[ii]);
        copyvar( buffer, 38, 46, varbuffer, 10);
        sscanf( varbuffer, "%8f", &molecule->coor.ycoor[ii]);
        copyvar( buffer, 46, 54, varbuffer, 10);
        sscanf( varbuffer, "%8f", &molecule->coor.zcoor[ii]);
        copyvar( buffer, 54, 60, varbuffer, 10);
        sscanf( varbuffer, "%6f", &molecule->rawmol.occup[ii]);
        copyvar( buffer, 60, 66, varbuffer, 10);
        sscanf( varbuffer, "%6f", &molecule->rawmol.bFac[ii]);
//       copyvar( buffer, 72, 76, molecule->rawmol.segId[ii]);     // works but gives trouble if field does not exists
        copyvar( buffer, 72, 76, varbuffer, 10);
        sscanf( varbuffer, "%4s", molecule->rawmol.segId[ii]);
//       copyvar( buffer, 76, 78, molecule->rawmol.element[ii]);   // works but gives trouble if field does not exists
        copyvar( buffer, 76, 78, varbuffer, 10);
        sscanf( varbuffer, "%2s", molecule->rawmol.element[ii]);
//       copyvar( buffer, 78, 80, molecule->rawmol.charge[ii]);    // works but gives trouble if field does not exists
        copyvar( buffer, 78, 80, varbuffer, 10);
        sscanf( varbuffer, "%2s", molecule->rawmol.charge[ii]);
      
        molecule->rawmol.atomn[ii] = ii+1;
        if ( strncmp(molecule->rawmol.atmId[ii], "*****", 5) && molecule->rawmol.atomn[ii] != atoi(molecule->rawmol.atmId[ii]) )
        {
          //fprintf( stderr, "Warning: skewed atom numbering (#%d - %d (%s))\n", molecule->rawmol.atomn[ii], atoi(molecule->rawmol.atmId[ii]), molecule->rawmol.atmId[ii] );
          skewed = 1;
        }
        if( molecule->rawmol.chainId[ii] == '\0' || molecule->rawmol.chainId[ii] == ' ' )
          molecule->rawmol.chainId[ii] = 'A';
        // taken from crdread: adapt!
        if( molecule->rawmol.segId[ii][0] == '\0' )
        {
          if( molecule->rawmol.chainId[ii] != '\0' )
            molecule->rawmol.segId[ii][0] = molecule->rawmol.chainId[ii];
          else
            molecule->rawmol.segId[ii][0] = 'A';
        }
        molecule->rawmol.weight[ii] = molecule->rawmol.bFac[ii];
        // end copied section
        
        if( molecule->altLocB )
        {
          molecule->coor2.xcoor[ii] = molecule->coor.xcoor[ii];
          molecule->coor2.ycoor[ii] = molecule->coor.ycoor[ii];
          molecule->coor2.zcoor[ii] = molecule->coor.zcoor[ii];
        }
        if( molecule->altLocC )
        {
          molecule->coor3.xcoor[ii] = molecule->coor.xcoor[ii];
          molecule->coor3.ycoor[ii] = molecule->coor.ycoor[ii];
          molecule->coor3.zcoor[ii] = molecule->coor.zcoor[ii];
        }
        
        ii++;
      }
      else if( buffer[16]=='B' )
      {
        copyvar( buffer, 30, 38, varbuffer, 10);
        sscanf( varbuffer, "%8f", &molecule->coor2.xcoor[ii-1]);
        copyvar( buffer, 38, 46, varbuffer, 10);
        sscanf( varbuffer, "%8f", &molecule->coor2.ycoor[ii-1]);
        copyvar( buffer, 46, 54, varbuffer, 10);
        sscanf( varbuffer, "%8f", &molecule->coor2.zcoor[ii-1]);
        
      }
      else if( buffer[16]=='C' )
      {
        copyvar( buffer, 30, 38, varbuffer, 10);
        sscanf( varbuffer, "%8f", &molecule->coor3.xcoor[ii-1]);
        copyvar( buffer, 38, 46, varbuffer, 10);
        sscanf( varbuffer, "%8f", &molecule->coor3.ycoor[ii-1]);
        copyvar( buffer, 46, 54, varbuffer, 10);
        sscanf( varbuffer, "%8f", &molecule->coor3.zcoor[ii-1]);
        
      }
    }
  }
  if( skewed )
    fprintf( stderr, "Warning: skewed atom numbering\n");
  
  tmp_presn = 1;
  molecule->rawmol.presn[1] = 1;
  for( ii=1 ; ii<molecule->nato; ii++ )
  {
    if( molecule->rawmol.resn[ii] != molecule->rawmol.resn[ii-1] )
      tmp_presn ++;
    else if( strncmp( molecule->rawmol.segId[ii],molecule->rawmol.segId[ii], 4 ) )
      tmp_presn ++;
    else if( molecule->rawmol.chainId[ii] != molecule->rawmol.chainId[ii] )
      tmp_presn ++;
    molecule->rawmol.presn[ii] = tmp_presn;
  }
  /*
  printf( "=== DEBUG: list right before Raw2Struct (%s) ===\n", molecule->rawmol.name );
  for( ii=1 ; ii<molecule->nato; ii++ )
  {
    printf( "DEBUG: atom %d has presn %d\n", ii+1, molecule->rawmol.presn[ii] );
  }
  printf( "==============end list before (%s)===============\n", molecule->rawmol.name );
  */
  
  Raw2Struct( molecule );
  
  molecule->coor.pbc_flag = 0;
  molecule->coor.pbc = NULL;
  if( molecule->rawmol.b_cryst == 1 )
  {
    molecule->a =       molecule->rawmol.a;
    molecule->b =       molecule->rawmol.b;
    molecule->c =       molecule->rawmol.c;
    molecule->alpha =   molecule->rawmol.alpha;
    molecule->beta =    molecule->rawmol.beta;
    molecule->gamma =   molecule->rawmol.gamma;
    strncpy( molecule->s_group, molecule->rawmol.s_group,11);
    molecule->z =       molecule->rawmol.z;
    
    molecule->coor.pbc_flag    = 1;
    molecule->coor.pbc = calloc( 1, sizeof( Pbc ));
    molecule->coor.pbc->a_size = molecule->rawmol.a;
    molecule->coor.pbc->b_size = molecule->rawmol.b;
    molecule->coor.pbc->c_size = molecule->rawmol.c;
    molecule->coor.pbc->angle1 = cos(molecule->rawmol.alpha * (PIE/180));
    molecule->coor.pbc->angle2 = cos(molecule->rawmol.beta  * (PIE/180));
    molecule->coor.pbc->angle3 = cos(molecule->rawmol.gamma * (PIE/180));
    strncpy( molecule->s_group, molecule->rawmol.s_group, 11);
    molecule->coor.pbc->zval       = molecule->rawmol.z;
  }
  
  molecule->origin = 1;
  molecule->filled = calloc( 4, sizeof(char)); // overalloc to avoid overflow
  sprintf( molecule->filled, "YES");
  
  fclose(pdb);
  
  if( tmpchainpop != NULL )
    free(tmpchainpop);
  
  
  return;
}
// ------------------------------------------------------------------
void WritePdb ( char *name, Molecule *molecule )
{
   int     ii, jj, kk, ll=0;
   FILE   *out;
   int     atomn = 0;
   
   out = fopen( name, "w");
   
   if( molecule->rawmol.b_cryst == 1 )
     fprintf( out, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %11s%4d\n", molecule->rawmol.a, molecule->rawmol.b, molecule->rawmol.c, molecule->rawmol.alpha, molecule->rawmol.beta, molecule->rawmol.gamma, molecule->rawmol.s_group, molecule->rawmol.z );
   
   ll=0;
   for ( ii = 0; ii < molecule->nSeg; ii++ )
    for ( jj = 0; jj < molecule->segment[ii].nRpS; jj++ )
     for ( kk = 0; kk < molecule->segment[ii].pRes[jj].nApR; kk++ )
     {
      if( molecule->segment[ii].pRes[jj].pAto[kk].atomn < 100000 )
        //atomn = molecule->segment[ii].pRes[jj].pAto[kk].atomn;
        atomn = ll+1;
      else
        //atomn = molecule->segment[ii].pRes[jj].pAto[kk].atomn % 100000;
        atomn = (ll+1) % 100000;
      if( molecule->segment[ii].pRes[jj].pAto[kk].hetatm == 0 )
        fprintf( out, "ATOM  ");
      else
        fprintf( out, "HETATM" );
      if ( (strlen (molecule->segment[ii].pRes[jj].pAto[kk].atmType)) < 4 ) 
       fprintf ( out, "%5d  %-4s%-4s%c%4d    %8.*f%8.*f%8.*f%6.2f%6.2f      %-4s%2s%2s\n",
                       atomn, 
                       molecule->segment[ii].pRes[jj].pAto[kk].atmType, 
                       molecule->segment[ii].pRes[jj].resType, 
                       molecule->segment[ii].pRes[jj].pAto[kk].chainId, 
                       molecule->segment[ii].pRes[jj].resn, 
                       //molecule->coor.xcoor[ll], 
                       //molecule->coor.ycoor[ll], 
                       //molecule->coor.zcoor[ll], 
                        /* The x,y,z coordinated are written normally with 3 figures after the comma. If less then -1000, only two figures are printed. */
                       molecule->segment[ii].pRes[jj].pAto[kk].xCoor <-1000 ? 2 : 3, molecule->segment[ii].pRes[jj].pAto[kk].xCoor,
                       molecule->segment[ii].pRes[jj].pAto[kk].yCoor <-1000 ? 2 : 3, molecule->segment[ii].pRes[jj].pAto[kk].yCoor,
                       molecule->segment[ii].pRes[jj].pAto[kk].zCoor <-1000 ? 2 : 3, molecule->segment[ii].pRes[jj].pAto[kk].zCoor,
                       molecule->segment[ii].pRes[jj].pAto[kk].occup, 
                       molecule->segment[ii].pRes[jj].pAto[kk].bFac , 
                       molecule->segment[ii].segName , 
                       molecule->segment[ii].pRes[jj].pAto[kk].element, 
                       molecule->segment[ii].pRes[jj].pAto[kk].charge);
      else
       fprintf ( out, "%5d %-4s %-4s%c%4d    %8.*f%8.*f%8.*f%6.2f%6.2f      %-4s%2s%2s\n", 
                       atomn,
                       molecule->segment[ii].pRes[jj].pAto[kk].atmType, 
                       molecule->segment[ii].pRes[jj].resType, 
                       molecule->segment[ii].pRes[jj].pAto[kk].chainId, 
                       molecule->segment[ii].pRes[jj].resn, 
                       //molecule->coor.xcoor[ll], 
                       //molecule->coor.ycoor[ll], 
                       //molecule->coor.zcoor[ll], 
                        /* The x,y,z coordinated are written normally with 3 figures after the comma. If less then -1000, only two figures are printed. */
                       molecule->segment[ii].pRes[jj].pAto[kk].xCoor <-1000 ? 2 : 3, molecule->segment[ii].pRes[jj].pAto[kk].xCoor,
                       molecule->segment[ii].pRes[jj].pAto[kk].yCoor <-1000 ? 2 : 3, molecule->segment[ii].pRes[jj].pAto[kk].yCoor,
                       molecule->segment[ii].pRes[jj].pAto[kk].zCoor <-1000 ? 2 : 3, molecule->segment[ii].pRes[jj].pAto[kk].zCoor,
                       molecule->segment[ii].pRes[jj].pAto[kk].occup, 
                       molecule->segment[ii].pRes[jj].pAto[kk].bFac , 
                       molecule->segment[ii].segName , 
                       molecule->segment[ii].pRes[jj].pAto[kk].element , 
                       molecule->segment[ii].pRes[jj].pAto[kk].charge);
      ll++;
     }
  
   fprintf(out, "END\n");
   fclose (out);
   return;
}
// ------------------------------------------------------------------
// ------------------------------------------------------------------
void WritePdb_unstr ( char *name, Molecule *molecule )
{
   int     ii;
   FILE   *out;
   int     atomn;
   
   out = fopen( name, "w");
   
   for ( ii = 0; ii < molecule->nato; ii++ )
   {
    if( molecule->rawmol.atomn[ii] < 100000 )
      atomn = molecule->rawmol.atomn[ii] ;
    else
      atomn = molecule->rawmol.atomn[ii]%100000;
    if ( molecule->rawmol.hetatm[ii] == 0 )
      fprintf( out, "ATOM  ");
    else
      fprintf( out, "HETATM" );
    if ( (strlen (molecule->rawmol.atmtype[ii])) < 4 ) 
     fprintf ( out, "%5d  %-4s%-3s %c%4d    %8.*f%8.*f%8.*f%6.2f%6.2f      %-4s%2s%2s\n",
                        atomn,
                        molecule->rawmol.atmtype[ii], 
                        molecule->rawmol.restype[ii], 
                        molecule->rawmol.chainId[ii],
                        molecule->rawmol.resn[ii],
                        /* The x,y,z coordinated are written normally with 3 figures after the comma. If less then -1000, only two figures are printed. */
                        molecule->coor.xcoor[ii] <-1000 ? 2 : 3, molecule->coor.xcoor[ii],
                        molecule->coor.ycoor[ii] <-1000 ? 2 : 3, molecule->coor.ycoor[ii],
                        molecule->coor.zcoor[ii] <-1000 ? 2 : 3, molecule->coor.zcoor[ii],
                        molecule->rawmol.occup[ii],
                        molecule->rawmol.bFac[ii],
                        molecule->rawmol.segId[ii],
                        molecule->rawmol.element[ii], 
                        molecule->rawmol.charge[ii]);
     else 
     //fprintf ( out, "ATOM  %5d  %-4s%-3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s\n",
     fprintf ( out, "%5d %-4s %-4s%c%4d    %8.*f%8.*f%8.*f%6.2f%6.2f      %-4s%2s%2s\n", 
                        atomn,
                        molecule->rawmol.atmtype[ii], 
                        molecule->rawmol.restype[ii], 
                        molecule->rawmol.chainId[ii],
                        molecule->rawmol.resn[ii],
                        /* The x,y,z coordinated are written normally with 3 figures after the comma. If less then -1000, only two figures are printed. */
                        molecule->coor.xcoor[ii] <-1000 ? 2 : 3, molecule->coor.xcoor[ii],
                        molecule->coor.ycoor[ii] <-1000 ? 2 : 3, molecule->coor.ycoor[ii],
                        molecule->coor.zcoor[ii] <-1000 ? 2 : 3, molecule->coor.zcoor[ii],
                        molecule->rawmol.occup[ii],
                        molecule->rawmol.bFac[ii],
                        molecule->rawmol.segId[ii],
                        molecule->rawmol.element[ii], 
                        molecule->rawmol.charge[ii]);
    
      }
  
   fprintf(out, "END\n");
   fclose (out);
   return;
}
// ------------------------------------------------------------------
// ------------------------------------------------------------------
void ReadCrd ( char *name, Molecule *molecule )
{
   FILE   *crd;
   int     ii;
   char    buffer[150], varbuffer[13];
   int     nato, nlines, crdnatofound = 0;
   
   crd = O_File( name, "r");
   
   molecule->rawmol.name = calloc( strlen(name)+1, sizeof(char));
   strcpy( molecule->rawmol.name, name );
   
   molecule->coor.pbc_flag = 0;
   molecule->coor.pbc = NULL;
   nlines=0;
   molecule->ext_flag = 0;
   molecule->nato = -1;
   memset( buffer, '\0', sizeof(buffer));
   if(fgets(buffer, 150, crd)==NULL && ( !feof(crd) || ferror(crd) ))
   {
     fprintf(stderr, "Warning! Premature end of file reached!\n");
   }
   while( !feof(crd))
   {
    nlines++;
    if( buffer[0] != '*' )
    {
     if( molecule->nato < 0 )
     {
       sscanf( buffer, "%d", &nato);
       if( strstr( buffer, "EXT") != NULL)
         molecule->ext_flag = 1;
     }
     ++molecule->nato;
    }
    else
      crdnatofound = 1;
    if(fgets(buffer, 150, crd)==NULL && ( !feof(crd) || ferror(crd) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
   }
   rewind(crd);
   
   if( crdnatofound == 0 )
   {
     fprintf( stderr, "Number of atoms (* field) not found! maybe %s not a crd?\n", name);
     exit(0);
   }
   
   if( molecule->nato != nato )
   {
     fprintf( stderr, "Problems while reading %s: check number of atoms and headers\n", name);
     exit(0);
   }
   
   if( molecule->ext_flag)
     printf("it's extended!\n");
   
   molecule->coor.nato = molecule->nato;
  // memalloc for coor and pdbs structures
   molecule->coor.cords =  walloc ( 3*molecule->nato, sizeof(float));
  // compatibility pointers
   molecule->coor.xcoor =  molecule->coor.cords;
   molecule->coor.ycoor = &molecule->coor.cords[molecule->nato];
   molecule->coor.zcoor = &molecule->coor.cords[2*molecule->nato];
   molecule->coor.coors =  calloc( 3, sizeof(float *));
   molecule->coor.coors[0] =  molecule->coor.cords;
   molecule->coor.coors[1] = &molecule->coor.cords[molecule->nato];
   molecule->coor.coors[2] = &molecule->coor.cords[2*molecule->nato];

   molecule->rawmol.hetatm  = calloc ( molecule->nato, sizeof(short) );
   molecule->rawmol.atomn   = calloc ( molecule->nato, sizeof( int ) );
   molecule->rawmol.chainId = calloc ( molecule->nato, sizeof( char) );
//   molecule->rawmol.cpresn  = calloc ( molecule->nato, sizeof( int ) );
//   molecule->rawmol.ipresn  = calloc ( molecule->nato, sizeof( int ) );
   molecule->rawmol.presn   = calloc ( molecule->nato, sizeof( int ) );
   molecule->rawmol.resn    = calloc ( molecule->nato, sizeof( int ) );
   molecule->rawmol.weight  = calloc ( molecule->nato, sizeof(float) );

   molecule->rawmol.altLoc  = calloc ( molecule->nato, sizeof( char) );
   molecule->rawmol.iCode   = calloc ( molecule->nato, sizeof( char) );
   molecule->rawmol.occup   = calloc ( molecule->nato, sizeof(float) );
   molecule->rawmol.bFac    = calloc ( molecule->nato, sizeof(float) );
   molecule->rawmol.segend  = calloc ( molecule->nato, sizeof(short) );
   molecule->rawmol.segbeg  = calloc ( molecule->nato, sizeof(short) );
   
   if( molecule->ext_flag == 0 )
   {
     molecule->rawmol.atmtype    = calloc ( molecule->nato, sizeof(char *));
     for ( ii=0; ii < molecule->nato ; ii++)
       molecule->rawmol.atmtype[ii] = calloc( 4, sizeof(char));

     molecule->rawmol.atmId      = calloc ( molecule->nato, sizeof(char *));
     for ( ii=0; ii < molecule->nato ; ii++)
       molecule->rawmol.atmId[ii]   = calloc( 5, sizeof(char));

     molecule->rawmol.restype    = calloc ( molecule->nato, sizeof(char *));
     for ( ii=0; ii < molecule->nato ; ii++)
       molecule->rawmol.restype[ii] = calloc( 4, sizeof(char));

     molecule->rawmol.segId      = calloc ( molecule->nato, sizeof(char *));
     for ( ii=0; ii < molecule->nato ; ii++)
       molecule->rawmol.segId[ii]   = calloc( 4, sizeof(char));

     molecule->rawmol.element    = calloc ( molecule->nato, sizeof(char *));
     for ( ii=0; ii < molecule->nato ; ii++)
       molecule->rawmol.element[ii] = calloc( 2, sizeof(char));

     molecule->rawmol.charge     = calloc ( molecule->nato, sizeof(char *));
     for ( ii=0; ii < molecule->nato ; ii++)
       molecule->rawmol.charge[ii] = calloc( 2, sizeof(char));
   }
   else
   {
     molecule->rawmol.atmtype    = calloc ( molecule->nato, sizeof(char *));
     for ( ii=0; ii < molecule->nato ; ii++)
       molecule->rawmol.atmtype[ii] = calloc( 10, sizeof(char));

     molecule->rawmol.atmId      = calloc ( molecule->nato, sizeof(char *));
     for ( ii=0; ii < molecule->nato ; ii++)
       molecule->rawmol.atmId[ii]   = calloc( 10, sizeof(char));

     molecule->rawmol.restype    = calloc ( molecule->nato, sizeof(char *));
     for ( ii=0; ii < molecule->nato ; ii++)
       molecule->rawmol.restype[ii] = calloc( 8, sizeof(char));

     molecule->rawmol.segId      = calloc ( molecule->nato, sizeof(char *));
     for ( ii=0; ii < molecule->nato ; ii++)
       molecule->rawmol.segId[ii]   = calloc( 8, sizeof(char));

     molecule->rawmol.element    = calloc ( molecule->nato, sizeof(char *));
     for ( ii=0; ii < molecule->nato ; ii++)
       molecule->rawmol.element[ii] = calloc( 2, sizeof(char));

     molecule->rawmol.charge     = calloc ( molecule->nato, sizeof(char *));
     for ( ii=0; ii < molecule->nato ; ii++)
       molecule->rawmol.charge[ii] = calloc( 2, sizeof(char));
   }

// fill up the pbds and coor structures (simple serial reading) - only time data is acquired from source (pdb file)
   
   for( ii=0; ii<(nlines-nato); ii++)
    if(fgets(buffer, 150, crd)==NULL && ( !feof(crd) || ferror(crd) ))
    {
       fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
   for ( ii=0; ii<molecule->nato; )
   {
    memset( buffer, '\0', sizeof(buffer));
    if(fgets(buffer, 150, crd)==NULL && ( !feof(crd) || ferror(crd) ))
    {
       fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
    if ( buffer[0] != '*' )
    {
      if( molecule->ext_flag == 0 )
      {
        molecule->rawmol.hetatm[ii] = 0;
        copyvar( buffer, 0, 5, varbuffer, 12);
        sscanf( varbuffer, "%5s", molecule->rawmol.atmId[ii]);
        copyvar( buffer, 5, 10, varbuffer, 12);
        sscanf( varbuffer, "%d", &molecule->rawmol.presn[ii]);
        copyvar( buffer, 11, 15, varbuffer, 12);
        sscanf( varbuffer, "%4s", molecule->rawmol.restype[ii]);
        copyvar( buffer, 16, 20, varbuffer, 12);
        sscanf( varbuffer, "%4s", molecule->rawmol.atmtype[ii]);
        copyvar( buffer, 20, 30, varbuffer, 12);
        sscanf( varbuffer, "%10f", &molecule->coor.cords[ii]);
        copyvar( buffer, 30, 40, varbuffer, 12);
        sscanf( varbuffer, "%10f", &molecule->coor.cords[ii+molecule->nato]);
        copyvar( buffer, 40, 50, varbuffer, 12);
        sscanf( varbuffer, "%10f", &molecule->coor.cords[ii+2*molecule->nato]);
        copyvar( buffer, 51, 55, varbuffer, 12);
        sscanf( varbuffer, "%4s", molecule->rawmol.segId[ii]);
//      copyvar( buffer, 51, 55, molecule->rawmol.segId[ii], 4);
        copyvar( buffer, 56, 60, varbuffer, 12);
        sscanf( varbuffer, "%4d", &molecule->rawmol.resn[ii]);
        copyvar( buffer, 60, 70, varbuffer, 12);
        sscanf( varbuffer, "%10f", &molecule->rawmol.weight[ii]);
      }
      else
      {
        copyvar( buffer,   0,  10, varbuffer, 12);
        sscanf( varbuffer, "%10s", molecule->rawmol.atmId[ii]);
        copyvar( buffer,  10,  20, varbuffer, 12);
        sscanf( varbuffer, "%d", &molecule->rawmol.presn[ii]);
        copyvar( buffer,  22,  30, varbuffer, 12);
        sscanf( varbuffer, "%8s", molecule->rawmol.restype[ii]);
        copyvar( buffer,  32,  40, varbuffer, 12);
        sscanf( varbuffer, "%8s", molecule->rawmol.atmtype[ii]);
        copyvar( buffer,  40,  60, varbuffer, 12);
/*      sscanf( varbuffer, "%f", &molecule->coor.xcoor[ii]);    */
        sscanf( varbuffer, "%f", &molecule->coor.cords[ii]);
        copyvar( buffer,  60,  80, varbuffer, 12);
/*      sscanf( varbuffer, "%f", &molecule->coor.ycoor[ii]);    */
        sscanf( varbuffer, "%f", &molecule->coor.cords[ii+molecule->nato]);
        copyvar( buffer,  80, 100, varbuffer, 12);
/*      sscanf( varbuffer, "%f", &molecule->coor.zcoor[ii]);    */
        sscanf( varbuffer, "%f", &molecule->coor.cords[ii+2*molecule->nato]);
        copyvar( buffer, 102, 110, varbuffer, 12);
        sscanf( varbuffer, "%8s", molecule->rawmol.segId[ii]);
//      copyvar( buffer, 51, 55, molecule->rawmol.segId[ii], 4);
        copyvar( buffer, 112, 120, varbuffer, 12);
        sscanf( varbuffer, "%d", &molecule->rawmol.resn[ii]);
        copyvar( buffer, 120, 140, varbuffer, 12);
        sscanf( varbuffer, "%f", &molecule->rawmol.weight[ii]);
      }
      
      molecule->rawmol.atomn[ii] = ii+1;
      
      if ( !strncmp(molecule->rawmol.atmId[ii], "*****", 5) && molecule->rawmol.atomn[ii] != atoi(molecule->rawmol.atmId[ii]) )
      {
        //fprintf( stderr, "Error in atom numbering (#%d - %d (%s))\n", molecule->rawmol.atomn[ii], atoi(molecule->rawmol.atmId[ii]), molecule->rawmol.atmId[ii] );
        //exit(99);
        fprintf( stderr, "Warning: mismatch in atom numbering (#%d - %d (%s))\n", molecule->rawmol.atomn[ii], atoi(molecule->rawmol.atmId[ii]), molecule->rawmol.atmId[ii] );
      }
      if( molecule->rawmol.segId[ii][0] != '\0' )
        molecule->rawmol.chainId[ii] = molecule->rawmol.segId[ii][0];
      else
      {
        molecule->rawmol.chainId[ii] = 'A';
        molecule->rawmol.segId[ii][0] = 'A';
      }
      molecule->rawmol.bFac[ii] = molecule->rawmol.weight[ii];
      ii++;
    }
   }

  Raw2Struct( molecule );
  
  molecule->coor.pbc_flag = 0;
  molecule->coor.pbc = NULL;
  if( molecule->rawmol.b_cryst == 1 )
  {
    molecule->a =       molecule->rawmol.a;
    molecule->b =       molecule->rawmol.b;
    molecule->c =       molecule->rawmol.c;
    molecule->alpha =   molecule->rawmol.alpha;
    molecule->beta =    molecule->rawmol.beta;
    molecule->gamma =   molecule->rawmol.gamma;
    strcpy( molecule->s_group, molecule->rawmol.s_group);
    molecule->z =       molecule->rawmol.z;
    
    molecule->coor.pbc_flag    = 1;
    molecule->coor.pbc = calloc( 1, sizeof( Pbc ));
    molecule->coor.pbc->a_size = molecule->rawmol.a;
    molecule->coor.pbc->b_size = molecule->rawmol.b;
    molecule->coor.pbc->c_size = molecule->rawmol.c;
    molecule->coor.pbc->angle1 = cos(molecule->rawmol.alpha);
    molecule->coor.pbc->angle2 = cos(molecule->rawmol.beta);
    molecule->coor.pbc->angle3 = cos(molecule->rawmol.gamma);
    strcpy( molecule->s_group, molecule->rawmol.s_group);
    molecule->coor.pbc->zval       = molecule->rawmol.z;
  }
  
  /*
// to fill up segment structures, gather data on the molecule(s)   

//  Count the number of segments
   molecule->nSeg=1;
   for ( ii=0; ii<molecule->nato-1; ii++)
   {
     if ( strncmp( molecule->rawmol.segId[ii], molecule->rawmol.segId[ii+1], 4) )
     {
       molecule->nSeg ++;
     }
   }

// Allocate space for the segments and set in each the number of atoms/residues to 1
   molecule->segment = NULL;
   molecule->segment = calloc ( molecule->nSeg, sizeof ( Segm ));
   if ( molecule->segment == NULL )
    printf(" ?!? Did not actually calloc!\n");
    
   molecule->nApS    = calloc ( molecule->nSeg, sizeof (  int  ));
   molecule->nRpS    = calloc ( molecule->nSeg, sizeof (  int  ));
   for ( ii=0; ii<molecule->nSeg; ii++)
   {
     molecule->segment[ii].nApS = 1;
     molecule->segment[ii].nRpS = 1;
     molecule->nApS[ii] = 1;
     molecule->nRpS[ii] = 1;
   }
   
//  Count the residues and atoms in each segment; get segment names
   jj = 0;
   StringCp ( molecule->rawmol.segId[0], molecule->segment[jj].segName, 4);
   for ( ii=0; ii<molecule->nato-1; ii++)
   {
     if ( !strncmp(molecule->rawmol.segId[ii], molecule->rawmol.segId[ii+1], 4) )
     {
       molecule->nApS[jj]++;
       molecule->segment[jj].nApS++;
     }
     if ( !strncmp(molecule->rawmol.segId[ii], molecule->rawmol.segId[ii+1], 4) && (molecule->rawmol.resn[ii] != molecule->rawmol.resn[ii+1]) )
     {
       molecule->nRpS[jj]++;
       molecule->segment[jj].nRpS++;
     }
     if (  strncmp(molecule->rawmol.segId[ii], molecule->rawmol.segId[ii+1], 4) )
     {
       jj++;
       StringCp ( molecule->rawmol.segId[ii+1], molecule->segment[jj].segName, 4);
     }
   }
   // count overall residue number and put it in molecule->nRes
   molecule->nRes=0;
   for( ii=0; ii<molecule->nSeg; ii++)
     molecule->nRes += molecule->nRpS[ii];
   
//  Allocate space for all residues and set in each the number of atoms to 1
   for ( ii=0; ii< molecule->nSeg; ii++)
   {
     molecule->segment[ii].pRes = calloc ( molecule->segment[ii].nRpS, sizeof( Res ) );
     for ( jj=0; jj<molecule->segment[ii].nRpS; jj++)
       molecule->segment[ii].pRes[jj].nApR = 1;
   }
   
//  Count atoms in each residue and set res.resType, res.resn and res.resId
   jj = 0;
   kk = 0;
   for ( ii=0; ii<molecule->nato-1; ii++)
    {
     if ( molecule->rawmol.resn[ii] == molecule->rawmol.resn[ii+1] && !strncmp(molecule->rawmol.segId[ii], molecule->rawmol.segId[ii+1], 4) )
      molecule->segment[jj].pRes[kk].nApR++;
     if ( !strncmp(molecule->rawmol.segId[ii], molecule->rawmol.segId[ii+1], 4) && (molecule->rawmol.resn[ii] != molecule->rawmol.resn[ii+1]) )
      {
       StringCp( molecule->rawmol.restype[ii], molecule->segment[jj].pRes[kk].resType, 4);
       molecule->segment[jj].pRes[kk].resn  = molecule->rawmol.resn[ii];
       kk++ ;
      }
     if (  strncmp(molecule->rawmol.segId[ii], molecule->rawmol.segId[ii+1], 4) )
      {
       if ( kk+1 != molecule->segment[jj].nRpS )
        {
         fprintf( stderr, "Error in nRpS for segment # %d (%s - %s)- aborting\n", jj+1, molecule->rawmol.segId[ii], molecule->rawmol.segId[ii+1] );
         fprintf( stderr, "kk+1 was %d while nRpS is %d\n", kk+1, molecule->segment[jj].nRpS );
         exit(98);
        }
       StringCp( molecule->rawmol.restype[ii], molecule->segment[jj].pRes[kk].resType, 4);
       molecule->segment[jj].pRes[kk].resn  = molecule->rawmol.resn[ii];
       jj++;
       kk = 0;
      }
    }
// Let's take care of the last residue of the last segment, too  
   StringCp( molecule->rawmol.restype[molecule->nato-1], molecule->segment[molecule->nSeg-1].pRes[molecule->segment[molecule->nSeg-1].nRpS-1].resType, 4);
   molecule->segment[molecule->nSeg-1].pRes[molecule->segment[molecule->nSeg-1].nRpS-1].resn  = molecule->rawmol.resn[molecule->nato-1];

//  Allocate space for all atoms
   
   for ( ii=0; ii< molecule->nSeg; ii++)
    for ( jj=0; jj< molecule->segment[ii].nRpS; jj++)
     molecule->segment[ii].pRes[jj].pAto = calloc ( molecule->segment[ii].pRes[jj].nApR, sizeof ( Atom ) );

   ll = 0;
   for ( ii=0; ii<molecule->nSeg; ii++ )
   {
    for ( jj=0; jj< molecule->segment[ii].nRpS; jj++)
    {
     for ( kk=0; kk< molecule->segment[ii].pRes[jj].nApR; kk++)
     {
      molecule->segment[ii].pRes[jj].pAto[kk].presn   = molecule->rawmol.presn[ll];
      molecule->segment[ii].pRes[jj].pAto[kk].hetatm   = molecule->rawmol.hetatm[ll];
      molecule->segment[ii].pRes[jj].pAto[kk].atomn   = molecule->rawmol.atomn[ll];
      StringCp( molecule->rawmol.atmId[ll],   molecule->segment[ii].pRes[jj].pAto[kk].atmId,   5);
      StringCp( molecule->rawmol.atmtype[ll], molecule->segment[ii].pRes[jj].pAto[kk].atmType, 5);
      StringCp( molecule->rawmol.segId[ll], molecule->segment[ii].pRes[jj].pAto[kk].segId, 4);
      molecule->segment[ii].pRes[jj].pAto[kk].altLoc  =                        ' ';
      molecule->segment[ii].pRes[jj].pAto[kk].chainId = molecule->rawmol.chainId[ll];
      molecule->segment[ii].pRes[jj].pAto[kk].iCode   =                        ' ';
      molecule->segment[ii].pRes[jj].pAto[kk].xCoor   =   molecule->coor.xcoor[ll];
      molecule->segment[ii].pRes[jj].pAto[kk].yCoor   =   molecule->coor.ycoor[ll];
      molecule->segment[ii].pRes[jj].pAto[kk].zCoor   =   molecule->coor.zcoor[ll];
      molecule->segment[ii].pRes[jj].pAto[kk].occup   =                          0;
      molecule->segment[ii].pRes[jj].pAto[kk].bFac    =                          0;
      molecule->segment[ii].pRes[jj].pAto[kk].weight  =  molecule->rawmol.weight[ll];
      StringCp( molecule->rawmol.element[ll], molecule->segment[ii].pRes[jj].pAto[kk].element, 2);
      StringCp( molecule->rawmol.charge[ll],  molecule->segment[ii].pRes[jj].pAto[kk].charge,  2);

      ll++;
     }
    }
   }
   */

   molecule->origin = 2;
   molecule->filled = calloc( 4, sizeof(char)); // overalloc to avoid overflow
   sprintf( molecule->filled, "YES");
   
   fclose(crd);
   return;
}
// ------------------------------------------------------------------
void WriteCrd ( char *name, Molecule *molecule )
{
   int ii, jj, kk, ll;
   
   FILE   *out;
   
   out = fopen( name, "w");
   fprintf(out, "* \n");
   fprintf(out, "%d\n", molecule->nato);
   
   ll=0;
   if( molecule->ext_flag == 1 )
     for ( ii = 0; ii < molecule->nSeg; ii++ )
       for ( jj = 0; jj < molecule->segment[ii].nRpS; jj++ )
         for ( kk = 0; kk < molecule->segment[ii].pRes[jj].nApR; kk++ )
         {
           fprintf(out, "%10s%10d  %-8s  %-8s%20.10f%20.10f%20.10f  %-8s  %-8d%20.10f\n",
                        molecule->segment[ii].pRes[jj].pAto[kk].atmId,
                        molecule->segment[ii].pRes[jj].pAto[kk].presn,
                        molecule->segment[ii].pRes[jj].resType,
                        molecule->segment[ii].pRes[jj].pAto[kk].atmType,
                        molecule->coor.xcoor[ll],
                        molecule->coor.ycoor[ll],
                        molecule->coor.zcoor[ll],
                        molecule->segment[ii].pRes[jj].pAto[kk].segId,
                        molecule->segment[ii].pRes[jj].resn,
                        molecule->segment[ii].pRes[jj].pAto[kk].weight  );
           ll++;
         }
   else
     for ( ii = 0; ii < molecule->nSeg; ii++ )
       for ( jj = 0; jj < molecule->segment[ii].nRpS; jj++ )
         for ( kk = 0; kk < molecule->segment[ii].pRes[jj].nApR; kk++ )
         {
           fprintf(out, "%5s%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
                        molecule->segment[ii].pRes[jj].pAto[kk].atmId,
                        molecule->segment[ii].pRes[jj].pAto[kk].presn,
                        molecule->segment[ii].pRes[jj].resType,
                        molecule->segment[ii].pRes[jj].pAto[kk].atmType,
                        molecule->coor.xcoor[ll],
                        molecule->coor.ycoor[ll],
                        molecule->coor.zcoor[ll],
                        molecule->segment[ii].pRes[jj].pAto[kk].segId,
                        molecule->segment[ii].pRes[jj].resn,
                        molecule->segment[ii].pRes[jj].pAto[kk].weight  );
           ll++;
         }
  
   fclose (out);
    
   return;
}
// ------------------------------------------------------------------
void WriteCrd_unstr ( char *name, Molecule *molecule )
{
  int     ii;
  FILE   *out;
  
  out = fopen( name, "w");
  fprintf(out, "* \n");
  fprintf(out, "%d\n", molecule->nato);
  
  if( molecule->ext_flag == 1 )
    for ( ii = 0; ii < molecule->nato; ii++ )
      fprintf(out, "%10s%10d  %-8s  %-8s%20.10f%20.10f%20.10f  %-8s  %-8d%20.10f\n",
                    molecule->rawmol.atmId[ii],
                    molecule->rawmol.presn[ii],
                    molecule->rawmol.restype[ii],
                    molecule->rawmol.atmtype[ii],
                    molecule->coor.xcoor[ii],
                    molecule->coor.ycoor[ii],
                    molecule->coor.zcoor[ii],
                    molecule->rawmol.segId[ii],
                    molecule->rawmol.resn[ii],
                    molecule->rawmol.weight[ii]  );
  else
    for ( ii = 0; ii < molecule->nato; ii++ )
      fprintf(out, "%5s%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
                    molecule->rawmol.atmId[ii],
                    molecule->rawmol.presn[ii],
                    molecule->rawmol.restype[ii],
                    molecule->rawmol.atmtype[ii],
                    molecule->coor.xcoor[ii],
                    molecule->coor.ycoor[ii],
                    molecule->coor.zcoor[ii],
                    molecule->rawmol.segId[ii],
                    molecule->rawmol.resn[ii],
                    molecule->rawmol.weight[ii]  );
  
   fclose (out);
   return;
}
// ------------------------------------------------------------------
int ShowMolInfo( Molecule *molecule)
{
  int       ii, jj, kk;
  
  fprintf( stdout, "Nato: %d\n", molecule->nato );
  fprintf( stdout, "nModels: %d\n", molecule->nModels );
  fprintf( stdout, "nChains: %d\n", molecule->nChains );
  fprintf( stdout, "nSegments: %d\n", molecule->nSeg );
  fprintf( stdout, "nRes: %d\n", molecule->nRes );
  fprintf( stdout, "Sequence: %s\n", molecule->seq1 );
  for( ii=0; ii<molecule->nChains; ii++ )
  {
    fprintf( stdout, "Chain #%d - ID: %c\n", ii, molecule->chains[ii].chainId);
    fprintf( stdout, " > number of segments: %d\n", molecule->chains[ii].nSeg );
    for( jj=0; jj<molecule->chains[ii].nSeg; jj++ )
    {
      fprintf( stdout, " >> segment #%d - ID %s\n", jj, molecule->chains[ii].pSegm[jj].segName );
      fprintf( stdout, " >> # of residues: %d\n", molecule->chains[ii].pSegm[jj].nRpS );
      for( kk=0; kk<molecule->chains[ii].pSegm[jj].nRpS; kk++ )
      {
        fprintf( stdout, " >>> residue #%d - ID %s\n", kk, molecule->chains[ii].pSegm[jj].pRes[kk].resType );
        fprintf( stdout, " >>> # of atoms: %d\n", molecule->chains[ii].pSegm[jj].pRes[kk].nApR );
      }
    }
  }
  fprintf( stdout, "total nSegments: %d\n", molecule->nSeg );
  
  return 0;
}
// ------------------------------------------------------------------
void ReadXyz ( char *name, Molecule *molecule )
{
   FILE   *xyzfile;
   int     ii;
   char    buffer[80];
   
   
   xyzfile = O_File( name, "r");
   
// first, get how many atoms are in the structure file
   molecule->nato = -1;
   molecule->rawmol.b_cryst=0;
   memset( buffer, '\0', sizeof(buffer));
   while ( !feof(xyzfile))
   {
    if ( strncmp(buffer, "XYZ", 3) )
       ++molecule->nato;
    if(fgets(buffer, 80, xyzfile)==NULL && ( !feof(xyzfile) || ferror(xyzfile) ))
    {
       fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
   }
   rewind(xyzfile);
   
// memalloc for coor and pdbs structures
   molecule->coor.cords =  walloc ( 3*molecule->nato, sizeof(float));
   molecule->coor.xcoor =  molecule->coor.cords;
   molecule->coor.ycoor = &molecule->coor.cords[molecule->nato];
   molecule->coor.zcoor = &molecule->coor.cords[2*molecule->nato];
   molecule->coor.coors =  calloc( 3, sizeof(float *));
   molecule->coor.coors[0] =  molecule->coor.cords;
   molecule->coor.coors[1] = &molecule->coor.cords[molecule->nato];
   molecule->coor.coors[2] = &molecule->coor.cords[2*molecule->nato];
   
// molecule->coor.xcoor = calloc ( molecule->nato, sizeof(float) );
// molecule->coor.ycoor = calloc ( molecule->nato, sizeof(float) );
// molecule->coor.zcoor = calloc ( molecule->nato, sizeof(float) );

   molecule->rawmol.atomn   = calloc ( molecule->nato, sizeof( int ) );
   molecule->rawmol.chainId = calloc ( molecule->nato, sizeof( char) );
   molecule->rawmol.altLoc  = calloc ( molecule->nato, sizeof( char) );
   molecule->rawmol.iCode   = calloc ( molecule->nato, sizeof( char) );
   molecule->rawmol.resn    = calloc ( molecule->nato, sizeof( int ) );
   molecule->rawmol.occup   = calloc ( molecule->nato, sizeof(float) );
   molecule->rawmol.bFac    = calloc ( molecule->nato, sizeof(float) );
   
   // this two just for crd-compatibility
   molecule->rawmol.presn   = calloc ( molecule->nato, sizeof(float) );
   molecule->rawmol.weight  = calloc ( molecule->nato, sizeof(float) );

   molecule->rawmol.atmtype    = calloc ( molecule->nato, sizeof(char *));
   molecule->rawmol.atmtype[0] = calloc ( 5*molecule->nato, sizeof(char) );
   for ( ii=0; ii < molecule->nato ; ii++)
    molecule->rawmol.atmtype[ii] = molecule->rawmol.atmtype[0] + 5*ii;

   molecule->rawmol.atmId      = calloc ( molecule->nato, sizeof(char *));
   molecule->rawmol.atmId[0]   = calloc ( 5*molecule->nato, sizeof(char) );
   for ( ii=0; ii < molecule->nato ; ii++)
    molecule->rawmol.atmId[ii]   = molecule->rawmol.atmId[0] + 5*ii;

   molecule->rawmol.restype    = calloc ( molecule->nato, sizeof(char *));
   molecule->rawmol.restype[0] = calloc ( 4*molecule->nato, sizeof(char) );
   for ( ii=0; ii < molecule->nato ; ii++)
    molecule->rawmol.restype[ii] = molecule->rawmol.restype[0] + 4*ii;

   molecule->rawmol.segId      = calloc ( molecule->nato, sizeof(char *));
   molecule->rawmol.segId[0]   = calloc ( 4*molecule->nato, sizeof(char) );
   for ( ii=0; ii < molecule->nato ; ii++)
    molecule->rawmol.segId[ii]   = molecule->rawmol.segId[0] + 4*ii;

   molecule->rawmol.element    = calloc ( molecule->nato, sizeof(char *));
   molecule->rawmol.element[0] = calloc ( 2*molecule->nato, sizeof(char) );
   for ( ii=0; ii < molecule->nato ; ii++)
    molecule->rawmol.element[ii] = molecule->rawmol.element[0] + 2*ii;

   molecule->rawmol.charge     = calloc ( molecule->nato, sizeof(char *));
   molecule->rawmol.charge[0]  = calloc ( 2*molecule->nato, sizeof(char) );
   for ( ii=0; ii < molecule->nato ; ii++)
    molecule->rawmol.charge[ii] = molecule->rawmol.charge[0] + 2*ii;

// no periodic boundary conditions   
   molecule->coor.pbc = 0;
   
// fill up the pbds and coor structures (simple serial reading) - only time data is acquired from source (pdb file)
   for ( ii=0; ii<molecule->nato; )
   {
     if(fgets(buffer, 80, xyzfile)==NULL && ( !feof(xyzfile) || ferror(xyzfile) ))
     {
        fprintf(stderr, "Warning! Premature end of file reached!\n");
     }
     if ( strncmp(buffer, "XYZ", 3) )
     {
       sscanf(buffer, "%f %f %f", &molecule->coor.xcoor[ii], &molecule->coor.ycoor[ii], &molecule->coor.zcoor[ii]);
       molecule->rawmol.atomn[ii] = ii+1;
       sprintf( molecule->rawmol.atmId[ii], "%5d", molecule->rawmol.atomn[ii]);
       ii++;
     }
   }
   
// build a totally fake molecular structure   

//  set the number of segments = 1
   molecule->nSeg=1;

   molecule->nApS    = calloc ( molecule->nSeg, sizeof (  int  ));
   molecule->nRpS    = calloc ( molecule->nSeg, sizeof (  int  ));
   molecule->nApS[0] = molecule->nato;
   molecule->nRpS[0] = 1;
   
// Allocate space for the segment and set the number of atoms = nato and residues to 1
   molecule->segment = NULL;
   molecule->segment = calloc ( molecule->nSeg, sizeof ( Segm ));
   if ( molecule->segment == NULL )
    printf(" ?!? Did not actually calloc!\n");
   sprintf( molecule->segment[0].segName, "DUMM" );
   molecule->segment[0].nApS = molecule->nato;
   molecule->segment[0].nRpS = 1;
   
   molecule->segment[0].pRes = calloc ( 1, sizeof( Res ) );
   molecule->segment[0].pRes[0].resn = 1;
   molecule->segment[0].pRes[0].nApR = molecule->nato;
   sprintf( molecule->segment[0].pRes[0].resType, "DUM");
   
   molecule->segment[0].pRes[0].pAto = calloc ( molecule->nato, sizeof ( Atom ) );
   for( ii=0; ii< molecule->nato; ii++ )
   {
     molecule->segment[0].pRes[0].pAto[ii].atomn = ii+1;
     sprintf( molecule->segment[0].pRes[0].pAto[ii].atmType, "Dumm");
     molecule->segment[0].pRes[0].pAto[ii].xCoor = molecule->coor.xcoor[ii];
     molecule->segment[0].pRes[0].pAto[ii].yCoor = molecule->coor.ycoor[ii];
     molecule->segment[0].pRes[0].pAto[ii].zCoor = molecule->coor.zcoor[ii];
     sprintf( molecule->segment[0].pRes[0].pAto[ii].atmId, "%6d", ii);
     sprintf( molecule->segment[0].pRes[0].pAto[ii].segId, "DUMM");
     sprintf( molecule->segment[0].pRes[0].pAto[ii].element, "DU");
     sprintf( molecule->segment[0].pRes[0].pAto[ii].charge, " 0");
     molecule->segment[0].pRes[0].pAto[ii].chainId = 'A';
     molecule->segment[0].pRes[0].pAto[ii].altLoc = ' ';
     molecule->segment[0].pRes[0].pAto[ii].iCode = ' ';
     molecule->segment[0].pRes[0].pAto[ii].bFac = 0.00;
     molecule->segment[0].pRes[0].pAto[ii].occup = 0.00;
     molecule->segment[0].pRes[0].pAto[ii].weight = 1.0;
     molecule->segment[0].pRes[0].pAto[ii].presn = ii+1;
     molecule->segment[0].pRes[0].pAto[ii].bb = 1;
   }

   molecule->origin = 2;
   
   fclose(xyzfile);
   
   return;
}
// ------------------------------------------------------------------
void WriteXyz (  char *name, Molecule *molecule )
{
   int   ii;
   FILE   *out;
   
   out = fopen( name, "w");
   fprintf ( out, "XYZ 0\n");
   for ( ii=0; ii<molecule->nato; ii++ )
     fprintf ( out, "%8.3f %8.3f %8.3f\n", molecule->coor.xcoor[ii], molecule->coor.ycoor[ii], molecule->coor.zcoor[ii]);
   
   return;
}
// ------------------------------------------------------------------
// ------------------------------------------------------------------
Molecule * ReadMolecule( char *name, moltype filetype )
{
   Molecule    *molecule;
   
   if( filetype == 0 )
     filetype = wrd_whichmoltype(name);
   
   molecule = calloc( 1, sizeof(Molecule));
   
   molecule->ext_flag = 0;
   molecule->coor.pbc = NULL;
   molecule->coor.pbc_flag = 0;
   molecule->filled = calloc( 4, sizeof(char));
   memset( molecule->filled, '\0', 4);

   switch ( filetype )
   {
    case 0:                     // get filetype from name; case should not be called
     fprintf( stderr, ".OOO files not supported (?)\n");
     exit(0);
    case 1:                     // read a pdb
     ReadPdb ( name, molecule );
     molecule->origin = 0;
     break;
    case 2:                     // read a crd
     ReadCrd ( name, molecule );
     molecule->origin = 1;
     break;
    case 3:                     // read a cor (crd)
     ReadCrd ( name, molecule );
     molecule->origin = 1;
     break;
    case 4:                     // read an zxyz ?
     ReadXyz ( name, molecule );
     break;
    case 5:                     // read a gro
     fprintf( stderr, "gro reading not implemented yet\n");
     exit(0);
     break;
    default:
     fprintf( stderr, "Error while calling ReadMolecule! check filetype\n");
     exit(78);
     break;
   }
   
   return molecule;
}
// ------------------------------------------------------------------
void WriteMolecule_str( Molecule *molecule, char *name, moltype filetype )
{
   if( filetype == 0)
     filetype = wrd_whichmoltype( name );

   switch ( filetype )
   {
    case 0:                     // get filetype from name; case should not be called
     fprintf( stderr, ".OOO files not supported (?)\n");
     exit(0);
    case 1:                     // write a pdb
     WritePdb ( name, molecule );
     break;
    case 2:                     // write a crd
     WriteCrd ( name, molecule );
     break;
    case 3:                     // write a cor (crd)
     WriteCrd ( name, molecule );
     break;
    case 4:                     // write a crd
     WriteXyz ( name, molecule );
     break;
    case 5:                     // write a gro
     fprintf( stderr, "gro writing not implemented yet\n");
     exit(0);
     break;
    default:
     fprintf( stderr, "Error while calling WriteMolecule! check filetype\n");
     exit(78);
     break;
   }
   
   return;
}
// ------------------------------------------------------------------
void WriteMolecule_unstr( Molecule *molecule, char *name, moltype filetype )
{
   if( filetype == 0)
     filetype = wrd_whichmoltype( name );

   switch ( filetype )
   {
    case 0:                     // get filetype from name; case should not be called
     fprintf( stderr, ".OOO files not supported (?)\n");
     exit(0);
    case 1:                     // write a pdb
     WritePdb_unstr ( name, molecule );
     break;
    case 2:                     // write a crd
     WriteCrd_unstr ( name, molecule );
     break;
    case 3:                     // write a cor (crd)
     WriteCrd_unstr ( name, molecule );
     break;
    case 4:                     // write a crd
     WriteXyz ( name, molecule );
     break;
    case 5:                     // write a gro
     fprintf( stderr, "gro writing not implemented yet\n");
     exit(0);
     break;
    default:
     fprintf( stderr, "Error while calling WriteMolecule_unstr! check filetype\n");
     exit(78);
     break;
   }
   
   return;
}
// ------------------------------------------------------------------
void WriteMolecule( Molecule *molecule, char *name, moltype filetype )
{
  WriteMolecule_unstr( molecule, name, filetype );
  return;
}
// ------------------------------------------------------------------
void GetMolecule( char *name, Molecule *molecule, moltype filetype )
{
   molecule->ext_flag = 0;
   /*
   molecule->coor.pbc = NULL;
   molecule->coor.pbc_flag = 0;
   if( molecule->filled = NULL )
   {
     molecule->filled = calloc( 4, sizeof(char));
     memset( molecule->filled, '\0', 4);
   }
   */
   if( filetype == 0 )
     filetype = wrd_whichmoltype(name);
   
   switch ( filetype )
   {
    case 0:                     // get filetype from name; case should not be called
     fprintf( stderr, ".OOO files not supported (?)\n");
     exit(0);
    case 1:                     // read a pdb
     ReadPdb ( name, molecule );
     molecule->origin = 0;
     break;
    case 2:                     // read a crd
     ReadCrd ( name, molecule );
     molecule->origin = 1;
     break;
    case 3:                     // read a cor (crd)
     ReadCrd ( name, molecule );
     molecule->origin = 1;
     break;
    case 4:                     // read an zxyz ?
     ReadXyz ( name, molecule );
     break;
    case 5:                     // read a gro
     fprintf( stderr, "gro reading not implemented yet\n");
     exit(0);
     break;
    default:
     fprintf( stderr, "Error while calling GetMolecule! check filetype\n");
     exit(78);
     break;
   }
   
   return;
}
// ------------------------------------------------------------------
void PutMolecule( Molecule *molecule , char *name, moltype filetype )
{
   if( filetype == 0)
     filetype = wrd_whichmoltype( name );
   switch ( filetype )
   {
    case 0:                     // get filetype from name; case should not be called
     fprintf( stderr, ".OOO files not supported (?)\n");
     exit(0);
    case 1:                     // write a pdb
     WritePdb ( name, molecule );
     break;
    case 2:                     // write a crd
     WriteCrd ( name, molecule );
     break;
    case 3:                     // write a cor (crd)
     WriteCrd ( name, molecule );
     break;
    case 4:                     // write a crd
     WriteXyz ( name, molecule );
     break;
    case 5:                     // write a gro
     fprintf( stderr, "gro writing not implemented yet\n");
     exit(0);
     break;
    default:
     fprintf( stderr, "Error while calling PutMolecule! check filetype\n");
     exit(78);
     break;
   }
   
   UpDateMol( molecule );
   
   return;
}
// ------------------------------------------------------------------
// ------------------------------------------------------------------
void UpDateMol( Molecule *molecule )
{
   int ii, jj, kk, presn;
   
   if( molecule->origin == 0 )
   {
    presn=0;
    for(ii=0; ii<molecule->nSeg; ii++ )
     for ( jj = 0; jj < molecule->segment[ii].nRpS; jj++ )
     {
      presn ++;
      for ( kk = 0; kk < molecule->segment[ii].pRes[jj].nApR; kk++ )
       molecule->segment[ii].pRes[jj].pAto[kk].presn = presn;
     }
   }
   
   return;
}
// ------------------------------------------------------------------
void CleanMolecule( Molecule *molecule )
{
   int  ii, jj;
   
   free( molecule->coor.cords );
   free( molecule->coor.coors );

   free( molecule->rawmol.atomn    );
   free( molecule->rawmol.chainId  );
   free( molecule->rawmol.presn    );
   free( molecule->rawmol.resn    );
   free( molecule->rawmol.weight  );

   free( molecule->rawmol.altLoc  );
   free( molecule->rawmol.iCode   );
   free( molecule->rawmol.occup   );
   free( molecule->rawmol.bFac   );
   
   //free( molecule->pRes );
   
   for ( ii=0; ii < molecule->nato ; ii++)
    free( molecule->rawmol.atmtype[ii] );
   free( molecule->rawmol.atmtype );

   for ( ii=0; ii < molecule->nato ; ii++)
    free( molecule->rawmol.atmId[ii] );
   free( molecule->rawmol.atmId );

   for ( ii=0; ii < molecule->nato ; ii++)
    free( molecule->rawmol.restype[ii] );
   free( molecule->rawmol.restype );

   for ( ii=0; ii < molecule->nato ; ii++)
    free( molecule->rawmol.segId[ii] );
   free( molecule->rawmol.segId );

   for ( ii=0; ii < molecule->nato ; ii++)
    free( molecule->rawmol.element[ii] );
   free( molecule->rawmol.element );

   for ( ii=0; ii < molecule->nato ; ii++)
    free( molecule->rawmol.charge[ii] );
   free( molecule->rawmol.charge );

   for ( ii=0; ii< molecule->nSeg; ii++)
     for ( jj=0; jj< molecule->segment[ii].nRpS; jj++)
       free( molecule->segment[ii].pRes[jj].pAto );
   
   for ( ii=0; ii< molecule->nSeg; ii++)
     free( molecule->segment[ii].pRes );
   
   free( molecule->segment );
   free( molecule->nApS );
   free( molecule->nRpS );
   
   free( molecule->nApR );
   free( molecule->pRes );
   
   if( molecule->coor.pbc_flag == 1 )
   {
     free(molecule->coor.pbc);
     molecule->coor.pbc_flag = 0;
     molecule->coor.pbc = NULL;
   }
   
   molecule->nato = 0;
   molecule->nSeg = 0;
   molecule->nChains = 0;
   molecule->nModels = 0;
   
   sprintf( molecule->filled, "NO ");
   free( molecule->filled );
   
   return;
}
// ------------------------------------------------------------------
void DelMolecule( Molecule *molecule )
{

   CleanMolecule( molecule );
   free(molecule);
   
   return;
}
// ------------------------------------------------------------------
Molecule * CopyMol( Molecule *mol1 )
{
   Molecule *mol2;
   int   ii, jj, kk;
   int   previouschainssegs;
   
   mol2 = calloc( 1, sizeof(Molecule));
   
   mol2->nato = mol1->nato;

   mol2->coor.cords = calloc ( mol2->nato*3, sizeof(float) );
   
   mol2->coor.xcoor = mol2->coor.cords;
   mol2->coor.ycoor = &mol2->coor.cords[mol2->nato];
   mol2->coor.zcoor = &mol2->coor.cords[mol2->nato * 2];
   
   mol2->rawmol.atomn   = calloc ( mol2->nato, sizeof( int ) );
   mol2->rawmol.chainId = calloc ( mol2->nato, sizeof( char) );
   mol2->rawmol.altLoc  = calloc ( mol2->nato, sizeof( char) );
   mol2->rawmol.iCode   = calloc ( mol2->nato, sizeof( char) );
   mol2->rawmol.resn    = calloc ( mol2->nato, sizeof( int ) );
   mol2->rawmol.occup   = calloc ( mol2->nato, sizeof(float) );
   mol2->rawmol.bFac    = calloc ( mol2->nato, sizeof(float) );
   mol2->rawmol.hetatm  = calloc ( mol2->nato, sizeof(short) );

   mol2->rawmol.atmtype    = calloc ( mol2->nato, sizeof(char *));
   mol2->rawmol.atmtype[0] = calloc ( 5*mol2->nato, sizeof(char) );
   for ( ii=0; ii < mol2->nato ; ii++)
    mol2->rawmol.atmtype[ii] = mol2->rawmol.atmtype[0] + 5*ii;

   mol2->rawmol.atmId      = calloc ( mol2->nato, sizeof(char *));
   mol2->rawmol.atmId[0]   = calloc ( 5*mol2->nato, sizeof(char) );
   for ( ii=0; ii < mol2->nato ; ii++)
    mol2->rawmol.atmId[ii]   = mol2->rawmol.atmId[0] + 5*ii;

   mol2->rawmol.restype    = calloc ( mol2->nato, sizeof(char *));
   mol2->rawmol.restype[0] = calloc ( 4*mol2->nato, sizeof(char) );
   for ( ii=0; ii < mol2->nato ; ii++)
    mol2->rawmol.restype[ii] = mol2->rawmol.restype[0] + 4*ii;
    

   mol2->rawmol.segId      = calloc ( mol2->nato, sizeof(char *));
   mol2->rawmol.segId[0]   = calloc ( 5*mol2->nato, sizeof(char) );
   for ( ii=0; ii < mol2->nato ; ii++)
    mol2->rawmol.segId[ii]   = mol2->rawmol.segId[0] + 5*ii;

   mol2->rawmol.element    = calloc ( mol2->nato, sizeof(char *));
   mol2->rawmol.element[0] = calloc ( 3*mol2->nato, sizeof(char) );
   for ( ii=0; ii < mol2->nato ; ii++)
    mol2->rawmol.element[ii] = mol2->rawmol.element[0] + 3*ii;

   mol2->rawmol.charge     = calloc ( mol2->nato, sizeof(char *));
   mol2->rawmol.charge[0]  = calloc ( 2*mol2->nato, sizeof(char) );
   for ( ii=0; ii < mol2->nato ; ii++)
    mol2->rawmol.charge[ii] = mol2->rawmol.charge[0] + 2*ii;

   for ( ii=0; ii<mol2->nato; ii++)
   {
     mol2->rawmol.atomn[ii] = mol1->rawmol.atomn[ii];
     mol2->rawmol.chainId[ii] = mol1->rawmol.chainId[ii];
     mol2->rawmol.altLoc[ii]  = mol2->rawmol.altLoc[ii];
     mol2->rawmol.iCode[ii]   = mol1->rawmol.iCode[ii];
     mol2->rawmol.resn[ii]    = mol1->rawmol.resn[ii];
     mol2->rawmol.occup[ii]   = mol1->rawmol.occup[ii];
     mol2->rawmol.bFac[ii]    = mol1->rawmol.bFac[ii];
     mol2->rawmol.hetatm[ii]  = mol1->rawmol.hetatm[ii];
     sprintf( mol2->rawmol.atmtype[ii], "%s",mol2->rawmol.atmtype[ii]);
     sprintf( mol2->rawmol.atmId[ii], "%s", mol1->rawmol.atmId[ii]);
     sprintf( mol2->rawmol.restype[ii], "%s",mol1->rawmol.restype[ii]);
     sprintf( mol2->rawmol.segId[ii], "%s",mol1->rawmol.segId[ii]);
     sprintf( mol2->rawmol.element[ii], "%s",mol1->rawmol.element[ii]);
     sprintf( mol2->rawmol.charge[ii], "%s",mol2->rawmol.charge[ii]);
   }
   
   
   mol2->nRes = mol1->nRes;
   mol2->nSeg = mol1->nSeg;
   mol2->segment = calloc ( mol2->nSeg, sizeof ( Segm ));
   mol2->nApS    = calloc ( mol2->nSeg, sizeof (  int  ));
   mol2->nRpS    = calloc ( mol2->nSeg, sizeof (  int  ));
   mol2->nApR    = calloc ( mol2->nRes, sizeof (  int  ));
   
   mol2->nChains = mol1->nChains;
   mol2->chains = calloc( mol2->nChains, sizeof( wrd_Chain ));
   mol2->nApC = calloc( mol2->nChains, sizeof( int )) ;
   mol2->nRpC = calloc( mol2->nChains, sizeof( int )) ;
   mol2->nSpC = calloc( mol2->nChains, sizeof( int )) ;
   
   previouschainssegs = 0;
   for( ii=0; ii<mol2->nChains; ii++ )
   {
     mol2->chains[ii].chainId = mol1->chains[ii].chainId;
     mol2->chains[ii].nApC = mol1->chains[ii].nApC;
     mol2->chains[ii].nRpC = mol1->chains[ii].nRpC;
     mol2->chains[ii].nSeg = mol1->chains[ii].nSeg;
     mol2->nApC[ii]   = mol1->nApC[ii];
     mol2->nRpC[ii]   = mol1->nApC[ii];
     mol2->nSpC[ii]   = mol1->nApC[ii];
     mol2->chains[ii].pSegm =  &mol2->segment[previouschainssegs];
     previouschainssegs += mol2->chains[ii].nSeg;

   }
   for( ii=0; ii<mol2->nSeg; ii++ )
   {
     mol2->nApS[ii] = mol1->nApS[ii];
     mol2->nRpS[ii] = mol1->nRpS[ii];
     sprintf( mol2->segment[ii].segName, "%s", mol1->segment[ii].segName);
     mol2->segment[ii].nApS = mol1->segment[ii].nApS;
     mol2->segment[ii].nRpS = mol1->segment[ii].nRpS;
     mol2->segment[ii].pRes = calloc( mol2->segment[ii].nRpS, sizeof( Res ));
     for( jj=0; jj<mol2->segment[ii].nRpS; jj++ )
     {
       mol2->segment[ii].pRes[jj].nApR = mol1->segment[ii].pRes[jj].nApR;
       mol2->segment[ii].pRes[jj].resn = mol1->segment[ii].pRes[jj].resn;
       sprintf( mol2->segment[ii].pRes[jj].resType, "%s", mol1->segment[ii].pRes[jj].resType);
       mol2->segment[ii].pRes[jj].pAto = calloc( mol2->segment[ii].pRes[jj].nApR, sizeof( Atom ) );
       for( kk=0; kk<mol2->segment[ii].pRes[jj].nApR; kk++ )
       {
         mol2->segment[ii].pRes[jj].pAto[kk].atomn = mol1->segment[ii].pRes[jj].pAto[kk].atomn ;
         mol2->segment[ii].pRes[jj].pAto[kk].xCoor = mol1->segment[ii].pRes[jj].pAto[kk].xCoor ;
         mol2->segment[ii].pRes[jj].pAto[kk].yCoor = mol1->segment[ii].pRes[jj].pAto[kk].yCoor ;
         mol2->segment[ii].pRes[jj].pAto[kk].zCoor = mol1->segment[ii].pRes[jj].pAto[kk].zCoor ;
         mol2->segment[ii].pRes[jj].pAto[kk].chainId = mol1->segment[ii].pRes[jj].pAto[kk].chainId ;
         mol2->segment[ii].pRes[jj].pAto[kk].altLoc = mol1->segment[ii].pRes[jj].pAto[kk].altLoc ;
         mol2->segment[ii].pRes[jj].pAto[kk].iCode = mol1->segment[ii].pRes[jj].pAto[kk].iCode ;
         mol2->segment[ii].pRes[jj].pAto[kk].bFac = mol1->segment[ii].pRes[jj].pAto[kk].bFac ;
         mol2->segment[ii].pRes[jj].pAto[kk].occup = mol1->segment[ii].pRes[jj].pAto[kk].occup ;
         mol2->segment[ii].pRes[jj].pAto[kk].weight = mol1->segment[ii].pRes[jj].pAto[kk].weight ;
         mol2->segment[ii].pRes[jj].pAto[kk].presn = mol1->segment[ii].pRes[jj].pAto[kk].presn ;
         mol2->segment[ii].pRes[jj].pAto[kk].bb = mol1->segment[ii].pRes[jj].pAto[kk].bb ;
         sprintf( mol2->segment[ii].pRes[jj].pAto[kk].atmType ,"%s", mol1->segment[ii].pRes[jj].pAto[kk].atmType );
         sprintf( mol2->segment[ii].pRes[jj].pAto[kk].atmId ,"%s", mol1->segment[ii].pRes[jj].pAto[kk].atmId );
         sprintf( mol2->segment[ii].pRes[jj].pAto[kk].segId ,"%s", mol1->segment[ii].pRes[jj].pAto[kk].segId );
         sprintf( mol2->segment[ii].pRes[jj].pAto[kk].element ,"%s", mol1->segment[ii].pRes[jj].pAto[kk].element );
         sprintf( mol2->segment[ii].pRes[jj].pAto[kk].charge ,"%s", mol1->segment[ii].pRes[jj].pAto[kk].charge );
       }
     }
   }
   return mol2;
}
// ------------------------------------------------------------------
//Selection * SeleDiff ( Selection *selection1,  Selection *selection2 )
//{
  //int         ii, jj;
  //Selection  *diffsele;
  //int         found=0;
  
  //diffsele = (Selection *)calloc( 1, sizeof(Selection));
  //diffsele->nselatm = 0;
  //diffsele->selatm = (int *)calloc( selection1->nselatm, sizeof(int));
  //for( ii=0; ii<selection1->nselatm; ii++ )
  //{
    //found = 0;
    //for( jj=0; jj<selection2->nselatm; jj++ )
      //if( selection1->selatm[ii] == selection2->selatm[jj] )
      //{
        //found = 1;
        //break;
      //}
    //if( !found )
    //{
      //diffsele->selatm[diffsele->nselatm] = selection1->selatm[ii];
      //diffsele->nselatm++;
    //}
  //}
  
  //diffsele->selatm = realloc( diffsele->selatm, diffsele->nselatm );
  
  //return diffsele;
//}
// ------------------------------------------------------------------
int SeleDiff ( Selection *selection1,  Selection *selection2, Selection *diffsele )
{
  int         ii, jj;
  int         found=0;
  int        *tmplist;
  
  diffsele->nselatm = 0;
  tmplist = (int *)calloc( selection1->nselatm, sizeof(int));
  for( ii=0; ii<selection1->nselatm; ii++ )
  {
    found = 0;
    for( jj=0; jj<selection2->nselatm; jj++ )
      if( selection1->selatm[ii] == selection2->selatm[jj] )
      {
        found = 1;
        break;
      }
    if( !found )
    {
      tmplist[diffsele->nselatm] = selection1->selatm[ii];
      diffsele->nselatm++;
    }
  }
  
  diffsele->selatm = (int *)calloc( diffsele->nselatm, sizeof(int));
  for( ii=0; ii<diffsele->nselatm; ii++ )
    diffsele->selatm[ii] = tmplist[ii];
  free(tmplist);
  sprintf( diffsele->selestring, "%s", selection1->selestring );
  
  return 0;
}
// ------------------------------------------------------------------
int GetSele ( char *selestring, Selection *selestr, Molecule *molecule )
{
  int          ii, jj;                                                  // Some Iterators
  int          iFirstSeleFlag = 1;                                      // First Selection Flag
  int          iSeleOp = 1;                                             // Selection Operator: 1 means ADD, 0 means DEL
  int          iNumOfSelectedAtoms = 0;                                 // Number of Atoms to Add to Calling Module's Selection
  int         *piSelectedAtoms;                                         // Selected Atoms
  
  char        *pcToken;                                                 // Selection Token
  char         pcSeleOp[5120];                                          // Selection Operator
  char         pcSeleStr[5120];                                         // Selection String
  char         pcOriginalSeleString[10240];                             // Original Module's Selection String

  Selection    seleTmp;                                                 // Temporary Selection Structure
  
  //fprintf( stdout, "DEBUG : %s\n", selestring ); fflush(stdout);
  strcpy(pcOriginalSeleString, selestring);
  
  piSelectedAtoms = (int *) calloc(molecule->nato, sizeof(int));
  
  pcToken = strtok(selestring, ";");
  while(pcToken != NULL)
  {
    if(iFirstSeleFlag == 1)
    {
      // For the first selection, 'add' operand is assumed
      sscanf(pcToken, "%s", pcSeleStr);
      iFirstSeleFlag = 0;
    }
    
    else
    {
      sscanf(pcToken, "%s %s", pcSeleOp, pcSeleStr);
    
      if(strcmp(pcSeleOp, "add") == 0)
        iSeleOp = 1;
        
      else if(strcmp(pcSeleOp, "del") == 0)
        iSeleOp = 0;
      
      else
      {
        printf("ERROR\nInvalid Selection Operator : %s => %s\n", pcToken, pcSeleOp);
        exit(1);
      }
    }
    
    GetSeleAtoms(pcSeleStr, &seleTmp, molecule);
    if(iSeleOp == 1)
    {
      // Add
      for(ii=0; ii<seleTmp.nselatm; ii++)
        piSelectedAtoms[seleTmp.selatm[ii]-1] = 1;
    }
    
    else
    {
      // Del
      for(ii=0; ii<seleTmp.nselatm; ii++)
        piSelectedAtoms[seleTmp.selatm[ii]-1] = 0;
    }
    
    pcToken = strtok(NULL, ";");
  }
  
  iNumOfSelectedAtoms = 0;
  for(ii=0; ii<molecule->nato; ii++)
  {
    if(piSelectedAtoms[ii] == 1)
    {
      iNumOfSelectedAtoms++;
    }
  }
  
  selestr->nselatm = iNumOfSelectedAtoms;
  selestr->selatm = (int *) calloc (selestr->nselatm, sizeof(int));
  jj = -1;
  
  for(ii=0; ii<molecule->nato; ii++)
  {
    if(piSelectedAtoms[ii] == 1)
    {
      jj++;
      selestr->selatm[jj] = (ii + 1);
    }
  }
  
  strcpy(selestr->selestring, pcOriginalSeleString);
  
  free(seleTmp.selatm);
  free(piSelectedAtoms);
  
  return iNumOfSelectedAtoms;
}
// ------------------------------------------------------------------
int GetSeleAtoms ( char *selestring, Selection *selestr, Molecule *molecule )
{
   int   ii, jj, kk, nn;
   char  chastring[1024];
   char  segstring[1024];
   char  resstring[10240];
   char *tmp_resstring;
   char  atmstring[1024];
   char  resnumbst[128];
   char *ext_selestring;
   int  *localsel;
   int   resbufferlentgh;
   
   // some variables to read the atom numbers straight from a file
   FILE *selef;
   int   nslashes;
   
   // === Within variables ===
   int    iBraceOpenFlag=0, iBraceCloseFlag=0;                          // Used to check within syntax
   int    iWithinFlag=0;                                                // 1 if Within selection used
   int    iWithinMode=0;                                                // 0 for inclusive selection, 1 for exclusive selection
   int    iNumOfAtomPairs=0;
   int   *piSele1, *piSele2;                                            // Used by DistanceVectCoor function
   int   *piSeleAtoms;
   int    iVectPos=-1;
   int    iNumOfSeleAtom=0;
   
   float  fShellLength=0.0;                                             // Length of within selection shell
   float *pfDistances;                                                  // Vector of distances
   // ========================
   
   // local variables used to expand selection (residue numbers)
   char *buffer;
   int   ca_found;
   
   // check if selestring is a selection or a file name for the atom list
   nslashes=0;
   ii = 0;
   while( selestring[ii]!='\0' )        // IS THIS RIGHT?!?
   {
    if(selestring[ii] == '/')
     nslashes++;
    ii++;
   }

   if( access( selestring, R_OK ) != -1 && nslashes != strlen(selestring) )
   {
     selef = O_File( selestring, "r");
     buffer=calloc(64, sizeof(char));
     selestr->nselatm = -1;
     while ( !feof(selef))
     {
      if (buffer[0]!='#')
         selestr->nselatm++;
      if(fgets(buffer, 64, selef)==NULL && ( !feof(selef) || ferror(selef) ))
      {
        fprintf(stderr, "Warning! Premature end of file reached!\n");
      }
     }
     rewind(selef);
     selestr->selatm = calloc ( selestr->nselatm, sizeof(int));
     for ( ii=0; ii<selestr->nselatm; )
     {
      if(fgets(buffer, 64, selef)==NULL && ( !feof(selef) || ferror(selef) ))
      {
        fprintf(stderr, "Warning! Premature end of file reached!\n");
      }
      if (buffer[0]!='#')
      {
       sscanf(buffer, "%d", &selestr->selatm[ii]);
       ii++;
      }
     }
     free(buffer);
     return(0);
   }
   
   for(ii=0; ii<strlen(selestring); ii++)
   {
     if(selestring[ii] == '[')
       iBraceOpenFlag = 1;

     if(selestring[ii] == ']')
       iBraceCloseFlag = 1;
   }
   
   memset ( chastring, '\0', sizeof(chastring));
   memset ( segstring, '\0', sizeof(segstring));
   memset ( resstring, '\0', sizeof(resstring));
   memset ( atmstring, '\0', sizeof(atmstring));
   
   ext_selestring = parseSlash( selestring );

   if(iBraceOpenFlag + iBraceCloseFlag == 0)
   {
     iWithinFlag = 0;
     
     sscanf(ext_selestring," /%[^/]/%[^/]/%[^/]/%s", chastring, segstring, resstring, atmstring);
   }
     
   else if(iBraceOpenFlag + iBraceCloseFlag == 2)
   {
     iWithinFlag = 1;
     
     sscanf(ext_selestring," /%[^/]/%[^/]/%[^/]/%[^\[][%f]", chastring, segstring, resstring, atmstring, &fShellLength);
     
     if(fShellLength >= 0.0)
       iWithinMode=0;
     else
     {
       iWithinMode=1;
       fShellLength = fShellLength * -1;
     }
   }
     
   else
   {
     printf("Within Selection is not valid\n");
     exit(1);
   }
   
   selestr->nselatm=0;
   localsel = calloc ( molecule->nato, sizeof(int));
   
// new selection; added srange ability
   
   if ( strchr( resstring, '-' ))
   {
     tmp_resstring = expandrange( resstring, &resbufferlentgh );
     sprintf( resstring, "%s", tmp_resstring);
     free(tmp_resstring);
   }

//  add robustness: look for weird chars in segstring and replace them
  for( ii=0; ii<strlen(segstring); ii++)
    if( segstring[ii] == ',' || segstring[ii] == ';' || segstring[ii] == ' ' )
      segstring[ii] = '|';

   //fprintf( stdout, "DEBUG nChains: %d\n", molecule->nChains ); fflush(stdout);
   
  /* selection loops proper */
  selestr->nselatm = 0;
  for(nn=0; nn<molecule->nChains; nn++)
  {
    if( !fnmatch(chastring, &molecule->chains[nn].chainId, FNM_EXTMATCH) )
      for(ii=0; ii<molecule->chains[nn].nSeg; ii++)
        if( !fnmatch(segstring, molecule->chains[nn].pSegm[ii].segName, FNM_EXTMATCH) )
          for(jj=0; jj<molecule->chains[nn].pSegm[ii].nRpS; jj++)
          {
            sprintf( resnumbst, "%d", molecule->chains[nn].pSegm[ii].pRes[jj].resn);
            if(!fnmatch(resstring, resnumbst, FNM_EXTMATCH) )
            {
              if( !strcmp( atmstring, "#" ) || !strcmp( atmstring, "XXXXX" ) )
              {
                // # calls for an atom per residue - possibly CA, failing that the middle (along the sequence) one
                selestr->nselatm++;
                ca_found = 0;
                for(kk=0; kk<molecule->chains[nn].pSegm[ii].pRes[jj].nApR; kk++)
                {
                  if ( !( fnmatch( "CA", molecule->chains[nn].pSegm[ii].pRes[jj].pAto[kk].atmType, FNM_EXTMATCH)) )
                  {
                   localsel[selestr->nselatm-1] = molecule->chains[nn].pSegm[ii].pRes[jj].pAto[kk].atomn;
                   ca_found = 1;
                   break;
                  }
                }
                if( ca_found == 0 )
                  localsel[selestr->nselatm-1] = molecule->chains[nn].pSegm[ii].pRes[jj].pAto[(int)molecule->chains[nn].pSegm[ii].pRes[jj].nApR/2].atomn;
              }
              else
                for(kk=0; kk<molecule->chains[nn].pSegm[ii].pRes[jj].nApR; kk++)
                {
                  if ( !( fnmatch(atmstring, molecule->chains[nn].pSegm[ii].pRes[jj].pAto[kk].atmType, FNM_EXTMATCH)) )
                  {
                   selestr->nselatm++;
                   localsel[selestr->nselatm-1] = molecule->chains[nn].pSegm[ii].pRes[jj].pAto[kk].atomn;
                  }
                }
            }
          }
  }
  /*
  for(nn=0; nn<molecule->nChains; nn++)
  {
    //for(ii=0; ii<molecule->nSeg; ii++)
    for(ii=0; ii<molecule->chains[nn].nSeg; ii++)
    {
      //for(jj=0; jj<molecule->segment[ii].nRpS; jj++)
      for(jj=0; jj<molecule->chains[nn].pSegm[ii].nRpS; jj++)
      {
        //sprintf( resnumbst, "%d", molecule->segment[ii].pRes[jj].resn);
        sprintf( resnumbst, "%d", molecule->chains[nn].pSegm[ii].pRes[jj].resn);
        //for(kk=0; kk<molecule->segment[ii].pRes[jj].nApR; kk++)
        for(kk=0; kk<molecule->chains[nn].pSegm[ii].pRes[jj].nApR; kk++)
        {
         if ( !( fnmatch(chastring, &molecule->chains[nn].chainId, FNM_EXTMATCH)) &&
              !( fnmatch(segstring, molecule->chains[nn].pSegm[ii].segName, FNM_EXTMATCH)) &&
              !( fnmatch(resstring, resnumbst, FNM_EXTMATCH)) &&
              !( fnmatch(atmstring, molecule->chains[nn].pSegm[ii].pRes[jj].pAto[kk].atmType, FNM_EXTMATCH)) )
          {
           selestr->nselatm++;
           localsel[selestr->nselatm-1] = molecule->chains[nn].pSegm[ii].pRes[jj].pAto[kk].atomn;
          }
        }
      }
    }
  }*/

  /* check whether you need to include atoms Within a radius */
  if(iWithinFlag != 1)
  {
    selestr->selatm = calloc ( selestr->nselatm, sizeof(int));
    for ( ii=0; ii<selestr->nselatm; ii++)
      selestr->selatm[ii] = localsel[ii];
  }
  else if(iWithinFlag == 1)
  {
   // =================== //
   // *-* within code *-* //
   // =================== //
    iNumOfAtomPairs = selestr->nselatm * molecule->nato;

    piSele1 = walloc(iNumOfAtomPairs, sizeof(int));
    piSele2 = walloc(iNumOfAtomPairs, sizeof(int));

    for(ii=0; ii<selestr->nselatm; ii++)
    {
      for(jj=0; jj<molecule->nato; jj++)
      {
        iVectPos++;
        piSele1[iVectPos] = localsel[ii];
        piSele2[iVectPos] = jj + 1;
      }
    }

    piSeleAtoms = calloc(molecule->nato,  sizeof(int));
    pfDistances = DistanceVectCoor(&molecule->coor, piSele1, &molecule->coor, piSele2, iNumOfAtomPairs);
    
    for(ii=0; ii<iNumOfAtomPairs; ii++)
    {
      if(iWithinMode == 0)
      {
        if(pfDistances[ii] <= fShellLength)
          piSeleAtoms[piSele2[ii]-1]++;
      }
      
      else
      {
        if(pfDistances[ii] >= fShellLength)
          piSeleAtoms[piSele2[ii]-1]++;
      }
    }
    
    if(iWithinMode == 0)
      for(ii=0; ii<selestr->nselatm; ii++)
        piSeleAtoms[localsel[ii]-1]++;
    
    for(ii=0; ii<molecule->nato; ii++)
      if(piSeleAtoms[ii] > 0)
        iNumOfSeleAtom++;
        
    selestr->nselatm = iNumOfSeleAtom;
    selestr->selatm = calloc(selestr->nselatm, sizeof(int));
    
    iVectPos = -1;
    for(ii=0; ii<molecule->nato; ii++)
    {
      if(piSeleAtoms[ii] > 0)
      {
        iVectPos++;
        selestr->selatm[iVectPos] = ii + 1;
      }
    }
    
    free(piSeleAtoms);
    free(pfDistances);
    free(piSele1);
    free(piSele2);
  }

  free(localsel);
  
  return(0);
}
// ------------------------------------------------------------------
Res * GetSeleRes ( char *selestring, Selection *selestr, Molecule *molecule, int *nres )
{
  Res    *seleres;
  
  // waiting for a real need ...
  // restype
  //switch( restype )
  //{
    //case 'A':
      ////dothis();
      //break;
    //default:
      //break;
  //}
  seleres = calloc( nres[0], sizeof(Res));
  
  return seleres;
}
// ------------------------------------------------------------------
void GetSeleCoor2Vecs (  CoorSet *source, float *xcoor, float *ycoor, float *zcoor, Selection *sele )
{
   int ii;
   
   for ( ii=0; ii<sele->nselatm; ii++)
    {
     xcoor[ii] = source->xcoor[sele->selatm[ii]-1];
     ycoor[ii] = source->ycoor[sele->selatm[ii]-1];
     zcoor[ii] = source->zcoor[sele->selatm[ii]-1];
    }

   return;
}
// ------------------------------------------------------------------
void GetSeleCoor ( CoorSet *source, CoorSet *target, Selection *sele )
{
   int ii;
   
   if( target->nato != sele->nselatm )
   {
      fprintf( stderr, "target coorset not the right size!");
      exit(0);
   }
   for ( ii=0; ii<sele->nselatm; ii++)
   {
     target->xcoor[ii] = source->xcoor[sele->selatm[ii]-1];
     target->ycoor[ii] = source->ycoor[sele->selatm[ii]-1];
     target->zcoor[ii] = source->zcoor[sele->selatm[ii]-1];
   }
   
   return;
}
// ------------------------------------------------------------------
void CopyCoor ( CoorSet *source, CoorSet *target )
{
   int ii;
   
   if( target->nato != source->nato )
   {
      fprintf( stderr, "target coorset not the right size!");
      exit(0);
   }
   for ( ii=0; ii<source->nato; ii++)
   {
     target->xcoor[ii] = source->xcoor[ii];
     target->ycoor[ii] = source->ycoor[ii];
     target->zcoor[ii] = source->zcoor[ii];
   }
   
   return;
}
// ------------------------------------------------------------------
void HdrCp ( Trjh  *hdr1, Trjh  *hdr2 )
{
  int ii, jj;
  
  for ( ii=0; ii<5; ii++ )
   hdr2->type[ii] = hdr1->type[ii];
  
  for ( ii=0; ii<3; ii++ )
   for ( jj=0; jj<80; jj++ )
    hdr2->comm[ii][jj] = hdr1->comm[ii][jj];
  
  hdr2->nframe   = hdr1->nframe;
  hdr2->begstp   = hdr1->begstp;
  hdr2->skpstp   = hdr1->skpstp;
  hdr2->numstp   = hdr1->numstp;
  hdr2->tmstp    = hdr1->tmstp;
  hdr2->varpbc   = hdr1->varpbc;
  hdr2->comnum   = hdr1->comnum;
  hdr2->nato     = hdr1->nato;
  hdr2->fixatm   = hdr1->fixatm;
  hdr2->extra1   = hdr1->extra1;
  hdr2->extra2   = hdr1->extra2;
  hdr2->extra3   = hdr1->extra3;
  hdr2->c_nframe = hdr1->c_nframe;
  hdr2->c_numstp = hdr1->c_numstp;

  hdr2->sex       = 0;
  hdr2->sixtyfour = 0;
  hdr2->size      = 0;
  
  return;
}
// ------------------------------------------------------------------
RadList * ReadRadius(char *filename)
{
   int          ii;
   int          nlines;
   FILE        *radfile;
   char         buffer[128];
   RadList     *radlist;
   
   radfile = O_File( filename, "r");
   nlines = 0;
   while ( !feof(radfile))
   {
    if(fgets(buffer, 90, radfile)==NULL && ( !feof(radfile) || ferror(radfile) ))
    {
       fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
     ++nlines;
   }
   rewind(radfile);
   radlist = (RadList *) calloc( 1, sizeof(RadList));
   radlist->nradii = nlines;
   radlist->resname = calloc( nlines, sizeof(char *));
   for( ii=0; ii<nlines; ii++)
     radlist->resname[ii] = calloc( 4, sizeof(char));
   radlist->atomname = calloc( nlines, sizeof(char *));
   for( ii=0; ii<nlines; ii++)
     radlist->atomname[ii] = calloc( 4, sizeof(char));
   radlist->radius = calloc( nlines, sizeof(float));
   
   for( ii=0; ii<nlines; ii++)
   {
     if(fgets(buffer, 90, radfile)==NULL && ( !feof(radfile) || ferror(radfile) ))
     {
       fprintf(stderr, "Warning! Premature end of file reached!\n");
     }
     if( buffer[0] != '*' && buffer[0] != '#' )
     {
       sscanf( buffer, "%s %s %f", radlist->atomname[ii], radlist->resname[ii], &radlist->radius[ii]);
     }
   }
   
   return radlist;
}
// ------------------------------------------------------------------
int RadAssign( RadList *radlist, Molecule *molecule, double *pdAtomsRadii)
{
   int          ii, jj, found;
   
   for( ii=0; ii<molecule->nato; ii++ )
   {
     found = 0;
     for( jj=0; jj<radlist->nradii; jj++ )
     {
       if( !strcmp(molecule->rawmol.restype[ii], radlist->resname[jj]))
       {
         if( !strcmp(molecule->rawmol.atmtype[ii], radlist->atomname[jj]))
         {
           pdAtomsRadii[ii] = (double)radlist->radius[jj];
           found = 1;
           break;
         } 
       }
     }
     if( !found)
       fprintf( stderr, "Did not find a radius for %s %s\n", molecule->rawmol.restype[ii], molecule->rawmol.atmtype[ii]);
   }

   return 1;
}
// ------------------------------------------------------------------
void * walloc( int nblock, int blocksize )
{
   void *pointer;
   
//  #ifdef MACOSX
//   pointer = malloc( nblock*blocksize );
//  #else
// pointer = memalign( 16, nblock*blocksize);
   pointer = calloc( nblock, blocksize);
//  #endif
   
   return pointer;
}
// ------------------------------------------------------------------
void ** wrd_matAlloc( int rows, int columns, int blocksize )
{
   void   **pointer;
   int     ii;
   
   pointer    = calloc( rows, blocksize);
   pointer[0] = calloc( rows*columns, blocksize);
   for( ii=0; ii<rows; ii++)
     pointer[ii] = pointer[0] + ii*columns;
   
   return pointer;
}
// ------------------------------------------------------------------
float ** wrd_fmatAlloc( int rows, int columns )
{
   float   **pointer;
   int       ii;
   
   pointer    = calloc( rows,         sizeof(float *));
   pointer[0] = calloc( rows*columns, sizeof( float ));
   for( ii=0; ii<rows; ii++)
     pointer[ii] = pointer[0] + ii*columns;
   
   return pointer;
}
// ------------------------------------------------------------------
void  wrd_fmatFree( float **pointer )
{
   free(pointer[0]);
   free(pointer);
   
   return;
}
// ------------------------------------------------------------------
void CalloCoor ( CoorSet *crd, int nato)
{
   crd->cords = calloc( 3*nato, sizeof(float));
   crd->coors = calloc( 3, sizeof(float *));
   crd->xcoor = crd->coors[0] = &crd->cords[0];
   crd->ycoor = crd->coors[1] = &crd->cords[nato];
   crd->zcoor = crd->coors[2] = &crd->cords[2*nato];

   crd->nato = nato;
   crd->pbc_flag = 0;
   crd->pbc = NULL;

   return;
}
// ------------------------------------------------------------------
CoorSet * InitCoor ( int nato )
{
   CoorSet *crd;
   
   crd = calloc( 1, sizeof(CoorSet));
   crd->nato = nato;
   CalloCoor( crd, nato );
   
   return crd;
}
// ---------------------------------------------------------------------
void PrintCoor( CoorSet *crd, int start, int end )
{
  int     ii;
  
  if( end == 0 )
    end = crd->nato;
  
  for( ii=start; ii<end; ii++ )
  {
    fprintf( stdout, "%4d %10.5f%10.5f%10.5f\n", ii+1, crd->xcoor[ii], crd->ycoor[ii], crd->zcoor[ii] );
  }
  
  return;
}
// ---------------------------------------------------------------------
void PrintSeleCoor( CoorSet *crd, Selection *sele )
{
  int     ii;
  
  for( ii=0; ii<sele->nselatm; ii++ )
    fprintf( stdout, "%4d %10.5f%10.5f%10.5f\n", ii+1, crd->xcoor[sele->selatm[ii]-1], crd->ycoor[sele->selatm[ii]-1], crd->zcoor[sele->selatm[ii]-1] );
  
  return;
}
// ---------------------------------------------------------------------
void CleanCoor ( CoorSet *crd)
{
   free( crd->cords );
   free( crd->coors );
   return;
}
// ---------------------------------------------------------------------
void DelCoor ( CoorSet *crd)
{
   CleanCoor( crd );
   free(crd);
}
// ---------------------------------------------------------------------
int copyvar( char *origin, int start, int end, char *copy , int csize)
{
   int           ii, jj;
   
   memset( copy, '\0', csize+1 );
   jj=0;
   for( ii=start; ii<end; ii++ )
   {
    copy[jj] = origin[ii];
    jj++;
   }
   return 1;
}
// ---------------------------------------------------------------------
char * expandrange( char *intstring, int *len )
{
   FILE    *beforefile, *afterfile;
   int      ii;
   char    *buffer, *newintstring;
   size_t   size;
   char     tmpchar, midchar;
   int      resnum, rangestart, rangeflag, tot_int=0;
   
   #ifdef MACOSX
    fprintf( stderr, "Error! Range-spanning selections not implemented for MACOSX\n");
    exit(47);
   #else
   
    beforefile = fmemopen ( intstring, 10240, "r");
    afterfile  = open_memstream ( &buffer, &size);
    
    //if ( (intstring[0]!='?' && intstring[0]!='*' && intstring[0]!='+' && intstring[0]!='@' && intstring[0]!='!' ) || intstring[1]!='(')
    //{
     //printf("Error! selection (%s) is not valid\n", intstring);
     //exit(0);
    //}
    if ( intstring[0]=='?' || intstring[0]=='*' || intstring[0]=='+' || intstring[0]=='@' || intstring[0]=='!'  )
    {
      // copy the starting "@(" string
      fscanf(beforefile, "%c", &tmpchar);
      fprintf(afterfile, "%c", tmpchar);
      fscanf(beforefile, "%c", &tmpchar);
      fprintf(afterfile, "%c", tmpchar);
    }
    
    midchar = 'x';
    rangeflag = 0;

    while( midchar != '\0' &&  midchar != ')')
    {
     fscanf( beforefile, "%d", &resnum  );
     fscanf( beforefile, "%c", &midchar );

     if ( rangeflag == 0)
     {
       if ( midchar == '|' || midchar == ')' || midchar == '\0' || midchar == ',' || midchar == ';' )
       {
         if( midchar == ',' || midchar == ';' )
           midchar = '|';
         fprintf( afterfile, "%d%c", resnum, midchar);
         tot_int += 1;
         rangeflag=0;
       }
       if ( midchar == '-' )
       {
         rangestart = resnum;
         rangeflag = 1;
       }
     }
     else if (rangeflag == 1)
     {
       if ( midchar == '|' || midchar == ')' || midchar == '\0' || midchar == ',' || midchar == ';' )
       {
         if( midchar == ',' || midchar == ';' )
           midchar = '|';
         if ( resnum < rangestart )
         {
           fprintf( stderr, "Error ! wrong range: %d is lower than %d\n", resnum, rangestart);
           exit(0);
         }
         for( ii=rangestart; ii<resnum; ii++ )
         {
           fprintf( afterfile, "%d|", ii);
           tot_int += 1;
         }
         fprintf( afterfile, "%d%c", resnum, midchar );
         tot_int += 1;
         rangeflag = 0;
       }
       if ( midchar == '-' )
       {
         fprintf(stderr, "Error! check sintax: double \"-\" not allowed\n");
         exit(0);
       }
     }
    }
    
    fclose (beforefile);
    fclose (afterfile);
    
    newintstring = calloc( strlen(buffer)+1, sizeof(char));
    memset( newintstring, '\0', strlen(buffer)+1);
    
    len[0] = tot_int;
    sprintf( newintstring, "%s", buffer );
    //printf( "DEBUG: %s\n", buffer );
    return newintstring;
   #endif
}
// ---------------------------------------------------------------------
int * readintstring( char *intstring, int tot_int )
{
   int  *intlist, ii;
   char delims[] = "|";
   char *result = NULL;
   
   if( tot_int == 0)
   {
     ii = 0;
     tot_int = 1;
     while( intstring[ii] !='\0')
     {
       if( intstring[ii] == '|' || intstring[ii] == ',')
         tot_int += 1;
       ii ++;
     }
   }
   
   intlist = calloc( tot_int, sizeof(int));
   
   ii = 0;
   result = strtok( intstring, delims );
   while( result != NULL )
   {
     sscanf( result, "%d", &intlist[ii] );
     ii++;
     result = strtok( NULL, delims );
   }
   
   return intlist;
}
//======================================================================
// -------------- GROMACS SECTION --------------------------------------
//======================================================================
// ------------------------------------------------------------------
void ReadXtcHeaderLITE ( Traj *trj )
  /* This function reads an XTC file and updates the header.
     Only the number of atoms, first step, skipping between steps 
     first time and dt are calculated. 
     
     A note to the programer:
     Use ReadXTC_HeaderFULL when more info is needed
   */
{

    char          *func="ReadXtcHeaderLITE";
    int           natom, step,first_step,skpstp;
    float         time,first_time,dt;
    rvec          *x;
    float         precision;
    int           result;
    matrix        box;

    xdrfile_close(trj->xtc);
    read_xtc_natoms(trj->filename,&natom);
    x = calloc(natom,sizeof(*x));
    trj->xtc=xdrfile_open(trj->filename,"r");
    result=read_xtc(trj->xtc,natom,&step,&time,box,x,&precision);
    if (exdrOK != result) {
      printf ("Wordom: %s\n\tError! could not read the first XTC coordinate\n",func);
      rfree(x);
      exit(99);
    }
    trj->hdr->begstp=step;
    trj->hdr->nato=natom;   
    first_step=step;
    first_time=time;
    // read next frame to calculate dt 
    result=read_xtc(trj->xtc,natom,&step,&time,box,x,&precision);
    if (exdrENDOFFILE == result) {
      // There's just one frame
      rfree(x);
      trj->hdr->nframe=1;
      trj->hdr->numstp=1;
      trj->hdr->tmstp=0;
      trj->hdr->skpstp=0;
      trj->hdr->varpbc=0;
      xdrfile_rewind(trj->xtc,trj->filename,"r");
      return;
    }
    skpstp=step-first_step;
    trj->hdr->skpstp=skpstp;
    dt=(time-first_time)/skpstp; // ps, between steps
    trj->hdr->tmstp=dt/AKMA2PS;
    rfree(x);
    xdrfile_rewind(trj->xtc,trj->filename,"r");
    return;
}
// ------------------------------------------------------------------
void ReadXtcHeaderFULL ( Traj *trj )
{
  /* This function reads an XTC file and updates the header,
     to include data on the number of frames and steps as well.
     It is needed due to difficulties with access to the termini of XTC files.
     
     A note to the programer:
     Since this function goes over the whole XTC it is relatively slow.
     Use ReadXTC_HeaderLITE instead whenever possible.
   */

    char          *func="ReadXtcHeaderFULL";
    int           natom, step,first_step,skpstp,numstp,nframe=0;
    float         time,first_time,dt;
    matrix        box,pbox;
    rvec          *x;
    float         precision;
    float         round_err;
    int           result;
    //bool          bVbox=0;
    int           bVbox=0;

    xdrfile_close(trj->xtc);
    read_xtc_natoms(trj->filename,&natom);
    x = calloc(natom,sizeof(*x));
    trj->xtc=xdrfile_open(trj->filename,"r");
    result=read_xtc(trj->xtc,natom,&step,&time,box,x,&precision);
    if (exdrOK != result) {
      printf ("Wordom: %s\n\tError! could not read the first XTC coordinate\n",func);
      rfree(x);
      exit(99);
    }
    round_err=1/precision;
    copymat(box,pbox);
    trj->hdr->begstp=step;
    trj->hdr->nato=natom;
    trj->hdr->fixatm=-1;
    trj->hdr->sixtyfour=-1;
    trj->hdr->c_nframe=-1;
    first_step=step;
    first_time=time;
    nframe++;
    // read next frame to calculate dt 
    result=read_xtc(trj->xtc,natom,&step,&time,box,x,&precision);
    if (exdrENDOFFILE == result) {
      // There's just one frame
      rfree(x);
      trj->hdr->nframe=1;
      trj->hdr->numstp=1;
      trj->hdr->tmstp=0;
      trj->hdr->skpstp=0;
      trj->hdr->varpbc=0;
      xdrfile_rewind(trj->xtc,trj->filename,"r");
      return;
    }
    skpstp=step-first_step;
    trj->hdr->skpstp=skpstp;
    dt=(time-first_time)/skpstp; // ps, between steps
    trj->hdr->tmstp=dt/AKMA2PS;
    bVbox= eqmat(box,pbox,round_err) ? 0 : 1;
    copymat(box,pbox);
    nframe++;
    while (exdrOK==result) {
      result=read_xtc(trj->xtc,natom,&step,&time,box,x,&precision);
      if (!bVbox) {
	bVbox=eqmat(box,pbox,round_err) ? 0 : 1;
	copymat(box,pbox);
      }
      nframe++;     
    }
    numstp=(time-first_time)/dt; // not counting the first step!
    numstp+=skpstp;
    nframe--;
    trj->hdr->nframe=nframe;
    trj->hdr->numstp=numstp;
    trj->hdr->varpbc=bVbox;
    rfree(x);
    xdrfile_rewind(trj->xtc,trj->filename,"r");
    return;
}
// ------------------------------------------------------------------
void ReadXtcCoor  ( Traj *trj, CoorSet *crd, int frn)
{
   int           natom,ii,step;
   float         time;
   matrix        box;
   rvec          *x;
   float         precision;
   int           result=exdrOK;

   /* unfortunately, no random access to xtc at the moment */
   natom = trj->hdr->nato;
   x = calloc(natom,sizeof(*x));
   if (trj->previousframe > frn-1)
   {
     xdrfile_rewind(trj->xtc,trj->filename,"r");
   }
   while ( trj->previousframe < frn && exdrOK == result)
   {
     result = read_xtc( trj->xtc, natom, &step, &time, box, x, &precision );
     trj->previousframe ++;
   }
   if (exdrOK != result)
   {
     fprintf(stderr,"\nWARNING: could not read frame %d from xtc file\n", frn);
     exit(0);
   }

   for ( ii=0; ii<natom; ii++)
   {
     crd->xcoor[ii]=x[ii][XX]/0.1;
     crd->ycoor[ii]=x[ii][YY]/0.1;
     crd->zcoor[ii]=x[ii][ZZ]/0.1;
   }

   if (x != NULL)
     rfree(x);

   if ( verbose )
     fprintf( stdout, "%8.3f%8.3f%8.3f\n%8.3f%8.3f%8.3f\n%8.3f%8.3f%8.3f\n", box[0][0], box[0][1], box[0][2], box[1][0], box[1][1], box[1][2], box[2][0], box[2][1], box[2][2]);
   
   if( crd->pbc_flag == 1 )
   {
     if( box[1][0] == 0 && box[2][0] == 0 && box[0][1] ==0 && box[0][2] == 0)
     {
       /* cubic box */
       crd->pbc->a_size = box[0][0]/0.1;
       crd->pbc->b_size = box[1][1]/0.1;
       crd->pbc->c_size = box[2][2]/0.1;
       crd->pbc->angle1 = 0.0;
       crd->pbc->angle2 = 0.0;
       crd->pbc->angle3 = 0.0;
     }
     else
     {
       fprintf( stderr, "non-cubic PBC found in xtc file: wordom can't actually handle those\n");
       crd->pbc->a_size = 0.0;
       crd->pbc->b_size = 0.0;
       crd->pbc->c_size = 0.0;
       crd->pbc->angle1 = 0.0;
       crd->pbc->angle2 = 0.0;
       crd->pbc->angle3 = 0.0;
     }
   }
   return;
}
// ------------------------------------------------------------------
void WriteXtcCoor (XDRFILE* xtc, CoorSet *crd, Trjh *hdr )
{

    char *func="WriteXTC_Coor";
    char *warn="Warning: - box matrix is not written\n";
    int ii,natom,step=0;
    float precision = 10000000; //precision=10000;

    float time=0;
    rvec *x;
    matrix box;
    
    natom=hdr->nato;    
    x=malloc(sizeof(rvec)*natom);
    if(x==NULL) {
      fprintf(stderr, "Could not allocate memory for %d atoms\n",natom);
      exit(99);
    }
//  clear_mat(box);
    box[XX][XX]=box[XX][YY]=box[XX][ZZ]=0.0;
    box[YY][XX]=box[YY][YY]=box[YY][ZZ]=0.0;
    box[ZZ][XX]=box[ZZ][YY]=box[ZZ][ZZ]=0.0;

    for ( ii=0; ii< natom; ii++) {
      x[ii][XX]=crd->xcoor[ii]*0.1;
      x[ii][YY]=crd->ycoor[ii]*0.1;
      x[ii][ZZ]=crd->zcoor[ii]*0.1;
    }

    write_xtc(xtc, natom, step, time, box, x, precision);
    printf ("%s: %s\n",func,warn);
    return;
}
// ------------------------------------------------------------------
bool xdrfile_rewind (XDRFILE *xtc, const char *filename, const char *mode)
{
    if (xtc!=NULL)
      xdrfile_close(xtc);
    xtc = xdrfile_open(filename, mode);
    if (xtc==NULL) {
      printf ("Could not reopen file %s\n",filename);
      return 0;
    }
    return 1;
}
// ------------------------------------------------------------------
void *ralloc (void *ptr,int nelem)
{
    char *funcname="ralloc";
    void *p;

    if (nelem<1) {
	fprintf(stderr, "%s: cannot allocate memory for %d elements\n",funcname, nelem);
	exit(99);
    }

    p=calloc((size_t) nelem,sizeof(*ptr));
    if (NULL==p) {
	fprintf(stderr, "Could not allocate memory in %s\n",funcname);
	exit(99);
    }
    return (p);
}
// ------------------------------------------------------------------
void rfree (void *ptr) 
{
    if (NULL!=ptr) 
      free(ptr);
}
// ------------------------------------------------------------------
void copymat ( matrix source,matrix dest)
{
    int i,j;
    for (i=XX;i<=ZZ;i++)
      for (j=XX;j<=ZZ;j++)
	dest[i][j]=source[i][j];
}
// ------------------------------------------------------------------
bool eqmat ( matrix m1,  matrix m2, float round_err)
{
    bool bEqual=1;
    int i,j;
    for (i=XX;i<=ZZ;i++)
      for (j=XX;j<=ZZ;j++) {
	if (abs(m1[i][j]-m2[i][j]) > round_err)
	  bEqual=0;
      }
    
    return bEqual;
}
// ------------------------------------------------------------------
void ShowXtcHeader( Trjh *hdr )
{
   float tmstp;
   // convert timestep from CHARMM units in fs
   // CHARMM units originate from the initial development
   // of WORDOM
   tmstp = floor ( (hdr -> tmstp * 488.8780249) + 0.5)/10;
   printf ( "XTC info\n");
   printf ( "Number frames\t%-20d\n",   hdr->nframe    );
   printf ( "Begin  step\t%-20d\t",     hdr->begstp    );
   printf ( "Skip   step\t\t%d\n",      hdr->skpstp    );
   printf ( "Number steps \t%-20d\t",   hdr->numstp    );
   printf ( "Time   step (fs)\t%2.1f\n",     tmstp     );
   printf ( "Number atoms\t%-20d\t",    hdr->nato      );
   printf ( "Variable PBC\t%d\t\t\n",   hdr->varpbc    );
   printf ( "END\n\n" );
   
   return;
}        
// ------------------------------------------------------------------
void AllocMol( int nato , Molecule *molecule )
{
  int ii;
  
  molecule->coor.cords = walloc ( 3*molecule->nato, sizeof(float));
  molecule->coor.xcoor = &molecule->coor.cords[0];
  molecule->coor.ycoor = &molecule->coor.cords[molecule->nato];
  molecule->coor.zcoor = &molecule->coor.cords[2*molecule->nato];
  molecule->coor.coors    =  calloc( 3, sizeof(float *));
  molecule->coor.coors[0] =  molecule->coor.cords;
  molecule->coor.coors[1] = &molecule->coor.cords[molecule->nato];
  molecule->coor.coors[2] = &molecule->coor.cords[2*molecule->nato];
  
  if( molecule->altLocB )
    CalloCoor( &molecule->coor2, molecule->nato );
  if( molecule->altLocC )
    CalloCoor( &molecule->coor3, molecule->nato );
  
  molecule->rawmol.hetatm  = calloc ( molecule->nato, sizeof(short) );
  molecule->rawmol.atomn   = calloc ( molecule->nato, sizeof( int ) );
  molecule->rawmol.chainId = calloc ( molecule->nato, sizeof( char) );
  molecule->rawmol.altLoc  = calloc ( molecule->nato, sizeof( char) );
  molecule->rawmol.iCode   = calloc ( molecule->nato, sizeof( char) );
//  molecule->rawmol.ipresn  = calloc ( molecule->nato, sizeof( int ) );
  molecule->rawmol.presn   = calloc ( molecule->nato, sizeof( int ) );
  molecule->rawmol.resn    = calloc ( molecule->nato, sizeof( int ) );
  molecule->rawmol.occup   = calloc ( molecule->nato, sizeof(float) );
  molecule->rawmol.bFac    = calloc ( molecule->nato, sizeof(float) );
  molecule->rawmol.segend  = calloc ( molecule->nato, sizeof(short) );
  molecule->rawmol.segbeg  = calloc ( molecule->nato, sizeof(short) );
  
  // this two just for crd-compatibility
//  molecule->rawmol.cpresn  = calloc ( molecule->nato, sizeof( int ) );
  molecule->rawmol.weight  = calloc ( molecule->nato, sizeof(float) );

  molecule->rawmol.atmtype    = calloc ( molecule->nato, sizeof(char *));
  for ( ii=0; ii < molecule->nato ; ii++)
   molecule->rawmol.atmtype[ii] = calloc( 4, sizeof(char));

  molecule->rawmol.atmId      = calloc ( molecule->nato, sizeof(char *));
  for ( ii=0; ii < molecule->nato ; ii++)
   molecule->rawmol.atmId[ii]   = calloc( 5, sizeof(char));

  molecule->rawmol.restype    = calloc ( molecule->nato, sizeof(char *));
  for ( ii=0; ii < molecule->nato ; ii++)
   molecule->rawmol.restype[ii] = calloc( 3, sizeof(char));

  molecule->rawmol.segId      = calloc ( molecule->nato, sizeof(char *));
  for ( ii=0; ii < molecule->nato ; ii++)
   molecule->rawmol.segId[ii]   = calloc( 4, sizeof(char));

  molecule->rawmol.element    = calloc ( molecule->nato, sizeof(char *));
  for ( ii=0; ii < molecule->nato ; ii++)
   molecule->rawmol.element[ii] = calloc( 2, sizeof(char));

  molecule->rawmol.charge     = calloc ( molecule->nato, sizeof(char *));
  for ( ii=0; ii < molecule->nato ; ii++)
   molecule->rawmol.charge[ii] = calloc( 2, sizeof(char));

/*
   molecule->coor.cords =  calloc ( 3*nato, sizeof(float));
   molecule->coor.xcoor = &molecule->coor.cords[0];
   molecule->coor.ycoor = &molecule->coor.cords[nato];
   molecule->coor.zcoor = &molecule->coor.cords[2*nato];
   molecule->coor.coors    =  calloc( 3, sizeof(float *));
   molecule->coor.coors[0] =  molecule->coor.cords;
   molecule->coor.coors[1] = &molecule->coor.cords[nato];
   molecule->coor.coors[2] = &molecule->coor.cords[2*nato];

   molecule->rawmol.hetatm  = calloc ( nato, sizeof(short) );
   molecule->rawmol.atomn   = calloc ( nato, sizeof( int ) );
   molecule->rawmol.chainId = calloc ( nato, sizeof( char) );
   molecule->rawmol.altLoc  = calloc ( nato, sizeof( char) );
   molecule->rawmol.iCode   = calloc ( nato, sizeof( char) );
   molecule->rawmol.resn    = calloc ( nato, sizeof( int ) );
   molecule->rawmol.occup   = calloc ( nato, sizeof(float) );
   molecule->rawmol.bFac    = calloc ( nato, sizeof(float) );

   // this two just for crd-compatibility
   molecule->rawmol.presn   = calloc ( nato, sizeof( int ) );
   molecule->rawmol.weight  = calloc ( nato, sizeof(float) );

   molecule->rawmol.atmtype    = calloc ( nato, sizeof(char *));
   for ( ii=0; ii < nato ; ii++)
    molecule->rawmol.atmtype[ii] = calloc( 4, sizeof(char));

   molecule->rawmol.atmId      = calloc ( nato, sizeof(char *));
   for ( ii=0; ii < nato ; ii++)
    molecule->rawmol.atmId[ii]   = calloc( 5, sizeof(char));

   molecule->rawmol.restype    = calloc ( nato, sizeof(char *));
   for ( ii=0; ii < nato ; ii++)
    molecule->rawmol.restype[ii] = calloc( 3, sizeof(char));

   molecule->rawmol.segId      = calloc ( nato, sizeof(char *));
   for ( ii=0; ii < nato ; ii++)
    molecule->rawmol.segId[ii]   = calloc( 4, sizeof(char));

   molecule->rawmol.element    = calloc ( nato, sizeof(char *));
   for ( ii=0; ii < nato ; ii++)
    molecule->rawmol.element[ii] = calloc( 2, sizeof(char));

   molecule->rawmol.charge     = calloc ( nato, sizeof(char *));
   for ( ii=0; ii < nato ; ii++)
    molecule->rawmol.charge[ii] = calloc( 2, sizeof(char));
*/   
   return;
}
// ---------------------------------------------------------------------
Molecule * MakeSeleMol ( Selection *selstr , Molecule *inmol ) //, Molecule **selmol )
{
  int ii;
  Molecule *selemol;

  selemol = calloc (1, sizeof(Molecule));
  selemol->nato = selstr->nselatm ;
  AllocMol(selstr->nselatm,selemol);
  
  /* fill coorset */
  for( ii=0; ii<selstr->nselatm; ii++)
  {
    selemol->coor.xcoor[ii] = inmol->coor.xcoor[selstr->selatm[ii]-1];
    selemol->coor.ycoor[ii] = inmol->coor.ycoor[selstr->selatm[ii]-1];
    selemol->coor.zcoor[ii] = inmol->coor.zcoor[selstr->selatm[ii]-1];
  }
  selemol->rawmol.name = calloc( strlen(inmol->rawmol.name)+1, sizeof(char));
  strcpy( selemol->rawmol.name, inmol->rawmol.name );
  for( ii=0; ii<selstr->nselatm; ii++)
  {
    selemol->rawmol.hetatm[ii]   = inmol->rawmol.hetatm[selstr->selatm[ii]-1];
    //selemol->rawmol.atomn[ii]   = inmol->rawmol.atomn[selstr->selatm[ii]-1];
    selemol->rawmol.atomn[ii]   = ii+1;
    selemol->rawmol.chainId[ii] = inmol->rawmol.chainId[selstr->selatm[ii]-1];
    selemol->rawmol.altLoc[ii]  = inmol->rawmol.altLoc[selstr->selatm[ii]-1];
    selemol->rawmol.iCode[ii]   = inmol->rawmol.iCode[selstr->selatm[ii]-1];
    selemol->rawmol.resn[ii]    = inmol->rawmol.resn[selstr->selatm[ii]-1];
    selemol->rawmol.occup[ii]   = inmol->rawmol.occup[selstr->selatm[ii]-1];
    selemol->rawmol.bFac[ii]    = inmol->rawmol.bFac[selstr->selatm[ii]-1];
    selemol->rawmol.presn[ii]   = inmol->rawmol.presn[selstr->selatm[ii]-1];
    selemol->rawmol.weight[ii]  = inmol->rawmol.weight[selstr->selatm[ii]-1];
    StringCp ( inmol->rawmol.atmtype[selstr->selatm[ii]-1] , selemol->rawmol.atmtype[ii], 4);
    StringCp ( inmol->rawmol.atmId[selstr->selatm[ii]-1] , selemol->rawmol.atmId[ii], 5);
    StringCp ( inmol->rawmol.restype[selstr->selatm[ii]-1] , selemol->rawmol.restype[ii], 3);
    StringCp ( inmol->rawmol.segId[selstr->selatm[ii]-1] , selemol->rawmol.segId[ii], 4);
    StringCp ( inmol->rawmol.element[selstr->selatm[ii]-1] , selemol->rawmol.element[ii], 2);
    StringCp ( inmol->rawmol.charge[selstr->selatm[ii]-1] , selemol->rawmol.charge[ii], 2);
  }
  /*
  fprintf( stdout, "DEBUG : called makeselemol !!!\n");
  printf( "=== DEBUG: list right before Raw2Struct (%s) ===\n", selemol->rawmol.name );
  for( ii=1 ; ii<selemol->nato; ii++ )
  {
    printf( "DEBUG: atom %d has presn %d\n", ii+1, selemol->rawmol.presn[ii] );
  }
  printf( "==============end list before (%s)===============\n", selemol->rawmol.name );
  */
  Raw2Struct( selemol );
    //fprintf( stdout, "Checking selection (in MKmol)\n");
    //for( ii=0; ii<selstr->nselatm; ii++)
      //fprintf( stdout, "%d\n", selstr->selatm[ii] );
  
  return selemol;
}
