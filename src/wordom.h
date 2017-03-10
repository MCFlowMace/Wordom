// ===================== wordom.h ====================================
// ------------------------------------------------------------------
// - Author: 	Michele Seeber , Francesco Rao and Marco Cecchini
//   		University of Zurich
// - Date:	21/03/2003
// - License:	GPL
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

/*! \file wordom.h
 \brief command-line options structure
 
 wordom.h lists the command-line options structure, which is used
 by a number of other files in he wordom program
*/

#ifndef WORDOM
#define WORDOM

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <fnmatch.h>
#define  FNM_EXTMATCH (1 << 5)
#include <getopt.h>
#ifdef THREADED
 #include <pthread.h>
#endif
/*#ifdef MACOSX
 #include "/usr/include/malloc/malloc.h"
#else
 #include <malloc.h>
#endif*/

#include "fileio.h"

// ------------------------------------------------------------------
extern int errno;
// ===================== STRUCTURES ==================================
//! options to be passed to iacalc_read
/*!
 when analysis is run from the command line rather than through a script
 this structure gathers the options to be passed to a parser
*/
typedef struct iaOpt            // options to be passed to iacalc_read
{
   char        *option;
   char        *value;
   int          noptions;
   struct iaOpt      *next;
} _iaOpt;        
// ------------------------------------------------------------------
// ** System options
//! Command line options - used program-wide
/*!
 Structure with system-wide options. Mainly filled in the main with data
 from the command line. Flags point out action to take.
*/
typedef struct sopt
{
   int     DEBUG_FLAG;  /*!< \brief Debug flag */

   int       FRN_FLAG;  /*!< \brief Frame number flag */
   int          FRN_N;  /*!< \brief Frame number */

   int      HEAD_FLAG;  /*!< \brief Print header flag */

   int      ITRJ_FLAG;  /*!< \brief Input trajectory flag */
   int      OTRJ_FLAG;  /*!< \brief Output trajectory flag */
   int      ATRJ_FLAG;  /*!< \brief Appendurum (to be appended) trajectory flag */
   int      IMOL_FLAG;  /*!< \brief Input molecule flag */
   int      OMOL_FLAG;  /*!< \brief Output molecule flag */
   int      AMOL_FLAG;  /*!< \brief Appendurum (to be appended) molecule flag */
   int       FRL_FLAG;  /*!< \brief Frame File flag */
   int      T2Ms_FLAG;  /*!< \brief Trj to many Mols flag */
   int      MONO_FLAG;  /*!< \brief Single mol to trj flag */
   int     TCONV_FLAG;  /*!< \brief Trajectory conversion flag */
   int       BEG_FLAG;  /*!< \brief Generic BEGIN flag */
   int       END_FLAG;  /*!< \brief Generic END flag */
   int      TRJL_FLAG;  /*!< \brief Trajectory List File flag */
   int      TRJS_FLAG;  /*!< \brief Trajectory Sum flag */
   int      TRJC_FLAG;  /*!< \brief Trajectory Combine flag */
   int      SKIP_FLAG;  /*!< \brief Skip flag */
   int        IA_FLAG;  /*!< \brief Generic Analysis Flag */
   int        OA_FLAG;  /*!< \brief Generic Analysis Output File Flag */
   int        Ia_FLAG;  /*!< \brief Command-line Analysis Flag */
   int        IE_FLAG;  /*!< \brief Generic non-trj Analysis Flag */
   int        Ie_FLAG;  /*!< \brief Generic non-trj Analysis Flag */
   int      CONV_FLAG;  /*!< \brief Conversion flag */
   int       PBC_FLAG;  /*!< \brief Periodic Boundary Condition Flag */
   int     NOPBC_FLAG;  /*!< \brief NEGLECT Periodic Boundary Condition Flag */
   int    ATMNUM_FLAG;  /*!< \brief Xyz atom number list Flag */
   int      HELP_FLAG;  /*!< \brief Man Help Flag */
   int   VERBOSE_FLAG;  /*!< \brief Verbose Flag */
   int       MOD_FLAG;  /*!< \brief DCD mod Flag */
   int       THR_FLAG;  /*!< \brief Threading Flag */
   int       AVG_FLAG;  /*!< \brief Average Extraction Flag */
   int   MERGEAN_FLAG;  /*!< \brief Merge files for analysis */
   int      SELE_FLAG;  /*!< \brief Selection flag */
   int     CSELE_FLAG;  /*!< \brief Check Selection flag */
   int      INFO_FLAG;  /*!< \brief Generic Info flag */

   char	   *ITRJ_FILE;  /*!< \brief Input trj filename */
   char	   *OTRJ_FILE;  /*!< \brief Input trj filename */
   char    *ATRJ_FILE;  /*!< \brief Addendurum trj filename */
   char    *IMOL_FILE;  /*!< \brief Input  molecule filename */
   char    *OMOL_FILE;  /*!< \brief Output molecule filename */
   char    *AMOL_FILE;  /*!< \brief Addendurum  molecule file name */
   char	   *TRJL_FILE;  /*!< \brief Trj List Filename */
   char	   *TRJS_FILE;  /*!< \brief Trj List (for Trj Sum) file name */
   char	   *TRJC_FILE;  /*!< \brief Trj List (for Trj Combine) file name */
   char	    *FRL_FILE;  /*!< \brief Frame List Filename */
   char	   *T2Ms_FILE;  /*!< \brief Frame List for Trj-to-many-Mols File Name */
   char	     *IA_FILE;  /*!< \brief Generic Analysis Input File Name */
   char	   *Ia_STRING;  /*!< \brief Generic Analysis Command Line Selection */
   char	     *OA_FILE;  /*!< \brief Generic Analysis Output File Name */
   char  *ATMNUM_FILE;  /*!< \brief Extraction Atom Number File Name */
   char	    *HELP_ARG;  /*!< \brief Argument of Help topic */
   char   *MOD_STRING;  /*!< \brief Field=Value for Head_Mod */
   char  *SELE_STRING;  /*!< \brief Selection String */
   char *MERGEAN_FILE;  /*!< \brief Merge files for analysis */
   char	   *DCDN_BASE;  /*!< \brief Trajectory Basename */
   char *SUBSELE_STRING;/*!< \brief subselection string */
   
   moltype  IMOL_TYPE;  /*!< \brief Molecule input  file type */
   moltype  OMOL_TYPE;  /*!< \brief Molecule output file type */
   moltype  AMOL_TYPE;  /*!< \brief Molecule append file type */
   trjtype  ITRJ_TYPE;  /*!< \brief Trajectory input  file type */
   trjtype  OTRJ_TYPE;  /*!< \brief Trajectory output file type */
   trjtype  ATRJ_TYPE;  /*!< \brief Trajectory append file type */
   
   int	     TRJN_BEG;  /*!< \brief First Trajectory Number */
   int	     TRJN_END;  /*!< \brief Last  Trajectory Number */

   int	     TYPE_NUM;  /*!< \brief Number of Atom Types to Filter */
   char     *TYPE_ATM;  /*!< \brief Atom Type Array */

   int	    SKIP_STEP;  /*!< \brief Skip Step */
   int	     BEG_STEP;  /*!< \brief First step */
   int	     END_STEP;  /*!< \brief Last Step  */
   
   int	     THR_NUM;   /*!< \brief Number of Required Threads  */
   
   float    PBCBOX[3]; 
   
   _iaOpt       *iaopt, *current;

} sopt;
void setOPT(sopt *);
// ------------------------------------------------------------------

#endif
