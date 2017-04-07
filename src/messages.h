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
// ===================== FUNCTIONS ===================================
/*! \file messages.h
 \brief messages printed by wordom
 
 Misc messages (help, about) printed by wordom in response to
 command line calls.
*/
void Print_About	( )     
{
   printf  ( "\n" );
   printf  ( "\t---------------------   WORDOM   ---------------------\n");
//   printf  ( "\t---------------- dev snapshot (121211) ---------------\n");
   printf  ( "\t----------------  version 0.23 -  rc1  ---------------\n");
   printf  ( "\t      A command-line based program to manipulate\n"       );
   printf  ( "\t        and analyze molecular structural data\n");
   printf  ( "\t- Authors:   M.Seeber, A.Felline, F.Raimondi, R.Friedman,\n"  );
   printf  ( "\t-            S.Muff, M.Cecchini, F.Rao and G.Settanni \n"  );
   printf  ( "\t- Copyright: University of Modena and Reggio Emilia (IT) and\n");
   printf  ( "\t-            University of Zurich (CH) \n");
   printf  ( "\t- License:   GPL\n" 				);
  #ifndef LAPACK
   printf  ( "\t*** This executable compiled without lapack support ***\n" );
   printf  ( "\t***         Some features might be missing          ***\n" );
  #endif
   printf  ( "\n" );
}
// ------------------------------------------------------------------
void Print_Help	( struct sopt *OPT  )
{

   printf  ( "\n" );
   printf  ( " Usage:\t wordom [options]\n"         );
   printf  ( "\n" );
   printf  ( " Parameters:\n"     );
   printf  ( " -omol  FILE \tName of the output molecule file\n"     );
   printf  ( " -imol  FILE \tName of the input  molecule file\n"     );
   printf  ( " -otrj  FILE \tName of the output trajectory file\n"     );
   printf  ( " -itrjc  FILE \tName of the input  trajectory file\n"     );
   printf  ( " -beg   INT  \tBegin - first frame to consider\n"    );
   printf  ( " -end   INT  \tEnd   - last frame to consider\n"     );
   printf  ( " -skip  INT  \tNumber of steps between frames\n"  );
   
   printf  ( " ==============================================\n"                );
   printf  ( " Functions (list not complete):\n"      );
   printf  ( " -F     FILE \tName of Frame numbers list file\n"    );
   printf  ( " -f     INT  \tTrajectory frame number\n"            );
   printf  ( " -M     FILE \tMerge - Name of TRJ names list file\n"              );
   printf  ( " -m     NAME \tMerge - Basename\n"             );
   printf  ( " -S     FILE \tSum - Name of DCD names list file\n"              );
   printf  ( " -mono       \tMol2trj: writes a single mol to a trj\n"           );
   printf  ( " -oxyz  NAME \tTrj2Xyz: prints coordinates of selected atoms\n"           );
   printf  ( " -head  FILE \tPrint the trj header - trj name\n"            );
   printf  ( " -iA    FILE \tAnalysis Calculation\n"            );
   printf  ( "\n" );
   printf  ( " -about\t\tPrint version and authors\n"      );
   printf  ( " -h\t\tPrint this screen\n"  );
   printf  ( "\n" );
   printf  ( " For full documentation see:\n"      );
   printf  ( " htpp://wordom.sf.net/\n"      );
   printf  ( " and the User Guide\n"      );
   printf  ( "\n" );
   
   if ( OPT->HELP_FLAG )
    {
       printf("Detailed info about \"%s\" are not available,sorry.\n\n", OPT->HELP_ARG);
    }
}
// ------------------------------------------------------------------
