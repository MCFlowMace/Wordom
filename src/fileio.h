/* ------------------------------------------------------------------ */
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
/* ------------------------------------------------------------------ */

/*! \file fileio.h
 \brief the input/output functions and structures
 
 fileio provides the low-level input/output functions
 to acces data in structure (crd/pdb) and trajectory(dcd/xtc) files
*/

#ifndef FILEIO
#define FILEIO

// ===   GROMACS section   ===

#ifndef XTC_MAGIC
 #define XTC_MAGIC 1995
#endif

#define DEG2RAD           (M_PI/180.0) // from GMX physics.h
#define AKMA2FS           48.88780249   
#define AKMA2PS         0.04888780249  //
#define XX 0
#define YY 1
#define ZZ 2

#include "xdrfile.h"
#include "xdrfile_xtc.h"

// from GMX vec.h //
/*
static inline void clear_mat(matrix a)
{
  a[XX][XX]=a[XX][YY]=a[XX][ZZ]=0.0;
  a[YY][XX]=a[YY][YY]=a[YY][ZZ]=0.0;
  a[ZZ][XX]=a[ZZ][YY]=a[ZZ][ZZ]=0.0;
}
*/
/////////////////////

typedef enum Moltype { OOO, PDB, CRD, COR, XYZ, GRO } moltype;       /*!< \brief enum for molecule file type (pdb/crd/txt/gro(?))*/
typedef enum Trjtype { ooo, dcd, xtc, txt, pdb, crd, cor } trjtype;  /*!< \brief enum for trajectory file type (dcd, xtc, txt, pdb, crd)*/
typedef enum Headertype {LITE, FULL } headertype;
/*typedef enum resType {ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL} restype;*/
/* ------------------------------------------------------------------ */
// ** Periodic Boundary Info
//! Periodic Boundary Structure
/*!
  Info about periodic boundary conditions - size, angles, CPT.
*/
typedef struct pbc
{
   //! Constant Pressure and Temperature flag (size may vary)
   int           cpt;
   //! PBC data: lenght along first axis
   double        a_size;
   //! PBC data: lenght along second axis
   double        b_size;
   //! PBC data: lenght along third axis
   double        c_size;
   //! PBC data: first angle
   double        angle1;
   //! PBC data: second angle
   double        angle2;
   //! PBC data: third angle
   double        angle3;
   //! crystal space group - Z value
   int           zval;
} Pbc   ;
/* ------------------------------------------------------------------ */
// ** Coordinates set
//! Coordinate Set Structure
/*!
  The coordinate set from a trj file, stored in three vectors(x, y and z)
*/
typedef struct coorset
{
   //! number of atoms
   int           nato;
   //! single array for all coordinates
   float        *cords;
   //! vector of x coordinates
   float        *xcoor;   // x coordinates (pointer to cords)
   //! vector of y coordinates
   float        *ycoor;   // y coordinates (pointer to cords+nato)
   //! vector of z coordinates
   float        *zcoor;   // z coordinates (pointer to cords+2*nato)
   //! pointers array for all coordinates (coors[0] = xcoor)
   float       **coors;
   //! PBC data: PBC flag
   int           pbc_flag;
   //! PBC data: PBC data structure
   Pbc          *pbc;
} CoorSet;
/* ------------------------------------------------------------------ */
// ** TRJ header - generic trajectory header
//! Trajectory (trj) info structure
/*!
  Structure with info about a trajectory, mostly collected from the
  headers in the case of dcd files.
*/
typedef struct trjh
{		
   char	      type[5];  /*!< \brief trj type can be CORD (coordinate trj) or VELO (velocities trj) */
   int         nframe;  /*!< \brief _real_ number of frames present in trajectories*/
   int         begstp;  /*!< \brief step number at beginning of trj*/
   int         skpstp;  /*!< \brief number of steps between coordinate frames*/
   unsigned int numstp; /*!< \brief total _real_ number of steps*/
   float        tmstp;  /*!< \brief time between steps (time between frames = tmstp*skpstp)*/
   int         varpbc;  /*!< \brief whether (0/1) variable Periodic Boundary Conditions are employed*/


   int           nato;  /*!< \brief number of atoms in the trj*/
   int         fixatm;  /*!< \brief number of (optionallly) fixed atoms in the trj*/
   int   *fixatm_list;  /*!< \brief list of (optionallly) fixed atoms in the trj*/
   CoorSet   *fixcoor;  /*!< \brief coordinates of (optionallly) fixed atoms in the trj*/

   int            sex;  /*!< \brief whether byte order of trj file is the same as the running machine (0=y, 1=n)*/
   int           size;  /*!< \brief size of headers: used to compute offset for random-access reading of frames*/
   int      sixtyfour;  /*!< \brief whether the file was made by a 64-bit machine*/

   
   int         comnum;  /*!< \brief CHARMM only - number of 80-bytes long comment string*/
   char  comm[12][91];  /*!< \brief CHARMM only - comment strings*/

   int       c_nframe;  /*!< \brief CHARMM only - _claimed_ number of frames; incorrect if trj is not complete*/
   int       c_numstp;  /*!< \brief CHARMM only - _claimed_ number of steps; incorrect if trj is not complete*/

   int	    extra1;     /*!< \brief CHARMM only - (unused) extra field in the header section*/
   int	    extra2;     /*!< \brief CHARMM only - (unused) extra field in the header section*/
   int	    extra3;     /*!< \brief CHARMM only - (unused) extra field in the header section*/
} Trjh;
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// ** general TRJ infos
//! Overall Trajectory (trj) structure
/*!
  Structure with overall info about a trajectory, pointers to the file,
  to the header info,  trjtype, etc.
*/
typedef struct traj
{
   //! pointer to the file (if open)
   FILE        *file;
   //! the file name
   char         filename[1024];
   //! whether it's open or not
   short        open;
   //! if it is an xtc (gromacs) file, this points to the file
   XDRFILE *    xtc;
   //! type (dcd/xtc) of the trj file
   trjtype      filetype;
   //! pointer to hdr structure
   struct trjh *hdr;
   //! the detail level at which info from an xtc file has to be gained
   headertype   hdetail;
   //! while reading a trj, this is the last frame read
   int          previousframe;
   //! a container for a coordinate set extracted from the trajectory
   CoorSet     *crdset;
} Traj	;
/* ------------------------------------------------------------------ */
// ** raw (or row) molecule structure
//! Molecule file (pdb/crd) file row-style structure
/*!
  The info from a molecule file (crd/pdb) are gather here in an
  unstructured form, ie in arrays of variables, each relative to a 
  in a row in the file. The info from each atom is thus contained
  across all the arrays at the relevant index. Easy and fast but crude
  access to data. To enable selections and structured access, this data
  is usually processed and organized in a Molecule structure
*/
typedef struct rawmol
{
   int          *presn;  /*!< presumed residue number (for CRD) */
   int        b_cryst;   /*!< bool: if have crystal unit */
   float            a;   /*!< crystal box element - lenght #1 */
   float            b;   /*!< crystal box element - lenght #2 */
   float            c;   /*!< crystal box element - lenght #3 */
   float        alpha;   /*!< crystal box element - angle #1 */
   float         beta;   /*!< crystal box element - angle #2 */
   float        gamma;   /*!< crystal box element - angle #3 */
   char   s_group[12];   /*!< crystal space group - space group */
   int              z;   /*!< crystal space group - Z value */

   int           nato;   /*!< number of atoms */
   char         *name;   /*!< molecule name */
   char        **atom;   /*!< first column label (in a pdb, usually ATOM) */
   int         *atomn;   /*!< atom number (similar to atmId, which is a char) */
   char        *chain;   /*!< chain name (a single character) */
   char      *chainId;   /*!< chain ID (PDB only) */
   char       *altLoc;   /*!< altLoc field (PDB only) */
   char        *iCode;   /*!< iCode field (PDB only) */
   float        *bFac;   /*!< beta factor (temperature) (PDb only) */
   float       *occup;   /*!< occupancy (PDb only) */
   int          *resn;   /*!< residue number */
   float        *weight; /*!< weight (CRD only) */
   
   short      *hetatm;   /*!< whether it is an hetatm or not */
   char       **atmId;   /*!< atom ID */
   char       **segId;   /*!< segment ID */
   char     **element;   /*!< element (PDB only) */
   char      **charge;   /*!< charge (PDB only) */
   char     **atmtype;   /*!< atom type */
   char     **restype;   /*!< residue type (amino acid) */
   char        **tail;   /*!< eventually, what's left after the known fields */
   short      *segend;   /*!< whether it is the last atom in a segment */
   short      *segbeg;   /*!< whether it is the first atom in a segment */

}  RawMol;
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// ** Atom
//! Single Atom structure
/*!
  All info regarding a single atom are stored in an instance of this
  structure. ato instances are tipically grouped in res instances.
  ato are allocated and filled by processing data in _pdbs instances
*/
typedef struct atom
{
    int        atomn;      // atom number
    char     atmType[5];   // atom type - actually char[4], but with 5 works

    float       xCoor;      // x coordinates
    float       yCoor;      // y coordinates
    float       zCoor;      // z coordinates

    short      hetatm;     // whether it is an hetatm or not
    char       atmId[6];   // atom ID number (char) new from PDB.org!!!
    char       segId[5];   // segment ID (char)     new from PDB.org!!!
    char     element[2];   // new from PDB.org!!!
    char      charge[2];   // new from PDB.org!!!
    char     chainId;      // new from PDB.org!!!
    char      altLoc;      // new from PDB.org!!!
    char       iCode;      // new from PDB.org!!!
    float       bFac;      // beta Factor, usable as index
    float      occup;      // occupancy

    float     weight;      // added for CRD
    int        presn;      // added for CRD
    int           bb;      // backbone (1) or side chain (0) - not widely used
} Atom;
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// ** Residue
//! Single Residue structure
/*!
  info regarding a single residue, like residue type, residue number,
  the number of atoms belonging to it and pointers to the relative 
  _ato instances. res are grouped in segm instances and are allocated
  and filled by processing data in _pdbs instances
*/
typedef struct res
{
    int           nApR;         // number of Atoms Per Residue (nAPR)
    char       resType[5];      // residue type
    int           resn;         // residue number
    int          presn;         // progressive residue number
    Atom         *pAto;         // Atoms 
    int           isAA;         // whether it is an aminoacid or not (has C, CA, N and O atoms)
    void         *pSeg;         // parent structure
    int       segIndex;         // index of segment in molecule->pSeg
    char       segName[5];      // parent segment name
} Res;
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// ** Segment
//! Single Segment structure
/*!
  info regarding a single segment, like segment name, number of atoms
  belonging to it, number of residues belonging to it and pointers to
  the relative _res instances
*/
typedef struct segm
{
    int           nApS;         // number of Atoms (per Segment)
    char       segName[5];      // segment name
    int           nRpS;         // n Residues per Segment (nRpS) 
    Res          *pRes;         // Residues 
} Segm;
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// ** Chain
//! Single Chain structure
/*!
  info regarding a single chain, like chain id, number of atoms
  belonging to it, number of segments belonging to it, number of 
  residues belonging to it and pointers to the relative Segm instances
*/
typedef struct wrd_chain
{
    int           nApC;         // number of Atoms (per Chain)
    int           nRpC;         // number of Atoms (per Chain)
    int           nSeg;         // n Segment per Chain (nRpS) 
    char       chainId;         // chain Id
    Segm        *pSegm;         // Segments 
} wrd_Chain;
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// ** Model
//! Single Model structure
/*!
  info regarding a single model, actually a whole molecule
  Number of atoms (total, per residue, per chain and per segment),
  number of residues (total, per chain and per segment), number
  of segments, number of chains and pointers to appropriate wrd_Chain
  instances, coordinates and a _coor instance and a _pdbs instance with 
  original data. Also a note about the origin and if it comes from an 
  extended format.
  * NOT USED
*/
typedef struct wrd_model
{
    int           nApM;         // number of Atoms (per Model)
    int        modelId;         // model Id (number)
    int           nCpM;         // n Chains per Model (nCpM) 
    wrd_Chain  *pChain;         // Chains 
} wrd_Model;
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// ** Overall Molecular Structure File
//! Single Molecule structure
/*!
  info regarding a molecule. Number of atoms (total, per residue and
  per segment), number of residues (total and per segment), number
  of segments and pointers to appropriate _segm instances, coordinates ad
  a _coor instance and a _pdbs instance with original data. Also a
  note about the origin and if comes from an extended format.
*/
typedef struct _molecule
{
    int             nato;       // number of atoms
    int          nModels;       // number of models
    int          nChains;       // number of chains
    int             nSeg;       // number of segments
    int             nRes;       // number of residues
    int            *nApC;       // number of Atoms per Chain
    int            *nRpC;       // number of Residues per Chain
    int            *nSpC;       // number of Segments per Chain
    int            *nApS;       // number of Atoms per Segment
    int            *nRpS;       // number of Residues per Segment
    int            *nApR;       // number of Atoms per Residue in Molecule - not filled

    wrd_Chain    *chains;       // array of chains
    Segm        *segment;       // array of segments
    Res           **pRes;       // array of residues (pointer to residues inside segments)
    CoorSet         coor;       // coordinates
    CoorSet         coor2;      // additional coordinates (see altLoc)
    CoorSet         coor3;      // additional coordinates (see altLoc)
    CoorSet        *cords;      // used if more than 1 model is present
    RawMol         rawmol;      // raw molecule structural data
    char            *seq1;      // sequence in 1-letter format
    char           **seq3;      // sequence in 3-letter format
    
    int           origin;       // 0 if from pdb, 1 if from crd, 2 if filled for both;
    int         ext_flag;       // whether a crd is EXT(ended)
    
    char         *filled;       // whether the structure has already been used (and thus allocated)
    
    short        altLocB;       // whether a second, alternative set of coordinates is present
    short        altLocC;       // whether a third, alternative set of coordinates is present
    
   /*taken as is from rawmol structure - better have it handy */
   int        b_cryst;   /*!< bool: if have crystal unit */
   float            a;   /*!< crystal box element */
   float            b;   /*!< crystal box element */
   float            c;   /*!< crystal box element */
   float        alpha;   /*!< crystal box element */
   float         beta;   /*!< crystal box element */
   float        gamma;   /*!< crystal box element */
   char   s_group[12];   /*!< crystal space group */
   int              z;   /*!< crystal Z */

} Molecule;
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// ** Selection Structure
//! Atom Selection structure
/*!
  info about an atom selection, such as the selection string which is
  used to scan the molecule file and the arrays where the atom numbers
  are stored once selected. Selestring can also be a file name, with the
  file being a list of atom numbers
*/
typedef struct selection
{
    int   ii;
    //! number of selected atoms
    int   nselatm;
    //! selected atoms' id numbers
    int  *selatm;
    //! selection string
    char  selestring[20480];
} Selection;
/* ------------------------------------------------------------------ */
// ** Frames to extract structure
//! Frame List Structure
/*!
  a list of frames, usually taken from a file
*/
typedef struct intlist
{
   int	   nframes;   // number of frames to process
   int	   *nframe;   // ids of frame to process
} IntList;
/* ------------------------------------------------------------------ */
// ** List of Trajectory files
//! List of Trajectory files
/*!
  a list of trajectory files, usually taken from a file
*/
typedef struct trjlist {
   //! number of trajectories to handle
   int	      ntrj;
   //! names of trajectories to handle
   char	 **trjname;   
} TrjList;
/* ------------------------------------------------------------------ */
// ** List of radii data
//! List of radii data
/*!
  a list of trajectory files, usually taken from a file
*/
typedef struct radlist {
   //! number of radii data
   int	        nradii;
   //! 
   char	      **resname;
   char       **atomname;
   float       *radius;
} RadList;
/* ------------------------------------------------------------------ */
//! substitute for fopen
/*!
  like fopen, opens a file according to mode, prints additional error
  message is anything goes wrong
*/
FILE *  O_File  ( char  *filename,      /*!< name of the file to open*/
                  char  *mode           /*!< file open mode (r/w)*/
                ) ;
/* ------------------------------------------------------------------ */
//! converts from different endian system
/*!
 * in case a different endian system (sex) is detected 
 * this function allows to convert the values to the running
 *  system's
*/
char *swap(int n, char *byte);
/* ------------------------------------------------------------------ */
//! prints info about trajectory file
/*!
  prints to stdout (usually terminal) basic info about a trajectory file
  mostly taken from the headers. Complete support for dcd files, more
  lacking for xtc files. Info have to be already been gathered.
*/
void ShowTrjHeader     (const Traj *traj) ;
/* ------------------------------------------------------------------ */
//! prints info about DCD trajectory file
/*!
  prints to stdout (usually terminal) basic info about a trajectory file
  mostly taken from the headers. Complete support for dcd files, more
  lacking for xtc files. Info have to be already been gathered.
*/
void ShowDcdHeader     ( 
                          //! dcd header
                          struct trjh   *trj
                        ) ;
/* ------------------------------------------------------------------ */
//! reads info about a trajectory (dcd) file
/*!
 reads all available info about a trajectory from a dcd file.
 Number of frames and number of steps are computed from the total
 file length and a singel coor set's length
*/
void ReadDcdHeader     ( FILE   *stream,        /*!< dcd file */  
                         Trjh   *trj            /*!< dcd header */
                        ) ;
/* ------------------------------------------------------------------ */
//! fills a trj header structure
/*!
 fills a trajectory header structure with data from a 
 molecule file structure
*/
void FillDcdHeader     ( Molecule      *,   /*!< molecule */
                          struct trjh   *    /*!< dcd header */
                        ) ;
/* ------------------------------------------------------------------ */
//! write content to a dcd file's header
/*!
 Writes a dcd header structure to a dcd file - used to build a new dcd
 or modify an existing one (vie the Head_Mod functions in datahandler
*/
void WriteDcdHeader     ( FILE          *,   /*!< dcd file */
                          struct trjh   *,   /*!< dcd header */
                          char          *    /*!< field to modify */
                        ) ;
/* ------------------------------------------------------------------ */
//! copies the content of a dcd header
/*!
 copies the content of the source trajectory header structure to the
 target one.
*/
void HdrCp ( struct trjh *source,   /*!< source header */
             struct trjh *target    /*!< target header */
           ) ;
/* ------------------------------------------------------------------ */
//! reads a coordinate set from a DCD trajectory file
/*!
 reads a coordinate set from a DCD trajectory file; called from ReadTRJ_Coor
*/
void ReadDcdCoor        ( FILE           *,   /*!< dcd file */
                         CoorSet        *,   /*!< coordinates */
                         Trjh           *,   /*!< dcd header */
                         int                /*!< frame number */
			) ;
/* ------------------------------------------------------------------ */
//! writes a coordinate set from a DCD trajectory file
/*!
 writes a coordinate set from a DCD trajectory file; called from WriteTRJ_Coor
*/
void WriteDcdCoor      ( FILE           *,   /*!< dcd file */
                         CoorSet        *,   /*!< coordinates */
                         Trjh           *    /*!< dcd header */
			) ;
/* ------------------------------------------------------------------ */
//! reads a pdb-format molecule file
/*!
 only called through GetMolecule, reads molecule file in the pdb format
*/
void ReadPdb            ( char          *,  /*!< pdb name */
                          Molecule      *   /*!< pdb data structure */
                        ) ;
/* ------------------------------------------------------------------ */
//! writes a pdb-format molecule file
/*!
 writes a pdb-format molecule file using the 
 strfile - segm - res - ato cascade
*/
void WritePdb	        ( char          *,   /*!< pdb name */
                          Molecule      *    /*!< pdb data structure */
                        ) ;
/* ------------------------------------------------------------------ */
//! writes a pdb-format molecule file (alt)
/*!
 writes a pdb-format molecule file using the RawMol structure, not the
 Molecule - Segm - Res - Atom cascade
*/
void WritePdb_unstr     ( char *, Molecule * );
/* ------------------------------------------------------------------ */
//! reads a crd-format molecule file
/*!
 only called through GetMolecule, reads molecule file in the pdb format
*/
void ReadCrd            ( char          *,  /*!< crd name */
                          Molecule      *   /*!< crd data structure */
                        ) ;
/* ------------------------------------------------------------------ */
//! writes a crd-format molecule file
/*!
 
*/
void WriteCrd           ( char          *,   /*!< crd name */
                          Molecule      *    /*!< mol data structure */
                        ) ;
/* ------------------------------------------------------------------ */
//! reads a xyz-format molecule file
/*!
 reads an xyz-format molecule in a molecule structure, even if it is
 not strictly one and does not contain all necessary informations.
 called by the GetMolecule function
*/
void ReadXyz            ( char *name,
                          Molecule *molecule
                        );
/* ------------------------------------------------------------------ */
//! writes a xyz-format molecule file
/*!
 writes an xyz-format molecule file. Good for script-based processing
 of coordinates sets
*/
void WriteXyz           ( char           *,   /*!< crd name */
                          Molecule *    /*!< molecule data structure */
                        ) ;

/* ------------------------------------------------------------------ */
//! detects molecule file type
/*!
  whichmoltype detects the molecule file type according to the suffix;
  currently crd (charmm molecule file), pdb (Protein Data Bank) and xyz
  (not really a molecule file) file type are supported
*/
moltype wrd_whichmoltype    ( char          *filename /*!< molecule name */);    
/* ------------------------------------------------------------------ */
//! detects trajectory file type
/*!
  whichtrjtype detects the trajectory file type according to the suffix;
  currently dcd (charmm trajectory file, almost complete support),
  xtc (Gromacs trj file, partial support) and txt (dcd/xtc/crd/pdb files
  list) file type are supported
*/
trjtype wrd_whichtrjtype    ( char          *filename/*!< trajectory name */);     
/* ------------------------------------------------------------------ */
//! reads a molecule
/*!
 general molecule-reading function - calls format-specific functions to
 read in different file formats
*/
Molecule * ReadMolecule ( char          * name, /*!< molecule file name */
                          moltype       filetype/*!< molecule file type */
                        ) ;
/* ------------------------------------------------------------------ */
//! reads a molecule
/*!
 general molecule-reading function - calls format-specific functions to
 read in different file formats
*/
void  GetMolecule       ( char          * name, /*!< molecule file name */
                          Molecule      * molecule,/*!< molecule data structure */
                          moltype       filetype/*!< molecule file type */
                        ) ;
/* ------------------------------------------------------------------ */
//! writes a molecule
/*!
 general molecule-writing function - calls format-specific functions to
 write in different file formats. Loops over chains, segments, residues.

*/
void WriteMolecule_str  ( Molecule      *molecule,      /*!< molecule data structure */
                          char          *name,          /*!< molecule file name */
                          moltype        filetype       /*!< molecule file type */
                        );
/* ------------------------------------------------------------------ */
//! writes a molecule
/*!
 general molecule-writing function - calls format-specific functions to
 write in different file formats. Uses rawmol (order as original file)

*/
void WriteMolecule_unstr( Molecule      *molecule,      /*!< molecule data structure */
                          char          *name,          /*!< molecule file name */
                          moltype        filetype       /*!< molecule file type */
                        );
/* ------------------------------------------------------------------ */
//! writes a molecule
/*!
 general molecule-writing function - calls format-specific functions to
 write in different file formats _ calls WriteMolecule_unstr

*/
void WriteMolecule      ( Molecule      *molecule,      /*!< molecule data structure */
                          char          *name,          /*!< molecule file name */
                          moltype        filetype       /*!< molecule file type */
                        );
/* ------------------------------------------------------------------ */
//! updates data inside molecule structure
/*!
 updates data inside molecule structure to complete the informations 
 about file origin
*/
void UpDateMol          ( Molecule      *molecule       /*!< molecule data structure */);  
/* ------------------------------------------------------------------ */
//! Allocate memory for a molecule
/*!
*/
void AllocMol           ( int           nato,           /*!< number of atoms */
                          Molecule      *molecule       /*! molecule data structure */
);
/* ------------------------------------------------------------------ */
//! Create a Molecule from a selection of another molecule
/*!
*/
Molecule * MakeSeleMol    ( Selection *selstr,  /*!< Selection structure to be used (input)*/
                            Molecule *inmol     /*!< Parent molecule (input) */
) ;
//void MakeSeleMol        ( Selection *selstr,  /*!< Selection structure to be used (input)*/
//                          Molecule *inmol,    /*!< Parent molecule (input) */
//                          Molecule **selemol  /*!< New Molecule made from selection (output) */
//) ;
/* ------------------------------------------------------------------ */
//! clean molecule
/*!
 frees all allocated memory inside a molecule structure
*/
void CleanMolecule      ( Molecule      *molecule );
/* ------------------------------------------------------------------ */
//! deletes molecule
/*!
 frees all allocated memory inside a molecule structure and frees molecule pointer
*/
void DelMolecule        ( Molecule      *molecule );
/* ------------------------------------------------------------------ */
//! opens a trajectory file for reading/writing
/*!
 opens the stream, allocates memory inside the traj structure (for trjh)
 annotates filename
*/
int OpenTrj            ( char  *filename,
                         Traj *traj,
                         char  *mode
                        ) ;
/* ------------------------------------------------------------------ */
//! opens a trajectory file for reading/writing
/*!
 opens the stream, allocates memory inside the traj structure (for trjh)
 annotates filename
*/
Traj * InitTrj          ( char  *filename,
                          char  *mode
                        ) ;
/* ------------------------------------------------------------------ */
//! closes a trajectory file
/*!
 free allocated memory (for trjh), erases filename and closes the stream
*/
void CloseTrj           ( Traj  *traj ) ;
/* ------------------------------------------------------------------ */
//! reads info about a trajectory file
/*!
 
*/
void ReadTrjHeader      ( Traj  *traj);
/* ------------------------------------------------------------------ */
//! reads a coordinate set from a trajectory file
/*!
 general header-reading functions; calls format-specific functions to 
 gather info about the trajectory
*/
void ReadTrjCoor        ( Traj          *traj,  /*!< trajectory to read */
                          CoorSet       *crd,   /*!< coordinates (to be filled) */
                          int            frn    /*!< frame to get */
                        );
/* ------------------------------------------------------------------ */
//! copy headers to a (new) trajectory header structure
/*!
 general trajectory headers copying function: calls format-specific
 functions to copy headers
*/
void CopyTrjHeader ( Trjh *source,         // source header
                     Trjh *target  );      // target header
/* ------------------------------------------------------------------ */
//! write trajectory headers to a trajectory file
/*!
 general trajectory headers writing function: calls format-specific
 functions to write headers to a trajectory file
*/
void WriteTrjHeader     ( Traj  *traj,
                          char  *field );
/* ------------------------------------------------------------------ */
//! write a coordinate set to a trajectory file
/*!
 general trajectory coordinates set writing function: calls
 format-specific functions to write a coordinates set to a trajectory
 file
*/
void WriteTrjCoor ( CoorSet     *crd,   /*!< coordinate set*/
                    Traj        *traj   /*!< trajectory file structure*/
                   );
/* ------------------------------------------------------------------ */
//! old pdb-reading function
/*!
 
*/
void oldReadPDB	        ( FILE          *,   /*!< pdb file */
			  RawMol        *,   /*!< pdb structure */
			  CoorSet       *    /*!< coordinates */
		 	) ;
/* ------------------------------------------------------------------ */
//! old pdb-writing function
/*!
 
*/
void oldWritePDB        ( FILE		*,   /*!< pdb file */
			  RawMol        *,   /*!< pdb structure */
			  CoorSet       *    /*!< coordinates */
                        ) ;
/* ------------------------------------------------------------------ */
//! reads a list of frames
/*!
 takes a text file (.txt) and reads a list of frames to be handled by
 the program (different procedures can use this)
*/
void ReadFramesList     ( char     *,      /*!< framelist filename */
                          IntList  *       /*!< int list structure */
		 	                  ) ;
/* ------------------------------------------------------------------ */
//! reads a list of trajectories
/*!
 takes a text file (.txt) and reads a list of trajectories to be handled
 by the program (different procedures can use this)
*/
void ReadTrjList        ( char          *,      /*!< list file */
			  TrjList        *       /*!< list of trj file structure */
		 	) ;
/* ------------------------------------------------------------------ */
//! extracts selected atoms' id numbers
/*!
 inflate selection string from possible shortened forms to complete
 string to be fed to the selection mechanism itself (range expansion
 not dealt with )
*/
char * inflateSeleString( char *datastring );
/* ------------------------------------------------------------------ */
//! selection difference
/*!
 fills a selection (third argument) with the atoms belonging in the first
 argument not present in the second
*/
int SeleDiff ( Selection *selection1,  Selection *selection2,  Selection *diffsele );
/* ------------------------------------------------------------------ */
//! extracts selected atoms' id numbers
/*!
 high-level selection function: 
*/
int GetSele             ( char          *,      /*!< selection string */
                          Selection     *,      /*!< selection structure */
                          Molecule      *       /*!< molecule file */
                        );
/* ------------------------------------------------------------------ */
//! extracts selected atoms' id numbers
/*!
 selection function: 
*/
int GetSeleAtoms        ( char          *,      /*!< selection string */
                          Selection     *,      /*!< selection structure */
                          Molecule      *       /*!< molecule file */
                        );
/* ------------------------------------------------------------------ */
//! reads a radius list file
/*!
 radius reading function: 
*/
RadList * ReadRadius    (char *filename);
/* ------------------------------------------------------------------ */
//! assign radius based on atomname and residue type
/*!
 radius assigning function: 
*/
int RadAssign           ( RadList *radlist,
                          Molecule *molecule, 
                          double *pdAtomsRadii
                        );
//==============================================================================
//   MEMORY ALLOCATION
//=============================================================================
//! aligned memory allocation
/*!
  walloc allocates num*size memory aligned at 16 bits, to enable
  automatic vectorialization (icc, gcc>4.2). Where available it uses
  memalign, on MacOX malloc. It does NOT set to 0 the allocated memory
*/
inline void * walloc    ( int nblock, int blocksize );
/* ------------------------------------------------------------------ */
//! memory allocation for matrix
/*!
*/
void ** wrd_matAlloc( int rows, int columns, int blocksize );
/* ------------------------------------------------------------------ */
//! memory allocation for float matrix
/*!
*/
float ** wrd_fmatAlloc( int rows, int columns );
/* ------------------------------------------------------------------ */
//! memory free for float matrix allocated with wrd_fmatAlloc
/*!
*/
void  wrd_fmatFree( float **pointer );
/* ------------------------------------------------------------------ */
//! allocates memory for nato atoms in a coorset
void CalloCoor          ( CoorSet * , int nato);
/* ------------------------------------------------------------------ */
//! creates a coorset and calls CalloCoor
CoorSet * InitCoor          ( int nato);
// -----------------------------------------------------------------
//! prints coordinates in a coorset (debugging tool)
void PrintCoor          ( CoorSet *crd, int start, int end );
// -----------------------------------------------------------------
//! prints selected coordinates in a coorset (debugging tool)
void PrintSeleCoor      ( CoorSet *crd, Selection *sele );
// -----------------------------------------------------------------
//! frees memory allocated inside a coorset
void CleanCoor           ( CoorSet *);
// -----------------------------------------------------------------
//! frees memory allocated inside a coorset and free the pointer
void DelCoor           ( CoorSet *);
//! Creates a Molecule and fills it with the target
Molecule * CopyMol( Molecule *mol1 );

int ShowMolInfo( Molecule *molecule);

// -----------------------------------------------------------------
// =============================================================================
//! copies data from a buffer to another
/*< 
 copies data from buffer1 to buffer2 ; from start position to end position ;
 buffer2 is csize long
*/
int copyvar             ( char *buffer1, int start, int end, char *buffer2, int csize );
/* ------------------------------------------------------------------ */
char * expandrange( char *intstring, int *len );
int * readintstring( char *intstring, int tot_int );
/* ------------------------------------------------------------------ */
int wrd_getsuffix (char *file_name, char **suffix);
int wrd_isdattxt(char *filename);
void CopyCoor ( CoorSet *source, CoorSet *target );
void GetSeleCoor2Vecs ( CoorSet *source, float *xcoor, float *ycoor, float *zcoor, Selection *sele );
void GetSeleCoor ( CoorSet *source, CoorSet *target, Selection *sele );
/* ------------------------------------------------------------------ */

/* === GROMACS section === */
/* ------------------------------------------------------------------ */
bool xdrfile_rewind ( XDRFILE *xtc, const char *filename, const char *mode);
/* ------------------------------------------------------------------ */
void ShowXtcHeader( struct trjh *hdr );
/* ------------------------------------------------------------------ */
void ReadXtcCoor  ( Traj *trj, CoorSet *crd, int frn);
/* ------------------------------------------------------------------ */
void WriteXtcCoor ( XDRFILE* xtc, CoorSet *crd, Trjh *hdr );
/* ------------------------------------------------------------------ */
void ReadXtcHeaderFULL ( Traj *trj);
/* ------------------------------------------------------------------ */
void ReadXtcHeaderLITE ( Traj *trj);
/* ------------------------------------------------------------------ */
void *ralloc (void *ptr,int nelem);
/* ------------------------------------------------------------------ */
void rfree (void *ptr);
/* ------------------------------------------------------------------ */
void copymat ( matrix source,matrix dest); 
bool eqmat ( matrix m1,  matrix m2, float round_err); 
/* ------------------------------------------------------------------ */
#endif
