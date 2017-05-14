// ------------------------------------------------------------------
// Copyright (C) 2009  Francesca Fanelli, Angelo Felline, Michele Seeber
// University of Modena and Reggio Emilia - ITALY
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
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

#define iMAX_DDP_NUMBER_OF_LINKS 50

#ifndef PSN
#define PSN

struct       simplecluster
{
  int        nclusters;
  int       *cluspop;
  int      **clusters;
};
// ------------------------------------------------------------------
//  Protein Structure Network Graph structure
struct inp_psn
{
  int        **ppiIntAtomPairs;
  int          tot_nresidues;
  int         *reslength;
  int         *reslength2;
  int        **atmList;
  int        **ppiInteractions;
  int         *piNumResInteractions;
  int          iNumOfIntMinSteps;
  int          iMaxResInteractions;
  int          iFrameNum;
  int          iHubContCutoff;
  int          iMaxNumOfLink;
  int         *piLargestClusterSize;
  int          iWarningFlag;
  int          iVerboseFlag;
  int          iTerminiFlag;
  int          iProxCutOff;
  int          iPBCFlag;
  int          iNumOfNoLinkPairs;
  int        **ppiNoLinkPairs;
  int          iNumOfForceLinkPairs;
  int          iNumOfNoProxSeg;
  int         *piPDBResNum;
  int        **ppiForceLinkPairs;
  int          iNumOfMerge;
  int          iMergeMaxResNum;
  int       ***pppiStableResInt;
  int       ***pppiStableResIntMergeClust;
  int      ****ppppiHubCorr;
  int       ***pppiClusters;
  int          iMaxClustNum;
  int        **ppiTmpStableHubs;
  int        **ppiStableHubs;
  int        **ppiStableHubsMergeClust;
  int        **ppiTempStableHubsMergeClust;
  int         *piNodeDegree;
  int          iNumOfSelRes;
  int          iNumOfParam;
  int          iNumOfFrames;
  int          iNumOfNodes;
  int          iGetIminFlag;
  int          iSecondRound;
  int          iDecNum;
  int          iBestIcValue;
  int        **ppiResResIntFreq;
  int          iIntType;
  int          iHubEqFlag;
  int          iWritePDBFlag;
  int          iParamDBFlag;
  int          iUsedParamDBFlag;

  int          iMergeClustFlag;
  int          iMergeMinPop;
  int         *piNodeClust;
  int         *piClustsPop;
  int          iNumOfPoxInt;
  int        **ppiFrameLinks;
  
  int         *piNodeClusters;
  int         *piClustSize;
  int        **ppiIntMatrix;
  int          iGetIminAlgo;
  int          iNumOfParamDBEntries;
  int         *piUsedDBParamIndices;
  
  char         title[999];
  char         cRawFileName[99];
  char       **ppcMergedResidues;
  char       **pcParamResVect;
  char         cIminRange[999];
  char      ***pppcIntAtomNamePairs;
  char       **pcNoProxSeg;
  char       **pcSeleResSegId;
  char         cParamDBFileName[999];
  char       **pcParamDBResVect;
  char       **ppcResId2Lab;
  
  float      **ppfIntStrength;
  float      **ppfHubsIntStrength;
  float       *res_norm;
  float        fDistCutoff;
  float        fIntMinStart, fIntMinStop, fIntMinStep;
  float      **ppfAvgResInt;
  float     ***pppfStableResIntStrength;
  float     ***pppfStableResIntStrengthMergeClust;
  float       *pfParamValVect;
  float       *pfParamDBValVect;
  float        fPreIcValue, fPostIcValue;
  float        fPreIcValueDelta, fPostIcValueDelta;
  float       *pfIminValues;
  
  double       fPCNSize;
  
  struct       simplecluster cluster;
  
  Selection    sele;

  FILE        *outfile;

  time_t       time_tToday;
  
};
// ==================================================================
/*!
 * linkwalk finds the clusters defined as bunch of nodes (residudues) connected among them
 * no spacial requirement nor min numbero f links. A long line of residue connected
 * one-by-one is good, as is a clot of residues each connected to a central hub only
 * arguments:
 * nres : number of total residues in molecule
 * linklist : array - row index (+1) is residues, row memebers (+1) are residues linked with it
 * nlinks : vector of number of links for each residue
 *cluster : simplecluster structure where results are stored - have to be already allocated
*/
int     linkwalk(int nres, int **linklist, int *nlinks, struct simplecluster *cluster);
// ---------------------------------------------------------------------
float   GetAANormFactor(char *resType);
// ---------------------------------------------------------------------
int     Read_iPSG ( char **input, int inp_index, struct inp_psn *inp_psn, char *printout, Molecule *molecule, int iNumOfFrames, struct sopt *OPT);
// ---------------------------------------------------------------------
int     Compute_PSG ( struct inp_psn *inp_psn, struct sopt *OPT, CoorSet *trj_crd, char *outprint, Molecule *molecule );
// ---------------------------------------------------------------------
int     Post_PSG(struct inp_psn *inp_psn, int iNumOfFrames, Molecule *molecule, struct sopt *OPT);
// ---------------------------------------------------------------------
float   GetUserNormFactor(struct inp_psn *inp_psn, char *resType);
// ---------------------------------------------------------------------
float   GetParamDBNormFactor(struct inp_psn *inp_psn, char *resType);
// ---------------------------------------------------------------------
float   GetNormFactor(struct inp_psn *inp_psn, char *resType, int isAAFlag);
// ---------------------------------------------------------------------
void    GetImin(struct inp_psn *inp_psn, struct sopt *OPT, Molecule *molecule);
// ---------------------------------------------------------------------
void    GetImin2(struct inp_psn *inp_psn, struct sopt *OPT, Molecule *molecule);
// ---------------------------------------------------------------------
double  GetBigClsSize(FILE *FRawFile, int *piNodeClusters, int *piClustSize, float fImin, int iNumOfNodes, int iNumOfSeleRes, int **ppiIntMatrix);
// ---------------------------------------------------------------------
double  GetBigClsSize2(struct inp_psn *inp_psn, float fImin);
// ---------------------------------------------------------------------
int     UpdateClusters(int iTmpRes1, int iTmpRes2, int *piNodeClusters, int iNumOfNodes, int iLastCluster);
// ---------------------------------------------------------------------
void    GetIminClearAll(int iNumOfNodes, int *piNodeClusters, int *piClustSize, int iLastCluster);
// ---------------------------------------------------------------------
void    MergeClusters(struct inp_psn *inp_psn, int iIntMinIterationNum, float fImin);
// ---------------------------------------------------------------------
int     SetNodeClustVect(struct inp_psn *inp_psn);
// ---------------------------------------------------------------------
void    LoadParamDB(struct inp_psn *inp_psn);
// ---------------------------------------------------------------------
void    ReportUsedParamDBEntries(struct inp_psn *inp_psn, FILE *FileHandler, int iReportMode);
// ---------------------------------------------------------------------

/* PSN Path Section */

// === Function Prototypes =============================================
int     PSNPATH(char **, int);                                          // Main function
int     GetPSNParam(char **, int);                                      // Gets some psn parameters
int     SequenceCheck();                                                // Check seq in PSN & Corr files
int     NodeLabToIndex(char *);                                         // Returns the node index associated to a node label
int     ArraySorter(const void *, const void *);                        // Used to sort an array of float values
int     ClustCheck();                                                   // Checks if iRes1 and iRes2 are in the same cluster
void    GetCorrData();                                                  // Gets Corr data
void    ClearAll();                                                     // Clears all matrices and vectors
void    SetIntVar();                                                    // Set internal variables
void    GetCorrRes();                                                   // Gets Correlated residues
void    GetPaths();                                                     // Calculates paths
void    FindPath();                                                     // Find Shortest path
void    GenPath(int, int);                                              // Generates last path string
void    SaveData();                                                     // Save Paths on file
void    SavePath();                                                     // Save last path
void    SavePathLite();                                                 // Used with --MODE LITE
void    GetFrameInt();                                                  // Load res-res int of next frame
void    ClearAllLite();                                                 // Clears all matrices and vectors in --MODE LITE
void    ClearIntMatrix();                                               // Clear IntMatrix in --MODE LITE
void    PSNPathError(int);                                              // Error Handler
void    GetPSGData();                                                   // Loads PSG data from an avg* file
void    LoadAvgLabIdx();                                                // Loads Node labels and indexes
void    LoadStableLinks();                                              // Loads Stable Links
void    GetUniqIntFreqList();                                           // Gets the list of interaction frequency values
void    IntMatrixFilter();                                              // Filters ppfIntMatrix for mix-paths-search
float   GetBadFrames();                                                 // Gets Bad frames
float   GetIcritFromFile(int);                                          // Returns the Icritic value from raw file
float   GetLinkWeight(float);                                           // Calculates the weight of passed link on the basis of iWeightFlag and res-res interaction strength
void    GetStrongestInteraction();                                      // Finds the strongest res-res interaction strength
void    ApplyMergeInfoInPathSearch();                                   // Get MergeClust info from rawFile and use 'em all
// =====================================================================


/* PSN Param Section */
// === PSNParam Structure ==============================================
struct inp_psnparam
{
  int                   iAvgMode;                                       // Average Mode
  int                   iNumOfIgnore;                                   // Res to ignore in links
  int                   iNumOfMol;                                      // Number of Molecules
  int                   iNumOfTrj;                                      // Number of Trajectories
  int                   iVerboseFlag;                                   // If 1 a detailed file is created
  int                   iNumOfWarning;                                  // The number of warning
  int                   iSetNumOfTargetAtomsFlag;                       // Flag to Set the Number of target residues
  int                   iNumOfTargetAtoms;                              // Number of target residues
  
  char                **ccIgnoreVect;                                   // Vect of --IGNORE
  char                **ccMolFileVect;                                  // Vect of --MOL file name
  char                **ccMolSeleVect;                                  // Vect of --MOL sele string
  char                **ccTrjFileVect;                                  // Vect of --TRJ trajectory file name
  char                **ccRefFileVect;                                  // Vect of --TRJ reference file name
  char                **ccTrjSeleVect;                                  // Vect of --TRJ sele string
  char                  cTarget[1024];                                  // "Residue" to be parametrized
  char                  cTitle[1024];                                   // Title, used for verbose file
  char                  cVerboseFileName[1024];                         // Verbose filename
  
  float                 fDistCutoff;                                    // Distance CutOff
  float                 fAvgMaxInteractions;                            // Averaged number of max interactions
  
  Selection             sele;                                           // selection structure
  Molecule             *molecule;                                       // molecule structure ptr
  Traj                 *xxx;                                           // trajectory structure ptr
  
  FILE                 *FVerboseFile;                                   // Verbose file handler
  
  time_t                time_tToday;                                    // Date and Time structure
};
// =====================================================================

// === Function Prototypes =============================================
int  InitPSNParam(char **, int);
int  CheckIgnore(struct inp_psnparam *, char *);
void CalcPSNParamFromMol(struct inp_psnparam *, int);
void CalcPSNParamFromTrj(struct inp_psnparam *, int);
// =====================================================================

#endif
