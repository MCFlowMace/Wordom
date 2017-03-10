/******************************************************************************
 *
 *  RING Module for WORDOM
 *  Ring perception software inplementing the method proposed by Hanser:
 *  J. Chem. Inf. Comput. Sci., 1996, 36 (6), pp 1146?1152 DOI:10.1021/ci960322f
 *
 *  Copyright (C) 2013 Simone Conti
 *
 *  RING Module for WORDOM is free software: you can redistribute it and/or 
 *  modify it under the terms of the GNU General Public License as published by 
 *  the Free Software Foundation, either version 3 of the License, or (at your 
 *  option) any later version.
 *
 *  RING Module for WORDOM is distributed in the hope that it will be useful, 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General 
 *  Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "ring.h"

typedef struct Edge {
	int ln;
	int *mb;
	int dl;
} edge;

typedef struct Graph {
	int v_num;
	int e_num;
	int max_l;
	edge *e;
} graph;

graph	ReadGraph(char *filename);
void	RunHanserGraph(char *graph, int max_l, FILE *rings, FILE *stats) ;
void	RunHanserList(char *list, int max_l, FILE *rings, FILE *stats);
graph	HanserSearch(graph *g);
void	PrintRings(FILE *out, graph g);
void	MakeRingStats(FILE *out,graph *g,int k);
void	DeleteGraph(graph *g);


void version() {
	printf("\n");
        printf("RING Module for WORDOM v1.4\n");
        printf("\n");
        printf("Copyright (C) 2013 Simone Conti.\n");
        printf("License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.\n");
        printf("This is free software: you are free to change and redistribute it.\n");
        printf("There is NO WARRANTY, to the extent permitted by law.\n");
        printf("\n");
        printf("Written by Simone Conti.\n");
        printf("\n");
}

void usage() {
        printf("Usage: wordom -ie RING [OPTION]...\n");
        printf("\n");
        printf("Options:\n");
	printf("\t--GRAPH  \t Input filename (if analysing a single graph)\n");
        printf("\t--LIST   \t Input filename (if analysing a list of graphs)\n");
        printf("\t--DIM    \t The maximum dimension of the rings; bigger rings are not perceived; default: 10\n");
        printf("\t--RINGS  \t The output filename for the list of perceived rings; default is 'rings.dat'\n");
        printf("\t--STATS  \t Output filename for ring dimension statistics; if not set no statistics are performed\n");
        printf("\t--HELP   \t Show this help and exit\n");
        printf("\t--VERSION\t Print version information and exit\n");
}

void help() {
        printf("\n");
	printf("\
  RING Module for WORDOM perform a ring searching using the ring perception \n\
algorithm proposed by Hanser et al. (J. Chem. Inf. Comput. Sci., 1996, 36 (6), \n\
pp 1146?1152 DOI:10.1021/ci960322f). The main idea is to find supramolecular \n\
cyclic structures, for example the formation of hexagons or other polygons, but \n\
the algorithm is quite generic and can be used for any ring perception problem. \n\
  The input file consist essentially in a graph in a very easy format: a set of \n\
lines each one containing two numbers. These numbers are the id of an object \n\
(could be an atom or a molecule, or whatever you want), and each line means \n\
that the two object identified by the two number are in contact. For example to \n\
graph file for a pentagon need only the following lines:\n\n\
  1 2\n\
  2 3\n\
  3 4\n\
  4 5\n\
  5 1\n\n\
Very important, each contact must be present only one time and must not be \n\
present duplicates like 1 2 and 2 1. \n\
If you want to analyze only one graph file your WORDOM input file is simply:\n\n\
BEGIN RING\n\
--GRAPH mygraph.inp\n\
END\n\n\
  You can specify two optional parameters: RINGS and DIM. RINGS allow you to \n\
specify an output file where the perceived rings are saved, the default name is \n\
rings.dat. The format of this output is quite simple: the first column contain \n\
the dimension of the ring and the following columns contain the id of each \n\
object that build that ring. The parameter DIM is very useful when you have a \n\
lot of condensed rings. In this case a full ring search will find a very big \n\
number of rings. With the DIM parameter only the rings with a dimension less \n\
than DIM will be saved and printed. The default value for DIM is 10. A complete \n\
WORDOM input file if we are interested only in hexagons could be:\n\n\
BEGIN RING\n\
--GRAPH mygraph.inp\n\
--RINGS myroutput.dat\n\
--DIM 6\n\
END\n\n\
or directly via command line:\n\n\
wordom -ie RING ?GRAPH mygraph.inp --RINGS myoutput.dat --DIM 6\n\n\
  A second way to use this module is via the LIST parameter which is useful if \n\
you need to analyze a lot of graph. In such case you only need to save each \n\
graph in a separate file, and then put all filenames in a list file one per \n\
line, for example:\n\n\
mygraph1.inp\n\
graph2.inp\n\
mythirdgraph.inp\n\
finalgraph.inp\n\n\
Then you can run WORDOM like:\n\n\
BEGIN RING\n\
--LIST mygraphlist.inp\n\
--RINGS myroutput.dat\n\
--DIM 6\n\
END\n\n\
The GRAPH and LIST options are obviously mutually exclusive.\n\n\
  There is a last parameter you can use: STATS. If you specify this parameter \n\
followed by an output filename, the module will make a statistics about ring \n\
dimension for each analyzed graph and will print the result in the specified \n\
output file. This is interesting, for example, when the graphs are representative \n\
of the frame of a trajectory, so you will obtain the statistics of the formed \n\
rings along the trajectory. The format of this file is: in the first column a \n\
number indicating the input file, and then a list of number indicating the \n\
frequency of cycles at each dimension: so in the second column you will find \n\
the number of cycle with dimension 1 (always zero), in the third column with \n\
dimension two (always zero), in the fourth column the triangles, and so on.\n\
  Finally, --HELP will print this help, --VERSION will print version number \n\
and copyright of the current RING Module for WORDOM.\n");
        printf("\n");
        usage();
        printf("\n");
        printf("Examples:\n");
	printf("\tAnalyse only a graph in 'mygraph.dat' file and outputs the perceived \n\trings in rings.dat. No statistics are performed on ring dimension. \n\tUse the default dimension of 10 for maximum ring dimension.\n");
        printf("\t\twordom -ie RING --GRAPH mygraph.dat --RINGS rings.dat\n");
	printf("\n");
	printf("\tAnalyse all graphs listed in 'mylist.dat', save all perceived rings \n\tin 'rings.dat' and save a statistic of rings dimension along all graphs \n\tin 'stats.dat'. Perceive only rings smaller than a pentagon.\n");
        printf("\t\twordom -ie RING --LIST mylist.dat --RINGS rings.dat --STATS stats.dat --DIM 5\n");
        printf("\n");
        printf("Report bugs to <simonecnt@gmail.com>.\n");
        printf("\n");
}




int Compute_Ring(char **input_wrd, int inp_index) {

	FILE *o_stats, *o_rings;
	int graph=0, list=0;
	int dim=10;
	int gotit=0;
	char rings[80]="rings.dat";
	char stats[80]="";
	char input[80]="";
	char *buffer;

	// Parsing command line options
	buffer=input_wrd[inp_index];
	while( strncmp (buffer, "END", 3)) {

		gotit = 0;

		if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#') {
			gotit = 1;
		} else if ( !strncmp(buffer, "--RINGS", 7)) {
			gotit=sscanf(buffer, "--RINGS %80s", rings);
			if (gotit != 1) {
				printf("\n\tERROR: You need to specify an output filename after --RINGS\n");
				exit(EXIT_FAILURE);
			}
     		} else if ( !strncmp(buffer, "--GRAPH", 7)) {
			gotit=sscanf(buffer, "--GRAPH %80s", input);
			if (gotit != 1) {
				printf("\n\tERROR: You need to specify an input filename after --GRAPH\n");
				exit(EXIT_FAILURE);
			}
			graph=1;
		} else if ( !strncmp(buffer, "--LIST", 6)) {
			gotit=sscanf(buffer, "--LIST %80s", input);
			if (gotit != 1) {
				printf("\n\tERROR: You need to specify ad input filename after --LIST\n");
				exit(EXIT_FAILURE);
			}
			list=1;
		} else if ( !strncmp(buffer, "--STATS", 7)) {
			gotit=sscanf(buffer, "--STATS %80s", stats);
			if (gotit != 1) {
				printf("\n\tERROR: You need to specify an output filename after --STAT\n");
				exit(EXIT_FAILURE);
			}
		} else if ( !strncmp(buffer, "--DIM", 5)) {
			gotit=sscanf(buffer, "--DIM %d", &dim);
			if (gotit != 1) {
				printf("\n\tERROR: You need to specify an integer number after --DIM\n");
				exit(EXIT_FAILURE);
			}
		} else if ( !strncmp(buffer, "--HELP", 6)) {
			version();
			help();
			exit(EXIT_SUCCESS);
		} else if ( !strncmp(buffer, "--VERSION", 9)) {
			version();
			exit(EXIT_SUCCESS);
		}
		if( gotit==0 ) {
			printf("\n\tERROR: Sorry, could not understand option: %s\n", buffer);
			exit(EXIT_FAILURE);
		}
		inp_index++;
		buffer=input_wrd[inp_index];
	}

	// Ok, now we have our input option, let check them, print what we are doing and compute
        printf("\n");
        version();
        printf("Parsed options:\n");
	if (list==1 && graph==1) {
		printf("\n\tERROR: You specified both '--GRAPH' and '--LIST' options which are mutually exclusive. If your input file is a graph use '--GRAPH', if is a file containing a list of graphs filenames use '--LIST'\n\n");
		exit(EXIT_FAILURE);
	}
	if (list==1) {
		printf("\tMode: list\n");
	} else if (graph==1) {
		printf("\tMode: graph\n");
	} else {
		printf("\n\tERROR: You have to specify the input filename via the '--GRAPH' or '--LIST' flags. You can chose '--GRAPH' if the input file contains only one graph, and '--LIST' if your input file contains a list of graph filenames.\n\n");
		exit(EXIT_FAILURE);
	}
        printf("\tInput filename: %s\n",input);
	if (input==NULL || input=='\0') {
		printf("\n\tERROR: You didn't specify an input filename. How should I know where to read input data?\n\n");
		exit(EXIT_FAILURE);
	}
        printf("\tOutput for rings: %s\n", rings);
	if (rings==NULL || rings[0]=='\0') {
		printf("\n\tERROR: You didn't specify an output filename with the '--RINGS' flag.\n\n");
	}
        printf("\tDimension: %d\n", dim);
	if (stats==NULL || stats[0]=='\0') {
		printf("\tNo statistics on rings dimension will be performed\n");
	} else {
		printf("\tOutput for statistics: %s\n", stats);
	}
        printf("\n");


	// Open output rings file
	o_rings=fopen(rings,"w");
	if (o_rings==NULL) {
                printf("\n\tERROR!! Is not possible to open %s in mode w. Quit.\n\n", rings);
		exit(EXIT_FAILURE);
	}


	// Open output stats file
	if (!(stats==NULL || stats[0]=='\0')) {
		o_stats=fopen(stats,"w");
		if (o_stats==NULL) {
	                printf("\n\tERROR!! Is not possible to open %s in mode w. Quit.\n", stats);
			exit(EXIT_FAILURE);
		}
	} else {
		o_stats=NULL;
	}

	if (graph) {
		printf("Processing the graph...\n");
		RunHanserGraph(input, dim, o_rings, o_stats);
		if (o_stats!=NULL) fclose(o_stats);
		fclose(o_rings);
		printf("Finish! Thanks for using HanserSearch!\n\n");
		return(EXIT_SUCCESS);

	} else if (list) {
		printf("Processing the list of graphes...\n");
		RunHanserList(input, dim, o_rings, o_stats);
		if (o_stats!=NULL) fclose(o_stats);
		fclose(o_rings);
		printf("Finish! Thanks for using HanserSearch!\n\n");
		return(EXIT_SUCCESS);

	} else {

		printf("\n\tERROR!! Some unknown error occurs... You are very unluky :(\n\n");
		return (EXIT_FAILURE);

	}
}
		
void RunHanserGraph(char *grp, int max_l, FILE *rings, FILE *stats) {

	graph g;

	// Read the input graph
	g=ReadGraph(grp);
	g.max_l=max_l;

	// Performe the search
	printf("Calculating...\n");
	g=HanserSearch(&g);

	// Print rings
	printf("Printing rings...\n"); 
	PrintRings(rings,g);

	// Print statistics on ring dimension
	if (stats!=NULL) {
		printf("Performing statistics...\n");
		MakeRingStats(stats,&g,1);
	}

	// Clean...
	DeleteGraph(&g);

	return;
}


void RunHanserList(char *list, int max_l, FILE *rings, FILE *stats) {

	FILE *fp;	
	char *tmp;
	char row[80];
	char filename[80];
	graph g;
	int i,k;

        // Open input list file
        fp = fopen(list, "r");
        if (fp==NULL) {
                printf("\n\tERROR!! Is not possible to open %s in mode r. Quit.\n", list);
                exit(EXIT_FAILURE);
	}

	// Iterate along the list file
	printf("Calculating...\n");
	k=0;
        while ((tmp=fgets(row, 80, fp))!=NULL) {

		// Get the filename of the graph
		i=0;
		while (row[i]!='\0') {
			if (row[i]!='\n') {
				filename[i]=row[i];
			}
			i++;
		}
		filename[i-1]='\0';
			
		// Read de graph
		g=ReadGraph(filename);
		g.max_l=max_l;

		// Performe the search
		g=HanserSearch(&g);

		// Print the found rings
		fprintf(rings,"%s\n",filename);
		PrintRings(rings,g);

		// Make statistics on ring dimension
		if (stats!=NULL) MakeRingStats(stats,&g,k);

		// Delete the graph
		DeleteGraph(&g);
		k++;
	}

	fclose(fp);

	return;
}

graph HanserSearch(graph *g) {

	int i,j,k;	// Some cycle indexes
	int buf;	// Buffer space for edge vector
	int type;	// Type of the intersection between two edges: 0 -> wrong intersection, >0 -> good intersection
	int v_ch;	// The chosen vertex
	int e_old;	// The number of edges in the current iteration
	int deb=0;	// 1 -> debug printf on, 0 -> debug printf off
	edge *tmp;
	edge *e;
	int v_num, e_num;
	graph out;
	int jj,ii;

	e=g->e;
	e_num=g->e_num;
	v_num=g->v_num;

	// Make some space for new edges
	buf=e_num*2;
	tmp=(edge*)realloc(e,buf*sizeof(edge));
	if (tmp!=NULL) {
		e=tmp;
	} else {
		printf("ERROR!! realloc can't reallocate memory. Exit.\n");
		exit(EXIT_FAILURE);
	}

	while (v_num>0) {
	
		// Choose a vertex (From the last to the firt in order. Not the best solution...)
		v_ch=v_num-1;
		if (deb) printf("v_ch=%d\n",v_ch);
	
		// Foreach couple of edges
		e_old=e_num;
		for (i=0;i<e_old;i++) {
			//if (deb) printf("\t%d\n",i);

			// Check if this edge is 'alive'
			if (e[i].dl==0) {

				for (j=i+1;j<e_old;j++) {

					// Check if also this edge is 'alive'
					if (e[j].dl==0) {

						// Check if the two edges share -exactly- the same chosen edge end then merge the tho edges
						// and also how they should be concatenated (head-head, head-tail, tail-head, tail-tail)

						type=0;
	
						if (e[i].mb[0]==v_ch) {
							if (e[j].mb[0]==v_ch) {
								type=1;
							} else if (e[j].mb[e[j].ln-1]==v_ch) {
								type=2;
							}
						} else if (e[i].mb[e[i].ln-1]==v_ch) {
							if (e[j].mb[0]==v_ch) {
								type=3;
							} else if (e[j].mb[e[j].ln-1]==v_ch) {
								type=4;
							}
						}

						for (ii=1;ii<e[i].ln-1;ii++) {
							for (jj=1;jj<e[j].ln-1;jj++) {
								if (e[i].mb[ii]==e[j].mb[jj]) {
									type=0;
								}
							}
						}

						if (deb) printf("type=%d\n",type);

						// If yes merge the two edge	
						if (type>0) {

							if (e_num>=buf) {
								if (deb) printf("realloc\n");
								buf*=2;
								tmp=(edge*)realloc(e,buf*sizeof(edge));
								if (tmp!=NULL) {
									e=tmp;
								} else {
									printf("ERROR!! realloc can't reallocate memory. Exit.\n");
									exit(EXIT_FAILURE);
								}
							}
							e[e_num].ln = e[i].ln + e[j].ln - 1;

							// If the lenght of the new path is less than the imposed maximum lenght, continue.
							//   Otherwise do nothing, that is erase this path

							if (e[e_num].ln <= g->max_l+1) {

							e[e_num].mb = (int*)malloc((e[e_num].ln+1)*sizeof(int));
							e[e_num].dl = 0;
					
							switch  (type) {

							case 2:
								if (deb) printf("\t\t\tIf2: %d, %d\n",i,j);
								for (k=0;k<e[i].ln;k++) {
									e[e_num].mb[k] = e[i].mb[e[i].ln-1-k];
								}
								for (k=0;k<e[j].ln;k++) {
									e[e_num].mb[e[i].ln-1+e[j].ln-1-k] = e[j].mb[k];
								}
								break;
			
							case 4:
								if (deb) printf("\t\t\tIf4: %d, %d\n",i,j);
								for (k=0;k<e[i].ln;k++) {
									e[e_num].mb[k] = e[i].mb[k];
								}
								for (k=0;k<e[j].ln;k++) {
									e[e_num].mb[e[i].ln-1+e[j].ln-1-k] = e[j].mb[k];
								}
								break;

							case 3: 
								if (deb) printf("\t\t\tIf3: %d, %d\n",i,j);
								for (k=0;k<e[i].ln;k++) {
									e[e_num].mb[k] = e[i].mb[k];
								}
								for (k=1;k<e[j].ln;k++) {
									e[e_num].mb[e[i].ln+k-1] = e[j].mb[k];
								}
								break; 

							case 1:
								if (deb) printf("\t\t\tIf1: %d, %d\n",i,j);
								for (k=0;k<e[i].ln;k++) {
									e[e_num].mb[k] = e[i].mb[e[i].ln-1-k];
								}
								for (k=0;k<e[j].ln;k++) {
									e[e_num].mb[e[i].ln-1+k] = e[j].mb[k];
								}
								break;
							}

							// If we formed a cicle set it as 'dead'
							if (e[e_num].mb[0]==e[e_num].mb[e[e_num].ln-1]) {
								if (deb) printf("CICLE!");
								e[e_num].dl=1;
							}
							e_num++;

							// And print some debug info...
							if (deb) {
								printf("\t\t\t");
								for (k=0;k<e[e_num-1].ln;k++) {
									printf("%d\t",e[e_num-1].mb[k]);
								}
								printf("\n");
							}

							} // end if lenght
				
						} // end if type>0
					} // end if j alive
				} // end for j
			} //end if i alive
		} // end for i

		// Delete merged edges
		for (i=0;i<e_old;i++) {
			if (e[i].mb[0]==v_ch || e[i].mb[e[i].ln-1]==v_ch) {
				e[i].dl=1;
			}
		}

		v_num--;
	}
	
	out.v_num=0;
	out.e_num=e_num;
	out.e=e;
	out.max_l=g->max_l;

	return out;
}

graph ReadGraph(char *filename) {

	FILE *fp;
        int n, m;
        int nr;
       // int val;
        char *tmp;
        char row[80];
	graph g;
	int buf=64;
	edge *mem;

        // Open input graph file
        fp = fopen(filename, "r");
        if (fp==NULL) {
                printf("\n\tERROR!! ReadGraph: is not possible to open %s in mode r. Quit.\n", filename);
                exit(EXIT_FAILURE);
	}

	g.v_num=0;
	g.e=(edge*)malloc(buf*sizeof(edge));
	g.e_num=0;

        // Read all elements
        while ((tmp=fgets(row, 80, fp))!=NULL) {
                if (row[0]!='#' && row[0]!='\n') {
                        if ((nr=sscanf(row,"%d %d %*d",&n,&m))!=2) {
                                printf("ERROR! ReadGraph: fscanf return %d instead of 3\n",nr);
                                exit(EXIT_FAILURE);
                        }
			if (n>g.v_num) {
				g.v_num=n;
			}
			if (m>g.v_num) {
				g.v_num=m;
			}
			if (g.e_num>=buf) {
				buf*=2;
				mem=(edge*)realloc(g.e,buf*sizeof(edge));
				if (tmp!=NULL) {
					g.e=mem;
				} else {
					printf("ERROR!! ReadGraph: realloc can't reallocate memory. Exit.\n");
					exit(EXIT_FAILURE);
				}
			}
			g.e[g.e_num].ln=2;
			g.e[g.e_num].mb=(int*)malloc(2*sizeof(int));
			g.e[g.e_num].mb[0]=n;
			g.e[g.e_num].mb[1]=m;
			g.e[g.e_num].dl=0;
			g.e_num++;
                }
        }

	fclose(fp);

	return g;
}

void PrintRings(FILE *out, graph g) {

	int i,j;

	for (i=0;i<g.e_num;i++) {
		if (g.e[i].mb[0]==g.e[i].mb[g.e[i].ln-1]) {
			fprintf(out,"%d\t",g.e[i].ln-1);
			for (j=0;j<g.e[i].ln-1;j++) {
				fprintf(out,"%d\t",g.e[i].mb[j]);
			}
			fprintf(out,"\n");
		}
	}

	return;
}

void MakeRingStats(FILE *out, graph *g, int k) {

	int i;
	int *v;

	v = (int*)malloc(g->max_l*sizeof(int));

	fprintf(out,"%d ",k);
	for (i=0;i<g->max_l;i++) {
		v[i]=0;
	}
	for (i=0;i<g->e_num;i++) {
		if (g->e[i].mb[0] == g->e[i].mb[g->e[i].ln-1]) {
			v[g->e[i].ln-2]++;
		}
	}
	for (i=0;i<g->max_l;i++) {
		fprintf(out,"%d ",v[i]);
	}
	fprintf(out,"\n");

	free(v);

	return;
}


void DeleteGraph(graph *g) {

	int i;

	for (i=0;i<g->e_num;i++) {
		free(g->e[i].mb);
	}

	free(g->e);
	
	return;

}
