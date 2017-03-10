#ifndef GCLUSTER_H
#define GLUCSTER_H

//Uses the gpu to calculate the DRMS of the distance matrices and assigns the frames to clusters
int gClusterDrms (struct inp_Cluster *inp_cluster,float *distance);

//Uses the gpu to calculate the RMSD and assigns the frames to clusters
int gClusterRmsd (struct inp_Cluster *inp_cluster,float *distance);

#endif /*GCLUSTER_H*/
