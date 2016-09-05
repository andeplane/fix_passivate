/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_passivate.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "domain.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include <iostream>
using namespace std;
using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;

#define SILICONTYPE 1
#define OXYGENTYPE 2
#define HYDROGENTYPE 3
#define THREE_BODY_CUTOFF 3.0
#define SI_O_DISTANCE 1.65
#define O_H_DISTANCE 0.95 // Filip Sund: https://www.duo.uio.no/handle/10852/41122 (page 68)
/* ---------------------------------------------------------------------- */
FixPassivate::FixPassivate(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  time_depend = 1;
  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
  newAtomsList = NULL;
  newAtomsList = (MyAtom *) memory->srealloc(newAtomsList, 1e6*sizeof(MyAtom),
                                                "passivate:newAtomList");
  numNewAtoms = 0;

  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use fix_deposit unless atoms have IDs");
  oxygenPositions = new double*[3];
  for(int i=0; i<3; i++) oxygenPositions[i] = new double[3];
  
  if(comm->nprocs > 1) {
    error->all(FLERR,"Slutt å være løk. Fix::Passivate kan ikke brukes på mer enn 1 kjerne.");
  }
}

/* ---------------------------------------------------------------------- */

FixPassivate::~FixPassivate()
{
  memory->destroy(newAtomsList);
}

int FixPassivate::setmask()
{
	int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

void FixPassivate::init()
{
	int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 0;
}

void FixPassivate::init_list(int id, NeighList *ptr)
{
	list = ptr;
}

// void FixPassivate::handle3CoordinatedSilicon(int atomIndex, int *oxygenNeighbors) {
void FixPassivate::handle3CoordinatedSilicon(int atomIndex, double *oxygen0, double *oxygen1, double *oxygen2) {
  double *silicon = atom->x[atomIndex];
  double v1[3];
  double v2[3];
  double normal[3];
  // double *oxygen0 = atom->x[oxygenNeighbors[0]];
  // double *oxygen1 = atom->x[oxygenNeighbors[1]];
  // double *oxygen2 = atom->x[oxygenNeighbors[2]];

  // First find the two vectors spaning the plane defined by the three oxygen atoms
  sub3(oxygen1, oxygen0, v1);
  sub3(oxygen2, oxygen0, v2);
  cross3(v1, v2, normal); // Normal vector of that plane
  norm3(normal); // Normalize

  // Find if normal vector points towards silicon atom or away from it
  double vectorFromOxygenToSilicon[3];
  sub3(silicon, oxygen0, vectorFromOxygenToSilicon);
  int sign = dot3(vectorFromOxygenToSilicon, normal) > 0 ? 1 : -1;

  double newOxygenPosition[3];
  newOxygenPosition[0] = sign*normal[0]*SI_O_DISTANCE + silicon[0];
  newOxygenPosition[1] = sign*normal[1]*SI_O_DISTANCE + silicon[1];
  newOxygenPosition[2] = sign*normal[2]*SI_O_DISTANCE + silicon[2];
  createOHPair(newOxygenPosition, silicon);
  num3CoordinatePassivated++;
}

void FixPassivate::handle2CoordinatedSilicon(int atomIndex, double *oxygen0, double *oxygen1) {
  double *silicon = atom->x[atomIndex];
  double v1[3];
  double v2[3];
  double normal[3];

  sub3(oxygen0, silicon, v1);
  sub3(oxygen1, silicon, v2);
  cross3(v1, v2, normal); // Normal vector of that plane
  norm3(normal); // Normalize

  double v1PlusV2[3];
  add3(v1,v2,v1PlusV2);
  norm3(v1PlusV2);
  double alpha = 109.5 / 180. * M_PI;
  double omegaN = SI_O_DISTANCE*sin(alpha * 0.5);
  double omegaO = SI_O_DISTANCE*cos(alpha * 0.5);

  double newOxygenPosition[3];
  newOxygenPosition[0] = silicon[0] + omegaN*normal[0] - omegaO*v1PlusV2[0];
  newOxygenPosition[1] = silicon[1] + omegaN*normal[1] - omegaO*v1PlusV2[1];
  newOxygenPosition[2] = silicon[2] + omegaN*normal[2] - omegaO*v1PlusV2[2];
  createOHPair(newOxygenPosition, silicon);

  newOxygenPosition[0] = silicon[0] - omegaN*normal[0] - omegaO*v1PlusV2[0];
  newOxygenPosition[1] = silicon[1] - omegaN*normal[1] - omegaO*v1PlusV2[1];
  newOxygenPosition[2] = silicon[2] - omegaN*normal[2] - omegaO*v1PlusV2[2];
  createOHPair(newOxygenPosition, silicon);
  num2CoordinatePassivated++;
}

void FixPassivate::handle1CoordinatedSilicon(int atomIndex, double *oxygen0) {
  double *silicon = atom->x[atomIndex];
  double normal[3];
  
  sub3(silicon, oxygen0, normal);
  norm3(normal); // Normalize

  double alpha = 109.5 / 180. * M_PI;
  double circleMidPoint[3];
  circleMidPoint[0] = silicon[0] + cos(alpha*0.5)*normal[0];
  circleMidPoint[1] = silicon[1] + cos(alpha*0.5)*normal[1];
  circleMidPoint[2] = silicon[2] + cos(alpha*0.5)*normal[2];
  double tangentInPlaneContaining3MissingOxygens[3];
  findVectorInPlaneDefinedByNormalVector(normal, tangentInPlaneContaining3MissingOxygens);
  norm3(tangentInPlaneContaining3MissingOxygens);

  double cosAlphaHalf = cos(alpha*0.5);
  double lengthOfVectorsInPlane = sqrt(SI_O_DISTANCE*SI_O_DISTANCE-cosAlphaHalf*cosAlphaHalf);
  double newOxygenPosition[3];
  newOxygenPosition[0] = circleMidPoint[0] + lengthOfVectorsInPlane*tangentInPlaneContaining3MissingOxygens[0];
  newOxygenPosition[1] = circleMidPoint[1] + lengthOfVectorsInPlane*tangentInPlaneContaining3MissingOxygens[1];
  newOxygenPosition[2] = circleMidPoint[2] + lengthOfVectorsInPlane*tangentInPlaneContaining3MissingOxygens[2];
  createOHPair(newOxygenPosition, silicon);
  
  handle2CoordinatedSilicon(atomIndex, oxygen0, newOxygenPosition);

  num1CoordinatePassivated++;
  num2CoordinatePassivated--;
}

void FixPassivate::createOHPair(double *oxygenPosition, double *siliconPosition) {
  // Create first OH pair
  MyAtom &newOxygenAtom = newAtomsList[numNewAtoms];
  newOxygenAtom.type = OXYGENTYPE;
  newOxygenAtom.position[0] = oxygenPosition[0];
  newOxygenAtom.position[1] = oxygenPosition[1];
  newOxygenAtom.position[2] = oxygenPosition[2];
  numNewAtoms++;

  double vectorFromSiliconToNewOxygen[3];
  sub3(oxygenPosition, siliconPosition, (double*)vectorFromSiliconToNewOxygen);
  norm3(vectorFromSiliconToNewOxygen);

  MyAtom &newHydrogenAtom = newAtomsList[numNewAtoms];
  newHydrogenAtom.type = HYDROGENTYPE;
  newHydrogenAtom.position[0] = newOxygenAtom.position[0] + vectorFromSiliconToNewOxygen[0]*O_H_DISTANCE;
  newHydrogenAtom.position[1] = newOxygenAtom.position[1] + vectorFromSiliconToNewOxygen[1]*O_H_DISTANCE;
  newHydrogenAtom.position[2] = newOxygenAtom.position[2] + vectorFromSiliconToNewOxygen[2]*O_H_DISTANCE;
  numNewAtoms++;
}

void FixPassivate::findVectorInPlaneDefinedByNormalVector(double *normal, double *result) {
    // Check if any of the components of the normal vector is close or equal to zero
    bool smallComponentInNormalVector[3];
    int numberOfSmallComponentsInNormalVector = 0;
    double epsilon = 1e-10;
    for (int i = 0; i < 3; i++) {
        if (fabs(normal[i]) < epsilon) {
            smallComponentInNormalVector[i] = true;
            numberOfSmallComponentsInNormalVector++;
        }
    }

    // Creating vector in the plane defined by normal
    switch(numberOfSmallComponentsInNormalVector) {
    case 0:
        result[0] = 1.0;
        result[1] = 1.0;
        result[2] = -((normal[0] + normal[1])/normal[2]);
        break;
    case 1:
        if (smallComponentInNormalVector[0]) {
            result[0] = 1.0;
            result[1] = 1.0;
            result[2] = 1.0 - ((normal[0] + normal[1])/normal[2]);
        } else {
            result[0] = 1.0 -((normal[1] + normal[2])/normal[0]);
            result[1] = 1.0;
            result[2] = 1.0;
        }
        break;
    case 2:
        result[0] = ( smallComponentInNormalVector[0] ? 1.0 : (normal[1] + normal[2])/normal[0] );
        result[1] = ( smallComponentInNormalVector[1] ? 1.0 : (normal[0] + normal[2])/normal[1] );
        result[2] = ( smallComponentInNormalVector[2] ? 1.0 : (normal[0] + normal[1])/normal[2] );
        break;
    }
}

void FixPassivate::handleSilicon(int atomIndex) {
  double **x = atom->x;
  int *type = atom->type;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int *jlist = firstneigh[atomIndex];
  int jnum = numneigh[atomIndex];

  double xi = x[atomIndex][0];
  double yi = x[atomIndex][1];
  double zi = x[atomIndex][2];

  int oxygenNeighbors[3];
  int oxygenNeighborCount = 0;

  double oxygen0[3];
  double oxygen1[3];
  double oxygen2[3];

  for (int jj = 0; jj < jnum; jj++) {
    int j = jlist[jj];
    j &= NEIGHMASK;
    double xj = x[j][0];
    double yj = x[j][1];
    double zj = x[j][2];
    double dx = xi - xj;
    double dy = yi - yj;
    double dz = zi - zj;
    double dr2 = dx*dx + dy*dy + dz*dz;
    double dr = sqrt(dr2);
    if(type[j] == OXYGENTYPE && dr < THREE_BODY_CUTOFF) {
      if(oxygenNeighborCount == 0) {
        oxygen0[0] = xj;
        oxygen0[1] = yj;
        oxygen0[2] = zj;
      } else if(oxygenNeighborCount == 1) {
        oxygen1[0] = xj;
        oxygen1[1] = yj;
        oxygen1[2] = zj;
      } else if(oxygenNeighborCount == 2) {
        oxygen2[0] = xj;
        oxygen2[1] = yj;
        oxygen2[2] = zj;
      }
      oxygenNeighbors[oxygenNeighborCount] = j;
      oxygenNeighborCount += 1;
    }
  }

  if(oxygenNeighborCount == 3) {
    handle3CoordinatedSilicon(atomIndex, oxygen0, oxygen1, oxygen2);
    // handle3CoordinatedSilicon(atomIndex, oxygenNeighbors);
  } else if(oxygenNeighborCount == 2) {
    handle2CoordinatedSilicon(atomIndex, oxygen0, oxygen1);
  } else if(oxygenNeighborCount == 1) {
    handle1CoordinatedSilicon(atomIndex, oxygen0);
  } 
}


void FixPassivate::handleOxygen(int atomIndex) {
  double **x = atom->x;
  int *type = atom->type;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int *jlist = firstneigh[atomIndex];
  int jnum = numneigh[atomIndex];

  double xi = x[atomIndex][0];
  double yi = x[atomIndex][1];
  double zi = x[atomIndex][2];

  int numSiliconNeighbors = 0;
  int siliconAtomIndex = -1;
  
  for (int jj = 0; jj < jnum; jj++) {
    int j = jlist[jj];
    j &= NEIGHMASK;
    double xj = x[j][0];
    double yj = x[j][1];
    double zj = x[j][2];
    double dx = xj - xi;
    double dy = yj - yi;
    double dz = zj - zi;
    double dr2 = dx*dx + dy*dy + dz*dz;
    double dr = sqrt(dr2);
    if(type[j] == SILICONTYPE && dr < THREE_BODY_CUTOFF) {
      siliconAtomIndex = j;
      numSiliconNeighbors++;
    }
  }

  if(numSiliconNeighbors == 1) {
    // Create hydrogen to passivate dangling oxygen bond
    double xj = x[siliconAtomIndex][0];
    double yj = x[siliconAtomIndex][1];
    double zj = x[siliconAtomIndex][2];
    double deltaR[3];
    deltaR[0] = xj - xi;
    deltaR[1] = yj - yi;
    deltaR[2] = zj - zi;
    norm3(deltaR);

    MyAtom &newAtom = newAtomsList[numNewAtoms];
    newAtom.type = HYDROGENTYPE;
    newAtom.position[0] = xi - deltaR[0]*O_H_DISTANCE;
    newAtom.position[1] = yi - deltaR[1]*O_H_DISTANCE;
    newAtom.position[2] = zi - deltaR[2]*O_H_DISTANCE;
    numNewAtoms++;
    numOxygenPassivated++;
  }
}

void FixPassivate::addNewAtoms() {
  for(int i = 0; i<numNewAtoms; i++) {
    MyAtom &newAtom = newAtomsList[i];
    double *position = (double*)&newAtom.position;
    domain->remap(position);
    atom->avec->create_atom(newAtom.type, position);
    int n = atom->nlocal - 1;
    atom->tag[n] = maxtag_all + i+1;
  }
}

void FixPassivate::end_of_step()
{
  // if (next_reneighbor != update->ntimestep) return;
  neighbor->build_one(list);
  num1CoordinatePassivated = 0;
  num2CoordinatePassivated = 0;
  num3CoordinatePassivated = 0;
  numOxygenPassivated = 0;
  find_maxid();

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int *type = atom->type;
  // loop over full neighbor list of my atoms
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    
    if(type[i] == SILICONTYPE) {
      handleSilicon(i);
    } else if(type[i] == OXYGENTYPE) {
      handleOxygen(i);
    }
  }


  char output[4096];
  sprintf(output, "Added %d new atoms. \n Oxygens passivated: %d. \n"
    "1 coordinated silicon passivated: %d. \n"
    "2 coordinated silicon passivated: %d. \n"
    "3 coordinated silicon passivated: %d. \n", numNewAtoms, numOxygenPassivated, num1CoordinatePassivated, num2CoordinatePassivated, num3CoordinatePassivated);

  if (screen)  fputs(output,screen);
  if (logfile) fputs(output,logfile);

  atom->nghost = 0;
  addNewAtoms();
  atom->natoms += numNewAtoms;
  next_reneighbor = 0;

  if (atom->map_style) {
    atom->map_init();
    atom->map_set();
  }
}

void FixPassivate::find_maxid()
{
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  tagint max = 0;
  for (int i = 0; i < nlocal; i++) max = MAX(max,tag[i]);
  MPI_Allreduce(&max,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
}








