/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2015 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'molecule_properties.c' is part of RASPA-2.0

    RASPA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RASPA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *************************************************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include "constants.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "simulation.h"
#include "ewald.h"
#include "potentials.h"
#include "utils.h"
#include "output.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "grids.h"
#include "spacegroup.h"
#include "spectra.h"
#include "statistics.h"
#include "input.h" // for the cutoff of order parameter
#include "movies.h" // try: for CutOffROP

REAL ComputeBondDistanceFramework(int index)
{
  int A,B;
  VECTOR dr;
  REAL r;

  if(index>=Framework[CurrentSystem].NumberOfBonds[CurrentFramework])
    fprintf(stderr, "Error: bond index too large\n");

  A=Framework[CurrentSystem].Bonds[CurrentFramework][index].A;
  B=Framework[CurrentSystem].Bonds[CurrentFramework][index].B;

  dr.x=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.x-
       Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.x;
  dr.y=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.y-
       Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.y;
  dr.z=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.z-
       Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.z;
  dr=ApplyBoundaryCondition(dr);
  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
  return r;
}

REAL ComputeBendAngleFramework(int index)
{
  int A,B,C;
  VECTOR posA,posB,posC,Rab,Rbc;
  REAL CosTheta,Theta,rab,rbc;

  if(index>=Framework[CurrentSystem].NumberOfBends[CurrentFramework])
    fprintf(stderr, "Error: framework bend index too large\n");

  A=Framework[CurrentSystem].Bends[CurrentFramework][index].A;
  B=Framework[CurrentSystem].Bends[CurrentFramework][index].B;
  C=Framework[CurrentSystem].Bends[CurrentFramework][index].C;

  posA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position;
  posB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Position;
  posC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position;

  Rab.x=posA.x-posB.x;
  Rab.y=posA.y-posB.y;
  Rab.z=posA.z-posB.z;
  Rab=ApplyBoundaryCondition(Rab);
  rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
  Rab.x/=rab;
  Rab.y/=rab;
  Rab.z/=rab;

  Rbc.x=posC.x-posB.x;
  Rbc.y=posC.y-posB.y;
  Rbc.z=posC.z-posB.z;
  Rbc=ApplyBoundaryCondition(Rbc);
  rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
  Rbc.x/=rbc;
  Rbc.y/=rbc;
  Rbc.z/=rbc;

  CosTheta=(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
  CosTheta=SIGN(MIN2(fabs(CosTheta),(REAL)1.0),CosTheta);
  Theta=acos(CosTheta);

  return (RAD2DEG*Theta);
}

REAL ComputeTorsionAngleFramework(int index)
{
  int A,B,C,D;
  VECTOR posA,posB,posC,posD;
  VECTOR Rab,Rbc,Rcd,Pb,Pc;
  REAL pb2,rpb1,pc2,rpc1,pbpc,rrbc;
  REAL cost,sint,theta;

  if(index>=Framework[CurrentSystem].NumberOfTorsions[CurrentFramework])
    fprintf(stderr, "Error: framework torsion index too large\n");

  A=Framework[CurrentSystem].Torsions[CurrentFramework][index].A;
  B=Framework[CurrentSystem].Torsions[CurrentFramework][index].B;
  C=Framework[CurrentSystem].Torsions[CurrentFramework][index].C;
  D=Framework[CurrentSystem].Torsions[CurrentFramework][index].D;

  posA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position;
  posB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Position;
  posC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position;
  posD=Framework[CurrentSystem].Atoms[CurrentFramework][D].Position;

  Rab.x=posA.x-posB.x;
  Rab.y=posA.y-posB.y;
  Rab.z=posA.z-posB.z;

  Rbc.x=posB.x-posC.x;
  Rbc.y=posB.y-posC.y;
  Rbc.z=posB.z-posC.z;

  Rcd.x=posC.x-posD.x;
  Rcd.y=posC.y-posD.y;
  Rcd.z=posC.z-posD.z;

  rrbc=1.0/sqrt(Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z);

  // construct first dihedral vector
  Pb.x=Rab.y*Rbc.z-Rab.z*Rbc.y;
  Pb.y=Rab.z*Rbc.x-Rab.x*Rbc.z;
  Pb.z=Rab.x*Rbc.y-Rab.y*Rbc.x;
  pb2=Pb.x*Pb.x+Pb.y*Pb.y+Pb.z*Pb.z;
  rpb1=1.0/sqrt(pb2);

  // construct second dihedral vector
  Pc.x=Rbc.y*Rcd.z-Rbc.z*Rcd.y;
  Pc.y=Rbc.z*Rcd.x-Rbc.x*Rcd.z;
  Pc.z=Rbc.x*Rcd.y-Rbc.y*Rcd.x;
  pc2=Pc.x*Pc.x+Pc.y*Pc.y+Pc.z*Pc.z;
  rpc1=1.0/sqrt(pc2);

  // determine dihedral angle
  pbpc=Pb.x*Pc.x+Pb.y*Pc.y+Pb.z*Pc.z;
  cost=pbpc*rpb1*rpc1;
  sint=(Rbc.x*(Pc.y*Pb.z-Pc.z*Pb.y)+Rbc.y*(Pb.x*Pc.z-Pb.z*Pc.x)
        +Rbc.z*(Pc.x*Pb.y-Pc.y*Pb.x))*(rpb1*rpc1*rrbc);
  theta=atan2(sint,cost);
  if(theta<0.0) theta+=2.0*M_PI;
  return theta*180.0/M_PI;
}

REAL ComputeThetaPhi(int system, int i, int j, VECTOR distance_com, FILE *OrderParamPtr)
{
// function the compute the two parameters
// used by ComputeAromaticOrderParameter function 
  int A,B,C,D,E,F;
  int k,l; // for looping
  VECTOR posA,posB,posC,posD,posE,posF;
  VECTOR Rab,Rbc,Rde,Ref;
  VECTOR CenterABC, CenterDEF; // get the center of ring 
  REAL theta,phi, phie;
  int carbon_count,count; 
  char atom_type[32];
  int *SelectedAtoms;
  SelectedAtoms=(int*)calloc(3,sizeof(int)); // just pick 3 atoms
  REAL pb2, rpb1, pe2, rpe1, rpbpe_cross, sint, cost;
  VECTOR Pb, Pe, pbpe_cross;
  VECTOR pbcom, pecom;
  REAL rr_com;
  REAL sinp, cosp, sinpe, cospe;
  
  // First compute theta, which is the angle between two benzene ring planes
  // first select 3 atoms from the ring
  // loop over all atoms in that molecule
  count=0;
  carbon_count=0;
  while (carbon_count < 3){
    //fprintf(stderr, "looping the atom number %d\n", count);
    //fprintf(stderr, "Selected Atom Type is: %d\n", Adsorbates[system][i].Atoms[count].Type);
    strcpy(atom_type, PseudoAtoms[Adsorbates[system][i].Atoms[count].Type].Name);
    //fprintf(stderr, "Selected Pseudo atom is: %s\n", atom_type);
    //One may use either UA model or OPLS detailed model
    if ((strncasecmp(atom_type, "C_xyl", strlen("C_xyl")) == 0) || (strncasecmp(atom_type, "C_OPLS_benz", strlen("C_OPLS_benz")) == 0)){
      //fprintf(stderr, "carbon_count is: %d\n", carbon_count);
      SelectedAtoms[carbon_count] = count;
      //fprintf(stderr, "Selected Atom is %d\n", SelectedAtoms[carbon_count]);
      carbon_count++;
    }
    count++;
  }
  //after the while loop, the selected atoms are ready.
  A = SelectedAtoms[0];
  B = SelectedAtoms[1];
  C = SelectedAtoms[2];
  posA=Adsorbates[system][i].Atoms[A].Position;
  posB=Adsorbates[system][i].Atoms[B].Position;
  posC=Adsorbates[system][i].Atoms[C].Position;
  // if EBCB&fictional first bead, use the fictional bead
  if ((EBCBMC && (!UseEBCBMCEnergyFirstBead)) || (Framework[system].FrameworkModel==NONE))
  {
    CenterABC = Adsorbates[system][i].Atoms[0].Position;
  }
  else
  {
    CenterABC = GetAdsorbateCenterOfMass(i);
  }

  //fprintf(stderr, "A is: %d	%lf %lf %lf\n", A, posA.x, posA.y, posA.z);
  //fprintf(stderr, "B is: %d	%lf %lf %lf\n", B, posB.x, posB.y, posB.z);
  // let's then get three atoms from molecule j
  count=0;
  carbon_count=0;
  while (carbon_count < 3){
    strcpy(atom_type, PseudoAtoms[Adsorbates[system][j].Atoms[count].Type].Name);
    if ((strncasecmp(atom_type, "C_xyl", strlen("C_xyl")) == 0) || (strncasecmp(atom_type, "C_OPLS_benz", strlen("C_OPLS_benz")) == 0)){
      SelectedAtoms[carbon_count] = count;
      carbon_count++;
    }
    count++;
  }
  D = SelectedAtoms[0];
  E = SelectedAtoms[1];
  F = SelectedAtoms[2];
  posD=Adsorbates[system][j].Atoms[D].Position;
  posE=Adsorbates[system][j].Atoms[E].Position;
  posF=Adsorbates[system][j].Atoms[F].Position;
  // if EBCB&fictional first bead, use the fictional bead
  if ((EBCBMC && (!UseEBCBMCEnergyFirstBead)) || (Framework[system].FrameworkModel==NONE))
  {
    CenterDEF = Adsorbates[system][j].Atoms[0].Position;
  }
  else
  {
    CenterDEF = GetAdsorbateCenterOfMass(j);
  }
  // then calculate the distances
  // distance between A, B
  Rab.x = posA.x - posB.x;
  Rab.y = posA.y - posB.y;
  Rab.z = posA.z - posB.z;
  // distance between B, C
  Rbc.x = posB.x - posC.x;
  Rbc.y = posB.y - posC.y;
  Rbc.z = posB.z - posC.z;
  // distance between D, E
  Rde.x = posD.x - posE.x;
  Rde.y = posD.y - posE.y;
  Rde.z = posD.z - posE.z;
  // distance between E, F
  Ref.x = posE.x - posF.x;
  Ref.y = posE.y - posF.y;
  Ref.z = posE.z - posF.z;
  // since Rab, Rbc, Rde, Ref are all 
  // within a molecule, no need for PBD
/*  Rab = ApplyBoundaryConditionUnitCell(Rab);
  Rbc = ApplyBoundaryConditionUnitCell(Rbc);
  Rde = ApplyBoundaryConditionUnitCell(Rde);
  Ref = ApplyBoundaryConditionUnitCell(Ref); */
  // compute normal for the first ring plane
  Pb.x=Rab.y*Rbc.z-Rab.z*Rbc.y;
  Pb.y=Rab.z*Rbc.x-Rab.x*Rbc.z;
  Pb.z=Rab.x*Rbc.y-Rab.y*Rbc.x;
  pb2=Pb.x*Pb.x+Pb.y*Pb.y+Pb.z*Pb.z;
  rpb1=1.0/sqrt(pb2);
  // then compute the normal for the second plane
  Pe.x=Rde.y*Ref.z-Rde.z*Ref.y;
  Pe.y=Rde.z*Ref.x-Rde.x*Ref.z;
  Pe.z=Rde.x*Ref.y-Rde.y*Ref.x;
  pe2=Pe.x*Pe.x+Pe.y*Pe.y+Pe.z*Pe.z;
  rpe1=1.0/sqrt(pe2);
  // theta is defined as atan2(norm of cross, dot)
  // first calculate the cross product of Pb and Pe
  pbpe_cross.x = Pb.y*Pe.z-Pb.z*Pe.y;
  pbpe_cross.y = Pb.z*Pe.x-Pb.x*Pe.z;
  pbpe_cross.z = Pb.x*Pe.y-Pb.y*Pe.x;
  // calculate the norm for pbpe_cross
  rpbpe_cross = 1.0/sqrt(pbpe_cross.x*pbpe_cross.x + pbpe_cross.y*pbpe_cross.y + pbpe_cross.z*pbpe_cross.z);
  // sin of theta is defined as 
  // dot of pbpe_cross and itself, then divide by the norms
  sint = sqrt(pbpe_cross.x*pbpe_cross.x + pbpe_cross.y*pbpe_cross.y + pbpe_cross.z*pbpe_cross.z)*rpb1*rpe1;
  // then calculate the dot product
  cost = (Pb.x*Pe.x + Pb.y*Pe.y + Pb.z*Pe.z)*rpb1*rpe1;
  //fprintf(stderr, "sint is: %lf\n", sint);
  //fprintf(stderr, "cost is: %lf\n", cost);
  theta = atan2(sint, cost);
  if(theta<0.0) theta+=2.0*M_PI;

  // pffff, then we can calculate phi, 
  // the angle between normal of one ring plane and the line of the two COMs.
  // the distance between COMs vector is passed down
  // choose to use Pb as the normal vector
  // first calculate the norm of distance_COMs
  rr_com=1.0/sqrt(distance_com.x*distance_com.x + distance_com.y*distance_com.y + distance_com.z*distance_com.z);
  // calculate the cross between Pb and distance_com
  pbcom.x = Pb.y*distance_com.z-Pb.z*distance_com.y;
  pbcom.y = Pb.z*distance_com.x-Pb.x*distance_com.z;
  pbcom.z = Pb.x*distance_com.y-Pb.y*distance_com.x;
  
  sinp = sqrt(pbcom.x*pbcom.x + pbcom.y*pbcom.y + pbcom.z*pbcom.z)*rr_com*rpb1;
  cosp = (Pb.x*distance_com.x + Pb.y*distance_com.y + Pb.z*distance_com.z)*rr_com*rpb1;
  phi = atan2(sinp, cosp);
  if(phi<0.0) phi+=2.0*M_PI;

  // now, we compute the third angle, which is the other phi angle
  pecom.x = Pe.y*distance_com.z-Pe.z*distance_com.y;
  pecom.y = Pe.z*distance_com.x-Pe.x*distance_com.z;
  pecom.z = Pe.x*distance_com.y-Pe.y*distance_com.x;
  sinpe = sqrt(pecom.x*pecom.x + pecom.y*pecom.y + pecom.z*pecom.z)*rr_com*rpe1;
  cospe = (Pe.x*distance_com.x + Pe.y*distance_com.y + Pe.z*distance_com.z)*rr_com*rpe1;
  phie = atan2(sinpe, cospe);
  if(phie<0.0) phie+=2.0*M_PI;
  // finally, compute the order parameters 
  // it will be the averaging of 2nd term of lagendre polynomials
  // comming from theta and phi
  // Why 2*phi? This is because of T-shape stacking, 
  // 2*phi will guarentee that which-ever the normal we choose,
  // we get the same order parameter.
  // actually, exp(cos(4phi)) is better... maybe more change later...
  // see powerpoint 090919 for reference of change
  
  // then, print these out to the file pointer
  // first print out the positions 
  fprintf(OrderParamPtr, "---------------------------------------\n");
  fprintf(OrderParamPtr, "A: %lf %lf %lf\n", posA.x, posA.y, posA.z);
  fprintf(OrderParamPtr, "B: %lf %lf %lf\n", posB.x, posB.y, posB.z);
  fprintf(OrderParamPtr, "C: %lf %lf %lf\n", posC.x, posC.y, posC.z);
  fprintf(OrderParamPtr, "D: %lf %lf %lf\n", posD.x, posD.y, posD.z);
  fprintf(OrderParamPtr, "E: %lf %lf %lf\n", posE.x, posE.y, posE.z);
  fprintf(OrderParamPtr, "F: %lf %lf %lf\n", posF.x, posF.y, posF.z);
  fprintf(OrderParamPtr, "Center of ABC: %lf %lf %lf\n", CenterABC.x, CenterABC.y, CenterABC.z);
  fprintf(OrderParamPtr, "Center of DEF: %lf %lf %lf\n", CenterDEF.x, CenterDEF.y, CenterDEF.z);
  //fprintf(OrderParamPtr, "Rab: %lf %lf %lf\n", Rab.x, Rab.y, Rab.z);
  //fprintf(OrderParamPtr, "Rbc: %lf %lf %lf\n", Rbc.x, Rbc.y, Rbc.z);
  //fprintf(OrderParamPtr, "Rde: %lf %lf %lf\n", Rde.x, Rde.y, Rde.z);
  //fprintf(OrderParamPtr, "Ref: %lf %lf %lf\n", Ref.x, Ref.y, Ref.z);
  //fprintf(OrderParamPtr, "Pb: %lf %lf %lf\n", Pb.x, Pb.y, Pb.z);
  //fprintf(OrderParamPtr, "Pe: %lf %lf %lf\n", Pe.x, Pe.y, Pe.z);
  //fprintf(OrderParamPtr, "pbpe_cross: %lf %lf %lf\n", pbpe_cross.x, pbpe_cross.y, pbpe_cross.z);
  fprintf(OrderParamPtr, "Center VECTOR: %lf %lf %lf\n", distance_com.x, distance_com.y, distance_com.z);
  fprintf(OrderParamPtr, "Molecule Pairs: %d	%d\n", i, j);
  fprintf(OrderParamPtr, "Angles are: %lf     %lf     %lf\n", theta*180.0/M_PI, phi*180.0/M_PI, phie*180.0/M_PI);
  free(SelectedAtoms);
  return theta*180.0/M_PI;

}

void ComputeAromaticOrderParameter(int system, int cycle, int Cycle_type)
{
  int i,j;
  int counter;
  VECTOR pos_dummy_i, pos_dummy_j;
  VECTOR distance_com;
  REAL rr_com;
  double angle; // angle calculated from ComputeThetaPhi function
  // for doing histograming
  int* bins;
  int num_bins = 180;
  double min_angle = 0;
  double max_angle = 180;
  double hist_width = (max_angle-min_angle)/num_bins;
  double* mid_points;
  mid_points = (double*)calloc(num_bins, sizeof(double)); // for performing integrations
  for (i = 0; i < num_bins; i++){
     mid_points[i] = min_angle + (i+0.5)*hist_width;
  }
  // allocate the histogram bins
  bins = (int*)calloc(num_bins, sizeof(int));
  // for calculating P2
  double* normfactor;
  double sum_factor = 0;
  double* unnormal;
  double sum_normal = 0;
  double P2;
  normfactor = (double*)calloc(num_bins, sizeof(double));
  unnormal = (double*)calloc(num_bins, sizeof(double));

  char *MoleculeName_i,*MoleculeName_j;
  char buffer[256];
  FILE *OrderParamPtr; // where information about theta, phi are dumped
  // Cycle_type indicates which kind of cycle you are in:
  // 0 = Initialization cycle
  // 1 = Equilibration cycle
  // 2 = Production cycle
  // It will change the name of the printed file
 
  fprintf(stderr, "CutOff ROP is %lf\n", CutOffRingOP);

  sprintf(buffer, "RingOrderParameter/System_%d", system);
  mkdir(buffer, S_IRWXU);
  if(Cycle_type == 0) // Initialization
  {
    sprintf(buffer,"RingOrderParameter/System_%d/Initial_Cycle_%d.txt",system,cycle);
  }
  else if(Cycle_type == 1) // Equilibration
  {
    sprintf(buffer,"RingOrderParameter/System_%d/Equil_Cycle_%d.txt",system,cycle);
  }
  else //Production 
  {
    sprintf(buffer,"RingOrderParameter/System_%d/Product_Cycle_%d.txt",system,cycle);
  }
  OrderParamPtr=fopen(buffer,"w");
  fprintf(OrderParamPtr, "Aromatic Order Parameter Analysis\n");
  fprintf(OrderParamPtr, "Molecule_1 Molecule_2 Theta Phi OrderParameter_S");

  // re-zero all the counters and parameters
  counter = 0;
  for(i = 0; i < NumberOfAdsorbateMolecules[system]; i++){
     for(j = (i+1); j < NumberOfAdsorbateMolecules[system]; j++){
        MoleculeName_i = Components[Adsorbates[system][i].Type].Name;
	MoleculeName_j = Components[Adsorbates[system][j].Type].Name;
        // we want the distance between the center of rings
        // for xylenes, center of rings are labeled as the 18th atom, 
        // which is dummy atom
        // if it is doing EBCBMC and Not using 1st bead energy, use the 1st bead position (ring center)
        // else, get the COM of the molecule (may not be the ring center)
        // sometimes for bulk phase OP, we still use the fictional 1st bead
        if ((EBCBMC && (!UseEBCBMCEnergyFirstBead)) || (Framework[system].FrameworkModel==NONE))
        {
          pos_dummy_i = Adsorbates[system][i].Atoms[0].Position;
          pos_dummy_j = Adsorbates[system][j].Atoms[0].Position;
        }
        else
        {
          pos_dummy_i = GetAdsorbateCenterOfMass(i);
          pos_dummy_j = GetAdsorbateCenterOfMass(j);
        }
	// get the distance between center of masses
	distance_com.x = pos_dummy_i.x - pos_dummy_j.x;
	distance_com.y = pos_dummy_i.y - pos_dummy_j.y;
	distance_com.z = pos_dummy_i.z - pos_dummy_j.z;
	distance_com=ApplyBoundaryConditionUnitCell(distance_com);
        rr_com=sqrt(distance_com.x*distance_com.x + distance_com.y*distance_com.y + distance_com.z*distance_com.z);
        if(rr_com <= CutOffRingOP){
	  // threshold set to be 6 angstrom, 5 is too small...
	  counter++;
	  angle = ComputeThetaPhi(system, i, j, distance_com, OrderParamPtr);
	  // assign the angle value to the corresponding bin
	  bins[(int)((angle - min_angle)/hist_width)] += 1;
	}
     } 
  }
  // finally, perform averaging on the order parameters
  for (i = 0; i < num_bins; i++){
        normfactor[i] = bins[i]*sin(mid_points[i]/(180/M_PI));
        unnormal[i] = bins[i]*sin(mid_points[i]/(180/M_PI))*SQR(cos(mid_points[i]/(180/M_PI)));
        sum_factor += normfactor[i];
        sum_normal += unnormal[i];
  }
  P2 = (3*sum_normal/sum_factor - 1)/2;

  fprintf(OrderParamPtr, "--------------------------------------------------------------------------------\n");
  fprintf(OrderParamPtr, "In the box we have %d pairs with Center Of Mass less than 5A\n", counter);
  fprintf(OrderParamPtr, "After averaging, the P2 is: %lf\n", P2);
  // record these two into output file
  RingOrderPairs=counter;
  RingOrderParam=(REAL)P2;

  fclose(OrderParamPtr);
  free(bins);
  free(normfactor);
  free(unnormal);
  free(mid_points);
}


REAL ComputeBondDistanceAdsorbate(int m,int index)
{
  int Type,A,B;
  VECTOR posA,posB,dr;
  REAL r;


  Type=Adsorbates[CurrentSystem][m].Type;
  if(index>=Components[Type].NumberOfBonds)
    fprintf(stderr, "Error: bond index too large\n");

  A=Components[Type].Bonds[index].A;
  B=Components[Type].Bonds[index].B;
  posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
  posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
  dr.x=posA.x-posB.x;
  dr.y=posA.y-posB.y;
  dr.z=posA.z-posB.z;
  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
  return r;
}

REAL ComputeBondDistanceCation(int m,int index)
{
  int Type,A,B;
  VECTOR posA,posB,dr;
  REAL r;


  Type=Cations[CurrentSystem][m].Type;
  if(index>=Components[Type].NumberOfBonds)
    fprintf(stderr, "Error: bond index too large\n");

  A=Components[Type].Bonds[index].A;
  B=Components[Type].Bonds[index].B;
  posA=Cations[CurrentSystem][m].Atoms[A].Position;
  posB=Cations[CurrentSystem][m].Atoms[B].Position;
  dr.x=posA.x-posB.x;
  dr.y=posA.y-posB.y;
  dr.z=posA.z-posB.z;
  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
  return r;
}

REAL ComputeUreyBradleyDistanceAdsorbate(int m,int index)
{
  int Type,A,C;
  VECTOR posA,posC,dr;
  REAL r;

  Type=Adsorbates[CurrentSystem][m].Type;
  if(index>=Components[Type].NumberOfUreyBradleys)
    fprintf(stderr, "Error: Urey-Brdley index too large\n");

  A=Components[Type].UreyBradleys[index].A;
  C=Components[Type].UreyBradleys[index].C;
  posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
  posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
  dr.x=posA.x-posC.x;
  dr.y=posA.y-posC.y;
  dr.z=posA.z-posC.z;
  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
  return r;
}

REAL ComputeUreyBradleyDistanceCation(int m,int index)
{
  int Type,A,C;
  VECTOR posA,posC,dr;
  REAL r;

  Type=Cations[CurrentSystem][m].Type;
  if(index>=Components[Type].NumberOfUreyBradleys)
    fprintf(stderr, "Error: Urey-Bradley index too large\n");

  A=Components[Type].UreyBradleys[index].A;
  C=Components[Type].UreyBradleys[index].C;
  posA=Cations[CurrentSystem][m].Atoms[A].Position;
  posC=Cations[CurrentSystem][m].Atoms[C].Position;
  dr.x=posA.x-posC.x;
  dr.y=posA.y-posC.y;
  dr.z=posA.z-posC.z;
  r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
  return r;
}

REAL ComputeBendAngleAdsorbate(int m,int index)
{
  int Type,A,B,C;
  VECTOR posA,posB,posC,Rab,Rbc;
  REAL rab,rbc,theta;

  Type=Adsorbates[CurrentSystem][m].Type;
  if(index>=Components[Type].NumberOfBends)
    fprintf(stderr, "Error: bend index too large\n");

  A=Components[Type].Bends[index].A;
  B=Components[Type].Bends[index].B;
  C=Components[Type].Bends[index].C;

  posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
  posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
  posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;

  Rab.x=posA.x-posB.x;
  Rab.y=posA.y-posB.y;
  Rab.z=posA.z-posB.z;
  rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
  Rab.x/=rab;
  Rab.y/=rab;
  Rab.z/=rab;

  Rbc.x=posC.x-posB.x;
  Rbc.y=posC.y-posB.y;
  Rbc.z=posC.z-posB.z;
  rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
  Rbc.x/=rbc;
  Rbc.y/=rbc;
  Rbc.z/=rbc;

  theta=acos(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
  return (RAD2DEG*theta);
}

REAL ComputeBendAngleCation(int m,int index)
{
  int Type,A,B,C;
  VECTOR posA,posB,posC,Rab,Rbc;
  REAL rab,rbc,theta;

  Type=Cations[CurrentSystem][m].Type;
  if(index>=Components[Type].NumberOfBends)
    fprintf(stderr, "Error: bend index too large\n");

  A=Components[Type].Bends[index].A;
  B=Components[Type].Bends[index].B;
  C=Components[Type].Bends[index].C;

  posA=Cations[CurrentSystem][m].Atoms[A].Position;
  posB=Cations[CurrentSystem][m].Atoms[B].Position;
  posC=Cations[CurrentSystem][m].Atoms[C].Position;

  Rab.x=posA.x-posB.x;
  Rab.y=posA.y-posB.y;
  Rab.z=posA.z-posB.z;
  rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
  Rab.x/=rab;
  Rab.y/=rab;
  Rab.z/=rab;

  Rbc.x=posC.x-posB.x;
  Rbc.y=posC.y-posB.y;
  Rbc.z=posC.z-posB.z;
  rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
  Rbc.x/=rbc;
  Rbc.y/=rbc;
  Rbc.z/=rbc;

  theta=acos(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
  return theta*180.0/M_PI;
}

REAL ComputeTorsionAngleAdsorbate(int m,int index)
{
  int Type,A,B,C,D;
  VECTOR posA,posB,posC,posD;
  VECTOR Rab,Rbc,Rcd,Pb,Pc;
  REAL pb2,rpb1,pc2,rpc1,pbpc,rrbc;
  REAL cost,sint,theta;

  Type=Adsorbates[CurrentSystem][m].Type;
  if(index>=Components[Type].NumberOfTorsions)
    fprintf(stderr, "Error: torsion index too large\n");

  A=Components[Type].Torsions[index].A;
  B=Components[Type].Torsions[index].B;
  C=Components[Type].Torsions[index].C;
  D=Components[Type].Torsions[index].D;

  posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
  posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
  posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
  posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

  Rab.x=posA.x-posB.x;
  Rab.y=posA.y-posB.y;
  Rab.z=posA.z-posB.z;

  Rbc.x=posB.x-posC.x;
  Rbc.y=posB.y-posC.y;
  Rbc.z=posB.z-posC.z;

  Rcd.x=posC.x-posD.x;
  Rcd.y=posC.y-posD.y;
  Rcd.z=posC.z-posD.z;

  rrbc=1.0/sqrt(Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z);

  // construct first dihedral vector
  Pb.x=Rab.y*Rbc.z-Rab.z*Rbc.y;
  Pb.y=Rab.z*Rbc.x-Rab.x*Rbc.z;
  Pb.z=Rab.x*Rbc.y-Rab.y*Rbc.x;
  pb2=Pb.x*Pb.x+Pb.y*Pb.y+Pb.z*Pb.z;
  rpb1=1.0/sqrt(pb2);

  // construct second dihedral vector
  Pc.x=Rbc.y*Rcd.z-Rbc.z*Rcd.y;
  Pc.y=Rbc.z*Rcd.x-Rbc.x*Rcd.z;
  Pc.z=Rbc.x*Rcd.y-Rbc.y*Rcd.x;
  pc2=Pc.x*Pc.x+Pc.y*Pc.y+Pc.z*Pc.z;
  rpc1=1.0/sqrt(pc2);

  // determine dihedral angle
  pbpc=Pb.x*Pc.x+Pb.y*Pc.y+Pb.z*Pc.z;
  cost=pbpc*rpb1*rpc1;
  sint=(Rbc.x*(Pc.y*Pb.z-Pc.z*Pb.y)+Rbc.y*(Pb.x*Pc.z-Pb.z*Pc.x)
        +Rbc.z*(Pc.x*Pb.y-Pc.y*Pb.x))*(rpb1*rpc1*rrbc);
  theta=atan2(sint,cost);
  if(theta<0.0) theta+=2.0*M_PI;
  return theta*180.0/M_PI;
}

REAL ComputeTorsionAngleCation(int m,int index)
{
  int Type,A,B,C,D;
  VECTOR posA,posB,posC,posD;
  VECTOR Rab,Rbc,Rcd,Pb,Pc;
  REAL pb2,rpb1,pc2,rpc1,pbpc,rrbc;
  REAL cost,sint,theta;

  Type=Cations[CurrentSystem][m].Type;
  if(index>=Components[Type].NumberOfTorsions)
    fprintf(stderr, "Error: torsion index too large\n");

  A=Components[Type].Torsions[index].A;
  B=Components[Type].Torsions[index].B;
  C=Components[Type].Torsions[index].C;
  D=Components[Type].Torsions[index].D;

  posA=Cations[CurrentSystem][m].Atoms[A].Position;
  posB=Cations[CurrentSystem][m].Atoms[B].Position;
  posC=Cations[CurrentSystem][m].Atoms[C].Position;
  posD=Cations[CurrentSystem][m].Atoms[D].Position;

  Rab.x=posA.x-posB.x;
  Rab.y=posA.y-posB.y;
  Rab.z=posA.z-posB.z;

  Rbc.x=posB.x-posC.x;
  Rbc.y=posB.y-posC.y;
  Rbc.z=posB.z-posC.z;

  Rcd.x=posC.x-posD.x;
  Rcd.y=posC.y-posD.y;
  Rcd.z=posC.z-posD.z;

  rrbc=1.0/sqrt(Rbc.x*Rbc.x+Rbc.y*Rbc.y+Rbc.z*Rbc.z);

  // construct first dihedral vector
  Pb.x=Rab.y*Rbc.z-Rab.z*Rbc.y;
  Pb.y=Rab.z*Rbc.x-Rab.x*Rbc.z;
  Pb.z=Rab.x*Rbc.y-Rab.y*Rbc.x;
  pb2=Pb.x*Pb.x+Pb.y*Pb.y+Pb.z*Pb.z;
  rpb1=1.0/sqrt(pb2);

  // construct second dihedral vector
  Pc.x=Rbc.y*Rcd.z-Rbc.z*Rcd.y;
  Pc.y=Rbc.z*Rcd.x-Rbc.x*Rcd.z;
  Pc.z=Rbc.x*Rcd.y-Rbc.y*Rcd.x;
  pc2=Pc.x*Pc.x+Pc.y*Pc.y+Pc.z*Pc.z;
  rpc1=1.0/sqrt(pc2);

  // determine dihedral angle
  pbpc=Pb.x*Pc.x+Pb.y*Pc.y+Pb.z*Pc.z;
  cost=pbpc*rpb1*rpc1;
  sint=(Rbc.x*(Pc.y*Pb.z-Pc.z*Pb.y)+Rbc.y*(Pb.x*Pc.z-Pb.z*Pc.x)
        +Rbc.z*(Pc.x*Pb.y-Pc.y*Pb.x))*(rpb1*rpc1*rrbc);
  theta=atan2(sint,cost);
  if(theta<0.0) theta+=2.0*M_PI;
  return theta*180.0/M_PI;
}

int ReturnTorsionConformation(REAL theta)
{
  int conf;

  conf=0;
  if(((theta>=0.0)&&(theta<30.0))||((theta>=330.0)&&(theta<=360.0)))
    conf=SYNPERIPLANAR;
  else if ((theta>=30.0)&&(theta<90.0))
    conf=SYNCLINAL_PLUS;
  else if ((theta>=90.0)&&(theta<150.0))
    conf=ANTICLINAL_PLUS;
  else if((theta>=150.0)&(theta<210.0))
    conf=ANTIPERIPLANAR_PLUS;
  else if((theta>=210.0)&&(theta<270.0))
    conf=ANTICLINAL_MIN;
  else if((theta>=270.0)&&(theta<330.0))
    conf=SYNCLINAL_MIN;

  return conf;
}
