void WriteLammpsdata(void)
{
char buffer[256];
FILE *LMPFILE;
sprintf(buffer,"molecule_%d%s.data",count,FileNameAppend);
LMPPtr=fopen(LMPFILE,"w");
for (int i = 0; i < NumberOfSystems; i++){
   int numberofatoms = 0;
   int numberofbonds = 0;
   int numberofbends = 0;
   int numberoftorsions = 0;
   numberofatoms = Framework[i].TotalNumberOfAtoms;
  for (int j = 0; j < NumberOfComponents; j++){
     // calculate the number of atoms
     int component_molecule = Components[j].NumberOfMolecules[i]-(Components[j].FractionalMolecule[i]>=0?1:0)-Components[j].NumberOfRXMCMoleculesPresent[i];
     numberofatoms += Components[j].NumberOfAtoms * component_molecule;
     // calculate the number of bonds
     numberofbonds += Components[j].NumberOfBonds * component_molecule;
     // then calculate the number of bends (angles)
     numberofbends += Components[j].NumberOfBends * component_molecule;
     // then calculate the number of torsions (dihedral)
     numberoftorsions += Components[j].NumberOfTorsions * component_molecule;
  }  
  fprintf(LMPPtr, "%i  atoms\n", numberofatoms);
  fprintf(LMPPtr, "%i  bonds\n", numberofbonds);
  fprintf(LMPPtr, "%i  angles\n", numberofbends);
  fprintf(LMPPtr, "%i  dihedrals\n", numberoftorsions);
  // print out the size of the unit cell and the shifts in xyz directions
  fprintf(LMPPtr"0   -%18.12f xlo xhi\n", UnitCellBox[i].ax * NumberOfUnitCells[i].x);
  fprintf(LMPPtr"0   -%18.12f ylo yhi\n", UnitCellBox[i].by * NumberOfUnitCells[i].y);
  fprintf(LMPPtr"0   -%18.12f zlo zhi\n", UnitCellBox[i].cz * NumberOfUnitCells[i].z);
  fprintf(LMPPtr"0   -%18.12f -%18.12f -%18.12f xy xz yz\n", UnitCellBox[i].bx * NumberOfUnitCells[i].y, UnitCellBox[i].cx * NumberOfUnitCells[i].z, UnitCellBox[i].by * NumberOfUnitCells[i].z);
 int counter = 0;
 // print out the masses of the atoms
 for (int j = 0; j < NumberOfPseudoAtoms; j++){
    if (NumberOfPseudoAtomsType[i][j] > 0){
	counter += 1;
        Mass=PseudoAtoms[Framework[i].Atoms[0][i].Type].Mass;	
    }
 }
}
return 0;
}
