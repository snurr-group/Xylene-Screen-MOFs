# Three-to-Full-Xylene
Script that takes in text restart file using three-site model xylene and converts it to text restart file of full model xylene
It requires RASPA.
The full model it generates is using OPLS

**Requirement**
1. Restart file using the three-site model of adsorption of xylene (or isomer mixture) in a structure (text restart file)
2. Pdb movie file of the framework, it has a general name: "Movies/System_0/Framework_final.pdb"
3. Need to check the directories
**Logic**
1. Reads restart file of three-site model
2. For the mX and oX, since the three-site model already defined the plane, rotate the bonds and creates the carbon and hydrogen atoms
3. For pX molecules, randomly select an angle and then create the molecule
4. calculate overlaps and optimize within the code (optional)
5. Generates full-model restart file (text file)
******************************************************************************************************************************************
