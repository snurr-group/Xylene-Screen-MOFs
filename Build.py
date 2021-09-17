import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import os
import os.path as path
import re

from itertools import repeat
from ase.io import read
import pymatgen as pm
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.io import ase as pm_ase
import math
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                      [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                      [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])
#ratios = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 
#          1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5]
ratios = [1.5] # define deformation ratio here
bridge = pm_ase.AseAtomsAdaptor()
file_dir = os.getcwd() + '/'
node_file = file_dir + 'JustNode.cif'
asemof = read(node_file)
bridge.get_structure(asemof)
mof = bridge.get_structure(asemof)
nn_object = CrystalNN(x_diff_weight=0) #can optionally set search_cutoff=None
linker_ends = []
node_metals= []
for atomidx in range(len(mof)):
    atom = mof[atomidx]
    if atom.species_string == "C" : # If it is a carbon atom
        local_env = nn_object.get_nn_info(mof,atomidx)
        coord_num = len(local_env)
        attached_oxygen = 0
        for a in range(0, coord_num):
            if (local_env[a]['site'].species_string == "O"):
                attached_oxygen += 1
        if attached_oxygen == 2:
            print(atomidx)
            linker_ends.append(atomidx)
    elif(atom.species_string == "V"): # If it is a metal atom
         node_metals.append(atomidx)
###############################################################
# before doing anything,                                     ##
# We need to scale the box to fit the size of the new linker ##
###############################################################

# need to calculate distance of center of two Vanadium atom to the linker end
# sqrt(((0.4394-0.5)*31.5657 )^2+((0.75-0.6944)*28.235)^2), using JustNode.cif
yz_distance = 2.47459
ratio = (5.63+yz_distance*2)/(5.647+yz_distance*2)
target_b = mof.lattice.b*ratio
target_c = mof.lattice.c*ratio
newla = pm.core.structure.Structure([[mof.lattice.a,0,0], [0,target_b,0], [0,0,target_c]], [], [])
for atomidx in range(len(mof)):
    newla.append(mof[atomidx].species_string, mof[atomidx].coords, True)
#######################################################         
# First, identify the nodes and create a node list#####
#######################################################
counter = -1
atom_moved = []
all_node_atoms = []
node_list = []
for start in node_metals:
    if(start in atom_moved):
        continue
    counter += 1
    print("Node: " + str(counter))
    node_atoms = [start]
    loop_list = [start]
    cycle = 0
    ###############################
    # identify all the node atoms##
    ###############################
    while ((len(loop_list) > 0) and (cycle < 6)):
        for b in loop_list:
            local_env = nn_object.get_nn_info(mof, b)
            atom = mof[b]
            for a in range(0,len(local_env)):
                if not(local_env[a]['site_index'] in node_atoms): # if it is a new atom
                    # calculate distance, check if it is within the cutoff
                    if(atom.distance(mof[local_env[a]['site_index']]) <= 4):
                    # add to the list
                    # also add to the loop list
                        # also check if the appended atom is the other terminal atom
                        # if so, we don't need to loop over the neighbors for the 
                        # terminal atoms 
                        ### here, we don't need to exclude the linker_end atoms
                        #if not(local_env[a]['site_index'] in linker_ends):
                        loop_list.append(local_env[a]['site_index'])
                        node_atoms.append(local_env[a]['site_index'])
            loop_list.remove(b)
    all_node_atoms.extend(node_atoms)
    node_list.append(node_atoms)
    # scale the positions
    new_first_location = mof[node_atoms[0]].frac_coords*[mof.lattice.a, target_b, target_c]
    for a in range(1,len(node_atoms)): # skip the first atom, it serves as an anchor
        old_vector = pm.util.coord.pbc_diff(mof[node_atoms[0]].frac_coords, mof[node_atoms[a]].frac_coords)*[mof.lattice.a, mof.lattice.b, mof.lattice.c]
        #old_vector = mof[node_atoms[0]].coords - mof[node_atoms[a]].coords # turn off pbc
        new_pos = new_first_location - old_vector 
        newla[node_atoms[a]].coords = new_pos
        print(node_atoms[a], new_pos)
        #newl.append(mof[node_atoms[a]].species_string, new_pos, True) # set "True", telling it that coords are cartesian
    # finally, move the position of the first atom
    newla[node_atoms[0]].coords = new_first_location
newla.to(filename = 'ScaledNode.cif')
##################################    
# Then, identify the end pairs ###
##################################
end_pairs = []
paired = []
counter = -1
# find the corresponding carbons
for end in linker_ends:
    if(end in paired):
        paired.append(end)
        continue
    counter += 1
    print("Looping pair: " + str(counter) + ', end is: ' + str(end))
    other_ends = [x for x in linker_ends if x not in paired]
    for other in other_ends:
        # calculate distance between each pair
        dist = pm.util.coord.all_distances([mof[end].coords], [mof[other].coords])
        if((dist<6) and (dist > 5.6)):
            end_pairs.append([end, other])
            paired.append(other)
            break
#######################################
# Read and Rotate the desired Linker###
#######################################
# Read the linker
Linker_file = file_dir + 'JustCubane_real2.cif'
aselinker = read(Linker_file)
bridge.get_structure(aselinker)
newLinker = bridge.get_structure(aselinker)
linker_joint = 8 # 0th atom in the linker
other_linker_end = 6 # 12th in the linker
# Loop over the linker pairs
#end_pairs = [end_pairs[0]]
for ends in end_pairs:
    joint = ends[0]
    addedID = []
    # Translate the linker to attach it on the joint
    joint_distance = pm.util.coord.pbc_diff(newLinker[linker_joint].coords/[newla.lattice.a, newla.lattice.b, newla.lattice.c], newla[joint].coords/[newla.lattice.a, newla.lattice.b, newla.lattice.c])
    # Get the old distances for atoms in the linker
    for a in range(0,len(newLinker)):
        if(a == linker_joint):
            continue
        old_vector = pm.util.coord.pbc_diff(newLinker[a].frac_coords, newLinker[linker_joint].frac_coords)*[newLinker.lattice.a, newLinker.lattice.b, newLinker.lattice.c]
        new_pos = newla[joint].coords + old_vector 
        newla.append(newLinker[a].species_string, new_pos, True)
        if(a == other_linker_end):
            oldaxis = pm.util.coord.pbc_diff(newla[len(newla)-1].frac_coords, newla[joint].frac_coords)*[newla.lattice.a, newla.lattice.b, newla.lattice.c]
            oldaxis = oldaxis/np.linalg.norm(oldaxis)
        addedID.append(len(newla)-1)
    #Rotate
    axis = pm.util.coord.pbc_diff(newla[ends[1]].frac_coords, newla[joint].frac_coords)*[newla.lattice.a, newla.lattice.b, newla.lattice.c]
    axis = axis/np.linalg.norm(axis)
    # get the angle between the axes
    axisangle = np.arccos(np.dot(axis,oldaxis))
    rotateaxis = np.cross(axis, oldaxis)/np.linalg.norm(np.cross(axis, oldaxis))
    newla.rotate_sites(addedID, (3.14*2-axisangle), rotateaxis, newla[joint].coords)
# save the structure
newla.to(filename = 'MIL-V-Cubane.cif')