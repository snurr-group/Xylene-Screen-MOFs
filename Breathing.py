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
bridge = pm_ase.AseAtomsAdaptor()
file_dir = os.getcwd() + '/'
mof_file = file_dir + 'v-bicycle_auto_remove_free_linker_P1.cif'
asemof = read(mof_file)
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
        if (coord_num == 3) and (attached_oxygen == 2):
            print(atomidx)
            linker_ends.append(atomidx)
    elif(atom.species_string == "V"): # If it is a metal atom
         node_metals.append(atomidx)
         
newMOF = mof
mof.to(filename = 'Old.cif')

# randomly select half of the linkers to rotate
all_linkers = linker_ends[:int(len(linker_ends)/2)] # assuming length is not singular
#https://www.geeksforgeeks.org/randomly-select-n-elements-from-list-in-python/
import random
counter = -1
#screw = [all_linkers[0]]
screw = all_linkers
all_linker_atoms = []
linker_list = []
for start in screw: # loop over ends
    counter += 1 
    print('Rotating Linker: ' + str(counter))
    # start is a carbon index, at one end of the linker      
    # TO rotate, get all the atom index that are from this linker
    linker_atoms = [start] # Start from the start, stop when it reaches the end
    # the end is another carbon atom that is connected to 2 oxygens
    # start from one end
    local_env_start = nn_object.get_nn_info(mof,start)
    # go to the carbon atom connected 
    for a in range(0, len(local_env_start)):
        if local_env_start[a]['site'].species_string == "C":
            linker_atoms.append(local_env_start[a]['site_index'])
            break
    
    # start searching
    loop_list = [linker_atoms[1]]
    cycle = 0
    while ((len(loop_list) > 0) and (cycle < 6)):
        for b in loop_list:
            local_env = nn_object.get_nn_info(mof, b)
            atom = mof[b]
            for a in range(0,len(local_env)):
                if not(local_env[a]['site_index'] in linker_atoms): # if it is a new atom
                    # calculate distance, check if it is within the cutoff
                    if(atom.distance(mof[local_env[a]['site_index']]) <= 1.7):
                    # add to the list
                    # also add to the loop list
                        linker_atoms.append(local_env[a]['site_index'])
                        # also check if the appended atom is the other terminal atom
                        # if so, we don't need to loop over the neighbors for the 
                        # terminal atoms 
                        if not(local_env[a]['site_index'] in linker_ends):
                            loop_list.append(local_env[a]['site_index'])
                        else:
                            end = local_env[a]['site_index']
            # as the loop finishes, remove the current atom from loop_list
            loop_list.remove(b)
        cycle += 1
    linker_atoms.sort(key = end.__eq__) # put the other end to the end of the list
    all_linker_atoms.extend(linker_atoms)
    linker_list.append(linker_atoms)
    # Finally, remove the linkers
##########################
# move all the nodes######
##########################
length = 10.58775
differ = 0.25 # fractional distance in y and z (y is not the same length as z)
Area = 0.5*differ*newMOF.lattice.b*differ*newMOF.lattice.c# area of the triangle 
# change ratio in Z axis
ratio = 1.2
# target Z length
target_c = newMOF.lattice.c*ratio # shrinks in z-axis
newdistance_z = differ*target_c
newdistance_y = np.sqrt(length**2-newdistance_z**2)
target_b = newdistance_y/differ # expand proportionally in y-axis
# make a new lattice
newla = pm.core.structure.Structure([[newMOF.lattice.a,0,0], [0,target_b,0], [0,0,target_c]], [], [])
# append all atoms to the new structure
for atomidx in range(len(mof)):
    newla.append(mof[atomidx].species_string, mof[atomidx].coords, True)
atom_moved = []
#node_metals = [440]
for start in node_metals:
    if(start in atom_moved):
        continue
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
                        if not((local_env[a]['site_index'] in linker_ends) or (local_env[a]['site_index'] in all_linker_atoms)):
                            loop_list.append(local_env[a]['site_index'])
                            node_atoms.append(local_env[a]['site_index'])
            # as the loop finishes, remove the current atom from loop_list
            loop_list.remove(b)
        cycle += 1
        # move the linker to the newly scaled location
        # get the new location of the first atom
        new_first_location = mof[node_atoms[0]].frac_coords*[mof.lattice.a, target_b, target_c]
    for a in range(1,len(node_atoms)): # skip the first atom, it serves as an anchor
        old_vector = pm.util.coord.pbc_diff(newMOF[node_atoms[0]].frac_coords, newMOF[node_atoms[a]].frac_coords)*[newMOF.lattice.a, newMOF.lattice.b, newMOF.lattice.c]
        #old_vector = newMOF[node_atoms[0]].coords - newMOF[node_atoms[a]].coords # turn off pbc
        new_pos = new_first_location - old_vector 
        newla[node_atoms[a]].coords = new_pos
        print(node_atoms[a], new_pos)
        #newl.append(newMOF[node_atoms[a]].species_string, new_pos, True) # set "True", telling it that coords are cartesian
    # finally, move the position of the first atom
    newla[node_atoms[0]].coords = new_first_location
    #newl.append(newMOF[node_atoms[0]].species_string, new_first_location, True)
    atom_moved.extend(node_atoms)
# remove all at once!

#newla.to(filename = 'latest_just_nodes.cif')
###########################
## Move back the linkers###
###########################
#linker_list = [linker_list[0],linker_list[5]]
#linker_list = [linker_list[0]]
for linker in linker_list:
    start = linker[0]
    end = linker[-1]
    start_environment = nn_object.get_nn_info(mof, start)
    end_environment = nn_object.get_nn_info(mof, end)
    # get center of the two oxygens at the start
    pos_start =  np.array([0.0,0.0,0.0])
    pos_start_new = np.array([0.0,0.0,0.0]) # the compressed positions
    # also get the rotation axis: difference between the two oxygens
    counter = 0
    for a in start_environment:
        if a['site'].species_string == 'O':
            pos_start += 0.5*a['site'].coords
            pos_start_new += 0.5*newla[a['site_index']].coords
            if(counter == 0):
                rotation_point = [a['site_index']]
            else:
                rotation_point.append(a['site_index'])
            counter+=1
    rotation_axis = pm.util.coord.pbc_diff(newla[rotation_point[0]].frac_coords, newla[rotation_point[1]].frac_coords)*[newla.lattice.a, newla.lattice.b, newla.lattice.c]
    rotation_axis = rotation_axis/np.linalg.norm(rotation_axis)
    print(rotation_point)
    print(rotation_axis)
    # get center of the two oxygens at the end
    pos_end =  np.array([0.0,0.0,0.0])
    pos_end_new =  np.array([0.0,0.0,0.0])
    for a in end_environment:
        if a['site'].species_string == 'O':
            pos_end += 0.5*a['site'].coords
            pos_end_new += 0.5*newla[a['site_index']].coords
    linker_orien = pm.util.coord.pbc_diff(pos_start/[mof.lattice.a, mof.lattice.b, mof.lattice.c], 
                                          pos_end/[mof.lattice.a, mof.lattice.b, mof.lattice.c])*[mof.lattice.a, mof.lattice.b, mof.lattice.c]
    linker_orien_new = pm.util.coord.pbc_diff(pos_start_new/[newla.lattice.a, newla.lattice.b, newla.lattice.c], 
                                          pos_end_new/[newla.lattice.a, newla.lattice.b, newla.lattice.c])*[newla.lattice.a, newla.lattice.b, newla.lattice.c]
    angle = np.arccos(np.dot(linker_orien/np.linalg.norm(linker_orien),
                                 linker_orien_new/np.linalg.norm(linker_orien_new)))
    # first translate
    distance = pos_start_new-pos_start
    for atom in linker:
        new_pos = newla[atom].coords + distance
        newla[atom].coords = new_pos
    # then rotate
    # since rotating angle or (-angle) is arbitrary, we can do a buffer
    # and at the end, do a distance check
    # if the check is failed, rotate -angle
    temp_pos = [] # list
    for point in linker:
        AB = pm.util.coord.pbc_diff(newla[rotation_point[1]].frac_coords, newla[point].frac_coords)*[newla.lattice.a, newla.lattice.b, newla.lattice.c]
        
        AC = pm.util.coord.pbc_diff(newla[rotation_point[1]].frac_coords, newla[rotation_point[0]].frac_coords)*[newla.lattice.a, newla.lattice.b, newla.lattice.c] # this is the rotation axis
        #print(newMOF[point].species_string)
        #print("AB: " + str(AB[0]), str(AB[1]), str(AB[2]))
        #print("AC: " + str(AC[0]), str(AC[1]), str(AC[2]))
        axisangle = np.arccos(np.dot(AB/np.linalg.norm(AB),
                                 AC/np.linalg.norm(AC)))
        project_line = AC/np.linalg.norm(AC)*(np.linalg.norm(AB)*np.sin(axisangle))
        #print(project_line)
        Anchor_point = np.asarray([newla[rotation_point[1]].x - project_line[0], newla[rotation_point[1]].y - project_line[1], newla[rotation_point[1]].z - project_line[2]])
        #print(Anchor_point)
        rotation_radians = angle
        #rotation_axis = pm.util.coord.pbc_diff(newMOF[end].frac_coords, newMOF[end2].frac_coords)*[newMOF.lattice.a, newMOF.lattice.b, newMOF.lattice.c]
        rotation_axis = AC/np.linalg.norm(AC)
        rotate_vec = pm.util.coord.pbc_diff(newla[point].frac_coords, Anchor_point/[newla.lattice.a, newla.lattice.b, newla.lattice.c])*[newla.lattice.a, newla.lattice.b, newla.lattice.c]
        
        coord1 = [[newla[point].x, newla[point].y, newla[point].z]]
        coord2 = [Anchor_point]
        length = float(pm.util.coord.all_distances(coord1, coord2))
        rotate_vec = rotate_vec/length
        rotated_vec = np.dot(rotation_matrix(rotation_axis,rotation_radians), rotate_vec)
        newpos = Anchor_point + rotated_vec*length
        #print('newpos:' + str(newpos))
        # newla[point].x = newpos[0]
        # newla[point].y = newpos[1]
        # newla[point].z = newpos[2]
        temp_pos.append(newpos)
    # before finally recording the locations, do a distance check
    # check the distance at the end of the linker with the mid-point 
    # of the two oxygens that is near it
    distance = np.linalg.norm(pm.util.coord.pbc_diff(newpos/[newla.lattice.a, newla.lattice.b, newla.lattice.c], pos_end_new/[newla.lattice.a, newla.lattice.b, newla.lattice.c])*[newla.lattice.a, newla.lattice.b, newla.lattice.c])
    if(distance < 3):
        counter = 0
        for point in linker:
            newla[point].coords = temp_pos[counter]
            counter += 1
    else: ### rotate a negative angle, can be optimized though ###
        for point in linker:
            AB = pm.util.coord.pbc_diff(newla[rotation_point[1]].frac_coords, newla[point].frac_coords)*[newla.lattice.a, newla.lattice.b, newla.lattice.c]
            
            AC = pm.util.coord.pbc_diff(newla[rotation_point[1]].frac_coords, newla[rotation_point[0]].frac_coords)*[newla.lattice.a, newla.lattice.b, newla.lattice.c] # this is the rotation axis
            #print(newMOF[point].species_string)
            #print("AB: " + str(AB[0]), str(AB[1]), str(AB[2]))
            #print("AC: " + str(AC[0]), str(AC[1]), str(AC[2]))
            axisangle = np.arccos(np.dot(AB/np.linalg.norm(AB),
                                     AC/np.linalg.norm(AC)))
            project_line = AC/np.linalg.norm(AC)*(np.linalg.norm(AB)*np.sin(axisangle))
            #print(project_line)
            Anchor_point = np.asarray([newla[rotation_point[1]].x - project_line[0], newla[rotation_point[1]].y - project_line[1], newla[rotation_point[1]].z - project_line[2]])
            #print(Anchor_point)
            rotation_radians = -angle
            #rotation_axis = pm.util.coord.pbc_diff(newMOF[end].frac_coords, newMOF[end2].frac_coords)*[newMOF.lattice.a, newMOF.lattice.b, newMOF.lattice.c]
            rotation_axis = AC/np.linalg.norm(AC)
            rotate_vec = pm.util.coord.pbc_diff(newla[point].frac_coords, Anchor_point/[newla.lattice.a, newla.lattice.b, newla.lattice.c])*[newla.lattice.a, newla.lattice.b, newla.lattice.c]
            
            coord1 = [[newla[point].x, newla[point].y, newla[point].z]]
            coord2 = [Anchor_point]
            length = float(pm.util.coord.all_distances(coord1, coord2))
            rotate_vec = rotate_vec/length
            rotated_vec = np.dot(rotation_matrix(rotation_axis,rotation_radians), rotate_vec)
            newpos = Anchor_point + rotated_vec*length
            newla[point].coords = newpos
newla.to(filename = 'latest_moved_' + str(ratio) + '.cif')