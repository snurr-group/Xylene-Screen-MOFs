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

bridge = pm_ase.AseAtomsAdaptor()
file_dir = os.getcwd() + '/'
mof_file = file_dir + '../' + 'UiO-66-Anthracene_clean_ML_chg.cif'
asemof = read(mof_file)
bridge.get_structure(asemof)
mof = bridge.get_structure(asemof)
nn_object = CrystalNN(x_diff_weight=0) #can optionally set search_cutoff=None
# for atomidx in range(len(mof)):
#     atom = mof[atomidx]
#     if atom.species_string == "C" : # If it is a carbon atom
#         local_env = nn_object.get_nn_info(mof,atomidx)
#         coord_num = len(local_env)
#         if coord_num == 3 and local_env[0]['site'].species_string == "O":
#             print(atomidx)
# Given a carbon atom, find the other end of the linker
# to go through a ring, it needs to pass through 2 carbon atoms that are 3-connected
# and both the carbon atoms are only connected to carbon
start = 152+48-1 # carbon index, at one end of the linker      
local_env = nn_object.get_nn_info(mof,start)
# go to the carbon atom connected 
for a in range(0, len(local_env)):
    if local_env[a]['site'].species_string == "C":
        step = local_env[a]['site_index']
        break
# Get local_env
local_env = nn_object.get_nn_info(mof,step)
# find the one carbon that is 3-connected and that all connected atoms are carbon
for a in range(0, len(local_env)):
    if (local_env[a]['site'].species_string == "C") and (not local_env[a]['site_index'] == start): # do not do the starting point
        new_env = nn_object.get_nn_info(mof,local_env[a]['site_index'])
        new_cord = len(new_env) # Get coord number of the local atoms
        if(new_cord == 3):
            # check if all connected atoms are carbon
            ALL_Carbon = True
            for b in range(0, new_cord):
                if(not new_env[b]['site'].species_string == "C"):
                    ALL_Carbon = False
            if(ALL_Carbon):
                steptwo = local_env[a]['site_index']
            break
# Find the next one
local_env = new_env
# find the one carbon that is 3-connected and that all connected atoms are carbon
for a in range(0, len(local_env)):
    if (local_env[a]['site'].species_string == "C") and (not local_env[a]['site_index'] == step): # do not do the starting point
        new_env = nn_object.get_nn_info(mof,local_env[a]['site_index'])
        new_cord = len(new_env) # Get coord number of the local atoms
        if(new_cord == 3):
            # check if all connected atoms are carbon
            ALL_Carbon = True
            for b in range(0, new_cord):
                if(not new_env[b]['site'].species_string == "C"):
                    ALL_Carbon = False
            if(ALL_Carbon):
                stepthree = local_env[a]['site_index']
            break
# Go to the next one
local_env = new_env
# find the one carbon that is 3-connected and that all connected atoms are carbon
for a in range(0, len(local_env)):
    if (local_env[a]['site'].species_string == "C") and (not local_env[a]['site_index'] == steptwo): # do not do the starting point
        new_env = nn_object.get_nn_info(mof,local_env[a]['site_index'])
        new_cord = len(new_env) # Get coord number of the local atoms
        if(new_cord == 3):
            # check if all connected atoms are carbon
            ALL_Carbon = True
            for b in range(0, new_cord):
                if(not new_env[b]['site'].species_string == "C"):
                    ALL_Carbon = False
            if(ALL_Carbon):
                stepfour = local_env[a]['site_index']
            break
# Go to the other end
local_env = new_env
# find the one carbon that is 3-connected and that all connected atoms are carbon
for a in range(0, len(local_env)):
    if (local_env[a]['site'].species_string == "C") and (not local_env[a]['site_index'] == stepthree): # do not do the starting point
        new_env = nn_object.get_nn_info(mof,local_env[a]['site_index'])
        new_cord = len(new_env) # Get coord number of the local atoms
        if(new_cord == 3):
            # check if all connected atoms are carbon
            ALL_Carbon = True
            for b in range(0, new_cord):
                if(not new_env[b]['site'].species_string == "C"):
                    ALL_Carbon = False
            if(not ALL_Carbon): # the other end should be connected to 2 oxygen atoms
                final = local_env[a]['site_index']
            break
# TO rotate, get all the atom index that are from this linker
linker_atoms = [start, final]
# start from one end
local_env_start = nn_object.get_nn_info(mof,start)
local_env_final = nn_object.get_nn_info(mof,final)
# go to the carbon atom connected 
for a in range(0, len(local_env_start)):
    if local_env_start[a]['site'].species_string == "C":
        linker_atoms.append(local_env_start[a]['site_index'])
        break
for a in range(0, len(local_env_final)):
    if local_env_final[a]['site'].species_string == "C":
        linker_atoms.append(local_env_final[a]['site_index'])
        break
# start searching

loop_list = [linker_atoms[3]]
cycle = 0
while ((len(loop_list) > 0) and (cycle < 5)):
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
                    loop_list.append(local_env[a]['site_index'])
        # as the loop finishes, remove the current atom from loop_list
        loop_list.remove(b)
    cycle += 1
########################################################
####################### Finally, Rotate#################
########################################################
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

# rotation_degrees = -60
# rotation_radians = np.radians(rotation_degrees)
newMOF = mof
print(mof[linker_atoms[5]])
mof.to(filename = 'Old_UiO-66-Anthracene.cif')
#for a in range(2,len(linker_atoms))
a = 2
end = linker_atoms[0]
end2 = linker_atoms[1]
linker_atoms.remove(end)
linker_atoms.remove(end2)
for point in linker_atoms:
    AB = [newMOF[end2].x - newMOF[point].x, newMOF[end2].y - newMOF[point].y, newMOF[end2].z - newMOF[point].z]
    AC = [newMOF[end2].x - newMOF[end].x, newMOF[end2].y - newMOF[end].y, newMOF[end2].z - newMOF[end].z]
    angle = np.arccos(np.dot(AB/np.linalg.norm(AB), AC/np.linalg.norm(AC)))
    project_line = AC/np.linalg.norm(AC)*(np.linalg.norm(AB)*np.sin(angle))
    Anchor_point = [newMOF[end2].x - project_line[0], newMOF[end2].y - project_line[1], newMOF[end2].z - project_line[2]]
    rotation_degrees = 180
    rotation_radians = np.radians(rotation_degrees)
    rotation_axis = np.asarray([mof[end].x-mof[end2].x, mof[end].y-mof[end2].y, mof[end].z-mof[end2].z])
    rotation_axis = rotation_axis/np.linalg.norm(rotation_axis)
    rotate_vec = np.asarray([newMOF[point].x - Anchor_point[0], newMOF[point].y - Anchor_point[1], newMOF[point].z - Anchor_point[2]])
    length = np.linalg.norm(rotate_vec)
    rotate_vec = rotate_vec/length
    rotated_vec = np.dot(rotation_matrix(rotation_axis,rotation_radians), rotate_vec)
    newpos = Anchor_point + rotated_vec*length
    newMOF[point].x = newpos[0]
    newMOF[point].y = newpos[1]
    newMOF[point].z = newpos[2]
    #newMOF.rotate_sites(indices = [point], theta = np.radians(60), axis = np.asarray([mof[end].x-mof[end2].x, mof[end].y-mof[end2].y, mof[end].z-mof[end2].z]), 
    #                    anchor=np.asarray(Anchor_point))
print(newMOF[linker_atoms[5]])

newMOF.to(filename = 'Rotated_UiO-66-Anthracene.cif')