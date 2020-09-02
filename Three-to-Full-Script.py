# change in file: change from counting overlap (hard sphere) to LJ energies
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import os
import os.path as path
import re

from itertools import repeat

temperature = 433
three_restart_dir = "F:/Grad School/Xylene_Project/1D Channels/Prelim_results/All_3_site_restart/"
# EBMC or NO-EB: determines the final restart file
EBMC = False
# if we want partial charges on the adsorbate atoms
charges = True
# list all files in the folders
restarts_dirs = os.listdir(three_restart_dir)
for filenumber in range(0, len(restarts_dirs)): #len(restarts_dirs)
    input_dir = three_restart_dir + restarts_dirs[filenumber] + '/'
    restart_dir = input_dir + 'Restart/System_0/'
    # extra step: get the restart file name
    filename = os.listdir(restart_dir)[0] # assuming there is only 1 restart file
    output_dir="Results/" + restarts_dirs[filenumber] + '/'
    if not path.exists(output_dir):
        os.mkdir(output_dir)
    file = restart_dir + filename
    abc = open(file)
    lines = abc.read()
    df = pd.DataFrame(columns = ['Num', 'Type', 'X', 'Y', 'Z'])
    a=0
    # header
    
    with open(file) as f:
        for line in f:       
            if("Adsorbate-atom-position:" in line):
                a+=1
                spline = line.split()
                num_molec = int(spline[1])
                type = int(spline[2])
                x = float(spline[3])
                y = float(spline[4])
                z = float(spline[5])
                pos = [num_molec, type, x, y, z]
                df.loc[a] = pos
                
    header = []
    with open(file) as f:
        for line in f:
            # get the headers
            header.append(line)
            if("Reactions:" in line):
                break
    # eliminate the fractional molecules
    a=0
    Eliminate_Scale = False
    Mol_Scale = []
    with open(file) as f:
        for line in f:
            if("Adsorbate-atom-scaling:" in line):
                a+=1
                spline = line.split()
                num_molec = int(spline[1])
                scale = float(spline[3])
                if(scale < 1.0000 and Eliminate_Scale):
                    df = df.drop(a)
                else:
                    Mol_Scale.append(scale)
    
    df = df.sort_values(by = ['Num', 'Type'])
    df.reset_index(drop = True, inplace = True)
    num_particles = int(df.shape[0])
    # identify xylene isomer
    def test_xy_type(vec_origin, vec_1, vec_2):
        a = vec_origin[2:5]
        b = vec_1[2:5]
        c = vec_2[2:5]
        Dist_1O = a-b
        Dist_2O = a-c
        dotted = np.dot(Dist_1O/np.linalg.norm(Dist_1O), Dist_2O/np.linalg.norm(Dist_2O))
        if (dotted <= -0.98):
            mol_type = 1 #px
        elif ((dotted < 0.6) and (dotted > 0.4)):
            mol_type=2
        elif ((dotted >= -0.6) and (dotted <= -0.4)):
            mol_type = 3
        else:
            print("We are in BIG Trouble!")
            mol_type = 4
        return mol_type
    # grow mx
    length_CH3_S=2.91
    length_S_H=2.48
    length_S_C=1.4
    # Euler-Rodrigues Formula
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
    
    
    # need to modify the ranking of the atoms, Sp, then carbon, then CH3 and hydro
    def grow_mx(a, vec_origin, vec_1, vec_2):
        full_pos.loc[13*a] = vec_origin
        vec_origin = vec_origin[2:5]
        temp_vec1 = vec_1
        temp_vec2 = vec_2
        vec_1 = vec_1[2:5]
        vec_2 = vec_2[2:5]
        Dist_1O = vec_1-vec_origin
        Dist_2O = vec_2-vec_origin
        Norm = np.cross(Dist_1O/np.linalg.norm(Dist_1O), Dist_2O/np.linalg.norm(Dist_2O))
        rotation_degrees = 60
        rotation_radians = np.radians(rotation_degrees)
        rotation_axis = np.array(Norm)
        rotate_vec = Dist_1O/np.linalg.norm(Dist_1O)
        rotated_vec = np.dot(rotation_matrix(rotation_axis,rotation_radians), rotate_vec)
        # check for overlap
        methyl_dist = rotated_vec - Dist_2O/np.linalg.norm(Dist_2O)
        if(np.linalg.norm(methyl_dist) < 1e-2):
            print("overlap")
            rotation_degrees = -60
            rotation_radians = np.radians(rotation_degrees)
            rotation_axis = np.array(Norm)
            rotate_vec = Dist_1O/np.linalg.norm(Dist_1O)
            rotated_vec = np.dot(rotation_matrix(rotation_axis,rotation_radians), rotate_vec)
        # grow the carbons
        # 5 and 6 are those in line with methyl groups
        carbon1 = [a, 1, Dist_1O[0]*length_S_C/length_CH3_S+vec_origin[0],
            Dist_1O[1]*length_S_C/length_CH3_S+vec_origin[1],
            Dist_1O[2]*length_S_C/length_CH3_S+vec_origin[2]]
        
        carbon2 = [a, 2, Dist_2O[0]*length_S_C/length_CH3_S+vec_origin[0],
            Dist_2O[1]*length_S_C/length_CH3_S+vec_origin[1],
            Dist_2O[2]*length_S_C/length_CH3_S+vec_origin[2]]
        
        carbon3 = [a, 3, -Dist_1O[0]*length_S_C/length_CH3_S+vec_origin[0],
            -Dist_1O[1]*length_S_C/length_CH3_S+vec_origin[1],
            -Dist_1O[2]*length_S_C/length_CH3_S+vec_origin[2]]
        
        carbon4 = [a, 4, -Dist_2O[0]*length_S_C/length_CH3_S+vec_origin[0],
            -Dist_2O[1]*length_S_C/length_CH3_S+vec_origin[1],
            -Dist_2O[2]*length_S_C/length_CH3_S+vec_origin[2]]
    
        carbon5 = [a, 5, length_S_C*rotated_vec[0]+vec_origin[0],
            length_S_C*rotated_vec[1]+vec_origin[1],
            length_S_C*rotated_vec[2]+vec_origin[2]]
        
        carbon6 = [a, 6, -length_S_C*rotated_vec[0]+vec_origin[0],
            -length_S_C*rotated_vec[1]+vec_origin[1],
            -length_S_C*rotated_vec[2]+vec_origin[2]]
    
        # then grow the hydrogen atoms
        hydro1 = [a, 9, -Dist_1O[0]*length_S_H/length_CH3_S+vec_origin[0],
            -Dist_1O[1]*length_S_H/length_CH3_S+vec_origin[1],
            -Dist_1O[2]*length_S_H/length_CH3_S+vec_origin[2]]
        
        hydro2 = [a, 10, -Dist_2O[0]*length_S_H/length_CH3_S+vec_origin[0],
            -Dist_2O[1]*length_S_H/length_CH3_S+vec_origin[1],
            -Dist_2O[2]*length_S_H/length_CH3_S+vec_origin[2]]
        
        hydro3 = [a, 11, length_S_H*rotated_vec[0]+vec_origin[0],
            length_S_H*rotated_vec[1]+vec_origin[1],
            length_S_H*rotated_vec[2]+vec_origin[2]]
        
        hydro4 = [a, 12, -length_S_H*rotated_vec[0]+vec_origin[0],
            -length_S_H*rotated_vec[1]+vec_origin[1],
            -length_S_H*rotated_vec[2]+vec_origin[2]]
        
        full_pos.loc[13*a+1] = carbon1
        full_pos.loc[13*a+2] = carbon2
        full_pos.loc[13*a+3] = carbon3
        full_pos.loc[13*a+4] = carbon4
        full_pos.loc[13*a+5] = carbon5
        full_pos.loc[13*a+6] = carbon6
        full_pos.loc[13*a+7] = temp_vec1
        full_pos.loc[13*a+8] = temp_vec2
        full_pos.loc[13*a+9] = hydro1
        full_pos.loc[13*a+10] = hydro2
        full_pos.loc[13*a+11] = hydro3
        full_pos.loc[13*a+12] = hydro4
    
    def grow_ox(a, vec_origin, vec_1, vec_2):
        full_pos.loc[13*a] = vec_origin
        temp_vec1 = vec_1
        temp_vec2 = vec_2
        vec_origin = vec_origin[2:5]
        vec_1 = vec_1[2:5]
        vec_2 = vec_2[2:5]
        Dist_1O = vec_1-vec_origin
        Dist_2O = vec_2-vec_origin
        Norm = np.cross(Dist_1O/np.linalg.norm(Dist_1O), Dist_2O/np.linalg.norm(Dist_2O))
        rotation_degrees = 60
        rotation_radians = np.radians(rotation_degrees)
        rotation_axis = np.array(Norm)
        rotate_vec = Dist_1O/np.linalg.norm(Dist_1O)
        rotated_vec = np.dot(rotation_matrix(rotation_axis,rotation_radians), rotate_vec)
        # check for overlap
        methyl_dist = rotated_vec - Dist_2O/np.linalg.norm(Dist_2O)
        if(np.linalg.norm(methyl_dist) < 1e-2):
            print("overlap")
            rotation_degrees = -60
            rotation_radians = np.radians(rotation_degrees)
            rotation_axis = np.array(Norm)
            rotate_vec = Dist_1O/np.linalg.norm(Dist_1O)
            rotated_vec = np.dot(rotation_matrix(rotation_axis,rotation_radians), rotate_vec)
        # grow the carbons
        # 5 and 6 are those in line with methyl groups
        carbon1 = [a, 1, Dist_1O[0]*length_S_C/length_CH3_S+vec_origin[0],
            Dist_1O[1]*length_S_C/length_CH3_S+vec_origin[1],
            Dist_1O[2]*length_S_C/length_CH3_S+vec_origin[2]]
        
        carbon2 = [a, 2, Dist_2O[0]*length_S_C/length_CH3_S+vec_origin[0],
            Dist_2O[1]*length_S_C/length_CH3_S+vec_origin[1],
            Dist_2O[2]*length_S_C/length_CH3_S+vec_origin[2]]
        
        carbon3 = [a, 3, -Dist_1O[0]*length_S_C/length_CH3_S+vec_origin[0],
            -Dist_1O[1]*length_S_C/length_CH3_S+vec_origin[1],
            -Dist_1O[2]*length_S_C/length_CH3_S+vec_origin[2]]
        
        carbon4 = [a, 4, -Dist_2O[0]*length_S_C/length_CH3_S+vec_origin[0],
            -Dist_2O[1]*length_S_C/length_CH3_S+vec_origin[1],
            -Dist_2O[2]*length_S_C/length_CH3_S+vec_origin[2]]
    
        carbon5 = [a, 5, length_S_C*rotated_vec[0]+vec_origin[0],
            length_S_C*rotated_vec[1]+vec_origin[1],
            length_S_C*rotated_vec[2]+vec_origin[2]]
        
        carbon6 = [a, 6, -length_S_C*rotated_vec[0]+vec_origin[0],
            -length_S_C*rotated_vec[1]+vec_origin[1],
            -length_S_C*rotated_vec[2]+vec_origin[2]]
    
        # then grow the hydrogen atoms
        hydro1 = [a, 9, -Dist_1O[0]*length_S_H/length_CH3_S+vec_origin[0],
            -Dist_1O[1]*length_S_H/length_CH3_S+vec_origin[1],
            -Dist_1O[2]*length_S_H/length_CH3_S+vec_origin[2]]
        
        hydro2 = [a, 10, -Dist_2O[0]*length_S_H/length_CH3_S+vec_origin[0],
            -Dist_2O[1]*length_S_H/length_CH3_S+vec_origin[1],
            -Dist_2O[2]*length_S_H/length_CH3_S+vec_origin[2]]
        
        hydro3 = [a, 11, length_S_H*rotated_vec[0]+vec_origin[0],
            length_S_H*rotated_vec[1]+vec_origin[1],
            length_S_H*rotated_vec[2]+vec_origin[2]]
        
        hydro4 = [a, 12, -length_S_H*rotated_vec[0]+vec_origin[0],
            -length_S_H*rotated_vec[1]+vec_origin[1],
            -length_S_H*rotated_vec[2]+vec_origin[2]]
        
        full_pos.loc[13*a+1] = carbon1
        full_pos.loc[13*a+2] = carbon2
        full_pos.loc[13*a+3] = carbon3
        full_pos.loc[13*a+4] = carbon4
        full_pos.loc[13*a+5] = carbon5
        full_pos.loc[13*a+6] = carbon6
        full_pos.loc[13*a+7] = temp_vec1
        full_pos.loc[13*a+8] = temp_vec2
        full_pos.loc[13*a+9] = hydro1
        full_pos.loc[13*a+10] = hydro2
        full_pos.loc[13*a+11] = hydro3
        full_pos.loc[13*a+12] = hydro4
    
    # grow px molecules
    from random import random
    circle_radius = length_S_H*np.sqrt(3)/2
    length_S_circle = 0.5 * length_S_H
    def grow_px(a, vec_origin, vec_1, vec_2):
        full_pos.loc[13*a] = vec_origin
        temp_vec1 = vec_1
        temp_vec2 = vec_2
        vec_1 = vec_1[2:5]
        vec_2 = vec_2[2:5]
        vec_origin = vec_origin[2:5]
        long_vect = vec_2 - vec_1
        mod_origin_1 = vec_origin + (vec_1-vec_origin)/(length_CH3_S)*length_S_circle
        theta = random()*2*np.pi
        pointtt = np.array([circle_radius*np.cos(theta), circle_radius*np.sin(theta), 0])
        #now, we need the cross of long-vector of the molecule and z-axis
        normal_rotate = np.cross(long_vect/(2*length_CH3_S), [0, 0, 1]) 
        angle_vect_z = math.atan2(np.linalg.norm(np.cross(long_vect/(2*length_CH3_S),[0, 0, 1])),
                                  np.dot(long_vect/(2*length_CH3_S),[0, 0, 1]))
        rotated_vec = np.dot(rotation_matrix(normal_rotate,-angle_vect_z), pointtt)
        hydro1 = np.append([a, 9], rotated_vec + mod_origin_1)
        hydro2 = np.append([a, 10], mod_origin_1-(hydro1[2:5] - mod_origin_1))
        hydro3 = np.append([a, 11], vec_origin - (hydro1[2:5] - vec_origin))
        hydro4 = np.append([a, 12], vec_origin - (hydro2[2:5] - vec_origin))
        # now, get the carbons
        carbon1 = np.append([a, 1], (hydro1[2:5]-vec_origin)/length_S_H*length_S_C+vec_origin)
        carbon2 = np.append([a, 2], -(hydro1[2:5]-vec_origin)/length_S_H*length_S_C+vec_origin)
        carbon3 = np.append([a, 3], (hydro2[2:5]-vec_origin)/length_S_H*length_S_C+vec_origin)
        carbon4 = np.append([a, 4], -(hydro2[2:5]-vec_origin)/length_S_H*length_S_C+vec_origin)
        carbon5 = np.append([a, 5], (vec_1-vec_origin)/length_CH3_S*length_S_C+vec_origin)
        carbon6 = np.append([a, 6], (vec_2-vec_origin)/length_CH3_S*length_S_C+vec_origin)
        full_pos.loc[13*a+1] = carbon1
        full_pos.loc[13*a+2] = carbon2
        full_pos.loc[13*a+3] = carbon3
        full_pos.loc[13*a+4] = carbon4
        full_pos.loc[13*a+5] = carbon5
        full_pos.loc[13*a+6] = carbon6
        full_pos.loc[13*a+7] = temp_vec1
        full_pos.loc[13*a+8] = temp_vec2
        full_pos.loc[13*a+9] = hydro1
        full_pos.loc[13*a+10] = hydro2
        full_pos.loc[13*a+11] = hydro3
        full_pos.loc[13*a+12] = hydro4
    
    def optimize_px(vec_origin, vec_1, vec_2):
        temp.loc[0] = vec_origin
        temp_vec1 = vec_1
        temp_vec2 = vec_2
        vec_1 = vec_1[2:5]
        vec_2 = vec_2[2:5]
        vec_origin = vec_origin[2:5]
        long_vect = vec_2 - vec_1
        mod_origin_1 = vec_origin + (vec_1-vec_origin)/(length_CH3_S)*length_S_circle
        theta = random()*2*np.pi
        pointtt = np.array([circle_radius*np.cos(theta), circle_radius*np.sin(theta), 0])
        #now, we need the cross of long-vector of the molecule and z-axis
        normal_rotate = np.cross(long_vect/(2*length_CH3_S), [0, 0, 1]) 
        angle_vect_z = math.atan2(np.linalg.norm(np.cross(long_vect/(2*length_CH3_S),[0, 0, 1])),
                                 np.dot(long_vect/(2*length_CH3_S),[0, 0, 1]))
        rotated_vec = np.dot(rotation_matrix(normal_rotate,-angle_vect_z), pointtt)
        hydro1 = np.append([a, 9], rotated_vec + mod_origin_1)
        hydro2 = np.append([a, 10], mod_origin_1-(hydro1[2:5] - mod_origin_1))
        hydro3 = np.append([a, 11], vec_origin - (hydro1[2:5] - vec_origin))
        hydro4 = np.append([a, 12], vec_origin - (hydro2[2:5] - vec_origin))
        # now, get the carbons
        carbon1 = np.append([a, 1], (hydro1[2:5]-vec_origin)/length_S_H*length_S_C+vec_origin)
        carbon2 = np.append([a, 2], -(hydro1[2:5]-vec_origin)/length_S_H*length_S_C+vec_origin)
        carbon3 = np.append([a, 3], (hydro2[2:5]-vec_origin)/length_S_H*length_S_C+vec_origin)
        carbon4 = np.append([a, 4], -(hydro2[2:5]-vec_origin)/length_S_H*length_S_C+vec_origin)
        carbon5 = np.append([a, 5], (vec_1-vec_origin)/length_CH3_S*length_S_C+vec_origin)
        carbon6 = np.append([a, 6], (vec_2-vec_origin)/length_CH3_S*length_S_C+vec_origin)
        temp.loc[1] = carbon1
        temp.loc[2] = carbon2
        temp.loc[3] = carbon3
        temp.loc[4] = carbon4
        temp.loc[5] = carbon5
        temp.loc[6] = carbon6
        temp.loc[7] = temp_vec1
        temp.loc[8] = temp_vec2
        temp.loc[9] = hydro1
        temp.loc[10] = hydro2
        temp.loc[11] = hydro3
        temp.loc[12] = hydro4
    
    # step 2: detect px/ox/mx by dot product
    Mol_Type = []
    full_pos = pd.DataFrame(columns = ['Num', 'Type', 'X', 'Y', 'Z'])
    epsilon_SB = 440
    sigma_SB = 5.27
    epsilon_CH3 = 85.51
    sigma_CH3 = 3.8
    epsilon_Hxyl = 15.03
    sigma_Hxyl = 2.42
    epsilon_Cxyl = 35.24
    sigma_Cxyl = 3.55
    epsilon_list = [epsilon_SB, 
                  epsilon_Cxyl, epsilon_Cxyl, epsilon_Cxyl, epsilon_Cxyl, epsilon_Cxyl, epsilon_Cxyl, 
                  epsilon_CH3, epsilon_CH3, epsilon_Hxyl, epsilon_Hxyl, epsilon_Hxyl, epsilon_Hxyl]
    sigma_list = [sigma_SB, 
                  sigma_Cxyl, sigma_Cxyl, sigma_Cxyl, sigma_Cxyl, sigma_Cxyl, sigma_Cxyl, 
                  sigma_CH3, sigma_CH3, sigma_Hxyl, sigma_Hxyl, sigma_Hxyl, sigma_Hxyl]
    Total_Type_P = 0
    Total_Type_O = 0
    Total_Type_M = 0
    Molecule_types = []
    for a in range(0,int(num_particles/3)):
        vec_origin = df.iloc[3*a]
        vec_1 = df.iloc[3*a+1]
        vec_2 = df.iloc[3*a+2]
        # change the types for CH3
        vec_1['Type'] = 7 # used the molecule definitions of Sp_p-xylene.def
        vec_2['Type'] = 8 
        mol_type = test_xy_type(vec_origin, vec_1, vec_2)
        if(mol_type == 3):
            grow_mx(a, vec_origin, vec_1, vec_2)
            Total_Type_M+=1
        elif(mol_type == 2):
            grow_ox(a, vec_origin, vec_1, vec_2)
            Total_Type_O+=1
        elif(mol_type == 1):
            grow_px(a, vec_origin, vec_1, vec_2)
            Total_Type_P+=1
        Mol_Type.append(mol_type)
    full_pos.to_csv(output_dir + "full_pos_Initial.csv")
    old_full = full_pos
    # now here, add framework atoms
    from pymatgen import Structure # keep these, we still needs to lattice params
    from mpmath import *
    file = input_dir + "Movies/System_0/Framework_final.pdb"
    a=0
    cart_abc = pd.DataFrame(columns = ['Num', 'Type', 'X', 'Y', 'Z'])
    with open(file) as f:
        for line in f:
            a+=1
            if("ATOM" in line):
                spline = line.split()
                num_molec = int(spline[1])
                type = spline[2]
                x = float(spline[4])
                y = float(spline[5])
                z = float(spline[6])
                pos = [num_molec, type, x, y, z]
                cart_abc.loc[a] = pos
    cart_abc = cart_abc.reset_index(drop = True)
    cart_abc.to_csv(output_dir + "final_MOF.csv")
    
    
    dfcart_abc=cart_abc # for recording
    frames = [full_pos, dfcart_abc]
    results = pd.concat(frames)
    results.reset_index(drop=True)
    results.to_csv(output_dir + "Initial_pos_with_frame.csv")
    
    # replace the types to numbers
    # get the sigmas from the forcefield file
    ff_file = "force_field_mixing_rules.def"
    a=0
    epsilons = []
    sigmas = []
    with open(ff_file) as f:
        for line in f:
            if(("LENNARD_JONES" in line) or "lennard-Jones" in line):
                print(line)
                spline = line.split()
                print(spline)
                type_atom = spline[0]
                type_epsilon = spline[2]
                type_sigma = spline[3]
                epsilons.append(float(type_epsilon))
                sigmas.append(float(type_sigma))
                # convert the types in type_abc to number types from the ff file
                cart_abc['Type'] = [a if (str(x) == type_atom or str(x) + '_' == type_atom) else x for x in cart_abc['Type']]
                a+=1
               
    
    def find(lst, b):
        result = []
        for i, x in enumerate(lst):
            if x == b:
                result.append(i)
        return result
    if(Total_Type_O > 0):
        df_Ox = pd.DataFrame(columns = ['Num', 'Type', 'X', 'Y', 'Z'])
        indexs = find(Mol_Type, 2)
        for s in range(0,len(indexs)):
            df_Ox = pd.concat([df_Ox, full_pos[full_pos['Num'] == indexs[s]]])
    if(Total_Type_M > 0):
        df_Mx = pd.DataFrame(columns = ['Num', 'Type', 'X', 'Y', 'Z'])
        indexs = find(Mol_Type, 3)
        for s in range(0,len(indexs)):
            df_Mx = pd.concat([df_Mx, full_pos[full_pos['Num'] == indexs[s]]])
    # Then get the Px df
    if(Total_Type_O > 0 and Total_Type_M > 0):
        drop_index = df_Mx.index.append(df_Ox.index)
        df_Ox = df_Ox.reset_index(drop = True)
        df_Mx = df_Mx.reset_index(drop = True)
    elif(Total_Type_O > 0):
        drop_index = df_Ox.index
        df_Ox = df_Ox.reset_index(drop = True)
    elif(Total_Type_M > 0):
        drop_index = df_Mx.index
        df_Mx = df_Mx.reset_index(drop = True)
    df_Px = full_pos.drop(index = drop_index)
    df_Px = df_Px.reset_index(drop = True)
    a=0
    scaling_count = 0
    Mol_Scale = Mol_Scale[0::3]# for each molecule, just keep one value
    #output_file = "restart_MTW_DD_P1_1.1.1_433.000000_1e+08"
    if(EBMC):
        output_file = filename
        if(charges):
            output_file = 'Charged_' + output_file
        with open(output_dir + output_file, 'w') as f:
            for m in range(0, len(header)):
                if("Lambda-factors" in header[m]):
                    placeholder = header[m].split()
                    placeholder.extend(repeat(header[m].split()[3], 10))
                    header[m] = ' '.join(placeholder)
                    header[m] = header[m] + '\n'
                if("three_site" in header[m]):
                    header[m] = header[m].replace("three_site", "Sp")
                f.write(header[m])
            f.write("\n")
            f.write("Component: 0     Adsorbate " + str(Total_Type_P) + " molecules of Sp_p-xylene\n")
            f.write("------------------------------------------------------------------------\n")
            # write the positions
            for m in range(0, Total_Type_P):
                for Atom in range(0,13):
                    slice = df_Px.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-position: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + str(slice['X']) + " " + str(slice['Y']) + " " + str(slice['Z']) + " " + "\n")
            # write velocity
            for m in range(0, Total_Type_P):
                for Atom in range(0,13):
                    slice = df_Px.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-velocity: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + str(0.0) + " " + str(0.0) + " " + str(0.0)+ "\n")
            # write force
            for m in range(0, Total_Type_P):
                for Atom in range(0,13):
                    slice = df_Px.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-force: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + str(0.0) + " " + str(0.0) + " " + str(0.0)+ "\n")
            # write charge: get the charges
            for m in range(0, Total_Type_P):
                for Atom in range(0,13):
                    slice = df_Px.iloc[13*m + Atom]
                    if(charges):
                      if(int(slice['Type']) == 0):
                          charge = 0.0
                      elif(int(slice['Type']) < 6):
                          charge = -0.115
                      else:
                          charge = 0.115
                    else:
                        charges = 0.0
                    f.write("Adsorbate-atom-charge: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + str(charge)+ "\n")
            # write the scaling
            for m in range(0, Total_Type_P):
                for Atom in range(0,13):
                    slice = df_Px.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-scaling: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + str(Mol_Scale[scaling_count])+ "\n")
                scaling_count+=1
            # write fix
            for m in range(0, Total_Type_P):
                for Atom in range(0,13):
                    slice = df_Px.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-fixed: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + "0 0 0"+ "\n")
        
            f.write("\nComponent: 1     Adsorbate " + str(Total_Type_O) + " molecules of Sp_o-xylene\n")
            f.write("------------------------------------------------------------------------\n")
            # write the positions
            for m in range(0, Total_Type_O):
                for Atom in range(0,13):
                    slice = df_Ox.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-position: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + str(slice['X']) + " " + str(slice['Y']) + " " + str(slice['Z']) + " " + "\n")
            # write velocity
            for m in range(0, Total_Type_O):
                for Atom in range(0,13):
                    slice = df_Ox.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-velocity: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + str(0.0) + " " + str(0.0) + " " + str(0.0)+ "\n")
            # write force
            for m in range(0, Total_Type_O):
                for Atom in range(0,13):
                    slice = df_Ox.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-force: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + str(0.0) + " " + str(0.0) + " " + str(0.0)+ "\n")
            # write charge
            for m in range(0, Total_Type_O):
                for Atom in range(0,13):
                    slice = df_Ox.iloc[13*m + Atom]
                    if(charges):
                      if(int(slice['Type']) == 0):
                          charge = 0.0
                      elif(int(slice['Type']) < 6):
                          charge = -0.115
                      else:
                          charge = 0.115
                    else:
                        charge = 0.0
                        
                    f.write("Adsorbate-atom-charge: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + str(charge)+ "\n")
            # write the scaling
            for m in range(0, Total_Type_O):
                for Atom in range(0,13):
                    slice = df_Ox.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-scaling: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + str(Mol_Scale[scaling_count])+ "\n")
                scaling_count+=1
            # write fix
            for m in range(0, Total_Type_O):
                for Atom in range(0,13):
                    slice = df_Ox.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-fixed: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + "0 0 0"+ "\n")
        
            f.write("\nComponent: 2     Adsorbate " + str(Total_Type_M) + " molecules of Sp_m-xylene\n")
            f.write("------------------------------------------------------------------------\n")
            # write the positions
            for m in range(0, Total_Type_M):
                for Atom in range(0,13):
                    slice = df_Mx.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-position: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + str(slice['X']) + " " + str(slice['Y']) + " " + str(slice['Z']) + " " + "\n")
            # write velocity
            for m in range(0, Total_Type_M):
                for Atom in range(0,13):
                    slice = df_Mx.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-velocity: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + str(0.0) + " " + str(0.0) + " " + str(0.0)+ "\n")
            # write force
            for m in range(0, Total_Type_M):
                for Atom in range(0,13):
                    slice = df_Mx.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-force: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + str(0.0) + " " + str(0.0) + " " + str(0.0)+ "\n")
            # write charge
            for m in range(0, Total_Type_M):
                for Atom in range(0,13):
                    slice = df_Mx.iloc[13*m + Atom]
                    if(charges):
                      if(int(slice['Type']) == 0):
                          charge = 0.0
                      elif(int(slice['Type']) < 6):
                          charge = -0.115
                      else:
                          charge = 0.115
                    else:
                        charge = 0.0
                    f.write("Adsorbate-atom-charge: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + str(charge)+ "\n")
            # write the scaling
            for m in range(0, Total_Type_M):
                for Atom in range(0,13):
                    slice = df_Mx.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-scaling: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + str(Mol_Scale[scaling_count])+ "\n")
                scaling_count+=1
            # write fix
            for m in range(0, Total_Type_M):
                for Atom in range(0,13):
                    slice = df_Mx.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-fixed: " + str(int(slice['Num'])) + " " + str(int(slice['Type'])) + " " + "0 0 0"+ "\n")
    else:
        # if Not EBMC, then we don't need to spherical benzene atom (1st atom)
        output_file = "NOEB_" + filename
        if(charges):
            output_file = 'Charged_' + output_file
        with open(output_dir + output_file, 'w') as f:
            for m in range(0, len(header)):
                if("Lambda-factors" in header[m]):
                    placeholder = header[m].split()
                    placeholder.extend(repeat(header[m].split()[3], 9)) # one less atom, one less lambda number needed
                    header[m] = ' '.join(placeholder)
                    header[m] = header[m] + '\n'
                if("three_site" in header[m]):
                    header[m] = header[m].replace("three_site_", "")
                f.write(header[m])
            f.write("\n")
            f.write("Component: 0     Adsorbate " + str(Total_Type_P) + " molecules of p-xylene\n")
            f.write("------------------------------------------------------------------------\n")
            # write the positions
            for m in range(0, Total_Type_P):
                for Atom in range(1,13):
                    slice = df_Px.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-position: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + str(slice['X']) + " " + str(slice['Y']) + " " + str(slice['Z']) + " " + "\n")
            # write velocity
            for m in range(0, Total_Type_P):
                for Atom in range(1,13):
                    slice = df_Px.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-velocity: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + str(0.0) + " " + str(0.0) + " " + str(0.0)+ "\n")
            # write force
            for m in range(0, Total_Type_P):
                for Atom in range(1,13):
                    slice = df_Px.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-force: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + str(0.0) + " " + str(0.0) + " " + str(0.0)+ "\n")
            # write charge
            for m in range(0, Total_Type_P):
                for Atom in range(1,13):
                    slice = df_Px.iloc[13*m + Atom]
                    if(charges):
                      if((int(slice['Type']) - 1) < 6):
                          charge = -0.115
                      else:
                          charge = 0.115
                    else:
                        charge = 0.0
                    f.write("Adsorbate-atom-charge: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + str(charge)+ "\n")
            # write the scaling
            for m in range(0, Total_Type_P):
                for Atom in range(1,13):
                    slice = df_Px.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-scaling: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + str(Mol_Scale[scaling_count])+ "\n")
                scaling_count+=1
            # write fix
            for m in range(0, Total_Type_P):
                for Atom in range(1,13):
                    slice = df_Px.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-fixed: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + "0 0 0"+ "\n")
        
            f.write("\nComponent: 1     Adsorbate " + str(Total_Type_O) + " molecules of o-xylene\n")
            f.write("------------------------------------------------------------------------\n")
            # write the positions
            for m in range(0, Total_Type_O):
                for Atom in range(1,13):
                    slice = df_Ox.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-position: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + str(slice['X']) + " " + str(slice['Y']) + " " + str(slice['Z']) + " " + "\n")
            # write velocity
            for m in range(0, Total_Type_O):
                for Atom in range(1,13):
                    slice = df_Ox.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-velocity: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + str(0.0) + " " + str(0.0) + " " + str(0.0)+ "\n")
            # write force
            for m in range(0, Total_Type_O):
                for Atom in range(1,13):
                    slice = df_Ox.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-force: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + str(0.0) + " " + str(0.0) + " " + str(0.0)+ "\n")
            # write charge: no central benzene atom
            for m in range(0, Total_Type_O):
                for Atom in range(1,13):
                    slice = df_Ox.iloc[13*m + Atom]
                    if(charges):
                      if((int(slice['Type']) - 1) < 6):
                          charge = -0.115
                      else:
                          charge = 0.115
                    else:
                        charge = 0.0
                    f.write("Adsorbate-atom-charge: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + str(charge)+ "\n")
            # write the scaling
            for m in range(0, Total_Type_O):
                for Atom in range(1,13):
                    slice = df_Ox.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-scaling: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + str(Mol_Scale[scaling_count])+ "\n")
                scaling_count+=1
            # write fix
            for m in range(0, Total_Type_O):
                for Atom in range(1,13):
                    slice = df_Ox.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-fixed: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + "0 0 0"+ "\n")
        
            f.write("\nComponent: 2     Adsorbate " + str(Total_Type_M) + " molecules of m-xylene\n")
            f.write("------------------------------------------------------------------------\n")
            # write the positions
            for m in range(0, Total_Type_M):
                for Atom in range(1,13):
                    slice = df_Mx.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-position: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + str(slice['X']) + " " + str(slice['Y']) + " " + str(slice['Z']) + " " + "\n")
            # write velocity
            for m in range(0, Total_Type_M):
                for Atom in range(1,13):
                    slice = df_Mx.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-velocity: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + str(0.0) + " " + str(0.0) + " " + str(0.0)+ "\n")
            # write force
            for m in range(0, Total_Type_M):
                for Atom in range(1,13):
                    slice = df_Mx.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-force: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + str(0.0) + " " + str(0.0) + " " + str(0.0)+ "\n")
            # write charge: no central benzene atom
            for m in range(0, Total_Type_M):
                for Atom in range(1,13):
                    slice = df_Mx.iloc[13*m + Atom]
                    if(charges):
                      if((int(slice['Type']) - 1) < 6):
                          charge = -0.115
                      else:
                          charge = 0.115
                    else:
                        charge = 0.0
                    f.write("Adsorbate-atom-charge: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + str(charge)+ "\n")
            # write the scaling
            for m in range(0, Total_Type_M):
                for Atom in range(1,13):
                    slice = df_Mx.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-scaling: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + str(Mol_Scale[scaling_count])+ "\n")
                scaling_count+=1
            # write fix
            for m in range(0, Total_Type_M):
                for Atom in range(1,13):
                    slice = df_Mx.iloc[13*m + Atom]
                    f.write("Adsorbate-atom-fixed: " + str(int(slice['Num'])) + " " + str(int(slice['Type'] - 1)) + " " + "0 0 0"+ "\n")
            