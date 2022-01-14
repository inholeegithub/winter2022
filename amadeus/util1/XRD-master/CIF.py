from math import acos, pi, ceil, sin,cos,sqrt
import numpy as np
import re

class cif(object):
    """a class of cif reader
    Attributes:
        wavelength: default: 1.54181a, namely Cu-Ka
        max2theta: the range of 2theta angle
        intensity: intensities for all hkl planes
        pxrd: powder diffraction data
    """

    def __init__(self, filename):
        """Return a XRD object with the proper info"""
        self.from_file(filename)
        self.parse_cell()
        self.parse_atom()
        self.apply_symops()

    def from_file(self, filename):
        cif = np.genfromtxt(filename, dtype=str, delimiter='\n')
        
        # 3 modes in each flag:  
        # 0: not started; 
        # 1: reading; 
        # 2: done
        flags = {'cell':0, 'symops':0, 'atom':0}

        atom = {}
        cell = {}
        symops = {'string':[], 'matrix':[]}

        for lines in cif:

            if 'loop_' in lines:  
                #if a _loop lines starts, the current reading flag switch to 0
                for item in flags.keys():
                    if flags[item] == 1:
                        flags[item] = 2

            elif '_cell_' in lines:
                #_cell_length_a          4.77985

                flags['cell'] = 1
                cell_str = lines.split()
                item = cell_str[0].replace(' ','')
                value = float(cell_str[1].split("(")[0])
                cell[item] = value

            elif '_symmetry_equiv_pos_as_xyz' in lines:
                #_symmetry_equiv_pos_as_xyz
                flags['symops'] = 1
      
            elif '_space_group_symop_operation_xyz' in lines:
                #_space_group_symop_operation_xyz
                flags['symops'] = 1
                
            elif flags['symops'] == 1:
                #1, 'x, y, z'
                #    x, -y, z
                raw_line = lines.strip().strip("'").split(' ', 1)
                if raw_line[0].isdigit():     
                    sym_str = raw_line[1].strip("'")
                else:
                    sym_str = lines.strip().strip("'").replace(' ', '')
                sym_str = sym_str.replace("'","")
                symops['string'].append(sym_str)
                symops['matrix'].append(self.xyz2sym_ops(sym_str))

            elif '_atom_site' in lines: 
                flags['atom'] = 1
                atom_str = lines.replace(' ','')
                item = atom_str
                atom[item] = []

            elif flags['atom'] == 1:
                raw_line = lines.split()
                for i, item in enumerate(atom.keys()):
                    raw_text = raw_line[i]
                    
                    if item.find('fract')>0:
                       value = float(raw_text.split("(")[0])
                    elif item.find('symbol')>0:
                       m_symbol = re.compile("([A-Z]+[a-z]*)")
                       value = m_symbol.findall(raw_text)
                       #print(raw_text, value)
                    else:
                       value = raw_text
                       
                    atom[item].append(value)

            elif flags['cell'] + flags['symops'] + flags['atom'] == 6:
                break

        self.cell = cell
        self.atom = atom
        self.symops = symops
   
    def parse_cell(self):
        cell_para = np.zeros(6)
        cell = self.cell
        for item in cell.keys():
            if item.find('_length_a') > 0:
                cell_para[0] = cell[item]
            elif item.find('_length_b') > 0:
                cell_para[1] = cell[item]
            elif item.find('_length_c') > 0:
                cell_para[2] = cell[item]
            elif item.find('_angle_alpha') > 0:
                cell_para[3] = cell[item]
            elif item.find('_angle_beta') > 0:
                cell_para[4] = cell[item]
            elif item.find('_angle_gamma') > 0:
                cell_para[5] = cell[item]
        self.cell_para = cell_para

    def parse_atom(self):
        atom = self.atom
        N_atom = len(atom['_atom_site_fract_x'])
        cif_xyz = np.zeros([N_atom, 3])

        for item in atom.keys():
            if item.find('_fract_x') > 0:
                cif_xyz[:,0] = np.array(atom[item])
            elif item.find('_fract_y') > 0:
                cif_xyz[:,1] = np.array(atom[item])
            elif item.find('_fract_z') > 0:
                cif_xyz[:,2] = np.array(atom[item])

        self.cif_xyz = cif_xyz

    #generates all coordinates from rotation matrices and translation vectors
    def apply_symops(self):
        fract_xyz = self.cif_xyz
        symops_matrix = self.symops['matrix']
        atom_type = self.atom['_atom_site_type_symbol']
        sym_coordinates = {}

        for ii,item in enumerate(atom_type):
            sym_coordinates[str(item)] = []
            mat_vec = symops_matrix[ii]
            sym_temp = np.dot(mat_vec[0], fract_xyz[ii].transpose()) + mat_vec[1]
            sym_coordinates[str(item)].append(sym_temp)

        self.coordinate, self.composition = self.remove_duplicate(sym_coordinates)

    #remove equivalent points and keep the unique ones
    #get the numbers of atoms per species
    @staticmethod
    def remove_duplicate(sym_coordinates):
        coordinate = []
        composition = []
        for item in sym_coordinates.keys():
            raw_equiv = np.array(sym_coordinates[item])
            raw_equiv = raw_equiv - np.floor(raw_equiv)
            raw_equiv = np.around(raw_equiv, 4)
            raw_equiv = np.unique(raw_equiv, axis=0)
            composition.append(len(raw_equiv))
            if coordinate == []:
                coordinate = raw_equiv
            else:
                coordinate = np.concatenate((coordinate,raw_equiv),axis=0)

        return coordinate, composition


    #function generates rotation matrices and translation vectors from equivalent points
    @staticmethod
    def xyz2sym_ops(string):
        #rotational matrix dictionary
        rot_dic = {}
        rot_dic['x'] = np.array([1.0,0,0])
        rot_dic['y'] = np.array([0,1.0,0])
        rot_dic['z'] = np.array([0,0,1.0])
        parts = string.strip().replace(' ','').lower().split(',')
        rot_mat = []
        rot_temp = np.array([0.,0.,0.])
        trans_vec = np.array([0.,0.,0.])
        #use re module to read xyz strings
        m_rot = re.compile(r"([+-]?)([\d\.]*)/?([\d\.]*)([x-z])")
        m_trans = re.compile(r"([+-]?)([\d\.]+)/?([\d\.]*)(?![x-z])")
        for jj,item in enumerate(parts):
            #rotation matrix
            for ii,m in enumerate(m_rot.finditer(item)):
                coef = -1 if m.group(1) == '-' else 1
                if m.group(2) != '':
                    if m.group(3) != '':
                        coef *= float(m.group(2))/float(m.group(3))
                    else:
                        coef *= float(m.group(2))
                if ii == 0:                  
                    rot_temp = rot_dic[m.group(4)]*coef
                else:
                    rot_temp += rot_dic[m.group(4)]*coef
            rot_mat.append(rot_temp)
            #translation vector
            for m in m_trans.finditer(item):
                coef = -1 if m.group(1) == '-' else 1
                if m.group(3) != '':
                    coef = float(m.group(2))/float(m.group(3))
                else:
                    coef = float(m.group(2))
                trans_vec[jj] = 1.0*coef
        return (rot_mat, trans_vec)
                   

from optparse import OptionParser
import pandas as pd
from tabulate import tabulate

if __name__ == "__main__":
    
    cif = cif('P3N5-alpha-70GPa.cif')
    print(cif.cell_para)
    print(cif.coordinate)
    print(cif.composition)
