#!/usr/bin/env python
# ***********************************************************
#   File: vasp2cif[.py]
#   Description: a tool to make CIF format files out of
#               VASP POSCAR+POTCAR/OUTCAR files.
#               Output files acquire a .cif extension
#               example: POSCAR --> POSCAR.cif
#
#   Copyright 2008-2013 Peter Larsson
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
#   Revision history:
#      2013-07-14  Torbjorn Bjorkman
#        - Runs Harold Stoke's FINDSYM program (if available) 
#          to determine the space group and standard setting.
#      2012-12-15  Torbjorn Bjorkman
#        - Extracts geometries also from OUTCAR files to a
#          set of blocks in a single CIF file. Useful for
#          visualization of a relaxation or MD run.
#      2011-05-06  Torbjorn Bjorkman
#        - You can now give many input POSCAR files and they
#          will be generated as separate data blocks in the
#          output CIF. Useful for visualization.
#        - Support for cartesian coordinates also for non-
#          orthorhombic systems.
#      2010-01-08  Peter Larsson
#        - Support for VASP 5 style CONTCAR files
#      2009-09-29  Peter Larsson
#        - More descriptive help for command line options
#      2009-04-17  Peter Larsson
#        - More robust error handling and support
#          for Cartesian coordinates in orthorhombic cells
#      2008-10-13  Peter Larsson
#        - Ported to Python and added support for
#          "Selective Dynamics" format in POSCAR and
#          volume scaling
#      2006-10-24  Peter Larsson
#        - Original version in Ruby.
# ***********************************************************

import os
import sys
import commands
import math
import re
from subprocess import Popen,call,PIPE
from optparse import OptionParser

# Input parser
parser = OptionParser()
parser.add_option("-v","--verbose",dest="verbose",help="Print CIF to stdout",action="store_true")
parser.add_option("-o","--output",dest="output",help="Save CIF to named file",metavar="FILE")
parser.add_option("-e","--elements",dest="elements",help="""Supply elements if no POTCAR is present. Example: --elements="Fe,Co,Ni" """,metavar="list of elements")
parser.add_option("--findsym-tolerance",dest="findsymtol",help="""Tolerance used for FINDSYM (default=0, minimal value).""")
parser.add_option("--no-findsym",dest="nofindsym",help="""Don't run FINDSYM to find symmetry of the crystal.""")
(options,args) = parser.parse_args()

####### FINDSYM ######
# Check if FINDSYM is up and running
global findsym
global findsymtolerance
full2shortHM = { '2/m 2/m 2/m' : 'mmm',
		 '4/m 2/m 2/m' : '4/mmm',
		 '-3 2/m' : '-3m',
		 '6/m 2/m 2/m' : '6/mmm',
		 '2/m -3' : 'm-3',
		 '4/m -3 2/m' : 'm-3m'}
fullHMlist = ['4/m -3 2/m',
	      '2/m -3',
	      '6/m 2/m 2/m',
	      '-3 2/m',
	      '4/m 2/m 2/m',
	      '2/m 2/m 2/m']
if options.findsymtol:
	findsymtolerance = float(options.findsymtol)
else:
	findsymtolerance = 0
try:
	p=Popen(["findsym"],stdin=PIPE,stdout=PIPE)
	p.stdin.write("blah\n0\n2\n1 1 1 90 90 90\n1\n1 0 0\n0 1 0\n0 0 1\n1\n1\n0 0 0")
	output = p.communicate()[0]
	if p.returncode != 0:
		findsym = False
	else:
		# Check if it worked
		if re.match("Error",output):
			findsym = False
		else:
			# If we get here everything should be OK
			findsym = True
except:
	findsym = False
if options.nofindsym:
	findsym = False

# Classes
class Cell:
	def __init__(self):
		label = ""
		latticevectors = None
		a = None
		b = None
		c = None
		alpha = None
		beta = None
		gamma = None
		HMSymbol = "'P 1'"
		sites = [] # list of (elementname, x, y, z) tuples

def ciffilestring(cell):
	# Make CIF header with cell parameters and coord record info
	outstring = ""
	nofindsym = not findsym
	if findsym:
		# Construct FINDSYM input
		findsymstring = " \n"
		findsymstring += str(findsymtolerance)+"\n"
#		findsymstring += "2\n"
#		findsymstring += str(cell.a)+" "+str(cell.b)+" "+str(cell.c)+" "+str(cell.alpha)+" "+str(cell.beta)+" "+str(cell.gamma)+"\n"
		findsymstring += "1\n"
		findsymstring += str(cell.latticevectors[0][0])+" "+str(cell.latticevectors[0][1])+" "+str(cell.latticevectors[0][2])+"\n"
		findsymstring += str(cell.latticevectors[1][0])+" "+str(cell.latticevectors[1][1])+" "+str(cell.latticevectors[1][2])+"\n"
		findsymstring += str(cell.latticevectors[2][0])+" "+str(cell.latticevectors[2][1])+" "+str(cell.latticevectors[2][2])+"\n"
		findsymstring += "1\n"
		findsymstring += "1 0 0\n"
		findsymstring += "0 1 0\n"
		findsymstring += "0 0 1\n"
		findsymstring += str(len(cell.sites))+"\n"
		for a in cell.sites:
			findsymstring += a[0]+" "
		findsymstring += "\n"
		for a in cell.sites:
			findsymstring += "%1.15f   %1.15f   %1.15f\n"%(a[1],a[2],a[3])
		# Call findsym with input and catch output
		p=Popen(["findsym"],stdin=PIPE,stdout=PIPE)
		p.stdin.write(findsymstring)
		findsymoutput=p.communicate()[0]
		findsymoutlines = findsymoutput.split("\n")
		# Delete log
		call(['rm','-f','findsym.log'])

		# Check that FINDSYM returned a cif at the end 
		cifout = False
		i = 0
		for line in findsymoutlines:
			i += 1
			if re.match('_audit_creation_method',line):
				cifout = True
				break
		if not cifout:
			print "***Error: FINDSYM failed to produce a CIF file, no symmetrization done."
			nofindsym = True

	outstring += "data_" + cell.label.strip(" ")+"\n"
	if nofindsym:
		outstring += "_audit_creation_method   'vasp2cif'\n"
	else:
		outstring += "_audit_creation_method   'vasp2cif/FINDSYM'\n"
		# Append FINDSYM cif to outstring.
		for line in findsymoutlines[i:]:
			# Convert to short version of H-M symbol.
			if re.match('_symmetry_space_group_name_H-M',line):
				tmp = line
				for k in fullHMlist:
					if re.search(k,tmp):
						tmp = tmp.replace(k,full2shortHM[k])
				outstring += tmp+"\n"
			else:
				outstring += line+"\n"

	if nofindsym:
		# Construct outstring
		outstring += "_cell_length_a    " + str(cell.a)+"\n"
		outstring += "_cell_length_b    " + str(cell.b)+"\n"
		outstring += "_cell_length_c    " + str(cell.c)+"\n"
		outstring += "_cell_angle_alpha    " + str(cell.alpha)+"\n"
		outstring += "_cell_angle_beta    " + str(cell.beta)+"\n"
		outstring += "_cell_angle_gamma    " + str(cell.gamma)+"\n"
		outstring += "\n"
		outstring += "_symmetry_space_group_name_H-M    "+cell.HMSymb+"\n"
		outstring += "loop_\n"
		outstring += "_atom_site_label\n"
		outstring += "_atom_site_type_symbol\n"
		outstring += "_atom_site_fract_x\n"
		outstring += "_atom_site_fract_y\n"
		outstring += "_atom_site_fract_z\n"
		outstring += "_atom_site_occupancy\n"
		i = 1
		for a in cell.sites:
			outstring += "%s%i   %s   %1.15f   %1.15f   %1.15f   1.0\n" % (a[0], i, a[0], a[1], a[2], a[3])
			i += 1		 
		outstring += "\n\n"
	return outstring

def gstrip(s):
	#Strip all whitespace in string s
	return re.sub("\s+" , "", s)

# Return OUTCAR or POSCAR depending on file type.
# If not identified, return UNKNOWN
def filetype(f):
	filetype = "UNKNOWN"
	# Something that should recognize an OUTCAR file
	sym = re.compile("Analysis of symmetry for")
	# POSCAR identified if line 2-5 can be interpreted as
	# a length scale followed by three lattice vectors
	lines = []
	for i in range(6):
		lines.append(f.readline())
	try:
		t = float(lines[1].split()[0])
		for i in range(2,5):
			t = [float(s) for s in lines[i].split()[:3]]
		filetype = "POSCAR"
	except:
		for line in f:
			if sym.match(line):
				filetype = "OUTCAR"
				break
	f.seek(0) # Rewind file
	return filetype
	
def mvmult3(mat,vec):
	# matrix-vector multiplication
	w = []
	for i in range(3):
		t = 0
		for j in range(3):
			t += mat[j][i]*vec[j]
		w.append(t)
	return w

def det3(m):
	# Determinant of 3x3 dimensional matrix
	a = m[1][1]*m[2][2]-m[1][2]*m[2][1]
	b = m[1][2]*m[2][0]-m[1][0]*m[2][2]
	c = m[1][0]*m[2][1]-m[1][1]*m[2][0]
	return m[0][0]*a + m[0][1]*b + m[0][2]*c

def minv3(m):
	# Inverse of 3x3 dimensional matrix
	det = det3(m)
	w = [[(m[1][1]*m[2][2]-m[1][2]*m[2][1])/det, (m[0][2]*m[2][1]-m[0][1]*m[2][2])/det, (m[0][1]*m[1][2]-m[0][2]*m[1][1])/det],
	     [(m[1][2]*m[2][0]-m[1][0]*m[2][2])/det, (m[0][0]*m[2][2]-m[0][2]*m[2][0])/det, (m[0][2]*m[1][0]-m[0][0]*m[1][2])/det],
	     [(m[1][0]*m[2][1]-m[1][1]*m[2][0])/det, (m[0][1]*m[2][0]-m[0][0]*m[2][1])/det, (m[0][0]*m[1][1]-m[0][1]*m[1][0])/det]]
	return w


#Let's get started, read input files.
# Store in lists as (inputfile,filename) tuples.
if len(args) == 0:
	#Pipe mode, read and write to stdin and stdout
	input_files = [(sys.stdin,None)]
	cif_file = sys.stdout
elif len(args) == 1:
	#Write to input.cif
	input_files = [(file(arg,'r'),arg) for arg in args]
	if options.output:
		cif_file = file(options.output,'w')
	else:
		cif_file = file(args[0] + ".cif",'w')
else:
	#Write to input.cif
	input_files = [(file(arg,'r'),arg) for arg in args[0:]]
	if options.output:
		cif_file = file(options.output,'w')
	else:
		cif_file = file(args[0] + "_etc.cif",'w')

# Initialize Cell object.
cell = Cell()
cell.HMSymb = "'P 1'"
# loop over input files
inputfilenr = 1
cifblocknr = 1
for input_file,filename in input_files:
	if filetype(input_file) == "POSCAR":
		poscar = input_file.readlines()

		# CIF block number
		cell.label = str(cifblocknr)

		#We need to determine the data format, VASP 5 stores element names in line 5
		if gstrip(poscar[5]).isdigit():
			#Old school format
			vasp5 = False
			offset = 0
		elif gstrip(poscar[5]).isalpha() and gstrip(poscar[6]).isdigit():
			#Looks like vasp5 like format
			vasp5 = True
			offset = 1

		#First deal with potential POTCAR problems
		atoms = []
		if options.elements:
			#Read atoms from supplied string, eg "Li,Fe,Si,O"
			atoms = options.elements.split(",")
			assert(len(atoms) > 0)
		else:		
			if vasp5:
				#Read elements from line 5
				words = poscar[5].split()
				atoms = [w.strip() for w in words]
			else:
				#Try to read atoms from POTCAR
				if not os.path.exists("POTCAR"):
					sys.stderr.write("ERROR: Cannot find POTCAR. Please supply atom labels with the -e flag.\n")
					sys.exit(1)

				potcar_lines = commands.getoutput("grep TITEL POTCAR").split("\n")
				if len(potcar_lines) == 0:
					sys.stderr.write("ERROR: POTCAR file exists, but is empty? Supply atom labels with the -e flag.\n")
					sys.exit(1)

				for line in potcar_lines:
					words = line.split()
					assert(words[0] == 'TITEL')

					#Note, we need the split _ to deal with names like "Li_sv"
					atoms.append(words[3].split("_")[0])

		#Lattice scaling factor
		lattice_constant = float(poscar[1].strip())

		#Dealing with volume scaling in POSCAR
		final_volume = -lattice_constant
		scale_volume = False
		if lattice_constant < 0.0:
			lattice_constant = 1.0
			scale_volume = True

		#Read cell vectors
		a = []
		b = []
		c = []
		a.append(lattice_constant*float(poscar[2].split()[0].strip()))
		a.append(lattice_constant*float(poscar[2].split()[1].strip()))
		a.append(lattice_constant*float(poscar[2].split()[2].strip()))

		b.append(lattice_constant*float(poscar[3].split()[0].strip()))
		b.append(lattice_constant*float(poscar[3].split()[1].strip()))
		b.append(lattice_constant*float(poscar[3].split()[2].strip()))

		c.append(lattice_constant*float(poscar[4].split()[0].strip()))
		c.append(lattice_constant*float(poscar[4].split()[1].strip()))
		c.append(lattice_constant*float(poscar[4].split()[2].strip()))

		unscaled_volume = a[0]*b[1]*c[2]-a[0]*b[2]*c[1]+a[1]*b[2]*c[0]-a[1]*b[0]*c[2]+a[2]*b[0]*c[1]-a[2]*b[1]*c[0]

		if scale_volume:
			lattice_constant = (final_volume/unscaled_volume)**(1.0/3.0)
			a = map(lambda x: lattice_constant*x,a)
			b = map(lambda x: lattice_constant*x,b)
			c = map(lambda x: lattice_constant*x,c)
			volume = a[0]*b[1]*c[2]-a[0]*b[2]*c[1]+a[1]*b[2]*c[0]-a[1]*b[0]*c[2]+a[2]*b[0]*c[1]-a[2]*b[1]*c[0]

		cell.a = math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
		cell.b = math.sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2])
		cell.c = math.sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2])

		cell.alpha = math.acos((b[0]*c[0]+b[1]*c[1]+b[2]*c[2])/(cell.b*cell.c))*180/math.pi
		cell.beta = math.acos((a[0]*c[0]+a[1]*c[1]+a[2]*c[2])/(cell.a*cell.c))*180/math.pi
		cell.gamma = math.acos((b[0]*a[0]+b[1]*a[1]+b[2]*a[2])/(cell.a*cell.b))*180/math.pi

		cell.latticevectors = [a,b,c]

		#Read atoms counts and make label array
		atomlabels = []

		atomcounts = poscar[5+offset].split()

		if len(atomcounts) != len(atoms):
			sys.stderr.write("ERROR: Not the same number of atom species in POTCAR and POSCAR. Please check.\n")
			sys.exit(1)

		n_atoms = 0
		for i in range(0,len(atomcounts)):
			n = int(atomcounts[i].strip())
			n_atoms = n_atoms + n
			for j in range(0,n):
				atomlabels.append(atoms[i])

		#Check for selective dynamics
		if poscar[6+offset].upper()[0] == 'S':
			offset = offset + 7
		else:
			offset = offset + 6

		#Check for direct coordinates
		direct_coordinates = True
		if poscar[offset].upper()[0] == 'D':
			direct_coordinates = True
		if poscar[offset].upper()[0] == 'C':
			direct_coordinates = False
			lattice_vectors = [a,b,c]
			inverse_lattice_vectors = minv3(lattice_vectors)

		#Scan and print atomic positions from offset
		if len(atomlabels) > (len(poscar)-offset):
			sys.stderr.write(("WARNING: vasp2cif expected to find %d coordinates, but there are only %d coordinate lines in the file!\n") % (len(atomlabels),len(poscar)-offset))
			atomlabels = atomlabels[0:len(poscar)-offset-1]

		cell.sites = []
		for i in range(0,len(atomlabels)):
			#extract first three fields in POSCAR line
			coords = map(float,poscar[i+offset+1].split()[0:3])

			if not direct_coordinates:
				coords = mvmult3(inverse_lattice_vectors, coords)

			cell.sites.append((atomlabels[i],coords[0],coords[1],coords[2]))

		# Print cell to cif files
		cifstring = ciffilestring(cell)
		cif_file.write(cifstring)
		if options.verbose:
			sys.stdout.write(cifstring)
		# increment cif block counter
		cifblocknr += 1
	
	elif filetype(input_file) == "OUTCAR":
		# First find elements and how many of each.
		atoms = []
		titellines = commands.getoutput("grep TITEL "+filename).split("\n")
		if len(titellines) == 0:
			sys.stderr.write("ERROR: Cannot read elements. Damaged OUTCAR file?\n")
			sys.exit(1)
		for line in titellines:
			words = line.split()
			assert(words[0] == 'TITEL')
			#Note, we need the split _ to deal with names like "Li_sv"
			atoms.append(words[3].split("_")[0])
		# How many of each?
		natoms = [int(s) for s in commands.getoutput("grep 'ions per type =' "+filename).split()[4:]]
		# Set up initial position array
		i = 0
		cell.sites = []
		for a in atoms:
			for j in range(natoms[i]):
				cell.sites.append((a,0.0,0.0,0.0))
			i += 1
		# Precompiled regular expressions.
		re_iter = re.compile("aborting loop because EDIFF is reached")
		re_lattice = re.compile("direct lattice vectors")
		re_positions = re.compile("POSITION")
		latticevectors = []
		introread = False
		vectorsread = False
		positionsread = False
		vecline = 5
		posline = sum(natoms)+1
		linenr = 0
		for line in input_file:
			linenr += 1
			# Start for lattices and positions after the end of first iteration.
			if not introread:
				if re_iter.search(line):
					introread = True
				continue
			# Get lattice vectors.
			if re_lattice.search(line):
				vecline = 1
				latticevectors = []
				continue
			if vecline <= 3:
				latticevectors.append([float(s) for s in line.split()[0:3]])
			if vecline == 4:
				# Lattice vectors read, set parameters
				vectorsread = True
				cell.latticevectors = latticevectors
				inverse_lattice_vectors = minv3(latticevectors)
				a,b,c = latticevectors[0],latticevectors[1],latticevectors[2]
				# crystallographic parameters
				cell.a = math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
				cell.b = math.sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2])
				cell.c = math.sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2])

				cell.alpha = math.acos((b[0]*c[0]+b[1]*c[1]+b[2]*c[2])/(cell.b*cell.c))*180/math.pi
				cell.beta = math.acos((a[0]*c[0]+a[1]*c[1]+a[2]*c[2])/(cell.a*cell.c))*180/math.pi
				cell.gamma = math.acos((b[0]*a[0]+b[1]*a[1]+b[2]*a[2])/(cell.a*cell.b))*180/math.pi
				
			# Read positions
			if re_positions.search(line):
				posline = -2
			if 0 <= posline < sum(natoms):
				c = [float(p) for p in line.split()[0:3]]
				c = mvmult3(inverse_lattice_vectors, c)
				cell.sites[posline] = (cell.sites[posline][0],c[0],c[1],c[2])
			if posline == sum(natoms) and introread:
				# Positions read, now print cell.
				cell.label = str(cifblocknr)
				cifstring = ciffilestring(cell)
				cif_file.write(cifstring)
				if options.verbose:
					sys.stdout.write(cifstring)
				# increment cif block counter
				cifblocknr += 1
			posline += 1
			vecline += 1
	else:
		sys.stderr.write("ERROR: Format of file %i not recognized.\n"%inputfilenr)
		sys.exit(1)
	# increment input file counter
	inputfilenr += 1

