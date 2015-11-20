#!/usr/bin/python
# -*- coding: utf-8 -*-
from argparse import ArgumentParser
import os, numpy
from Bio.PDB import PDBParser, PDBIO
import matplotlib.pyplot as plt

"""

Usage: python ./contact_map.py -f file.pdb -o output_name [-d directory, -cmax cutoff_max, -cmin cutoff_min, -fc first_chain, -sc second_chain]

"""

######################## Parsing stuff ########################

# Defaults
cutoff_max = False
cutoff_min = False
directory = ""
chainA = "A"
chainB = "B"

parser = ArgumentParser(description=""" Makes a contact map between two chains (Defaults A and B) from a single PDB. \n
The contact map is white if there's no contact within the specified cutoff. \n
Also outputs a modified .pdb file, in case the submitted one has problems. \n
Outputs the separate models too""")

# Named arguments
parser.add_argument("-f", "--file", help="The name of the input pdb file")
parser.add_argument("-o", "--output", help="The generic name of the output file(s)")

# Optional arguments
parser.add_argument("-d", "--directory", help="Output directory ? Default is the current directory.")
parser.add_argument("-cmax", "--cutoff_max", help="Optional maximum distance cutoff")
parser.add_argument("-cmin", "--cutoff_min", help="Optional minimum distance cutoff")
parser.add_argument("-fc", "--first_chain", help="Optional first chain name")
parser.add_argument("-sc", "--second_chain", help="Optional second chain name")

args = parser.parse_args()

######################## Directory stuff ########################

if args.directory:
    # Checks if the directory has a /
    if args.directory[-1] == "/":
        directory += args.directory
    else:
        directory += args.directory + "/"
    os.system("mkdir " + directory)

######################## Parser checking stuff ########################

# Check if the cutoff(s) are numbers
if args.cutoff_max:
    try:
        cutoff_max=float(args.cutoff_max)
    except ValueError:
        print "Oops! That's not a valid number. Please enter a valid number for the maximum cutoff.\n"
        while type(cutoff_max != float):
            try:
                cutoff_max=float(raw_input("Enter the maximum cutoff: \n"))
            except ValueError:
                print "Oops! That's not a valid number. Please enter a valid number for the maximum cutoff.\n"

if args.cutoff_min:
    try:
        cutoff_min=float(args.cutoff_min)
    except ValueError:
        print "Oops! That's not a valid number. Please enter a valid number for the minimum cutoff.\n"
        while type(cutoff_min != float):
            try:
                cutoff_min=float(raw_input("Enter the minimum cutoff: \n"))
            except ValueError:
                print "Oops! That's not a valid number. Please enter a valid number for the minimum cutoff.\n"


# Asks if you want cutoffs, if they aren't specified. Only the maximum cutoff seems relevant
while cutoff_max == False:
    try:
        cutoff_max=raw_input("Enter the maximum cutoff (Na for no maximum cutoff): \n")
        if cutoff_max == "Na":
            cutoff_max = False
            break
        cutoff_max = float(cutoff_max)
    except ValueError:
        print "Oops! That's not a valid number. Please enter a valid number for the maximum cutoff.\n"
        cutoff_max = False

# Minimum cutoff not necessarily useful (may actually make the resulting contact map ambiguous by having both close and far points whitened)
#~ while cutoff_min == False:
    #~ try:
        #~ cutoff_min=raw_input("Enter the minimum cutoff (Enter 'Na' for no minimum cutoff).\nWARNING !! Minimum cutoff might be useless, use it wisely:\n")
        #~ if cutoff_min == "Na":
            #~ cutoff_min = False
            #~ break
        #~ cutoff_min = float(cutoff_min)
    #~ except ValueError:
        #~ print "Oops! That's not a valid number. Please enter a valid number for the minimum cutoff.\n"
        #~ cutoff_min = False

if args.first_chain:
    chainA = args.first_chain

if args.second_chain:
    chainB = args.second_chain

pdb_code = args.output
pdb_filename = args.file

######################## Functions ########################


def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""

    # Gets the vector between two residues
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord

    # Gets the length of that vector
    return numpy.sqrt(numpy.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""

    # First create  an empty (with zeros) 2d matrix
    answer = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)

    # For every atom pair (C-alpha), gets the residue distance and put it into the matrix
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist(residue_one, residue_two)

    return answer

def get_residue_names(chain):
    """Returns a list of residues"""

    # Make an empty list
    res = []

    # Loop through all the residues in the chain and append to the list
    for residue in chain:
        res.append(residue.get_resname())

    return res

def fix_pdb(pdb_code, pdb_filename):
    """Fixes the pdb models by removing CONECT and END statements that biopython doesn't like"""

    global directory

    # Make sure you want to create the pdb_code directory with the right name:
    while pdb_code == "" or pdb_code == None:
        pdb_code = raw_input("What is the wanted output name ? \n")
    pdb_code = pdb_code.rstrip("/")

    # Makes a directory with the pdb_code
    if directory != "":
        os.system("mkdir " + directory + pdb_code)
    else:
        os.system("mkdir " + pdb_code)

    directory += pdb_code + "/"

    # Create the fixed pdb file
    out_name = open(directory + pdb_code + "_fixed.pdb", "w")

    # Now that we have a working output name
    # Try the input name
    try:
        in_name = open(pdb_filename, "r")
    except:
        print "%s does not exist or missing pdb name argument. Did you forget to use the -f option ?\nPlease choose your pdb filename" % (pdb_filename)
        from Tkinter import Tk
        from tkFileDialog import askopenfilename
        while 'in_name' not in locals():
            try:
                Tk().withdraw()
                pdb_filename = askopenfilename()
                in_name = open(pdb_filename, "r")
            except not KeyboardInterrupt:
                print "Still not a valid input pdb file ...\n"

    # Reads everything in the input into memory (a list of each line)
    in_name_content = in_name.readlines()
    in_name.close()

    # For every line, checks if the model number is at the right spot (else biopython won't read it)
    # Also removes the "CONECT" lines, useless for biopython
    for line in in_name_content:
        if line[0:5] == "MODEL":
            #Check the place of the model number
            rank = 6
            for i in line[6:]:
                if i != " ":
                    if rank == 10:
                        break
                    else:
                        out_name.write("MODEL     "+i)
                        break
                rank += 1
        if line[0:6] != "CONECT" and line[0:3] != "END":
            out_name.write(line)
        else:
            pass

    out_name.close()

    return pdb_code, directory + pdb_code + "_fixed.pdb"

def sliced_max(matrix):
    """Only numbers in the 2d matrix, replaces the numbers above the cutoff by the cutoff (max)"""

    # First create  an empty (with zeros) 2d matrix
    tmp = numpy.zeros((len(matrix), len(matrix[0])), numpy.float)

    # Now replaces every value above cutoff_min with the maximum cutoff, which will show white spots
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] > cutoff_max:
                tmp[i][j] += cutoff_max
            else:
                tmp[i][j] += matrix[i][j]

    return tmp

def sliced_min(matrix):
    """Only numbers in the 2d matrix, replaces the numbers below the cutoff by the maximum (to whiten those)"""

    # First create an empty (with zeros) 2d matrix
    tmp = numpy.zeros((len(matrix), len(matrix[0])), numpy.float)

    # If there is no maximum cutoff, find the maximum value (as it will be white, to assign it to the values under cutoff_min)
    if 'cutoff_max' not in locals():
        cutoff_max = 0
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                if matrix[i][j] > cutoff_max:
                    cutoff_max = matrix[i][j]

    # Now replaces every value under cutoff_min with the maximum cutoff, which will show white spots
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] < cutoff_min:
                tmp[i][j] += cutoff_max
            else:
                tmp[i][j] += matrix[i][j]

    return tmp

######################## Main ########################

# In case the pdb comes from preddimer, fixes the pdb for use with biopython
pdb_code, pdb_filename = fix_pdb(pdb_code, pdb_filename)

# Loading the pdb into a biopython module
structure = PDBParser(QUIET=True).get_structure(pdb_code, pdb_filename)

# Loading the biopython module for writing pdbs
io = PDBIO()

# Loop through all the models
for i in range(len(structure)):
    # Selecting the ith model
    model = structure[i]

    io.set_structure(model)

    # Writing a separate model file
    io.save(directory + pdb_code + "_model_" + str(i+1) + ".pdb")

    # Calculates the distance matrix between chain A and B
    dist_matrix = calc_dist_matrix(model["A"], model["B"])

    # Gets the residue names
    try:
        chainA_names = get_residue_names(model[chainA])
    except:
        print "No chain named %s" % (chainA)
        exit()
    try:
        chainB_names = get_residue_names(model[chainB])
    except:
        print "No chain named %s" % (chainB)
        exit()

    # Cutoff handling
    if cutoff_max != False:
        dist_matrix = sliced_max(dist_matrix)
    if cutoff_min != False:
        dist_matrix = sliced_min(dist_matrix)

    # Matplotlib part
    plt.figure(str(i))
    # flipud lets us flip the matrix up and down (as we will reverse the y axis)
    #~ plt.imshow(numpy.flipud(dist_matrix), interpolation='none', cmap=plt.get_cmap('hot'))
    # No need to flip anything here, we'll flip later
    plt.imshow(dist_matrix, interpolation='none', cmap=plt.get_cmap('hot'))
    ax = plt.gca()
    ax.set_xticks(range(len(chainA_names)))
    ax.set_xticklabels([str(x+1) + " " + chainA_names[x] for x in range(len(chainA_names))])
    ax.set_yticks(range(len(chainB_names)))
    # Get the numbering right
    numbers = [y+1 for y in list(reversed(range(len(chainB_names))))]
    # [::-1] reverses a list
    ax.set_yticklabels([str(numbers[x]) + " " + chainB_names[::-1][x] for x in range(len(chainB_names))])
    # Now we're flippin'
    ax.invert_yaxis()
    # Rotates the xticks labels
    plt.xticks(rotation=90)
    plt.title(pdb_code + "_" + str(i+1))
    # Add the colored scale
    plt.colorbar()
    # Saving
    plt.savefig(directory + pdb_code + "_contactmap_" + str(i+1) + ".png", bbox_inches='tight')
