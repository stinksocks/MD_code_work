import sys
import random
import math
import MDAnalysis as mda

# WELCOME TO CG BILAYER
# To use this program please invoke the following command: python CG_bilayer.py file.txt
# In file.txt should be a list of the lipids that are present in a single leaflet of the bilayer you are trying to build

"""EXAMPLE .TXT FOLLOWS:
CHL     800
PC_16:0_16:0   800
"""

#Based on this .txt file you will create a randomised 50:50 ratio Cholestrol and PC 16:0/16:0 bilayer, where there are 1600 lipids present in each bilayer#

#This code reads the .txt file to get the lipid types and how many you want.
def structure(inlist):
    fil = []
    with open(inlist, "r") as infile:
        for line in infile:
            line = line.strip()
            columns = line.split()
            fil.append(columns)

    coords, copied = [], []

#num lipids should always square root into an integer, i.e. 100, 10000, and so on. You will provide this value interactively at the command line

    num_lipids = int(raw_input("How many lipids per bilayer?: "))
    grid_size = int(math.sqrt(num_lipids))
    running_total = 0
    string = ""

#Here a square grid of coordinates is created for lipid placement

    for x in range(0, grid_size):
        for y in range(0, grid_size):
            string = str(x) + ":" + str(y)
            coords.append(string)

#These coordinates are shuffled to allow for random lipid placement

    random.shuffle(coords)

#Lipids are inserted in this code using a library.

#IMPORTANT THIS LIBRARY MUST CONTAIN A STRUCTURES DIRECTORY AND A TOPOLOGY DIRECTORY, WHERE THE .gro STRUCTURE FILE and the .itp TOPOLOGY FILE have the same prefix before the file type. Example, CHL.gro and CHL.itp for the cholestrol structure and topology.

#Here the directory is labelled as lib/structures in the same directory as the script, however you are free to edit this to provide your own path to a library.

#This code makes a single leaflet

    for item in fil:
        for insert in range(int(item[1])):
            val = coords.pop()
            xy = val.split(":")
            x = int(xy[0])
            y = int(xy[1])
            u = mda.Universe("../lib/annealed_structures/fixed/processed/"+item[0]+".gro")
            atoms = u.select_atoms("all")
            u2 = u.copy()
            move_by = [(x*7), (y*7), 0]
            u2.atoms.translate(move_by)
            copied.append(u2.atoms)

    new_u = mda.Merge(*copied)
    output = new_u.select_atoms("all")
    ps = output.positions           

    #This code creates the box size

    def box_value(inputlist, pos):
        return max([sublist[pos] for sublist in inputlist]) - min([sublist[pos] for sublist in inputlist])

    x_box, y_box, z_box = box_value(ps, 0), box_value(ps, 1), box_value(ps, 2)

    #This code finds the centre of mass of the leaflet and places it in the middle of the box

    u_2_dim = [x_box, y_box, z_box*3, 90, 90, 90]

    new_u.dimensions = u_2_dim

    box_center = [x_box/2, y_box/2, z_box*1.5]
    com = output.center_of_mass()
    center = (box_center-com)

    new_u.atoms.translate(center)

    leaflet2 = new_u.copy()
    box_center2 = [x_box/2, y_box/2, (z_box*1.5)-20]

    #This code creates a clone of the lipid leaflet and flips it, arranging such that a bilayer is formed.

    output2 = leaflet2.select_atoms("all")
    com2 = output2.center_of_mass()
    center2 = (box_center2-com2)
    leaflet2.atoms.rotateby(180, [1, 0, 0])
    leaflet2.atoms.translate(center2)
    output2 = leaflet2.select_atoms("all")

    layers = []
    layers.append(new_u.atoms)
    layers.append(leaflet2.atoms)

    bilayer = mda.Merge(*layers)
    bilayer.dimensions = u_2_dim
    bl = bilayer.select_atoms("all")

    #AT THIS STAGE THE BILAYER IS WRITTEN OUT AS A .GRO FILE, THIS STRUCTURE FILE IS THE START POINT FOR EQUILIBRATION

    bl.write("bilayer.gro")



#This code writes the topology
def topology(inputlst):
    file_lst, top_lst = [], []
    with open(inputlst, "r") as infile:
        for line in infile:
            line = line.strip()
            file_lst.append(line)


    #The file template.top, or an equivalent blank topology can be used, where template.top is included as a reference

    with open("template.top", "r") as infile2:
        for line in infile2:
            line = line.strip()
            top_lst.append(line)


    include_lst = []
    res_names = []
    sc = 0


    #This code uses the names of your topolgies in you topology library to insert the name of the molecule into the .top file

    for x in file_lst:
        columns = x.split()
        name = columns[0]
        num = columns[1]
        include = '#include "lib/topologies/' + name + '.itp"'
        include_lst.append(include)
        with open("lib/topologies/"+name+".itp", "r") as itpfile:
            sc = 0
            for line in itpfile:
                line = line.strip()
                columns = line.split()
                if line.startswith("["):
                    sc += 1
                if sc == 1 and "moleculetype" not in line and not line.startswith(";"):
                    res = line[:4]
                    if len(res) > 0:
                        res_names.append(res + "  " + num)

    #This includes the martini force field
    ind = top_lst.index('#include "../martini2_2.itp"')
    for x in include_lst:
        top_lst.insert((ind+1), x)
        #top_lst.insert((ind+2), '#include "str.itp"')

    top_lst = top_lst + res_names*2

    #This produces your .top file used for minimization and equilibration etc, make sure you use it in the right place to avoid missing .itp errors

    with open("auto_bilayer.top", "w") as outfile:
            for item in top_lst:
                print >> outfile, item

structure(sys.argv[1])
topology(sys.argv[1])
#After this script runs you will have a bilayer.gro and auto_bilayer.top files, from these you can run your equilibration and production.
