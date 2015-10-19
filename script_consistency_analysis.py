#!/usr/bin/python
# Python script to perform consistency analysis on metabolic models
# Copyright (C) 2015  Miguel Ponce de Leon
# Contact: miguelponcedeleon@gmail.com

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

import sys,re,os
import csv
from cobra.io import read_sbml_model
from modules.fastcore import fast_find_blocked
from modules import utils
import modules.modular_gap_find as mgf
import networkx as nx

def usage():
    print "USAGE: %s <sbml_file> [output_folde]" % sys.argv[0]

def save_gap_graph_gml(G,fname,compounds_table={},reaction_table={}):
    cmps = [n for n in G.nodes() if G.node[n]['node_class']=='metabolite']
    rxns = [n for n in G.nodes() if G.node[n]['node_class']=='reactions']
    labels = {}
    for n in cmps:
        if n not in compounds_table:
            labels[n] = n
            continue
        labels[n] = compounds_table[n]['abbrev']

    for n in rxns:
        label = ''
        if n not in reaction_table:
            labels[n] = n
            continue
        if re.search('[0-9]\.[0-9]',reaction_table[n]['enzyme']):
            EC = reaction_table[n]['enzyme']
            EC = EC.split("|")[0]
            labels[n] = reaction_table[n]['enzyme']
        else:
            labels[n] = reaction_table[n]['abbrev']

    G = utils.decorate_graph(G,labels=labels)
    nx.write_gml(G,fname)

def read_seed_table(fname,delimiter="\t"):
    reader = csv.reader(open(fname),delimiter=delimiter)
    table_dict = {}
    tab_header = reader.next()
    _id = tab_header.pop(0)
    n = len(tab_header)
    for row in reader:
        item_id = row.pop(0)
        table_dict[item_id] = {tab_header[i]:row[i] for i in range(n)}

    return table_dict



project_root = "/home/mponce/Research/Data/GSR/THE_SEED_MODELS/Large_Scale_Analysis_of_UM"
models_folder = project_root+"/Models/*.xml"

reactions_tab_file = "Data/" + "seed_reactions_tab.csv"
compounds_tab_file = "Data/" + "seed_compounds_tab.csv"

reaction_table = read_seed_table(reactions_tab_file)
compounds_table = read_seed_table(compounds_tab_file)

######################################################
# PREPARE THE MODEL                                  #
######################################################

if len(sys.argv) < 2:
    print "Error: no SBML file name supplied"
    sys.exit(1)
else:
    model_file = sys.argv[1]

try:
    model = read_sbml_model(model_file)
except:
    print "Invalid SBML for file %s" % model_file
    usage()
    sys.exit(1)
    
    
if len(sys.argv) >= 3:
    output_folder = sys.argv[2]
    if not os.path.isdir(output_folder):
        print "Invalid output folder %s" % output_folder
        usage()
        sys.exit(1)
else:
    output_folder = "./"


print "========================================================"
print "Model %s loaded succesfully" % model.id
print
print "Begining consistency analysis on model: %s" % model.id
print "Finding blocked reactions..."
blocked_reactions = fast_find_blocked(model)
print "Finding gap metabolites..."
gap_metabs = mgf.find_gap_metabolites(model,blocked_reactions=blocked_reactions)
print "Computing Unconnected Modules..."
unnconected_modules = mgf.find_unconnected_modules(model,blocked_reactions=blocked_reactions)
print "Consistency analysis finished"

print "========================================================"
print "Consistency report:"
print "- Total blocked reactions: %i" % len(blocked_reactions)
print "- Total gap metabolites %i" % len(blocked_reactions)
print "- Total Unconected Modules %i" % len(unnconected_modules)
print "========================================================"

blocked_reactions_fname = os.path.join(output_folder,model.id + "_BLK_RXNS.csv")
gaps_fname = os.path.join(output_folder,model.id + "_GAPS_METAB.csv")
ums_fname = os.path.join(output_folder,model.id + "_UMS.csv")
gap_graph_fname = os.path.join(output_folder,model.id + "_Gap_graph.gml")

blocked_reactions.sort()
print "Save blocked reactions as %s" % blocked_reactions_fname
utils.csv_save(blocked_reactions,blocked_reactions_fname)

gap_metabs.sort()
print "Save blocked reactions as %s" % gaps_fname
utils.csv_save(gap_metabs,gaps_fname)

print "Saving unconnected modules as %s" % ums_fname
um_ids = sorted(unnconected_modules.keys(),key=lambda x: int(x[2:]))
UMs = [k+": "+" ".join(unnconected_modules[k]) for k in um_ids] 
utils.csv_save(UMs,ums_fname)

print "Saving gap graph as %s" % gap_graph_fname
G = mgf.create_gap_graph(model,gap_metabs,blocked_reactions)
save_gap_graph_gml( G,gap_graph_fname,
                    compounds_table=compounds_table,
                    reaction_table=reaction_table
                )

print "========================================================\n"



