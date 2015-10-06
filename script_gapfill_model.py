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

import sys, csv, re, os, time
import json
from cobra.flux_analysis import find_blocked_reactions
from cobra.io import read_sbml_model, write_sbml_model
from cobra_extensions.model_reconfiguration import add_expantion_fluxes
from modules import utils
from modules.fastcore import fast_find_blocked, fast_core
#import warnings
#warnings.filterwarnings("ignore")

######################################################
# GLOBAL PARAMETERS                                  #
######################################################

METAMODEL_PATH = "Data/MM197x_consistent.xml"
EXPAND_METAMODEL = False
RXN2ECS_PATH = "Data/rxn2ec.csv"
OUTPUT_FOLDER = "./"
OUTPUT_SUFIX = "gapfilled"
REACTION_PREFIX = "rxn|new"
BIOMASS_PREFIX = "bio"
EXCHANGE_PREFIX = "EX|DM|EFF"
TRANSPORT_WEIGHT = 50
REACTION_WEIGHT = 10
ADD_TRANSPORTS = False
EPSILON = 1e-4
DEBUG = False
SOLVER = ""
SPONTANEOUS = "SPONTANEOUS"

def load_parameters(fname):
    param_dict = {}
    with open(fname) as f:
        param_dict = json.load(f)

    global ADD_TRANSPORTS
    global BIOMASS_PREFIX
    global EXCHANGE_PREFIX
    global METAMODEL_PATH
    global EXPAND_METAMODEL
    global OUTPUT_FOLDER
    global OUTPUT_SUFIX
    global REACTION_PREFIX
    global REACTION_WEIGHT
    global RXN2ECS_PATH
    global TRANSPORT_WEIGHT
    global EPSILON
    global DEBUG
    global SOLVER
    global SPONTANEOUS
    
    METAMODEL_PATH = param_dict["metamodel_path"].encode()
    RXN2ECS_PATH = param_dict["rxn2ecs_path"].encode()
    
    OUTPUT_FOLDER = param_dict["output_folder"].encode()
    OUTPUT_SUFIX = param_dict["output_model_suffix"].encode()
    
    EXPAND_METAMODEL = bool(param_dict["expand_metamodel"])
    REACTION_PREFIX = param_dict["reaction_prefix"].encode() 
    BIOMASS_PREFIX = param_dict["biomass_prefix"].encode() 
    EXCHANGE_PREFIX = param_dict["exchange_prefix"].encode()
    
    REACTION_WEIGHT = param_dict["reaction_weight"]
    TRANSPORT_WEIGHT = param_dict["transport_weight"]
    EPSILON_WEIGHT = param_dict["epsilon"]
    DEBUG = bool(param_dict["debug"])
    ADD_TRANSPORTS = param_dict["add_transporters"]

    SOLVER = bool(param_dict["solver"])
    SPONTANEOUS = param_dict["spontaneous"]

def f_rxn(x):
    return re.search(REACTION_PREFIX,x) 

def f_ex(x):
    return re.search(EXCHANGE_PREFIX,x)

def f_flux(x):
    return f_ex(x) or f_rxn(x) 



######################################################
# INIT REQUIRED DATA STRUCTURES                      #
######################################################

param_fname = './parameters.json'
try:
    load_parameters(param_fname)
    print "Parameters loaded from %s " % param_fname
except:
    print "The parameter file %s not found, running with defaul parameter" % param_fname
    


print "Reading Metamodel", METAMODEL_PATH,
metamodel = read_sbml_model(METAMODEL_PATH)
print " - loaded!"

ec = '^[1-6]\.[0-9][0-9]*\.[0-9][0-9]*'
ECs_rxns = utils.read_ec_numbers(RXN2ECS_PATH)
rxn2ec = {r.id:r.annotation['ec_number'] for r in metamodel.reactions 
                            if re.search(ec,r.annotation['ec_number'])}


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
    print "Model %s loaded" % model.id
except:
    print "Invalid SBML file"
    sys.exit(1)


initial_reactions = {r.id for r in model.reactions}

for r in model.reactions:
    if r.upper_bound > 0:
        r.upper_bound = 999999.0
    if r.lower_bound < 0:
        r.lower_bound = -999999.0

print "Finding blocked reactions"
blocked_reactions = find_blocked_reactions(model)

if len(blocked_reactions) == 0:
    print "No blocked reactions found, the model is consistent, nothing to be done"
    sys.exit(0)
else:
    print "%i blocked reactions found" % len(blocked_reactions)


print "Model prepared for Gapfilling"
metablocked_reactions =  [r for r in blocked_reactions if r not in metamodel.reactions]
model = utils.prepare_model(model,metamodel,reactions_to_remove=metablocked_reactions,correct_seed=True)


biomass = [r for r in model.reactions if r.startswith(BIOMASS_PREFIX)][0]
metamodel.add_reactions([biomass.copy()])
for r in metamodel.reactions:
    r.objective_coefficient = 0

r = metamodel.reactions.get_by_id(biomass.id)
r.objective_coefficient = 1
metamodel.optimize()
if metamodel.solution.f < 1e-7:
    print "ERROR: metamodel %s cannot produce biomass %s from model %s" % (metamodel.id,biomass.id,model.id)
    sys.exit(0)


if EXPAND_METAMODEL:
    new_blocked_reactions = fast_find_blocked(model)
    metablocked_reactions += [r for r in new_blocked_reactions if r not in metamodel.reactions]
    model.remove_reactions([r for r in new_blocked_reactions if r not in metamodel.reactions])
    metamodel.add_reactions([r.copy() for r in model.reactions if r.id not in metamodel.reactions])
else:
    metablocked_reactions += [r.id for r in model.reactions if r not in metamodel.reactions]
    model.remove_reactions([r for r in metablocked_reactions if r in model.reactions])

print "%i metablocked reactions found" % len(metablocked_reactions)
print "Removing the following metablocked reaction from the model:"
print "==============================================================="
print " ".join(metablocked_reactions)
print "==============================================================="



######################################################
# INIT GAPFILLING                                    #
######################################################

if ADD_TRANSPORTS:
    # Exclude metabolites which already have an exchange flux
    efflux_metabolites = [r.metabolites.keys()[0].id for r in metamodel.reactions if f_Ex(r.id)]
    metabolites = [m.id for m in metamodel.metabolites if m.endswith('_c')]
    metabolites_to_expand = set(metabolites) - set(efflux_metabolites)

    SUX = utils.add_expantion_fluxes(metamodel,metabolites={metamodel.metabolites.get_by_id(m):-1 for m in metabolites_to_expand}, prefix='X_')
    for r in SUX.reactions:
        if 'gene_locus_tag' not in r.annotation:
            r.annotation['gene_locus_tag'] = ''
        if not r.startswith('X_'):
            continue
        r.upper_bound = 100
        r.lower_bound = -100
else:
    SUX = metamodel

model_reactions = [r.id for r in model.reactions]
consistent_reactions = list(set(model_reactions) - set(blocked_reactions))


model_ECs = {rxn2ec[r.id] for r in model.reactions if r.id in rxn2ec}
unweigthed_rxn = []
if SPONTANEOUS in metamodel.genes:
    spontaneous_rxns = metamodel.genes.get_by_id(SPONTANEOUS).reactions
    unweigthed_rxn += [r.id for r in spontaneous_rxns]

unweigthed_rxn += [r for ec in model_ECs if ec in ECs_rxns for r in ECs_rxns[ec] if r not in model_reactions]
unweigthed_rxn = set(unweigthed_rxn)

weights = {}
for r in SUX.reactions:
    weights[r.id] = 0
    if r.id in unweigthed_rxn:
        weights[r.id] = 0
    elif r.startswith('X_'):
        weights[r.id] = TRANSPORT_WEIGHT
    elif f_rxn(r.id):
        weights[r.id] = REACTION_WEIGHT



print len(blocked_reactions),"Blocked reactions to solve"
print "Running FastCore on",model.id


computing_time = -1

for r in SUX.reactions:
    if r.upper_bound > 0:
        r.upper_bound = 999999.0
    if r.lower_bound < 0:
        r.lower_bound = -999999.0
try:
    start = time.clock()
    consistent_reactions = fast_core( 
                                       SUX,
                                       model_reactions,
                                       weights=weights,
                                       epsilon=EPSILON,
                                       debug = DEBUG
                                    )  
    computing_time = (time.clock() - start)
except Exception, e:
    print "ERROR: while runnning fastcore on",model.id
    print "Exception msg:",str(e)
    sys.exit(1)


print "FASTCORE COMPUTING TIME: %s SECONDS" % computing_time
gapfilling_reactions = set(consistent_reactions) - set(model_reactions)
print "FastCore Finish Succsesfully!",len(gapfilling_reactions),"gapfilling reactions added to model"


consistent_model = utils.create_consisten_model(model,metamodel,consistent_reactions)

new_blocked_reactions = find_blocked_reactions(consistent_model)

if len(new_blocked_reactions) > 0:
    print "ERROR: %i blocked reactions found in %s after running Fastcore" % (len(new_blocked_reactions),model.id)
    sys.exit(1)

final_model_reactions = {r.id for r in consistent_model.reactions}

sbml_out = '.'.join([model.id,OUTPUT_SUFIX,"xml"])
sbml_out = os.path.join(OUTPUT_FOLDER,sbml_out)
write_sbml_model(consistent_model,sbml_out,use_fbc_package=False)

print "Curated model %s saved as %s" % (model.id,sbml_out)
print "====================================================="
print "MODEL_ID\tRi\tBlk\tMetaBlk\tGF\tRf"
print "%s\t%i\t%i\t%i\t%i\t%i" % (  consistent_model.id,
                                    len(initial_reactions),
                                    len(blocked_reactions),
                                    len(initial_reactions - final_model_reactions),
                                    len(gapfilling_reactions),
                                    len(final_model_reactions)
                                )

print "====================================================="


    

