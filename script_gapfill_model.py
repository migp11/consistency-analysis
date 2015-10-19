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


import sys, csv, re, os, time
import json
from cobra.flux_analysis import find_blocked_reactions
from cobra.io import read_sbml_model, write_sbml_model
from cobra_extensions.model_reconfiguration import add_expantion_fluxes
from modules import utils,settings
from modules.fastcore import fast_find_blocked, fast_core

#import warnings
#warnings.filterwarnings("ignore")


######################################################
# INIT REQUIRED DATA STRUCTURES                      #
######################################################

if len(sys.argv) < 2:
    print "Error: no SBML file name supplied"
    sys.exit(1)
else:
    model_file = sys.argv[1]

if len(sys.argv) >= 3:
    param_fname = sys.argv[2]
else:
    param_fname = './parameters.json'

try:
    settings.load_parameters(param_fname)
    print "Parameters loaded from %s " % param_fname
except Exception, e:
    print "The parameter file %s not found, running with defaul parameter" % param_fname
    print str(e)
    sys.exit(0)


print "Reading Metamodel", settings.METAMODEL_PATH,
metamodel = read_sbml_model(settings.METAMODEL_PATH)
print " - loaded!"

ec = '^[1-6]\.[0-9][0-9]*\.[0-9][0-9]*'
ECs_rxns = utils.read_ec_numbers(settings.RXN2ECS_PATH)
rxn2ec = {r.id:r.annotation['ec_number'] for r in metamodel.reactions 
                            if re.search(ec,r.annotation['ec_number'])}


######################################################
# PREPARE THE MODEL                                  #
######################################################

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
metablocked_reactions =  [r for r in blocked_reactions if r not in metamodel.reactions ]
model = utils.prepare_model(model,metamodel,reactions_to_remove=metablocked_reactions,correct_seed=True)


biomass = [r for r in model.reactions if r.startswith(settings.BIOMASS_PREFIX)][0]
metamodel.add_reactions([biomass.copy()])
for r in metamodel.reactions:
    r.objective_coefficient = 0

r = metamodel.reactions.get_by_id(biomass.id)
r.objective_coefficient = 1
metamodel.optimize()
if metamodel.solution.f < 1e-7:
    print "ERROR: metamodel %s cannot produce biomass %s from model %s" % (metamodel.id,biomass.id,model.id)
    sys.exit(0)


if settings.EXPAND_METAMODEL:
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

if settings.ADD_TRANSPORTS:
    # Exclude metabolites which already have an exchange flux
    efflux_metabolites = [r.metabolites.keys()[0].id for r in metamodel.reactions if utils.f_Ex(r.id)]
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
if settings.SPONTANEOUS in metamodel.genes:
    spontaneous_rxns = metamodel.genes.get_by_id(settings.SPONTANEOUS).reactions
    unweigthed_rxn += [r.id for r in spontaneous_rxns]

unweigthed_rxn += [r for ec in model_ECs if ec in ECs_rxns for r in ECs_rxns[ec] if r not in model_reactions]
unweigthed_rxn = set(unweigthed_rxn)

weights = {}
for r in SUX.reactions:
    weights[r.id] = 0
    if r.id in unweigthed_rxn:
        weights[r.id] = 0
    if r in model_reactions or r.startswith('EFF') or r.startswith('EX'):
        weights[r.id] = 0
    elif r.startswith('X_'):
        weights[r.id] = settings.TRANSPORT_WEIGHT
    elif utils.f_rxn(r.id):
        weights[r.id] = settings.REACTION_WEIGHT
        


print len([r for r in blocked_reactions if r in model_reactions]),"Blocked reactions to solve"
print "Running FastCore on %s solver=%s using milp=%s" % (model.id,settings.SOLVER,settings.USE_MILP)



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
                                       epsilon=settings.EPSILON,
                                       debug = settings.DEBUG,
                                       solver=settings.SOLVER,
                                       use_milp=settings.USE_MILP
                                    )  
    computing_time = (time.clock() - start)
except Exception, e:
    print "ERROR: while runnning fastcore on",model.id
    print "Exception msg:",str(e)
    sys.exit(1)


blocked_curated = set(consistent_reactions).intersection(blocked_reactions)
print "FASTCORE COMPUTING TIME: %s SECONDS" % computing_time
gapfilling_reactions = set(consistent_reactions) - set(model_reactions)
print "FastCore Finish Succsesfully!",len(gapfilling_reactions),"gapfilling reactions added to model"

consistent_model = utils.create_consisten_model(model,metamodel,consistent_reactions)
new_blocked_reactions = find_blocked_reactions(consistent_model)
if len(new_blocked_reactions) > 0:
    print "ERROR: %i blocked reactions found in %s after running Fastcore" % (len(new_blocked_reactions),model.id)
    sys.exit(1)

clean = True
if clean:
    aux = consistent_model.copy()
    biomass = aux.reactions.get_by_id(biomass.id)
    aux.remove_reactions([r for r in aux.reactions if r.startswith('EFF') and r._metabolites.keys()[0] in biomass.metabolites])
    new_blocked_reactions = find_blocked_reactions(aux)
    if len(new_blocked_reactions) == 0:
        print "Cleaning model"
        consistent_model = aux
        gapfilling_reactions = [r for r in gapfilling_reactions if r in consistent_model.reactions]

final_model_reactions = {r.id for r in consistent_model.reactions}

sbml_out = '.'.join([model.id,settings.OUTPUT_SUFIX,"xml"])
sbml_out = os.path.join(settings.OUTPUT_FOLDER,sbml_out)
write_sbml_model(consistent_model,sbml_out,use_fbc_package=False)

print "Curated model %s saved as %s" % (model.id,sbml_out)
print "====================================================="
print "MODEL_ID\tRXN_initial\tBLK_curated\tBLK_Removed\tGF\tRXN_final"
print "%s\t%i\t%i\t%i\t%i\t%i" % (  consistent_model.id,
                                    len(initial_reactions),
                                    len(blocked_curated),
                                    len(initial_reactions - final_model_reactions),
                                    len(gapfilling_reactions),
                                    len(final_model_reactions)
                                )

print "====================================================="


    

