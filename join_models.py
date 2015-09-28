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

import sys
import re
import cobra
from cobra.io import read_sbml_model



def create_reactions_dictionary(model):
    reactions_dict = {}
    for r in model.reactions:
        reaction_dict = dict([(m.id,r.get_coefficient(m)) for m in r.get_reactants()] +[(m.id,r.get_coefficient(m)) for m in r.get_products()])
        gene_rule = r.gene_reaction_rule
        reactions_dict[r.id] = { 
                                  'id':r.id,
                                  'name':r.name,
                                  'lower_bound':r.lower_bound,
                                  'upper_bound':r.upper_bound,
                                  'reversibility':r.reversibility,
                                  'gene_rule':gene_rule,
                                  'reaction_dict':reaction_dict
                                }
    
    return reactions_dict

def create_metabolites_dictionary(model):
    metabolites = { m.id:
                        { 
                          'id': m.id,
                          'name':m.name,
                          'charge':m.charge,
                          'compartment':m.compartment
                        } 
                    for m in model.metabolites 
                  }
    
    return metabolites
    
def create_metabolite_from_dict(m_dict):
    metabolite = cobra.Metabolite(m_dict['id'])
    metabolite.name = m_dict['name']
    metabolite.charge = m_dict['charge']
    metabolite.compartment = m_dict['compartment']
    return metabolite

def create_reaction_from_dict(r_dict,metabolites_dict):
    reaction = cobra.Reaction(r_dict['id'])
    reaction.name = r_dict['name']
    reaction.lower_bound = r_dict['lower_bound']
    reaction.upper_bound = r_dict['upper_bound']
    reaction.reversibility = r_dict['reversibility']
    for m,coeff in r_dict['reaction_dict'].items():
        reaction.add_metabolites({metabolites_dict[m]:coeff})
    
    return reaction

def create_metamodel(list_of_models,model_id='metamodel'):
 
    metamodel = cobra.Model(model_id)
    
    for model in list_of_models:
        metabolites_dictionary = create_metabolites_dictionary(model)
        for m_id,m_dict in metabolites_dictionary.items():
            if m_id in metamodel.metabolites:
                continue
            metabolite = create_metabolite_from_dict(m_dict)
            metamodel.add_metabolites([metabolite])
            
    all_metabolites_dict = {m.id:m for m in metamodel.metabolites}
            
    for model in list_of_models:
        reactions_dictionary = create_reactions_dictionary(model)
        for r_id,r_dict in reactions_dictionary.items():
            if r_id in metamodel.reactions:
                r = metamodel.reactions.get_by_id(r_id)
                r.lower_bound = min(r_dict['lower_bound'],r.lower_bound)
                r.upper_bound = max(r_dict['upper_bound'],r.upper_bound)
                r.reversibility = max(r_dict['reversibility'],r.reversibility)
                continue
                
            reaction = create_reaction_from_dict(r_dict,all_metabolites_dict)
            metamodel.add_reaction(reaction)
    
    
    return metamodel


def main():
    metamodel=create_metamodel(list_of_models)
    G=create_metabolic_graph(metamodel)    
    gap_graph=G.subgraph(metabolites_always_gap+reactions_always_blocked.tolist())
    i = 1
    for um_graph in nx.connected_component_subgraphs(gap_graph.to_undirected()):
        um_graph = decorate_graph(um_graph,labels=nodes_label_dict)
        nx.write_gml(um_graph,'Metamodel_UM'+str(i)+'.gml')
        i += 1
    
    

reactions_copies = []
rxns_subsys = {}
for model in models:

    for r in model.reactions:
        if r.id not in rxns_subsys:
            rxns_subsys[r.id] = []
            reactions_copies.append(r.copy())

        rxns_subsys[r.id].append(r.subsystem)

RXN2SUBSYS = {}
for k,v in rxns_subsys.items():
    RXN2SUBSYS[k] = Counter(v).most_common(1)[0][0]
    if len(RXN2SUBSYS[k]) < 3:
        RXN2SUBSYS[k] = 'Unclassified'
