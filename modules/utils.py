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

import csv, re
from cobra.core import Gene,Model
from settings import REACTION_PREFIX


def read_ec_numbers(fname):
    rxn2ec = {row[0]:row[1] for row in csv.reader(open(fname))}
    ECs_rxns = {}
    for rxn,ec in rxn2ec.items():
        if not re.search('^[1-6]\.[0-9][0-9]*\.[0-9][0-9]*',ec):
            continue
        elif ec not in ECs_rxns:
            ECs_rxns[ec] = []
        ECs_rxns[ec].append(rxn)
    
    return ECs_rxns

def add_expantion_fluxes(model,metabolites={},prefix='EFFLUX_',set_as_objective=False,copy_model=True):
    if copy_model:
        model = model.copy()
    
    for metab,coef in metabolites.items():
        if not hasattr(metab, 'id'):
            metab = model.metabolites.get_by_id(metab)
        reaction = cobra.Reaction(prefix+metab.id)
        reaction.add_metabolites({metab : coef})
        if set_as_objective:
            reaction.objective_coefficient = 1
        model.add_reaction(reaction)

    return model

def correct_seed_model(model,metamodel):
    if 'rxn02297' in model.reactions and 'rxn02296' not in model.reactions:
        new_rxn = metamodel.reactions.new_rxn00021.copy()
        old_rxn = model.reactions.rxn02297
        model.add_reaction(new_rxn)
        new_rxn.gene_reaction_rule = old_rxn.gene_reaction_rule
        model.remove_reactions([old_rxn])

    # the following metabolites should be removed because they appear as conserved pools
    # on the 
    biomass = [r for r in model.reactions if r.startswith('bio')][0]
    conflictive_metabolties = ['cpd01997_c','cpd03422_c','cpd11416_c']
    #conflictive_metabolties = ['cpd01997_c','cpd03422_c']
    for m in conflictive_metabolties:
        #print m
        if m not in model.metabolites:
            continue
        metabolite = model.metabolites.get_by_id(m)
        if metabolite not in biomass.products:
            continue
        s_coeff = biomass.get_coefficient(metabolite.id)
        biomass.add_metabolites({metabolite:-s_coeff})
    
    if 'EX_cpd11416_c' in model.reactions:
        model.remove_reactions(['EX_cpd11416_c'])
    
    if 'rxn05029' in model.reactions:
        model.remove_reactions(['rxn05029'])
    
    return model

def prepare_model(model,metamodel,reactions_to_remove=[],correct_seed=False):
    model = model.copy()
    if len(reactions_to_remove) > 0:
        model.remove_reactions(reactions_to_remove)

    if correct_seed:
        correct_seed_model(model,metamodel)

    list_of_reactions = [r for r in model.reactions]
    
    for r in list_of_reactions:
        if r not in metamodel.reactions:
            print "%s not in metamodel %s" % (r.id,metamodel.id)
            continue
        reaction_reference = metamodel.reactions.get_by_id(r.id)
        result = r - reaction_reference
        if len(result.metabolites) > 0:
            genes = r.genes
            gene_rule = r.gene_reaction_rule
            model.remove_reactions([r])
            model.add_reaction(reaction_reference.copy())
            new_reaction = model.reactions.get_by_id(r.id)
            new_reaction._genes = genes
            new_reaction.gene_reaction_rule = gene_rule
            for g in genes:
                g._reaction.add(new_reaction)
        
        model.reactions.get_by_id(r.id).lower_bound = reaction_reference.lower_bound
        
    metabolites_to_remove = [m.id for m in model.metabolites if len(m.reactions) == 0]
    for m in metabolites_to_remove:
        if m not in model.metabolites:
            continue
        model.metabolites.get_by_id(m).remove_from_model()
  
    return model

def create_consisten_model(model,metamodel,consistent_reactions):
    consistent_model = Model()
    consistent_model.id = model.id
    consistent_model.description = model.id
    
    auxiliar_gene = Gene('MODULAR_GAPFILLING')
    auxiliar_gene._model = consistent_model
    consistent_model.genes.append(auxiliar_gene)

    for reaction_id in consistent_reactions:
        new_reaction = metamodel.reactions.get_by_id(reaction_id).copy()
        
        if reaction_id in model.reactions:
            reaction_reference = model.reactions.get_by_id(reaction_id)
            gene_list = []
            
            for gene in reaction_reference.genes:
                if gene.id in consistent_model.genes:
                    gene_list.append(consistent_model.genes.get_by_id(gene.id))
                else:
                    new_gene = Gene(gene.id)
                    new_gene._model = consistent_model
                    consistent_model.genes.append(new_gene)
                    gene_list.append(new_gene)
            
            for gene in gene_list:
                gene._reaction.add(new_reaction)
            new_reaction._genes = gene_list
            new_reaction.gene_reaction_rule = reaction_reference.gene_reaction_rule

        else:
            new_reaction.gene_reaction_rule = auxiliar_gene.name
            auxiliar_gene._reaction.add(new_reaction)

        consistent_model.add_reaction(new_reaction)

    return consistent_model
    
def get_full_coupled_sets(reaction_names,fctab, exclude_preffix=None):

    """ Interpretation for element (i, j):
          1 - fully coupled <=>
          2 - partially coupled <->
          3 - reaction i is directionally coupled to j ( v[i]<>0 -> v[j]<>0 )
          4 - reaction j is directionally coupled to i ( v[j]<>0 -> v[i]<>0 )
          5 - uncoupled
    """

    assert fctab.shape[0] == len(reaction_names)
    
    already_coupled = set()
    coupling_sets = []
    for i in np.arange(fctab.shape[0]):
        if i in already_coupled:
            continue
        indexes = np.where(fctab[i,:]==1)[0]
        if len(indexes) < 2:
            continue
        coupling_sets.append(indexes)
        already_coupled = already_coupled.union(indexes)

    #coupling_sets = np.array([np.array(b) for b in set([tuple(a) for a in fctab])])
    #coupling_sets = [subset for subset in [np.where(es==1)[0] for es in coupling_sets] if len(subset)>1]
    result = {}
    counter = 1
    for subset in coupling_sets:
        rs_id = 'RS_'+str(counter)
        if exclude_preffix:
            reaction_ids = [reaction_names[i] for i in subset if not re.search(exclude_preffix,reaction_names[i])]
        else:
            reaction_ids = [reaction_names[i] for i in subset]
        
        if len(reaction_ids) > 1:
            result[rs_id] = reaction_ids
            counter += 1
        
    return result

def decorate_graph(G,labels={},colors={}):

    for tail,head in G.edges():
        G.edge[tail][head]['graphics'] = {}

        G.edge[tail][head]['graphics']['targetArrow'] = "standard"

        try:
            width = float(G[tail][head]['label'])
            G.edge[tail][head]['graphics']['width'] = width
        except:
            G.edge[tail][head]['graphics']['width'] = 1.0
            

    for n in G.nodes():
        label = n
        if n in labels:
            G.node[n]['label'] = labels[n]
            label = labels[n]

        graphics = {}

        if n in colors:
            color = colors[n]
        else:
            color = None
            

        if G.node[n]['node_class'] == 'reaction':
            outline = "#000000"
            if not color:
                color = "#c0c0c0"
            height = 16.0
            width = max((len(label) * 8.0),85.0)

            
                            
            graphics = {"w":width, "h":height, "type":"roundrectangle", "fill":color, "outline":outline}

        elif G.node[n]['node_class'] == 'metabolite':
            outline = "#ffffff"
            if not color:
                color = "#ffffff"
            height = 15.0
            width = max((len(label) * 8.0),60.0)
            if n in colors:
                color = colors[n]
                outline = "#000000"
               
            graphics = {"w":width, "h":height, "type":"rectangle", "fill":color, "outline":outline}


        G.node[n]['graphics'] = graphics

    return G

def csv_save(a_list,fname):
    if not isinstance(a_list,list):
        return
    elif len(a_list)<1:
        return
    elif not isinstance(a_list[0],list):
        a_list = [[e] for e in a_list]
    f = open(fname,'w')
    w = csv.writer(f)
    x = w.writerows(a_list)
    f.close()
    return x



def f_rxn(x):
    return re.search(REACTION_PREFIX,x) 

def f_ex(x):
    return re.search(EXCHANGE_PREFIX,x)

def f_flux(x):
    return f_ex(x) or f_rxn(x) 

