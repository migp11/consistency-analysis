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

import re
import csv
import numpy as np
import networkx as nx
from cobra import Model, Reaction, Metabolite
from cobra.core import Object
from cobra.solvers import solver_dict, get_solver_name
from cobra.flux_analysis.variability import find_blocked_reactions




def fast_core( model,C,epsilon=1e-4, zero_tolerance=1e-9, 
               solver=None, weights={}, debug=False, 
               allow_dilution=False, use_milp=False):

    model = model.copy()
    if allow_dilution:
        for m in model.metabolites:
            m._constraint_sense = 'G'

    if hasattr(C[0], 'id'):
        C = [r.id for r in C]

    C = set(C)

    all_reactions = {r.id for r in model.reactions}
    irrev_reactions = {r.id for r in model.reactions if r.lower_bound >= 0}
    rev_reactions = {r.id for r in model.reactions if r.lower_bound < 0}

    A = []
    flipped = False
    singleton = False


    J = C.intersection(irrev_reactions)
    P = all_reactions - C
    if len(weights) == 0:
        weights = {r:10.0 for r in P}

    Supp = find_sparse_mode(model,J,P,singleton,weights=weights,
                            epsilon=epsilon,solver=solver,use_milp=use_milp)
    

    if len(J - set(Supp)) > 0:
        print J - set(Supp)
        raise Exception("Error: Inconsistent irreversible core reactions")

    A = Supp
    J = C - A
    if debug:
        print "|A|=",len(A)
        print "|J|=",len(J)

    while len(J) > 0:
        P = P - A
        if singleton:
            Supp = find_sparse_mode(model,[list(J)[0]],P,singleton,weights=weights,epsilon=epsilon,solver=solver)
        else:
            Supp = find_sparse_mode(model,J,P,singleton,weights=weights,epsilon=epsilon,solver=solver)

        A = A.union(Supp)
        if debug:
            print "|A|=",len(A)

        if len(J.intersection(A)) > 0:
            J = J - A
            if debug:
                print "|J|=",len(J)
            flipped = False
            singleton = False
        else:
            if singleton:
                JiRev = {list(J)[0]} - irrev_reactions
            else:
                JiRev = J - irrev_reactions
            
            if flipped or len(JiRev) == 0:
                if singleton:
                    print J
                if singleton:
                    raise Exception("Global network not consisten")
                else:
                    flipped = False
                    singleton = True
            else:
                if debug:
                    print "Flipp"
                flipped = True
                for r in JiRev:
                    reaction = model.reactions.get_by_id(r)
                    reaction_dict = dict([(k, -2*v) for k, v in reaction._metabolites.items()])
                    reaction.add_metabolites(reaction_dict)
                    lb = reaction.lower_bound
                    reaction.lower_bound = -reaction.upper_bound
                    reaction.upper_bound = -lb

    return A

def find_sparse_mode( model,J,P,singleton, weights={}, epsilon=1e-4, 
                      zero_tolerance=1e-9, solver=None, use_milp=False):
    
    if len(J) == 0:
        return []
    solution_dict = LP7(model,J,epsilon=epsilon,solver=solver)

    Supp = {k for k,v in solution_dict.items() if (abs(v) - 0.99 * epsilon) > zero_tolerance}
    K = set(J).intersection(Supp)

    if len(K) == 0:
        return []
    
    solution_dict = LP9(model,K,P,weights=weights,epsilon=1.0,solver=solver,use_milp=use_milp)
    Supp = {k for k,v in solution_dict.items() if (abs(v) - 0.99 * epsilon) > zero_tolerance}

    return Supp
   
def fast_find_blocked(model,the_reactions=None,epsilon=1e-3,zero_tolerance=1e-7):

    for r in model.reactions:
        r.objective_coefficient = 0

    if the_reactions is None:
        the_reactions = model.reactions

    if hasattr(the_reactions[0], 'id'):
        the_reactions = [r.id for r in the_reactions]
    
    consistent = fast_consistency_check(
                                        model,
                                        the_reactions=the_reactions,
                                        epsilon=epsilon,zero_tolerance=zero_tolerance
                                       )
    
    blocked_reactions = set(the_reactions) - set(consistent)

    return list(blocked_reactions)

def fast_consistency_check(model,the_reactions=None, epsilon=1e-4, zero_tolerance=1e-7, debug=False):

    solver = solver_dict[get_solver_name()]    

    if the_reactions is None:
        the_reactions = model.reactions

    if hasattr(the_reactions[0], 'id'):
        the_reactions = [r.id for r in the_reactions]


    the_reactions = set(the_reactions)
    irrev_reactions = {r for r in  the_reactions if model.reactions.get_by_id(r).lower_bound >= 0}


    solution_dict = LP7(model,list(irrev_reactions),epsilon=epsilon)
    consistent = {r for r,v in solution_dict.items() if abs(v) > zero_tolerance }

    inconsistent_irrev = irrev_reactions - consistent
    the_reactions = the_reactions - consistent


    
    J = the_reactions - irrev_reactions
    flipped = False
    singleton = False
    FVA_Flag = False
    while len(J) > 0 and not FVA_Flag:
        if singleton:
            blocked = find_blocked_reactions(model,reaction_list=map(model.reactions.get_by_id,J))
            consistent = consistent.union(J - set(blocked))
            FVA_Flag = True
        else:
            Ji = J
            solution_dict = LP7(model,list(Ji),epsilon=epsilon)
            consistent = consistent.union({r for r,v in solution_dict.items()
                                        if abs(v) > zero_tolerance and not r.startswith('dummy')})

        if len(J.intersection(consistent)) > 0:
            J -= consistent
            flipped = False
        else:
            if flipped or len(Ji) == 0:
                flipped = False
                if singleton:
                    J -= Ji
                    if debug:
                        print "Inconsistent reversible reactions detected:",len(Ji)
                else:
                    singleton = True
            else:
                for r in Ji:
                    reaction = model.reactions.get_by_id(r)
                    reaction_dict = dict([(k, -2*v)
                                  for k, v in reaction._metabolites.items()])

                    reaction.add_metabolites(reaction_dict)
                    lb = reaction.lower_bound
                    reaction.lower_bound = -reaction.upper_bound
                    reaction.upper_bound = -lb

                flipped = True

    return list(consistent)


def LP7(model,the_reactions=None,epsilon=1e-3,solver=None):
    model_lp7 = model.copy()
    
    for reaction in model_lp7.reactions:
        reaction.objective_coefficient = 0

    if the_reactions is None:
        the_reactions = [r.id for r in model_lp7.reactions]

    if not hasattr(list(the_reactions)[0], 'id'):
        the_reactions = map(model_lp7.reactions.get_by_id,the_reactions)

    for reaction in the_reactions:
        dummy_reaction = Reaction("dummy_rxn_" + reaction.id)
        dummy_reaction.lower_bound = 0
        dummy_reaction.upper_bound = epsilon
        dummy_reaction.objective_coefficient = 1

        model_lp7.add_reaction(dummy_reaction)

        dummy_metabolite = Metabolite("dummy_met_" + reaction.id)
        dummy_metabolite._constraint_sense = "L"
        dummy_metabolite._bound = 0

        reaction.add_metabolites({dummy_metabolite: -1})
        dummy_reaction.add_metabolites({dummy_metabolite: 1})

    model_lp7.optimize(solver=solver)
    if model_lp7.solution is None or model_lp7.solution.x_dict is None:
        print "INFEASIBLE LP7"
        return {}

    
    return dict([(k,v) for k,v in model_lp7.solution.x_dict.items() if not k.startswith('dummy_')])

def LP9(model,K,P,weights={},epsilon=1e-4,solver=None,scaling_factor=1e3,use_milp=False):
    model = model.copy()

    scaling_factor = 1/epsilon
    for reaction in model.reactions:
        reaction.objective_coefficient = 0
        reaction.lower_bound *= scaling_factor
        reaction.upper_bound *= scaling_factor

    for reaction in P:
        reaction = model.reactions.get_by_id(reaction)
        dummy_reaction = Reaction("dummy_rxn_" + reaction.id)
        dummy_reaction.lower_bound = 0.
        if use_milp:
            dummy_reaction.upper_bound = 1.
            dummy_reaction.variable_kind = 'integer'
        else:
            dummy_reaction.upper_bound = 1000.
        
        dummy_reaction.objective_coefficient = weights[reaction.id]
        model.add_reaction(dummy_reaction)

        dummy_metabolite = Metabolite("dummy_met_ub_" + reaction.id)
        dummy_metabolite._constraint_sense = "L"
        dummy_metabolite._bound = 0.
        reaction.add_metabolites({dummy_metabolite: 1.})
        if use_milp:
            dummy_reaction.add_metabolites({dummy_metabolite: -1000})
        else:
            dummy_reaction.add_metabolites({dummy_metabolite: -1})

        dummy_metabolite = Metabolite("dummy_met_lb_" + reaction.id)
        dummy_metabolite._constraint_sense = "L"
        dummy_metabolite._bound = 0.
        reaction.add_metabolites({dummy_metabolite: -1.})
        if use_milp:
            dummy_reaction.add_metabolites({dummy_metabolite: -1000.0})
        else:
            dummy_reaction.add_metabolites({dummy_metabolite: -1.})

    for reaction in K:
        reaction = model.reactions.get_by_id(reaction)
        dummy_metabolite = Metabolite("dummy_met_lb_" + reaction.id)
        dummy_metabolite._constraint_sense = "G"
        dummy_metabolite._bound = 1.
        reaction.add_metabolites({dummy_metabolite: 1.})

    
    model.optimize(objective_sense='minimize',solver=solver)
    if model.solution is None or model.solution.x_dict is None:
        print "INFEASIBLE LP9"
        return {}

    return dict([(k,v) for k,v in model.solution.x_dict.items() if not k.startswith('dummy_')])




