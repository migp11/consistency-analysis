import networkx as nx


def find_gap_metabolites(cobra_model, blocked_reactions):

    gap_metabolites = []

    if len(blocked_reactions)>0 and hasattr(blocked_reactions[0], 'id'):
        blocked_reactions = set([r.id for r in blocked_reactions])
    
    blocked_reactions_set = set(blocked_reactions)

    for m in cobra_model.metabolites:
        reactions = set([r.id for r in m.reactions])
        if reactions.issubset(blocked_reactions_set):
            gap_metabolites.append(m.id)

    return gap_metabolites

def create_metabolic_graph(cobra_model,directed=True,reactions=None,edges_labels=False):

    if directed:
        G = nx.DiGraph()
    else:
        G = nx.Graph()

    if reactions:
        if not hasattr(reactions[0], 'id'):
            reactions = [cobra_model.reactions.get_by_id(r) for r in reactions]

        metabolites = [m for r in reactions for m in r.reactants + r.products]
    else:
        reactions = [r for r in cobra_model.reactions]
        metabolites = [m for m in cobra_model.metabolites]
        
        
    for r in reactions:
        G.add_node(r.id, label=r.id, node_class="reaction", node_id=r.id)
    for m in metabolites:
        G.add_node(m.id, label=m.id, node_class="metabolite", node_id=r.id)

    for r in reactions:
        the_metabolites = r.reactants + r.products
        for m in the_metabolites:
            # evaluating if the metabolite has been defined as a reactant or product
            if r.get_coefficient(m) < 0 :
                (tail,head) = (m.id,r.id)
            else:
                (tail,head) = (r.id,m.id)

            # adding an arc between a metabolite and a reactions
            G.add_edge(tail,head)
            label = 'irrev'
            if r.lower_bound<0:
                label = 'rev'
            if edges_labels:
                G[tail][head]['label'] = label

            G[tail][head]['reversibile'] = r.lower_bound<0

    return G

def create_gap_graph(cobra_model,gap_metabolites,blocked_reactions):

    if len(gap_metabolites)>0 and hasattr(gap_metabolites[0], 'id'):
        gap_metabolites = [m.id for m in gap_metabolites]

    if len(blocked_reactions)>0 and hasattr(blocked_reactions[0], 'id'):
        blocked_reactions = [r.id for r in blocked_reactions]

    G = create_metabolic_graph(cobra_model)

    return G.subgraph(gap_metabolites + blocked_reactions)

def find_unconnected_modules(cobra_model,blocked_reactions):
    result_dict = {}

    if len(blocked_reactions)>0 and hasattr(blocked_reactions[0], 'id'):
        blocked_reactions = [r.id for r in blocked_reactions]

    gap_metabolites = find_gap_metabolites(cobra_model, blocked_reactions)

    if len(gap_metabolites)>0 and hasattr(gap_metabolites[0], 'id'):
        gap_metabolites = [m.id for m in model.metabolites]
    G = create_gap_graph(cobra_model,gap_metabolites,blocked_reactions)
    unconnected_modules = nx.connected_components(G.to_undirected())
    unconnected_modules = sorted(unconnected_modules,key=lambda x: len(x),reverse=True)
    
    for i,um in enumerate(unconnected_modules):
        um_id = "UM"+str(i+1)
        result_dict[um_id] = um

    return result_dict
    
