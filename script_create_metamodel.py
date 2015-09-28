import sys, re, os, glob
from cobra.core import Model
from cobra.io import read_sbml_model,write_sbml_model



folder = sys.argv[1]
metamodel_id = sys.argv[2]

assert os.path.isdir(folder)

metamodel = Model(metamodel_id)
metamodel.description = metamodel_id

reactions = set()
models = []
for fname in glob.glob(os.path.join(folder,"*.xml")):
    model = read_sbml_model(fname)
    models.append(model)
    print "%s loaded" % model.id
    for r in model.reactions:
        r.id = re.sub('_[ec][0-9]','',r.id)
        if r.id in reactions:
            continue
        metamodel.add_reaction(r.copy())
        reactions.add(r.id) 


write_sbml_model(metamodel_id,metamodel)
