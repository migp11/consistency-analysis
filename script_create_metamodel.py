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
