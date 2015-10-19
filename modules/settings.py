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

#import warnings
#warnings.filterwarnings("ignore")

######################################################
# GLOBAL PARAMETERS                                  #
######################################################

METAMODEL_PATH = "Data/MM166.1.xml"
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
SOLVER = None
USE_MILP = False
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
    global USE_MILP
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
    USE_MILP = bool(param_dict["use_milp"])
    SOLVER = param_dict["solver"]
    SPONTANEOUS = param_dict["spontaneous"]
