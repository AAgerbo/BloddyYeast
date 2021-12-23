import sys
import csv
from cobra.io import read_sbml_model, write_sbml_model
from copy import copy
from dotenv import find_dotenv
import os
from os.path import dirname
from functions import *
MODEL_PATH = "model/yeast-GEM.xml"
model = read_yeast_model(MODEL_PATH) # loading
write_yeast_model(model, "model/model_yeast8.xml")   # saving

from cobra import Reaction, Metabolite
add_rhb(model, "model/model_yeast8_rhb.xml")
MODEL_PATH = "model/model_yeast8_rhb.xml"
model = read_yeast_model(MODEL_PATH) # loading
