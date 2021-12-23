import sys
import csv
from cobra.io import read_sbml_model, write_sbml_model
from copy import copy
from dotenv import find_dotenv
import os
from os.path import dirname
from functions import *
from cameo import models
from cameo.visualization.plotting.with_plotly import PlotlyPlotter
import numpy as np
from cameo.strain_design.deterministic.flux_variability_based import FSEOF
import matplotlib.pyplot as plt

MODEL_PATH = "model/model_yeast8_rhb.xml"
model = read_yeast_model(MODEL_PATH) # loading

fseof = FSEOF(model)
result = fseof.run(target=model.reactions.EX_rHb) # Recombinant haemoglobin as the target reaction.
df = result.data_frame
change = []
for i in range(len(df)):
    change.append(df.iloc[i,9]-df.iloc[i,0])
df['Change']=change

percent_change = []
for i in range(len(df)):
    if df.iloc[i,0] != 0:
        percent_change.append(((df.iloc[i,9]-df.iloc[i,0])/df.iloc[i,0]))
    else:
        percent_change.append(np.nan)
df['Fold_change']=percent_change
df_sorted = df.reindex(df.Fold_change.sort_values(ascending=False).index) 

knockout_targets = []
file = open("fseof_reactions.txt", "w")
file.write("Fold change;Reaction/enzyme;Short enzyme name""\n")
for i in range(len(df_sorted)):
    temp_fold_change = str(df_sorted.iloc[i:i+1, [-1]].values[0][0])[:4]
    temp_reaction = model.reactions.get_by_id(df_sorted.index.values[i]).name
    temp_short_id = model.reactions.get_by_id(df_sorted.index.values[i]).id
    file.write(temp_fold_change+";")
    file.write(temp_reaction+";")
    file.write(temp_short_id)
    file.write("\n")
    if df_sorted.iloc[i:i+1, [-1]].values[0][0] == -1:
        knockout_targets.append(model.reactions.get_by_id(df_sorted.index.values[i]))


file = open("knockout_reactions.txt", "w")
file.write("Reaction name;Reaction;Involved genes""\n")
for reaction in knockout_targets:
    #descriptive_targets = print(reaction.name, reaction, reaction.gene_reaction_rule)
    file.write(str(reaction.name)+";")
    file.write(str(reaction)+";")
    file.write(str(reaction.gene_reaction_rule))
    file.write("\n")
file.close() 