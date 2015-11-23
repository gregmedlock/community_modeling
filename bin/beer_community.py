# beer_community.py
# scripts and function for creating community metabolic model describing
# S. cerivisiae and T. maritima community

import cobra
import cobra.flux_analysis.reaction
import cobra.flux_analysis.essentiality
import cobra.flux_analysis.variability
import libsbml
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('ggplot')

from model_stitcher import stitch_models

#def combine_models(model1,model2):

def essential_exchanges(model, exchange_list):
    essential = cobra.flux_analysis.essentiality.assess_medium_component_essentiality(\
                    model, the_components = exchange_list, solver = 'cglpk')
    return essential

def get_exchange_reaction_list(bigg_model):
    exchange_reaction_ids = []
    for reaction in bigg_model.reactions:
        if reaction.id.startswith('EX_'):
            exchange_reaction_ids.append(reaction.id)
    return exchange_reaction_ids

def get_transport_reaction_list(bigg_model,suffix="t"):
    transport_reaction_ids = []
    for reaction in bigg_model.reactions:
        if reaction.id.endswith(suffix):
            transport_reaction_ids.append(reaction.id)
    return transport_reaction_ids



def open_exchanges(model):
    # enable flux through all exchange reactions
    for reaction in model.reactions:
        if reaction.id.startswith('EX_'):
            reaction.lower_bound = -1000
            reaction.upper_bound = 1000

def simulate_variable_abundance(stitched_model,number_of_objective_ratios=20):
    objective_step = float(1)/number_of_objective_ratios
    yeast_objective = np.zeros(number_of_objective_ratios+1)
    thermotoga_objective = np.zeros(number_of_objective_ratios+1)
    yeast_objective[0] = 1
    thermotoga_objective[0] = 0
    for i in range(1,number_of_objective_ratios+1):
        yeast_objective[i] = yeast_objective[i-1] - objective_step
        thermotoga_objective[i] = thermotoga_objective[i-1] + objective_step

    yeast_biomass = np.zeros(21)
    thermotoga_biomass = np.zeros(21)
    for i in range(0,len(yeast_objective)):
        stitched_model.reactions.get_by_id('BIOMASS_SC5_notrace_mod1').\
            objective_coefficient = yeast_objective[i]
        stitched_model.reactions.get_by_id('BIOMASS_Ecoli_TM_mod2').\
            objective_coefficient = thermotoga_objective[i]
        solution = stitched_model.optimize()
        yeast_biomass[i] = solution.x_dict['BIOMASS_SC5_notrace_mod1']
        thermotoga_biomass[i] = solution.x_dict['BIOMASS_Ecoli_TM_mod2']

    biomass_dict = {'Yeast biomass':yeast_biomass,'Thermotoga biomass':thermotoga_biomass}

    for species in biomass_dict.keys():
        biomass_dict[species][biomass_dict[species] < 0] = 0.000001
    yeast_objective[yeast_objective < 0] = 0

    biomass_dataframe = pd.DataFrame(biomass_dict)
    biomass_dataframe['Yeast objective coefficient'] = pd.Series(yeast_objective)
    return biomass_dataframe


os.chdir('../models')
yeast_file = 'iMM904.xml'
maritima_file = 'iLJ478.xml'
yeast_path = '/'.join([os.getcwd(),yeast_file])
maritima_path = '/'.join([os.getcwd(),maritima_file])
yeast_model = cobra.io.read_sbml_model(yeast_path)
maritima_model = cobra.io.read_sbml_model(maritima_path)
os.chdir('..')

print("--------------------------------------------------------------")
print("Determining solutions without changing parameters")
print(yeast_model.optimize())
print(maritima_model.optimize())

print("--------------------------------------------------------------")
print("Finding transport reactions")
yeast_transport_reactions = get_transport_reaction_list(yeast_model)
maritima_transport_reactions = get_transport_reaction_list(maritima_model)

print("--------------------------------------------------------------")
print("Finding exchange reactions")
yeast_exchange_reactions = get_exchange_reaction_list(yeast_model)
maritima_exchange_reactions = get_exchange_reaction_list(maritima_model)

print("--------------------------------------------------------------")
print("Stitching models")
stitched_model = stitch_models(yeast_model,maritima_model)

print("--------------------------------------------------------------")
print("Performing FBA with a range of biomass ratios")
number_of_objective_ratios = 20
fig = plt.figure()
stitched_model.reactions.get_by_id('ETOHt_mod1').lower_bound = -10000
stitched_model.reactions.get_by_id('ETOHt_mod1').upper_bound = 10000
biomass_dataframe = simulate_variable_abundance(stitched_model,number_of_objective_ratios)
ax1 = fig.add_subplot(1,5,1)
biomass_dataframe.plot(kind='area',stacked=True,ax=ax1,sharey=True,\
        y=['Yeast biomass','Thermotoga biomass'])
ax1.xaxis.set_visible(False)

stitched_model.reactions.get_by_id('ETOHt_mod1').upper_bound = 0
biomass_dataframe = simulate_variable_abundance(stitched_model,number_of_objective_ratios)
ax2 = fig.add_subplot(1,5,2)
biomass_dataframe.plot(kind='area',stacked=True,ax=ax2,sharey=True,\
        y=['Yeast biomass','Thermotoga biomass'])
ax2.legend().set_visible(False)
ax2.xaxis.set_visible(False)

stitched_model.reactions.get_by_id('ETOHt_mod1').upper_bound = -1000
biomass_dataframe = simulate_variable_abundance(stitched_model,number_of_objective_ratios)
ax3 = fig.add_subplot(1,5,3)
biomass_dataframe.plot(kind='area',stacked=True,ax=ax3,sharey=True,\
        x='Yeast objective coefficient',y=['Yeast biomass','Thermotoga biomass'])
ax3.legend().set_visible(False)

stitched_model.reactions.get_by_id('ETOHt_mod1').upper_bound = -2000
biomass_dataframe = simulate_variable_abundance(stitched_model,number_of_objective_ratios)
ax4 = fig.add_subplot(1,5,4)
biomass_dataframe.plot(kind='area',stacked=True,ax=ax4,sharey=True,\
        y=['Yeast biomass','Thermotoga biomass'])
ax4.legend().set_visible(False)
ax4.xaxis.set_visible(False)

stitched_model.reactions.get_by_id('ETOHt_mod1').upper_bound = -3000
biomass_dataframe = simulate_variable_abundance(stitched_model,number_of_objective_ratios)
ax5 = fig.add_subplot(1,5,5)
biomass_dataframe.plot(kind='area',stacked=True,ax=ax5,sharey=True,\
        y=['Yeast biomass','Thermotoga biomass'])
ax5.legend().set_visible(False)
ax5.xaxis.set_visible(False)


plt.show()
# yeast_biomass = np.zeros(21)
# thermotoga_biomass = np.zeros(21)
# stitched_model.reactions.get_by_id('ETOHt_mod1').lower_bound = -10000
# stitched_model.reactions.get_by_id('ETOHt_mod1').upper_bound = -3000
# for i in range(0,len(yeast_objective)):
#     stitched_model.reactions.get_by_id('BIOMASS_SC5_notrace_mod1').\
#         objective_coefficient = yeast_objective[i]
#     stitched_model.reactions.get_by_id('BIOMASS_Ecoli_TM_mod2').\
#         objective_coefficient = thermotoga_objective[i]
#     solution = stitched_model.optimize()
#     yeast_biomass[i] = solution.x_dict['BIOMASS_SC5_notrace_mod1']
#     thermotoga_biomass[i] = solution.x_dict['BIOMASS_Ecoli_TM_mod2']
#
# biomass_dict = {'Yeast biomass':yeast_biomass,'Thermotoga biomass':thermotoga_biomass}
# for species in biomass_dict.keys():
#     biomass_dict[species][biomass_dict[species] < 0] = 0
#
#
# biomass_dataframe = pd.DataFrame(biomass_dict)
# biomass_dataframe['Yeast objective coefficient'] = pd.Series(yeast_objective)

# ax2 = fig.add_subplot(1,2,2)
# biomass_dataframe.plot(kind='area',stacked=True,ax=ax2,\
#         x='Yeast objective coefficient',y=['Yeast biomass','Thermotoga biomass'])
# plt.show()
# print biomass_dataframe
# print yeast_biomass
# print thermotoga_biomass
# stitched_model.reactions.get_by_id('BIOMASS_SC5_notrace_mod1').objective_coefficient = .6
# stitched_model.reactions.get_by_id('BIOMASS_Ecoli_TM_mod2').objective_coefficient = .4
# stitched_model_solution = stitched_model.optimize()
# print 'objective is ',stitched_model.objective
# print 'Flux through yeast biomass is ',stitched_model_solution.x_dict['BIOMASS_SC5_notrace_mod1']
# print 'Flux through thermotoga biomass is ',stitched_model_solution.x_dict['BIOMASS_Ecoli_TM_mod2']

# print("--------------------------------------------------------------")
# print("Opening all exchange reactions, simulating abundance ratios")
# open_exchanges(stitched_model)
# stitched_model_solution = stitched_model.optimize()
# print 'Flux through yeast biomass is ',stitched_model_solution.x_dict['BIOMASS_SC5_notrace_mod1']
# print 'Flux through thermotoga biomass is ',stitched_model_solution.x_dict['BIOMASS_Ecoli_TM_mod2']

# Open exchange reactions necessary for yeast, then open all yeast essential
# exchanges in the maritima model
#open_exchanges(yeast_model)
#open_exchanges(maritima_model)

# yeast_model.reactions.get_by_id('BIOMASS_SC5_notrace').objective_coefficient = 1
# single_objective_solution = yeast_model.optimize()
# print(single_objective_solution)
# yeast_model.reactions.get_by_id('ETOHt').objective_coefficient = 0
# ethanol_transport = yeast_model.reactions.get_by_id('ETOHt')
