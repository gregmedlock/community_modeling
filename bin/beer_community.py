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
import matplotlib.pyplot as plt

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

def plot_biomass_scatter(biomass_dict,xlabel):


def open_exchanges(model):
    for reaction in model.reactions:
        if reaction.id.startswith('EX_'):
            reaction.lower_bound = -10
            #reaction.upper_bound = 1000


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
print maritima_transport_reactions
print yeast_transport_reactions

print("--------------------------------------------------------------")
print("Finding exchange reactions")
yeast_exchange_reactions = get_exchange_reaction_list(yeast_model)
maritima_exchange_reactions = get_exchange_reaction_list(maritima_model)
print maritima_exchange_reactions
print yeast_exchange_reactions

# Open exchange reactions necessary for yeast, then open all yeast essential
# exchanges in the maritima model
#open_exchanges(yeast_model)
#open_exchanges(maritima_model)

# yeast_model.reactions.get_by_id('BIOMASS_SC5_notrace').objective_coefficient = 1
# single_objective_solution = yeast_model.optimize()
# print(single_objective_solution)
# yeast_model.reactions.get_by_id('ETOHt').objective_coefficient = 0
# ethanol_transport = yeast_model.reactions.get_by_id('ETOHt')
