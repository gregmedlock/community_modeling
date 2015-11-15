# beer_community.py
# scripts and function for creating community metabolic model describing
# S. cerivisiae and T. maritima community

import cobra
import cobra.flux_analysis.reaction
import cobra.flux_analysis.essentiality
import cobra.flux_analysis.variability
import libsbml
import os

#def combine_models(model1,model2):

def essential_exchanges(model, exchange_list):
    essential = cobra.flux_analysis.essentiality.assess_medium_component_essentiality(\
                    model, the_components = exchange_list, solver = 'cglpk')
    return essential

def get_exchange_list(bigg_model):
    exchanges = []
    for reaction in bigg_model.reactions:
        #print reaction.id.endswith('t')
        if reaction.id.endswith('t'):
            print reaction.id
            exchanges.append(reaction.id)
    return exchanges

def open_exchanges(model, exchange_list):
    for reaction in model.reactions:
        if reaction.id in exchange_list:
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

print(yeast_model.optimize())
print(maritima_model.optimize())
print("--------------------------------------------------------------")
print("Finding exchange reactions")
yeast_exchanges = get_exchange_list(yeast_model)
maritima_exchanges = get_exchange_list(maritima_model)

print maritima_exchanges
print yeast_exchanges

# print("--------------------------------------------------------------")
# print("Determining essential exchange reactions")
# yeast_essential = essential_exchanges(yeast_model, yeast_exchanges)
# maritima_essential = essential_exchanges(maritima_model, maritima_exchanges)
#
# print("yeast essential exchanges")
# print yeast_essential
# print("maritima essential exchanges")
# print maritima_essential

# Open exchange reactions necessary for yeast, then open all yeast essential
# exchanges in the maritima model
open_exchanges(yeast_model, yeast_exchanges)
open_exchanges(maritima_model, maritima_exchanges)
#print(yeast_model.optimize())
#print( [reaction.objective_coefficient for reaction in yeast_model.reactions] )
#print (cobra.flux_analysis.reaction.assess(yeast_model,yeast_model.reactions.get_by_id('BIOMASS_SC5_notrace')))

yeast_model.reactions.get_by_id('BIOMASS_SC5_notrace').objective_coefficient = 1
single_objective_solution = yeast_model.optimize()
print(single_objective_solution)
yeast_model.reactions.get_by_id('ETOHt').objective_coefficient = 1
print( [reaction.objective_coefficient for reaction in yeast_model.reactions] )
multi_objective_solution = yeast_model.optimize()
print(multi_objective_solution)
