# model_stitcher.py
# Combines two genome-scale metabolic network reconstructions

import cobra
import cobra.flux_analysis.reaction
import cobra.flux_analysis.essentiality
import cobra.flux_analysis.variability
import libsbml
import os


def stitch_models(model_1,model_2,model_suffixes=['_mod1','_mod2']):
    # get compartments and abbreviations for each model
    model_1_compartments = model_1.compartments
    model_2_compartments = model_1.compartments
    # Create new dict for new compartments
    model_1_new_compartments = {}
    model_1_abbreviation_mapping = {} # mapping from old to new abbreviations
    for abbreviation,full_compartment_name in model_1_compartments.iteritems():
        if abbreviation != 'e':
            model_1_new_compartments[abbreviation + model_suffixes[0]] = \
                full_compartment_name + model_suffixes[0]
            model_1_abbreviation_mapping[abbreviation] = abbreviation + model_suffixes[0]
    model_2_new_compartments = {}
    model_2_abbreviation_mapping = {} # mapping from old to new abbreviations
    for abbreviation,full_compartment_name in model_2_compartments.iteritems():
        if abbreviation != 'e':
            model_2_new_compartments[abbreviation + model_suffixes[1]] = \
                full_compartment_name + model_suffixes[1]
            model_2_abbreviation_mapping[abbreviation] = abbreviation + model_suffixes[1]

    model_1_new_compartments['e'] = 'extracellular space'
    model_1.compartments = model_1_new_compartments
    model_2_new_compartments['e'] = 'extracellular space'
    model_2.compartments = model_2_new_compartments
    # Loop through model 1 to rename compartments for every met, from old to new dict
    for metabolite in model_1.metabolites:
        if metabolite.compartment in model_1_abbreviation_mapping.keys():
            metabolite.compartment = model_1_abbreviation_mapping[metabolite.compartment]
    # Loop through model 2 to rename compartments for every met, from old to new dict
    for metabolite in model_2.metabolites:
        if metabolite.compartment in model_2_abbreviation_mapping.keys():
            metabolite.compartment = model_2_abbreviation_mapping[metabolite.compartment]
    # Add model suffix to every reaction in model 1 and 2
    model_1_exchange_reactions = []
    for reaction in model_1.reactions:
        if reaction.id.startswith('EX_') == False: # Don't rename exchange reactions
            reaction.id = reaction.id + model_suffixes[0]
        else:
            model_1_exchange_reactions.append(reaction.id)
    print model_1_exchange_reactions
    for reaction in model_2.reactions:
        if reaction.id.startswith('EX_') == False: # Don't rename exchange reactions
            reaction.id = reaction.id + model_suffixes[1]
        else:
            if reaction.id in model_1_exchange_reactions: # remove exchanges present in other model
                print reaction.id
                model_2.remove_reactions([reaction.id])
    # Add models
    stitched_model = model_1 + model_2
    return stitched_model


os.chdir('../models')
yeast_file = 'iMM904.xml'
maritima_file = 'iLJ478.xml'
yeast_path = '/'.join([os.getcwd(),yeast_file])
maritima_path = '/'.join([os.getcwd(),maritima_file])
yeast_model = cobra.io.read_sbml_model(yeast_path)
maritima_model = cobra.io.read_sbml_model(maritima_path)
os.chdir('..')

stitched_model = stitch_models(yeast_model,maritima_model)
