# model_stitcher.py
# Combines two genome-scale metabolic network reconstructions

import cobra
import libsbml
import os

def stitch_models(base_model,add_model,model_suffixes=['_mod1','_mod2']):
    model_1 = base_model.copy()
    model_2 = add_model.copy()
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
    #print model_2.reactions.get_by_id('EX_ala__L_e').id
    # Add model suffix to every reaction in model 1 and 2
    for reaction in model_1.reactions:
        if reaction.id.startswith('EX_') == False: # Don't rename exchange reactions
            reaction.id = reaction.id + model_suffixes[0]

    for reaction in model_2.reactions:
        if reaction.id.startswith('EX_') == False: # Don't rename exchange reactions
            new_reaction = reaction.copy()
            new_reaction.id = reaction.id + model_suffixes[1]
            model_1.add_reaction(new_reaction)
        else: # search for non-duplicate exchange reactions and add to model 1
            if reaction.id not in [i.id for i in model_1.reactions]:
                new_reaction = reaction.copy()
                model_1.add_reaction(new_reaction)
    model_1.repair()
    return model_1

if __name__ == '__main__':
    os.chdir('../models')
    yeast_file = 'iMM904.xml'
    maritima_file = 'iLJ478.xml'
    yeast_path = '/'.join([os.getcwd(),yeast_file])
    maritima_path = '/'.join([os.getcwd(),maritima_file])
    yeast_model = cobra.io.sbml3.read_sbml_model(yeast_path)
    maritima_model = cobra.io.sbml3.read_sbml_model(maritima_path)
    os.chdir('..')
    stitched_model = stitch_models(yeast_model,maritima_model)
else:
    print 'loading model_stitcher.py'
