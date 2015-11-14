# beer_community.py
# scripts and function for creating community metabolic model describing
# S. cerivisiae and T. maritima community

import cobra
import cobra.flux_analysis.reaction
import cobra.flux_analysis.essentiality
import cobra.flux_analysis.variability
import libsbml
import os

def combine_models(model1,model2):
