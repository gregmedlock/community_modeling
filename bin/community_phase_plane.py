# community_phase_plane.py
import os
import cobra
from cobra.flux_analysis import calculate_phenotype_phase_plane
import model_stitcher
import matplotlib.pyplot as plt



if __name__ == '__main__':
    os.chdir('../models')
    yeast_file = 'iMM904.xml'
    maritima_file = 'iLJ478.xml'
    yeast_path = '/'.join([os.getcwd(),yeast_file])
    maritima_path = '/'.join([os.getcwd(),maritima_file])
    yeast_model = cobra.io.sbml3.read_sbml_model(yeast_path)
    maritima_model = cobra.io.sbml3.read_sbml_model(maritima_path)
    os.chdir('..')
    stitched_model = model_stitcher.stitch_models(yeast_model,maritima_model)

    data = calculate_phenotype_phase_plane(yeast_model, "EX_glc__D_e", "EX_o2_e")
    print data
    fig = plt.figure()
    ax = data.plot_matplotlib();
    plt.show()


else:
    print 'loading community_phase_plane.py'
