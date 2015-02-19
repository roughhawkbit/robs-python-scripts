#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
from optparse import OptionParser
import toolbox_idynomics
import toolbox_plotting
import toolbox_results


parser = OptionParser()
parser.add_option("-r", "--ResultsDir", dest="results_dir", default=os.getcwd(),
                                      help="path to results directory")
parser.add_option("-i", "--IterateNum", dest="iter_num", default=0,
                    type="int", help="number of the iterate to be plotted")
parser.add_option("-a", "--PlotAll", dest="plot_all", default=False,
                    action="store_true", help="plot all iterates, ignoring -i")
parser.add_option("-s", "--SoluteName", dest="solute_name", default="none",
                    help="name of the solute to be plotted behind cells")
parser.add_option("-b", "--ColorBar", dest="color_bar", default=True,
                    action="store_false", help="include a colorbar")
''' TODO
parser.add_option("-c", "--ColorDict", dest="color_dict", default="none",
                    help="path to file containing a color dictionary")
parser.add_option("-f", "--FigureType", dest="fig_type", default="SlideFigure",
                    )
'''

(options, args) = parser.parse_args()

sim = toolbox_idynomics.SimulationDirectory(options.results_dir)

''' TODO replace with a better way of coloring cells '''
colors = [(0, 0, 1), (0, 1, 0), (1, 0, 0), (0, 1, 1), (1, 1, 0), (1, 0, 1)]
species_color_dict = {}
for species_name in sim.get_species_names():
    species_color_dict[species_name] = colors[0]
    colors.pop(0)


def plot(iter_info, min_max_concns):
    nI = iter_info.agent_output.grid_nI
    nJ = iter_info.agent_output.grid_nJ
    res = iter_info.agent_output.grid_res
    width = toolbox_plotting.mm2inch(nJ*res)
    height = toolbox_plotting.mm2inch(nI*res)
    figure = toolbox_plotting.SlideFigure(width=width, height=height)
    axis = figure.add_subplot('', 111, frameon=False)
    toolbox_idynomics.color_cells_by_species(iter_info.agent_output,
                                                            species_color_dict)
    toolbox_idynomics.plot_cells_2d(axis, iter_info.agent_output)
    if not options.solute_name == "none":
        solute_output = toolbox_results.SoluteOutput(iter_info.env_output,
                                                     name=options.solute_name)
        cs = toolbox_idynomics.solute_contour(axis, solute_output,
                            concn_range=min_max_concns[options.solute_name],
                                                    interpolation='bicubic')
        if options.color_bar:
            toolbox_plotting.make_colorbar(axis, cs)
    #axis.set_title(r'Biofilm (%s g.L$^{-1}$)'%(solute_name))
    lb, rt = 0.01, 0.99
    figure.inset_axes()
    figure.subplots_adjust(left=lb, bottom=lb, right=rt, top=rt)
    figure.save(os.path.join(sim.figures_dir,
                'biofilm_%s(%d).png'%(options.solute_name, iter_info.number)))


if options.plot_all:
    min_max_concns = sim.get_min_max_concns()
    for i in sim.get_iterate_numbers():
        iter_info = sim.get_single_iterate(i)
        plot(iter_info, min_max_concns)
else:
    iter_info = sim.get_single_iterate(options.iter_num)
    min_max_concns = iter_info.get_min_max_concns()
    plot(iter_info, min_max_concns)
