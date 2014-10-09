#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import numpy
import os
import sys
import toolbox_basic
import toolbox_idynomics
import toolbox_plotting
import toolbox_results


path = sys.argv[1]
sim_dir = toolbox_idynomics.SimulationDirectory(path)

colors = [(0, 0, 1), (0, 1, 0), (1, 0, 0), (0, 1, 1), (1, 1, 0), (1, 0, 1)]

species_color_dict = {}
for species_name in sim_dir.get_species_names():
    species_color_dict[species_name] = colors[0]
    colors.pop(0)

min_max_concns = sim_dir.get_min_max_concns()

solute_name = sim_dir.get_solute_names()[0]

figures_dir = os.path.join(path, 'figures')
toolbox_basic.make_dir(figures_dir)

for iter_info in sim_dir.get_iterate_information():
    #if not iter_info.number%12 == 0:
    #    continue
    solute_output = toolbox_results.SoluteOutput(iter_info.env_output, name=solute_name)
    fig = toolbox_plotting.PlosFigure()
    axis = fig.add_subplot('', 111)
    toolbox_idynomics.color_cells_by_species(iter_info.agent_output, species_color_dict)
    toolbox_idynomics.plot_cells_2d(axis, iter_info.agent_output)
    cs = toolbox_idynomics.solute_contour(axis, solute_output,
                    concn_range=min_max_concns[solute_name], interpolation='bicubic')
    axis.contour(x, y, array, colors='k')
    toolbox_plotting.make_colorbar(axis, cs)
    axis.set_title(r'Biofilm (%s g.L$^{-1}$)'%(solute_name))
    #fig.process_subplots()
    fig.process_lines()
    lb, rt = 0.12, 0.82
    fig.subplots_adjust(left=lb, bottom=lb, right=rt, top=rt)
    fig.save(os.path.join(figures_dir, 'biofilm(%d).png'%iter_info.number))
