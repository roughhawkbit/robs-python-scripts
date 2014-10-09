#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
import sys
import toolbox_basic
import toolbox_idynomics
import toolbox_plotting
import toolbox_results


path = sys.argv[1]
sim_dir = toolbox_idynomics.SimulationDirectory(path)

species_color_dict = {'MyHeterotroph' : (1, 0, 0)}#,
#                      'MyAutotrophs' : (0, 1, 0),
#                      'MySwitchHeterotroph' : (0, 0, 1)}

fig = toolbox_plotting.PlosFigure()
axis = fig.add_subplot('', 111)

plot_left, plot_right, text_left = 30, 35, 38
y, y_diff = 100, -60

for species_name in species_color_dict.keys():
    time = []
    abundance = []
    for iter_info in sim_dir.get_iterate_information():
        time.append(iter_info.time)
        species = iter_info.agent_output.get_species_by_name(species_name)
        abundance.append(species.population())
    axis.plot(time, abundance, color=species_color_dict[species_name])
    #axis.plot([plot_left, plot_right], [y]*2, color=species_color_dict[species_name])
    #axis.text(text_left, y, species_name, va='center', ha='left')
    #y += y_diff


axis.set_xlabel('Time (h)')
axis.set_ylabel('Number of cells')

fig.process_subplots()
fig.subplots_adjust(left=0.2, bottom=0.15, right=0.98, top=0.98)
fig.save(os.path.join(path, 'figures', 'total_species_abundance.png'))
