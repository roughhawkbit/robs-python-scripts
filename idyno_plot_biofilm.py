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
parser.add_option("-i", "--IterateNum", dest="iter_num", default=-1,
                    type="int", help="number of the iterate to be plotted")
parser.add_option("-a", "--PlotAll", dest="plot_all", default=False,
                    action="store_true", help="plot all iterates, ignoring -i")
parser.add_option("-s", "--SoluteName", dest="solute_name", default="none",
                    help="name of the solute to be plotted behind cells")
parser.add_option("-b", "--ColorBar", dest="color_bar", default=True,
                    action="store_false", help="include a colorbar")
parser.add_option("-m", "--Movie", dest="make_movie", default=False,
                    action="store_true", help="make a movie using the images")
''' TODO
parser.add_option("-c", "--ColorDict", dest="color_dict", default="none",
                    help="path to file containing a color dictionary")
parser.add_option("-f", "--FigureType", dest="fig_type", default="SlideFigure",
                    )
'''

(options, args) = parser.parse_args()

sim = toolbox_idynomics.SimulationDirectory(options.results_dir)

save_name = 'biofilm_'+options.solute_name

num_digits = len(str(sim.get_last_iterate_number()))

''' TODO replace with a better way of coloring cells '''
colors = [(0, 0, 1), (0, 1, 0), (1, 0, 0), (0, 1, 1), (1, 1, 0), (1, 0, 1)]
species_color_dict = {}
script = 'Species name\t\tRGB color\n'
for species_name in sim.get_species_names():
    script += species_name+'\t\t'+str(colors[0])+'\n'
    species_color_dict[species_name] = colors[0]
    colors.pop(0)
p = os.path.join(sim.figures_dir, save_name+'_info.txt')
with open(p, 'w') as f:
    f.write(script)

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
    axis.set_title(r'Biofilm (%s g L$^{-1}$)'%(options.solute_name))
    axis.fill_between([0, nJ*res], [0]*2, y2=[-res]*2, color='k', zorder=-1)
    #lb, rt = 0.05, 0.95
    figure.inset_axes()
    #figure.subplots_adjust(left=lb, bottom=lb, right=rt, top=rt)
    figure.subplots_adjust(left=0.01, bottom=0.01, right=0.9, top=0.9)
    save_num = (num_digits - len(str(iter_info.number)))*'0' + str(iter_info.number)
    figure.save(os.path.join(sim.figures_dir, '%s_%s.png'%(save_name, save_num)))


if options.plot_all:
    min_max_concns = sim.get_min_max_concns()
    for i in sim.get_iterate_numbers():
        iter_info = sim.get_single_iterate(i)
        plot(iter_info, min_max_concns)
elif options.iter_num >= 0:
    iter_info = sim.get_single_iterate(options.iter_num)
    min_max_concns = iter_info.get_min_max_concns()
    plot(iter_info, min_max_concns)

if options.make_movie:
    cmd = "ffmpeg -framerate 2  -i '"
    cmd += os.path.join(os.path.abspath(sim.figures_dir), save_name)
    cmd+= "_%"+str(num_digits)+"d.png'"
    cmd += " -pix_fmt yuv420p -r 24  "
    cmd += os.path.join(os.path.abspath(sim.movies_dir), save_name+".mp4")
    print cmd
    os.system(cmd)
