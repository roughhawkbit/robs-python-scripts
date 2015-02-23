#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
from optparse import OptionParser
import toolbox_idynomics
import toolbox_plotting
import toolbox_results


parser = OptionParser()
parser.add_option("-a", "--PlotAll", dest="plot_all", default=False,
                    action="store_true", help="plot all iterates, ignoring -i")
parser.add_option("-b", "--ColorBar", dest="color_bar", default=False,
                            action="store_true", help="include a colorbar")
parser.add_option("-f", "--FrameOn", dest="frameon", default=False,
                        action="store_true", help="turn the figure frame on")
parser.add_option("-F", "--FigureType", dest="figure_type", default=None,
                         help="type of figure to use. Default is 'Slide'")
parser.add_option("-H", "--Height", dest="height", default=0,
                    type="int", help="figure height in inches")
parser.add_option("-i", "--IterateNum", dest="iter_num", default=0,
                    type="int", help="number of the iterate to be plotted")
parser.add_option("-I", "--IMax", dest="i_max", default=0,
                        type="int", help="maximum height to plot")
parser.add_option("-r", "--ResultsDir", dest="results_dir",
                      default=os.getcwd(), help="path to results directory")
parser.add_option("-s", "--SoluteName", dest="solute_name", default="none",
                        help="name of the solute to be plotted behind cells")
parser.add_option("-t", "--TitleOn", dest="titleon", default=False,
                        action="store_true", help="turn the figure title on")
parser.add_option("-W", "--Width", dest="width", default=0,
                    type="int", help="figure width in inches")
parser.add_option("-z", "--ZeroColorBar", dest="zero_color", default=False,
                    action="store_true",
                    help="forces the lower limit of the color bar to zero")
(options, args) = parser.parse_args()


sim = toolbox_idynomics.SimulationDirectory(options.results_dir)

save_name = 'biofilm_'+options.solute_name

num_digits = len(str(sim.get_last_iterate_number()))
    
species_color_dict = toolbox_idynomics.get_default_species_colors(sim)
toolbox_idynomics.save_color_dict(species_color_dict, sim.figures_dir)

nI, nJ, nK, res = sim.find_domain_dimensions()
if options.i_max > 0:
    nI = options.i_max

if options.figure_type == None:
    if options.height > 0: height = options.height
    else: height = toolbox_plotting.mm2inch(nI * res)
    if options.width > 0: width = options.width
    else: width = toolbox_plotting.mm2inch(nJ * res)
    figure = toolbox_plotting.SlideFigure(width=width,
                                          height=height, projection='3d')
else:
    script = "figure = toolbox_plotting."+options.figure_type+"Figure("
    if nI > 2*nJ:
        script += "height='double'"
    elif nJ > 2*nI:
        script += "double_column=True, height='single'"
    else:
        script += "double_column=True, height='double'"
    script += ")"
    try:
        exec(script)
    except:
        print 'Could not make figure!'
        print script


def plot(iter_info, min_max_concns):
    axis = figure.add_subplot('', 111, frameon=options.frameon, projection='3d')
    toolbox_idynomics.color_cells_by_species(iter_info.agent_output,
                                                            species_color_dict)
    if options.three_dim:
        toolbox_idynomics.plot_cells_3d(axis, iter_info.agent_output)
    else:
        toolbox_idynomics.plot_cells_2d(axis, iter_info.agent_output)
        axis.fill_between([0, nJ*res], [0]*2, y2=[-res]*2, color='k', zorder=-1)
        figure.subplots_adjust(left=0.01, bottom=0.01, right=0.9, top=0.9)
        figure.inset_axes()
    if not options.solute_name == "none":
        solute_output = toolbox_results.SoluteOutput(iter_info.env_output,
                                                     name=options.solute_name)
        cs = toolbox_idynomics.solute_contour(axis, solute_output,
                            concn_range=min_max_concns[options.solute_name],
                                                    interpolation='bicubic')
        if options.color_bar:
            toolbox_plotting.make_colorbar(axis, cs)
    axis.set_title(r'Biofilm (%s g L$^{-1}$)'%(options.solute_name))
    save_num = (num_digits - len(str(iter_info.number)))*'0' + str(iter_info.number)
    figure.save(os.path.join(sim.figures_dir, '%s_%s.png'%(save_name, save_num)))


if options.plot_all:
    min_max_concns = sim.get_min_max_concns()
    if options.zero_color:
        min_max_concns[0] = 0
    for i in sim.get_iterate_numbers():
        if i == 0: continue
        iter_info = sim.get_single_iterate(i)
        plot(iter_info, min_max_concns)
elif options.iter_num >= 0:
    iter_info = sim.get_single_iterate(options.iter_num)
    min_max_concns = iter_info.get_min_max_concns()
    plot(iter_info, min_max_concns)

if options.make_movie:
    cmd = "ffmpeg -framerate 8  -i '"
    cmd += os.path.join(os.path.abspath(sim.figures_dir), save_name)
    cmd+= "_%"+str(num_digits)+"d.png'"
    cmd += " -pix_fmt yuv420p -r 24  '"
    cmd += os.path.join(os.path.abspath(sim.movies_dir), save_name+".mp4'")
    print cmd
    os.system(cmd)
