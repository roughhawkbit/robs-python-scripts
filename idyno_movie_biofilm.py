#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
from optparse import OptionParser
import toolbox_idynomics


parser = OptionParser()
parser.add_option("-f", "--FrameRate", dest="frame_rate", default=24,
                        type="int", help="number of images per second")
parser.add_option("-r", "--ResultsDir", dest="results_dir",
                      default=os.getcwd(), help="path to results directory")
parser.add_option("-s", "--SoluteName", dest="solute_name", default="none",
                        help="name of the solute to be plotted behind cells")
(options, args) = parser.parse_args()

sim = toolbox_idynomics.SimulationDirectory(options.results_dir)

num_digits = len(str(sim.get_last_iterate_number()))

save_name = 'biofilm_'+options.solute_name

cmd = "ffmpeg -framerate "+str(options.frame_rate)+"  -i '"
cmd += os.path.join(os.path.abspath(sim.figures_dir), save_name)
cmd += "_%"+str(num_digits)+"d.png'"
cmd += " -pix_fmt yuv420p -r 24  '"
cmd += os.path.join(os.path.abspath(sim.movies_dir), save_name+".mp4'")
print cmd
os.system(cmd)