#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import matplotlib
import numpy
import os
import sys
import toolbox_basic
import toolbox_plotting
import toolbox_results
import toolbox_schematic


class SimulationDirectory:
    def __init__(self, path):
        self.path = toolbox_basic.check_path(path)
        self.iterate_numbers = []
        self.iterate_information = []
        self.min_max_concns = {}

    def get_iterate_numbers(self):
        if not self.iterate_numbers == []:
            return self.iterate_numbers
        agent_dir = os.path.join(self.path, 'agent_Sum')
        if not os.path.isdir(agent_dir):
            toolbox_basic.unzip_files(agent_dir+'.zip')
        for f in toolbox_basic.file_list(agent_dir, filetype='*.xml'):
            output = toolbox_results.Output(path=f)
            self.iterate_numbers.append(output.iterate)
        return self.iterate_numbers

    def get_iterate_information(self):
        if not self.iterate_information == []:
            return self.iterate_information
        agent_dir = os.path.join(self.path, 'agent_State')
        if not os.path.isdir(agent_dir):
            toolbox_basic.unzip_files(agent_dir+'.zip')
        env_dir = os.path.join(self.path, 'env_State')
        if not os.path.isdir(env_dir):
            toolbox_basic.unzip_files(env_dir+'.zip')
        for i in self.get_iterate_numbers():
            self.iterate_information.append(IterateInformation(self, i))
        self.iterate_information.sort(key=lambda x: x.number)
        return self.iterate_information

    def get_single_iterate(self, number):
        try:
            return [i for i in self.get_iterate_information() \
                                                    if i.number == number][0]
        except Error:
            toolbox_basic.error_message('No iterate %d in'%(number), self.path)

    def get_min_max_concns(self):
        if not self.min_max_concns == {}:
            return self.min_max_concns
        for solute_name in self.get_solute_names():
            self.min_max_concns[solute_name] = [sys.float_info.max, 0.0]
        for iter_info in self.get_iterate_information():
            for solute_name in self.min_max_concns.keys():
                solute_output = toolbox_results.SoluteOutput( \
                                        iter_info.env_output, name=solute_name)
                self.min_max_concns[solute_name] = \
                        [min(self.min_max_concns[solute_name][0],
                                            min(solute_output.values)),
                         max(self.min_max_concns[solute_name][1],
                                            max(solute_output.values))]
        return self.min_max_concns

    def get_solute_names(self):
        return self.get_iterate_information()[0].env_output.get_solute_names()

    def get_species_names(self):
        return self.get_iterate_information()[0].agent_output.get_species_names()

    def find_protocol_file_xml_tree(self, filename=None):
        if filename is None:
            filename = toolbox_basic.find_protocol_file_path(self.path)
        self.protocol_file_xml_tree = toolbox_basic.get_xml_tree(filename)

    def find_domain_dimensions(self):
        try:
            pfxt = self.protocol_file_xml_tree
        except Error:
            self.find_protocol_xml_tree()


class ProtocolFile:
    def __init__(self, path):
        


class IterateInformation:
    def __init__(self, simulation_directory, iterate_number):
        self.number = iterate_number
        agent_dir = os.path.join(simulation_directory.path, 'agent_State')
        if not os.path.isdir(agent_dir):
            toolbox_basic.unzip_files(agent_dir+'.zip')
        agent_path = os.path.join(agent_dir, 'agent_State(%d).xml'%(iterate_number))
        agent_path = toolbox_basic.check_path(agent_path)
        self.agent_output = toolbox_results.AgentOutput(path=agent_path)
        self.time = self.agent_output.time
        env_dir = os.path.join(simulation_directory.path, 'env_State')
        if not os.path.isdir(env_dir):
            toolbox_basic.unzip_files(env_dir+'.zip')
        env_path = os.path.join(env_dir, 'env_State(%d).xml'%(iterate_number))
        env_path = toolbox_basic.check_path(env_path)
        self.env_output = toolbox_results.EnvOutput(path=env_path)


def draw_cell_2d(axis, cell_output, total_radius=False, zorder=0):
    (x, y, z) = cell_output.get_location()
    rad = cell_output.get_radius(total_radius=total_radius)
    col = (0, 1, 0) if cell_output.color == None else cell_output.color
    #col = cell_output.color
    circle = toolbox_schematic.Circle()
    circle.set_defaults(edgecolor='none', facecolor=col, zorder=zorder)
    circle.set_points((y, x), rad)
    circle.draw(axis)



def plot_cells_2d(axis, agent_output, zorder=0):
    for cell in agent_output.get_all_cells():
        draw_cell_2d(axis, cell, zorder=zorder)
    width = agent_output.grid_nJ * agent_output.grid_res
    height = agent_output.grid_nI * agent_output.grid_res
    axis.set_xlim(0, width)
    axis.set_ylim(0, height)


def color_cells_by_species(agent_output, species_color_dict):
    for species in agent_output.species_outputs:
        for cell in species.members:
            cell.color = species_color_dict[species.name]


# Find a list of standard colormaps (cmap) at
# http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
# It is also possible to define your own
def solute_contour(axis, solute_output, interpolation='nearest', zorder=-10,
                    cmap='gray', concn_range=[None]*2, array_multiplier=1):
    width = solute_output.grid_nJ * solute_output.grid_res
    height = solute_output.grid_nI * solute_output.grid_res
    extent = [0, width, 0, height]
    array = solute_output.concentration_array()
    if not array_multiplier == 1:
        array = numpy.multiply(array, array_multiplier)
    cs = axis.imshow(array,
            interpolation=interpolation, origin='lower', cmap=cmap,
            extent=extent, zorder=zorder, vmin=concn_range[0], vmax=concn_range[1])
    return cs
