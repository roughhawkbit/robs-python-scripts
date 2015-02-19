#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import numpy
import os
import sys
import toolbox_basic
import toolbox_results
import toolbox_schematic


class SimulationDirectory:
    def __init__(self, path):
        self.path = toolbox_basic.check_path(path)
        self.iterate_numbers = []
        self.iterate_information = []
        self.min_max_concns = {}
        # agent_Sum
        try:
            self.agent_Sum = os.path.join(self.path, 'agent_Sum')
            if not os.path.isdir( self.agent_Sum ):
                toolbox_basic.unzip_files(self.agent_Sum + '.zip')
            self.agent_Sum = toolbox_basic.check_path(self.agent_Sum)
        except TypeError:
            print('Could not find agent_Sum info! '+self.path)
        # agent_State
        try:
            self.agent_State = os.path.join(self.path, 'agent_State')
            if not os.path.isdir( self.agent_State ):
                toolbox_basic.unzip_files(self.agent_State + '.zip')
            self.agent_State = toolbox_basic.check_path(self.agent_State)
        except TypeError:
            print('Could not find agent_State info! '+self.path)
        # env_Sum
        try:
            self.env_Sum = os.path.join(self.path, 'env_Sum')
            if not os.path.isdir( self.env_Sum ):
                toolbox_basic.unzip_files(self.env_Sum + '.zip')
            self.env_Sum = toolbox_basic.check_path(self.env_Sum)
        except TypeError:
            print('Could not find env_Sum info! '+self.path)
        # env_State
        try:
            self.env_State = os.path.join(self.path, 'env_State')
            if not os.path.isdir( self.env_State ):
                toolbox_basic.unzip_files(self.env_State + '.zip')
            self.env_State = toolbox_basic.check_path(self.env_State)
        except TypeError:
            print('Could not find env_State info! '+self.path)
        # Figures directory
        self.figures_dir = os.path.join(self.path, 'figures')
        if not os.path.isdir(self.figures_dir):
            toolbox_basic.make_dir(self.figures_dir)
        self.movies_dir = os.path.join(self.path, 'movies')
        if not os.path.isdir(self.movies_dir):
            toolbox_basic.make_dir(self.movies_dir)


    def get_iterate_numbers(self):
        """
        Returns a (sorted) list of the iterate numbers, from agent_Sum
        """
        if not self.iterate_numbers == []:
            return self.iterate_numbers
        for f in toolbox_basic.file_list(self.agent_Sum, filetype='*.xml'):
            output = toolbox_results.Output(path=f)
            self.iterate_numbers.append(output.iterate)
        self.iterate_numbers.sort()
        return self.iterate_numbers

    def get_iterate_information(self):
        """
        Tries to read in all of the iterates for this simulation. Can be
        time-consuming for large or long simulations.
        """
        #if self.iterate_information == []:
        for i in self.get_iterate_numbers():
            self.iterate_information.append(IterateInformation(self, i))
        return self.iterate_information
    
    def get_last_iterate_number(self):
        """
        
        """
        return max(self.get_iterate_numbers())
    
    def get_single_iterate(self, number):
        """
        Tries to get information for a single iteration, first by checking the
        list of iterates already read in, then by reading in the output files.
        """
        for i in self.iterate_information:
            if i.number == number:
                return i
        i = IterateInformation(self, number)
        self.iterate_information.append(i)
        return i

    def get_min_max_concns(self):
        """

        """
        if self.min_max_concns == {}:        
            for solute_name in self.get_solute_names():
                self.min_max_concns[solute_name] = [sys.float_info.max, 0.0]
            for i in self.get_iterate_information():
                iter_min_max = i.get_min_max_concns()
                for solute_name in self.min_max_concns.keys():
                    self.min_max_concns[solute_name] = \
                        [min(self.min_max_concns[solute_name][0],
                                    iter_min_max[solute_name][0]),
                         max(self.min_max_concns[solute_name][1],
                                    iter_min_max[solute_name][1])]
        return self.min_max_concns

    def get_solute_names(self):
        """

        """
        return self.get_iterate_information()[0].env_output.get_solute_names()

    def get_species_names(self):
        """

        """
        return self.get_single_iterate(0).agent_output.get_species_names()

    def find_protocol_file_xml_tree(self, filename=None):
        """

        """
        if filename is None:
            filename = toolbox_basic.find_protocol_file_path(self.path)
        self.protocol_file_xml_tree = toolbox_basic.get_xml_tree(filename)

    def find_domain_dimensions(self):
        """

        """
        try:
            pfxt = self.protocol_file_xml_tree
        except Error:
            self.find_protocol_xml_tree()

    def clean_up(self):
        """
        Deletes all unzipped output folders TODO
        """
        pass



class ProtocolFile:
    def __init__(self, path):
        pass


class IterateInformation:
    def __init__(self, simulation_directory, iterate_number):
        self.number = iterate_number
        self.min_max_concns = {}
        agent_path = os.path.join(simulation_directory.agent_State,
                                    'agent_State(%d).xml'%(iterate_number))
        agent_path = toolbox_basic.check_path(agent_path)
        self.agent_output = toolbox_results.AgentOutput(path=agent_path)
        self.time = self.agent_output.time
        env_path = os.path.join(simulation_directory.env_State,
                                        'env_State(%d).xml'%(iterate_number))
        env_path = toolbox_basic.check_path(env_path)
        self.env_output = toolbox_results.EnvOutput(path=env_path)

    def get_min_max_concns(self):
        if self.min_max_concns == {}:
            for solute_name in self.env_output.get_solute_names():
                solute_output = toolbox_results.SoluteOutput(self.env_output,
                                                            name=solute_name)
                self.min_max_concns[solute_name] = [min(solute_output.values),
                                                    max(solute_output.values)]
        return self.min_max_concns



def draw_cell_2d(axis, cell_output, total_radius=False, zorder=0):
    """

    """
    (x, y, z) = cell_output.get_location()
    rad = cell_output.get_radius(total_radius=total_radius)
    if cell_output.color == None:
        print 'Cell has no defined color!'
        col = (0, 1, 0)
    else:
        col = cell_output.color
    #col = (0, 1, 0) if cell_output.color == None else cell_output.color
    #col = cell_output.color
    circle = toolbox_schematic.Circle()
    circle.set_defaults(edgecolor='none', facecolor=col, zorder=zorder)
    circle.set_points((y, x), rad)
    circle.draw(axis)



def plot_cells_2d(axis, agent_output, zorder=0):
    """

    """
    for cell in agent_output.get_all_cells():
        draw_cell_2d(axis, cell, zorder=zorder)
    width = agent_output.grid_nJ * agent_output.grid_res
    height = agent_output.grid_nI * agent_output.grid_res
    axis.set_xlim(0, width)
    axis.set_ylim(0, height)


def color_cells_by_species(agent_output, species_color_dict):
    """

    """
    for species in agent_output.species_outputs:
        for cell in species.members:
            cell.color = species_color_dict[species.name]


# Find a list of standard colormaps (cmap) at
# http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
# It is also possible to define your own
def solute_contour(axis, solute_output, interpolation='nearest', zorder=-10,
                    cmap='gray', concn_range=[None]*2, array_multiplier=1):
    """

    """
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
