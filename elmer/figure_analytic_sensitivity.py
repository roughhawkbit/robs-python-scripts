#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import math
import numpy
import toolbox_elmer
import toolbox_plotting

qmaxA = 1
pmax = 10
kA = 10
qmaxB = 5
pmin = 0.04
kB = 30

fig = toolbox_plotting.ThesisFigure()
sub = fig.add_subplot('', 111)

analytic = toolbox_elmer.AnalyticApproach()
analytic.set_parameters(A=1, qmaxA=qmaxA, pmax=pmax, kA=kA,
                        B=1, qmaxB=qmaxB, pmin=pmin, kB=kB)

solution = analytic.calc_equal_concn_rate()

sensitivities = numpy.linspace(0, 50, 51)
variation = [analytic.sensitivity_analysis(cv=sensitivity/100, \
        return_rate=True, return_diffs=False) for sensitivity in sensitivities]
minus_var = [i[0]/solution for i in variation]
plus_var = [i[1]/solution for i in variation]
#minus_var = [math.log(i[0]/solution, 2) for i in variation]
#plus_var = [math.log(i[1]/solution, 2) for i in variation]

sub.plot(sensitivities, minus_var, 'r')
sub.plot(sensitivities, plus_var, 'b')

#sub.text(20, 35, 'Underestimate', color='r', va='center', ha='left')
#sub.text(30, 150, 'Overestimate', color='b', va='center', ha='right')
sub.text(20, 0.7, 'Underestimate', color='r', va='center', ha='left')
sub.text(30, 2.6, 'Overestimate', color='b', va='center', ha='right')

#sub.set_xlim([0, 50])
#sub.set_ylim([0, 300])

sub.set_xlabel('Uncertainity in parameters (%)')
sub.set_ylabel('Worst-case estimate/exact estimate')
fig.process_subplots()
fig.subplots_adjust(left=0.2, right=0.97, top=0.97, bottom=0.16)
fig.save('figure_analytic_sensitivity.pdf')
