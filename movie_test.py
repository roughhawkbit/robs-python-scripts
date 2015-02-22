#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
# This script is adapted from moviewriter.py example given at
# http://matplotlib.org/examples/animation/moviewriter.html

# Need to do this first so that we have the correct back end
import matplotlib
matplotlib.use("Agg")
# Now import the matplotlib packages we need
import matplotlib.pyplot
import matplotlib.animation
# Finally, everything else
import numpy


try:
    movieWriter = matplotlib.animation.writers['ffmpeg']
except KeyError:
    print('Could not find ffmpeg, attempting to use imageMagick instead.')
    try:
        movieWriter = matplotlib.animation.writers['imagemagick']
        print('Found! But imageMagick is still problematic...')
    except KeyError:
        print('Could not find imageMagick either!')
        if matplotlib.animation.writers.list() == []:
            print('Try installing ffmpeg, imageMagick, etc')
        else:
            print('Consider trying writers in this list:')
            for item in matplotlib.animation.writers.list():
                print('\t'+item)


metadata = dict(title='Movie Test', artist='Matplotlib',
        comment='Movie support!')
writer = movieWriter(fps=15, metadata=metadata)

fig = matplotlib.pyplot.figure()
l, = matplotlib.pyplot.plot([], [], 'k-o')

matplotlib.pyplot.xlim(-5, 5)
matplotlib.pyplot.ylim(-5, 5)

x0,y0 = 0, 0

with writer.saving(fig, "writer_test.mp4", 100):
    for i in range(100):
        x0 += 0.1 * numpy.random.randn()
        y0 += 0.1 * numpy.random.randn()
        l.set_data(x0, y0)
        writer.grab_frame()