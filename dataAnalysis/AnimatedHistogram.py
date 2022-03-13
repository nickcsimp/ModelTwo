"""
   ==================
   Animated histogram
   ==================
   """

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation

HIST_BINS = np.linspace(0, 10, 10)

def animate_histogram(inputFileName):
    f = open(inputFileName, 'r')
    line = eval(f.readline())
    mylist = []
    for i in range(len(line)):
        if i != 1:
            for j in range(int(line[i])):
                mylist.append(i+1)

    n, _ = np.histogram(mylist, HIST_BINS)


    ###############################################################################
    # To animate the histogram, we need an ``animate`` function, which generates
    # a random set of numbers and updates the heights of rectangles. We utilize a
    # python closure to track an instance of `.BarContainer` whose `.Rectangle`
    # patches we shall update.


    ###############################################################################
    # Using :func:`~matplotlib.pyplot.hist` allows us to get an instance of
    # `.BarContainer`, which is a collection of `.Rectangle` instances. Calling
    # ``prepare_animation`` will define ``animate`` function working with supplied
    # `.BarContainer`, all this is used to setup `.FuncAnimation`.
    def prepare_animation(bar_container):
        def animate(frame_number):
            next_line = eval(f.readline())
            next_mylist = []
            for i in range(len(next_line)):
                if i != 1:
                    for j in range(int(next_line[i])):
                        next_mylist.append(i + 1)
            next_n, _ = np.histogram(next_mylist, HIST_BINS)
            for count, rect in zip(next_n, bar_container.patches):
                rect.set_height(count)
            return bar_container.patches

        return animate

    fig, ax = plt.subplots()
    _, _, bar_container = ax.hist(line, HIST_BINS, lw=1,
                                  ec="yellow", fc="green", alpha=0.5)
    ax.set_ylim(top=100)  # set safe limit to ensure that all data is visible.

    ani = animation.FuncAnimation(fig, prepare_animation(bar_container), 500,
                                  repeat=False, blit=True)
    ani.save('figures/myAnimation.gif', fps=30)