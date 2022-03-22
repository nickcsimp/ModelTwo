import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
import matplotlib.cm as cm  # matplotlib's color map library
from mpl_toolkits.axes_grid1 import make_axes_locatable


def create_contour_plot(inputFileName):
    contour_data = pd.read_csv(inputFileName)
    length = contour_data.pivot_table(index='Ggen', columns='Gbb', values='Length')
    x_unique = np.sort(contour_data.Ggen.unique())
    y_unique = np.sort(contour_data.Gbb.unique())

    Ggen, Gbb = np.meshgrid(x_unique, y_unique)
    # Initialize plot objects
    rcParams['figure.figsize'] = 16, 16 # sets plot size
    fig = plt.figure()

    # Generate a color mapping of the levels we've specified
    plt.contourf(Ggen, Gbb, length, cmap=cm.Reds)
    plt.colorbar()
    plt.savefig("figures/distributions.png")
