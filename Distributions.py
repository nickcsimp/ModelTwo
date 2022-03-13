import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
import matplotlib.cm as cm  # matplotlib's color map library
from mpl_toolkits.axes_grid1 import make_axes_locatable


def create_contour_plot(inputFileName):
    contour_data = pd.read_csv(inputFileName)
    z = contour_data.pivot_table(index='x', columns='y', values='z')
    x_unique = np.sort(contour_data.x.unique())
    y_unique = np.sort(contour_data.y.unique())

    x, y = np.meshgrid(x_unique, y_unique)
    # Initialize plot objects
    rcParams['figure.figsize'] = 5, 5 # sets plot size
    fig = plt.figure()

    # Generate a color mapping of the levels we've specified
    plt.contourf(y, x, z, cmap=cm.Reds)
    plt.colorbar()
    plt.show()
