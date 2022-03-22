if __name__ == '__main__':
    import sys
    sys.path.append('dataAnalysis')
    import Distributions
    Distributions.create_contour_plot('LengthDist.csv')