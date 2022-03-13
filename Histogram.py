import matplotlib.pyplot as plt

def create_histogram(inputFileName):
    # Create graph for visualising the average length of polymers in the system over time
    f = open(inputFileName, 'r')

    line = f.readline()
    while line != '':
        final_line = eval(line)
        line = f.readline()

    mylist = []
    for i in range(len(final_line)):
        if i != 0:
            for j in range(int(final_line[i])):
                mylist.append(i+1)

    plt.hist(mylist, len(final_line), color="green", ec="yellow")
    plt.show()


def create_average_length_graph(inputFileName):
    # Create graph for visualising the average length of polymers in the system over time
    f = open(inputFileName, 'r')
    line = f.readline()
    average_length = []

    while line != '':
        histo = eval(line)
        count = 0
        ave_count = 0
        for i in range(len(histo)):
            count = count + histo[i]
            ave_count = ave_count + i * histo[i]
        average_length.append(ave_count / count)
        line = f.readline()

    plt.plot(average_length)
    plt.show()
