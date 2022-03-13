import matplotlib.pyplot as plt
import networkx as nx
import os
import glob

def create_plots(inputFileName, outputFolderName):
    # Create images for visualising the system
    files = glob.glob(outputFolderName+"/*")
    for f in files:
        os.remove(f)

    G = nx.Graph()
    f = open(inputFileName, 'r')
    line = f.readline()
    count = 0

    while line != '':
        G.clear()
        plt.clf()
        count = count + 1
        Polymers = eval(line)
        line = f.readline()
        family = eval(line)
        line = f.readline()
        polymer_connections = eval(line)
        line = f.readline()
        index_connections = eval(line)
        line = f.readline()
        val_map = {}
        print_me = False

        if type(Polymers[0]) == list:
            print_me = True
            for pol in range(len(Polymers)):
                if type(Polymers[pol]) == list:
                    for mon in range(len(Polymers[pol]) - 1):
                        G.add_edge(Polymers[pol][mon], Polymers[pol][mon + 1])
                        val_map[Polymers[pol][mon]] = family[pol]
                        val_map[Polymers[pol][mon + 1]] = family[pol]
                else:
                    G.add_node(Polymers[pol])
                    if type(family) == int:
                        val_map[Polymers[pol]] = family
                    else:
                        val_map[Polymers[pol]] = family[pol]

            for con in range(len(polymer_connections)):
                G.add_edge(polymer_connections[con][0][index_connections[con][0]], polymer_connections[con][1][index_connections[con][1]])
        else:
            if len(Polymers) > 1:
                print_me = True
                for mon in range(len(Polymers) - 1):
                    G.add_edge(Polymers[mon], Polymers[mon + 1])
                    val_map[Polymers[mon]] = family
                    val_map[Polymers[mon + 1]] = family
            else:
                G.add_node(Polymers[0])
                val_map[Polymers[0]] = family

        values = [val_map.get(node, 0.25) for node in G.nodes()]
        nx.draw_planar(G, cmap=plt.get_cmap('viridis'), node_color=values, node_size=10)
        plt.savefig(outputFolderName + "/image" + str(count) + ".png")

    f.close()
    plt.close()