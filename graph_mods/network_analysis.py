import numpy as np
import h5py
import networkx as nx
import matplotlib.pyplot as plt
import os
import scipy.io 

path = os.getcwd()

def save_data(rmax):
    '''
    Save network data corresponding to parameter 'mc' used to denote the number of mitral cells 
    in the file. By convention, h5py transposes MATLAB arrays.
    '''

    f = h5py.File(f'{path}/data{rmax}/fullNetwork.mat', 'r')
    data = np.array(f.get("network"))
    G = nx.Graph()

    for i in range(data.shape[0]):
        G.add_node(f"Granule Cell {i}", bipartite=0)
        for j in range(data.shape[1]):
            G.add_node(f"Mitral Cell {j}", bipartite=1)
            if data[i][j] == 1:
                G.add_edge(f"Granule Cell {i}", f"Mitral Cell {j}")
    nx.write_graphml(G, path + f"/data{rmax}/graph{rmax}.graphml")


def plot_bipartition(G, l, r):
    '''
    Plots mitral and granule cells as a bipartite graph, with nodes in each
    independent set on either half of the plot.
    '''

    pos = {}
    pos.update((node, (1, index)) for index, node in enumerate(l))
    pos.update((node, (2, index)) for index, node in enumerate(r))

    # Pretty large picture drawn, as graph usually has a lot of nodes. Adjust accordingly

    plt.figure(figsize=(80, 80))
    nx.draw(G, pos=pos, node_size=60,font_size=8)
    plt.savefig(os.getcwd() + "/graph_plot.png")
    
def compute_and_save_centralities(G, k):
    ''' 
    Computes the centralities of the nodes and saves them in an easily readable format
    This saves compute time, as computing betweenness and closeness centralities is 
    expensive.
    '''

    eig = sorted(nx.eigenvector_centrality_numpy(G).items(), reverse=True,key = lambda x: x[1])
    bet = sorted(nx.betweenness_centrality(G).items(), reverse=True,key = lambda x: x[1])
    clo = sorted(nx.closeness_centrality(G).items(), reverse=True,key = lambda x: x[1])

    np.save(path + f'/data{k}/closeness.npy', np.array(clo))
    np.save(path + f'/data{k}/betweenness.npy', np.array(bet))
    np.save(path + f'/data{k}/eigenvector.npy', np.array(eig))

def plot_degree_dist(G, mitral, granule):
    '''
    Plot degree (number of neighbors) distributions of mitral and granule cells.
    '''

    mitral_dist = {k : G.degree[k] for k in mitral}
    gc_dist = {k : G.degree[k] for k in granule}

    plt.hist(mitral_dist.values())
    plt.xlabel("Degree")
    plt.ylabel("Number of Mitral Cells")

    plt.figure()

    plt.hist(gc_dist.values())
    plt.xlabel("Degree")
    plt.ylabel("Number of Granule Cells")
    plt.show()


def read_graph(k):
    '''
    Read graph and output graph and both independent sets
    '''

    G = nx.read_graphml(path + f"/data{k}/graph{k}.graphml")
    granule, mitral = [[n for n in G.nodes if G.nodes[n]['bipartite'] == i] for i in range(2)]
    return G, granule, mitral


def load_and_save_centralities(k):
    '''
    Load centralities from numpy format, parse them, and convert to MATLAB readable data (.mat). 
    This can be merged with with 'compute_and_save_centralities()'.
    '''
    
    eigenvector = np.load(path + f'/data{k}/eigenvector.npy')
    closeness = np.load(path + f'/data{k}/closeness.npy')
    betweenness = np.load(path + f'/data{k}/betweenness.npy')

    centralities = {"Eigenvector" : eigenvector, "Closeness" : closeness, "Betweenness": betweenness}

    for m, c in centralities.items():
        c = list(filter(lambda x : "Granule" not in x[0], c))
        c = list(map(lambda x : int(x[0][11:].strip()), c))
        centralities[m] = c

    scipy.io.savemat(path + f"/data{k}/centralities.mat",  centralities)


