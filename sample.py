import network_analysis as na

k = 50 # k here is rmax in NetworkGenMod

na.save_data(k)
G, granule, mitral = na.read_graph(k)
na.compute_and_save_centralities(G, k)
na.load_and_save_centralities(k)

