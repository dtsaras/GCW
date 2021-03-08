import numpy as np
import networkx as nx 
# from scipy.cluster.vq import kmeans2, whiten
from sklearn.cluster import DBSCAN, KMeans
from array import array
import argparse
import random
from DGCD import DGCD
import sys
# import snap
import forest_fire


#Need to read
#clust_labels have the id of users starting from 0
#clusters is sets of users with ids starting from 1 

parser = argparse.ArgumentParser(description='Create isolated communities with either kmeans or dbscan')
parser.add_argument('-alg', '--algorithm', type=str, default="kmeans", help='The algorithm to perform clustering (dbscan, kmeans)')
parser.add_argument('-dt', '--dataset', type=str, default="gowalla", help='The dataset to transform')
parser.add_argument('-k', '--num_clusts', type=int, default=32, help='Number of clusters for kmeans')
parser.add_argument('-e_db', '--epsilon_db', type=float, default=2, help='Epsilon for DBSCAN')
parser.add_argument('-e_dgcd', '--epsilon_dgcd', type=float, default=0.45, help='Epsilon for DGCD')
parser.add_argument('-g', '--gamma', type=float, default=300, help='Gamma for DGCD')
parser.add_argument('-m', '--min_samples', type=int, default=3, help='Min number of neighbors for a node to be core for DBSCAN or DGCD')
parser.add_argument('-p', '--prob', type=float, default=0, help='Percent of edges to discard')
parser.add_argument('-w', '--weight', type=float, default=0.001, help='Constant weight for the bridges')
parser.add_argument('-f', '--factor', type=float, default=1, help='Factot to increase the weights by')
parser.add_argument('-fk', '--forrest_fire_k', type=int, default=1000, help='The number of output nodes of the graph')
parser.add_argument('-stats', '--graph_statistics', action='store_true', default=False, help='Print Graph statistics')


args = parser.parse_args()
if args.dataset.lower() == "gowalla" or args.dataset.lower() == "g":
	path_in = "./data/Gowalla/"
	path_out = "./data/Gowalla/" + args.algorithm.lower() + "/"
	dataset = "gowalla"
	user_locs_name = path_in + "users_locations_" + dataset + ".txt"
	user_locs_f = open(user_locs_name, "r")
	user_locs_f_out = open(path_out + "users_locations_" + dataset + ".txt", "w")
elif args.dataset.lower() == "flixster" or args.dataset.lower() == "x":
	path_in = "./data/Flixster/"
	path_out = "./data/Flixster/" + args.algorithm.lower() + "/"
	dataset = "flixster"
	similarities_in = open("./data/Flixster/flixster_similarities.txt", "r")
	similarities_out = open("./data/Flixster/" + args.algorithm.lower() + "/flixster_similarities.txt", "w")
elif args.dataset.lower() == "diggs" or args.dataset.lower() == "d":
	path_in = "./data/Diggs/"
	path_out = "./data/Diggs/" + args.algorithm.lower() + "/"
	dataset = "diggs"
	similarities_in = open("./data/Diggs/diggs_similarities.txt", "r")
	similarities_out = open("./data/Diggs/" + args.algorithm.lower() + "/diggs_similarities.txt", "w")
elif args.dataset.lower() == "lastfm" or args.dataset.lower() == "l":
	path_in = "./data/Lastfm/"
	path_out = "./data/Lastfm/" + args.algorithm.lower() + "/"
	dataset = "lastfm"
	similarities_in = open("./data/Lastfm/lastfm_similarities.txt", "r")
	similarities_out = open("./data/Lastfm/" + args.algorithm.lower() + "/lastfm_similarities.txt", "w")
elif args.dataset.lower() == "goodreads" or args.dataset.lower() == "r":
	path_in = "./data/Goodreads/"
	path_out = "./data/Goodreads/" + args.algorithm.lower() + "/"
	dataset = "goodreads"
	similarities_in = open("./data/Goodreads/goodreads_similarities.txt", "r")
	similarities_out = open("./data/Goodreads/" + args.algorithm.lower() + "/goodreads_similarities.txt", "w")
elif args.dataset.lower() == "foursquare" or args.dataset.lower() == "f":
	path_in = "./data/Foursquare/"
	path_out = "./data/Foursquare/" + args.algorithm.lower() + "/"
	dataset = "foursquare"
	user_locs_name = path_in + "users_locations_" + dataset + ".txt"
	user_locs_f = open(user_locs_name, "r")
	user_locs_f_out = open(path_out + "users_locations_" + dataset + ".txt", "w")



f_name = path_in + dataset + ".txt"
f = open(f_name, "r")

if not args.graph_statistics:
	fout = open(path_out + dataset + ".bin", "wb")

if dataset == "gowalla" or dataset == "foursquare":
	user_locs = []
	V = int(user_locs_f.readline())
	print("Original num of nodes:", V)
	for i in range(V):
		l = user_locs_f.readline().split()
		user_locs.append([float(l[1]), float(l[2])])

	user_locs = np.array(user_locs)
	user_locs_f.close()

	# Clustering
	kms_per_radian = 6371.0088
	epsilon = args.epsilon_db / kms_per_radian
	ff = False
	#Create graph for algs
	G = nx.DiGraph()
	G.add_nodes_from(list(range(1, V+1)))
	E = int(f.readline().split()[1])
	for i in range(E):
		u, v = list(map(int, f.readline().split()))
		G.add_edge(u, v)

	print("Working on Clustering\r")
	clusters = {}
	if args.algorithm == "dbscan":
		print("Algorithm:" , args.algorithm, "| Epsilon:", args.epsilon_db, "| Mu:", args.min_samples)
		clust = DBSCAN(eps=epsilon, min_samples=args.min_samples, algorithm='ball_tree', metric='haversine').fit(np.radians(user_locs))
		cluster_labels = clust.labels_
		for i in range(len(cluster_labels)):
			if cluster_labels[i] in clusters:
				clusters[cluster_labels[i]].add(i+1)
			else:
				clusters[cluster_labels[i]] = {i+1}
	elif args.algorithm == "dgcd":
		print("Algorithm:" , args.algorithm, "| Gamma:", args.gamma, "| Epsilon:", args.epsilon_dgcd, "| Mu:", args.min_samples)
		# clust = DGCD(G, gamma=300, epsilon=0.45, mu=3, user_locs=user_locs)
		cluster_labels = [0] * G.number_of_nodes() #This is clust_labels the ids start from 0

		clust = DGCD(G, gamma=args.gamma, epsilon=args.epsilon_dgcd, mu=args.min_samples, user_locs=user_locs)
		clusters[0] = set()
		c = 1
		for i in clust.keys():
			if len(clust[i]) == 1:
				clusters[0].add(clust[i][0])
			else:
				for j in clust[i]:
					cluster_labels[j - 1] = c
				clusters[c] = set(clust[i])
				c += 1
		for i in range(len(cluster_labels)):
			if cluster_labels[i] == 0:
				clusters[0].add(i + 1)
	elif args.algorithm == "ff":
		print("Donkey")
		ff = True
	elif args.algorithm == "kmeans":
		#These use as ids starting from 0 instead of 1
		print("Algorithm:" , args.algorithm, "| k:", args.num_clusts)
		clust = KMeans(n_clusters=args.num_clusts, random_state=0, n_jobs=-1).fit(np.radians(user_locs))
		cluster_labels = clust.labels_
		#This part fixes the ids to start from 1 instead
		for i in range(len(cluster_labels)):
			if cluster_labels[i] in clusters:
				clusters[cluster_labels[i]].add(i+1)
			else:
				clusters[cluster_labels[i]] = {i+1}

	if not ff:
		print("Number of clusters:", len(set(cluster_labels)))
	# print(set(cluster_labels))


	#Recreate the graph
	bridges = []
	G = nx.DiGraph()
	G.clear()
	G.add_nodes_from(list(range(1, V+1)), in_deg=0)
	# nx.set_node_attributes(G, 0, 'in_deg')
	f.seek(0)
	E = int(f.readline().split()[1])

	for i in range(E):
		u, v = list(map(int, f.readline().split()))
		# label_u = cluster_labels[u - 1]
		# if ff or v in clusters[label_u]:
		if ff or v in clusters[cluster_labels[u - 1]]:
			G.add_edge(u, v)
			G.nodes[v]["in_deg"] += 1
		elif args.prob == 0 or random.random() >= args.prob:
			G.add_edge(u, v)
			bridges.append((u, v))
	print('Number of bridges:', len(bridges))

	if ff:
		# G_snap = snap.TNGraph.New()
		# for i in range(1, V+1):
		# 	G_snap.AddNode(i)
		# for i in range(1, V+1):
		# 	for j in G.edges(i):
		# 		G_snap.AddEdge(i, h)
		# FF = snap.TForestFire(G_snap, 0.7, 0.2)
		# H_snap = FF.GenGraph(100, 0.7, 0.2)
		# for EI in H_snap.Edges():
		# 	print("edge: (%d, %d)" % (EI.GetSrcNId(), EI.GetDstNId()))
		# sys.exit()


		# G = nx.DiGraph()
		# G.add_nodes_from(list(range(1, V+1)))

		# object4=Graph_Sampling.ForestFire()
		# G = object4.forestfire(G,args.forrest_fire_k)
		for i in [500000, 0]:
			if args.forrest_fire_k != 0:
				forest_fire.forest_fire(G, 0.7, 0.2, i)
		print(G.number_of_edges())
		sys.exit()
	#Relabel the graph
	new_map = {}
	new_node = 1
	new_nodes = sorted(list(G.nodes()))
	for i in new_nodes:
		if G.in_degree(i) == 0 and G.out_degree(i) == 0:
			G.remove_node(i)
		else:
			new_map[i] = new_node
			new_node += 1
	G = nx.relabel_nodes(G, new_map)
	relabel_bridges = set()
	for e in bridges:
		relabel_bridges.add((new_map[e[0]], new_map[e[1]]))

	#Set the weights of the edges
	tot_weight = 0
	num_ed = 0
	for e in list(G.edges()):
		if args.weight != 0 and (e[0], e[1]) in relabel_bridges:
			G[e[0]][e[1]]['weight'] = args.weight
		else:
			# G[e[0]][e[1]]['weight'] = min(args.factor/G.in_degree[e[1]], 1)
			G[e[0]][e[1]]['weight'] = min(args.factor/G.in_degree(e[1]), 1)
			# tot_weight += min(args.factor/G.in_degree[e[1]], 1)
			tot_weight += min(args.factor/G.in_degree(e[1]), 1)
			num_ed += 1

	if args.graph_statistics:
		print("Number of Nodes:", G.number_of_nodes(), "Number of edges:", G.number_of_edges())
		sys.exit()


	user_locs_f_out.write(str(G.number_of_nodes()) + "\n")
	for i in range(V):
		if i + 1 in new_map:
			user_locs_f_out.write(str(new_map[i+1]) + " " + str(user_locs[i][0]) + " " + str(user_locs[i][1]) + "\n")
	user_locs_f_out.close()
else:
	#Create graph
	f.seek(0)
	V, E = list(map(int, f.readline().split()))
	G = nx.DiGraph()
	G.clear()
	G.add_nodes_from(list(range(1, V+1)), in_deg=0)
	# nx.set_node_attributes(G, 0, 'in_deg')
	edgs = {}
	for i in range(E):
		l = f.readline().split()
		u, v, w = int(l[0]), int(l[1]), float(l[2])
		if w != 0:
			G.add_edge(u + 1, v + 1)
			edgs[(u + 1, v + 1)] = w

	#Relabel the graph
	new_map = {}
	new_node = 1
	for i in range(1, V + 1):
		if G.in_degree(i) == 0 and G.out_degree(i) == 0:
			G.remove_node(i)
		else:
			new_map[i] = new_node
			new_node += 1
	G = nx.relabel_nodes(G, new_map)
	for e in edgs.keys():
		G[new_map[e[0]]][new_map[e[1]]]['weight'] = edgs[e]
		# G.add_edge(new_map[e[0]], new_map[e[1]], weight=edgs[e])
		# relabel_bridges.add((new_map[e[0]], new_map[e[1]]))

	if args.graph_statistics:
		print("Number of Nodes:", G.number_of_nodes(), "Number of edges:", G.number_of_edges())
		sys.exit()

	sims = []
	num_events = similarities_in.readline().split()[1]
	for i in range(V):
		sims.append(similarities_in.readline().split(" ", 1)[1])

	similarities_out.write(str(G.number_of_nodes()) + " " + str(num_events) + "\n")
	for i in range(V):
		if i + 1 in new_map:
			similarities_out.write(str(new_map[i+1]) + " " + sims[i])
	similarities_out.close()





#Create a reverse vesion of the graph for DSSA
H = G.reverse(copy=True)

print('New Number of Nodes:', G.number_of_nodes(), "edges:", G.number_of_edges())

degs = []
reverse_degs = []
zero_deg = 0
for i in range(1, G.number_of_nodes() + 1):
	degs.append(G.out_degree(i))
	reverse_degs.append(H.out_degree(i))
	if G.in_degree(i) == 0 and G.out_degree(i) == 0:
		zero_deg += 1

print('nodes with zero deg', zero_deg)

n = array('i', [G.number_of_nodes()])
m = array('q', [G.number_of_edges()])
n.tofile(fout)
m.tofile(fout)

degs = array('i', degs)
reverse_degs = array('i', reverse_degs)

reverse_degs.tofile(fout)
degs.tofile(fout)

weights = []
reverse_weights = []

for i in range(1, G.number_of_nodes() + 1):
	tmp = []
	rev_tmp = []
	a = array('i', list(G[i]))
	ra = array('i', list(H[i]))
	for j in a:
		tmp.append(G[i][j]['weight'])

	for j in ra:
		rev_tmp.append(H[i][j]['weight'])

	weights.append(array('f', tmp))
	reverse_weights.append(array('f', rev_tmp))

	ra.tofile(fout)
	a.tofile(fout)

for i in range(G.number_of_nodes()):
	reverse_weights[i].tofile(fout)
	weights[i].tofile(fout)





f.close()
fout.close()
# user_locs_f_out.close()







