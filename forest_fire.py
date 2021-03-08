import networkx as nx
from queue import Queue
import random
import numpy as np

def forest_fire(G, forward_prob, backward_prob, target_nodes):
	H = nx.DiGraph()
	#Get a random node
	#Add it to BurnedNIdH
	#NBurnedTmV, NBurningTmV, NewBurnedTmV clear these total burned, currently burning, newly burned in current time step
	visited_nodes = [False]*(G.number_of_nodes()+1)
	unvisited_set = list(G.nodes())
	random.shuffle(unvisited_set)
	q = Queue()

	while (H.number_of_nodes() < target_nodes):
		print(H.number_of_nodes())
		while (visited_nodes[unvisited_set[-1]]):
			unvisited_set.pop()
		rand_node = unvisited_set.pop()
		q.put(rand_node)
		while not q.empty():
			cur_node = q.get()
			H.add_node(cur_node)
			if not visited_nodes[cur_node]:
				visited_nodes[cur_node] = True
				for out_node in G.successors(cur_node):
					if H.number_of_nodes() + q.qsize() < target_nodes and random.random() < forward_prob:
						if not visited_nodes[out_node]:
							q.put(out_node)
						H.add_edge(cur_node, out_node)
				for in_node in G.successors(cur_node):
					if H.number_of_nodes() + q.qsize() < target_nodes and random.random() < backward_prob:
						if not visited_nodes[in_node]:
							q.put(in_node)
						H.add_edge(in_node, cur_node)
	print("Forest Fire (Number of nodes, Number of edges):", H.number_of_nodes(), H.number_of_edges())
	return H


	# for (int NNodes = Graph->GetNodes()+1; NNodes <= GraphNodes; NNodes++){}
	#loop
	#	clear NewBurnedNIdV the nodes burned at current step
	#	for the currently burning nodes:
	#		cur_node = i
	#		has_live_neighs=false
	#		NDiedFire = 0;
	#		loop outgoing neighs
	#			if neihb is not in burned node vec BurnedNIdH:
	#				has_live_neighs=true
	#				if flip coin succeeds Rnd.GetUniDev() < FwdBurnProb:
	#					add neigh to BurnedNIdH and NewBurnedNIdV and NBurned++
	#		same thing for back prop
	#
	#	`	if not has_live_neihs then NDiedFire++
	#	NBurnedTmV.Add(NBurned);
    #	NBurningTmV.Add(BurningNIdV.Len() - NDiedFire);
    #	NewBurnedTmV.Add(NewBurnedNIdV.Len());
    #	if (BurningNIdV.Empty()) break;


