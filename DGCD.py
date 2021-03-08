from queue import Queue 
import networkx as nx
from haversine import haversine
import math

def DGCD(G, gamma, epsilon, mu, user_locs): 
	# mu is the minimum number of users to be core
	# gamma is the radius of spatial neighbors
	# epsilon is the minimum allowed similarity
	q = Queue()
	community_id = 1 
	communities = {}
	H = {}
	processed_nodes = [False] * (G.number_of_nodes() + 1)
	for v in G.nodes():
		if not processed_nodes[v]:
			N_e, H = e_Neighbor(G, v, gamma, epsilon, H, user_locs)
			if len(N_e) >= mu:
				#assign community_id to v
				processed_nodes[v] = True #added line
				if community_id not in communities:
					communities[community_id] = []
				communities[community_id].append(v)
				for u in N_e:
					#assign community_id to u #Check here if he should be added 		
					if not processed_nodes[u]:
						#mark as processed the node u
						communities[community_id].append(u) #Maybe move this line after the if statement
						processed_nodes[u] = True
						q.put(u)
			while not q.empty():
				i = q.get()
				N_e, H = e_Neighbor(G, i, gamma, epsilon, H, user_locs)
				if len(N_e) >= mu:
					for j in N_e:
						#assign community_id to j
						if not processed_nodes[j]:
							#mark as processed the node u
							communities[community_id].append(j)
							processed_nodes[j] = True
							q.put(j)
			community_id += 1
	return communities

def DGCD1(G, gamma, epsilon, mu, user_locs): 
	# mu is the minimum number of users to be core
	# gamma is the radius of spatial neighbors
	# epsilon is the minimum allowed similarity
	q = Queue()
	community_id = 1 
	communities = {}
	H = {}
	processed_nodes = [False] * (G.number_of_nodes() + 1)
	for v in G.nodes():
		if not processed_nodes[v]:
			if community_id not in communities:
				communities[community_id] = []
			N_e, H = e_Neighbor(G, v, gamma, epsilon, H, user_locs)
			if len(N_e) >= mu:
				#assign community_id to v
				processed_nodes[v] = True #added line
				communities[community_id].append(v)
				for u in N_e:
					#assign community_id to u #Check here if he should be added 
					communities[community_id].append(u) #Maybe move this line after the if statement
					if not processed_nodes[u]:
						#mark as processed the node u
						q.put(u)
			while not q.empty():
				i = q.get()
				if not processed_nodes[i]:
					N_e, H = e_Neighbor(G, i, gamma, epsilon, H, user_locs)
					processed_nodes[i] = True #added line
					if len(N_e) >= mu:
						for j in N_e:
							#assign community_id to j
							communities[community_id].append(j)
							if not processed_nodes[j]:
								#mark as processed the node u
								q.put(j)
			community_id += 1
	return communities

def e_Neighbor(G, v, gamma, epsilon, H, user_locs):
	S = []
	N_gs_v = gs_Neighnorhood(G, v, gamma, user_locs)
	for u in N_gs_v:
		if (v,u) in H:
			if H[(v,u)] == True:
				S.append(u)
		else:
			N_gs_u = gs_Neighnorhood(G, u, gamma, user_locs)
			sigma = len(N_gs_v.intersection(N_gs_u))/math.sqrt(len(N_gs_v)*len(N_gs_u))
			if sigma >= epsilon:
				S.append(u)
				H[(v, u)] = True
	return S, H


def gs_Neighnorhood(G, v, gamma, user_locs):
	N_gs = {v}
	for u in G[v]:
		loc1 = (user_locs[v-1][0], user_locs[v-1][1])
		loc2 = (user_locs[u-1][0], user_locs[u-1][1])
		if haversine(loc1, loc2) <= gamma:
			N_gs.add(u)
	return N_gs


