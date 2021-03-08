import sys
import copy
import os
import numpy as np


import matplotlib
import matplotlib.pyplot as plt
import matplotlib.style as style
import seaborn as sns
import math

style.use('seaborn-pastel')
style.use('seaborn-paper')
matplotlib.rcParams.update({'font.size': 30})
f_size = 20
matplotlib.rc('xtick', labelsize=f_size) 
matplotlib.rc('ytick', labelsize=f_size)
matplotlib.rc('legend', fontsize=f_size)
matplotlib.rc('axes', labelsize=f_size)
matplotlib.rc('axes', titlesize=f_size)
plt.rcParams["font.family"] = "stix"
matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams['text.usetex'] = True

markers = "os^Dxh"
colors = "bgrycm"


############# Common Parameters ###############
dataset = "g" #(f,g,x)

############### Alg Parameters ################

#File name
clustalg = "original"
model = "lt" #(ic, lt)
similarity = True
rand_edges = False
sigma = 0 
mu = 0
variable = "k" #(c,k,"")
budget = "10"
competitors = "16"

#Execution
algs = "EGNA"
mc1 = 100 #For updating AP sets
mc2 = 1000 #For total similarity
out = "def"
rnd = 2
box = True

#box 
loc = "mean"
decimal_spots = 0

################ Plot Type ####################
scatter = True
bar = False

############ Cluster Parameters ###############
calg = "dgcd" #(kmeans, dbscan)
k = 512 #number of clusters for kmeans
eps_db = 2.0 #Epsilon for DBSCAN (default = 2.0)
p = 0 #percent of edges to discard (default = 1)
w = 0.0001 #constant weight for the bridges (default = 0)
eps_dgcd = 0.4 #Epsilon for DGCD (default = 2.0)
gamma = 300 #Gamma for DGCD (default = 300)
min_samples = 3 #'Min number of neighbors for a node to be core for DBSCAN or DGCD'


############ Running Parameters ###############
# run_clust = True
run_box = False
run_clust = False
run_alg = True
exp_type = 4 #(1, 2, 3, 4) 
save = True
show = False

def bound(rho):
	return (math.e - 1)/(2*math.e - 1 + math.e * rho)

def scalability_exp(dt, r, sizes):
	icfile = open("results/scalability/" + dt + "_ic_scalability.dat", "r")
	ltfile = open("results/scalability/" + dt + "_lt_scalability.dat", "r")

	data_ic_gcw = np.empty([r, len(sizes)])
	data_lt_gcw = np.empty([r, len(sizes)])
	data_ic_na = np.empty([r, len(sizes)])
	data_lt_na = np.empty([r, len(sizes)])

	for i in range(len(sizes)):
		for j in range(r):
			for k in range(3):
				# icfile.readline()
				ltfile.readline()
			# data_ic_gcw[j][i] = float(icfile.readline().split()[0])
			data_lt_gcw[j][i] = float(ltfile.readline().split()[0])

		# l1 = float(icfile.readline().split()[0])
		# l2 = float(ltfile.readline().split()[0])
			for k in range(2):
				# icfile.readline()
				ltfile.readline()
		# data_ic_na[j][i] = float(icfile.readline().split()[0])
			data_lt_na[j][i] = float(ltfile.readline().split()[0])

	print(np.mean(data_lt_gcw, axis=0))
	print(np.mean(data_lt_na, axis=0))
	sys.exit()
	plt.style.use('grayscale')
	names = [str(i) for i in sizes]
	w = 0.4
	x1 = [i - w/2 for i in range(1, len(names) + 1)]
	x2 = [i + w/2 for i in range(1, len(names) + 1)]

	#IC
	plt.bar(x1, np.mean(data_ic, axis=0), label="GCW IC", yerr=np.std(data_ic, axis=0), width=0.4, capsize=0.8, capstyle='projecting', color='lightgray')

	#LT
	plt.bar(x2, np.mean(data_lt, axis=0), label="GCW LT", yerr=np.std(data_lt, axis=0), width=0.4, capsize=0.8, capstyle='projecting', color='dimgray')

	plt.ylim(ymin=0)
	plt.xticks(list(range(1, len(names) + 1)), names)
	plt.ylabel(r'Time$(s)$')
	plt.legend()
	plt.grid(axis='y', ls=":", alpha=0.5)
	plt.xlabel(r'$\# of Nodes$')
	

	if save:
		plt.savefig("./plots/scalability/" + dt + "_" + scalability, dpi=300, bbox_inches='tight') 
	if show:
		plt.show()
	plt.clf()

def scalability_table(dt):
	# icfile = open("results/scalability/" + dt + "_ic_scalability.dat", "r")
	# ltfile = open("results/scalability/" + dt + "_lt_scalability.dat", "r")
	data = {}
	data["IC"] = {}
	data["LT"] = {}

	data["IC"]['NA'] = {}
	data["LT"]['NA'] = {}
	data["IC"]['EG'] = {}
	data["LT"]['EG'] = {}	

	c = 0
	for m in ["IC", "LT"]:
		f = open("results/scalability/" + dt + "_" + m.lower() + "_scalability.dat", "r")
		line = f.readline()
		while line != "":
			c = int(line)
			if c not in data[m]['NA']:
				data[m]['NA'][c] = []
				data[m]['EG'][c] = []
			for i in range(2):
				name = f.readline()[:-1]
				f.readline()
				time = float(f.readline().split()[0])
				data[m][name][c].append(time)
			line = f.readline()

	s = "|"
	z = " | "
	print("{0:^17}".format("Number of Nodes"), s, "{0:^8}".format("GCW_IC"), s, "{0:^8}".format("NA_IC"), s, "{0:^8}".format("GCW_LT"), s, "{0:^8}".format("NA_LT"))
	for i in data["IC"]['NA'].keys():
		if i == 0:
			out = "{0:^17}".format("2.1M")
		else:
			out = "{0:^17}".format(str(i))
		for k in ["IC", "LT"]:
			for j in ["EG", "NA"]:
				out += z + "{0:^8}".format(round(sum(data[k][j][i])/len(data[k][j][i]),5))
		print(out)

	# print("{0:^17}".format("Number of Nodes"), s, "{0:^8}".format("GCW_IC"), s, "{0:^8}".format("GCW_LT"))
	# for i in data["IC"]['NA'].keys():
	# 	if i == 0:
	# 		out = "{0:^17}".format("2.1M")
	# 	else:
	# 		out = "{0:^17}".format(str(i))
	# 	for k in ["IC", "LT"]:
	# 		j = "EG"
	# 		out += z + "{0:^8}".format(round(sum(data[k][j][i])/len(data[k][j][i]),5))
	# 	print(out)

def box_command(locb, dec_s, dt):
	return "python3 box.py -loc " + locb + " -dec " + str(dec_s) + " -dt " + dt

def alg_command(al, dt, mod, var, e_type, m1, m2, o, rounds, cl, b, sim, re, m, s, bud, comp):
	if e_type == 2:
		al = "EGNA"
		# rounds = 10
		rounds = 1
		var = ""
	if e_type == 3:
		al = "EG"
		rounds = 1
	if e_type == 5:
		al = "EGNA"
		rounds = 1
	if dt == "test":
		al = "EGNA"
		rounds = 2
		var = ""
		cl = False
		b = False
		sim = True
		e_type = 4
		bud = 1
		comp = 2

	algorithms = " -alg " + al
	dtset = " -dt " + dt 
	md = " -m 1"
	if mod != "ic":
		md = " -m 0"
	v = " -v " + var
	if var == "":
		v = ""
	default = ""
	if var == "" or (var == "c" and bud != ""):
		default = " -k " + str(bud)
	if var == "" or (var == "k" and comp != ""):
		default += " -c " + str(comp)

	exp = " -exp " + str(e_type)
	monte_carlo1 = " -mc1 " + str(m1)
	monte_carlo2 = " -mc2 " + str(m2)
	rnds = " -rnd " + str(rounds)

	if cl != "":
		clst = " -clust " + cl

	bx = " -box F"
	if b:
		bx = " -box T"
	outp = " -o def"
	if not o == "def":
		outp = ""
	similar = " -sim yes"
	if not sim: 
		similar = " -sim no"
	r_edges = ""
	if re:
		r_edges = " -rand_edges T -mu " + str(m) + " -sigma " + str(s)

	cmd = "./GCW" + algorithms + dtset + md + v + exp + monte_carlo1 + monte_carlo2 + rnds + clst + bx + outp + similar + r_edges + default

	return cmd

def clust_command(al, dt, num_clust, prc, wght, e_db, e_dgcd, gam, core_num, ff_k):
	cmd = "python3 community_creation.py"
	d = " -dt " + dt
	perc = " -p " + str(prc)
	weight = " -w " + str(wght)
	n_clust = " -k " + str(num_clust)
	epsilon_db = " -e_db " + str(e_db)
	epsilon_dgcd = " -e_dgcd " + str(e_dgcd)
	forest_k = " -fk " + str(ff_k)
	m = " -m" + str(core_num)
	g = " -g " + str(gam)
	if al == "kmeans":
		cmd += " -alg kmeans" + d + n_clust + perc + weight
	elif al == "dgcd":
		cmd += " -alg dgcd" + d + epsilon_dgcd + m + perc + weight + g
	elif al == "ff":
		cmd += " -alg ff" + d + forest_k
	else:
		cmd += d + epsilon_db + perc + weight


	return cmd

def clust_exper(al, dt, prc, clusts=[8,16, 32, 64, 128]):
	fout = open(dt + "_" + al + ".dat", "w")
	s = ""
	for i in clusts:
		s += str(i) + " "
	fout.write(s[:-1] + "\n")
	for md in ['ic', 'lt']:
		fout.write(md + "\n")
		for c in clusts:
			fout.write(str(c)+"\n")
			cmd = clust_command(al, dt, c, prc, 0.0001, 0)
			os.system(cmd)
			if dt=="g":
				m1 = 1000
				m2 = 10000
			else:
				m1 = 1000
				m2 = 10000
			cmd = alg_command("", dt, md, "", 5, m1, m2, "def", 1, True, True, True, False, 0, 0, 10, 16)
			os.system(cmd)
			fin = open("tmp_clust_exper.dat", "r")
			fout.write(fin.read())
			fin.close()
	fout.close()

def create_fname(dt, cl, mod, sim, rnd_e, m, s, var, exp):
	fname = "results/"
	
	if cl == "":
		fname += "original/"
	else:
		fname += cl + "/"
	fname += dt
	fname += "_" + mod
	if sim:
		fname += "_ysim"
	else:
		fname += "_nsim"
	if rnd_e:
		fname += "_rand_" + str(m) + "_" + str(s)
	if exp != 2 and var != "":
		fname += "_" + var
	if exp == 1:
		fname += ".dat"
	elif exp == 2:
		fname += "_rounds.dat"
	elif exp == 3:
		fname += "_utility.dat"
	elif exp == 4:
		fname += "_mult_rounds.dat"
	return fname

def plot_exp1_new(dt, var, alg):
	icfile = open("results/gain/" + alg + "/" + dt + "_ic_" + var + "_gain.dat", "r")
	ltfile = open("results/gain/" + alg + "/" + dt + "_lt_" + var + "_gain.dat", "r")

	icfile.seek(0)
	icfile.readline()
	tmp = icfile.readline()[:-1].split()[1:]
	x = [int(c) for c in tmp]

	reps = 0
	for line in icfile:
		if line[:-1] == "total similarity":
			reps += 1
	icfile.seek(0)
	icfile.readline()
	tmp = icfile.readline()[:-1].split()[1:]
	x = [int(c) for c in tmp]

	

	data_sim_ic = {}
	data_tim_ic = {}
	data_mem_ic = {}
	algs_ic = {}

	data_sim_lt = {}
	data_tim_lt = {}
	data_mem_lt = {}
	algs_lt = {}


	icfile.seek(0)
	icfile.readline()
	for r in range(reps):
		icfile.readline()
		line = icfile.readline()[:-1]
		while line != "time" and line != "":
			l = line.split()
			if l[0] not in data_sim_ic:
				data_sim_ic[l[0]] = np.empty([reps, len(x)])
			for i in range(len(l) - 1):
				data_sim_ic[l[0]][r][i] = float(l[i + 1])
			if r == 0:
				if l[0][:l[0].find(".")] in algs_ic:
					algs_ic[l[0][:l[0].find(".")]] += 1
				else:
					algs_ic[l[0][:l[0].find(".")]] = 1
			line = icfile.readline()[:-1]

		icfile.readline()
		line = icfile.readline()[:-1]
		while line != "memory" and line != "":
			l = line.split()
			if l[0] not in data_tim_ic:
				data_tim_ic[l[0]] = np.empty([reps, len(x)])
			for i in range(len(l) - 1):
				data_tim_ic[l[0]][r][i] = float(l[i + 1])
			line = icfile.readline()[:-1]

		icfile.readline()
		line = icfile.readline()[:-1]
		while line != "total similarity" and line != "":
			l = line.split()
			if l[0] not in data_mem_ic:
				data_mem_ic[l[0]] = np.empty([reps, len(x)])
			for i in range(len(l) - 1):
				data_mem_ic[l[0]][r][i] = float(l[i + 1])
			line = icfile.readline()[:-1]


			################################
	ltfile.seek(0)
	ltfile.readline()
	for r in range(reps):
		ltfile.readline()
		line = ltfile.readline()[:-1]
		while line != "time" and line != "":
			l = line.split()
			if l[0] not in data_sim_lt:
				data_sim_lt[l[0]] = np.empty([reps, len(x)])
			for i in range(len(l) - 1):
				data_sim_lt[l[0]][r][i] = float(l[i + 1])
			if r == 0:
				if l[0][:l[0].find(".")] in algs_lt:
					algs_lt[l[0][:l[0].find(".")]] += 1
				else:
					algs_lt[l[0][:l[0].find(".")]] = 1
			line = ltfile.readline()[:-1]

		ltfile.readline()
		line = ltfile.readline()[:-1]
		while line != "memory" and line != "":
			l = line.split()
			if l[0] not in data_tim_lt:
				data_tim_lt[l[0]] = np.empty([reps, len(x)])
			for i in range(len(l) - 1):
				data_tim_lt[l[0]][r][i] = float(l[i + 1])
			line = ltfile.readline()[:-1]

		ltfile.readline()
		line = ltfile.readline()[:-1]
		while line != "total similarity" and line != "":
			l = line.split()
			if l[0] not in data_mem_lt:
				data_mem_lt[l[0]] = np.empty([reps, len(x)])
			for i in range(len(l) - 1):
				data_mem_lt[l[0]][r][i] = float(l[i + 1])
			line = ltfile.readline()[:-1]

	plt.style.use('grayscale')
	names = [str(i) for i in x]
	w = 0.4
	x1 = [i - w/2 for i in range(1, len(names) + 1)]
	x2 = [i + w/2 for i in range(1, len(names) + 1)]

	
	############# Gain ################

	#IC
	tmp = 100*np.subtract(data_sim_ic["GCW.1"], data_sim_ic["NA.1"])/data_sim_ic["NA.1"]
	plt.bar(x1, np.mean(tmp, axis=0), label="GCW IC", yerr=np.std(tmp, axis=0), width=0.4, capsize=0.8, capstyle='projecting', color='lightgray')

	#LT
	tmp2 = 100*np.subtract(data_sim_lt["GCW.1"], data_sim_lt["NA.1"])/data_sim_lt["NA.1"]
	plt.bar(x2, np.mean(tmp2, axis=0), label="GCW LT", yerr=np.std(tmp2, axis=0), width=0.4, capsize=0.8, capstyle='projecting', color='dimgray')


	plt.ylim(ymin=0)
	plt.xticks(list(range(1, len(names) + 1)), names)
	plt.ylabel(r'Gain$\%$')
	plt.legend()
	plt.grid(axis='y', ls=":", alpha=0.5)
	if var == "k":
		plt.xlabel(r'$k$')
	else:
		plt.xlabel(r'$|C|$')
	

	if save:
		plt.savefig("./plots/gain/" + alg + "/new_gain_" + dt + "_" + var, dpi=300, bbox_inches='tight') 
	if show:
		plt.show()
	plt.clf()


	# w = 0.2
	# x1 = [i - (3/2)*w for i in range(1, len(names) + 1)]
	# x2 = [i - (1/2)*w for i in range(1, len(names) + 1)]
	# x3 = [i + (1/2)*w for i in range(1, len(names) + 1)]
	# x4 = [i + (3/2)*w for i in range(1, len(names) + 1)]
	# ############# Time ################
	# #NA IC
	# plt.bar(x1, np.mean(data_tim_ic['NA.1'], axis=0), label="NA IC", yerr=np.std(data_tim_ic['NA.1'], axis=0), width=w, capsize=0.8, capstyle='projecting', color='whitesmoke', linewidth=0.2, edgecolor='black')

	# #IC
	# plt.bar(x2, np.mean(data_tim_ic['GCW.1'], axis=0), label="GCW IC", yerr=np.std(data_tim_ic['GCW.1'], axis=0), width=w, capsize=0.8, capstyle='projecting', color='lightgray')

	# #NA LT
	# plt.bar(x3, np.mean(data_tim_lt['NA.1'], axis=0), label="NA LT", yerr=np.std(data_tim_lt['NA.1'], axis=0), width=w, capsize=0.8, capstyle='projecting', color='grey')

	# #LT
	# plt.bar(x4, np.mean(data_tim_lt['GCW.1'], axis=0), label="GCW LT", yerr=np.std(data_tim_lt['GCW.1'], axis=0), width=w, capsize=0.8, capstyle='projecting', color='dimgray')

	#print times normalized
	#IC
	tmp = data_tim_ic["GCW.1"]/data_tim_ic["NA.1"]

	plt.bar(x1, np.mean(tmp, axis=0), label="GCW IC", yerr=np.std(tmp, axis=0), width=0.4, capsize=0.8, capstyle='projecting', color='lightgray')
	tmp_abs = np.mean(data_tim_ic["GCW.1"])
	s = f_ics
	y = f_ic
	# s = list(np.std(tmp, axis=0))
	# y = list(np.mean(tmp, axis=0))
	v = list(np.mean(data_tim_ic["GCW.1"], axis=0))
	for i in range(len(names)):
		plt.text(x=(x1[i]-w/2 + 0.1*(4 - len(str(round(v[i], 1))))) , y=(y[i]+(s[i]*1.1)) , s=str(round(v[i], 1)) , fontdict=dict(fontsize=16))

	#LT
	tmp2 = data_tim_lt["GCW.1"]/data_tim_lt["NA.1"]
	plt.bar(x2, np.mean(tmp2, axis=0), label="GCW LT", yerr=np.std(tmp2, axis=0), width=0.4, capsize=0.8, capstyle='projecting', color='dimgray')
	# sys.exit()
	s2 = f_lts
	y2 = f_lt
	# s2 = list(np.std(tmp2, axis=0))
	# y2 = list(np.mean(tmp2, axis=0))
	v2 = list(np.mean(data_tim_lt["GCW.1"], axis=0))
	for i in range(len(names)):
		plt.text(x=(x2[i]-w/2 + 0.05*(4 - len(str(round(v2[i], 1))))) , y=(y2[i]+(s2[i]*1.1)) , s=str(round(v2[i], 1)) , fontdict=dict(fontsize=16))

	plt.ylim(ymin=0)
	plt.xticks(list(range(1, len(names) + 1)), names)
	plt.ylabel(r'Normalized Time')
	plt.legend()
	plt.grid(axis='y', ls=":", alpha=0.5)
	if var == "k":
		plt.xlabel(r'$k$')
	else:
		plt.xlabel(r'$|C|$')
	

	if save:
		plt.savefig("./plots/gain/" + alg + "/new_time_" + dt + "_" + var, dpi=300, bbox_inches='tight') 
	if show:
		plt.show()
	plt.clf()

def plot_exp2_new(dt, alg):
	icfile = open("results/rounds/" + alg + "/" + dt + "_ic_all_rounds.dat", "r")
	ltfile = open("results/rounds/" + alg + "/" + dt + "_lt_all_rounds.dat", "r")

	reps_ic = 0
	reps1_ic = 0
	rounds_ic = 0
	for line in icfile:
		if line[0].isdigit():
			rounds_ic += 1
		if line[:-1] == "EG" or line[:-1] == "NA":
			reps1_ic += 1
		if line[:-1] == "EG":
			reps_ic += 1
	print(rounds_ic)
	rounds_ic = rounds_ic//reps1_ic
	print(rounds_ic, rounds_ic, reps1_ic)

	rounds_ic = 10
	reps_lt = 0
	reps1_lt = 0
	rounds_lt = 0
	for line in ltfile:
		if line[0].isdigit():
			rounds_lt += 1
		if line[:-1] == "EG" or line[:-1] == "NA":
			reps1_lt += 1
		if line[:-1] == "EG":
			reps_lt += 1

	rounds_lt = rounds_lt//reps1_lt
	rounds_lt = 10

	eg_gains_lt = np.empty([reps_lt, rounds_lt-1])
	eg_time_lt = np.empty([reps_lt, rounds_lt-1])
	na_gains_lt = np.empty([reps_lt, rounds_lt-1])
	na_time_lt = np.empty([reps_lt, rounds_lt-1])

	eg_gains_ic = np.empty([reps_ic, rounds_ic-1])
	eg_time_ic = np.empty([reps_ic, rounds_ic-1])
	na_gains_ic = np.empty([reps_ic, rounds_ic-1])
	na_time_ic = np.empty([reps_ic, rounds_ic-1])


	icfile.seek(0)
	for i in range(reps_ic):
		for k in range(3):
			icfile.readline()
		for j in range(rounds_ic - 1):
			l = icfile.readline().split()
			eg_time_ic[i][j] = float(l[0])
			eg_gains_ic[i][j] = float(l[1])
		icfile.readline()
		icfile.readline()
		l = icfile.readline().split()
		na_gains_ic[i][0] = float(l[1])
		na_time_ic[i][0] = float(l[0])
		for j in range(1, rounds_ic - 1):
			na_gains_ic[i][j] = na_gains_ic[i][0]
			na_time_ic[i][j] = na_time_ic[i][0]
			icfile.readline()

	ltfile.seek(0)
	for i in range(reps_lt):
		for k in range(3):
			ltfile.readline()
		for j in range(rounds_lt - 1):
			l = ltfile.readline().split()
			eg_time_lt[i][j] = float(l[0])
			eg_gains_lt[i][j] = float(l[1])
		ltfile.readline()
		ltfile.readline()
		l = ltfile.readline().split()
		na_gains_lt[i][0] = float(l[1])
		na_time_lt[i][0] = float(l[0])
		for j in range(1, rounds_lt - 1):
			na_gains_lt[i][j] = na_gains_lt[i][0]
			na_time_lt[i][j] = na_time_lt[i][0]
			ltfile.readline()

	final_gain_ic = 100*np.subtract(eg_gains_ic, na_gains_ic)/na_gains_ic
	final_gain_lt = 100*np.subtract(eg_gains_lt, na_gains_lt)/na_gains_lt

	names = []
	for i in range(rounds_lt - 1):
		names.append(str(i + 1))


	w = 0.4
	x1 = [i - w/2 for i in range(1, len(names) + 1)]
	x2 = [i + w/2 for i in range(1, len(names) + 1)]
	plt.style.use('grayscale')
	plt.ylabel(r'Gain$\%$')
	plt.xlabel('Rounds')
	plt.grid(axis='y', ls=":", alpha=0.5)
	# plt.bar(names, np.mean(final, axis=0), capsize=0.8, yerr=np.std(final, axis=0), hatch='xxx', fill=False, edgecolor="steelblue", ls='-', lw='0.6', width=0.4)
	plt.bar(x1, np.mean(final_gain_ic, axis=0), yerr=np.std(final_gain_ic, axis=0), width=w, capsize=0.8, capstyle='projecting', label="GCW IC", color='lightgray')
	plt.bar(x2, np.mean(final_gain_lt, axis=0), yerr=np.std(final_gain_lt, axis=0), width=w, capsize=0.8, capstyle='projecting', label="GCW LT", color='dimgray')

	plt.legend()
	plt.xticks(list(range(1, len(names) + 1)), names)
	if save:
		plt.savefig("./plots/rounds/" + alg + "/new_all_rounds_gain" + "_" + dt, dpi=300, bbox_inches='tight') 
	if show:
		plt.show()
	plt.clf()


	w = 0.4
	x1 = [i - w/2 for i in range(1, len(names) + 1)]
	x2 = [i + w/2 for i in range(1, len(names) + 1)]
	############# Time ################
	#IC
	plt.bar(x1, np.mean(eg_time_ic, axis=0), label="GCW IC", yerr=np.std(eg_time_ic, axis=0), width=w, capsize=0.8, capstyle='projecting', color='lightgray')

	#LT
	plt.bar(x2, np.mean(eg_time_lt, axis=0), label="GCW LT", yerr=np.std(eg_time_lt, axis=0), width=w, capsize=0.8, capstyle='projecting', color='dimgray')

	plt.ylabel(r'Time$(s)$')
	plt.xlabel('Rounds')
	plt.grid(axis='y', ls=":", alpha=0.5)
	plt.legend()
	plt.xticks(list(range(1, len(names) + 1)), names)
	if save:
		plt.savefig("./plots/rounds/" + alg + "/new_all_rounds_time" + "_" + dt, dpi=300, bbox_inches='tight') 
	if show:
		plt.show()
	plt.clf()

	icfile.close()
	ltfile.close()

def plot_exp3_new(dt, var, alg, typ):
	icfile = open("results/utility/" + alg + "/" + dt + "_ic_ysim_" + var + "_utility.dat", "r")
	ltfile = open("results/utility/" + alg + "/" + dt + "_lt_ysim_" + var + "_utility.dat", "r")

	icfile.seek(0)

	xaxis_ic = []
	beta_ic = []
	util_ic = []
	DSS_ic = []
	tot_sim_ic = []

	xaxis_lt = []
	beta_lt = []
	util_lt = []
	DSS_lt = []
	tot_sim_lt = []

	icfile.readline()
	for line in icfile:
		l = line.split()
		print(l)
		xaxis_ic.append(int(l[0]))
		util_ic.append(float(l[1]))
		beta_ic.append(float(l[2]))
		DSS_ic.append(float(l[3]))
		tot_sim_ic.append(float(l[4]))

	ltfile.readline()
	for line in ltfile:
		l = line.split()
		xaxis_lt.append(int(l[0]))
		util_lt.append(float(l[1]))
		beta_lt.append(float(l[2]))
		DSS_lt.append(float(l[3]))
		tot_sim_lt.append(float(l[4]))

	names = [str(i) for i in xaxis_ic]

	plt.style.use('grayscale')
	w = 0.4
	x1 = [i - w/2 for i in range(1, len(names) + 1)]
	x2 = [i + w/2 for i in range(1, len(names) + 1)]

	
	if typ == "r" or typ == "rho":
		
		rho_ic = []
		rho_lt = []
		bound_ic = []
		bound_lt = []

		for i in range(len(xaxis_ic)):
			rho_ic.append(abs(DSS_ic[i])/(tot_sim_ic[i]))
			rho_lt.append(abs(DSS_lt[i])/(tot_sim_lt[i]))
			bound_ic.append(max((bound(abs(DSS_ic[i])/(tot_sim_ic[i]))), 0.2))
			bound_lt.append(max((bound(abs(DSS_lt[i])/(tot_sim_lt[i]))), 0.2))

	
		f = plt.figure()
		# ax1 = f.add_subplot(111)
		# #IC
		# ax1.bar(x1, rho_ic, label="GCW IC", width=0.4, capsize=0.8, capstyle='projecting', color='lightgray')

		# # #LT
		# ax1.bar(x2, rho_lt, label="GCW LT", width=0.4, capsize=0.8, capstyle='projecting', color='dimgray')
		# plt.ylabel(r'$\rho(\mathcal{S})$')
		# plt.legend()
		# plt.ylim(ymin=0)
		#Bounds
		# ax2 = f.add_subplot(111, sharex=ax1, frameon=False)
		ax2 = f.add_subplot(111)
		ax2.plot(list(range(1, len(names) + 1)), bound_ic, ls="-", marker='*', color='black', label='GCW IC Bounds')
		ax2.plot(list(range(1, len(names) + 1)), bound_lt, ls="-", marker='o', color='k', label='GCW LT Bounds')
		ax2.yaxis.tick_right()
		ax2.set_ylim([0,bound(0)])
		ax2.yaxis.set_label_position("right")
		plt.legend()
		plt.ylabel(r'$b$')

	else:
		w = 0.2
		x1 = [i - (3/2)*w for i in range(1, len(names) + 1)]
		x2 = [i - (1/2)*w for i in range(1, len(names) + 1)]
		x3 = [i + (1/2)*w for i in range(1, len(names) + 1)]
		x4 = [i + (3/2)*w for i in range(1, len(names) + 1)]
		############# Time ################
		#NA IC
		plt.bar(x1, util_ic, label=r'avg[$\sigma_j (\mathcal{S}^j)$] IC', width=w, capsize=0.8, capstyle='projecting', color='whitesmoke', linewidth=0.2, edgecolor='black')

		#IC
		plt.bar(x2, beta_ic, label=r'avg[$\beta_j (\mathcal{S}^j)$] IC', width=w, capsize=0.8, capstyle='projecting', color='lightgray')

		#NA LT
		plt.bar(x3, util_lt, label=r'avg[$\sigma_j (\mathcal{S}^j)$] LT', width=w, capsize=0.8, capstyle='projecting', color='grey')

		#LT
		plt.bar(x4, beta_lt, label=r'avg[$\beta_j (\mathcal{S}^j)$] LT', width=w, capsize=0.8, capstyle='projecting', color='dimgray')
		plt.legend()
		plt.ylabel(r'$\sigma(\mathcal{S})$')


	
	plt.xticks(list(range(1, len(names) + 1)), names)
	plt.grid(axis='y', ls=":", alpha=0.5)
	if var == "k":
		plt.xlabel(r'$k$')
	else:
		plt.xlabel(r'$|C|$')
	

	if save:
		plt.savefig("./plots/utility/" + alg + "/utility_" + dt + "_" + var, dpi=300, bbox_inches='tight') 
	if show:
		plt.show()
	plt.clf()

def plot_exp5_new(dt, alg):
	infile = open(dt+"_"+alg+".dat")
	clusts = infile.readline().split()
	data = {}
	for d1 in ["ic", "lt"]:
		data[d1] = {}
		for d2 in ["GCW", "NA"]:
			data[d1][d2] = {}
			for d3 in ["sim", "time"]:
				data[d1][d2][d3] = []
	c_model = ""
	for i in range(2):
		c_model = infile.readline()[:-1]
		for j in range(len(clusts)):
			infile.readline()
			l1 = infile.readline().split()
			l2 = infile.readline().split()
			data[c_model][l1[0]]["sim"].append(float(l1[1]))
			data[c_model][l1[0]]["time"].append(float(l1[2]))
			data[c_model][l2[0]]["sim"].append(float(l2[1]))
			data[c_model][l2[0]]["time"].append(float(l2[2]))
	temp_out = []
	for i in range(len(clusts)):
		temp_out.append((data["lt"]["GCW"]["sim"][i] - data["lt"]["NA"]["sim"][i])/data["lt"]["NA"]["sim"][i])
	print(temp_out)






	# if scatter:
	# 	xaxis = [str(i) for i in xaxis]
	# 	col = ["darkorange", 'dimgrey', "steelblue"]
	# 	plt.plot(xaxis, util, ls="-", marker=markers[0], color=col[0], label=r'avg[$\sigma_j (\mathcal{S}^j)$]')
	# 	plt.plot(xaxis, beta, ls="-", marker=markers[1], color=col[2], label=r'avg[$\beta_j (\mathcal{S}^j)$]')
	# 	plt.ylabel(r'$\sigma (S)$')
	# 	plt.xticks(xaxis)


	################################# Percent #######################################
	# plt.ylabel(r'Violation of basicness$\%$')
	# diff_ic = []
	# diff_lt = []
	# for i in range(len(xaxis_ic)):
	# 	diff_ic.append(100*(util_ic[i] - beta_ic[i])/util_ic[i])
	# 	diff_lt.append(100*(util_lt[i] - beta_lt[i])/util_lt[i])

	# w = 0.4
	# x1 = [i - w/2 for i in range(1, len(names) + 1)]
	# x2 = [i + w/2 for i in range(1, len(names) + 1)]
	# ############# Time ################
	# #IC
	# plt.bar(x1, diff_ic, label="GCW IC", width=w, capsize=0.8, capstyle='projecting', color='lightgray')

	# #LT
	# plt.bar(x2, diff_lt, label="GCW LT", width=w, capsize=0.8, capstyle='projecting', color='dimgray')
	# 	# plt.bar(names, diff)

	################################ Flat ###########################################
	plt.ylabel(r'$\sigma (S)$')

	w = 0.2
	x1 = [i - (3/2)*w for i in range(1, len(names) + 1)]
	x2 = [i - (1/2)*w for i in range(1, len(names) + 1)]
	x3 = [i + (1/2)*w for i in range(1, len(names) + 1)]
	x4 = [i + (3/2)*w for i in range(1, len(names) + 1)]
	############# Time ################
	#NA IC
	plt.bar(x1, util_ic, label=r'avg[$\sigma_j (\mathcal{S}^j)$]' + " IC", width=w, capsize=0.8, capstyle='projecting', color='whitesmoke', linewidth=0.2, edgecolor='black')

	#IC
	plt.bar(x2, beta_ic, label=r'avg[$\beta_j (\mathcal{S}^j)$]' + " IC", width=w, capsize=0.8, capstyle='projecting', color='lightgray')

	#NA LT
	plt.bar(x3, util_lt, label=r'avg[$\sigma_j (\mathcal{S}^j)$]' + " LT", width=w, capsize=0.8, capstyle='projecting', color='grey')

	#LT
	plt.bar(x4, beta_lt, label=r'avg[$\beta_j (\mathcal{S}^j)$]' + " LT", width=w, capsize=0.8, capstyle='projecting', color='dimgray')

	# plt.plot(xaxis, util, ls="-", marker=markers[0], color=col[0], label=r'avg[$\sigma_j (\mathcal{S}^j)$]')
	# plt.plot(xaxis, beta, ls="-", marker=markers[1], color=col[2], label=r'avg[$\beta_j (\mathcal{S}^j)$]')
	









	plt.xticks(list(range(1, len(names) + 1)), names)
	# plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", fancybox=True)
	plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=2, prop={'size': 19})
	plt.grid(axis='y', ls=":", alpha=0.5)
	if var == "k":
		plt.xlabel(r'$k$')
		# plt.legend()
	else:
		plt.xlabel(r'$|C|$')
		# plt.legend(loc=4, prop={'size': 6})

	if save:
		plt.savefig("new_utility_" + dt + "_" + var, dpi=300, bbox_inches='tight') 
	if show:
		plt.show()
	plt.clf()
	icfile.close()
	ltfile.close()

def plot_exp1(fname):
	# var = "c"
	# if "_k_" in fname:
	# 	var = "k"
	# infile = open(fname, "r")
	# ts_data = {}
	# t_data = {}
	# m_data = {}
	# line = infile.readline()[:-1]
	# tmp = infile.readline()[:-1].split()[1:]
	# line = infile.readline()[:-1]
	# x = [int(c) for c in tmp]
	# max_ts = 0
	# while line != "time":
	# 	l = line.split()
	# 	ts_data[l[0]] = [float(n) for n in l[1:]]
	# 	if max(ts_data[l[0]]) > max_ts:
	# 		max_ts = max(ts_data[l[0]])
	# 	line = infile.readline()[:-1]
	# line = infile.readline()[:-1]
	# line = infile.readline()[:-1]
	# while line != "memory":
	# 	l = line.split()
	# 	t_data[l[0]] = [float(n) for n in l[1:]]
	# 	line = infile.readline()[:-1]
	# infile.readline()
	# for line in infile:
	# 	l = line.split()
	# 	m_data[l[0]] = [float(n) for n in l[1:]]

	# if scatter:
	# 	i = 0
	# 	for j in ts_data.keys():
	# 		print(j, ts_data[j])
	# 		plt.plot(x, ts_data[j], markers[i] + "-" + colors[i], label=j)
	# 		i += 1
	# 	plt.legend()
	# 	plt.grid(True)
	# 	plt.axis([min(x), max(x), 0, max_ts * 1.1])
	# 	plt.ylabel(r'$\sigma (S)$')
	# 	plt.xticks(x)

	# if len(ts_data.keys()) == 2:
	# 	diff = []
	# 	k1 = list(ts_data.keys())[0]
	# 	k2 = list(ts_data.keys())[1]
	# 	for i in range(len(ts_data[k1])):
	# 		diff.append(round(abs(ts_data[k1][i] - ts_data[k2][i])/min(ts_data[k1][i], ts_data[k2][i]), 3)*100)
	# 	print(diff)
	# 	print("average gain:", sum(diff)/len(diff))
	var = "c"
	if "_k_" in fname:
		var = "k"
	infile = open(fname, "r")
	reps = 0
	for line in infile:
		if line[:-1] == "total similarity":
			reps += 1
	infile.seek(0)
	infile.readline()
	tmp = infile.readline()[:-1].split()[1:]
	x = [int(c) for c in tmp]

	infile.seek(0)
	data_sim = {}
	data_tim = {}
	data_mem = {}
	algs = {}

	infile.readline()
	for r in range(reps):
		infile.readline()
		line = infile.readline()[:-1]
		while line != "time" and line != "":
			l = line.split()
			if l[0] not in data_sim:
				data_sim[l[0]] = np.empty([reps, len(x)])
			for i in range(len(l) - 1):
				data_sim[l[0]][r][i] = float(l[i + 1])
			if r == 0:
				if l[0][:l[0].find(".")] in algs:
					algs[l[0][:l[0].find(".")]] += 1
				else:
					algs[l[0][:l[0].find(".")]] = 1
			line = infile.readline()[:-1]

		infile.readline()
		line = infile.readline()[:-1]
		while line != "memory" and line != "":
			l = line.split()
			if l[0] not in data_tim:
				data_tim[l[0]] = np.empty([reps, len(x)])
			for i in range(len(l) - 1):
				data_tim[l[0]][r][i] = float(l[i + 1])
			line = infile.readline()[:-1]

		infile.readline()
		line = infile.readline()[:-1]
		while line != "total similarity" and line != "":
			l = line.split()
			if l[0] not in data_mem:
				data_mem[l[0]] = np.empty([reps, len(x)])
			for i in range(len(l) - 1):
				data_mem[l[0]][r][i] = float(l[i + 1])
			line = infile.readline()[:-1]

	gs = []
	col = ["steelblue", "darkorange", 'dimgrey']
	col1 = ['dimgrey', 'black']
	names = [str(i) for i in x]
	error_kw=dict(lw=5, capsize=5, capthick=3)
	f = 0
	for k in algs.keys():
		if k != "NA":
			for i in range(algs[k], 0, -1):
				new_key = k + "." + str(i)
				tmp = 100*np.subtract(data_sim[new_key], data_sim["NA.1"])/data_sim["NA.1"]
				lab = new_key
				gs.append(round(min(np.mean(tmp, axis=0)), 2))
				gs.append(round(max(np.mean(tmp, axis=0)), 2))
				plt.bar(names, np.mean(tmp, axis=0), label=lab, yerr=np.std(tmp, axis=0), color=col[f], width=0.4, capsize=0.8, capstyle='projecting', ecolor=col1[f])
				# plt.bar(names, np.mean(tmp, axis=0), label=lab, yerr=np.std(tmp, axis=0), hatch='xxx', edgecolor=col[f], fill=False, ls='-', lw='0.6', width=0.4, capsize=0.8, capstyle='projecting')
				# plt.bar(names, np.mean(tmp, axis=0), label=new_key, yerr=np.std(tmp, axis=0), color=col[f])
				f += 1

	print("File name:",fname, "| min:", min(gs), "| max:", max(gs))
	plt.ylabel(r'Gain$\%$')
	plt.legend()
	plt.grid(axis='y', ls=":", alpha=0.5)
	if var == "k":
		plt.xlabel(r'$k$')
	else:
		plt.xlabel(r'$|C|$')
	

	if save:
		plt.savefig(fname[:-4], dpi=300, bbox_inches='tight') 
	if show:
		plt.show()
	plt.clf()

	j = 0
	for k in algs.keys():
		for i in range(algs[k], 0, -1):
			new_key = k + "." + str(i)
			lab = new_key
			if new_key == "NA.1":
				lab = "NA"
			# plt.errorbar(x, np.mean(data_tim[j], axis=0), ls="-", marker=markers[i], c=colors[i], label=j, yerr=np.std(data_tim[j], axis=0))
			# plt.errorbar(x, np.mean(data_tim[j], axis=0), ls="-", marker=markers[i], label=j, yerr=np.std(data_tim[j], axis=0))
			plt.errorbar(x, np.mean(data_tim[new_key], axis=0), ls="-", label=lab, marker=markers[j], yerr=np.std(data_tim[new_key], axis=0))
			# sns.pointplot(x, np.mean(data_tim[j], axis=0), label=j)
			# plt.plot(names, np.mean(data_tim[j], axis=0), markers[i] + "-" + colors[i], label=j)
			j += 1

	# plt.axis([min(x), max(x), 0, max_ts * 1.1])
	plt.ylabel(r'Time$(s)$')
	plt.legend()
	plt.grid(ls=":", alpha=0.5)
	if var == "k":
		plt.xlabel(r'$k$')
	else:
		plt.xlabel(r'$|C|$')


	if save:
		plt.savefig(fname[:-4] + "_time", dpi=300, bbox_inches='tight') 
	if show:
		plt.show()
	plt.clf()

	infile.close()

def plot_exp2(fname):
	infile = open(fname, "r")
	reps = 0
	reps1 = 0
	rounds = 0
	for line in infile:
		if line[0].isdigit():
			rounds += 1
		if line[:-1] == "EG" or line[:-1] == "NA":
			reps1 += 1
		if line[:-1] == "EG":
			reps += 1

	rounds = rounds//reps1

	eg_gains = np.empty([reps, rounds-1])
	eg_time = np.empty([reps, rounds-1])
	na_gains = np.empty([reps, rounds-1])
	na_time = np.empty([reps, rounds-1])

	infile.seek(0)
	for i in range(reps):
		for k in range(3):
			infile.readline()
		for j in range(rounds - 1):
			eg_gains[i][j] = float(infile.readline().split()[1])
		infile.readline()
		infile.readline()
		na_gains[i][0] = float(infile.readline().split()[1])
		for j in range(1, rounds - 1):
			na_gains[i][j] = na_gains[i][0]
			infile.readline()

	final = 100*np.subtract(eg_gains, na_gains)/na_gains
	names = []
	for i in range(rounds - 1):
		names.append(str(i + 1))

	plt.ylabel(r'Gain$\%$')
	plt.xlabel('Rounds')
	plt.grid(axis='y', ls=":", alpha=0.5)
	# plt.bar(names, np.mean(final, axis=0), capsize=0.8, yerr=np.std(final, axis=0), hatch='xxx', fill=False, edgecolor="steelblue", ls='-', lw='0.6', width=0.4)
	plt.bar(names, np.mean(final, axis=0), yerr=np.std(final, axis=0), color="steelblue", width=0.4, capsize=0.8, capstyle='projecting')
	if save:
		plt.savefig(fname[:-4], dpi=300, bbox_inches='tight') 
	if show:
		plt.show()
	plt.clf()
	infile.close()

def plot_exp3(fname):
	infile = open(fname, "r")
	var = "c"
	if "_k_" in fname:
		var = "k"
	xaxis = []
	beta = []
	util = []
	DSS = []
	infile.readline()
	for line in infile:
		l = line.split()
		xaxis.append(int(l[0]))
		util.append(float(l[1]))
		beta.append(float(l[2]))
		DSS.append(float(l[3]))

	if scatter:
		xaxis = [str(i) for i in xaxis]
		col = ["darkorange", 'dimgrey', "steelblue"]
		plt.plot(xaxis, util, ls="-", marker=markers[0], color=col[0], label=r'avg[$\sigma_j (\mathcal{S}^j)$]')
		plt.plot(xaxis, beta, ls="-", marker=markers[1], color=col[2], label=r'avg[$\beta_j (\mathcal{S}^j)$]')
		plt.ylabel(r'$\sigma (S)$')
		plt.xticks(xaxis)

	if bar:
		plt.ylabel(r'Violation of basicness$\%$')
		diff = []
		names = [str(i) for i in xaxis]
		for i in range(len(xaxis)):
			diff.append(100*(util[i] - beta[i])/util[i])
		plt.bar(names, diff)

	plt.legend()
	plt.grid(True)
	if var == "k":
		plt.xlabel(r'$k$')
	else:
		plt.xlabel(r'$|C|$')

	if save:
		plt.savefig(fname[:-4], dpi=300, bbox_inches='tight') 
	if show:
		plt.show()
	plt.clf()
	infile.close()

def plot(fname):
	if "utility" in fname:
		plot_exp3(fname)
	elif "gain" in fname or "ic_all_rounds" not in fname and ("c_all_rounds" in fname or "k_all_rounds" in fname or "mult_rounds" in fname):
		plot_exp1(fname)
	else:
		plot_exp2(fname)

def plot_old_exp2(fname):
	infile = open(fname, "r")
	reps = 0
	rounds = 0
	for line in infile:
		if line[0].isdigit():
			rounds += 1
		if line[:-1] == "EG":
			reps += 1
	rounds = rounds//reps

	gains = np.empty([reps, rounds-2])
	time = np.empty([reps, rounds-2])

	infile.seek(0)
	for i in range(reps):
		for k in range(3):
			infile.readline()[:-1]
		r1 = float(infile.readline().split()[1])
		for j in range(rounds - 2):
			r2 = float(infile.readline().split()[1])
			gains[i][j] = 100*(r2 - r1)/r1
			r1 = r2

	names = []
	for i in range(rounds - 2):
		names.append(str(i + 1) + "-" + str(i + 2))

	plt.ylabel(r'Gain$\%$')
	plt.xlabel('Rounds')
	plt.xticks(rotation=30)
	plt.grid(axis='y')
	plt.bar(names, np.mean(gains, axis=0), yerr=np.std(gains, axis=0))
	if save:
		plt.savefig(fname[:-4], dpi=300, bbox_inches='tight') 
	if show:
		plt.show()
	plt.clf()
	infile.close()

def copyfile(infile_name, outfile):
	infile = open(infile_name, "r")
	outfile.write(infile.read())
	infile.close()

def main():
	# f = create_fname(dataset, clust, model, similarity, rand_edges, mu, sigma, variable, exp_type)
	# if run_box:
	# 	os.system(box_command(loc))
	# cmd = alg_command(algs, dataset, model, variable, exp_type, mc1, mc2, out, rnd, clust, box, similarity, rand_edges, mu, sigma, budget, competitors)
	# print(cmd)
	if run_clust:
		cmd = clust_command(calg, dataset, k, p, w, eps_db, eps_dgcd, gamma, min_samples)
		os.system(cmd)

	# clust_exper("kmeans", "g", 0, clusts=[8,16, 32, 64, 128])
	# plot_exp5_new("g", "kmeans")
	

	# print(f)
	# if not os.path.exists(f) or run_alg:
	# 	cmd = alg_command(algs, dataset, model, variable, exp_type, mc1, mc2, out, rnd, clust, box, similarity, rand_edges, mu, sigma)
	# 	print(cmd)
	# 	os.system("make")
	# 	os.system(cmd)

	# if exp_type == 1:
	# 	plot_exp1(f)
	# elif exp_type == 2:
	# 	plot_exp2(f)

	# ans = input("Would you like to open the data file? (y,n): ")
	# if ans == "y":
	# 	os.system("open " + f)



	###############################################
	################ Experiments ##################
	###############################################
	# scalability_exp('f', 5, [100000, 500000, 1000000, 1500000])
	# for v in ["c"]:
	# plot_exp1_new("g", "k", "kmeans")
	# plot_exp1_new("g", "c", "kmeans")
	# plot_exp1_new("f", "k", "kmeans")
	# plot_exp1_new("f", "c", "kmeans")
	# plot_exp1_new("g", "k", "dgcd")
	# plot_exp1_new("g", "c", "dgcd")
	# plot_exp1_new("f", "k", "dgcd")
	# plot_exp1_new("f", "c", "dgcd")
	# plot_exp1_new("x", "k", "original")
	# plot_exp1_new("x", "c", "original")
	# plot_exp1_new("l", "k", "original")
	# plot_exp1_new("l", "c", "original")

	# plot_exp2_new("l", "original")
	# plot_exp2_new("g", "kmeans")

	# plot_exp3_new("g", "k", "kmeans", 'r')
	# plot_exp3_new("g", "c", "kmeans", 'r')
	# plot_exp3_new("l", "k", "original", 'r')
	# plot_exp3_new("l", "c", "original", 'r')
	# plot_exp3_new("g", "k", "dgcd", 'r')
	# plot_exp3_new("g", "c", "dgcd", 'r')
	# plot_exp3_new("x", "k", "original", 'r')
	# plot_exp3_new("x", "c", "original", 'r')
	# plot_exp3_new("g", "k", "kmeans", 'b')
	# plot_exp3_new("g", "c", "kmeans", 'b')
	# plot_exp3_new("l", "k", "original", 'b')
	# plot_exp3_new("l", "c", "original", 'b')
	# plot_exp3_new("g", "k", "dgcd", 'b')
	# plot_exp3_new("g", "c", "dgcd", 'b')
	# plot_exp3_new("x", "k", "original", 'b')
	# plot_exp3_new("x", "c", "original", 'b')

	scalability_table("f")

	

	# Exper 1	
	if False: 
		for d in ['g']:
			for model in ["ic", "lt"]:
				for v in ["c", "k"]:
					f_gain =  open("results/gain/" + clustalg + "/" + d + "_" + model + "_" + v + "_gain.dat", "w")
					for i in range(r):
						f = create_fname(d, clustalg, model, similarity, rand_edges, mu, sigma, v, 4)
						cmd = alg_command(algs, d, model, v, 4, mc1, mc2, out, rnd, clustalg, box, similarity, rand_edges, mu, sigma, budget, competitors)
						print(cmd)
						os.system(cmd)
						copyfile(f, f_gain)
					f_gain.close()
						
	
	#Exper 2
	clustalg = 'original'
	if False:
		comps = [16]
		budjs = [10]
		for d in ['f', 'g']:
			for model in ["ic", "lt"]:
				all_rounds = open("results/rounds/" + clustalg + "/" + d + "_" + model + "_all_rounds.dat", "w")
				for i in range(r):
					for c in comps:
						for b in budjs:
							f = create_fname(d, clustalg, model, similarity, rand_edges, mu, sigma, variable, 2)
							cmd = alg_command(algs, d, model, variable, 2, mc1, mc2, out, rnd, clustalg, box, similarity, rand_edges, mu, sigma, budget, competitors)
							os.system(cmd)
							ifile = open(f, "r")
							all_rounds.write("c=" + str(c) + " k=" + str(b) + "\n")
							all_rounds.write(ifile.read())
							ifile.close()
				all_rounds.close()

	
	#Exper 3
	clustalg = 'kmeans'
	if False:
		for d in ['g']:
			for model in ["ic", "lt"]:
				for v in ["k", "c"]: 
					cmd = alg_command(algs, d, model, v, 3, mc1, mc2, out, rnd, clustalg, box, similarity, rand_edges, mu, sigma, budget, competitors)
					os.system(cmd)
					f = create_fname(d, clustalg, model, similarity, rand_edges, mu, sigma, v, 3)
					os.system("mv " + f + " " + f[:f.find("/")] + "/utility" + f[f.find("/"):])


	#Dataset expers
	if False:
		for d in ['x']:
			for model in ["ic", 'lt']:
				cmd = alg_command("NA", d, model, "k", 1, mc1, mc2, out, rnd, clust, box, similarity, rand_edges, mu, sigma, budget, "1")
				os.system(cmd)


	#Scalability exper
	if False:
		clustalg = 'ff'
		for d in ['f']:
			for siz in [500000, 1000000, 1500000, 0]:
				cmd = clust_command(clustalg, d, k, p, w, eps_db, eps_dgcd, gamma, min_samples, siz)
				os.system(cmd)
				for model in ["ic", "lt"]:
					outf = open("results/scalability/" + d + "_" + model + "_scalability.dat", "a")
					cmd = clust_command(clustalg, d, k, p, w, eps_db, eps_dgcd, gamma, min_samples, siz)
					os.system(cmd)
					for i in range(r):
						f = create_fname(d, clustalg, model, similarity, rand_edges, mu, sigma, variable, 2)
						cmd = alg_command(algs, d, model, variable, 2, mc1, mc2, out, rnd, clustalg, box, similarity, rand_edges, mu, sigma, budget, competitors)
						os.system(cmd)
						ifile = open(f, "r")
						outf.write(str(siz) + "\n")
						outf.write(ifile.read())
						outf.flush()
						ifile.close()




	###############################################
	################## Plots ######################
	###############################################



	################exp1###########################
	# plot_all_rounds("g_all_rounds.dat", 10)
	if False:
		for d in ["r"]:
			for q in ["lt", "ic"]:
				for v in ["k", "c"]: 
					plot("results/gain/" + clustalg + "/" + d + "_" + q + "_" + v + "_gain.dat")



	###############exp2############################
	if False:
		for q in ["lt", "ic"]:
			# for d in ["g", "f"]:
			for d in ["x"]:
				print("results/rounds/" + clustalg + "/" + d + "_" + q + "_all_rounds.dat")
				plot("results/rounds/" + clustalg + "/" + d + "_" + q + "_all_rounds.dat")


	#################exp3##########################
	if False:
		for q in ["lt", "ic"]:
			for d in ["g", "f"]:
				plot("results/utility/" + clustalg + "/" + d + "_clust_" + q + "_ysim_k_utility.dat")
				plot("results/utility/" + clustalg + "/" + d + "_clust_" + q + "_ysim_c_utility.dat")


	if False:
		for d in ['f', 'g']:
			# plot_exp1_new(d, "k")
			# plot_exp1_new(d, "c")
			# plot_exp2_new(d)
			plot_exp3_new(d, "k")
			plot_exp3_new(d, "c")



main()
# ./CGW -alg EGNA -dt test -m 1 -c 2 -k 2 -mc1 100 -mc2 1000 -rnd 2 -o def -exp 4






