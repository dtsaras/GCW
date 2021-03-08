import random
import argparse


# loc_box = "tl" #(rand, median, mean, tr, tl, br, bl)

parser = argparse.ArgumentParser(description='Create events in a tight box area')
parser.add_argument('-loc', '--location', type=str, default="mean", help='The location of the box')
parser.add_argument('-dt', '--dataset', type=str, default="f", help='Which dataset')
parser.add_argument('-dec', '--decimal_spots', type=int, default=2, help='The number of constant decimal spots')

args = parser.parse_args()
loc_box = args.location

r = "data/"

if args.dataset == "g" or args.dataset == "Gowalla" or args.dataset == "gowalla":
	r += "Gowalla/"
	c = "gowalla"
elif args.dataset == "f" or args.dataset == "Foursquare" or args.dataset == "foursquare":
	r += "Foursquare/"
	c = "foursquare"

inf = r + "events_locations_" + c + ".txt"
outf = r + "box_events_locations_" + c + ".txt"
fin = open(inf, 'r')
fout = open(outf, 'w')

fout.write(fin.readline())
p_lat = 0
p_lon = 0
if loc_box == "rand":
	i = int(random.random() * 128)
	for i in range(i):
		fin.readline()
	l = fin.readline().split()
	p_lat = l[1]
	p_lon = l[2]

else:
	lat = []
	lon = []
	for l in fin:
		k = l.split()
		lat.append(float(k[1]))
		lon.append(float(k[2]))

	if loc_box == "mean": 
		p_lat = str(sum(lat)/len(lat))
		p_lon = str(sum(lon)/len(lon))

	elif loc_box == "median": 
		lat.sort()
		lon.sort()
		p_lat = str((lat[63] + lat[64])/2)
		p_lon = str((lon[63] + lon[64])/2)

	elif loc_box == "tr": 
		p_lat = str(max(lat))
		p_lon = str(max(lon))

	elif loc_box == "tl": 
		p_lat = str(min(lat))
		p_lon = str(max(lon))

	elif loc_box == "br": 
		p_lat = str(max(lat))
		p_lon = str(min(lon))

	elif loc_box == "bl": 
		p_lat = str(min(lat))
		p_lon = str(min(lon))

base_lat = p_lat[:p_lat.find(".") + 1 + args.decimal_spots]
base_lon = p_lon[:p_lon.find(".") + 1 + args.decimal_spots]

for i in range(128):
	fout.write(str(i) + " " + base_lat + str(int(random.random() * 1000)) + " " + base_lon + str(int(random.random() * 1000)) + "\n")

fin.close()
fout.close()
