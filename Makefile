PARA = -std=c++11 -Wall -O3 -w

GCW: GLib.hpp mappedHeap.hpp HeapData.hpp hypergraph.hpp GCW.cpp Utilities.o Locations.o BinaryHeap.o rwgraph.o option.o sfmt/SFMT.c
	g++ -o GCW $(PARA) GCW.cpp Utilities.o Locations.o BinaryHeap.o rwgraph.o option.o sfmt/SFMT.c
# 	g++ -o GCW -fopenmp $(PARA) GCW.cpp Utilities.o Locations.o BinaryHeap.o rwgraph.o option.o sfmt/SFMT.c

el2bin: el2bin.cpp
	g++ el2bin.cpp -o el2bin $(PARA)

option.o: option.cpp option.h
	g++ -c option.cpp -o option.o $(PARA)

BinaryHeap.o: BinaryHeap.cc BinaryHeap.h
	g++ -c BinaryHeap.cc -o BinaryHeap.o $(PARA)

Utilities.o: Utilities.cpp Utilities.h
	g++ -c Utilities.cpp -o Utilities.o $(PARA)

Locations.o: Locations.cpp Locations.h
	g++ -c Locations.cpp -o Locations.o $(PARA)

rwgraph.o: rwgraph.cpp	rwgraph.h
	g++ -c rwgraph.cpp -o rwgraph.o $(PARA)


# com: rwgraph.cpp GCW.cpp el2bin.cpp option.cpp GLib.hpp
# 	g++ -c GLib.hpp -o GLib.o $(PARA)
# 	g++ -c mappedHeap.hpp -o mappedHeap.o $(PARA)
# 	g++ -c HeapData.hpp -o HeadData.o $(PARA)
# 	g++ -c option.cpp -o option.o $(PARA)
# 	g++ -c BinaryHeap.cc -o BinaryHeap.o $(PARA)
# 	g++ -c Utilities.cpp -o Utilities.o $(PARA)
# 	g++ -c Locations.cpp -o Locations.o $(PARA)
# 	g++ -c rwgraph.cpp -o rwgraph.o $(PARA)
# 	g++ GCW.cpp Utilities.o Locations.o BinaryHeap.o rwgraph.o option.o -o GCW -fopenmp $(PARA) sfmt/SFMT.c
# 	g++ nDSSA.cpp Utilities.o Locations.o BinaryHeap.o rwgraph.o option.o -o nDSSA -fopenmp $(PARA) sfmt/SFMT.c
	

clean:
	rm *.o GCW el2bin
