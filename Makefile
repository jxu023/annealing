CFLAGS=-pedantic -Wall -g -std=c++14
SOURCES=anneal.cpp

a: $(SOURCES)
	g++ -o anneal $(CFLAGS) $(SOURCES)
clean:
	rm -rf anneal *.o
run:
	#./anneal 4-Insertions_3.col 13 -fixed_k
	#./anneal input1.txt 13 -fixed_k
	#./anneal input1.txt 13 -penalty -st
	#./anneal 2-FullIns_5.col 13 -fixed_k -st
	./anneal 4-Insertions_3.col 100 -penalty
