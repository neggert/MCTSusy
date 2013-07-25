INCLUDES=-I/usr/include/python2.7


all:  example

example: mt2_bisect.cpp example.cpp mt2_bisect.h 
	g++ -o example mt2_bisect.cpp example.cpp 

tar:
	tar -cvzf mt2-1.01a.tar.gz mt2_bisect.cpp mt2_bisect.h Makefile example.cpp README update_notes

mt2_bisect_wrap.cxx: mt2_bisect.i
	swig -c++ -python mt2_bisect.i

mt2_bisect.o: mt2_bisect.cpp mt2_bisect.h
	g++ -fpic -c mt2_bisect.cpp

mt2_bisect_wrap.o: mt2_bisect_wrap.cxx
	g++ -fpic -c mt2_bisect_wrap.cxx $(INCLUDES)

_mt2_bisect.so: mt2_bisect.o mt2_bisect_wrap.o
	g++ -shared mt2_bisect_wrap.o mt2_bisect.o -o _mt2_bisect.so -lpython2.7

clean:
	rm -f *.so *.o

.PHONY: clean