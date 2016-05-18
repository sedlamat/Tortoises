
CFLAGS = $(shell pkg-config --cflags opencv)
LIBS = $(shell pkg-config --libs opencv)
OBJS = main.o #tortoise.o generalHough.o

tortoise : $(OBJS)
	g++ $(OBJS) -o tortoise $(LIBS) $(CFLAGS) 

main.o : main.cpp
	g++ -c main.cpp

#tortoise.o : tortoise.cpp
#	g++ -c tortoise.cpp

#generalHough.o : generalHough.cpp
#	g++ -c generalHough.cpp

clean:
	rm tortose $(OBJS)
