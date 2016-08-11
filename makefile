
CFLAGS = `pkg-config --cflags --libs opencv`
TARGET = Tortoise
#~ OBJS = main.o tortoise.o
SRC = main.cpp Tortoise.cpp GeneralHough.cpp
DEPS = GeneralHough.hpp Tortoise.hpp

$(TARGET) : $(SRC)
	g++ -o $(TARGET) $(SRC) -Wall -g -std=c++11 $(CFLAGS) -lX11

#~ $(TARGET) : $(OBJS)
#~ 	g++ -o $(TARGET) $(OBJS) -Wall -g -std=c++11 $(CFLAGS) -lX11

#~ main.o: main.cpp Tortoise.hpp
#~ 	g++ -c main.cpp -o main.o -std=c++11

#~ tortoise.o : Tortoise.cpp Tortoise.hpp
#~ 	g++ -c Tortoise.cpp -o Tortoise.o -std=c++11

#~ tortoise.o : GeneralHough.cpp GeneralHough.hpp
#~ 	g++ -c GeneralHough.cpp -o GeneralHough.o -std=c++11



