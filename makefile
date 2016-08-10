
CFLAGS = `pkg-config --cflags --libs opencv`
TARGET = tortoise
OBJS = main.o tortoise.o
SRC = main.cpp tortoise.cpp
#~ DEPS = my_img_proc.hpp general_hough.hpp tortoise.hpp

#~ $(TARGET) : $(OBJS)
#~ 	g++ -o $(TARGET) $(OBJS) -Wall -g -std=c++11 $(CFLAGS) -lX11

$(TARGET) : $(OBJS)
	g++ -o $(TARGET) $(OBJS) -Wall -g -std=c++11 $(CFLAGS) -lX11

main.o: main.cpp tortoise.hpp
	g++ -c main.cpp -o main.o -std=c++11

tortoise.o : tortoise.cpp tortoise.hpp
	g++ -c tortoise.cpp -o tortoise.o -std=c++11



