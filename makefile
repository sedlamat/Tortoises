
CFLAGS = `pkg-config --cflags --libs opencv`
TARGET = tortoise
SRC = *.cpp
DEPS = my_img_proc.hpp general_hough.hpp tortoise.hpp


$(TARGET) : $(SRC) $(DEPS)
	g++ -Wall -g -std=c++11 $(SRC) -o $(TARGET) $(CFLAGS) -lX11


