
CFLAGS = `pkg-config --cflags --libs opencv`
TARGET = tortoise
SRC = main.cpp
DEPS = my_img_proc.hpp general_hough.hpp

$(TARGET) : $(SRC) $(DEPS)
	g++ -Wall -g $(SRC) -o $(TARGET) $(CFLAGS)


