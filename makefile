
CFLAGS = `pkg-config --cflags --libs opencv`
TARGET = tortoise
SRC = main.cpp
DEPS = my_img_proc.hpp my_general_hough.hpp

$(TARGET) : $(SRC) $(DEPS)
	g++ -Wall $(SRC) -o $(TARGET) $(CFLAGS) 


