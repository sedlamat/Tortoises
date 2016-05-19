
CFLAGS = `pkg-config --cflags --libs opencv`
TARGET = tortoise
SRC = main.cpp
DEPS = my_img_proc.hpp

$(TARGET) : $(SRC) $(DEPS)
	g++ -Wall $(SRC) $(DEPS) -o $(TARGET) $(CFLAGS) 


