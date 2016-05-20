
CFLAGS = `pkg-config --cflags --libs opencv`
TARGET = tortoise
SRC = main.cpp

$(TARGET) : $(SRC) $(DEPS)
	g++ -Wall $(SRC) -o $(TARGET) $(CFLAGS) 


