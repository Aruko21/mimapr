# Exec file name
TARGET = lab_pa9

# Plot file
PLOT = results.gp

# cc invoking
CC = g++

# object file names
OBJECTS = src/main.o

# source files directory
SRCS = src/

# extern libraries
LIBS = -lm

STD = c++14

# project dependence
all: $(TARGET)
$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGET) $(OBJECTS) $(LIBS)

# $< means first dependent argument. Example: all: obj obj2 - $< = obj
# $@ means target. For previous example this is "all"
# rule for implicit handling of all object files
%.o: %.c
	$(CC) -c $< -o $@ -std=${STD}
plot:
	gnuplot ./$(PLOT)
clean:
	rm -rf $(TARGET) $(SRCS)*.o