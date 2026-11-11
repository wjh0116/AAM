CXX = g++
CXXFLAGS = -O3 -std=c++17
TARGET = AAM

all: $(TARGET)

MODE_HDRS = HybridSample.h HybridGen.h DecGreedy.h MODE.h AAMS.h

$(TARGET): main.cpp head.h graph.h alias.h AdoptionFunction.h $(MODE_HDRS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) main.cpp

clean:
	rm -f $(TARGET)

.PHONY: all clean
