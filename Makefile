SRC=main.cpp CSIDH.cpp Montgomery.cpp fp.cpp mcl.cpp
OBJ=$(SRC:.cpp=.o)
DEP=$(SRC:.cpp=.d)


TARGET=csidh
CFLAGS+=-O3 -Wall -Wextra -DNDEUBG -g -I ./mcl/include
CFLAGS+= -I ./mcl/src
LDFLAGS=-L ./mcl/lib -lmcl -lgmpxx -lgmp

all: $(TARGET)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CFLAGS) -MMD -MP -MF $(@:.o=.d)

-include $(DEP)

$(TARGET): $(OBJ)
	$(CXX) -o $@ $(OBJ) $(LDFLAGS)

clean:
	$(RM) -rf $(OBJ) $(DEP) $(TARGET)

.PHONY: clean

.SECONDARY: $(DEP)
