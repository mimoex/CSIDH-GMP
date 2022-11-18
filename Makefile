SRC=main.cpp CSIDH.cpp Montgomery.cpp
OBJ=$(SRC:.cpp=.o)
DEP=$(SRC:.cpp=.d)


TARGET=csidh
CFLAGS+=-O3 -Wall -Wextra -DNDEUBG -g
LDFLAGS=-lgmp -lgmpxx

all: $(TARGET)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CFLAGS) -MMD -MP -MF $(@:.o=.d)

-include $(DEP)

$(TARGET): $(OBJ)
	$(CXX) -o $@ $? $(LDFLAGS)

clean:
	$(RM) -rf $(OBJ) $(DEP) $(TARGET)

.PHONY: clean

.SECONDARY: $(DEP)
