FLAGS=-O3 -std=c++11 -mcmodel=large -L /lib64 
EXEC=closest-pair
LIBS     = -lusb-1.0 -l pthread
LDFLAGS="-Wl,--copy-dt-needed-entries"

all: $(EXEC)

$(EXEC):
	$(CXX) $(FLAGS) $(EXEC).cpp -c -o $(EXEC).o
	$(CXX) $(FLAGS) $(EXEC).o -o $(EXEC)

clean:
	rm -rf $(EXEC) *.o
