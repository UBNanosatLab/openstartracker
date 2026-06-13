CC = g++
PYTHON ?= python3
PYTHONHEADERS ?= $(shell $(PYTHON) -c 'import sysconfig; print(sysconfig.get_paths()["include"])')

all: $(OBJ)
	swig -python -c++ beast.i
	$(CC) -g -std=c++11 -Ofast -fPIC -c beast_wrap.cxx -o beast_wrap.o -lstdc++ -I$(PYTHONHEADERS)
	$(CC) -g -shared -fPIC beast_wrap.o -o _beast.so
	
clean:
	rm -f _beast.so beast.so beast.o beast.py beast.pyc beast_wrap.cxx beast_wrap.o
