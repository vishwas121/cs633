CC=mpicc
EXE=code
OUT=output.txt
HOST=hostfile
DUMP=dump
PLOTS=plot_*.jpg

.PHONY=all clean

all : src.c
	$(CC) src.c -o $(EXE)

clean :
	rm -f $(EXE) $(OUT) $(HOST) $(DUMP) $(PLOTS)

