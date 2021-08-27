CC=mpicc
EXE=collectives
OUT=data*.txt
PLOTS=plot*.jpg
HOST=hostfile

.PHONY=all clean

all : src.c
	$(CC) src.c -o $(EXE)

clean :
	rm -f $(EXE) $(OUT) $(HOST) $(PLOTS)