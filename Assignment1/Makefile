CC=mpicc
EXE=halo_exchange
OUT=data*.txt
HOST=hostfile

.PHONY=all clean

all : src.c
	$(CC) src.c -o $(EXE) -lm

clean :
	rm -f $(EXE) $(OUT) $(HOST)
