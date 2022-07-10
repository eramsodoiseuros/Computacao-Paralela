CC=gcc
CFLAGS= -g -O3
DEPS =
OBJ = sequential.o
LIB = -fopenmp

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< $(LIB)

sequential: $(OBJ)
	gcc $(CFLAGS) -o $@ $^ $(LIB)

clean:
	rm *.o sequential