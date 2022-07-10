CC=gcc
CFLAGS= -g -O3
DEPS =
OBJ = teste.o
LIB = -fopenmp

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< $(LIB)

teste: $(OBJ)
	gcc $(CFLAGS) -o $@ $^ $(LIB)

clean:
	rm *.o teste