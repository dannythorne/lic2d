
CC=g++
RM=/bin/rm -f

CFLAGS=
CFLAGS+=-g

lic: ./src/bmp_class.o ./src/lic_class.o
	$(CC) $(CFLAGS) -o lic ./src/lic.cc ./src/bmp_class.o ./src/lic_class.o

apple: ./src/bmp_class.o
	$(CC) $(CFLAGS) ./src/apple.cc ./src/bmp_class.o

lic_class.o: ./src/lic_class.h ./src/lic_class.cc
	$(CC) -c $(CFLAGS) ./src/lic_class.cc

bmp_class.o: ./src/bmp_class.h ./src/bmp_class.cc
	$(CC) -c $(CFLAGS) ./src/bmp_class.cc

clean:
	$(RM) lic.exe *.o *.bmp
