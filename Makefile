all: abm.so

abm.so: abm.c
	gcc -fPIC -shared -O2 -lm -o abm.so abm.c

clean:
	rm abm.so