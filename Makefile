all: abm.so

abm.so: abm.c
	cc -fPIC -shared -O2 -lm -o abm.so abm.c