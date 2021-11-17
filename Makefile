all: abm.so

abm.so: abm.c
	cc -fPIC -shared -Ofast -lm -o abm.so abm.c