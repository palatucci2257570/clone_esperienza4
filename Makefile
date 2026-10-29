# Compilare:                 make
# Cancellare gli eseguibili: make clean

CFLAGS = -std=c17 -Wall -Wextra -Wno-unused-parameter -O2

all: oscillatore pendolo_smorzato

oscillatore: oscillatore.c
	cc $(CFLAGS) oscillatore.c -o oscillatore -lm

pendolo_smorzato: pendolo_smorzato.c
	cc $(CFLAGS) pendolo_smorzato.c -o pendolo_smorzato -lm

clean:
	rm -f oscillatore pendolo_smorzato
