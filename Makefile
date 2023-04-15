

PROGS = basic-least-squre gram-schmidt-qr qr-for-ls
ROOT = .

CC = gcc
LDFLAGS=
LDDIR=-L$(ROOT)/lib
LDLIBS=$(LDDIR) -lm $(EXTRALIBS)
CFLAGS= -std=c99 -g -I$(ROOT)/include -Wall -DLINUX -D_GNU_SOURCE $(EXTRA)
AR=ar

all:	$(PROGS)

%:	%.c
	$(CC) $(CFLAGS) $@.c -o $@ $(LDFLAGS) $(LDLIBS)

naive-gauss:    basic-least-square.c
		$(CC) $(CFLAGS) basic-least-square.c -o basic-least-square $(LDFLAGS) $(LDLIBS)

gram-schmidt-qr:    gram-schmidt-qr.c
		$(CC) $(CFLAGS) gram-schmidt-qr.c -o gram-schmidt-qr $(LDFLAGS) $(LDLIBS)

qr-for-ls:    qr-for-ls.c
		$(CC) $(CFLAGS) qr-for-ls.c -o qr-for-ls $(LDFLAGS) $(LDLIBS)

clean:
		rm -f $(PROGS)

