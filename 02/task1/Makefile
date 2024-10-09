CC = gcc
CFLAGS = -Wall -Werror -std=gnu11 -O3
TARGETS = pi_seq

.PHONY: all clean

all: $(TARGETS)

clean:
	$(RM) $(TARGETS)

%: %.c
	$(CC) $(CFLAGS) $< -o $@
