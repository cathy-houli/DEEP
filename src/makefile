all: main

Compiler	= gcc
FLAGS		= -g -Wall -O0
LIB		= -lz -lm -lpthread
HEADER		= main.h
SOURCE		= main.c index.c bwa_index.c stat.c NW.c variation_hash.c variation_align.c tool.c map_align.c
PROG		= DEEP
OBJECT		= $(SOURCE:%.c=%.o) 

all:		$(PROG)

main:		$(OBJECT)
			$(Compiler) $(FLAGS) $(OBJECT) -o $(PROG) $(LIB)
                        
%.o:		%.c $(HEADER)
			$(Compiler) $(FLAGS) -c $<

clean:
		rm -f *.o
