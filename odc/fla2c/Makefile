
CC = gcc

LEXSRC = fla2c.l
PROGRAM = fla2c

all: lex.yy.o
	$(CC) -o $(PROGRAM) lex.yy.o -ll

lex.yy.o: lex.yy.c

lex.yy.c: $(LEXSRC)
	flex $(LEXSRC)
