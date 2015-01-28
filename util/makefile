CC     := gcc
CFLAGS := -I../util -Wall -g
VPATH  := ../util

OBJ := doubly_linked_list.o \
       epsilon.o            \
       memfunc.o            \
       quantifier.o

all: $(OBJ)

doubly_linked_list.o: doubly_linked_list.h \
                      all.h

epsilon.o: epsilon.h

memfunc.o: memfunc.h \
           all.h

quantifier.o: quantifier.h \
              all.h

clean:
	rm $(OBJ)
