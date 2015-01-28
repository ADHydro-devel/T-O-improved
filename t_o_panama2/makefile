CC      := gcc
CFLAGS  := -I../util -Wall -g
LDFLAGS := -lm -lX11
VPATH   := ../util

EXE := test_panama
OBJ := test_panama.o           \
       t_o.o                \
       doubly_linked_list.o \
       epsilon.o            \
       memfunc.o           

$(EXE): $(OBJ)

test_exp2.o: t_o.h     \
            epsilon.h \
            all.h

t_o.o: t_o.h                \
       doubly_linked_list.h \
       epsilon.h            \
       memfunc.h            \
       all.h

doubly_linked_list.o: doubly_linked_list.h \
                      all.h

epsilon.o: epsilon.h

memfunc.o: memfunc.h \
           all.h

clean:
	rm $(EXE) $(OBJ)
