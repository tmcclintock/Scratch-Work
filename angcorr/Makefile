CC = gcc
ifdef ALONE
ifeq ($(ALONE),yes)
$(info Building executable)
OBJS = angcorr.o
EXEC = angcorr.exe
CFLAGS = -g
OFLAGS = 
endif
else
$(info Building shared library)
OBJS = angcorr.o
EXEC = _angcorr.so
CFLAGS = -fPIC
OFLAGS = -shared
endif

INCL = -I/${GSLI} -O3

LIBS = -lm -lgsl -lgslcblas -L/${GSLL} -O3

.SUFFIXES:.c.o

%.o: %.c
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@

$(EXEC): $(OBJS)
	$(CC) $(OFLAGS) $(OBJS) $(LIBS) -o $(EXEC)
	rm $(OBJS)

.PHONY: clean

clean:
	rm -f $(OBJS) $(EXEC) *~ *.exe *.so
