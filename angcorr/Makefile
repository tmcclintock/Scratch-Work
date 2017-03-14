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
EXEC = c_angcorr.so
CFLAGS = -fPIC
OFLAGS = -shared
endif

INCL = -I/${GSLI} -fopenmp -O3

LIBS = -lm -lgsl -L/${GSLL} -O2

.SUFFIXES:.c.o

%.o: %.c
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@

$(EXEC): $(OBJS)
	$(CC) $(OFLAGS) $(OBJS) $(LIBS) -o $(EXEC)
	rm $(OBJS)

.PHONY: clean

clean:
	rm -f $(OBJS) $(EXEC) *~ *.exe *.so
