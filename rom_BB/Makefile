
CC = mpic++
CFLAGS = -Wall -Wextra -pedantic -O3
LDFLAGS =
INCLUDES = -I/usr/local/include
LIBS = -L/usr/local/lib -lm

TARGET_LIB = librom_BB.a
OBJS = read_BB.o write_BB.o std.o

$(TARGET_LIB): $(OBJS)
	ar rcs $@ $?
	rm -f $(OBJS)

clean:
	-rm -f $(OBJS) $(TARGET_LIB)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(INCLUDES)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(INCLUDES)

.PHONY: clean;
	$(RM) -f $(OBJS) $(TARGET_LIB)

clean_data:
	$(RM) -f $(OBJS)
