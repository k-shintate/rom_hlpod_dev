CC = mpic++
FC = mpif90

CFLAGS   := -g -O3
LDFLAGS  :=
INCLUDES := -I/usr/local/include -I./include -I./../../../test_thermal/submodule/monolis/include
LIBS     := -L/usr/local/lib -lm -lstdc++ \
            -L../../../test_thermal/submodule/monolis/lib -lmonolis_solver \
            -lgedatsu -lmonolis_utils -lmetis -lscalapack -llapack -lblas

INCLUDES_BB = -I../../../test_thermal/include
LIBS_BB     = -L../../../test_thermal/lib -lBB -lBBFE_std -lBBFE_sys -lBBFE_elemmat -lBBFE_manusol

INCLUDES_ROM = -I./../../include/rom_BB

LIBS_ROM     = -L./../../lib -lrom_BB

SRC_DIR := ./src
SRC     := $(wildcard $(addsuffix /*.c,$(SRC_DIR)))
OBJ1    := $(patsubst %.c,%.o,$(SRC)) ./src_main/main.o
OBJ2    := $(patsubst %.c,%.o,$(SRC)) ./src_main/main_bc.o

TARGET1 := merge_graph
TARGET2 := merge_graph_bc

# Default make target
all: build clean_data

# Build targets
build: $(TARGET1) $(TARGET2)

$(TARGET1): $(OBJ1)
	$(FC) $(CFLAGS) $(INCLUDES) $(INCLUDES_BB) $(INCLUDES_ROM) $(LDFLAGS) -o $@ $^ -Wl,--start-group $(LIBS_BB) $(LIBS_ROM) $(LIBS) -Wl,--end-group

$(TARGET2): $(OBJ2)
	$(FC) $(CFLAGS) $(INCLUDES) $(INCLUDES_BB) $(INCLUDES_ROM) $(LDFLAGS) -o $@ $^ -Wl,--start-group $(LIBS_BB) $(LIBS_ROM) $(LIBS) -Wl,--end-group

# Compilation rule for .c to .o
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) $(INCLUDES_BB) $(INCLUDES_ROM) -c $< -o $@

.PHONY clean:;
	$(RM) -f $(OBJ1) $(TARGET1) $(OBJ2) $(TARGET2)

clean_data:
	$(RM) -f $(OBJ1) $(OBJ2)
