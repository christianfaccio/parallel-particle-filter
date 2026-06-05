# Makefile for the parallel particle filter
#
# Targets:
#   make            build the serial filter (bin/pf_serial)
#   make run        build then run on the data in ./data
#   make debug      build with -O0 -g for gdb / sanitizers
#   make asan       build with AddressSanitizer (leak + OOB checks)
#   make clean      remove build artifacts
#
# Layout:
#   include/   public headers
#   src/       translation units
#   build/     object files (generated)
#   bin/       executables (generated)
#
# The serial build is the correctness baseline; the OpenMP / CUDA backends
# will be added as additional targets below (see "Parallel backends").

# ---- Toolchain & flags -------------------------------------------------------

CC       := cc
CSTD     := -std=c11
WARN     := -Wall -Wextra
OPT      := -O3 -march=native
CPPFLAGS := -Iinclude -Ithird_party/lidarsim
CFLAGS   := $(CSTD) $(WARN) $(OPT)
LDLIBS   := -lm

# ---- Sources & objects -------------------------------------------------------

SRC_DIR   := src
BUILD_DIR := build
BIN_DIR   := bin

SOURCES := $(wildcard $(SRC_DIR)/*.c)
OBJECTS := $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SOURCES))

TARGET := $(BIN_DIR)/pf_serial

# ---- Default target ----------------------------------------------------------

.PHONY: all
all: $(TARGET)

$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(CC) $(CFLAGS) $(OBJECTS) $(LDLIBS) -o $@

# Compile each .c -> .o, regenerating when its headers change (see deps below).
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR) $(BIN_DIR):
	mkdir -p $@

# ---- Convenience targets -----------------------------------------------------

.PHONY: run
run: $(TARGET)
	./$(TARGET)

# Unoptimized build with symbols, for gdb.
.PHONY: debug
debug: CFLAGS := $(CSTD) $(WARN) -O0 -g
debug: clean $(TARGET)

# AddressSanitizer build: catches leaks, use-after-free, out-of-bounds.
.PHONY: asan
asan: CFLAGS := $(CSTD) $(WARN) -O1 -g -fsanitize=address,undefined
asan: LDLIBS := -lm -fsanitize=address,undefined
asan: clean $(TARGET)

# ---- Parallel backends (to be added) -----------------------------------------
# When the OpenMP backend lands, add e.g.:
#   bin/pf_omp: CFLAGS += -fopenmp
#   bin/pf_omp: LDLIBS += -fopenmp
# guarded behind its own object dir so serial and OpenMP objects don't collide.

# ---- Housekeeping ------------------------------------------------------------

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

# Auto-generated header dependencies (-MMD). Silently ignored on first build.
-include $(OBJECTS:.o=.d)
