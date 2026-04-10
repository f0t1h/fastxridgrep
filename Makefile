CC ?= cc

CPPFLAGS := -D_POSIX_C_SOURCE=200809L
SRC_CFLAGS := -std=c11 -O3 -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes -Wmissing-prototypes
EXT_CFLAGS := -std=c11 -O3 -Iext/aho-corasick

BIN := fastxridgrep
BUILD_DIR := build

SRC_OBJS := $(BUILD_DIR)/main.o
EXT_OBJS := $(BUILD_DIR)/acism.o $(BUILD_DIR)/acism_create.o $(BUILD_DIR)/acism_file.o

.PHONY: all clean

all: $(BIN)

$(BIN): $(SRC_OBJS) $(EXT_OBJS)
	$(CC) $(SRC_OBJS) $(EXT_OBJS) -o $@ -lz

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/main.o: src/main.c | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(SRC_CFLAGS) -c $< -o $@

$(BUILD_DIR)/acism.o: ext/aho-corasick/acism.c | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(EXT_CFLAGS) -c $< -o $@

$(BUILD_DIR)/acism_create.o: ext/aho-corasick/acism_create.c | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(EXT_CFLAGS) -c $< -o $@

$(BUILD_DIR)/acism_file.o: ext/aho-corasick/acism_file.c | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(EXT_CFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(BIN)
