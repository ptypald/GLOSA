CC := gcc
CFLAGS := -O2 -std=c11
LDLIBS := -lm
BUILD_DIR := build

.PHONY: all clean

all: \
	$(BUILD_DIR)/glosa_dp \
	$(BUILD_DIR)/glosa_sdp \
	$(BUILD_DIR)/glosa_dddp \
	$(BUILD_DIR)/glosa_ddp

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/glosa_dp: red2green/dp/main.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

$(BUILD_DIR)/glosa_sdp: red2green/journal\ TRC/sdp.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) "$<" $(LDLIBS) -o $@

$(BUILD_DIR)/glosa_dddp: red2green/journal\ TRC/dddp.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) "$<" $(LDLIBS) -o $@

$(BUILD_DIR)/glosa_ddp: red2green/journal\ TRC/ddp.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) "$<" $(LDLIBS) -o $@

clean:
	rm -rf $(BUILD_DIR)
