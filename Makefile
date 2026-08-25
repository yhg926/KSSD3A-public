CC ?= gcc
CSTD ?= -std=gnu11
WARNFLAGS ?= -Wno-format-overflow -Wno-unused-result
OPTFLAGS ?= -O3 -flto
ARCH_FLAGS ?=
OMPFLAGS ?= -fopenmp
CFLAGS ?= $(CSTD) $(WARNFLAGS) $(OPTFLAGS) $(ARCH_FLAGS) $(OMPFLAGS)
LDFLAGS ?=
LDLIBS ?= -lz -lm

ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
SRCDIR := $(ROOT)/src
OBJDIR := $(ROOT)/obj
BINDIR := $(ROOT)/bin
PRONAME := kssd3a
TARGET := $(BINDIR)/$(PRONAME)
PREFIX := /usr/local
COMPLETION_FILE := $(ROOT)/etc/$(PRONAME).bash
COMPLETIONDIR ?= $(PREFIX)/share/bash-completion/completions

all: $(TARGET)
	@echo "Build completed."

native: clean
	$(MAKE) ARCH_FLAGS="-march=native" all

avx2: clean
	$(MAKE) ARCH_FLAGS="-march=native -mavx2 -mbmi2" all

SRCS := $(wildcard $(SRCDIR)/*.c)
KLIB_SRCS := $(SRCDIR)/command_ani.c \
             $(SRCDIR)/command_composite.c \
             $(SRCDIR)/command_matrix.c \
             $(SRCDIR)/command_operate.c \
             $(SRCDIR)/command_sketch.c \
             $(SRCDIR)/kssdlib_sort.c
KLIB_OBJS := $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(KLIB_SRCS))
KLIB_LIB_SRCS := klib/kstring.c
KLIB_LIB_OBJS := $(patsubst klib/%.c,$(OBJDIR)/%.o,$(KLIB_LIB_SRCS))
GEN_SRCS := $(filter-out $(KLIB_SRCS),$(SRCS))
GEN_OBJS := $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(GEN_SRCS))
OBJS := $(GEN_OBJS) $(KLIB_OBJS) $(KLIB_LIB_OBJS)
DEPS := $(OBJS:.o=.d)

$(TARGET): $(OBJS) | $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

$(KLIB_OBJS): CFLAGS += -Iklib

$(OBJDIR) $(BINDIR):
	mkdir -p $@

$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

$(OBJDIR)/%.o: klib/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) -Iklib -MMD -MP -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS) $(DEPS)

install:
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 $(TARGET) $(DESTDIR)$(PREFIX)/bin/$(PRONAME)
	@echo "Installed $(PRONAME) to $(PREFIX)/bin"

install_completion:
	install -d $(DESTDIR)$(COMPLETIONDIR)
	install -m 644 $(COMPLETION_FILE) $(DESTDIR)$(COMPLETIONDIR)/$(PRONAME)
	@echo "Installed bash completion to $(COMPLETIONDIR)/$(PRONAME)"

uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/$(PRONAME)
	@echo "Removed $(PRONAME) from $(PREFIX)/bin"

install_env:
	@printf '%s\n' 'export PATH="$(BINDIR):$$PATH"'
	@printf '%s\n' 'Add the line above to your shell profile if you want to run $(PRONAME) directly.'

completion:
	@printf '%s\n' 'source "$(COMPLETION_FILE)"'

test-smoke: all
	"$(ROOT)/tests/smoke_main.sh"

test: test-smoke

-include $(DEPS)

.PHONY: all native avx2 clean install install_completion uninstall install_env completion test-smoke test
