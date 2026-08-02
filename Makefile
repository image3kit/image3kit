-include .env

VENV_DIR ?= $(shell \
	if [ -d .venv ]; then echo .venv; \
	elif [ -d ../.venv ]; then echo ../.venv; \
	elif [ -d ../../.venv ]; then echo ../../.venv; \
	else echo .venv; fi)

VENV_BIN := $(VENV_DIR)/bin
PYTHON   := $(VENV_BIN)/python
PIP      := $(VENV_BIN)/pip

help:
	@echo "Install prerequisite using commands:"
	@echo "  python -m venv ./.venv"
	@echo "  source .venv/bin/activate"
	@echo "  pip install pybind11-stubgen numpy pytest"
	@echo
	@echo "Available target commands:"
	@echo "  make all          # Install and update pybind11 stubs"
	@echo "  make install      # Install the package"
	@echo "  make test         # Run tests"

all: setup_venv
	$(PIP) uninstall -y image3kit || true
	$(PIP) install -e .
	@$(MAKE) stubgen
	@$(MAKE) pre-commit
	@$(MAKE) test

stubgen:
	@[ -f $(VENV_BIN)/pybind11-stubgen ] || (set -x && $(PIP) install pybind11-stubgen)
	$(VENV_BIN)/pybind11-stubgen image3kit --output-dir src
	sed -i 's/sirun.dbl3/sirun.dbl3\|tuple\[float, float, float\]/g' src/image3kit/_core/voxlib.pyi
	sed -i 's/sirun.int3/sirun.int3\|tuple\[int, int, int\]/g'       src/image3kit/_core/voxlib.pyi
	sed -i 's/ = None,/\|None = None,/g'       src/image3kit/_core/voxlib.pyi


install_global:
	pip install .

test:
	@[ -f $(VENV_BIN)/pytest ] || (set -x && $(PIP) install pytest)
	$(PYTHON) -m pytest
	@if [ -d build ]; then $(VENV_BIN)/cmake --build build --target common_tests voxlib_tests && $(VENV_BIN)/ctest --test-dir build --output-on-failure; fi

.PHONY: setup_venv clean build tests all help stubgen install_global test pre-commit

setup_venv:
	@if [ ! -f $(PYTHON) ]; then \
		echo "Creating virtual environment at .venv..."; \
		python3 -m venv .venv; \
	fi
	$(PIP) install cmake pre-commit
	$(VENV_BIN)/cmake -S . -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DBUILD_TESTING=ON
	ln -sf build/compile_commands.json ./

clean:
	rm -rf build compile_commands.json

pre-commit:
	$(VENV_BIN)/pre-commit run -a || $(VENV_BIN)/pre-commit run -a
