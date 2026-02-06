
help:
	@echo "Install prerequisite using commands:"
	@echo "  python -m venv ./.venv"
	@echo "  source .venv/bin/activate"
	@echo "  pip install pybind11-stubgen numpy  pytest"
	@echo
	@echo "Available target commands:"
	@echo "  make all # Install and update pybind11 stubs"
	@echo "  make install      # Install the package"
	@echo "  make test         # Run tests"

all: install test
	@[ -f .venv/bin/pybind11-stubgen ] || (set -x && .venv/bin/pip install pybind11-stubgen)
	.venv/bin/pybind11-stubgen image3kit --output-dir src

install:
	pip install .

test:
	@[ -f .venv/bin/pytest ] || (set -x && .venv/bin/pip install pytest)
	.venv/bin/pytest

.PHONY: setup-ide clean build tests


setup-ide:
	python3 -m venv .venv
	./.venv/bin/pip install cmake pybind11 pre-commit
	cmake -S . -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -Dpybind11_DIR=$(shell .venv/bin/python -m pybind11 --cmakedir)
	ln -sf build/compile_commands.json .

clean:
	rm -rf build compile_commands.json

pre-commit:
	pre-commit run -a