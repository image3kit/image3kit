
PY=
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

all: setup_venv
	.venv/bin/pip install .
	@$(MAKE) test
	@$(MAKE) stubgen
	@$(MAKE) pre-commit

stubgen:
	@[ -f .venv/bin/pybind11-stubgen ] || (set -x && .venv/bin/pip install pybind11-stubgen)
	.venv/bin/pybind11-stubgen image3kit --output-dir src

install_global:
	pip install .

test:
	@[ -f .venv/bin/pytest ] || (set -x && .venv/bin/pip install pytest)
	.venv/bin/pytest

.PHONY: setup_venv clean build tests


setup_venv:
	python3 -m venv .venv
	.venv/bin/pip install cmake pre-commit
	.venv/bin/cmake -S . -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
	ln -sf build/compile_commands.json ./

clean:
	rm -rf build compile_commands.json

pre-commit:
	.venv/bin/pre-commit run -a || .venv/bin/pre-commit run -a
