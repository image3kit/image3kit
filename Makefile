
build_pyi:
	pip install .
	pybind11-stubgen voxlib --output-dir src

install:
	pip install .

test:
	# .venv/bin/pip install pytest numpy
	python -m pytest

.PHONY: setup-ide clean build tests


setup-ide:
	python3 -m venv .venv
	./.venv/bin/pip install cmake pybind11 pre-commit
	cmake -S . -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -Dpybind11_DIR=$(shell ./.venv/bin/python -m pybind11 --cmakedir)
	ln -sf build/compile_commands.json .

clean:
	rm -rf build compile_commands.json

pre-commit:
	pre-commit run -a