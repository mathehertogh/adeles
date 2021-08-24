# Makefile for Sphinx documentation
#

# Check if SageMath is installed.
ifeq ($(shell which sage >/dev/null 2>&1; echo $$?), 1)
$(error The 'sage' command was not found. Make sure you have SageMath installed (cf. https://www.sagemath.org/).)
endif

.PHONY: all install docs clean

all: install docs

install:
	sage -pip install --upgrade .

docs:
	sage -sh -c "sphinx-build -b html docs/src docs"

clean:
	rm -rf docs/html

