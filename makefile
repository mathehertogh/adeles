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
	rm -f docs/*.html
	rm -f docs/searchindex.js
	rm -f docs/objects.inv
	rm -f docs/.buildinfo
	rm -rf docs/_static
	rm -rf docs/_images
	rm -rf docs/_sources
	rm -rf docs/.doctrees

