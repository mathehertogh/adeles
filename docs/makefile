# Makefile for Sphinx documentation
#

# Check if SageMath is installed
ifeq ($(shell which sage >/dev/null 2>&1; echo $$?), 1)
$(error The 'sage' command was not found. Make sure you have SageMath installed (cf. https://www.sagemath.org/).)
endif

# Check if Sphinx is installed
ifeq ($(shell which sphinx-build >/dev/null 2>&1; echo $$?), 1)
$(error The 'sphinx-build' command was not found. Make sure you have Sphinx installed (cf. https://www.sphinx-doc.org/).)
endif

.PHONY: build-html clean

build-html:
	sage -sh -c "sphinx-build -b html src html"

clean:
	rm -rf html

