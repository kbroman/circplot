all: doc
.PHONY: doc

# build package documentation
doc:
	R -e 'devtools::document()'

# run tests
test:
	R -e 'devtools::test()'
