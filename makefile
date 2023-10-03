.PHONY: install

install:
	sage -pip install --upgrade --no-index .

test:
	sage -t sign_vector_conditions/

doc:
	cd docs && make html

doc-pdf:
	cd docs && make latexpdf
