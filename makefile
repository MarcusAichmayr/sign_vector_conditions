.PHONY: install

install:
	sage -pip install --upgrade --no-index .

test:
	sage -t sign_vector_conditions/

	sage tests/test_conditions_injectivity.py -v # remove this later
	sage tests/test_conditions_bijectivity_robust.py -v # remove this later

doc:
	cd docs && make html

doc-pdf:
	cd docs && make latexpdf
