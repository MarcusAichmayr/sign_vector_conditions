.PHONY: install

install:
	sage -pip install --upgrade --no-index .

test:
	sage -t sign_vector_conditions/

	sage tests/test_uniqueness.py -v # TODO remove this later
	sage tests/test_robustness.py -v # TODO remove this later

doc:
	cd docs && make html

doc-pdf:
	cd docs && make latexpdf
