.PHONY: install

install:
	sage -pip install --upgrade --no-index .

test:
	sage -t bijectivity_exponential_maps/

	sage tests/test_conditions_injectivity.py -v # remove this later
	sage tests/test_conditions_bijectivity_robust.py -v # remove this later

