.PHONY: install

install:
	sage -pip install --upgrade --no-index .

test:
	sage -t bijectivity_exponential_maps/functions.py
	sage -t bijectivity_exponential_maps/conditions_injectivity.py
	sage -t bijectivity_exponential_maps/conditions_bijectivity.py
	sage -t bijectivity_exponential_maps/conditions_bijectivity_robust.py
	sage -t bijectivity_exponential_maps/utility.py

	sage tests/test_conditions_injectivity.py -v # remove this later
	sage tests/test_conditions_bijectivity_robust.py -v # remove this later

