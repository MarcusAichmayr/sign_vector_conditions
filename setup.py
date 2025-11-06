from setuptools import setup


# Get information from separate files (README)
def readfile(filename):
    with open(filename, encoding="utf-8") as f:
        return f.read()


setup(
    name="sign_vector_conditions",
    version="2.1",
    description="a SageMath package to work with sign vector conditions for chemical reaction networks",
    long_description=readfile("README.md"),
    long_description_content_type="text/markdown",
    url="https://github.com/MarcusAichmayr/sign_vector_conditions",
    author="Marcus S. Aichmayr",
    author_email="aichmayr@mathematik.uni-kassel.de",
    license="GPL-3.0-or-later",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],  # classifiers list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords=[
        "reaction networks",
        "equilibrium",
        "robustness",
        "generalized mass-action",
        "sign vector conditions",
        "oriented matroids",
    ],
    packages=["sign_vector_conditions", "examples"],
    install_requires=["elementary_vectors>=2.0", "sign_vectors>=1.0", "certlin>=1.0"],
    extras_require={
        "passagemath": [
            "passagemath-symbolics",
            "passagemath-flint",
            "passagemath-graphs",
            "passagemath-repl",
        ],
    },
)
