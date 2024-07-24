from setuptools import setup


# Get information from separate files (README)
def readfile(filename):
    with open(filename, encoding="utf-8") as f:
        return f.read()


setup(
    name="sign_vector_conditions",
    version="v1.0",
    description="a SageMath package to work with sign vector conditions for chemical reaction networks",
    long_description=readfile("README.md"),
    long_description_content_type="text/markdown",
    url="https://github.com/MarcusAichmayr/sign_vector_conditions",
    author="Marcus S. Aichmayr",
    author_email="aichmayr@mathematik.uni-kassel.de",
    license="GPLv3",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],  # classifiers list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords=[
        "chemical reaction networks",
        "sign vector conditions, oriented matroids",
    ],
    packages=["sign_vector_conditions", "examples"],
    setup_requires=["elementary_vectors>=1.0"],
)
