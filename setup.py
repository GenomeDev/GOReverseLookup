from setuptools import find_packages, setup
import os

# Get the absolute path to the Readme.md file
readme_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "Readme.md")
with open(readme_path, "r") as f:
    long_description = f.read()

setup(
    name="GOReverseLookup",
    version="1.0.4",
    description="Python library for Gene Ontology Reverse Lookup",
    package_dir={"": "app"},
    packages=find_packages(where="app"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MediWizards/GOReverseLookup",
    author="Aljoša Škorjanc, Vladimir Smrkolj",
    author_email="skorjanc.aljosa@gmail.com, vladimir.smrkolj@gmail.com",
    license="Apache License 2.0",
    install_requires=[],
    python_requires=">=3.10"
)

