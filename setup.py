from setuptools import setup, find_packages

setup(
    name="complax", 
    version="1.0.0",
    author="Federica Lauria",
    author_email="federica.lauria95@gmail.com",
    description="A Python tool for automatic microsolvation and geometry optimization using xTB.",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Fedelau/complax",
    license="MIT",
    keywords=['compchem', 'microsolvation', 'xtb', 'geometry optimization'],
    packages=find_packages(),
    install_requires=[
        "numpy",
        "ase",
        "colorama",
        "tqdm",
        "tabulate"
    ],
    python_requires=">=3.8",
    classifiers=[
        "Programming Language :: Python3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={
    "console_scripts": [
        "complax=complax.complax:main",
    ],
    },
    
)