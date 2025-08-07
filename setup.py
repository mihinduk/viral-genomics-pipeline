from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="viral-genomics-pipeline",
    version="2.0.0",
    author="Kathie Mihindu",
    author_email="mihindu@wustl.edu",
    description="A modular pipeline for viral genomics analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mihinduk/viral-genomics-pipeline-v2",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "pandas>=1.3.0",
        "numpy>=1.21.0",
        "matplotlib>=3.4.0",
        "pyyaml>=6.0",
        "biopython>=1.79",
        "plotly>=5.0.0",
    ],
    entry_points={
        "console_scripts": [
            "viral-pipeline=viral_pipeline.cli:main",
        ],
    },
)
