#\!/usr/bin/env python3
"""
Viral Genomics Pipeline Setup
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="viral-genomics-pipeline",
    version="1.0.0",
    author="Viral Genomics Team",
    author_email="noreply@example.com",
    description="End-to-end viral genomics pipeline: from raw reads to publication-ready visualizations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mihinduk/viral-genomics-pipeline",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.8",
    install_requires=[
        "pandas>=1.3.0",
        "numpy>=1.20.0",
        "matplotlib>=3.5.0",
        "biopython>=1.79",
    ],
    entry_points={
        "console_scripts": [
            "viral-pipeline=scripts.run_viral_pipeline:main",
        ],
    },
)
