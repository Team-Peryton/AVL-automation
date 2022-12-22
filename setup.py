from setuptools import setup, find_packages
import codecs
import os

VERSION = '0.1'
DESCRIPTION = "Automation tools for vortex lattice method, AVL"

# setup
setup(
    name="avl_automation",
    version=VERSION,
    author="jj-foster",
    url="https://github.com/Team-Peryton/AVL-automation",
    packages=find_packages("avl_automation"),
    install_requires=["numpy","scipy","pandas","matplotlib","tqdm"],
    python_requires='>=3.9'
)