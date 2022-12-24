from setuptools import setup, find_packages
from pathlib import Path

VERSION = '1.0.2'
DESCRIPTION = "Automation tools for vortex lattice method, AVL"
LONG_DESCRIPTION = (Path(__file__).parent / "README.md").read_text()

# setup
setup(
    name="avlautomation",
    version=VERSION,
    author="jj-foster",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/Team-Peryton/AVL-automation",
    packages=find_packages(),
    install_requires=["numpy","scipy","pandas","matplotlib","tqdm"],
    python_requires='>=3.9'
)
