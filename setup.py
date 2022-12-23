from setuptools import setup, find_packages

VERSION = '0.2'
DESCRIPTION = "Automation tools for vortex lattice method, AVL"
LONG_DESCRIPTION = "Automatic tail sizing and dihedral angle investigation made possible through lose wrapping of AVL."

# setup
setup(
    name="avl_automation",
    version=VERSION,
    author="jj-foster",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    url="https://github.com/Team-Peryton/AVL-automation",
    packages=find_packages("avl_automation"),
    install_requires=["numpy","scipy","pandas","matplotlib","tqdm"],
    python_requires='>=3.9'
)
