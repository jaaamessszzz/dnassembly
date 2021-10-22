try:
    from setuptools import setup, find_packages
    import os.path
except ImportError:
    from distutils.core import setup, find_packages

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
PKG_DIR = "dnassembly"

version = {}

with open("README.md", "r") as f:
    long_description = f.read()

with open("requirements.txt", "r") as f:
    REQUIRED = f.readlines()

with open(os.path.join(BASE_DIR, PKG_DIR, "version.py")) as version_file:
    exec(version_file.read(), version)

# Setup
setup(
    name='dnassembly',
    version=version['__version__'],
    author='James Lucas',
    author_email='james.lucas@berkeley.edu',
    description='Digest and assemble DNA',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/outpace-bio/dnassembly',
    keywords=[
        'DNA',
        'assembly',
        'cloning',
        'plasmid',
        'digest',
        'restriction',
        'enzyme',
        'molecular',
        'biology'
    ],
    classifiers=[
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
    packages=find_packages(),
    install_requires=REQUIRED,
    python_requires=">=3.7",
    include_package_data=True,
    zip_safe=False,
)