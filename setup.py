from setuptools import (setup, find_packages)

setup(
    name='MultimediaBNS',
    version='0.0.1',
    author='E. Payne',
    packages=find_packages(),
    scripts=['bin/testing.py'],
    license='LICENSE.txt',
    description='Code for the analysis of BNS merger \
     signals in GW and EM bands',
    long_description=open('README.md').read(),
)
