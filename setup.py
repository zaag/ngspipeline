#!/usr/bin/env python

from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()


setup(name='NGSPIPE',
      version='1.0',
      description='Pipeline for NGS data-analysis',
      long_description=readme(),
      long_description_content_type='text/markdown',
      author='Martin Haagmans',
      author_email='mahaagmans@gmail.com',
      url='https://www.github.com/zaag/ngspipeline',
      license='MIT',
      scripts=['NGSPIPE'],
      install_requires=['ngsscriptlibrary>=1.0', 
                        'pycnv>=1.0',
                        'openpyxl>=2.4.10']
      )
