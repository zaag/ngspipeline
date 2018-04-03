#!/usr/bin/env python

from distutils.core import setup

setup(name='NGSPIPE',
      version='1.0',
      description='Pipeline for NGS data-analysis',
      author='Martin Haagmans',
      author_email='mahaagmans@gmail.com',
      packages=['diagnostiek'],
      url='www.github.com/zaag/ngspipeline',
      license='MIT',
      scripts=['NGSPIPE'],
      # install_requires=['diagnostiek'],
      )
