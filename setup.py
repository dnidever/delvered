#!/usr/bin/env python

from distutils.core import setup

setup(name='delvered',
      version='1.0',
      description='David Nidever DELVE reduction code'
      author='David Nidever',
      author_email='dnidever@montana.edu',
      url='https://github.com/dnidever/delvered',
      packages=['delvered'],
      scripts=['bin/make_delvered_table'],
      requires=['numpy','astropy','scipy','dlnpyutils']
)
