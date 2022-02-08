***************
Getting Started
***************

Quickstart
==========

There are two main steps to |delvered|, the `exposure` level processing and the `brick` level forced-photometry processing.  Most of the
software is IDL.

Exposure
--------

The exposure-level processing is grouped on a nightly basis.  The "driver" program is called `delvered_exposures` and can be given a list
or range of nights to process.

.. code-block:: idl

	IDL>delvered_exposures,'20160101'


Brick
-----

The brick-level processing is grouped by `bricks` that cover 0.25x0.25 deg on the sky.  The driver program is called `delvered_bricks` and
can be given a list of brick names.

.. code-block:: idl

	IDL>delvered_bricks,'1887m827'

