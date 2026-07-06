Configuration Guide
===================

This page walks through ``config.txt`` section by section. Every key has a
sensible default; you only need to set what differs from those defaults.

``config.txt`` combines what used to be three separate files in the old
PyStructure (``master_key.txt``, ``data_key.txt``, ``config_key.txt``) into
a single INI-style file with named sections.


[meta]
------

Metadata stored in the output ``.ecsv`` table header for provenance.

.. code-block:: ini

   [meta]
   user     = Dr. Blocksberg
   comments = Example HexMaps run


[paths]
-------

All file and directory paths. Relative paths are resolved relative to the
location of ``config.txt``.

.. code-block:: ini

   [paths]
   data_dir        = data/
   out_dir         = output/
   folder_savefits = ./saved_fits_files/

   # Defaults to keys/target_definitions.txt next to config.txt
   # geom_file = keys/target_definitions.txt

   # Only needed if use_hfs_lines = true
   # hfs_file  = keys/hfs_lines.txt


[sources] and [overlay]
------------------------

.. code-block:: ini

   [sources]
   sources = ngc5194

   [overlay]
   overlay_file = _12co21.fits

``sources`` is a comma-separated list of source names. Each name is
prepended to the file extensions in the map and cube tables to form full
filenames. ``overlay_file`` defines the 3D spectral cube that sets the
spatial extent and spectral axis of the hexagonal grid.


Map and Cube Tables
-------------------

Maps (2D) and cubes (3D) are defined as comma-separated table rows
immediately after their comment markers inside the ``[sources]`` section.

.. code-block:: ini

   # ---- maps ----
   # name,  description,  unit,  file_extension,  directory,  [uc_extension]
   spire250,  SPIRE 250 um,  MJy/sr,  _spire250_gauss21.fits,  data/,  _spire250_gauss21_unc.fits

   # ---- cubes ----
   # name,  description,  unit,  file_extension,  directory,  [map_ext],  [map_uc_ext]
   12co21,  12CO(2-1),  K,  _12co21.fits,  data/
   12co10,  12CO(1-0),  K,  _12co10.fits,  data/

.. IMPORTANT::

   By default, the **first cube in the list is used as the reference line**
   for mask construction. Put your brightest, highest-SNR line first.


Mask Table
----------

An optional mask section supports three types of entries:

.. code-block:: ini

   # ---- mask ----

   # External FITS mask (sampled onto the hex grid):
   # name,  description,  file_extension,  directory
   # co_mask,  CO signal mask,  _co_mask.fits,  data/

   # Fixed velocity window (signal integration range):
   # vel_mask,  description,  v_start,  v_end,  unit
   # vel_mask,  Signal window,  -200,  200,  km/s

   # Noise velocity windows (for RMS estimation):
   # noise_mask,  description,  v_start,  v_end,  unit
   # noise_mask,  Noise blue,  -300,  -150,  km/s
   # noise_mask,  Noise red,    150,   300,  km/s


[resolution]
------------

.. code-block:: ini

   [resolution]
   target_res      = 27.0      # arcsec (for resolution = angular)
   resolution      = angular   # angular | physical | native
   pixels_per_beam = 2         # hex grid spacing in units of beam FWHM
   max_rad         = auto      # map radius in degrees, or "auto"
   NAXIS_shuff     = 200       # channels in shuffled spectrum
   CDELT_SHUFF     = 4000.0    # channel width of shuffled spectrum [m/s]

``resolution`` controls how ``target_res`` is interpreted:

* ``angular`` — use ``target_res`` directly in arcseconds
* ``physical`` — convert ``target_res`` (in pc) to arcseconds using the
  source distance from ``target_definitions.txt``
* ``native`` — use the native beam of the overlay cube


.. _step2:

[masking]
---------

.. code-block:: ini

   [masking]
   ref_line             = first
   SN_processing        = 2, 4
   strict_mask          = false
   use_input_mask       = false
   use_fixed_vel_mask   = false
   use_fixed_noise_mask = false
   use_hfs_lines        = false
   fov_erosion_beams    = 0.5
   mom_thresh           = 5
   conseq_channels      = 3
   mom2_method          = fwhm

``ref_line`` selects which lines build the signal mask:

* ``first`` — first cube in the cube list (default)
* ``<LINE_NAME>`` — a specific named line, e.g. ``12co21``
* ``all`` — OR-combine all cubes
* ``n`` — OR-combine the first *n* cubes
* ``individual`` — build and apply one mask per line independently

**Combination tokens** can be appended to any keyword above:

* ``AND(input)`` — AND-combine with the external input mask
* ``OR(input)``  — OR-combine with the external input mask
* ``AND(fixed)`` — AND-combine with the fixed velocity-window mask
* ``OR(fixed)``  — OR-combine with the fixed velocity-window mask

Examples:

.. code-block:: ini

   ref_line = 12co21, AND(input)       # 12co21 mask AND external mask
   ref_line = first, OR(fixed)         # first-line mask OR fixed window
   ref_line = all, AND(input)          # union of all lines, limited by input


[spectral]
----------

.. code-block:: ini

   [spectral]
   spec_smooth        = default   # default | overlay | <float km/s>
   spec_smooth_method = binned    # binned | gauss | combined

.. WARNING::

   Gaussian spectral smoothing (``spec_smooth_method = gauss``) can
   underestimate RMS by 10–15%. Use ``binned`` or ``combined`` instead.


[output]
--------

.. code-block:: ini

   [output]
   save_cubes    = false   # save convolved PPV cubes (fits stage)
   save_mom_maps = true    # save moment map FITS files
   save_maps     = true    # save 2D band map FITS files
   save_mask     = false   # save the PPV mask as a FITS cube


[structure]
-----------

.. code-block:: ini

   [structure]
   structure_creation = default   # default | fill | archive
   # fname_fill = ngc5194_hexmaps_27p0as_2025_01_01.ecsv

* ``default`` — create or overwrite the ``.ecsv`` file each run
* ``fill`` — open an existing file and add only missing maps/cubes
* ``archive`` — create a new versioned file each run


.. _geomFile:

target_definitions.txt
-----------------------

The ``keys/target_definitions.txt`` file lists source geometry for all
targets that may ever be processed. Add sources here once; only those
listed in ``config.txt [sources]`` will be processed on any given run.

.. code-block:: text

   # source, ra_ctr, dec_ctr, dist_mpc, e_dist_mpc,
   #         incl_deg, e_incl_deg, posang_deg, e_posang_deg, r25, e_r25
   ngc5194, 202.4696, 47.1952, 8.58, 0.10, 22.0, 3.0, 173.0, 3.0, 3.54, 0.05

Galaxy geometry columns (``incl_deg``, ``posang_deg``, ``r25``) are
optional. Leave them blank or omit them for non-galaxy targets such as
Milky Way molecular clouds — the pipeline will skip deprojection and
print a warning.
