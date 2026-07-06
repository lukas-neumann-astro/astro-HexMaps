Advanced Configuration
======================

This page documents options that are less commonly changed but useful for
specific use cases.


Resolution Modes
----------------

.. code-block:: ini

   [resolution]
   resolution = angular    # arcseconds (default)
   # resolution = physical # target_res in parsecs, converted via distance
   # resolution = native   # use the overlay beam as-is

For ``resolution = physical``, set ``target_res`` in **parsecs**. HexMaps
converts to arcseconds per source using its distance from
``target_definitions.txt``.


Grid Parameters
---------------

.. code-block:: ini

   [resolution]
   pixels_per_beam = 2       # spacing = target_res / pixels_per_beam
   max_rad         = auto    # "auto" derives radius from overlay footprint
   NAXIS_shuff     = 200     # channels in the shuffled spectrum
   CDELT_SHUFF     = 4000.0  # channel width of shuffled spectrum [m/s]


FOV Edge Erosion
----------------

Pixels near the map boundary are computed from a partial kernel and are
biased. HexMaps trims these by eroding the footprint:

.. code-block:: ini

   [masking]
   fov_erosion_beams = 0.5   # trim 0.5 × beam FWHM (default)
   # fov_erosion_beams = 0   # disable — keep full overlay footprint
   # fov_erosion_beams = 1.0 # conservative — trim one full beam

The same erosion is applied to the hex-grid, moment maps, and FITS outputs.


Reference Line and Mask Combinations
--------------------------------------

The ``ref_line`` key controls both which lines define the primary S/N mask
and how that mask is combined with any external masks.

**Line selection:**

.. code-block:: ini

   ref_line = first         # first cube (default)
   ref_line = 12co21        # specific named line
   ref_line = all           # OR-combine all lines
   ref_line = 2             # first 2 lines
   ref_line = individual    # one mask per line, applied independently

**Combination tokens** (append to any line selection, comma-separated):

.. code-block:: ini

   AND(input)   # AND with external input mask (use_input_mask must be true)
   OR(input)    # OR  with external input mask
   AND(fixed)   # AND with fixed velocity-window mask (use_fixed_vel_mask = true)
   OR(fixed)    # OR  with fixed velocity-window mask

Examples:

.. code-block:: ini

   ref_line = 12co21, AND(input)
   ref_line = first, OR(fixed)
   ref_line = all, AND(input), AND(fixed)

The external masks are defined in the ``# ---- mask ----`` table and
enabled by ``use_input_mask`` / ``use_fixed_vel_mask``.

**Two-level S/N mask:**

.. code-block:: ini

   SN_processing   = 2, 4   # [low_SN, high_SN] thresholds
   strict_mask     = false  # if true, remove spatially isolated detections
   conseq_channels = 3      # min consecutive channels for a valid mask signal


Velocity Windows
----------------

Explicit velocity windows for signal integration and noise estimation:

.. code-block:: ini

   # ---- mask ----
   vel_mask,   Signal window,  -200,  200,   km/s
   noise_mask, Noise blue,    -300, -150,   km/s
   noise_mask, Noise red,      150,  300,   km/s

Enable with:

.. code-block:: ini

   [masking]
   use_fixed_vel_mask   = true
   use_fixed_noise_mask = true


Spectral Smoothing
------------------

.. code-block:: ini

   [spectral]
   spec_smooth        = default   # no smoothing
   # spec_smooth      = overlay   # smooth to overlay spectral resolution
   # spec_smooth      = 5.0       # convolve to 5.0 km/s

   spec_smooth_method = binned    # recommended
   # spec_smooth_method = gauss       # ±10-15% RMS bias; avoid for science
   # spec_smooth_method = combined    # bin first, then Gaussian residual


Hyperfine Structure Correction
--------------------------------

For lines with hyperfine structure (HCN, N₂H⁺, CN, CCH), HexMaps can
extend the signal mask to satellite components:

.. code-block:: ini

   [paths]
   hfs_file = keys/hfs_lines.txt

   [masking]
   use_hfs_lines = true

Add entries for each line with hyperfine structure to ``hfs_lines.txt``.


Database Fill Mode
------------------

To add new maps or cubes to an existing ``.ecsv`` without re-running the
full pipeline:

.. code-block:: ini

   [structure]
   structure_creation = fill
   fname_fill = ngc5194_hexmaps_27p0as_2025_01_01.ecsv

HexMaps opens the existing file and adds only the maps/cubes that are not
yet present.
