"""
stage_products.py — spectral processing, masking, and moment computation.

This stage reads the .ecsv file written by stage_regrid, processes the
spectra, and writes the result back to the same file.

Processing steps
----------------
1. Determine the reference line for mask construction (first cube by default,
   or the line named in config_key ref_line).
2. Build or load the velocity-integration mask.
3. Optionally combine masks from multiple lines or from an HI map.
4. Optionally apply a spatial connectivity (strict) filter.
5. For each spectral line:
   a. Compute moment maps (mom0/1/2, Tpeak, rms, EW) within the mask.
   b. Compute shuffled spectra (shifted so emission is centred at v=0).
6. Write the enriched table back to disk.

Mask construction (construct_mask)
-----------------------------------
The mask uses a two-level S/N approach to capture both bright cores and faint
line wings:

  high_thresh = SN_processing[1] × per-spectrum MAD
  low_thresh  = SN_processing[0] × per-spectrum MAD

A channel is masked if it exceeds high_thresh AND has at least one adjacent
channel above high_thresh.  This core mask is then dilated outward to include
all adjacent channels above low_thresh (up to 5 dilation passes), followed by
two additional edge-growing passes.  This naturally captures asymmetric line
profiles without requiring a fixed velocity window.

Output column naming convention
---------------------------------
MOM0_<LINE>   : integrated intensity (K km/s or equiv.)
EMOM0_<LINE>  : propagated uncertainty on MOM0
MOM1_<LINE>   : intensity-weighted mean velocity
EMOM1_<LINE>  : uncertainty on MOM1
MOM2_<LINE>   : velocity dispersion (FWHM, sqrt(mom2), or math. def.)
EMOM2_<LINE>  : uncertainty on MOM2
TPEAK_<LINE>  : peak brightness temperature within the mask
RMS_<LINE>    : noise rms outside the mask
EW_<LINE>     : equivalent width (∑ T dv / Tpeak / sqrt(2π))
EEW_<LINE>    : uncertainty on EW
SPEC_SHUFF_<LINE>  : shuffled spectrum (n_pts × n_shuff_chan)
SPEC_MASK          : combined velocity-integration mask (n_pts × n_chan)
SPEC_MASK_<LINE>   : per-line mask (stored if ref_line != "first")
SPEC_VAXIS         : velocity axis in km/s (n_pts × n_chan)
SPEC_VAXIS_SHUFF   : shuffled velocity axis in km/s (n_pts × n_shuff_chan)
"""

import numpy as np
import pandas as pd
from astropy import units as au
from astropy.stats import median_absolute_deviation
from astropy.table import Table, Column

from hexmaps.utils_table import shuffle, get_mom_maps, build_noise_mask, parse_ref_line

from hexmaps.logger import get_logger

LOG = get_logger("Products")


# ============================================================================
# Mask construction
# ============================================================================


def construct_mask(ref_line, this_data, SN_processing):
    """
    Build a two-level S/N velocity-integration mask from *ref_line*.

    The algorithm operates spectrum by spectrum, using per-spectrum noise
    estimated via the median absolute deviation (MAD) of emission-free channels.
    Two threshold levels (low and high S/N) are used to ensure that the mask
    captures faint wings attached to high-S/N cores while rejecting isolated
    noise spikes.

    Parameters
    ----------
    ref_line      : str         — column name suffix, e.g. "12CO21" (without "SPEC_")
    this_data     : Table       — the HexMaps table (must contain SPEC_<ref_line>)
    SN_processing : list[float] — [low_SN_thresh, high_SN_thresh]

    Returns
    -------
    mask_q        : astropy Quantity (n_pts × n_chan) — 0/1 mask
    line_vmean    : astropy Quantity (n_pts,)         — intensity-weighted mean velocity
    line_vaxis    : astropy Quantity (n_chan,)         — velocity axis in km/s

    Notes
    -----
    The per-spectrum MAD is estimated using only channels below 3× the global
    MAD (a two-pass approach to avoid contamination from strong emission).
    """
    ref_line_data = this_data[f"SPEC_{ref_line.upper()}"]
    n_pts = np.shape(ref_line_data)[0]
    n_chan = np.shape(ref_line_data)[1]

    # Reconstruct velocity axis from table metadata
    line_vaxis = (
        this_data.meta["SPEC_VCHAN0"]
        + (np.arange(n_chan) - (this_data.meta["SPEC_CRPIX"] - 1))
        * this_data.meta["SPEC_DELTAV"]
    )
    line_vaxis = line_vaxis.to(au.km / au.s)

    # Two-pass global MAD to estimate the noise floor
    rms = median_absolute_deviation(ref_line_data, axis=None, ignore_nan=True)
    rms = median_absolute_deviation(
        ref_line_data[np.where(ref_line_data < 3 * rms)], ignore_nan=True
    )

    # Per-spectrum noise: MAD of channels below the global 3-sigma threshold
    mask_rough = ref_line_data < 3 * rms
    masked_cube = np.where(mask_rough, ref_line_data, np.nan)
    med_mask = np.nanmedian(masked_cube, axis=1)
    mad_mask = np.nanmedian(np.abs(masked_cube - med_mask[:, None]), axis=1)

    low_thresh = SN_processing[0] * mad_mask[:, None]
    high_thresh = SN_processing[1] * mad_mask[:, None]

    # Initial high-S/N mask: channel above high_thresh with adjacent support
    mask = (ref_line_data > high_thresh).astype(int)
    low_mask = (ref_line_data > low_thresh).astype(int)
    mask = mask & (np.roll(mask, 1, 1) | np.roll(mask, -1, 1))

    # Require ≥3 of 3 consecutive channels to suppress single-channel spikes
    mask = ((mask + np.roll(mask, 1, 1) + np.roll(mask, -1, 1)) >= 3).astype(int)
    low_mask = (
        (low_mask + np.roll(low_mask, 1, 1) + np.roll(low_mask, -1, 1)) >= 3
    ).astype(int)

    # Dilate high-S/N core into low-S/N wings (5 passes)
    for _ in range(5):
        mask = ((mask + np.roll(mask, 1, 1) + np.roll(mask, -1, 1)) >= 1).astype(
            int
        ) * low_mask

    # Grow mask edge by 2 channels to ensure full line coverage
    for _ in range(2):
        mask = ((mask + np.roll(mask, 1, 1) + np.roll(mask, -1, 1)) >= 1).astype(int)

    # Compute intensity-weighted mean velocity per sampling point
    mask_q = mask * au.dimensionless_unscaled
    line_vmean = np.zeros(n_pts) * np.nan * au.km / au.s
    for jj in range(n_pts):
        denom = np.nansum(ref_line_data[jj, :] * mask_q[jj, :])
        if denom != 0:
            line_vmean[jj] = (
                np.nansum(line_vaxis * ref_line_data[jj, :] * mask_q[jj, :]) / denom
            )

    return mask_q, line_vmean, line_vaxis


# ============================================================================
# Mask coherence filters  (strict and broad modes)
# ============================================================================


def _coords_and_beam(this_data):
    """Return (ra, dec, beam_deg) from table metadata and coordinate columns."""
    _skip = {"incl_deg", "posang_deg"}
    _coord_cols = [c for c in this_data.colnames
                   if c.lower().endswith("_deg") and c not in _skip]
    if len(_coord_cols) >= 2:
        ra  = np.asarray(this_data[_coord_cols[0]])
        dec = np.asarray(this_data[_coord_cols[1]])
    else:
        ra  = np.asarray(this_data["RA"])
        dec = np.asarray(this_data["DEC"])
    beam_deg = float(this_data.meta["beam_as"].to(au.deg).value)
    return ra, dec, beam_deg


def _apply_strict_mask(mask, this_data):
    """
    Remove spatially isolated mask features (strict mask mode).

    For each spectral channel, label spatially connected groups of masked
    sightlines using pairwise angular distances on the hex grid.  Groups
    smaller than one beam area (estimated as pi*(beam/2)^2 on the hex grid,
    which works out to ~5 sightlines at half-beam spacing) are removed.

    This is the hex-grid analogue of ACES beam-area pruning: any detection
    that does not have at least one beam's worth of spatially connected
    sightlines in the same channel is discarded as a noise spike.

    Parameters
    ----------
    mask      : np.ndarray (n_pts × n_chan) int
    this_data : Table

    Returns
    -------
    mask : np.ndarray — filtered mask (same shape)
    """
    ra, dec, beam_deg = _coords_and_beam(this_data)
    n_chan = mask.shape[1]
    # At half-beam hex spacing, one beam area ~ pi*(0.5)^2 / (sqrt(3)/4) ≈ 3.6
    # sightlines; we round up to 5 to match the ACES ≥3-beam-area criterion.
    half_sep  = beam_deg / 2.0
    tolerance = 0.15 * beam_deg   # ±15% of beam to account for grid irregularity
    min_group = max(5, 1)         # minimum sightlines = ~1 beam area

    mask = mask.copy()
    for jj in range(n_chan):
        mask_spec  = mask[:, jj]
        if not mask_spec.any():
            continue
        labels = np.zeros(len(mask_spec), dtype=int)
        current_label = 1

        for n in range(len(labels)):
            if labels[n] != 0 or mask_spec[n] == 0:
                if mask_spec[n] == 0:
                    labels[n] = -1
                continue
            dist = np.sqrt((ra - ra[n])**2 + (dec - dec[n])**2)
            neigh = np.where(np.abs(dist - half_sep) < tolerance)[0]
            neigh_labels = np.unique(labels[neigh])
            pos_labels   = neigh_labels[neigh_labels > 0]
            if len(pos_labels) > 0:
                labels[n] = pos_labels[0]
                for extra in pos_labels[1:]:
                    labels[labels == extra] = pos_labels[0]
            else:
                labels[n] = current_label
                current_label += 1

        for lab in np.unique(labels):
            if lab <= 0:
                continue
            members = np.where(labels == lab)[0]
            if len(members) < min_group:
                mask[members, jj] = 0

    return mask


def _apply_broad_mask(spec_cube, SN_processing, this_data):
    """
    Build an inclusive (broad) mask using spatial smoothing + two-level S/N
    dilation, following the PHANGS-ALMA broad-masking strategy.

    Steps
    -----
    1. Smooth each channel spatially to a scale of ~2× the beam using a
       Gaussian-weighted average of hex-grid neighbours.  This improves
       sensitivity to faint, spatially coherent emission at the cost of
       spatial resolution.
    2. Estimate per-sightline noise (MAD) on the *smoothed* cube.
    3. Identify a high-S/N core (≥ high_thresh, with consecutive-channel
       support) in the smoothed cube.
    4. Dilate the core into a low-S/N wing mask (≥ low_thresh) in the
       smoothed cube — captures faint line wings attached to bright cores.
    5. Return the dilated mask to be applied to the *original* (unsmoothed)
       spectra when computing moments.

    Parameters
    ----------
    spec_cube     : np.ndarray (n_pts × n_chan) — raw spectra
    SN_processing : list[float] — [low_SN, high_SN] thresholds
    this_data     : Table — for coordinate columns and beam metadata

    Returns
    -------
    mask : np.ndarray (n_pts × n_chan) int — broad integration mask
    """
    ra, dec, beam_deg = _coords_and_beam(this_data)
    n_pts, n_chan = spec_cube.shape

    # ------------------------------------------------------------------
    # Step 1: smooth spectra spatially to ~2× beam using Gaussian weights
    # ------------------------------------------------------------------
    smooth_sigma_deg = beam_deg  # target sigma = 1 beam FWHM ≈ 2.35 sigma
    smoothed = np.full_like(spec_cube, np.nan)

    for n in range(n_pts):
        dist_sq = (ra - ra[n])**2 + (dec - dec[n])**2
        weights  = np.exp(-0.5 * dist_sq / smooth_sigma_deg**2)
        weights[np.all(np.isnan(spec_cube), axis=1)] = 0.0
        w_sum = weights.sum()
        if w_sum > 0:
            smoothed[n, :] = np.nansum(
                weights[:, None] * spec_cube, axis=0
            ) / w_sum

    # ------------------------------------------------------------------
    # Step 2: per-sightline noise on the smoothed cube
    # ------------------------------------------------------------------
    rms = median_absolute_deviation(smoothed, axis=None, ignore_nan=True)
    rms = median_absolute_deviation(
        smoothed[smoothed < 3 * rms], ignore_nan=True
    )
    mask_rough  = smoothed < 3 * rms
    masked_sm   = np.where(mask_rough, smoothed, np.nan)
    med_sm      = np.nanmedian(masked_sm, axis=1)
    mad_sm      = np.nanmedian(np.abs(masked_sm - med_sm[:, None]), axis=1)

    low_thresh  = SN_processing[0] * mad_sm[:, None]
    high_thresh = SN_processing[1] * mad_sm[:, None]

    # ------------------------------------------------------------------
    # Step 3: core mask on smoothed cube (consecutive-channel requirement)
    # ------------------------------------------------------------------
    core = (smoothed > high_thresh).astype(int)
    core = ((core + np.roll(core, 1, 1) + np.roll(core, -1, 1)) >= 3).astype(int)

    # ------------------------------------------------------------------
    # Step 4: wing mask; dilate core into wing (5 passes)
    # ------------------------------------------------------------------
    wing = (smoothed > low_thresh).astype(int)
    wing = ((wing + np.roll(wing, 1, 1) + np.roll(wing, -1, 1)) >= 3).astype(int)

    mask = core.copy()
    for _ in range(5):
        mask = ((mask + np.roll(mask, 1, 1) + np.roll(mask, -1, 1)) >= 1
                ).astype(int) * wing

    # Grow edge by 2 channels for completeness
    for _ in range(2):
        mask = ((mask + np.roll(mask, 1, 1) + np.roll(mask, -1, 1)) >= 1).astype(int)

    return mask


# ============================================================================
# Hyperfine structure mask
# ============================================================================


def _build_hfs_mask(mask, line_name, hfs_data, this_data):
    """
    Extend the mask to cover hyperfine satellite lines.

    Shifts the existing mask by the velocity offset of each satellite line
    relative to the main component.  The union of all shifted masks forms the
    HFS mask, ensuring that all spectral components of the line are included in
    the integration window.

    Parameters
    ----------
    mask      : np.ndarray (n_pts × n_chan) — existing 0/1 mask
    line_name : str        — name of the line to look up in hfs_data
    hfs_data  : pd.DataFrame — hyperfine structure table from handler_keys
    this_data : Table      — used for SPEC_DELTAV metadata

    Returns
    -------
    mask_hfs : astropy Quantity (n_pts × n_chan) — extended mask, or None if
               line_name is not in the HFS table.
    """
    lines_hfs = list(set(hfs_data["hfs_name"]))
    if line_name not in lines_hfs:
        return None

    idx_cols = hfs_data["hfs_name"] == line_name
    restfreqs = [
        f * au.Unit(str(u))
        for f, u in zip(hfs_data["hfs_ref_freq"][idx_cols], hfs_data["unit"][idx_cols])
    ]
    hfs_freqs = [
        f * au.Unit(str(u))
        for f, u in zip(hfs_data["hfs_freq"][idx_cols], hfs_data["unit"][idx_cols])
    ]

    v_ch = this_data.meta["SPEC_DELTAV"].to(au.km / au.s)
    mask_hfs = np.copy(mask)

    for freq, restfreq in zip(hfs_freqs, restfreqs):
        v_shift = freq.to(au.km / au.s, equivalencies=au.doppler_radio(restfreq))
        shift_ch = int(np.rint(v_shift.value / v_ch.value))

        mask_shift = np.zeros_like(mask, dtype=float)
        if shift_ch > 0:
            mask_shift[:, shift_ch:] = mask[:, :-shift_ch]
        elif shift_ch < 0:
            mask_shift[:, :shift_ch] = mask[:, -shift_ch:]
        else:
            mask_shift = mask.copy()

        mask_hfs[mask_shift == 1] = 1

    return mask_hfs * au.dimensionless_unscaled


# ============================================================================
# Stage entry point
# ============================================================================


def run_products(target, fname, meta, cubes, input_mask, hfs_data,
                 window_mask=None, noise_mask_df=None):
    """
    Process all spectra for *target*: mask, moments, shuffle.

    This is the entry point for the "products" pipeline stage.

    Reads the .ecsv file written by stage_regrid, enriches it with the
    columns listed in the module docstring, and overwrites the file.

    Parameters
    ----------
    target     : str
    fname      : str          — path to the .ecsv file from stage_regrid
    meta       : dict         — from KeyHandler.meta
    cubes      : pd.DataFrame — cube definitions from KeyHandler
    input_mask : pd.DataFrame — mask definition from KeyHandler
    hfs_data   : pd.DataFrame or None — hyperfine data from KeyHandler
    """
    # ------------------------------------------------------------------
    # Unpack settings from meta
    # ------------------------------------------------------------------
    ref_line_method  = meta.get("ref_line", "first")
    SN_processing    = meta.get("SN_processing", [2, 4])
    mask_mode        = str(meta.get("strict_mask", "false")).lower()
    use_hfs_lines    = meta.get("use_hfs_lines", False)
    velocity_window  = meta.get("velocity_window", None)
    mom_calc = [
        meta.get("mom_thresh", 5),
        meta.get("conseq_channels", 3),
        meta.get("mom2_method", "fwhm"),
    ]
    shuff_axis = [meta.get("NAXIS_shuff", 200), meta.get("CDELT_SHUFF", 4000.0)]

    this_data  = Table.read(fname)
    line_names = [str(l) for l in cubes["line_name"]]
    n_lines    = len(line_names)

    # ref_line: the first line in the list (used for vmean / shuffling fallback)
    ref_line = (
        ref_line_method
        if isinstance(ref_line_method, str)
        and ref_line_method in line_names
        else line_names[0]
    )

    n_chan = np.shape(this_data[f"SPEC_{ref_line.upper()}"])[1]
    #   mask_lines  — cube lines for S/N masking
    #   use_input   — whether to include the external input mask
    #   use_window  — whether to include the fixed velocity-window mask
    #   combinator  — "AND" or "OR" (default "OR")
    #
    # All requested masks are collected and combined with the combinator.
    # ------------------------------------------------------------------
    mask_lines, use_individual, use_input, use_window, combinator = parse_ref_line(
        ref_line_method, line_names
    )

    # Helper: fetch the external mask array from the table column
    def _get_ext_col(tag):
        """Fetch and remove a pre-sampled mask column from the table."""
        if tag not in this_data.colnames:
            LOG.warning(f"External mask column {tag} not found in table; skipping.")
            return None
        ext = this_data[tag]
        del this_data[tag]
        return np.asarray(ext.value if hasattr(ext, "value") else ext).astype(int)

    def _input_tag():
        return (
            f'SPEC_{str(input_mask["mask_name"].iloc[0]).upper()}'
            if input_mask is not None and len(input_mask) > 0 else None
        )

    def _window_tag():
        return (
            f'SPEC_{str(window_mask["mask_name"].iloc[0]).upper()}'
            if window_mask is not None and len(window_mask) > 0 else None
        )

    # ---- individual mode ------------------------------------------------
    if use_individual:
        LOG.info(f"Building individual masks for {', '.join(mask_lines)}.")
        ref_line_vmean = None
        for line in mask_lines:
            mask_line, vmean_line, _ = construct_mask(
                line, this_data, SN_processing
            )
            this_data[f"SPEC_MASK_{line.upper()}"] = Column(
                mask_line,
                unit=au.dimensionless_unscaled,
                description=f"Velocity-integration mask for {line}",
            )
            if ref_line_vmean is None:
                ref_line_vmean = vmean_line

        # Apply external masks per-line if requested.
        # When use_hfs_lines is True and a line has HFS satellite entries,
        # the external mask is first duplicated and shifted to each satellite
        # frequency before combining, so the external mask covers the same
        # spectral extent as the per-line S/N mask.
        if use_input:
            ext_tag = _input_tag()
            if ext_tag is None:
                LOG.error("ref_line contains 'input' but no input_mask is defined.")
            else:
                ext_arr = _get_ext_col(ext_tag)
                if ext_arr is not None:
                    for line in mask_lines:
                        ext_line = ext_arr
                        if use_hfs_lines and hfs_data is not None:
                            ext_hfs = _build_hfs_mask(
                                ext_arr.astype(float), line, hfs_data, this_data
                            )
                            if ext_hfs is not None:
                                ext_line = np.asarray(
                                    ext_hfs.value if hasattr(ext_hfs, "value") else ext_hfs
                                ).astype(int)
                                LOG.info(f"External input mask shifted to HFS frequencies for {line}.")
                        mask_line = this_data[f"SPEC_MASK_{line.upper()}"].astype(int)
                        mask_line = (
                            (mask_line & ext_line) if combinator == "AND" else (mask_line | ext_line)
                        ) * au.dimensionless_unscaled
                        this_data[f"SPEC_MASK_{line.upper()}"] = Column(
                            mask_line,
                            unit=au.dimensionless_unscaled,
                            description=f"Velocity-integration mask for {line}",
                        )
                    LOG.info(f"Individual masks {combinator} input mask.")
        if use_window:
            ext_tag = _window_tag()
            if ext_tag is None:
                LOG.error("ref_line contains 'window' but no window_mask is defined.")
            else:
                ext_arr = _get_ext_col(ext_tag)
                if ext_arr is not None:
                    for line in mask_lines:
                        ext_line = ext_arr
                        if use_hfs_lines and hfs_data is not None:
                            ext_hfs = _build_hfs_mask(
                                ext_arr.astype(float), line, hfs_data, this_data
                            )
                            if ext_hfs is not None:
                                ext_line = np.asarray(
                                    ext_hfs.value if hasattr(ext_hfs, "value") else ext_hfs
                                ).astype(int)
                                LOG.info(f"External window mask shifted to HFS frequencies for {line}.")
                        mask_line = this_data[f"SPEC_MASK_{line.upper()}"].astype(int)
                        mask_line = (
                            (mask_line & ext_line) if combinator == "AND" else (mask_line | ext_line)
                        ) * au.dimensionless_unscaled
                        this_data[f"SPEC_MASK_{line.upper()}"] = Column(
                            mask_line,
                            unit=au.dimensionless_unscaled,
                            description=f"Velocity-integration mask for {line}",
                        )
                    LOG.info(f"Individual masks {combinator} velocity-window mask.")


    # ---- combined mask mode ---------------------------------------------
    else:
        mask_parts = []

        # S/N masks from cube lines
        ref_line_vmean = None
        if mask_lines:
            LOG.info(f"Building velocity mask from: {', '.join(mask_lines)}.")
            for line in mask_lines:
                mask_line, vmean_line, _ = construct_mask(
                    line, this_data, SN_processing
                )
                this_data[f"SPEC_MASK_{line.upper()}"] = Column(
                    mask_line,
                    unit=au.dimensionless_unscaled,
                    description=f"Velocity-integration mask for {line}",
                )
                if ref_line_vmean is None:
                    ref_line_vmean = vmean_line
                mask_parts.append(
                    np.asarray(
                        mask_line.value
                        if hasattr(mask_line, "value") else mask_line
                    ).astype(int)
                )

        # External FITS input mask (pre-sampled onto hex grid by stage_regrid).
        # When use_hfs_lines is active, the external mask is first extended to
        # cover the HFS satellite frequencies of every HFS-capable line, then
        # combined into mask_parts so the master mask covers all components.
        if use_input:
            ext_tag = _input_tag()
            if ext_tag is None:
                LOG.error("ref_line contains 'input' but no input_mask is defined.")
            else:
                ext_arr = _get_ext_col(ext_tag)
                if ext_arr is not None:
                    mask_parts.append(ext_arr)
                    LOG.info("Input mask included.")

        # Fixed velocity-window mask (pre-sampled onto hex grid by stage_regrid).
        # Same HFS extension logic as for the input mask above.
        if use_window:
            ext_tag = _window_tag()
            if ext_tag is None:
                LOG.error("ref_line contains 'window' but no window_mask is defined.")
            else:
                ext_arr = _get_ext_col(ext_tag)
                if ext_arr is not None:
                    mask_parts.append(ext_arr)
                    LOG.info("Velocity-window mask included.")
                    # Use ref_line vmean for shuffling when only window provided
                    if ref_line_vmean is None:
                        _, ref_line_vmean, _ = construct_mask(
                            line_names[0], this_data, SN_processing
                        )

        # Combine all parts with the combinator
        if not mask_parts:
            LOG.warning("No mask parts resolved; using empty mask.")
            n_pts, n_chan = np.shape(this_data[f"SPEC_{ref_line.upper()}"])
            mask = np.zeros((n_pts, n_chan), dtype=int) * au.dimensionless_unscaled
        else:
            combined = mask_parts[0].copy()
            for part in mask_parts[1:]:
                combined = (combined & part) if combinator == "AND" else (combined | part)
            mask = combined * au.dimensionless_unscaled
            if len(mask_parts) > 1:
                LOG.info(
                    f"Combined {len(mask_parts)} mask(s) with {combinator}."
                )

        # Apply spatial coherence filter according to mask_mode:
        #   'strict' — remove features smaller than ~1 beam area per channel
        #   'broad'  — re-derive mask from spatially smoothed cube
        #   'false'  — no additional filtering
        if mask_mode == "strict":
            LOG.info("Applying strict spatial mask filter (beam-area pruning).")
            mask = _apply_strict_mask(
                np.asarray(mask.value if hasattr(mask, "value") else mask).astype(int),
                this_data,
            ) * au.dimensionless_unscaled
        elif mask_mode == "broad":
            LOG.info("Applying broad spatial mask (smoothed cube + two-level S/N).")
            # Recompute from the primary S/N reference line(s) on the
            # spatially smoothed cube; external masks are re-combined below.
            if mask_lines:
                ref_name = mask_lines[0]
                spec_data = np.array(this_data[f"SPEC_{ref_name.upper()}"])
                broad = _apply_broad_mask(spec_data, SN_processing, this_data)
                for extra_line in mask_lines[1:]:
                    spec_ex = np.array(this_data[f"SPEC_{extra_line.upper()}"])
                    broad_ex = _apply_broad_mask(spec_ex, SN_processing, this_data)
                    broad = broad | broad_ex
                mask = broad * au.dimensionless_unscaled
            else:
                LOG.warning(
                    "broad mask mode requested but no S/N lines resolved; "
                    "falling back to current mask."
                )

        # Store the combined mask
        this_data["SPEC_MASK"] = Column(
            mask,
            unit=au.dimensionless_unscaled,
            description="Velocity-integration mask",
        )

        # HFS mask extension
        # if not use_input and not use_window:
        lines_hfs = (
            list(set(hfs_data["hfs_name"]))
            if (use_hfs_lines and hfs_data is not None)
            else []
        )
        if use_hfs_lines and hfs_data is not None:
            for jj in range(n_lines):
                if line_names[jj] in lines_hfs:
                    LOG.info(f"Building HFS mask for {line_names[jj]}.")
                    mask_hfs = _build_hfs_mask(
                        mask.value, line_names[jj], hfs_data, this_data
                    )
                    if mask_hfs is not None:
                        this_data[f"SPEC_MASK_{line_names[jj].upper()}"] = Column(
                            mask_hfs,
                            unit=au.dimensionless_unscaled,
                            description=f"HFS mask for {line_names[jj].upper()}",
                        )

    LOG.info(f"Mask(s) complete. Computing moments.")

    # ------------------------------------------------------------------
    # Velocity axis columns
    #
    # SPEC_VAXIS and SPEC_VAXIS_SHUFF are written ONCE here, before the loop
    # over lines, because:
    #   (a) they are identical for all lines (both derived from the overlay header)
    #   (b) Astropy Table raises ValueError if you try to overwrite an existing
    #       column by direct assignment — doing this inside the loop causes the
    #       second line onward to fail silently.
    # ------------------------------------------------------------------
    cdelt = shuff_axis[1] * au.m / au.s
    naxis_shuff = int(shuff_axis[0])
    new_vaxis = (cdelt * (np.arange(naxis_shuff) - naxis_shuff / 2)).to(au.km / au.s)

    n_pts_total = len(this_data)
    _v0 = this_data.meta["SPEC_VCHAN0"]
    _dv = this_data.meta["SPEC_DELTAV"]
    _crpix = this_data.meta["SPEC_CRPIX"]
    _n_chan = np.shape(this_data["SPEC_" + line_names[0].upper()])[1]
    _vaxis = (_v0 + (np.arange(_n_chan) - (_crpix - 1)) * _dv).to(au.km / au.s)

    this_data["SPEC_VAXIS"] = Column(
        np.array([_vaxis] * n_pts_total),
        unit=au.km / au.s,
        description="Velocity axis (km/s)",
    )
    this_data["SPEC_VAXIS_SHUFF"] = Column(
        np.array([new_vaxis] * n_pts_total),
        unit=au.km / au.s,
        description="Shuffled velocity axis (km/s)",
    )
    this_data.meta["SPEC_VCHAN0_SHUFF"] = new_vaxis[0]
    this_data.meta["SPEC_DELTAV_SHUFF"] = new_vaxis[1] - new_vaxis[0]

    # ------------------------------------------------------------------
    # Build the noise channel mask (hex-grid path).
    # If use_fixed_noise_mask is True and noise_mask_df is non-empty,
    # build a (n_pts, n_chan) boolean array selecting the channels to use
    # for noise (RMS) estimation. Passed to get_mom_maps as noise_mask.
    # ------------------------------------------------------------------
    use_fixed_noise_mask = meta.get("use_fixed_noise_mask", False)
    hex_noise_mask = None
    if use_fixed_noise_mask:
        if noise_mask_df is not None and len(noise_mask_df) > 0:
            hex_noise_mask = build_noise_mask(
                noise_mask_df, _vaxis, (n_pts_total, n_chan)
            )
            if hex_noise_mask is not None:
                LOG.info(
                    f"Noise RMS will be estimated from {len(noise_mask_df)} "
                    "fixed velocity window(s)."
                )
        else:
            LOG.warning(
                "use_fixed_noise_mask is True but no noise_mask rows found "
                "in the [mask] table. Falling back to mask-inverted noise."
            )

    # ------------------------------------------------------------------
    # Loop over spectral lines: compute moments and shuffled spectra
    # ------------------------------------------------------------------
    for jj in range(n_lines):
        line_name = line_names[jj]
        tag_spec = f"SPEC_{line_name.upper()}"

        if tag_spec not in this_data.keys():
            LOG.error(f"{tag_spec} not found in table; skipping {line_name}.")
            continue

        this_spec = this_data[tag_spec]

        if np.nansum(this_spec, axis=None) == 0:
            LOG.error(f"{line_name} spectrum is all zeros; skipping.")
            continue

        dim_sz = np.shape(this_spec)
        n_pts_l = dim_sz[0]
        n_chan_l = dim_sz[1]

        # Reconstruct the velocity axis for this line's channel grid
        this_v0 = this_data.meta["SPEC_VCHAN0"]
        this_deltav = this_data.meta["SPEC_DELTAV"]
        this_crpix = this_data.meta["SPEC_CRPIX"]
        this_vaxis = (
            this_v0 + (np.arange(n_chan_l) - (this_crpix - 1)) * this_deltav
        ).to(au.km / au.s)

        # Choose the appropriate mask for this line:
        # use_individual → each line has its own mask built from that line's data
        # otherwise      → all lines share the combined master mask
        if use_individual:
            active_mask = this_data[f"SPEC_MASK_{line_name.upper()}"].astype(int)
        else:
            if use_hfs_lines and line_name in lines_hfs:
                hfs_tag = f"SPEC_MASK_{line_name.upper()}"
                active_mask = (
                    this_data[hfs_tag] * au.Unit(1) if hfs_tag in this_data.keys() else mask
                )
            else:
                active_mask = mask
        # LN: can think about implementing option to use the line-specific mean velocities
        line_vmean = ref_line_vmean
            
        # Compute moment maps; use the noise channel mask if available,
        # otherwise fall back to inverting the integration mask.
        line_noise_mask = (
            hex_noise_mask[:n_pts_l] if hex_noise_mask is not None else None
        )
        mom_maps = get_mom_maps(
            this_spec, 
            active_mask, 
            this_vaxis, 
            mom_calc,
            noise_mask=line_noise_mask,
        )
        line_desc = str(cubes["line_desc"].iloc[jj])

        # Only derive moments from the cube if no pre-computed 2D map was
        # provided in data_key for this line (indicated by a non-empty map_ext)
        band_ext_val = str(cubes["map_ext"].iloc[jj]).strip()
        if band_ext_val in ("", "nan"):
            this_data["MOM0_" + line_name.upper()] = Column(
                mom_maps["mom0"], description=f"{line_desc} integrated intensity (mom0)"
            )
            this_data["EMOM0_" + line_name.upper()] = Column(
                mom_maps["mom0_err"], description=f"Uncertainty: {line_desc} mom0"
            )
            this_data["TPEAK_" + line_name.upper()] = Column(
                mom_maps["tpeak"],
                description=f"{line_desc} peak brightness temperature",
            )
            this_data["RMS_" + line_name.upper()] = Column(
                mom_maps["rms"], description=f"{line_desc} rms noise"
            )
            this_data["MOM1_" + line_name.upper()] = Column(
                mom_maps["mom1"],
                description=f"{line_desc} intensity-weighted mean velocity (mom1)",
            )
            this_data["EMOM1_" + line_name.upper()] = Column(
                mom_maps["mom1_err"], description=f"Uncertainty: {line_desc} mom1"
            )
            this_data["MOM2_" + line_name.upper()] = Column(
                mom_maps["mom2"],
                description=f"{line_desc} velocity dispersion (mom2; {mom_calc[2]})",
            )
            this_data["EMOM2_" + line_name.upper()] = Column(
                mom_maps["mom2_err"], description=f"Uncertainty: {line_desc} mom2"
            )
            this_data["EW_" + line_name.upper()] = Column(
                mom_maps["ew"], description=f"{line_desc} equivalent width"
            )
            this_data["EEW_" + line_name.upper()] = Column(
                mom_maps["ew_err"], description=f"Uncertainty: {line_desc} EW"
            )
            LOG.info(f"Moments computed for {line_name}.")
        else:
            LOG.info(
                f"Pre-computed 2D map found for {line_name}; skipping moment computation."
            )

        # Compute shuffled spectrum
        shuffled = shuffle(
            spec=this_spec,
            vaxis=this_vaxis,
            zero=line_vmean, #ref_line_vmean,
            new_vaxis=new_vaxis,
            interp=0,  # nearest-neighbour to preserve noise statistics
        )
        this_data[f"SPEC_SHUFF_{line_name.upper()}"] = Column(
            shuffled,
            unit=this_spec.unit,
            description=f"Velocity-shuffled {line_desc} brightness temperature",
        )
        LOG.info(f"Shuffled spectrum computed for {line_name}.")

    # Shuffled spectral axis reference pixel (always channel 1 after shuffling)
    this_data.meta["SPEC_CRPIX_SHUFF"] = 1

    # ------------------------------------------------------------------
    # Write enriched table back to disk
    # ------------------------------------------------------------------
    # Update the embedded pipeline log so it covers both the regrid and
    # products stages when they are run together.
    from hexmaps.logger import logger as _logger
    try:
        this_data.meta["pipeline_log"] = _logger.as_text().replace("\n", "\\n")
    except Exception as e:
        LOG.warning(f"Could not update pipeline log in metadata: {e}")
   
    this_data.write(fname, format="ascii.ecsv", overwrite=True)
    LOG.info(f"Spectra processing complete for {target}.")
    LOG.info(f"Database written to: {fname}")
