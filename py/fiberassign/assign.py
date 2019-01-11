# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.assign
=====================

Functions and classes for computing and writing fiber assignment.

"""
from __future__ import absolute_import, division, print_function

import os

import re

import numpy as np

import multiprocessing as mp
from multiprocessing.sharedctypes import RawArray

from functools import partial

from collections import OrderedDict

import fitsio

import desimodel.focalplane

from ._version import __version__

from .utils import Logger, Timer, default_mp_proc

from .tiles import load_tiles

from .targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                      TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE,
                      TargetTree, TargetsAvailable, FibersAvailable,
                      load_target_file, desi_target_type)

from .hardware import (load_hardware, FIBER_STATE_STUCK, FIBER_STATE_BROKEN)

from ._internal import Assignment


# The columns and types that are written to FITS format.  The raw data has
# 3 HDUs, and these are designed to efficiently contain all information used
# internally to fiber assignment without duplicating information:
#
# FASSIGN
# - Per-fiber information, including assignment.  Sorted by fiber ID.
#
# FTARGETS
# - All available targets and their properties relevant to assignment.  Sorted
#   by target ID.
#
# FAVAIL
# - Target IDs available to each fiber.  Sorted by fiber ID and then target ID.
#

results_assign_columns = OrderedDict([
    ("FIBER", "i4"),
    ("TARGETID", "i8"),
    ("LOCATION", "i4"),
    ("FIBERSTATUS", "i4"),
    ("LAMBDA_REF", "f4"),
    ("PETAL_LOC", "i2"),
    ("DEVICE_LOC", "i4"),
    ("DEVICE_TYPE", "a3"),
    ("TARGET_RA", "f8"),
    ("TARGET_DEC", "f8"),
    ("DESIGN_X", "f4"),
    ("DESIGN_Y", "f4"),
    ("DESIGN_Q", "f4"),
    ("DESIGN_S", "f4"),
])

results_targets_columns = OrderedDict([
    ("TARGETID", "i8"),
    ("TARGET_RA", "f8"),
    ("TARGET_DEC", "f8"),
    ("FBATYPE", "u1"),
    ("PRIORITY", "i4"),
    ("SUBPRIORITY", "f8"),
    ("OBSCONDITIONS", "i4"),
    ("NUMOBS_MORE", "i4")
])

results_avail_columns = OrderedDict([
    ("FIBER", "i4"),
    ("TARGETID", "i8")
])


def result_tiles(dir=".", prefix="fiberassign_"):
    # Find all the per-tile files and get the tile IDs
    tiles = list()
    for root, dirs, files in os.walk(dir):
        for f in files:
            mat = re.match(r"{}(\d+).fits".format(prefix), f)
            if mat is not None:
                # Matches the prefix
                tiles.append(int(mat.group(1)))
    return tiles


def result_path(tile_id, dir=".", prefix="fiberassign_",
                ext="fits", create=False, split=False):
    tiledir = dir
    if split:
        tilegroup = tile_id // 1000
        tiledir = os.path.join(dir, "{:02d}".format(tilegroup))
    if create:
        os.makedirs(tiledir, exist_ok=True)
    path = os.path.join(tiledir,
                        "{}{:05d}.{}".format(prefix, tile_id, ext))
    return path


def write_assignment_fits_tile(asgn, fulltarget, params):
    """Write a single tile assignment to a FITS file.

    Args:
        outroot (str):  full path of the output root file name.
        asgn (Assignment):  the assignment class instance.
        fulltarget (bool):  if True, dump the target information for all
            available targets, not just the ones that are assigned.
        params (tuple):  tuple containing the tile ID, RA, DEC,
            output path, and GFA targets

    """
    tm = Timer()
    tm.start()
    tile_id, tile_ra, tile_dec, tile_file, gfa_targets = params
    log = Logger.get()

    # Hardware properties
    hw = asgn.hardware()

    # Targets available for all tile / fibers
    tgsavail = asgn.targets_avail()

    # Target properties
    tgs = asgn.targets()

    # Data for this tile
    tdata = asgn.tile_fiber_target(tile_id)
    avail = tgsavail.tile_data(tile_id)

    # The recarray dtypes
    assign_dtype = np.dtype([(x, y) for x, y in
                            results_assign_columns.items()])
    targets_dtype = np.dtype([(x, y) for x, y in
                             results_targets_columns.items()])
    avail_dtype = np.dtype([(x, y) for x, y in
                            results_avail_columns.items()])
    tm.stop()
    tm.report("  data pointers for tile {}".format(tile_id))
    tm.clear()
    tm.start()

    log.debug("Write:  indexing tile {}".format(tile_id))

    fibers = np.array(hw.fiber_id)

    # Compute the total list of targets
    navail = np.sum([len(avail[x]) for x in avail.keys()])
    tgids_avail = np.unique(np.concatenate([np.array(avail[x], dtype=np.int64)
                                            for x in avail.keys()]))

    tm.stop()
    tm.report("  available targets tile {}".format(tile_id))

    if len(tgids_avail) > 0:
        # We have some available targets.
        # tm.clear()
        # tm.start()

        if os.path.isfile(tile_file):
            raise RuntimeError("output file {} already exists"
                               .format(tile_file))
        # This tile has some available targets
        log.info("Writing tile {}".format(tile_id))
        fd = fitsio.FITS(tile_file, "rw")

        tm.stop()
        tm.report("  opening file for tile {}".format(tile_id))
        tm.clear()
        tm.start()

        # Unpack all our target properties from C++ objects into numpy arrays

        tgids = None
        if fulltarget:
            # We are dumping all targets
            tgids = tgids_avail
        else:
            # We are dumping only assigned targets
            tgids = np.array(sorted([y for x, y in tdata.items()]),
                             dtype=np.int64)

        tg_ra = np.empty(len(tgids), dtype=np.float64)
        tg_dec = np.empty(len(tgids), dtype=np.float64)
        tg_type = np.empty(len(tgids), dtype=np.uint8)
        tg_priority = np.empty(len(tgids), dtype=np.int32)
        tg_subpriority = np.empty(len(tgids), dtype=np.float64)
        tg_obscond = np.empty(len(tgids), dtype=np.int32)
        tg_obsrem = np.empty(len(tgids), dtype=np.int32)
        off = 0
        for tg in tgids:
            props = tgs.get(tg)
            tg_obsrem[off] = 0
            if props.obs_remain > 0:
                tg_obsrem[off] = props.obs_remain
            tg_ra[off] = props.ra
            tg_dec[off] = props.dec
            tg_type[off] = props.type
            tg_priority[off] = props.priority
            tg_subpriority[off] = props.subpriority
            tg_obscond[off] = props.obscond
            off += 1
        tg_indx = {y: x for x, y in enumerate(tgids)}

        tm.stop()
        tm.report("  extract target props for {}".format(tile_id))
        tm.clear()
        tm.start()

        header = dict()
        header["TILEID"] = tile_id
        header["TILERA"] = tile_ra
        header["TILEDEC"] = tile_dec
        header["FBAVER"] = __version__

        # FIXME:  write "-1" for unassigned targets.  Write all other fiber
        # status and other properties to this table.

        log.debug("Write:  copying assignment data for tile {}"
                  .format(tile_id))

        fdata = np.empty(len(fibers), dtype=assign_dtype)
        fdata["FIBER"] = fibers

        # For unassigned fibers, we give each fiber a unique negative
        # number based on the tile and fiber.
        unassign_offset = tile_id * len(fibers)
        assigned_tgids = np.array([tdata[x] if x in tdata.keys()
                                  else -(unassign_offset + x)
                                  for x in fibers], dtype=np.int64)
        fdata["TARGETID"] = assigned_tgids

        # Rows containing assigned fibers
        assigned_valid = np.where(assigned_tgids >= 0)[0]
        assigned_invalid = np.where(assigned_tgids < 0)[0]

        # Buffers for X/Y/RA/DEC
        assigned_tgx = np.empty(len(assigned_tgids), dtype=np.float64)
        assigned_tgy = np.empty(len(assigned_tgids), dtype=np.float64)
        assigned_tgra = np.empty(len(assigned_tgids), dtype=np.float64)
        assigned_tgdec = np.empty(len(assigned_tgids), dtype=np.float64)

        if (len(assigned_invalid) > 0):
            # Fill our unassigned fiber X/Y coordinates with the central
            # positioner locations.  Then convert these to RA/DEC.
            empty_fibers = fibers[assigned_invalid]
            fpos_xy_mm = dict(hw.fiber_pos_xy_mm)
            empty_x = np.array(
                [fpos_xy_mm[f][0] for f in empty_fibers], dtype=np.float64)
            empty_y = np.array(
                [fpos_xy_mm[f][1] for f in empty_fibers], dtype=np.float64)
            radec = hw.xy2radec_multi(tile_ra, tile_dec, empty_x, empty_y, 0)
            assigned_tgx[assigned_invalid] = empty_x
            assigned_tgy[assigned_invalid] = empty_y
            assigned_tgra[assigned_invalid] = [x for x, y in radec]
            assigned_tgdec[assigned_invalid] = [y for x, y in radec]

        if (len(assigned_valid) > 0):
            # The target IDs assigned to fibers (note- NOT sorted)
            assigned_real = np.copy(assigned_tgids[assigned_valid])

            # Mapping of our assigned target rows into the target properties.
            target_rows = [tg_indx[x] for x in assigned_real]

            real_ra = np.array(tg_ra[target_rows])
            real_dec = np.array(tg_dec[target_rows])

            # Now fill in our real targets
            xy = hw.radec2xy_multi(tile_ra, tile_dec, real_ra, real_dec, 0)
            assigned_tgra[assigned_valid] = real_ra
            assigned_tgdec[assigned_valid] = real_dec
            assigned_tgx[assigned_valid] = [x for x, y in xy]
            assigned_tgy[assigned_valid] = [y for x, y in xy]

        fdata["TARGET_RA"] = assigned_tgra
        fdata["TARGET_DEC"] = assigned_tgdec
        fdata["DESIGN_X"] = assigned_tgx
        fdata["DESIGN_Y"] = assigned_tgy

        assigned_q, assigned_s = desimodel.focalplane.xy2qs(
            assigned_tgx, assigned_tgy)

        fdata["DESIGN_Q"] = assigned_q
        fdata["DESIGN_S"] = assigned_s

        location = dict(hw.fiber_location)
        device = dict(hw.fiber_device)
        petal = dict(hw.fiber_petal)
        device_type = dict(hw.fiber_device_type)

        fdata["LOCATION"] = np.array(
            [location[x] for x in fibers]).astype(np.int32)
        fdata["DEVICE_LOC"] = np.array(
            [device[x] for x in fibers]).astype(np.int32)
        fdata["PETAL_LOC"] = np.array(
            [petal[x] for x in fibers]).astype(np.int16)
        fdata["DEVICE_TYPE"] = np.array(
            [device_type[x] for x in fibers]).astype(np.dtype("a3"))

        # This hard-coded value copied from the original code...
        lambda_ref = np.ones(len(fibers), dtype=np.float32) * 5400.0
        fdata["LAMBDA_REF"] = lambda_ref

        # Fiber status
        fstate = dict(hw.state)
        fstatus = np.zeros(len(fibers), dtype=np.int32)
        # Set unused bit
        fstatus |= [0 if x in tdata.keys() else 1 for x in fibers]
        # Set stuck / broken bits
        fstatus |= [2 if (fstate[x] & FIBER_STATE_STUCK) else 0
                    for x in fibers]
        fstatus |= [4 if (fstate[x] & FIBER_STATE_BROKEN) else 0
                    for x in fibers]
        fdata["FIBERSTATUS"] = fstatus

        tm.stop()
        tm.report("  copy and compute assign data tile {}".format(tile_id))
        tm.clear()
        tm.start()

        log.debug("Write:  writing assignment data for tile {}"
                  .format(tile_id))

        fd.write(fdata, header=header, extname="FASSIGN")
        del fdata

        tm.stop()
        tm.report("  write assign data tile {}".format(tile_id))
        tm.clear()
        tm.start()

        log.debug("Write:  writing target data for tile {}"
                  .format(tile_id))

        fdata = np.empty(len(tgids), dtype=targets_dtype)
        fdata["TARGETID"] = tgids
        fdata["TARGET_RA"] = tg_ra
        fdata["TARGET_DEC"] = tg_dec
        fdata["FBATYPE"] = tg_type
        fdata["PRIORITY"] = tg_priority
        fdata["SUBPRIORITY"] = tg_subpriority
        fdata["OBSCONDITIONS"] = tg_obscond
        fdata["NUMOBS_MORE"] = tg_obsrem

        tm.stop()
        tm.report("  copy targets data tile {}".format(tile_id))
        tm.clear()
        tm.start()

        fd.write(fdata, header=header, extname="FTARGETS")
        del fdata

        tm.stop()
        tm.report("  write targets data tile {}".format(tile_id))
        tm.clear()
        tm.start()

        log.debug("Write:  writing avail data for tile {}"
                  .format(tile_id))

        fdata = np.empty(navail, dtype=avail_dtype)
        off = 0
        for fid in sorted(avail.keys()):
            for tg in avail[fid]:
                fdata[off] = (fid, tg)
                off += 1

        tm.stop()
        tm.report("  copy avail data tile {}".format(tile_id))
        tm.clear()
        tm.start()

        fd.write(fdata, header=header, extname="FAVAIL")
        del fdata

        if gfa_targets is not None:
            try:
                #- Astropy Table
                fd.write(gfa_targets.as_array(), extname='GFA_TARGETS')
            except AttributeError:
                #- numpy structured array
                fd.write(gfa_targets, extname='GFA_TARGETS')

        fd.close()

        tm.stop()
        tm.report("  write avail data tile {}".format(tile_id))
    return


def write_assignment_fits(tiles, asgn, out_dir=".", out_prefix="fiberassign",
                          split_dir=False, all_targets=False, gfa_targets=None):
    """Write out assignment results in FITS format.

    For each tile, all available targets (not only the assigned targets) and
    their properties are written to the first HDU.  The second HDU contains
    one row for every fiber and target available to that fiber.

    Args:
        tiles (Tiles):  The Tiles object containing the properties of each
            tile.
        asgn (Assignment):  The assignment object.
        out_dir (str):  The output directory for writing per-tile files.
        out_prefix (str):  The output file name prefix.
        split_dir (bool):  Optionally split files by tile ID prefix.
        all_targets (bool):  Optionally dump the target properties of all
            available targets for each tile.  If False, only dump the target
            properties of assigned targets.
        gfa_targets (list of numpy arrays): Include these as GFA_TARGETS HDUs
    """
    tm = Timer()
    tm.start()

    # The tile IDs that were assigned
    tileids = asgn.tiles_assigned()
    tileorder = tiles.order
    tilera = tiles.ra
    tiledec = tiles.dec

    write_tile = partial(write_assignment_fits_tile, asgn, all_targets)

    for i, tid in enumerate(tileids):
        tra = tilera[tileorder[tid]]
        tdec = tiledec[tileorder[tid]]
        outfile = result_path(tid, dir=out_dir, prefix=out_prefix,
                              create=True, split=split_dir)
        if gfa_targets is None:
            params = (tid, tra, tdec, outfile, None)
        else:
            params = (tid, tra, tdec, outfile, gfa_targets[i])

        write_tile(params)

    tm.stop()
    tm.report("Write output files")

    return


def write_assignment_ascii(tiles, asgn, out_dir=".", out_prefix="fiberassign",
                           split_dir=False):
    """Write out assignment results in ASCII format.

    For each tile, only the final assignment to each tile is written out.  For
    the full information, including available targets, use the FITS format.

    Args:
        tiles (Tiles):  The Tiles object containing the properties of each
            tile.
        asgn (Assignment):  The assignment object.
        out_dir (str):  The output directory for writing per-tile files.
        out_prefix (str):  The output file name prefix.
        split_dir (bool):  Optionally split files by tile ID prefix.

    """
    log = Logger.get()
    # Go through the assignment, one tile at a time.  For each tile, get the
    # best assignment and potential targets.

    tileids = asgn.tiles_assigned()
    tileorder = tiles.order
    tilera = tiles.ra
    tiledec = tiles.dec
    # Target properties
    tgs = asgn.targets()

    for t in tileids:
        tdata = asgn.tile_fiber_target(t)
        nfiber = len(tdata)
        tfile = result_path(t, dir=out_dir, prefix=out_prefix,
                            create=True, split=split_dir)
        if nfiber > 0:
            log.debug("Writing tile {}".format(t))
            with open(tfile, "w") as f:
                f.write("# TILE_RA = {}\n".format(tilera[tileorder[t]]))
                f.write("# TILE_DEC = {}\n".format(tiledec[tileorder[t]]))
                f.write("#\n")
                f.write("# FIBER  TARGETID  RA  DEC  PRIORITY  "
                        "SUBPRIORITY  OBSCONDITIONS  NUMOBS_MORE  FBATYPE\n")
                for fid in sorted(tdata.keys()):
                    tgid = tdata[fid]
                    tg = tgs.get(tgid)
                    f.write("{:d} {:d} {:.6f} {:.6f}\n"
                            .format(fid, tgid, tg.ra, tg.dec, tg.priority,
                                    tg.subpriority, tg.obscond, tg.obs_remain,
                                    tg.type))
    return


def avail_table_to_dict(avail_data):
    """Convert a recarray of available targets into a dictionary.

    Args:
        avail_data (array):  The available targets read from a FITS HDU
            (for example).

    Returns:
        (dict):  A dictionary of numpy arrays, one per fiber, containing the
            available target IDs for the fiber.

    """
    avail_target = avail_data["TARGETID"]
    avail_fiber = None
    if "LOCATION" in avail_data.dtype.names:
        avail_fiber = avail_data["LOCATION"]
    else:
        avail_fiber = avail_data["FIBER"]
    avail = dict()
    for fid, tgid in zip(avail_fiber, avail_target):
        if fid in avail:
            avail[fid].append(tgid)
        else:
            avail[fid] = list([tgid])
    avail = {f: np.array(av) for f, av in avail.items()}
    return avail


def read_assignment_fits_tile(params):
    """Read in results.

    This reads in only the result information that was originally written by
    the fiber assignment.  This function is used internally for reading data
    for use in plotting and for reading the original data for merging with
    input target catalogs.

    Args:
        params (tuple):  The tile ID and tile file path packed as a tuple for
            use in multiprocessing.

    Returns:
        (tuple):  The FITS header, assignment recarray, target property
            recarray, the available targets recarray, and the GFA targets

    """
    (tile_id, tile_file) = params
    log = Logger.get()
    # Output tile file
    if not os.path.isfile(tile_file):
        raise RuntimeError("input file {} does not exist".format(tile_file))

    log.debug("Reading tile data {}".format(tile_file))

    # Open the file
    fd = fitsio.FITS(tile_file, "r")

    header = None
    fiber_data = None
    targets_data = None
    avail_data = None
    if "FASSIGN" in fd:
        # We have a new-style file
        header = fd["FASSIGN"].read_header()
        fiber_data = fd["FASSIGN"].read()
        targets_data = fd["FTARGETS"].read()
        avail_data = fd["FAVAIL"].read()
    else:
        # We have old-style data.  Build new recarrays from the FIBERASSIGN
        # and POTENTIAL_ASSIGNMENTS extensions.
        header = fd["FIBERASSIGN"].read_header()
        nfiber = fd["FIBERASSIGN"].get_nrows()
        fbassign = fd["FIBERASSIGN"].read()
        assign_dtype = np.dtype([(x, y) for x, y in
                                results_assign_columns.items()])
        fiber_data = np.empty(nfiber, dtype=assign_dtype)
        fiber_data[:] = [(x["FIBER"], x["TARGETID"]) for x in fbassign]

        tgrows = np.where(fbassign["TARGETID"] >= 0)[0]
        ntarget = np.sum([1 for x in tgrows])
        targets_dtype = np.dtype([(x, y) for x, y in
                                 results_targets_columns.items()])
        targets_data = np.empty(ntarget, dtype=targets_dtype)
        targets_data["TARGETID"] = fbassign["TARGETID"][tgrows]
        targets_data["RA"] = fbassign["TARGET_RA"][tgrows]
        targets_data["DEC"] = fbassign["TARGET_DEC"][tgrows]
        targets_data["PRIORITY"] = fbassign["PRIORITY"][tgrows]
        targets_data["SUBPRIORITY"] = fbassign["SUBPRIORITY"][tgrows]
        targets_data["OBSCONDITIONS"] = fbassign["OBSCONDITIONS"][tgrows]
        targets_data["NUMOBS_MORE"] = fbassign["NUMOBS_MORE"][tgrows]
        targets_data["FBATYPE"] = [desi_target_type(x) for x in
                                   fbassign["DESI_TARGET"][tgrows]]
        del fbassign
        if "POTENTIALTARGETID" in fd:
            log.warning("File {} is an old format with an incompatible "
                        "packing of the potential targets HDU."
                        .format(tile_file))
            avail_data = fd["POTENTIALTARGETID"].read()
        else:
            avail_data = fd["POTENTIAL_ASSIGNMENTS"].read()

    if "GFA_TARGETS" in fd:
        gfa_targets = fd["GFA_TARGETS"].read()
    else:
        gfa_targets = None

    return header, fiber_data, targets_data, avail_data, gfa_targets


merge_results_tile_tgbuffers = None
merge_results_tile_tgdtypes = None
merge_results_tile_tgshapes = None


def merge_results_tile_initialize(bufs, dtypes, shapes):
    global merge_results_tile_tgbuffers
    global merge_results_tile_tgdtypes
    global merge_results_tile_tgshapes
    merge_results_tile_tgbuffers = bufs
    merge_results_tile_tgdtypes = dtypes
    merge_results_tile_tgshapes = shapes
    return


merged_fiberassign_swap = {
    "RA": "TARGET_RA",
    "DEC": "TARGET_DEC",
    "RA_IVAR": "TARGET_RA_IVAR",
    "DEC_IVAR": "TARGET_DEC_IVAR"
}

merged_fiberassign_req_columns = OrderedDict([
    ("FIBER", "i4"),
    ("LOCATION", "i4"),
    ("NUMTARGET", "i2"),
    ("TARGETID", "i8"),
    ("DESI_TARGET", "i8"),
    ("BGS_TARGET", "i8"),
    ("MWS_TARGET", "i8"),
    ("TARGET_RA", "f8"),
    ("TARGET_DEC", "f8"),
    ("DESIGN_X", "f4"),
    ("DESIGN_Y", "f4"),
    ("BRICKNAME", "a8"),
    ("FIBERSTATUS", "i4"),
    ("DESIGN_Q", "f4"),
    ("DESIGN_S", "f4"),
    ("LAMBDA_REF", "f4"),
    ("OBJTYPE", "a3"),
    ("PETAL_LOC", "i2"),
    ("DEVICE_LOC", "i4"),
    ("PRIORITY", "i4"),
    ("SUBPRIORITY", "f8"),
    ("OBSCONDITIONS", "i4"),
    ("NUMOBS_MORE", "i4")
])

merged_skymon_columns = OrderedDict([
    ("FIBER", "i4"),
    ("LOCATION", "i4"),
    ("NUMTARGET", "i2"),
    ("TARGETID", "i8"),
    ("DESI_TARGET", "i8"),
    ("BGS_TARGET", "i8"),
    ("MWS_TARGET", "i8"),
    ("TARGET_RA", "f8"),
    ("TARGET_DEC", "f8"),
    ("DESIGN_X", "f4"),
    ("DESIGN_Y", "f4"),
    ("BRICKNAME", "a8"),
    ("FIBERSTATUS", "i4"),
    ("DESIGN_Q", "f4"),
    ("DESIGN_S", "f4"),
    ("PETAL_LOC", "i2"),
    ("DEVICE_LOC", "i4"),
    ("PRIORITY", "i4"),
])

merged_potential_columns = OrderedDict([
    ("TARGETID", "i8"),
    ("FIBER", "i4"),
    ("LOCATION", "i4")
])


def merge_results_tile(out_dtype, params):
    """Merge results for one tile.

    This uses target catalog data which has been pre-staged to shared memory.
    This function should not be called directly, but only from the
    merge_results() function.

    Args:
        out_dtype (np.dtype):  The output recarray dtype for the merged target
            HDU.  This is the union of columns chosen from the input catalogs.
        params (tuple):  The tile ID and input / output files.  Set by
            multiprocessing call.
    Returns:
        None.

    """
    (tile_id, infile, outfile) = params
    log = Logger.get()

    log.info("Reading raw tile data {}".format(infile))

    tm = Timer()
    tm.start()

    inhead, fiber_data, targets_data, avail_data, gfa_targets = \
            read_assignment_fits_tile((tile_id, infile))

    tm.stop()
    tm.report("  read input data {}".format(tile_id))
    tm.clear()
    tm.start()

    # Get the list of all targets we are considering for this tile.  This
    # will be either just the assigned targets or all available targets
    # depending on how the user wrote the file.

    tile_tgids = np.copy(targets_data["TARGETID"])
    tile_tgindx = {y: x for x, y in enumerate(tile_tgids)}

    # The row indices of our assigned targets.  These indices are valid for
    # both the original input FTARGETS data and also our per-tile copy of the
    # target catalog data.  These indices are essentially random access into
    # the target table.
    target_rows = np.array([tile_tgindx[x] for x in tile_tgids])

    # Extract just these targets from the full set of catalogs

    tile_targets_dtype = np.dtype(out_dtype.fields)
    tile_targets = np.empty(len(tile_tgids), dtype=tile_targets_dtype)

    # Loop over input target files and copy data.  Note: these are guaranteed
    # to be sorted by TARGETID, since that is done during staging to shared
    # memory.
    targetfiles = list(merge_results_tile_tgbuffers.keys())
    for tf in targetfiles:
        tgview = np.frombuffer(merge_results_tile_tgbuffers[tf],
                               dtype=merge_results_tile_tgdtypes[tf])\
                               .reshape(merge_results_tile_tgshapes[tf])
        # Some columns may not exist in all target files (e.g. PRIORITY),
        # So we select the valid columns for this file and only copy those.
        # The ordering of the targets in the catalog is not guaranteed to be
        # sorted, and the output table is sorted by fiber ID, not target.  So
        # must build an explicit mapping from target catalog rows to output
        # table rows.
        tgids = tgview["TARGETID"]
        inrows = np.where(np.isin(tgids, tile_tgids))[0]
        outrows = np.where(np.isin(tile_tgids, tgids))[0]
        tfcolsin = list()
        tfcolsout = list()
        for c in tgview.dtype.names:
            nm = c
            if c in merged_fiberassign_swap:
                nm = merged_fiberassign_swap[c]
            if nm in out_dtype.names:
                tfcolsin.append(c)
                tfcolsout.append(nm)

        tm.stop()
        tm.report("  indexing {} for {}".format(tf, tile_id))
        tm.clear()
        tm.start()

        if len(outrows) > 0:
            for irw, orw in zip(inrows, outrows):
                for c, nm in zip(tfcolsin, tfcolsout):
                    tile_targets[nm][orw] = tgview[c][irw]
        del tgids
        del inrows
        del outrows
        del tgview

        tm.stop()
        tm.report("  copy targets from {} for {}".format(tf, tile_id))
        tm.clear()
        tm.start()

    # Now we have a reduced set of target data including only the targets
    # relevant for this tile, and with data merged from all the files.
    # Next, merge this with the assignment information.

    # Determine the rows of the assignment that are for science and sky
    # monitor positioners.
    science_rows = np.where(fiber_data["DEVICE_TYPE"] == b"POS")[0]
    sky_rows = np.where(fiber_data["DEVICE_TYPE"] == b"ETC")[0]

    # Construct output recarray
    outdata = np.zeros(len(fiber_data), dtype=out_dtype)

    # Build mapping from FTARGETS rows (sorted by target) to FASSIGN rows
    # (sorted by fiber).

    # Rows containing assigned fibers
    fassign_valid = np.where(fiber_data["TARGETID"] >= 0)[0]

    tm.stop()
    tm.report("  fiber / target index mapping {}".format(tile_id))
    tm.clear()
    tm.start()

    # Copy original data and also determine which of our required columns
    # will come from external catalogs
    external_cols = list()
    for field in out_dtype.names:
        if field in results_assign_columns:
            # Copy assignment and fiber property columns directly.
            outdata[field] = fiber_data[field]
        else:
            if field in results_targets_columns:
                # This is a column we are copying from our raw output.
                if (len(target_rows) > 0):
                    outdata[field][fassign_valid] = \
                        targets_data[field][target_rows]
            else:
                # This required column is coming from external catalogs
                external_cols.append(field)

    tm.stop()
    tm.report("  copy raw data to output {}".format(tile_id))
    tm.clear()
    tm.start()

    # Now copy external target properties for the remaining columns.
    if (len(target_rows) > 0):
        # # Looping over rows and then columns is faster, likely because there
        # # are so many columns and the data is stored row-major.
        # for irw, orw in zip(target_rows, fassign_valid):
        #     for c in external_cols:
        #         realname = c
        #         if c in merged_fiberassign_swap:
        #             realname = merged_fiberassign_swap[c]
        #         outdata[realname][orw] = tile_targets[c][irw]
        for c in external_cols:
            outdata[c][fassign_valid] = tile_targets[c][target_rows]

    tm.stop()
    tm.report("  copy external data to output {}".format(tile_id))
    tm.clear()
    tm.start()

    # Create the file
    if os.path.isfile(outfile):
        os.remove(outfile)
    fd = fitsio.FITS(outfile, "rw")

    # Write the main HDU- only the data for science positioners
    log.info("Writing new data {}".format(outfile))
    fd.write(outdata[science_rows], header=inhead, extname="FIBERASSIGN")

    # Now write out the sky monitor fibers.  We extract the rows and columns
    # from the already-computed recarray.

    skymon_dtype = np.dtype([(x, y) for x, y in merged_skymon_columns.items()])
    skymon = np.zeros(len(sky_rows), dtype=skymon_dtype)
    # Sky monitor fake FIBER column with values 0-19.  The fake FIBER value in
    # the raw data is already based on increasing LOCATION value.
    skymon_fiber = np.arange(len(sky_rows), dtype=np.int32)
    for field in skymon_dtype.names:
        if field == "FIBER":
            skymon["FIBER"] = skymon_fiber
        else:
            skymon[field] = outdata[field][sky_rows]
    fd.write(skymon, header=inhead, extname="SKY_MONITOR")

    if gfa_targets is not None:
        fd.write(gfa_targets, extname="GFA_TARGETS")

    # Write the per-tile catalog information also.  Sadly, this HDU is
    # expected to have the original column names.  We swap them back, which
    # is fine since we are going to delete this data after writing anyway.
    backswap = {y: x for x, y in merged_fiberassign_swap.items()}
    curnames = np.copy(tile_targets.dtype.names)
    newnames = list()
    for nm in curnames:
        if nm in backswap:
            newnames.append(backswap[nm])
        else:
            newnames.append(nm)
    tile_targets.dtype.names = newnames
    fd.write(tile_targets, header=inhead, extname="TARGETS")
    del tile_targets

    # Write the "POTENTIAL_ASSIGNMENTS" HDU
    potential_dtype = np.dtype([(x, y) for x, y
                                in merged_potential_columns.items()])
    potential = np.zeros(len(avail_data), dtype=potential_dtype)

    fiberloc = {x: y for x, y in
                zip(fiber_data["FIBER"], fiber_data["LOCATION"])}
    potential["FIBER"] = avail_data["FIBER"]
    potential["TARGETID"] = avail_data["TARGETID"]
    potential["LOCATION"] = [fiberloc[x] for x in avail_data["FIBER"]]
    fd.write(potential, header=inhead, extname="POTENTIAL_ASSIGNMENTS")

    # Now copy the original HDUs
    fd.write(fiber_data, header=inhead, extname="FASSIGN")
    fd.write(targets_data, header=inhead, extname="FTARGETS")
    fd.write(avail_data, header=inhead, extname="FAVAIL")
    fd.close()

    tm.stop()
    tm.report("  write data to file {}".format(tile_id))
    tm.clear()
    tm.start()

    # Try to encourage python to free some memory...
    del fd

    del avail_data
    del targets_data
    del fiber_data
    return


def merge_results(targetfiles, tiles, result_dir=".",
                  result_prefix="fiberassign", result_split_dir=False,
                  out_dir=None, out_prefix="tile", out_split_dir=False,
                  columns=None):
    """Merge target files and assignment output.

    Full target data is stored in shared memory and then multiple processes
    copy this data into the per-tile files.

    Args:
        targetfiles (list):  List of pathnames containing the original input
            target files.  The rows of the file MUST be sorted by TARGETID.
        tiles (list):  List of tile IDs to process.
        result_dir (str):  Top-level directory of fiberassign results.
        result_prefix (str):  Prefix of each per-tile file name.
        result_split_dir (bool):  Results are in split tile directories.
        out_dir (str):  Top-level directory for merged outputs.
        out_prefix (str):  Prefix of each per-tile output file name.
        out_split_dir (bool):  Write outputs in split tile directories.
        columns (list):  List of column names to propagate from the input
            target files (default is all).

    Returns:
        None.

    """
    # Load the full set of target files into memory.  Also build a mapping of
    # target ID to row index.  We assume that the result columns have the same
    # dtype in any of the target files.  We take the first target file and
    # construct the output recarray dtype from the columns in that file.
    out_dtype = None
    dcols = [(x, y) for x, y in merged_fiberassign_req_columns.items()]
    dcolnames = [x for x in merged_fiberassign_req_columns.keys()]

    tgdata = dict()
    tgdtype = dict()
    tgshape = dict()
    tghead = dict()

    for tf in targetfiles:
        tm = Timer()
        tm.start()
        fd = fitsio.FITS(tf)
        tghead[tf] = fd[1].read_header()
        # Allocate a shared memory buffer for the target data
        tglen = fd[1].get_nrows()
        tgshape[tf] = (tglen,)
        tgdtype[tf], tempoff, tempisvararray = fd[1].get_rec_dtype()
        tgbytes = tglen * tgdtype[tf].itemsize
        tgdata[tf] = RawArray("B", tgbytes)
        tgview = np.frombuffer(tgdata[tf],
                               dtype=tgdtype[tf]).reshape(tgshape[tf])
        # Read data directly into shared buffer
        tgview[:] = fd[1].read()

        # Sort rows by TARGETID if not already done
        tgviewids = tgview["TARGETID"]
        if not np.all(tgviewids[:-1] <= tgviewids[1:]):
            tgview.sort(order="TARGETID", kind="heapsort")

        tm.stop()
        tm.report("Read {} into shared memory".format(tf))

        # Add any missing columns to our output dtype record format.
        tfcols = list(tgview.dtype.names)
        if columns is not None:
            tfcols = [x for x in tfcols if x in columns]
        for col in tfcols:
            subd = tgview.dtype[col].subdtype
            colname = col
            if col in merged_fiberassign_swap:
                colname = merged_fiberassign_swap[col]
            if colname not in dcolnames:
                if subd is None:
                    dcols.extend([(colname, tgview.dtype[col].str)])
                else:
                    dcols.extend([(colname, subd[0], subd[1])])
                dcolnames.append(colname)

    out_dtype = np.dtype(dcols)

    # For each tile, find the target IDs used.  Construct the output recarray
    # and copy data into place.

    merge_tile = partial(merge_results_tile, out_dtype)

    if out_dir is None:
        out_dir = result_dir

    tile_map_list = [(x, result_path(x, dir=result_dir, prefix=result_prefix,
                                     split=result_split_dir),
                      result_path(x, dir=out_dir, prefix=out_prefix,
                                  create=True, split=out_split_dir))
                     for x in tiles]

    with mp.Pool(processes=default_mp_proc,
                 initializer=merge_results_tile_initialize,
                 initargs=(tgdata, tgdtype, tgshape)) as pool:
        results = pool.map(merge_tile, tile_map_list)

    return
