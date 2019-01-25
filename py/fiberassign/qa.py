# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.qa
=======================

Quality Assurance tools.

"""
from __future__ import absolute_import, division, print_function

import os
import re
import numpy as np

import multiprocessing as mp
from functools import partial

import json

from .utils import Logger, default_mp_proc

from .targets import (Targets, append_target_table)

from .assign import (result_tiles, result_path, avail_table_to_dict,
                     read_assignment_fits_tile)


def qa_parse_table(tgdata):
    """Extract target info from a table.
    """
    tgs = Targets()
    typecol = None
    if "FBATYPE" not in tgdata.dtype.names:
        typecol = "DESI_TARGET"
    append_target_table(tgs, tgdata, typecol=typecol)
    return tgs


def qa_tile(hw, tile_id, tgs, tile_assign, tile_avail):
    props = dict()
    # tile_idx = tiles.order[tile_id]
    # props["tile_ra"] = tiles.ra[tile_idx]
    # props["tile_dec"] = tiles.dec[tile_idx]
    # props["tile_obscond"] = tiles.obscond[tile_idx]
    fibers = np.array(hw.fiber_id)
    nassign = 0
    nscience = 0
    nstd = 0
    nsky = 0
    unassigned = list()
    for fid in fibers:
        if fid not in tile_assign:
            unassigned.append(int(fid))
            continue
        tgid = tile_assign[fid]
        if tgid < 0:
            unassigned.append(int(fid))
            continue
        nassign += 1
        tg = tgs.get(tgid)
        if tg.is_science():
            nscience += 1
        if tg.is_standard():
            nstd += 1
        if tg.is_sky():
            nsky += 1
    props["unassigned"] = unassigned
    props["assign_total"] = nassign
    props["assign_science"] = nscience
    props["assign_std"] = nstd
    props["assign_sky"] = nsky
    return props


def qa_tile_file(hw, params):
    (tile_id, tile_file) = params
    log = Logger.get()

    log.info("Processing tile {}".format(tile_id))

    header, fiber_data, targets_data, avail_data, gfa_data = \
        read_assignment_fits_tile((tile_id, tile_file))

    # Target properties
    tgs = qa_parse_table(targets_data)

    # Target assignment
    tassign = {x["FIBER"]: x["TARGETID"] for x in fiber_data
               if (x["FIBER"] >= 0)}

    tavail = avail_table_to_dict(avail_data)

    qa_data = qa_tile(hw, tile_id, tgs, tassign, tavail)

    return qa_data


def qa_tiles(hw, tiles, result_dir=".", result_prefix="fiberassign_",
             result_split_dir=False, qa_out=None):
    """Plot assignment output.

    Args:
        hw (Hardware):  the hardware description.
        tiles (list):  List of tile IDs to process.
        result_dir (str):  Top-level directory of fiberassign results.
        result_prefix (str):  Prefix of each per-tile file name.
        result_split_dir (bool):  Results are in split tile directories.
        qa_out (str):  Override the output path.

    Returns:
        None.

    """
    log = Logger.get()

    foundtiles = result_tiles(dir=result_dir, prefix=result_prefix)
    log.info("Found {} fiberassign tile files".format(len(foundtiles)))

    if qa_out is None:
        qa_out = os.path.join(result_dir, "qa.json")

    qa_tile = partial(qa_tile_file, hw)

    avail_tiles = np.array(tiles.id)
    select_tiles = [x for x in foundtiles if x in avail_tiles]

    tile_map_list = [(x, result_path(x, dir=result_dir, prefix=result_prefix,
                                     split=result_split_dir))
                     for x in select_tiles]

    log.info("Selecting {} fiberassign tile files".format(len(tile_map_list)))
    with mp.Pool(processes=default_mp_proc) as pool:
        qa_result = pool.map(qa_tile, tile_map_list)

    qa_full = dict()
    for tid, props in zip(foundtiles, qa_result):
        tidx = tiles.order[tid]
        props["tile_ra"] = float(tiles.ra[tidx])
        props["tile_dec"] = float(tiles.dec[tidx])
        props["tile_obscond"] = int(tiles.obscond[tidx])
        qa_full[tid] = props

    with open(qa_out, "w") as f:
        json.dump(qa_full, f, indent=4, sort_keys=True)

    return
