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

from desitarget.targetmask import desi_mask

from .utils import Logger, default_mp_proc

from .targets import (Targets, load_target_table)

from .assign import (result_tiles, result_path, avail_table_to_dict,
                     read_assignment_fits_tile)


def qa_parse_table(header, tgdata):
    """Extract target info from a table.
    """
    tgs = Targets()
    if "FA_SURV" in header:
        load_target_table(tgs, tgdata,
                          survey=str(header["FA_SURV"]).rstrip(),
                          typecol="FA_TYPE")
    else:
        load_target_table(tgs, tgdata)
    survey = tgs.survey()

    tgprops = dict()
    if survey == "main":
        lrgmask = int(desi_mask["LRG"].mask)
        elgmask = int(desi_mask["ELG"].mask)
        qsomask = int(desi_mask["QSO"].mask)
        badmask = int(desi_mask["BAD_SKY"].mask)
        for row in range(len(tgdata)):
            tgid = tgdata["TARGETID"][row]
            tg = tgs.get(tgid)
            dt = tg.bits
            tgprops[tgid] = dict()
            if dt & lrgmask:
                tgprops[tgid]["type"] = "LRG"
            elif dt & elgmask:
                tgprops[tgid]["type"] = "ELG"
            elif dt & qsomask:
                tgprops[tgid]["type"] = "QSO"
            elif dt & badmask:
                tgprops[tgid]["type"] = "BAD"
            else:
                tgprops[tgid]["type"] = "NA"
    else:
        # Could define similar things for other surveys here...
        for row in range(len(tgdata)):
            tgid = tgdata["TARGETID"][row]
            tgprops[tgid] = {"type": "NA"}
    return tgs, tgprops


def qa_tile(hw, tile_id, tgs, tgprops, tile_assign, tile_avail):
    props = dict()
    # tile_idx = tiles.order[tile_id]
    # props["tile_ra"] = tiles.ra[tile_idx]
    # props["tile_dec"] = tiles.dec[tile_idx]
    # props["tile_obscond"] = tiles.obscond[tile_idx]
    locs = np.array(hw.device_locations("POS"))
    nassign = 0
    nscience = 0
    nstd = 0
    nsky = 0
    nsafe = 0
    unassigned = list()
    objtypes = dict()
    for lid in locs:
        if lid not in tile_assign:
            unassigned.append(int(lid))
            continue
        tgid = tile_assign[lid]
        if tgid < 0:
            unassigned.append(int(lid))
            continue
        nassign += 1
        tg = tgs.get(tgid)
        if tg.is_science():
            nscience += 1
            ot = tgprops[tgid]["type"]
            if ot in objtypes:
                objtypes[ot] += 1
            else:
                objtypes[ot] = 1
        if tg.is_standard():
            nstd += 1
        if tg.is_sky():
            nsky += 1
        if tg.is_safe():
            nsafe += 1
    props["assign_total"] = nassign
    props["assign_science"] = nscience
    props["assign_std"] = nstd
    props["assign_sky"] = nsky
    props["assign_safe"] = nsafe
    for ot, cnt in objtypes.items():
        props["assign_obj_{}".format(ot)] = cnt
    props["unassigned"] = unassigned
    return props


def qa_tile_file(hw, params):
    (tile_id, tile_file) = params
    log = Logger.get()

    log.info("Processing tile {}".format(tile_id))

    header, fiber_data, targets_data, avail_data, gfa_data = \
        read_assignment_fits_tile((tile_id, tile_file))

    # Target properties
    tgs, tgprops = qa_parse_table(header, targets_data)

    # Only do QA on positioners.
    pos_rows = np.where(fiber_data["DEVICE_TYPE"] == b"POS")[0]

    # Target assignment
    tassign = {x["FIBER"]: x["TARGETID"] for x in fiber_data[pos_rows]
               if (x["FIBER"] >= 0)}

    tavail = avail_table_to_dict(avail_data)

    qa_data = qa_tile(hw, tile_id, tgs, tgprops, tassign, tavail)

    return qa_data


def qa_tiles(hw, tiles, result_dir=".", result_prefix="fiberassign_",
             result_split_dir=False, qa_out=None):
    """Run QA on a set of tiles.

    This will run QA on a set of output files.

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
