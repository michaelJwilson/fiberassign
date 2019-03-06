# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.scripts.assign
==============================

High-level functions for running assignment.

"""
from __future__ import absolute_import, division, print_function

import os
import sys
import argparse
from datetime import datetime

from desitarget.targetmask import desi_mask

from fiberassign.utils import GlobalTimers, Logger

from fiberassign.hardware import load_hardware

from fiberassign.tiles import load_tiles

from fiberassign.gfa import get_gfa_targets

from fiberassign.targets import (str_to_target_type, TARGET_TYPE_SCIENCE,
                                 TARGET_TYPE_SKY, TARGET_TYPE_STANDARD,
                                 TARGET_TYPE_SAFE, Targets, TargetsAvailable,
                                 TargetTree, FibersAvailable,
                                 load_target_file,
                                 get_sciencemask, get_stdmask, get_skymask,
                                 get_safemask, get_excludemask
                                 )

from fiberassign.assign import (Assignment, write_assignment_fits,
                                result_path)


def parse_assign(optlist=None):
    """Parse assignment options.

    This parses either sys.argv or a list of strings passed in.  If passing
    an option list, you can create that more easily using the
    :func:`option_list` function.

    Args:
        optlist (list, optional): Optional list of arguments to parse instead
            of using sys.argv.

    Returns:
        (namespace):  an ArgumentParser namespace.

    """
    log = Logger.get()
    parser = argparse.ArgumentParser()

    parser.add_argument("--targets", type=str, required=True, nargs="+",
                        help="Input file with targets of any type.  This "
                        "argument can be specified multiple times (for "
                        "example if standards / skies / science targets are "
                        "in different files).  By default, the "
                        "'--mask_column' (default DESI_TARGET)"
                        "column and bitfield values defined in desitarget "
                        "are used to determine the type of each target.  "
                        "Each filename may be optionally followed by comma "
                        "and then one of the strings 'science', 'standard', "
                        "'sky' or 'safe' to force all targets in that file "
                        "to be treated as a fixed target type.")

    parser.add_argument("--gfafile", type=str, required=False, default=None,
                        help="Optional GFA targets FITS file")

    parser.add_argument("--footprint", type=str, required=False, default=None,
                        help="Optional FITS file defining the footprint.  If"
                        " not specified, the default footprint from desimodel"
                        " is used.")

    parser.add_argument("--tiles", type=str, required=False, default=None,
                        help="Optional text file containing a subset of the"
                        " tile IDs to use in the footprint, one ID per line."
                        " Default uses all tiles in the footprint.")

    parser.add_argument("--positioners", type=str, required=False,
                        default=None,
                        help="Optional FITS file describing the fiber "
                        "positioner locations.  Default uses the file from "
                        "desimodel.")

    parser.add_argument("--status", type=str, required=False, default=None,
                        help="Optional fiber status file in astropy ECSV "
                        "format.  Default treats all fibers as good.")

    parser.add_argument("--rundate", type=str, required=False, default=None,
                        help="Optional date to simulate for this run of "
                        "fiber assignment, used with the fiber status file "
                        "to determine which fibers currently have problems.  "
                        "Default uses the current date.  Format is "
                        "YYYY-MM-DDTHH:mm:ss in UTC time.")

    parser.add_argument("--dir", type=str, required=False, default=None,
                        help="Output directory.")

    parser.add_argument("--prefix", type=str, required=False,
                        default="fiberassign_",
                        help="Prefix of each file (before the <tile>.fits).")

    parser.add_argument("--split", required=False, default=False,
                        action="store_true",
                        help="Split output into tile prefix directories.")

    parser.add_argument("--standards_per_petal", type=int, required=False,
                        default=10, help="Required number of standards per"
                        " petal")

    parser.add_argument("--sky_per_petal", type=int, required=False,
                        default=40, help="Required number of sky targets per"
                        " petal")

    parser.add_argument("--write_all_targets", required=False, default=False,
                        action="store_true",
                        help="When writing target properties, write data "
                        "for all available targets, not just those which are "
                        "assigned.  This is convenient, but increases the "
                        "write time and the file size.")

    parser.add_argument("--overwrite", required=False, default=False,
                        action="store_true",
                        help="Overwrite any pre-existing output files")

    parser.add_argument("--mask_column", required=False, default="DESI_TARGET",
                        help="Default FITS column to use for applying target "
                             "masks")

    parser.add_argument("--sciencemask", required=False,
                        default=get_sciencemask(),
                        help="Default DESI_TARGET mask to use for science "
                             "targets")

    parser.add_argument("--stdmask", required=False,
                        default=get_stdmask(),
                        help="Default DESI_TARGET mask to use for stdstar "
                             "targets")

    parser.add_argument("--skymask", required=False,
                        default=get_skymask(),
                        help="Default DESI_TARGET mask to use for sky targets")

    parser.add_argument("--safemask", required=False,
                        default=get_safemask(),
                        help="Default DESI_TARGET mask to use for safe "
                        "backup targets")

    parser.add_argument("--excludemask", required=False,
                        default=get_excludemask(),
                        help="Default DESI_TARGET mask to exclude from "
                        "any assignments")

    parser.add_argument("--by_tile", required=False, default=False,
                        action="store_true",
                        help="Do assignment one tile at a time.  This disables"
                        " redistribution.")

    args = None
    if optlist is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(optlist)

    # Allow sciencemask, stdmask, etc. to be int or string
    if isinstance(args.sciencemask, str):
        args.sciencemask = desi_mask.mask(args.sciencemask.replace(",", "|"))

    if isinstance(args.stdmask, str):
        args.stdmask = desi_mask.mask(args.stdmask.replace(",", "|"))

    if isinstance(args.skymask, str):
        args.skymask = desi_mask.mask(args.skymask.replace(",", "|"))

    if isinstance(args.safemask, str):
        args.safemask = desi_mask.mask(args.safemask.replace(",", "|"))

    if isinstance(args.excludemask, str):
        args.excludemask = desi_mask.mask(args.excludemask.replace(",", "|"))

    log.info("sciencemask {}".format(
        " ".join(desi_mask.names(args.sciencemask))))
    log.info("stdmask     {}".format(" ".join(desi_mask.names(args.stdmask))))
    log.info("skymask     {}".format(" ".join(desi_mask.names(args.skymask))))
    log.info("safemask    {}".format(" ".join(desi_mask.names(args.safemask))))
    log.info("excludemask {}".format(
        " ".join(desi_mask.names(args.excludemask))))

    # Set output directory
    if args.dir is None:
        args.dir = "out_fiberassign_{}".format(args.rundate)

    return args


def run_assign_init(args):
    """Initialize assignment inputs.

    This uses the previously parsed options to load the input files needed.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        (tuple):  The (Hardware, Tiles, Targets) needed to run assignment.

    """
    log = Logger.get()
    # Read hardware properties
    hw = load_hardware(fiberpos_file=args.positioners, rundate=args.rundate,
                       status_file=args.status)

    # Read tiles we are using
    tileselect = None
    if args.tiles is not None:
        tileselect = list()
        with open(args.tiles, "r") as f:
            for line in f:
                # Try to convert the first column to an integer.
                try:
                    tileselect.append(int(line.split()[0]))
                except ValueError:
                    pass
    tiles = load_tiles(tiles_file=args.footprint, select=tileselect)

    # Before doing significant calculations, check for pre-existing files
    if not args.overwrite:
        for tileid in tiles.id:
            outfile = result_path(tileid, dir=args.dir,
                                  prefix=args.prefix, split=args.split)
            if os.path.exists(outfile):
                outdir = os.path.split(outfile)[0]
                log.error("Output files already exist in {}".format(outdir))
                log.error("either remove them or use --overwrite")
                sys.exit(1)

    # Create empty target list
    tgs = Targets()

    # Append each input target file
    for tgarg in args.targets:
        tgprops = tgarg.split(",")
        tgfile = tgprops[0]
        typeforce = None
        if len(tgprops) > 1:
            # we are forcing the target type for this file
            typeforce = str_to_target_type(tgprops[1])
        load_target_file(tgs, tgfile, typeforce=typeforce,
                         typecol=args.mask_column,
                         sciencemask=args.sciencemask,
                         stdmask=args.stdmask,
                         skymask=args.skymask,
                         safemask=args.safemask,
                         excludemask=args.excludemask)
    return (hw, tiles, tgs)


def run_assign_full(args):
    """Run fiber assignment over all tiles simultaneously.

    This uses the previously parsed options to read input data and run through
    the typical assignment sequence doing one step at a time over all tiles.
    It then writes to the outputs to disk.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        None

    """
    gt = GlobalTimers.get()
    gt.start("run_assign_full calculation")

    # Load data
    hw, tiles, tgs = run_assign_init(args)

    # Create a hierarchical triangle mesh lookup of the targets positions
    tree = TargetTree(tgs)

    # Compute the targets available to each fiber for each tile.
    tgsavail = TargetsAvailable(hw, tgs, tiles, tree)

    # Free the tree
    del tree

    # Compute the fibers on all tiles available for each target and sky
    favail = FibersAvailable(tgsavail)

    # Create assignment object
    asgn = Assignment(tgs, tgsavail, favail)

    # First-pass assignment of science targets
    asgn.assign_unused(TARGET_TYPE_SCIENCE)

    # Redistribute science targets across available petals
    asgn.redistribute_science()

    # Assign standards, up to some limit
    asgn.assign_unused(TARGET_TYPE_STANDARD, args.standards_per_petal)
    asgn.assign_force(TARGET_TYPE_STANDARD, args.standards_per_petal)

    # Assign sky to unused fibers, up to some limit
    asgn.assign_unused(TARGET_TYPE_SKY, args.sky_per_petal)
    asgn.assign_force(TARGET_TYPE_SKY, args.sky_per_petal)

    # If there are any unassigned fibers, try to place them somewhere.
    asgn.assign_unused(TARGET_TYPE_SCIENCE)
    asgn.assign_unused(TARGET_TYPE_STANDARD)
    asgn.assign_unused(TARGET_TYPE_SKY)

    # NOTE:  This was removed since we are treating BAD_SKY as science targets
    # with very low priority.
    #
    # # Assign safe location to unused fibers (no maximum).  There should
    # # always be at least one safe location (i.e. "BAD_SKY") for each fiber.
    # # So after this is run every fiber should be assigned to something.
    # asgn.assign_unused(TARGET_TYPE_SAFE)

    # Assign sky monitor fibers
    asgn.assign_unused(TARGET_TYPE_SKY, -1, "ETC")

    gt.stop("run_assign_full calculation")
    gt.start("run_assign_full write output")

    # Make sure that output directory exists
    if not os.path.isdir(args.dir):
        os.makedirs(args.dir)

    # Optionally get GFA targets
    gfa_targets = None
    if args.gfafile is not None:
        gfa_targets = get_gfa_targets(tiles, args.gfafile)

    # Write output
    write_assignment_fits(tiles, asgn, out_dir=args.dir,
                          out_prefix=args.prefix, split_dir=args.split,
                          all_targets=args.write_all_targets,
                          gfa_targets=gfa_targets, overwrite=args.overwrite)

    gt.stop("run_assign_full write output")

    gt.report()

    return


def run_assign_bytile(args):
    """Run fiber assignment tile-by-tile.

    This uses the previously parsed options to read input data and run through
    the typical assignment sequence on a single tile before moving on to the
    next.  It then writes to the outputs to disk.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        None

    """
    gt = GlobalTimers.get()
    gt.start("run_assign_bytile calculation")

    # Load data
    hw, tiles, tgs = run_assign_init(args)

    # Create a hierarchical triangle mesh lookup of the targets positions
    tree = TargetTree(tgs)

    # Compute the targets available to each fiber for each tile.
    tgsavail = TargetsAvailable(hw, tgs, tiles, tree)

    # Free the tree
    del tree

    # Compute the fibers on all tiles available for each target and sky
    favail = FibersAvailable(tgsavail)

    # Create assignment object
    asgn = Assignment(tgs, tgsavail, favail)

    # We are now going to loop over tiles and assign each one fully before
    # moving on to the next

    for tile_id in tiles.id:

        # First-pass assignment of science targets
        asgn.assign_unused(TARGET_TYPE_SCIENCE, -1, "POS", tile_id, tile_id)

        # Assign standards, up to some limit
        asgn.assign_unused(TARGET_TYPE_STANDARD, args.standards_per_petal,
                           "POS", tile_id, tile_id)
        asgn.assign_force(TARGET_TYPE_STANDARD, args.standards_per_petal,
                          tile_id, tile_id)

        # Assign sky to unused fibers, up to some limit
        asgn.assign_unused(TARGET_TYPE_SKY, args.sky_per_petal, "POS",
                           tile_id, tile_id)
        asgn.assign_force(TARGET_TYPE_SKY, args.sky_per_petal,
                          tile_id, tile_id)

        asgn.assign_unused(TARGET_TYPE_STANDARD, -1, "POS", tile_id, tile_id)
        asgn.assign_unused(TARGET_TYPE_SKY, -1, "POS", tile_id, tile_id)

        # NOTE:  This was removed since we are treating BAD_SKY as science
        # targets with very low priority.
        #
        # # Assign safe location to unused fibers (no maximum).  There should
        # # always be at least one safe location (i.e. "BAD_SKY") for each
        # # fiber.  So after this is run every fiber should be assigned to
        # # something.
        # asgn.assign_unused(TARGET_TYPE_SAFE)

        # Assign sky monitor fibers
        asgn.assign_unused(TARGET_TYPE_SKY, -1, "ETC", tile_id, tile_id)

    gt.stop("run_assign_bytile calculation")
    gt.start("run_assign_bytile write output")

    # Make sure that output directory exists
    if not os.path.isdir(args.dir):
        os.makedirs(args.dir)

    # Optionally get GFA targets
    gfa_targets = None
    if args.gfafile is not None:
        gfa_targets = get_gfa_targets(tiles, args.gfafile)

    # Write output
    write_assignment_fits(tiles, asgn, out_dir=args.dir,
                          out_prefix=args.prefix, split_dir=args.split,
                          all_targets=args.write_all_targets,
                          gfa_targets=gfa_targets, overwrite=args.overwrite)

    gt.stop("run_assign_bytile write output")

    gt.report()

    return
