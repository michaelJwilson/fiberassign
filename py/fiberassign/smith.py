import time
import desimodel.io
import pylab as pl
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from   desitarget.targetmask import desi_mask
from   desimodel.focalplane.geometry import xy2radec
from   astropy.table import Table


def get_fao(tids, basedir='/project/projectdirs/desi/datachallenge/reference_runs/19.2/'):
    ##  fao_fibs, fao_tiles = get_fao(tids)

    fao_fibs   = []
    fao_tiles  = []

    for tid in tids:
      fname    = basedir + '/fiberassign/fiberassign_{0:05d}.fits'.format(tid)
      fao_fib  = Table.read(fname, 'FASSIGN')
      fao_fib.sort('FIBER')
      fao_fibs.append(fao_fib)

      fname    = basedir+ '/fiberassign/tile-{0:05d}.fits'.format(tid)
      fao_tile = Table.read(fname, 'POTENTIAL_ASSIGNMENTS')
      fao_tiles.append(fao_tile)

    return  fao_fibs, fao_tiles
      
def fa_plot(axarr, i, j, tids, fao_fibs, fao_tiles, tiles, targets, ttypes=['QSO'], notin=True, s=7):
    '''
    Plot the distribution of assigned and unassigned targets for a list of tiles. 

    Input:

    tids      --  Tile IDs to plot.
    num       --  Index of Tile to plot in tids.
    fao_fibs  --  Loaded fiberassign fiber-* files. 
    fao_tiles --  Loaded fiberassign tile-* files. 
    tiles     --  (Associated) Nominal tile info. 
    targets   --  instance of the targets class.
    ttypes    --  Target types to plot ['ELG', 'LRG', 'QSO', 'BGS']
    notin     --  Add markers for unassigned (more busy).
    '''

    colordict    = {'ELG': 'b', 'LRG': 'r', 'QSO': 'orange', 'BGS': 'y'}        
    tile_radius  = 1.65

    if isinstance(axarr, list) | isinstance(axarr, np.ndarray):
      ax         = axarr[i,j]

    else:
      ax         = axarr

    ##  Plot boundaries of all in tids. 
    for count, tid in enumerate(tids):
        tileloc  = np.where(tiles['TILEID'].quantity == tid)[0][0]
        ra, dec  = tiles['RA'][tileloc], tiles['DEC'][tileloc]
        c        = plt.Circle((ra, dec), tile_radius, facecolor='None', edgecolor='k')
        ax.add_artist(c)

    ##  Plot for main tile.
    tid          = tids[i]

    fao_fib      = fao_fibs[i]
    fao_tile     = fao_tiles[i]

    tileloc      = np.where(tiles['TILEID'].quantity == tid)[0][0]
    prog         = tiles['PROGRAM'][tileloc]

    ##  Centre plot on first tile in tids. 
    ctileloc     = np.where(tiles['TILEID'].quantity == tids[0])[0][0]

    tmp          = fao_tile['TARGETID']
    
    not_assigned = (~np.in1d(tmp, fao_fib['TARGETID']))
    assigned     = ( np.in1d(tmp, fao_fib['TARGETID']))

    dsorted      = np.argsort(targets['TARGETID'])
        
    ypos         = np.searchsorted(targets['TARGETID'][dsorted], tmp)
    inds         = np.zeros_like(ypos)

    inds[tmp <= np.max(targets['TARGETID'])] = dsorted[ypos[tmp <= np.max(targets['TARGETID'])]]

    for ttype in ttypes:
        ##  Enforce BGS on bright.  MWS?
        if prog == 'BRIGHT':
            ttype = 'BGS'

        istype    = (targets[inds]['DESI_TARGET'] & desi_mask.mask(ttype) != 0)

        ax.scatter(targets[inds]['RA'][assigned & istype & (tmp <= np.max(targets['TARGETID']))],
                   targets[inds]['DEC'][assigned & istype & (tmp <= np.max(targets['TARGETID']))],
                   s=s, edgecolor="none", c=colordict[ttype], zorder=10, label=ttype)
        
        if notin:
            ax.scatter(targets[inds]['RA'][not_assigned & istype & (tmp <= np.max(targets['TARGETID']))],
                       targets[inds]['DEC'][not_assigned & istype & (tmp <= np.max(targets['TARGETID']))],
                       s=s, edgecolor="none", c=colordict[ttype], marker="x", zorder=10, label='Unassigned', alpha=0.6)
        
        ax.set_xlim(tiles['RA'][ctileloc]  - tile_radius, tiles['RA'][ctileloc]  + tile_radius)
        ax.set_ylim(tiles['DEC'][ctileloc] - tile_radius, tiles['DEC'][ctileloc] + tile_radius)
        
        ax.legend(frameon=False, title='TILE %d:  %s' % (tid, prog), loc=1)

        ax.set_aspect('equal', 'datalim')
        ax.invert_xaxis()


if __name__ == '__main__':
    ##  Hardware.
    fiberpos = Table(desimodel.io.load_fiberpos())
    fiberpos.sort('FIBER')

    fiber_locations = sorted(zip(fiberpos['FIBER'], fiberpos['LOCATION']))

    ## Nominal tiling. 
    tiles     = Table(desimodel.io.load_tiles())
    
    ## Targets.
    targets   = Table.read('/project/projectdirs/desi/datachallenge/reference_runs/19.2/targets/mtl.fits')

    '''
    ## One tile, multiple targets. 
    ttypes   = ['LRG', 'QSO']
    tids     = np.array([1165], dtype=np.int32)

    fao_fibs, fao_tiles = get_fao(tids) 

    fig, axarr = plt.subplots(figsize=(12, 14))
    mpl.rcParams['figure.dpi'] = 500

    fa_plot(axarr, 0, 0, tids, fao_fibs, fao_tiles, tiles, targets, ttypes=ttypes)

    axarr.set_xlabel('Right ascension [deg.]')
    axarr.set_ylabel('Declination [deg.]')

    pl.show()
    '''

    ##  All tiles and targets. 
    ttypes   = ['LRG', 'QSO', 'ELG']
    tids     = np.array([1165, 6927], dtype=np.int32)

    fao_fibs, fao_tiles = get_fao(tids)

    fig, axarr = plt.subplots(len(tids), len(ttypes), figsize=(12, 14))
    mpl.rcParams['figure.dpi'] = 100

    for i, _ in enumerate(tids):
      for j, ttype in enumerate(ttypes):
        fa_plot(axarr, i, j, tids, fao_fibs, fao_tiles, tiles, targets, ttypes=[ttype])

    axarr[1,0].set_xlabel('Right ascension [deg.]')
    axarr[1,0].set_xlabel('Right ascension [deg.]')

    axarr[0,0].set_ylabel('Declination [deg.]')
    axarr[1,0].set_ylabel('Declination [deg.]')

    plt.tight_layout()
    pl.show()

    print('\n\nDone.\n\n')
