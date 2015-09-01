#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script provides a command line interface
for plotting data files with matplotlib
"""

from __future__ import division

import sys
from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def main(*sources, title='None', xlabel='None', ylabel='None', savefig='None',
         vertical_lines='None', labels='None', xscale='linear', yscale='linear',
         subplots='1,1', markers='osdv*', column_factors='None', columns='all',
         file_factors='None', sharex=False, sharey=False, figsize='None'):
    sp_nrows, sp_ncols = map(int, subplots.split(','))
    fig_kw = {

    }
    figsize = None if figsize == 'None' else tuple(map(float, figsize.split(',')))
    fig, axes = plt.subplots(sp_nrows, sp_ncols, squeeze=False,
                             sharex=sharex, sharey=sharey, figsize=figsize)
    for ax in axes.flatten():
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
    all_data = OrderedDict()
    all_interp = OrderedDict()
    prev_ncols = 0
    master_x = None
    for source in sources:
        all_data[source] = A = np.loadtxt(source)
        if master_x is None:
            master_x = A[:, 0]
        ncols = A.shape[1]
        if prev_ncols == 0:
            prev_ncols = ncols
        else:
            if ncols != prev_ncols:
                raise ValueError("Mismatching number of columns")
        all_interp[source] = interp1d(A[:, 0], A[:, 1:], axis=0)

    if file_factors != 'None':
        file_factors = [float(_) for _ in file_factors.split(',')]
        A = np.zeros((master_x.size, prev_ncols))
        A[:, 0] = master_x
        for factor, data in zip(file_factors, all_interp.values()):
            A[:, 1:] += factor*data(master_x)
        all_data = {'//all': A}
        all_interp = {'//all': interp1d(A[:, 0], A[:, 1:], axis=0)}

    if columns == 'all':
        columns = list(range(1, prev_ncols))
    else:
        columns = [int(_) for _ in columns.split(',')]

    if column_factors != 'None':
        column_factors = [float(_) for _ in column_factors.split(',')]
        if len(column_factors) != len(columns):
            raise ValueError("Mismatching number of factors and columns")
        if len(columns) > 1:
            if len(axes) != len(all_interp):
                raise ValueError("Mismatching number of subplots and files")
        for col_idx, ax in zip(columns, axes.flatten()):
            y = np.zeros((1, 1))
            for source in sources:
                raise NotImplementedError
                master = None
            ax.plot()
        # assert len(sources) == 2
        # data1 = data[sources[0]]
        # data2 = data[sources[1]]
        # assert data1.shape == data2.shape
        # data = data1.copy()
        # data[:, 1:] -= interp1d(data2[:, 0], data2[:, 1:], axis=0)(data1[:, 0])
        # for i in range(1, data.shape[1]):
        #     plt.plot(data[:, 0], np.abs(data[:, i]))
    else:
        for ax, (source, data) in zip(axes.flatten(), all_data.items()):
            for col_idx in columns:
                ax.plot(data[:, 0], data[:, col_idx])

    if vertical_lines != 'None':
        for vline in [float(_) for _ in vertical_lines.split(',')]:
            for ax in axes:
                ax.plot([vline]*2, ax.get_ylim())

    if xlabel != 'None':
        plt.xlabel(xlabel)

    if ylabel != 'None':
        plt.ylabel(ylabel)

    if labels != 'None':
        plt.legend(loc='best')

    plt.draw()

    if len(axes) > 1:
        plt.tight_layout()

    if savefig == 'None':
        plt.show()
    else:
        plt.savefig(savefig)
    return 0


if __name__ == '__main__':
    import argh
    argh.dispatch_command(main)
