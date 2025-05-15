#!/usr/bin/env python3
#
# setup.py: this file is part of the EZmock python package.
#
# EZmock: Effective Zel'dovich approximation mock generator.
#
# Github repository:
#       https://github.com/cheng-zhao/EZmock
#
# Copyright (c) 2023-2025 Cheng Zhao <zhaocheng03@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

import os


class pyclustering:
  """A class to call clustering measurement codes for EZmock calibration."""

  def __init__(self, workdir='.', odir=None,
               pk=True, xi=True, bk=True,
               pk_exe=None, xi_exe=None, bk_exe=None,
               pk_conf=None, xi_conf=None, bk_conf=None,
               pk_ref=None, xi_ref=None, bk_ref=None,
               nthread=1, verbose=False):
    """
    Initialize the pyclustering class.

    Parameters
    ----------
    workdir : str, default: '.' (current directory)
        The working directory for the clustering measurements.
    odir : str, default: f'{workdir}/clustering'
        The directory where the clustering measurements are stored.
    pk : bool, default: True
        If True, the power spectrum will be measured.
    xi : bool, default: True
        If True, the 2-point correlation function will be measured.
    bk : bool, default: True
        If True, the bispectrum will be measured.
    pk_exe : str, optional
        The path to the executable for power spectrum measurement.
        Available at https://github.com/cheng-zhao/powspec.
    xi_exe : str, optional
        The path to the executable for 2-point correlation function estimator.
        Available at https://github.com/cheng-zhao/FCFC.
    bk_exe : str, optional
        The path to the executable for bispectrum measurement.
        Available upon request to Cheng Zhao (zhaocheng03@gmail.com).
    pk_conf : str, default: f'{workdir}/conf/powspec.conf'
        The configuration file for power spectrum measurement.
    xi_conf : str, default: f'{workdir}/conf/fcfc_2pt_box.conf'
        The configuration file for 2-point correlation function measurement.
    bk_conf : str, default: f'{workdir}/conf/bispec.conf'
        The configuration file for bispectrum measurement.
    pk_ref : str, optional
        The reference power spectrum.
    xi_ref : str, optional
        The reference 2-point correlation function.
    bk_ref : str, optional
        The reference bispectrum.
    nthread : int, default: 1
        The number of threads to use for the measurements.
    verbose : bool, default: False
        If True, verbose output will be printed.
    """
    self.workdir = workdir
    self.pk = pk
    self.xi = xi
    self.bk = bk
    self.pk_exe = pk_exe
    self.xi_exe = xi_exe
    self.bk_exe = bk_exe
    self.pk_ref = pk_ref
    self.xi_ref = xi_ref
    self.bk_ref = bk_ref

    try:
      self.nthread = int(nthread)
    except ValueError:
      raise ValueError("nthread must be an integer.")
    if self.nthread <= 0:
      raise ValueError("nthread must be positive.")
    self.verbose = verbose

    self.pk_conf = pk_conf if pk_conf else f'{self.workdir}/conf/powspec.conf'
    self.xi_conf = xi_conf if xi_conf else f'{self.workdir}/conf/fcfc_2pt_box.conf'
    self.bk_conf = bk_conf if bk_conf else f'{self.workdir}/conf/bispec.conf'

    if self.pk:
      if self.pk_exe is None or not os.path.isfile(self.pk_exe):
        raise ValueError("Power spectrum executable not found.")
      if not os.access(self.pk_exe, os.X_OK):
        raise ValueError("Power spectrum executable is not executable.")
      if self.pk_conf is None or not os.path.isfile(self.pk_conf):
        raise ValueError("Power spectrum configuration file not found.")
      if self.pk_ref is None or not os.path.isfile(self.pk_ref):
        raise ValueError("Power spectrum reference file not found.")
      
    if self.xi:
      if self.xi_exe is None or not os.path.isfile(self.xi_exe):
        raise ValueError("2-point correlation function executable not found.")
      if not os.access(self.xi_exe, os.X_OK):
        raise ValueError("2-point correlation function executable is not executable.")
      if self.xi_conf is None or not os.path.isfile(self.xi_conf):
        raise ValueError("2-point correlation function configuration file not found.")
      if self.xi_ref is None or not os.path.isfile(self.xi_ref):
        raise ValueError("2-point correlation function reference file not found.")
      
    if self.bk:
      if self.bk_exe is None or not os.path.isfile(self.bk_exe):
        raise ValueError("Bispectrum executable not found.")
      if not os.access(self.bk_exe, os.X_OK):
        raise ValueError("Bispectrum executable is not executable.")
      if self.bk_conf is None or not os.path.isfile(self.bk_conf):
        raise ValueError("Bispectrum configuration file not found.")
      if self.bk_ref is None or not os.path.isfile(self.bk_ref):
        raise ValueError("Bispectrum reference file not found.")

    self.odir = odir if odir else f'{self.workdir}/clustering'

    if not os.path.exists(self.odir):
      try:
        os.makedirs(self.odir)
      except OSError as e:
        raise OSError(f"Failed to create directory {self.odir}: {e}")


  def run(self, fname, Lbox):
    """
    Run the clustering measurements.

    Parameters
    ----------
    fname : str
        The name of the input catalogue.
    Lbox : float
        Side length of the simulation box.
    """
    from subprocess import Popen, PIPE

    if fname is None or not os.path.isfile(fname):
      raise ValueError("Input catalogue file not found.")
    bname = os.path.basename(fname)
    
    try:
      Lbox = float(Lbox)
    except ValueError:
      raise ValueError("Lbox must be a real number.")
    if Lbox <= 0:
      raise ValueError("Lbox must be positive.")
    
    omp_threads = os.environ.get('OMP_NUM_THREADS')
    os.environ['OMP_NUM_THREADS'] = str(self.nthread)
    
    if self.pk:
      pk_ofile = f'{self.odir}/PK_{bname}'
      if not os.path.exists(pk_ofile):
        pk_coord = f'[($1+{Lbox:g})%{Lbox:g},($2+{Lbox:g})%{Lbox:g},($3+{Lbox:g})%{Lbox:g}]'
        pk_cmd = [self.pk_exe, '-c', self.pk_conf, '-d', fname, '-B', f'{Lbox:g}', '-p', pk_coord, '-a', pk_ofile]
        if self.verbose: pk_cmd.append('-v')

        process = Popen(pk_cmd, stdout=PIPE, stderr=PIPE, text=True)
        sts = process.wait()
        for line in process.stdout: print(line, end='')
        for line in process.stderr: print(line, end='')
        if sts != 0:
          raise RuntimeError(f"Power spectrum measurement failed with exit code {sts}.")
      
    if self.xi:
      xi_ofile = f'{self.odir}/2PCF_{bname}'
      if not os.path.exists(xi_ofile):
        xi_cmd = [self.xi_exe, '-c', self.xi_conf, '-i', fname, '-b', str(Lbox), '-P', f'{xi_ofile}.dd', '-E', f'{xi_ofile}.xi2d', '-M', xi_ofile]
        if self.verbose: xi_cmd.append('-v')

        process = Popen(xi_cmd, stdout=PIPE, stderr=PIPE, text=True)
        sts = process.wait()
        for line in process.stdout: print(line, end='')
        for line in process.stderr: print(line, end='')
        if sts != 0:
          raise RuntimeError(f"2-point correlation function measurement failed with exit code {sts}.")
      
    if self.bk:
      bk_ofile = f'{self.odir}/BK_{bname}'
      if not os.path.exists(bk_ofile):
        bk_cmd = [self.bk_exe, '-c', self.bk_conf, '-i', fname, '-B', str(Lbox), '-o', bk_ofile]
        if self.verbose: bk_cmd.append('--verbose')

        process = Popen(bk_cmd, stdout=PIPE, stderr=PIPE, text=True)
        sts = process.wait()
        for line in process.stdout: print(line, end='')
        for line in process.stderr: print(line, end='')
        if sts != 0:
          raise RuntimeError(f"Bispectrum measurement failed with exit code {sts}.")

    if omp_threads is not None:
      os.environ['OMP_NUM_THREADS'] = omp_threads


  def plot(self, fnames):
    """
    Plot the clustering measurements.

    Parameters
    ----------
    fnames : list of str
        The names of the input files to plot.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 3, figsize=(10, 6))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)

    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # Plot EZmock clustering
    for i, fname in enumerate(fnames):
      if i == len(fnames) - 1:
        color = 'k'
        label = 'current'
        zorder = i + 2
      else:
        color = color_cycle[i % len(color_cycle)]
        label = str(i)
        zorder = i + 1

      bname = os.path.basename(fname)
      if self.pk:
        pk_fname = f'{self.odir}/PK_{bname}'
        try:
          pk = np.loadtxt(pk_fname, unpack=True)
          ax[0, 0].plot(pk[0], pk[0]**1.5 * pk[5], c=color, label=label, zorder=zorder)
          ax[1, 0].plot(pk[0], pk[0]**1.5 * pk[6], c=color, label=label, zorder=zorder)
        except Exception:
          print(f"Failed to load power spectrum data from {pk_fname}.")

      if self.xi:
        xi_fname = f'{self.odir}/2PCF_{bname}'
        try:
          xi = np.loadtxt(xi_fname, unpack=True)
          ax[0, 1].plot(xi[0], xi[0]**2 * xi[3], c=color, label=label, zorder=zorder)
          ax[1, 1].plot(xi[0], xi[0]**2 * xi[4], c=color, label=label, zorder=zorder)
        except Exception:
          print(f"Failed to load 2-point correlation function data from {xi_fname}.")

      if self.bk:
        bk_fname = f'{self.odir}/BK_{bname}'
        try:
          bk = np.loadtxt(bk_fname, unpack=True)
          ax[0, 2].plot(bk[0] / np.pi, bk[4], c=color, label=label, zorder=zorder)
          ax[1, 2].plot(bk[0] / np.pi, bk[5], c=color, label=label, zorder=zorder)
        except Exception:
          print(f"Failed to load bispectrum data from {bk_fname}.")

    # Plot reference clustering
    label = 'ref'
    zorder = len(fnames)
    if self.pk:
      pk = np.loadtxt(self.pk_ref, unpack=True)
      ax[0, 0].plot(pk[0], pk[0]**1.5 * pk[5], 'k--', label=label, zorder=zorder)
      ax[1, 0].plot(pk[0], pk[0]**1.5 * pk[6], 'k--', label=label, zorder=zorder)

    if self.xi:
      xi = np.loadtxt(self.xi_ref, unpack=True)
      ax[0, 1].plot(xi[0], xi[0]**2 * xi[3], 'k--', label=label, zorder=zorder)
      ax[1, 1].plot(xi[0], xi[0]**2 * xi[4], 'k--', label=label, zorder=zorder)

    if self.bk:
      bk = np.loadtxt(self.bk_ref, unpack=True)
      ax[0, 2].plot(bk[0] / np.pi, bk[4], 'k--', label=label, zorder=zorder)
      ax[1, 2].plot(bk[0] / np.pi, bk[5], 'k--', label=label, zorder=zorder)

    # Set axis labels and titles
    ax[0, 0].set_xlabel(r'$k$ [h/Mpc]')
    ax[0, 0].set_ylabel(r'$k^{1.5} P_0 (k)$')
    ax[1, 0].set_xlabel(r'$k$ [h/Mpc]')
    ax[1, 0].set_ylabel(r'$k^{1.5} P_2 (k)$')
    ax[0, 1].set_xlabel(r'$s$ [Mpc/h]')
    ax[0, 1].set_ylabel(r'$s^2 \xi_0 (s)$')
    ax[1, 1].set_xlabel(r'$s$ [Mpc/h]')
    ax[1, 1].set_ylabel(r'$s^2 \xi_2 (s)$')
    ax[0, 2].set_xlabel(r'$\theta_{12} / \pi$')
    ax[0, 2].set_ylabel(r'$B (\theta_{12})$')
    ax[1, 2].set_xlabel(r'$\theta_{12} / \pi$')
    ax[1, 2].set_ylabel(r'$Q (\theta_{12})$')

    ax[1, 2].legend(loc='upper left')
    for i in range(2):
      for j in range(3):
        ax[i, j].grid(ls=':', c='dimgray', alpha=0.6)
