"""
Author: Anirban Mukhopadhyay, Prof. Anil's magnonics group, IITM, India
Created on: April 9, 2021
Version: 0.2
"""

import numpy as np
import matplotlib.pyplot as plt
import os as os
import fnmatch as fn
from scipy.signal import windows as win


class SWLib:
    class CalcSWProp:
        def __init__(self, h_0, k, film_thickness, saturation_magnetization, sw_type, gamma_0=0.0352, lambda_ex=None,
                     exchange=False):
            self.h_0 = h_0
            self.k = k
            self.dk = np.average(np.diff(k))
            self.film_thickness = film_thickness
            self.ms = saturation_magnetization
            self.sw_type = sw_type
            self.gamma_0 = gamma_0
            self.exch = exchange
            self.f_0_bar = self.h_0 / self.ms
            if self.exch is True:
                self.lambda_ex = lambda_ex
                self.f_ex_bar = self.f_0_bar + self.lambda_ex * self.k ** 2
            else:
                pass

        def bvsw_disp(self):
            with np.errstate(divide='ignore', invalid='ignore'):
                if self.exch is True:
                    buffer = np.sqrt(self.f_ex_bar ** 2 + self.f_ex_bar * (
                            1 - np.exp(-np.abs(self.k) * self.film_thickness)) / (np.abs(self.k) * self.film_thickness))
                else:
                    buffer = np.sqrt(self.f_0_bar ** 2 + self.f_0_bar * (
                            1 - np.exp(-np.abs(self.k) * self.film_thickness)) / (np.abs(self.k) * self.film_thickness))
            buffer[np.where(np.isnan(buffer))] = np.sqrt(self.f_0_bar ** 2 + self.f_0_bar)
            return buffer

        def fvsw_disp(self):
            with np.errstate(divide='ignore', invalid='ignore'):
                if self.exch is True:
                    buffer = np.sqrt(self.f_ex_bar ** 2 + self.f_ex_bar - self.f_ex_bar * (
                            1 - np.exp(-np.abs(self.k) * self.film_thickness)) / (np.abs(self.k) * self.film_thickness))
                else:
                    buffer = np.sqrt(self.f_0_bar ** 2 + self.f_0_bar - self.f_0_bar * (
                            1 - np.exp(-np.abs(self.k) * self.film_thickness)) / (np.abs(self.k) * self.film_thickness))
            buffer[np.where(np.isnan(buffer))] = self.f_0_bar
            return buffer

        def ssw_disp(self):
            with np.errstate(divide='ignore', invalid='ignore'):
                if self.exch is True:
                    buffer = np.sqrt(
                        self.f_ex_bar ** 2 + self.f_ex_bar + (
                                    1 - np.exp(-2 * np.abs(self.k) * self.film_thickness)) / 4)
                else:
                    buffer = np.sqrt(
                        self.f_0_bar ** 2 + self.f_0_bar + (1 - np.exp(-2 * np.abs(self.k) * self.film_thickness)) / 4)
            buffer[np.where(np.isnan(buffer))] = np.sqrt(self.f_0_bar ** 2 + self.f_0_bar)
            return buffer

        def calc_dispersion_curve_lowest_order_mode(self):
            dict = {"bvsw": self.bvsw_disp, "fvsw": self.fvsw_disp, "ssw": self.ssw_disp}
            func = dict.get(self.sw_type, "Invalid spin-wave class")
            fig, ax = plt.subplots(1, 1, squeeze=True)
            ax.plot(self.k, func(), "b--")
            ax.set_xlabel(r"$k$ (rad/nm)", fontsize=15)
            ax.set_ylabel(r"$f/f_{\rm{M}}$", fontsize=15)
            ax.tick_params(labelsize=12)

        def calc_grvel_lowest_order_mode(self):
            dict = {"bvsw": self.bvsw_disp, "fvsw": self.fvsw_disp, "ssw": self.ssw_disp}
            func = dict.get(self.sw_type, "Invalid spin-wave class")
            group_velocity = np.diff(func()) / self.dk
            fig, ax = plt.subplots(1, 1, squeeze=True)
            ax.plot(self.k[:-1], group_velocity, "b--")
            ax.set_xlabel(r"$k$ (rad/nm)", fontsize=15)
            ax.set_ylabel(r"$v_{\rm{g}}/f_{\rm{M}}$ (nm/rad)", fontsize=15)
            ax.tick_params(labelsize=12)

    class MagOvfPostProcessor:
        def __init__(self, folder_path, run_time, dt, x_in, length, dx, window_choice,
                     zero_pad=None, temporal_limits=None, spatial_limits=None):
            # Optional arguments value assignment if None
            if zero_pad is None:
                zero_pad = [2048, 2048]

            if temporal_limits is None:
                temporal_limits = [0, run_time]
            if spatial_limits is None:
                spatial_limits = [x_in, x_in + length]

            self.path = folder_path
            self.run_time, self.dt, self.n_t = run_time, dt, int(run_time / dt + 1)
            self.x_in, self.length, self.dx, self.n_x = x_in, length, dx, int(length / dx)
            self.window_choice = window_choice
            self.zero_pad = zero_pad
            self.temporal_limits = temporal_limits
            self.spatial_limits = spatial_limits

        def importer(self):
            my_file_list = []  # creating a empty list to store paths
            with os.scandir(self.path) as entries:
                for entry in entries:
                    if fn.fnmatch(entry.name, 'm*.ovf'):
                        my_file_list.append(entry.name)

            my_file_list = sorted(my_file_list)  # sorting filenames in increasing order of current value
            mag_data_spacetime = []
            [mag_data_spacetime.append(np.loadtxt(self.path + file, usecols=1)) for file in
             my_file_list]  # choosing y-component

            # reshape(no. if time instances (axis 0), no. of cells along y-axis (axis 1), no. of cells along x-axis (axis
            # 2))
            mag_data_spacetime = np.array(mag_data_spacetime).reshape(self.n_t, -1, self.n_x)
            return mag_data_spacetime

        def slicer(self):
            mag_data_spacetime = self.importer()
            i_upr, i_lwr = int(self.temporal_limits[1] / self.dt), int(self.temporal_limits[0] / self.dt)
            j_upr, j_lwr = int((self.spatial_limits[1] - self.x_in) / self.dx), int(
                (self.spatial_limits[0] - self.x_in) / self.dx)
            return mag_data_spacetime[i_lwr: i_upr + 1, :, j_lwr: j_upr + 1], mag_data_spacetime[0, :, j_lwr: j_upr + 1]

        def kaiser2d(self, shape, beta=14):
            return np.outer(win.kaiser(shape[0], beta), win.kaiser(shape[1], beta))

        def chebwin2d(self, shape, att=95):
            return np.outer(win.chebwin(shape[0], att), win.chebwin(shape[1], att))

        def window2d(self, shape):
            dict = {"kai": self.kaiser2d, "che": self.chebwin2d}
            func = dict.get(self.window_choice, "Invalid window method")
            return func(shape)

        def pre_conditioner(self):
            mag_data_spacetime, init_mag_data_spacetime = self.slicer()
            # removing the intial magnetization state
            mag_data_spacetime = np.array([mag_data_spacetime[i, :, :] - init_mag_data_spacetime for i in
                                           range(mag_data_spacetime.shape[0])])
            # averaging over y-axis & row no. = no. of time instances
            # column no. = no. of cells along x-axis
            mag_data_spacetime = np.average(mag_data_spacetime, axis=1)
            mag_data_spacetime *= self.window2d(mag_data_spacetime.shape)
            mag_pwr = 20 * np.log10(np.abs(np.fft.fftshift(np.fft.fft2(
                mag_data_spacetime, s=self.zero_pad, norm="ortho"))))  # dBm
            fmin, fmax = np.amin(np.fft.fftshift(np.fft.fftfreq(mag_pwr.shape[0], self.dt))), \
                         np.amax(np.fft.fftshift(np.fft.fftfreq(mag_pwr.shape[0], self.dt)))  # GHz
            kmin, kmax = np.amin(np.fft.fftshift(np.fft.fftfreq(mag_pwr.shape[1], self.dx) * 2 * np.pi)), \
                         np.amax(np.fft.fftshift(np.fft.fftfreq(mag_pwr.shape[1], self.dx) * 2 * np.pi))  # rad/nm
            return [kmin, kmax], [fmin, fmax], mag_pwr

        def spatio_temporal_plottr(self, axs_label, axs_range, cmap_choice="bwr", cbar_limits=None):
            mag_data_spacetime, init_mag_data_spacetime = self.slicer()
            # removing the intial magnetization state
            mag_data_spacetime = np.array([mag_data_spacetime[i, :, :] - init_mag_data_spacetime for i in
                                           range(mag_data_spacetime.shape[0])])
            # averaging over y-axis & row no. = no. of time instances
            # column no. = no. of cells along x-axis
            mag_data_spacetime = np.average(mag_data_spacetime, axis=1)
            boundary = [self.spatial_limits, self.temporal_limits]
            boundary = list(np.concatenate(boundary).flat)
            fig, ax = plt.subplots(1, 1)
            if cbar_limits is None:
                im = ax.imshow(mag_data_spacetime, cmap=cmap_choice, vmin=np.amin(mag_data_spacetime),
                               vmax=np.amax(mag_data_spacetime), interpolation='antialiased',
                               extent=boundary, origin='lower', aspect='auto')
            else:
                im = ax.imshow(mag_data_spacetime, cmap=cmap_choice, vmin=cbar_limits[0],
                               vmax=cbar_limits[1], interpolation='antialiased',
                               extent=boundary, origin='lower', aspect='auto')
            ax.set_xlabel(r"$" + axs_label[0] + "$ (nm)", fontsize=15)
            ax.set_ylabel(r"$" + axs_label[1] + "$ (ns)", fontsize=15)
            ax.set_xlim(axs_range[:2])
            ax.set_ylim(axs_range[2:])
            cbar = plt.colorbar(im, shrink=0.95)
            cbar.set_label(label=r"Magnetization (1)", size=12)
            ax.tick_params(labelsize=12)

        def dispersion_plottr(self, axs_label, axs_range, cmap_choice="inferno", cbar_limits=None):
            k, f, mag_pwr = self.pre_conditioner()
            boundary = [k, f]
            boundary = list(np.concatenate(boundary).flat)
            fig, ax = plt.subplots(1, 1)
            if cbar_limits is None:
                im = ax.imshow(mag_pwr, cmap=cmap_choice, vmin=np.amin(mag_pwr),
                                vmax=np.amax(mag_pwr), interpolation='antialiased',
                                extent=boundary, origin='lower', aspect='auto')
            else:
                im = ax.imshow(mag_pwr, cmap=cmap_choice, vmin=cbar_limits[0],
                                vmax=cbar_limits[1], interpolation='antialiased',
                                extent=boundary, origin='lower', aspect='auto')
            ax.set_xlabel(r"$" + axs_label[0] + "$ (rad/nm)", fontsize=15)
            ax.set_ylabel(r"$" + axs_label[1] + "$ (GHz)", fontsize=15)
            ax.set_xlim(axs_range[:2])
            ax.set_ylim(axs_range[2:])
            cbar = plt.colorbar(im, shrink=0.95)
            cbar.set_label(label=r"Power (dB)", size=12)
            ax.tick_params(labelsize=12)


