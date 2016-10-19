import numpy as np
import matplotlib.pyplot as plt
from collections import Iterable
from mpl_toolkits.axes_grid1 import make_axes_locatable


class Plotter(object):
    def __init__(self, prot_list, scores_list, df_list):

        self.scores = scores_list[0:-1]
        self.prot_list = prot_list[0:-1]
        self.df_list = df_list
        self.plot_data = self._get_plot_data()

    def plotting_function(self):

        plt.figure(figsize=(15, 25))
        ax = plt.gca()
        ax.set_ylim([-1, len(self.prot_list) + 1])
        ax.invert_yaxis()
        plt.title('High Affinity Score Locations')
        # colors = cm.rainbow(np.linspace(0, 1, len(plot_data)))
        colors = self._get_color_map()

        fig, ax = plt.subplots()
        cax = fig.add_axes([-100, max(self.plot_data[0][0]) + 500, -1, len(self.prot_list)])
        im = ax.imshow(colors, cmap='rainbow_r')
        fig.colorbar(im, cax=cax, orientation='horizontal')

        for i, coordinate in enumerate(self.plot_data):
            plt.scatter(coordinate[0], coordinate[1], c=colors[i],
                        linewidth='0', s=50, vmin=-1, vmax=1, label=self.prot_list[i], cmap='rainbow_r')
            # plt.scatter(coordinate[0], coordinate[1], color=colors, label=prot_list[i])


            plt.legend(self.scores, bbox_to_anchor=(1.5, 1), loc=1, ncol=1)
            # plt.legend

    def plotting_function_updated(self):

        plt.figure(figsize=(15, 25))
        ax = plt.gca()

        ax.set_ylim([-1, len(self.plot_data) + 1])
        ax.invert_yaxis()
        plt.title('High Affinity Score Locations')

        new_scores = self.scores
        new_scores[0] = 0.0

        colors = self._get_color_map(new_scores)
        legend_data = list(zip(self.prot_list, new_scores))
        legend_data = ["".join([str(i[0]), '; Score: ', str(colors[idx][0])]) for idx, i in enumerate(legend_data)]

        for i, coordinate in enumerate(self.plot_data[0:-1]):
            plt.scatter(coordinate[0], coordinate[1], c=colors[i],
                        linewidth='0', s=50, vmin=0, vmax=1, cmap='rainbow_r', label=legend_data[i])
            plt.legend(bbox_to_anchor=(1.5, 1), loc=1, ncol=1)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        bar = plt.colorbar(cax=cax)
        bar.ax.invert_yaxis()

    # ('Difference between peptide normalized on a 0-1 scale')


    def _get_color_map(self, scores):

        color_list = []

        for i, color in enumerate(scores):
            interpolated_col = np.interp(color, [0, max(scores)], [0, 1])
            full_list = np.full(len(self.plot_data[i][0]), interpolated_col)
            color_list.append(full_list)

        return color_list

    def _get_plot_data(self):

        plotting_data = []

        for idx, pep in enumerate(self.prot_list):
            for i in self.df_list:

                df_protein = i.ID.unique()[0]
                if pep == df_protein:
                    ranges_list = self._get_ranges_list(i)
                    to_plot = self._add_x_axis_data(ranges_list, idx)
                    plotting_data.append(to_plot)

        return plotting_data

    @staticmethod
    def _add_x_axis_data(ranges_list, idx):
        return [ranges_list, [idx + 1] * len(ranges_list)]

    def _get_ranges_list(self, df):

        from_df = list(df.Range.values)
        flattened = list(self.flatten(from_df))
        concatd_ranges_list = list(set(flattened))

        return concatd_ranges_list

    def flatten(self, l):
        for el in l:
            if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
                yield from self.flatten(el)
            else:
                yield el

