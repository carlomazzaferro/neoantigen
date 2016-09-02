class BlastConservationCalculator:

    def __init__(self, matrix_name="blosum62"):
        self._subs_mat = getattr(MatrixInfo, matrix_name)
        self._no_use_thresh = 0.95

    def conservation_dict(self, blast_rec):
        cons_dict = {}
        rec_size = int(blast_rec.query_letters)
        for base_index in range(rec_size):
            cons_dict[base_index] = []
            for align in blast_rec.alignments:
                for hsp in align.hsps:
                    if (float(hsp.identities) / float(rec_size) <= self._no_use_thresh):
                        cons_dict = self._add_hsp_conservation(hsp, cons_dict)

        return cons_dict

    def _add_hsp_conservation(self, hsp, cons_dict):
        start_index = int(hsp.query_start) - 1
        hsp_index = 0
        for q_index in range(len(hsp.query)):
            if hsp.query[q_index] != '-':
                if hsp.sbjct[q_index] != '-':
                    try:
                        sub_val = self._subs_mat[(hsp.query[q_index], hsp.sbjct[q_index])]
                    except KeyError:
                        sub_val = self._subs_mat[(hsp.sbjct[q_index],  hsp.query[q_index])]
                        cons_dict[start_index + hsp_index].append(sub_val)
            hsp_index += 1
        return cons_dict

        """
        The result of this work is a dictionary of score conservation at each position.
        If you plot the average of these scores directly, it results in a very choppy graph
        which is hard to interpret. Luckily, Andrew Dalke has tackled this problem and presented
        a detailed writeup of smoothing scores for plotting. Using the Savitzky-Golay technique
        described there, the smoothed average of the results are plotted using matplotlib:

        """

import pylab

window_size = 29
data_smoother = SavitzkyGolayDataSmoother(window_size)
pos_data = []
cons_data = []
for pos in indexes:
    pos_data.append(pos + 1)
if len(cons_dict[pos]) > 0:
    cons_data.append(numpy.median(cons_dict[pos]))
else:
    cons_dict.append(0)


smooth_data = data_smoother.smooth_values(cons_data)
smooth_pos_data = pos_data[data_smoother.half_window():
len(pos_data) - data_smoother.half_window()]
pylab.plot(smooth_pos_data, smooth_data)
pylab.axis(xmin=min(pos_data), xmax=max(pos_data))
pylab.xlabel("Amino acid position")
pylab.ylabel("Conservation")
pylab.savefig('%s_conservation.png' % accession.replace(".", "_"))
