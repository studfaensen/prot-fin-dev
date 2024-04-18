from tools import *
from .algorithm import create_constellation, get_aa_vector
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as Cmap, Normalize
from matplotlib.cm import ScalarMappable
from multiprocessing import Pool


def plot_frequencies(prot_file: str, out_file: str, cpu_count=1):
    fasta = Fasta(prot_file)
    protein_count = len(fasta)

    with Pool(cpu_count) as p:
        subproc_slices = divide_evenly(protein_count, cpu_count)
        main_slice = next(subproc_slices)
        subprocesses = p.map_async(_process, ((fasta, slc) for slc in subproc_slices))
        sel_freqs = _process((fasta, main_slice))
        for res_sel_freqs in subprocesses.get():
            for res_freq, (res_freq_count, res_prot_count) in res_sel_freqs.items():
                freq_count, prot_count = sel_freqs.get(res_freq, (0, 0))
                sel_freqs[res_freq] = (freq_count + res_freq_count, prot_count + res_prot_count)

    if sel_freqs:
        cmap = Cmap.from_list(
            'custom_cmap',
            [(0, 'orange'), (1/3, 'red'), (1, 'black')]
        )
        norm = Normalize(vmin=0, vmax=protein_count)

        freqs, counts = zip(*sorted(sel_freqs.items()))
        freq_counts, prot_counts = zip(*counts)

        fig = plt.figure(figsize=(10, 5))
        plt.scatter(freqs, freq_counts, marker='o', c=list(map(lambda c: cmap(norm(c)), prot_counts)))
        for outlier_freq, (outlier_freq_count, _) in filter(lambda i: i[1][0] > max(freq_counts) / 2, sel_freqs.items()):
            plt.text(
                outlier_freq,
                outlier_freq_count - max(freq_counts) / 100 * 1.5,
                round(outlier_freq, 2),
                ha="center",
                va="top"
            )
        plt.colorbar(ScalarMappable(cmap=cmap, norm=norm)).set_label("Sequences", rotation=-90, va="bottom")
        plt.xlabel("Frequencies")
        plt.ylabel("Absolute count in all sequences")
        plt.title("Occurences of selected STFT frequences")
        plt.savefig(out_file)


def _process(args) -> Dict[float, Tuple[int, int]]:
    fasta, prot_slice = args
    selected_freqs = {}

    for _, _, seq in fasta[prot_slice]:
        constellation_map = create_constellation(get_aa_vector(seq))
        all_freqs = set() 
        for window in constellation_map:
            for freq, _ in window:
                freq_count, prot_count = selected_freqs.get(freq, (0, 0))
                selected_freqs[freq] = (freq_count + 1, prot_count)
                all_freqs.add(freq)

        for freq in all_freqs:
            freq_count, prot_count = selected_freqs[freq]
            selected_freqs[freq] = (freq_count, prot_count + 1)

    return selected_freqs
