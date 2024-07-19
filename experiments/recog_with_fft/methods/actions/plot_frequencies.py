from tools import *
from .algorithm import create_constellation, get_aa_vector
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as Cmap, Normalize
from matplotlib.cm import ScalarMappable
from multiprocessing import Pool
from scipy.fft import fftfreq


def plot_frequencies(prot_file: str, out_file: str, cpu_count=1, **kwargs):
    fasta = Fasta(prot_file)
    protein_count = len(fasta)
    db_config = DBConfig(**{k: v for k, v in kwargs.items() if k in DBConfig._fields})

    if cpu_count > 1:
        with Pool(cpu_count - 1) as p:
            subprocesses = p.map_async(_process, ((fasta, slice(i, None, cpu_count), db_config) for i in range(1, cpu_count)))
            sel_freqs, freqs_per_win, quantiles_per_win_without_first_ones = _process((fasta, slice(0, None, cpu_count), db_config))
            for res_sel_freqs, res_freqs_per_win, res_quantiles_per_win_without_first_ones in subprocesses.get():
                freqs_per_win.extend(res_freqs_per_win)
                quantiles_per_win_without_first_ones.extend(res_quantiles_per_win_without_first_ones)
                for res_freq, (res_freq_count, res_prot_count) in res_sel_freqs.items():
                    freq_count, prot_count = sel_freqs[res_freq]
                    sel_freqs[res_freq] = (freq_count + res_freq_count, prot_count + res_prot_count)
    else:
        sel_freqs, freqs_per_win, quantiles_per_win_without_first_ones = _process((fasta, slice(None), db_config))

    if sel_freqs:
        cmap = Cmap.from_list(
            'custom_cmap',
            [(0, 'orange'), (1/3, 'red'), (1, 'black')]
        )
        norm = Normalize(vmin=0, vmax=protein_count)

        freqs, counts = zip(*sorted((freq_count, prot_count) for freq_count, prot_count in sel_freqs.items() if freq_count > 0))
        freq_counts, prot_counts = zip(*counts)

        fig = plt.figure(figsize=(10, 5))
        plt.scatter(freqs, freq_counts, marker='o', c=list(map(lambda c: cmap(norm(c)), prot_counts)))
        for f, f_count in zip(freqs, freq_counts):
            plt.text(
                f,
                f_count - (max(freq_counts)-min(freq_counts)) / 100 * 1.5,
                round(f, 2),
                ha="center",
                va="top"
            )

        freqs_per_win_without_first_quantiles = [len(quantiles) for freq_count, quantiles in zip(freqs_per_win, quantiles_per_win_without_first_ones)]
        freqs_per_win = round(np.mean(freqs_per_win), 2)
        freqs_per_win_without_first_quantiles = round(np.mean(freqs_per_win_without_first_quantiles), 2)

        print(freqs_per_win, freqs_per_win_without_first_quantiles, round(np.concatenate(quantiles_per_win_without_first_ones).mean(), 2), sep=",")
        plt.colorbar(ScalarMappable(cmap=cmap, norm=norm), ax=plt.gca()).set_label("Sequences", rotation=-90, va="bottom")
        plt.xlabel("Frequencies")
        plt.ylabel("Absolute count in all sequences")
        plt.title("Occurences of selected STFT frequencies with average %g frequencies per window" % freqs_per_win)
        plt.savefig(out_file, bbox_inches='tight')


def _process(args) -> Dict[float, Tuple[int, int]]:
    fasta, prot_slice, config = args

    FREQS = sorted(abs(i) for i in fftfreq(config.window_size) if i >= 0 or i == -.5)
    selected_freqs = dict.fromkeys(FREQS, (0, 0))
    freqs_per_win = []
    quantiles_per_win_without_first_ones = []

    for _, _, seq in fasta[prot_slice]:
        constellation_map = create_constellation(get_aa_vector(seq), config)
        all_freqs = set()
        for window in constellation_map:
            freqs_per_win.append(len(window))
            win_quantiles = []
            for freq_idx, _, quantile in window:
                if freq_idx > config.n_peaks - 1:
                    win_quantiles.append(quantile)
                freq = FREQS[freq_idx]
                freq_count, prot_count = selected_freqs[freq]
                selected_freqs[freq] = (freq_count + 1, prot_count)
                all_freqs.add(freq)
            quantiles_per_win_without_first_ones.append(win_quantiles)

        for freq in all_freqs:
            freq_count, prot_count = selected_freqs[freq]
            selected_freqs[freq] = (freq_count, prot_count + 1)

    return selected_freqs, freqs_per_win, quantiles_per_win_without_first_ones
