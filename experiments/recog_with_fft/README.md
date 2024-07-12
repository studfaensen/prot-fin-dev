# Recognition with Fast Fourier Transformation
`NOTE: This README and the code are updated synchronously.`

The concept of this experiment is to generate a vector of values from an amino acid sequence, using one factor of [the table of Kidera factors](../../materials/Amino_Acid_Kidera_Factors.csv), for now. These factors represent the 10 most identifying features of amino acids.<br>
With doing a Short Time Fourier Transformation (STFT) on the created vector, it is now possible to find structual features in the sequence to create a specific constellation of them.

Now, instead of looking on text-based likelihood, the goal is to identify a protein by its own sequence from this constellation, so it needs to be unique like a fingerprint, but still similar to related sequences of the same family and less similar to others.

---

## Usage
1. go to `./methods`
2. create a database of reference proteins: `python3 protfin.py create-db <ref-fasta>`
3. find best scored matches for protein sequence samples: `python3 protfin.py find-matches <samples-fasta>`

### Tools
```sh
cd methods

# create a boxplot for protein counts for a hash, grouped by the hash's window distance
python3 evaluation.py plot-prots-per-windist database.pickle plot.png

# select sample proteins from different mapman bins
python3 evaluation.py select-samples mapmanreferencebins.results.txt protein.fa > samples.fa

# summarize protfin output
python3 evaluation.py eval protfin_out.csv > protfin_out.summary.csv

# summarize the *.summary.csv
python3 summary.py *.summary.csv

# generate a plot for counts of calculated hashes
TITLE="Distribution of sequences' hash counts" X_LABEL="Hash counts" \
Rscript raincloud_plot.R normal <(python3 evaluation.py print-hash-counts database.pickle) plot.png

# generate a plot for counts of proteins per hash with a log10 transformation
TITLE="Distribution of proteins per hash" X_LABEL="Protein counts" \
Rscript raincloud_plot_log10.R normal <(python3 evaluation.py print-prots-per-hash database.pickle) plot.png
```

### Unit Testing
To run the unit tests, just run the following:
```sh
cd methods
TQDM_DISABLE=1 python3 test.py
```


## Experiments
<ul>
    <li><b>Current</b>
        <ul>
            <li>
                <details>
                    <summary><code>v0.3-exp-stft_params</code> - Redoing <code>v0.1-exp-stft_params</code> because of a previous bug: <a href="https://github.com/usadellab/prot-fin/blob/v0.3-exp-stft_params/experiments/recog_with_fft">go to branch</a></summary>
                    The creation of the constellation map is based on the STFT.<br>
                    To increase the accuracy of the recognition algorithm, it is very important to optimize the parameters to generate the most effective constellation map for a protein.
                    <br><br>
                    Therefore, window size, overlap and number of selected peaks are passed to <code>prot-fin</code>.<br>
                    As every configuration of parameters needs a custom database, this procedure is done in parallel on a compute cluster.<br>
                    The results of each recognition process are summarized in <a href="https://github.com/usadellab/prot-fin/blob/v0.3-exp-stft_params/experiments/recog_with_fft/results/summary.csv">summary.csv</a>.
                    <br><br>
                    The results don't differ too much in their F1-scoring, but when treating it as significant, higher overlap and window sizes lead to better results.
                </details>
            </li>
            <li>
                <details>
                    <summary><code>v0.3-exp-hashed_amplitudes</code> - Redoing <code>v0.2-exp-hashed_amplitudes</code> because of a previous bug: <a href="https://github.com/usadellab/prot-fin/blob/v0.3-exp-hashed_amplitudes/experiments/recog_with_fft">go to branch</a></summary>
                    The recognition of proteins in protfin is based on hashes.<br>
                    Currently, hashes are created of STFT frequency pairs and the distance between them in the constelllation map.<br>
                    Including the STFT amplitudes could increase the hashes' quality, as they store more information then.
                    <br><br>
                    Therefore, the amplitudes will be included in hash generation as the result of the comparisons between the amplitudes of the frequency pairs that are included in a hash.<br>
                    So the amplitudes take only a few bits, as the full values may lead to overfitting.
                    <br><br>
                    Looking at <a href="https://github.com/usadellab/prot-fin/blob/v0.3-exp-hashed_amplitudes/experiments/recog_with_fft/results/summary.csv">summary.csv</a>, it seems that the amplitude information really improves the accuracy.<br>
                    The more bits for amplitude are used, the better the F1-Score. May be something to keep in mind. Currently, the performance is bad.
                </details>
            </li>
            <li>
                <details>
                    <summary><code>v0.3-exp-peak_selection</code> - Redoing <code>v0.2-exp-peak_selection</code> because of a previous bug: <a href="https://github.com/usadellab/prot-fin/blob/v0.3-exp-peak_selection/experiments/recog_with_fft">go to branch</a></summary>
                    To create hashes from a STFT, only a subset of frequencies is selected to be included in the hash generation.<br>
                    Currently, the <code>scipy.signal.find_peaks</code> function is used to select the local maxima only. But in case of protein sequences instead of music, this may not be that useful, as the maxima's neighbors could be still relevant for identification of familiar proteins.
                    <br><br>
                    Therefore, an alternative way of selection is going to be developed. The current approach is to just sort the frequencies by their amplitudes descending and select the first ones.
                    <br><br>
                    The results don't differ that much from previous selection, as <a href="https://github.com/usadellab/prot-fin/blob/v0.3-exp-peak_selection/experiments/recog_with_fft/results/summary.csv">summary.csv</a> shows.<br>
                    The average match count is lower, but not that much. Further analysises necessary.
                </details>
            </li>
        </ul>
    </li>
    <li><details><summary><b>Previous</b></summary>
        <ul>
            <li>
                <details>
                    <summary><code>v0.1-exp-stft_params</code> - Trying different parameters for the STFT to fit the best: <a href="https://github.com/usadellab/prot-fin/blob/v0.1-exp-stft_params/experiments/recog_with_fft">go to branch</a></summary>
                    The creation of the constellation map is based on the STFT.<br>
                    To increase the accuracy of the recognition algorithm, it is very important to optimize the parameters to generate the most effective constellation map for a protein.
                    <br><br>
                    Therefore, window size, overlap and number of selected peaks are passed to <code>prot-fin</code>.<br>
                    As every configuration of parameters needs a custom database, this procedure is done in parallel on a compute cluster.<br>
                    The results of each recognition process are summarized in <a href="https://github.com/usadellab/prot-fin/blob/v0.1-exp-stft_params/experiments/recog_with_fft/results/stft_param_exp.summary.csv">stft_param_exp.summary.csv</a>.
                    <br><br>
                    It looks like that the maximum overlap (so hop size of 1) is the best option for accuracy.<br>
                    Currently, for window size and selected peaks are further analyses necessary.
                </details>
            </li>
            <li>
                <details>
                    <summary><code>v0.2-exp-hash_analysises</code> - Analyzing the generated hashes: <a href="https://github.com/usadellab/prot-fin/blob/v0.2-exp-hash_analysises/experiments/recog_with_fft">go to branch</a></summary>
                    The recognition of proteins in protfin is based on hashes.<br>
                    To increase the accuracy of the recognition algorithm, a high quantity and quality of hashes is of interest.<br>
                    To understand how to improve both efficiently, is the purpose of this experiment.
                    <br><br>
                    Therefore, hash counts and their components will be analyzed.
                    <br><br>
                    Currently, there are very many unused hashes that are just ignored, as <a href="https://github.com/usadellab/prot-fin/blob/v0.2-exp-hash_analysises/experiments/recog_with_fft/results/potential_hashes.png">potential_hashes.png</a> implies.
                </details>
            </li>
            <li>
                <details>
                    <summary><code>v0.2-exp-hashed_amplitudes</code> - Analyzing the influence of STFT amplitudes included in hashes: <a href="https://github.com/usadellab/prot-fin/blob/v0.2-exp-hashed_amplitudes/experiments/recog_with_fft">go to branch</a></summary>
                    The recognition of proteins in protfin is based on hashes.<br>
                    Currently, hashes are created of STFT frequency pairs and the distance between them in the constelllation map.<br>
                    Including the STFT amplitudes could increase the hashes' quality, as they store more information then.
                    <br><br>
                    Therefore, the amplitudes will be included in hash generation as the result of the comparisons between the amplitudes of the frequency pairs that are included in a hash.<br>
                    So the amplitudes take only a few bits, as the full values may lead to overfitting.
                    <br><br>
                    Using 1 or 2 bits seems to work good enough, as <a href="https://github.com/usadellab/prot-fin/blob/v0.2-exp-hashed_amplitudes/experiments/recog_with_fft/results/summary.csv">summary.csv</a> implies.<br>
                    Currently, only the first rank of matches is analyzed to see if the original match was identified. The other related matches need to be checked on familiarity concerning their mapman bins.
                </details>
            </li>
            <li>
                <details>
                    <summary><code>v0.2-exp-peak_selection</code> - Analyzing the peak selection method in STFT: <a href="https://github.com/usadellab/prot-fin/blob/v0.2-exp-peak_selection/experiments/recog_with_fft">go to branch</a></summary>
                    To create hashes from a STFT, only a subset of frequencies is selected to be included in the hash generation.<br>
                    Currently, the <code>scipy.signal.find_peaks</code> function is used to select the local maxima only. But in case of protein sequences instead of music, this may not be that useful, as the maxima's neighbors could be still relevant for identification of familiar proteins.
                    <br><br>
                    Therefore, an alternative way of selection is going to be developed. The current approach is to just sort the frequencies by their amplitudes descending and select the first ones.
                    <br><br>
                    The difference for 5 selected peaks doesn't seem that big, as <a href="https://github.com/usadellab/prot-fin/blob/v0.2-exp-peak_selection/experiments/recog_with_fft/results/summary.csv">summary.csv</a> shows.<br>
                    The average match count is lower, but not that much. Further analysises necessary.
                </details>
            </li>
        </ul>
    </details></li>
</ul>

## Methods (`methods/*`)
<ul>
    <li>
        <details>
            <summary><code>protfin.py</code> - The tool that is going to be developed</summary>
            <table>
                <th>method</th><th>steps</th>
                <tr>
                    <td>actions.algorithm.kidera:<br><code>get_aa_vector(seq, factor, normalize, file)</code></td>
                    <td>
                        <ul><li>defaults: <code>normalize=True</code>, <code>file="../../../materials/Amino_Acid_Kidera_Factors.csv"</code></li></ul>
                        <ol type="1">
                            <li>normalize values by adding the global table mean if <code>normalize</code> is <code>True</code></li>
                            <li>extend value table with columns for symbols representing multiple amino acids, by forming the mean of the corresponding amino acids' vectors</li>
                            <li>extend value table with columns for non-valued amino acids 'O' and 'U', by treating their value as zero</li>
                            <li>transform the sequence and return it</li>
                        </ol>
                    </td>
                </tr>
                <tr>
                    <td>actions.algorithm.constellation:<br><code>create_constellation(aa_vec, window_size, n_peaks, window, **kwargs)</code></td>
                    <td>
                        <ul><li>defaults: <code>n_peaks=0</code>, <code>window="boxcar"</code>, <code>overlap@kwargs=window_size//2</code></li></ul>
                        <ol type="1">
                            <li>Initialize values: set <code>overlap=window_size-1</code> if it is bigger than window size</li>
                            <li>If input sequence is shorter than window size, return empty map</li>
                            <li>Do a STFT on <code>aa_vec</code> with the given parameters</li>
                            <li>for each STF-transformed window filter the amplitudes by the quantiles calculated in the sampling experiment</li>
                            <li>for each filtered amplitudes, get the n most prominent peaks as set by <code>n_peaks</code> or select all if <code>n_peaks=0</code></li>
                            <li>append all triples of peak (frequency index), its amplitude and quantile as one whole n-tuple to the constellation map, so one n-tuple per window with all its frequencies</li>
                        </ol>
                    </td>
                </tr>
                <tr>
                    <td>actions.algorithm.hash_gen:<br><code>create_hashes(constellation_map, prot_id, kidera_factor)</code></td>
                    <td>
                        <ol type="1">
                            <li>
                                for each frequency and its quantile in each window in the map create combinatorial hashes (anker points) with all upcoming frequencies in the next 2<sup>12</sup> windows:<br>
                                as frequencies use a max. of 5 bits each and the quantiles 1 bit each and the kidera factor 4 bits, the hashes are generated by combining them into a 32-bit int like: <br>
                                <code>(zeros)-(kidera_factor)-(quantile)-(other_quantile)-(index_diff)-(freq_of_other_pair)-(frequency)</code><br>
                                So currently there are 4 unused bits of zeros that can be assigned in further experiments.
                            </li>
                            <li>also, as the frequencies in the last window in the map doesn't have any upcoming frequencies to pair up with, they are combined with a dummy frequency that never exists (2<sup>5</sup>-1)
                            <li>save index and protein id for each hash</li>
                        </ol>
                    </td>
                </tr>
                <tr>
                    <td>actions.find_matches:<br><code>score_prots(hashes, database, protein_lookup)</code></td>
                    <td>
                        <ol type="1">
                            <li>for each hash, collect for each protein its offsets to its occurences in the protein sequence</li>
                            <li>for each protein, calculate its Jaccard Similarity Index (JSI)</li>
                            <li>the offset having the most matching occurences and the JSI form the score for a protein, as it is the best fitting constellation of the hashes</li>
                            <li>return the scores as Dictionary of protein identifiers pointing to their scores</li>
                        </ol>
                    </td>
                </tr>
                <tr>
                    <td>actions.create_db:<br><code>create_db(prot_file, db_out)</code></td>
                    <td>
                        <ol type="1">
                            <li>create a database for all proteins in the file by joining the results of <code>create_hashes</code></li>
                            <li>create a protein-lookup as well to get to the hash count for each protein</li>
                            <li>dump both into <code>db_out</code></li>
                        </ol>
                    </td>
                </tr>
                <tr>
                    <td>actions.find_matches:<br><code>find_matches(family_file, db_in, filter_quantile)</code></td>
                    <td>
                        <ol type="1">
                            <li>filter the database hashes by <code>filter_quantile</code></li>
                            <li>for each protein in the file, find all match(es), using the database in <code>db_in</code>, and print them to stdout. The score consists of the custom score multiplied with the JSI</li>
                        </ol>
                    </td>
                </tr>
                <tr>
                    <td>actions.match_family:<br><code>match_family(fasta_file, db_in, filter_quantile)</code></td>
                    <td>
                        <ol type="1">
                            <li><code>family_file</code> is csv with header: <code>Family_ID,Protein_ID</code></li>
                            <li>filter the database hashes by <code>filter_quantile</code></li>
                            <li>for each hash, count how many proteins of a family share this hash</li>
                            <li>for each family, take all hashes shared by all its members and look for matches in database</li>
                            <li>calculate the F-Score for the result and print everything as csv</li>
                        </ol>
                    </td>
                </tr>
            </table>
            <h3>Convenience</h3>
            <code>actions.algorithm.hashes_from_seq(seq, prot_id)</code>
            <ul>
                <li>just the workflow <code>seq_to_vectors</code> $\rightarrow$ <code>create_constellation</code> $\rightarrow$ <code>create_hashes</code> for all kidera factors</li>
            </ul>
            <code>tools.Fasta(fasta_file)</code>
            <ul>
                <li>a class to iterate easily through the fasta file's contents with support of slicing, adding also a progress bar to indicate processed proteins</li>
                <li>currently not validating the file</li>
            </ul>
            <code>tools.count_appearances_in_file(pattern, file)</code>
            <ul>
                <li>used to count fastly e.g. the number of proteins in a file, which is necessary to create an appropriate progress bar</li>
            </ul>
            <code>tools.verify_type(val, ty)</code>
            <ul>
                <li>used in unit tests to easily and deeply verify a value's data type</li>
            </ul>
            <code>tools.pd_read_chunkwise(csv_file, chunksize)</code>
            <ul>
                <li>used for chunkwise iteration over the protfin output csv to reduce memory usage</li>
                <li>a returned item stores all matches of one input protein</li>
            </ul>
        </details>
    </li>
    <li>
        <details>
            <summary><code>evaluation.py</code> - Test prot-fin with training data and evaluate the results</summary>
            <table>
                <th>method</th><th>steps</th>
                <tr>
                    <td><code>evaluate_protfin(protfin_out_file)</code></td>
                    <td>
                        <ol type="1">
                            <li>for each output in <code>protfin_out_file</code>, extract the matches' data and count them</li>
                            <li>collect the input specific data from below the output</li>
                            <li>store everything into a dataframe and write it as csv to stdout</li>
                        </ol>
                    </td>
                </tr>
                <tr>
                    <td><code>select_samples(mapman, protein_file, samples_per_family)</code></td>
                    <td>
                        <ol type="1">
                            <li>identify the protein families in <code>mapman</code> file</li>
                            <li>for each family, select randomly <code>samples_per_family</code> proteins</li>
                            <li>find the selected proteins in <code>protein_file</code> and write them as new FASTA formatted output to stdout</li>
                        </ol>
                    </td>
                </tr>
                <tr>
                    <td><code>print_hash_counts(database)</code></td>
                    <td>
                        <ol type="1">
                            <li>Extract the hash counts from the protein lookup in <code>database</code></li>
                            <li>Print the extracted values comma separated to stdout</li>
                        </ol>
                    </td>
                </tr>
                <tr>
                    <td><code>print_prots_per_hash(database)</code></td>
                    <td>
                        <ol type="1">
                            <li>Extract the counts of proteins per hash from the <code>database</code></li>
                            <li>Print the extracted values comma separated to stdout</li>
                        </ol>
                    </td>
                </tr>
                <tr>
                    <td><code>plot_frequencies(prot_file, out_file, cpu_count)</code></td>
                    <td>
                        <ol type="1">
                            <li>Create the constellation maps of all sequences and collect the selected frequencies</li>
                            <li>Plot the frequences' rates and indicate how many sequences share a frequence</li>
                        </ol>
                    </td>
                </tr>
                <tr>
                    <td><code>plot_prots_per_windist(database, out_file)</code></td>
                    <td>
                        <ol type="1">
                            <li>Collect the protein counts per hash, grouped by the hash's window distance</li>
                            <li>Plot boxes per window distance</li>
                        </ol>
                    </td>
                </tr>
            </table>
        </details>
    </li>
    <li><code>raincloud_plot.R</code> - A script to plot groups of values into a raincloud plot</li>
    <li><code>summary.py</code> - A script to summarize the output of <code>evaluation.py eval</code></li>
    <li><code>eval_times.py</code> - A script to fetch the durations for database creation</li>
</ul>

---
## Results (`results/*`)
|                          file                            |     content
|----------------------------------------------------------|------------------

### Reproduce
In this repository, `protein.fa` is used to generate the database. You can extract the file from [this archive](https://github.com/usadellab/prot-fin/raw/5be77c4247327e3958c89200c03a938ec4734834/material/Mapman_reference_DB_202310.tar.bz2). The archive also includes `mapmanreferencebins.results.txt` which maps the proteins to their families.

The used table of Kidera factors is located in [../../materials/Amino_Acid_Kidera_Factors.csv](../../materials/Amino_Acid_Kidera_Factors.csv) and was generated by [this R-script](https://github.com/usadellab/prot-fin/blob/5be77c4247327e3958c89200c03a938ec4734834/methods/Amino_Acid_Kidera_Factors.R). (Read the [root-README](../../README.md) for more details).

---
## Discussion/Brainstorming

---
## Environment
<ul>
    <li><b>Personal</b><br>
        System: <code>Ubuntu 20.04.6 LTS</code>
        Shell: <code>zsh 5.8</code><br>
        <br>
        <table>
            <th>dependency</th><th>version</th>
            <tr><td>python3</td><td>3.8.10</td></tr>
            <tr><td>scipy</td><td>1.10.1</td></tr>
            <tr><td>numpy</td><td>1.23.0</td></tr>
            <tr><td>pandas</td><td>2.0.1</td></tr>
            <tr><td>tqdm</td><td>4.66.2</td></tr>
            <tr><td>matplotlib</td><td>3.5.2</td></tr>
            <tr><td>R</td><td>3.6.3</td></tr>
            <tr><td>tibble</td><td>3.2.1</td></tr>
            <tr><td>ggplot2</td><td>3.5.0</td></tr>
            <tr><td>ggdist</td><td>3.3.2</td></tr>
            <tr><td>dplyr</td><td>1.1.4</td></tr>
        </table>
    </li>
    <li><b>Slurm</b><br>
        System: <code>Ubuntu 22.04.4 LTS</code>
        Shell: <code>zsh 5.8.1</code><br>
        <br>
        <table>
            <th>dependency</th><th>version</th>
            <tr><td>slurm-wlm</td><td>21.08.5</td></tr>
            <tr><td>python3</td><td>3.10.12</td></tr>
            <tr><td>scipy</td><td>1.11.4</td></tr>
            <tr><td>numpy</td><td>1.26.2</td></tr>
            <tr><td>pandas</td><td>2.1.3</td></tr>
            <tr><td>tqdm</td><td>4.66.2</td></tr>
            <tr><td>matplotlib</td><td>3.8.2</td></tr>
            <tr><td>R</td><td>4.4.0</td></tr>
            <tr><td>tibble</td><td>3.2.1</td></tr>
            <tr><td>ggplot2</td><td>3.5.0</td></tr>
            <tr><td>ggdist</td><td>3.3.2</td></tr>
            <tr><td>dplyr</td><td>1.1.4</td></tr>
        </table>
    </li>
</ul>
