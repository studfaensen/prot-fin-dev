from typing import Dict, List, Tuple, TextIO
from tools import ProteinID, JSI, Score
from tools import count_appearances_in_file as count
from tqdm import tqdm
import pandas as pd
import re


Matches = Dict[ProteinID, Tuple[JSI, Score]]


def evaluate_protfin(protfin_out_file: str):
    """
    Summarizes the output of prot-fin's find-match output.
    For each sample protein, it collects the following information into csv:
      - sample protein id                         -> Sample_ID
      - number of found matches                   -> First_Match_Count
      - occurs the sample in the found matches    -> Sample_In_First_Matches
      - sequence length of sample                 -> Sequence_Length
      - created hashes for sample                 -> Sample_Hashes
      - Jaccard Similarity Score of sample        -> Sample_JSI
      - Score of sample                           -> Sample_Score
      - Jaccard Similarity Score of found matches -> Sample_JSI
      - Score of found matches                    -> Sample_Score

    The data is printed to stdout

    ...

    Parameters
    ----------
    protfin_out_file : str
        Path to the file the output of protfin was written to
    """
    with open(protfin_out_file) as f:

        # just counting number proteins used as protfin input
        result_count = count("Seems to be", f)

        # data will be collected into a dataframe
        evaluation = pd.DataFrame(columns=["Sample_ID", "First_Match_Count", "Sample_In_First_Matches", "Sequence_Length", "Sample_Hashes", "Sample_JSI", "Sample_Score", "Top_JSI", "Top_Score"])
        for _ in tqdm(range(result_count)):
            matches: Matches = get_match_info(f)

            # now crawl the few lines below a protfin result to collect data
            # about the input protein
            input_sample, seq, jsi, score, hash_count = get_appendix_info(f)

            # insert the data into the dataframe
            evaluation.loc[len(evaluation.index)] = (
                input_sample,  # Sample_ID
                len(matches),  # First_Match_Count
                input_sample in matches,  # Sample_In_First_Matches
                len(seq),  # Sequence_Length
                hash_count,  # Sample_Hashes
                jsi,  # Sample_JSI
                score,  # Sample_Score
                list(matches.values())[0][0] if matches else JSI(0),  # Top_JSI
                list(matches.values())[0][1] if matches else Score(0)  # Top_Score
            )

        # write to stdout
        print(evaluation.to_csv(index=False), end="")


def get_match_info(file: TextIO) -> Matches:
    matches: Matches = {}

    # iterate through the lines of matches and extract their information
    while (match := file.readline()) != "\n":
        prot_id, jacc_index, score = re.findall("(.*) - .*Jaccard Index of (\d+\.\d+) : Score of (\d+)", match)[0]
        matches[prot_id] = (jacc_index, score)

    return matches


def get_appendix_info(file: TextIO) -> Tuple[ProteinID, str, JSI, Score, int]:
    file.readline()  # Seems to be...
    file.readline()  #
    input_sample: ProteinID = re.findall(".*Input: +(.*) - .*", file.readline())[0]  # Input:  ...
    seq: str = file.readline()[:-1]  # sequence line
    jsi, score = re.findall("JSI: (\d+\.\d).*Score: (\d+)", file.readline())[0]  # Input-JSI ...
    file.readline()  #
    hashes: int = re.findall("Found hashes: (\d+)", file.readline())[0]  # Found hashes: ...

    jsi = JSI(jsi)
    score = Score(score)
    hashes = int(hashes)

    return input_sample, seq, jsi, score, hashes
