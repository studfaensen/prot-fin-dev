from .find_matches import get_filtered_db
from tools import *
import pandas as pd
from tqdm import tqdm


def match_family(
        family_file: str,
        db_in: str,
        filter_quantile=1.0
        ):

    db, hash_blacklist = get_filtered_db(db_in, filter_quantile)
    fam_data = pd.read_csv(family_file, sep=",", index_col="Protein_ID").squeeze()

    fam_hashes = {f: {} for f in fam_data.unique()}
    for hash_, prots in tqdm(db.db.items()):
        fam = fam_data.loc[[prot for _, prot in prots]]
        for f in fam:
            fam_hashes[f][hash_] = fam_hashes[f].get(hash_, 0) + 1

    print("Family_ID", "F_Score", "Precision", "Recall", "Sharpness", "Member_Count", "Match_Count", "Hash_Intersec_Size", sep=",")
    for fam in tqdm(fam_hashes):
        fam_member_count = (fam_data == fam).sum()
        if fam_member_count > 1:
            hashes = set(hash_ for hash_, count in fam_hashes[fam].items() if count == fam_member_count)
            if len(hashes):
                match_prots = {}
                for hash_ in hashes:
                    for _, match_prot in db.db.get(hash_):
                        match_prots[match_prot] = match_prots.get(match_prot, 0) + 1

                match_fams = fam_data.loc[sorted(match_prots, key=lambda x: match_prots[x] / len(hashes))]

                fam_mask = (match_fams == fam).groupby(match_fams.index).max()
                true_pos_cumsum = fam_mask.cumsum()
                true_pos = true_pos_cumsum[-1]
                false_neg = fam_member_count - true_pos
                assert false_neg == 0
                positives = len(match_prots) - true_pos_cumsum.value_counts().get(0, 0)
                false_pos = positives - true_pos

                precision = true_pos / positives
                assert precision >= 0, (true_pos, positives)
                recall = true_pos / (true_pos + false_neg)

                tp_max_score_prot = fam_mask.idxmax()
                fp_max_score_prot = (~fam_mask).idxmax()
                tp_max_score = match_prots[tp_max_score_prot] if fam_mask[tp_max_score_prot] else 0  # true score is divided by len(hashes), but for sharpness this is unnecessary float calculation
                fp_max_score = match_prots[fp_max_score_prot] if not fam_mask[fp_max_score_prot] else 0

                sharpness = (tp_max_score - fp_max_score) / tp_max_score if tp_max_score > 0 else -1

                f1_score = (2 * precision * recall) / (precision + recall)
                print(fam, round(f1_score, 2), round(precision, 2), round(recall, 2), round(sharpness, 2), fam_member_count, len(match_prots), len(hashes), sep=",")

            else:
                eprint("No intersection hashes for", fam)
        else:
            eprint("Ignored family", fam, "with only", fam_member_count, "member")
