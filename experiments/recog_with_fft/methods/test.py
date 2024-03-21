from protfin import *
from tools import *
import numpy as np


def tests():  # TODO
    def test_get_aa_vector():
        kidera = pd.read_csv("../../../materials/Amino_Acid_Kidera_Factors.csv")

        for factor in range(len(kidera)):
            factor_name = kidera["Kidera_Factor_Description"][factor]

            seq = "".join(kidera.loc[:, "A":].keys())
            expected = kidera.loc[factor, "A":].to_numpy().astype(np.float64)
            returned = get_aa_vector(seq, factor, normalize=False)[:len(expected)]

            # test unnormalized values
            assert np.allclose(expected, returned),\
                f"test_get_aa_vector: Unnormalized feature vector for Kidera factor {factor} ({factor_name}) not equal\n" +\
                 "    orig: {expected}" +\
                 "    got:  {returned}"

            # test normalization
            vecs = get_aa_vector(seq, factor, normalize=False)[:len(expected)]
            norm_diff = vecs - expected.astype(np.float32)
            assert np.allclose(norm_diff, norm_diff[0]),\
                f"test_get_aa_vector: Values for Kidera factor {factor} ({factor_name}) aren't normalized equally:\n" +\
                f"    {norm_diff}"

    return locals().values()


def main():
    for test in tests():
        try:
            test()
            eprint("%s: passed" % test.__name__)
        except Exception as err:
            eprint(f"{err}")


if __name__ == '__main__':
    main()