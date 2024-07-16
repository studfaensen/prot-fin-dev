from scipy import signal
import numpy as np
from tools import *
from os import environ as env
import pandas as pd
from operator import gt, lt

# parameters for the STFT
WINDOW_SIZE = int(env.get("WINDOW_SIZE", 30))
WINDOW_TYPE = env.get("WINDOW_TYPE", ["boxcar", "triang", "blackman", "hamming", "hann", "bartlett", "flattop", "parzen", "bohman", "blackmanharris", "nuttall", "barthann", "cosine", "exponential", "tukey", "taylor", "lanczos"]\
                                     [0])
OVERLAP = int(env.get("OVERLAP", 15))
N_PEAKS = int(env.get("N_PEAKS", 0))  # 0 means all

SIGNIFICANCE = round(float(env.get("SIGNIFICANCE", 5)), 5)  # for quantile selection in percent


AMPL_QUANTILES = {
    10: {
        5: [ 0.844, 0.54312, 0.0482, 0.1085, 0.04832, 0.025],
        95: [ 1.576, 1.0551, 0.36524, 0.58613, 0.35894, 0.59],
        0.1: [ 0.917, 0.5915, 0.06941, 0.15183, 0.06946, 0.05],
        99.9: [ 2.1, 1.35915, 0.5566, 0.86315, 0.64118, 1.108],
        0.01: [ 0.722, 0.4618, 0.02007, 0.04817, 0.01999, 0.005],
        99.99: [ 2.1, 1.35915, 0.5566, 0.86315, 0.64118, 1.108],
        0.001: [ 0.607, 0.3849, 0.00222, 0.01523, 0.00222, 0.0],
        99.999: [ 2.1, 1.35915, 0.5566, 0.86315, 0.64118, 1.108],
        "std_dev": [ 0.22193, 0.15551, 0.09711, 0.144, 0.09513, 0.17633],
        "mean": [ 1.20153, 0.79177, 0.19102, 0.33865, 0.18884, 0.26475]
    },
    20: {
        5: [ 0.9265, 0.58593, 0.03295, 0.12298, 0.0336, 0.0614, 0.03466, 0.05002, 0.0334, 0.04578, 0.009],
        95: [ 1.47, 0.96363, 0.25633, 0.45731, 0.25985, 0.37825, 0.2688, 0.33932, 0.25623, 0.32337, 0.2875],
        0.1: [ 0.9905, 0.62804, 0.04747, 0.1563, 0.04833, 0.0865, 0.04992, 0.07119, 0.04804, 0.06537, 0.018],
        99.9: [ 2.1, 1.34242, 0.54509, 0.80931, 0.54436, 0.75714, 0.5304, 0.66531, 0.53459, 0.63679, 0.86],
        0.01: [ 0.813, 0.51232, 0.01426, 0.06644, 0.01463, 0.0275, 0.015, 0.02228, 0.01454, 0.02032, 0.0015],
        99.99: [ 2.1, 1.34242, 0.54509, 0.80931, 0.54436, 0.75714, 0.5304, 0.66531, 0.53459, 0.63679, 0.86],
        0.001: [ 0.708, 0.44335, 0.0034, 0.02318, 0.00383, 0.00863, 0.00357, 0.00705, 0.0038, 0.00641, 0.0],
        99.999: [ 2.1, 1.34242, 0.54509, 0.80931, 0.54436, 0.75714, 0.5304, 0.66531, 0.53459, 0.63679, 0.86],
        "std_dev": [ 0.16448, 0.11433, 0.06861, 0.10129, 0.06943, 0.09637, 0.07178, 0.0883, 0.06838, 0.08489, 0.08828],
        "mean": [ 1.20059, 0.77457, 0.13121, 0.28611, 0.13315, 0.20942, 0.13794, 0.18184, 0.13195, 0.17117, 0.1166]
    },
    30: {
        5: [ 0.96233, 0.60741, 0.02684, 0.13735, 0.02703, 0.0618, 0.0273, 0.04408, 0.02884, 0.03842, 0.02779, 0.03539, 0.02723, 0.03403, 0.02677, 0.00933],
        95: [ 1.426, 0.92711, 0.20942, 0.41527, 0.20886, 0.31861, 0.21161, 0.28917, 0.22348, 0.27425, 0.21433, 0.25526, 0.20979, 0.24952, 0.2072, 0.278],
        0.1: [ 1.02267, 0.64754, 0.03862, 0.16597, 0.03885, 0.08498, 0.03924, 0.06262, 0.04148, 0.05489, 0.03994, 0.05059, 0.03912, 0.04869, 0.03847, 0.01867],
        99.9: [ 2.1, 1.33935, 0.49725, 0.74672, 0.47119, 0.65328, 0.48677, 0.57474, 0.47018, 0.57461, 0.56203, 0.56122, 0.50387, 0.54378, 0.54849, 0.94267],
        0.01: [ 0.85467, 0.5364, 0.0117, 0.08587, 0.0118, 0.02858, 0.01196, 0.01967, 0.01259, 0.01704, 0.01217, 0.01571, 0.01193, 0.01508, 0.01171, 0.002],
        99.99: [ 2.1, 1.33935, 0.49725, 0.74672, 0.47119, 0.65328, 0.48677, 0.57474, 0.47018, 0.57461, 0.56203, 0.56122, 0.50387, 0.54378, 0.54849, 0.94267],
        0.001: [ 0.75967, 0.4736, 0.00326, 0.03649, 0.00328, 0.00913, 0.00343, 0.00623, 0.00351, 0.00539, 0.00346, 0.00497, 0.00342, 0.00477, 0.00327, 0.00033],
        99.999: [ 2.1, 1.33935, 0.49725, 0.74672, 0.47119, 0.65328, 0.48677, 0.57474, 0.47018, 0.57461, 0.56203, 0.56122, 0.50387, 0.54378, 0.54849, 0.94267],
        "std_dev": [ 0.14008, 0.09658, 0.05621, 0.08436, 0.05593, 0.07776, 0.05667, 0.07473, 0.05976, 0.07216, 0.05733, 0.06738, 0.05611, 0.06608, 0.05548, 0.08505],
        "mean": [ 1.2002, 0.77018, 0.10673, 0.27363, 0.10688, 0.18493, 0.10806, 0.15685, 0.11443, 0.1444, 0.10973, 0.13396, 0.10745, 0.13011, 0.10593, 0.11537]
    },
    40: {
        5: [ 0.985, 0.62111, 0.02311, 0.14885, 0.0235, 0.06486, 0.02337, 0.04227, 0.0236, 0.03498, 0.02493, 0.03231, 0.02477, 0.02989, 0.02375, 0.02861, 0.02354, 0.02804, 0.02329, 0.02736, 0.00625],
        95: [ 1.40075, 0.90677, 0.18135, 0.39088, 0.18179, 0.29322, 0.18058, 0.25347, 0.18294, 0.23882, 0.19231, 0.23472, 0.1918, 0.22034, 0.18369, 0.21159, 0.18169, 0.20922, 0.18049, 0.20535, 0.20225],
        0.1: [ 1.04125, 0.65952, 0.03324, 0.17437, 0.03377, 0.08684, 0.03358, 0.0595, 0.03389, 0.04985, 0.03582, 0.0462, 0.0356, 0.04279, 0.03414, 0.04095, 0.0338, 0.04015, 0.03349, 0.03919, 0.01275],
        99.9: [ 2.1, 1.33828, 0.47572, 0.72965, 0.42329, 0.57971, 0.44903, 0.55927, 0.51544, 0.5019, 0.49672, 0.51945, 0.46734, 0.57032, 0.4699, 0.49472, 0.52439, 0.50631, 0.44334, 0.55157, 0.86],
        0.01: [ 0.88525, 0.55334, 0.01007, 0.10197, 0.01028, 0.03118, 0.01021, 0.01902, 0.01035, 0.01556, 0.01085, 0.01433, 0.01084, 0.01324, 0.01037, 0.01267, 0.01032, 0.01243, 0.01016, 0.01212, 0.00125],
        99.99: [ 2.1, 1.33828, 0.47572, 0.72965, 0.42329, 0.57971, 0.44903, 0.55927, 0.51544, 0.5019, 0.49672, 0.51945, 0.46734, 0.57032, 0.4699, 0.49472, 0.52439, 0.50631, 0.44334, 0.55157, 0.86],
        0.001: [ 0.796, 0.49532, 0.00274, 0.05267, 0.00287, 0.01009, 0.00277, 0.00602, 0.00297, 0.00493, 0.003, 0.00452, 0.00303, 0.00418, 0.00282, 0.00399, 0.00295, 0.00392, 0.00276, 0.00382, 0.0],
        99.999: [ 2.1, 1.33828, 0.47572, 0.72965, 0.42329, 0.57971, 0.44903, 0.55927, 0.51544, 0.5019, 0.49672, 0.51945, 0.46734, 0.57032, 0.4699, 0.49472, 0.52439, 0.50631, 0.44334, 0.55157, 0.86],
        "std_dev": [ 0.12568, 0.08624, 0.0488, 0.07352, 0.04871, 0.06906, 0.04838, 0.06422, 0.04905, 0.06231, 0.05145, 0.06199, 0.05134, 0.05844, 0.04921, 0.05618, 0.04866, 0.05565, 0.04836, 0.05469, 0.06228],
        "mean": [ 1.2, 0.76827, 0.09204, 0.26809, 0.09289, 0.17561, 0.09228, 0.14141, 0.0933, 0.12782, 0.09845, 0.12268, 0.09802, 0.11466, 0.09388, 0.10993, 0.09288, 0.10828, 0.09217, 0.10606, 0.08188]
    },
    50: {
        5: [ 1.0006, 0.63076, 0.02058, 0.15777, 0.02102, 0.06872, 0.02094, 0.04221, 0.02089, 0.03309, 0.02107, 0.0294, 0.02206, 0.02803, 0.02252, 0.02664, 0.02165, 0.02532, 0.02116, 0.02457, 0.02104, 0.02429, 0.021, 0.02382, 0.02053, 0.0066],
        95: [ 1.3842, 0.89365, 0.16329, 0.37428, 0.16251, 0.27711, 0.16205, 0.23588, 0.16137, 0.21476, 0.16336, 0.20636, 0.17051, 0.20539, 0.17526, 0.19843, 0.16774, 0.18882, 0.16347, 0.18405, 0.16241, 0.18309, 0.16234, 0.18014, 0.15974, 0.2026],
        0.1: [ 1.0536, 0.66751, 0.02961, 0.1811, 0.0302, 0.08959, 0.03007, 0.05888, 0.03, 0.04698, 0.03026, 0.04198, 0.03168, 0.0401, 0.03235, 0.03813, 0.0311, 0.03625, 0.03037, 0.03518, 0.0302, 0.03479, 0.03014, 0.03412, 0.02949, 0.0132],
        99.9: [ 2.1, 1.33778, 0.47351, 0.74305, 0.41786, 0.57513, 0.42881, 0.53669, 0.44679, 0.49059, 0.51544, 0.51766, 0.46139, 0.53103, 0.49072, 0.48614, 0.48685, 0.5944, 0.39577, 0.5187, 0.48607, 0.45298, 0.44729, 0.43524, 0.52556, 0.876],
        0.01: [ 0.909, 0.56699, 0.00901, 0.11422, 0.00923, 0.0348, 0.0092, 0.01918, 0.00918, 0.01477, 0.00928, 0.01306, 0.00969, 0.01243, 0.0099, 0.0118, 0.00952, 0.01122, 0.0093, 0.01088, 0.00927, 0.01075, 0.00924, 0.01055, 0.00902, 0.0014],
        99.99: [ 2.1, 1.33778, 0.47351, 0.74305, 0.41786, 0.57513, 0.42881, 0.53669, 0.44679, 0.49059, 0.51544, 0.51766, 0.46139, 0.53103, 0.49072, 0.48614, 0.48685, 0.5944, 0.39577, 0.5187, 0.48607, 0.45298, 0.44729, 0.43524, 0.52556, 0.876],
        0.001: [ 0.8246, 0.51204, 0.00268, 0.06627, 0.00275, 0.01157, 0.00274, 0.00609, 0.00274, 0.00467, 0.00279, 0.00413, 0.00289, 0.00392, 0.00295, 0.00372, 0.00283, 0.00353, 0.00277, 0.00343, 0.00279, 0.0034, 0.00275, 0.00333, 0.00268, 0.0002],
        99.999: [ 2.1, 1.33778, 0.47351, 0.74305, 0.41786, 0.57513, 0.42881, 0.53669, 0.44679, 0.49059, 0.51544, 0.51766, 0.46139, 0.53103, 0.49072, 0.48614, 0.48685, 0.5944, 0.39577, 0.5187, 0.48607, 0.45298, 0.44729, 0.43524, 0.52556, 0.876],
        "std_dev": [ 0.11586, 0.07931, 0.04411, 0.0658, 0.04356, 0.06303, 0.04344, 0.05877, 0.04325, 0.05543, 0.04382, 0.05419, 0.04568, 0.05439, 0.04696, 0.05274, 0.04497, 0.05025, 0.04381, 0.04901, 0.04354, 0.04883, 0.04351, 0.04807, 0.04288, 0.06224],
        "mean": [ 1.19991, 0.76722, 0.08237, 0.26497, 0.083, 0.17054, 0.0827, 0.13419, 0.08239, 0.11697, 0.08327, 0.10935, 0.08713, 0.10708, 0.0893, 0.10272, 0.08559, 0.0977, 0.08348, 0.09509, 0.08296, 0.09432, 0.08285, 0.09267, 0.08134, 0.08294]
    },
    60: {
        5: [ 1.0125, 0.63852, 0.01874, 0.16428, 0.019, 0.07258, 0.01915, 0.04301, 0.01908, 0.03235, 0.01904, 0.02775, 0.01922, 0.02562, 0.01995, 0.02488, 0.02062, 0.0242, 0.02023, 0.02301, 0.01955, 0.02226, 0.01922, 0.02182, 0.01917, 0.02168, 0.01917, 0.02139, 0.01879, 0.02104, 0.00517],
        95: [ 1.37183, 0.88416, 0.15074, 0.36246, 0.14746, 0.26519, 0.14853, 0.22382, 0.14754, 0.20151, 0.1473, 0.18865, 0.14897, 0.18342, 0.15454, 0.18351, 0.1599, 0.18209, 0.15699, 0.17321, 0.1515, 0.16778, 0.14879, 0.16478, 0.14824, 0.16449, 0.14877, 0.16251, 0.14613, 0.16061, 0.16483],
        0.1: [ 1.06267, 0.67353, 0.02696, 0.18612, 0.02731, 0.09236, 0.0275, 0.05928, 0.0274, 0.04576, 0.02734, 0.03955, 0.02758, 0.03662, 0.02866, 0.0356, 0.02963, 0.03464, 0.02907, 0.03295, 0.02809, 0.03189, 0.02761, 0.03126, 0.02752, 0.03106, 0.02754, 0.03065, 0.027, 0.03017, 0.0105],
        99.9: [ 2.1, 1.33751, 0.44079, 0.67723, 0.39627, 0.53053, 0.39825, 0.48413, 0.41697, 0.50125, 0.48108, 0.47664, 0.47828, 0.45953, 0.41263, 0.51732, 0.40633, 0.4974, 0.41199, 0.43519, 0.56203, 0.42898, 0.44082, 0.4813, 0.48557, 0.44283, 0.42607, 0.42437, 0.43693, 0.52717, 0.825],
        0.01: [ 0.92683, 0.57913, 0.00818, 0.12278, 0.00833, 0.03898, 0.0084, 0.01981, 0.00837, 0.01449, 0.00835, 0.01234, 0.00845, 0.01137, 0.00874, 0.01102, 0.00904, 0.01071, 0.00887, 0.01019, 0.00858, 0.00985, 0.00841, 0.00966, 0.00842, 0.0096, 0.00839, 0.00947, 0.00823, 0.00932, 0.001],
        99.99: [ 2.1, 1.33751, 0.44079, 0.67723, 0.39627, 0.53053, 0.39825, 0.48413, 0.41697, 0.50125, 0.48108, 0.47664, 0.47828, 0.45953, 0.41263, 0.51732, 0.40633, 0.4974, 0.41199, 0.43519, 0.56203, 0.42898, 0.44082, 0.4813, 0.48557, 0.44283, 0.42607, 0.42437, 0.43693, 0.52717, 0.825],
        0.001: [ 0.84667, 0.52688, 0.00231, 0.07551, 0.0024, 0.01354, 0.00241, 0.00631, 0.00241, 0.0046, 0.00236, 0.00389, 0.00249, 0.00359, 0.00248, 0.0035, 0.0026, 0.00339, 0.00255, 0.00321, 0.0025, 0.00311, 0.00239, 0.00305, 0.00249, 0.00303, 0.00238, 0.00299, 0.00236, 0.00294, 0.00017],
        99.999: [ 2.1, 1.33751, 0.44079, 0.67723, 0.39627, 0.53053, 0.39825, 0.48413, 0.41697, 0.50125, 0.48108, 0.47664, 0.47828, 0.45953, 0.41263, 0.51732, 0.40633, 0.4974, 0.41199, 0.43519, 0.56203, 0.42898, 0.44082, 0.4813, 0.48557, 0.44283, 0.42607, 0.42437, 0.43693, 0.52717, 0.825],
        "std_dev": [ 0.10848, 0.07404, 0.04089, 0.06028, 0.03957, 0.05832, 0.03984, 0.05476, 0.03957, 0.05153, 0.0395, 0.04921, 0.03998, 0.04839, 0.04144, 0.0487, 0.04284, 0.04851, 0.04208, 0.04617, 0.04065, 0.04474, 0.0399, 0.04396, 0.03978, 0.04393, 0.03993, 0.04342, 0.03922, 0.04296, 0.05085],
        "mean": [ 1.1998, 0.76656, 0.07549, 0.26288, 0.07518, 0.16717, 0.07572, 0.12963, 0.07529, 0.11118, 0.07515, 0.10117, 0.07591, 0.09647, 0.07884, 0.09541, 0.08161, 0.0938, 0.08007, 0.08928, 0.07726, 0.08644, 0.07593, 0.08485, 0.07567, 0.08451, 0.07583, 0.08344, 0.07443, 0.08233, 0.06677]
    },
    70: {
        5: [ 1.02114, 0.64475, 0.01728, 0.16884, 0.01742, 0.07655, 0.01777, 0.04418, 0.01767, 0.03215, 0.01765, 0.0269, 0.01762, 0.02423, 0.01777, 0.02291, 0.01833, 0.02247, 0.01899, 0.02215, 0.01899, 0.02138, 0.01837, 0.02056, 0.01798, 0.02008, 0.01778, 0.01982, 0.01774, 0.01974, 0.01778, 0.01955, 0.01752, 0.01923, 0.0173, 0.00529],
        95: [ 1.362, 0.87674, 0.14103, 0.35387, 0.13535, 0.2555, 0.13765, 0.2146, 0.13711, 0.19186, 0.1364, 0.17814, 0.13628, 0.16964, 0.13774, 0.1663, 0.14221, 0.16682, 0.14723, 0.16737, 0.14869, 0.162, 0.14257, 0.156, 0.13926, 0.15234, 0.13746, 0.15041, 0.13721, 0.15035, 0.13804, 0.14923, 0.13593, 0.147, 0.13481, 0.16614],
        0.1: [ 1.07, 0.67824, 0.02486, 0.18969, 0.02503, 0.09525, 0.02552, 0.06011, 0.02538, 0.04525, 0.02533, 0.03825, 0.02529, 0.03459, 0.02551, 0.03276, 0.02632, 0.03217, 0.02727, 0.03172, 0.02727, 0.03063, 0.02638, 0.02946, 0.02582, 0.02877, 0.02553, 0.0284, 0.02546, 0.02829, 0.02553, 0.02802, 0.02515, 0.02756, 0.02484, 0.01071],
        99.9: [ 2.1, 1.33735, 0.48007, 0.67338, 0.36076, 0.57199, 0.43599, 0.48666, 0.40007, 0.45565, 0.4564, 0.45292, 0.41361, 0.44242, 0.47828, 0.42646, 0.41665, 0.51105, 0.45718, 0.4197, 0.45269, 0.42405, 0.39823, 0.55711, 0.46476, 0.38358, 0.46457, 0.41452, 0.48494, 0.4302, 0.45439, 0.41766, 0.45318, 0.43214, 0.52539, 0.86143],
        0.01: [ 0.93871, 0.58886, 0.00757, 0.12849, 0.00765, 0.04367, 0.00781, 0.02071, 0.00777, 0.01446, 0.00777, 0.01198, 0.00775, 0.01075, 0.00782, 0.01015, 0.00806, 0.00996, 0.00835, 0.0098, 0.00836, 0.00946, 0.00809, 0.0091, 0.00791, 0.00889, 0.00783, 0.00878, 0.00781, 0.00875, 0.00784, 0.00866, 0.00771, 0.00852, 0.00761, 0.001],
        99.99: [ 2.1, 1.33735, 0.48007, 0.67338, 0.36076, 0.57199, 0.43599, 0.48666, 0.40007, 0.45565, 0.4564, 0.45292, 0.41361, 0.44242, 0.47828, 0.42646, 0.41665, 0.51105, 0.45718, 0.4197, 0.45269, 0.42405, 0.39823, 0.55711, 0.46476, 0.38358, 0.46457, 0.41452, 0.48494, 0.4302, 0.45439, 0.41766, 0.45318, 0.43214, 0.52539, 0.86143],
        0.001: [ 0.86229, 0.53843, 0.00227, 0.08172, 0.0023, 0.0161, 0.00235, 0.00664, 0.00234, 0.00459, 0.00237, 0.00378, 0.00233, 0.00339, 0.00237, 0.00321, 0.00242, 0.00314, 0.00251, 0.0031, 0.00254, 0.00299, 0.00243, 0.00288, 0.00238, 0.00281, 0.00236, 0.00277, 0.00237, 0.00277, 0.00238, 0.00273, 0.00232, 0.00269, 0.00229, 0.00014],
        99.999: [ 2.1, 1.33735, 0.48007, 0.67338, 0.36076, 0.57199, 0.43599, 0.48666, 0.40007, 0.45565, 0.4564, 0.45292, 0.41361, 0.44242, 0.47828, 0.42646, 0.41665, 0.51105, 0.45718, 0.4197, 0.45269, 0.42405, 0.39823, 0.55711, 0.46476, 0.38358, 0.46457, 0.41452, 0.48494, 0.4302, 0.45439, 0.41766, 0.45318, 0.43214, 0.52539, 0.86143],
        "std_dev": [ 0.10285, 0.06991, 0.03841, 0.05631, 0.03634, 0.05423, 0.03693, 0.05156, 0.03679, 0.04857, 0.0366, 0.0462, 0.03656, 0.04455, 0.03699, 0.04402, 0.03816, 0.04434, 0.03948, 0.04463, 0.03994, 0.04324, 0.03825, 0.04169, 0.03737, 0.04069, 0.03686, 0.04017, 0.03685, 0.04019, 0.03707, 0.03991, 0.03648, 0.03933, 0.03623, 0.05115],
        "mean": [ 1.19952, 0.76603, 0.07009, 0.26135, 0.06892, 0.16469, 0.07022, 0.12641, 0.06987, 0.10717, 0.06959, 0.09638, 0.06952, 0.09001, 0.0702, 0.08698, 0.07246, 0.08649, 0.07511, 0.08609, 0.07545, 0.08326, 0.07261, 0.08015, 0.07102, 0.07826, 0.07015, 0.07727, 0.07001, 0.07711, 0.07029, 0.07646, 0.06924, 0.07528, 0.06855, 0.06772]
    },
    80: {
        5: [ 1.02737, 0.64944, 0.01609, 0.17235, 0.01619, 0.08037, 0.01657, 0.04552, 0.01656, 0.03229, 0.0165, 0.02641, 0.01647, 0.0234, 0.01645, 0.0217, 0.01661, 0.02083, 0.01703, 0.02058, 0.01764, 0.02045, 0.01781, 0.02002, 0.0175, 0.01929, 0.01698, 0.01879, 0.01672, 0.01844, 0.01661, 0.01827, 0.01658, 0.01822, 0.01662, 0.01812, 0.01643, 0.01784, 0.01619, 0.01769, 0.0045],
        95: [ 1.35438, 0.87097, 0.13288, 0.34738, 0.12608, 0.2473, 0.12854, 0.20713, 0.12867, 0.18442, 0.1279, 0.17026, 0.12757, 0.16102, 0.12746, 0.15508, 0.12876, 0.15292, 0.13241, 0.15354, 0.13695, 0.15457, 0.13977, 0.15327, 0.13611, 0.14683, 0.13196, 0.143, 0.12966, 0.14046, 0.12845, 0.13916, 0.12836, 0.13922, 0.12922, 0.13867, 0.12762, 0.13672, 0.12605, 0.1359, 0.14263],
        0.1: [ 1.07562, 0.68187, 0.02314, 0.19251, 0.02325, 0.09804, 0.0238, 0.06109, 0.02378, 0.04516, 0.02369, 0.03747, 0.02364, 0.03337, 0.02362, 0.03102, 0.02384, 0.02981, 0.02446, 0.02946, 0.02533, 0.0293, 0.0256, 0.02869, 0.02512, 0.02764, 0.02438, 0.02693, 0.02401, 0.02644, 0.02384, 0.02618, 0.0238, 0.02612, 0.02386, 0.02597, 0.0236, 0.02557, 0.02326, 0.02536, 0.009],
        99.9: [ 2.1, 1.33725, 0.45016, 0.63042, 0.4358, 0.54621, 0.39476, 0.48861, 0.39089, 0.4695, 0.39222, 0.49419, 0.42353, 0.46703, 0.38656, 0.44181, 0.46739, 0.38499, 0.3913, 0.42046, 0.49672, 0.39708, 0.38536, 0.46232, 0.4055, 0.36364, 0.48449, 0.51119, 0.34165, 0.38548, 0.46674, 0.40446, 0.47092, 0.41181, 0.43038, 0.42924, 0.42081, 0.50911, 0.37938, 0.5475, 0.82287],
        0.01: [ 0.94725, 0.59617, 0.00704, 0.13271, 0.0071, 0.04861, 0.00727, 0.02179, 0.00728, 0.01462, 0.00725, 0.01179, 0.00724, 0.01041, 0.00722, 0.00963, 0.00731, 0.00923, 0.00747, 0.00911, 0.00776, 0.00906, 0.00782, 0.00887, 0.00769, 0.00854, 0.00746, 0.00832, 0.00734, 0.00816, 0.00729, 0.00809, 0.00729, 0.00807, 0.0073, 0.00803, 0.00722, 0.0079, 0.0071, 0.00783, 0.00088],
        99.99: [ 2.1, 1.33725, 0.45016, 0.63042, 0.4358, 0.54621, 0.39476, 0.48861, 0.39089, 0.4695, 0.39222, 0.49419, 0.42353, 0.46703, 0.38656, 0.44181, 0.46739, 0.38499, 0.3913, 0.42046, 0.49672, 0.39708, 0.38536, 0.46232, 0.4055, 0.36364, 0.48449, 0.51119, 0.34165, 0.38548, 0.46674, 0.40446, 0.47092, 0.41181, 0.43038, 0.42924, 0.42081, 0.50911, 0.37938, 0.5475, 0.82287],
        0.001: [ 0.87475, 0.54784, 0.00206, 0.0864, 0.00209, 0.01939, 0.00213, 0.00704, 0.00217, 0.00464, 0.00214, 0.00372, 0.00213, 0.00329, 0.00212, 0.00304, 0.00219, 0.00292, 0.00218, 0.00288, 0.0023, 0.00286, 0.00229, 0.0028, 0.00229, 0.00269, 0.00218, 0.00262, 0.00216, 0.00258, 0.00215, 0.00255, 0.00219, 0.00255, 0.00214, 0.00253, 0.00213, 0.0025, 0.00208, 0.00247, 0.00012],
        99.999: [ 2.1, 1.33725, 0.45016, 0.63042, 0.4358, 0.54621, 0.39476, 0.48861, 0.39089, 0.4695, 0.39222, 0.49419, 0.42353, 0.46703, 0.38656, 0.44181, 0.46739, 0.38499, 0.3913, 0.42046, 0.49672, 0.39708, 0.38536, 0.46232, 0.4055, 0.36364, 0.48449, 0.51119, 0.34165, 0.38548, 0.46674, 0.40446, 0.47092, 0.41181, 0.43038, 0.42924, 0.42081, 0.50911, 0.37938, 0.5475, 0.82287],
        "std_dev": [ 0.09853, 0.06669, 0.03633, 0.0533, 0.03388, 0.05063, 0.0345, 0.04888, 0.03454, 0.0462, 0.03432, 0.0439, 0.03424, 0.04212, 0.03421, 0.04091, 0.03458, 0.04057, 0.03556, 0.04086, 0.03676, 0.04123, 0.03756, 0.04101, 0.03654, 0.03924, 0.03549, 0.03825, 0.03479, 0.03755, 0.03446, 0.03721, 0.03449, 0.03725, 0.03472, 0.03712, 0.03426, 0.03661, 0.03385, 0.03643, 0.04404],
        "mean": [ 1.19923, 0.76558, 0.0656, 0.26026, 0.06412, 0.16279, 0.06553, 0.12398, 0.06552, 0.10423, 0.06518, 0.09293, 0.06503, 0.08598, 0.065, 0.08166, 0.0656, 0.07964, 0.0674, 0.07942, 0.06979, 0.07951, 0.07091, 0.07837, 0.06926, 0.07529, 0.06719, 0.07336, 0.06607, 0.07205, 0.06552, 0.07137, 0.06545, 0.0713, 0.06576, 0.07096, 0.06499, 0.06992, 0.06413, 0.06944, 0.05778]
    },
    90: {
        5: [ 1.03189, 0.65308, 0.01509, 0.17548, 0.01525, 0.08371, 0.01551, 0.04708, 0.01567, 0.03259, 0.01559, 0.02616, 0.01555, 0.02285, 0.01551, 0.02096, 0.01551, 0.01981, 0.01565, 0.01921, 0.016, 0.01904, 0.01654, 0.01902, 0.01679, 0.01882, 0.01666, 0.0183, 0.01621, 0.01774, 0.01592, 0.01739, 0.01569, 0.01717, 0.01564, 0.01703, 0.01563, 0.017, 0.01567, 0.01694, 0.01556, 0.01672, 0.01535, 0.01651, 0.01524, 0.00456],
        95: [ 1.348, 0.86621, 0.12572, 0.34209, 0.11924, 0.24054, 0.12052, 0.20078, 0.12149, 0.17829, 0.12099, 0.16394, 0.12034, 0.1545, 0.12019, 0.14786, 0.1202, 0.14354, 0.12132, 0.14214, 0.12442, 0.14272, 0.12836, 0.14401, 0.13115, 0.14451, 0.13114, 0.13993, 0.1262, 0.13556, 0.12363, 0.13282, 0.12184, 0.13097, 0.12103, 0.13004, 0.12102, 0.1302, 0.12184, 0.13007, 0.12067, 0.1283, 0.11918, 0.12703, 0.1188, 0.14389],
        0.1: [ 1.08011, 0.68479, 0.02171, 0.19501, 0.02191, 0.10049, 0.02226, 0.06228, 0.0225, 0.0453, 0.02238, 0.03702, 0.02231, 0.03256, 0.02227, 0.02993, 0.02227, 0.02833, 0.02247, 0.0275, 0.02297, 0.02726, 0.02374, 0.02725, 0.02411, 0.02696, 0.02392, 0.02623, 0.02327, 0.02542, 0.02285, 0.02493, 0.02253, 0.0246, 0.02244, 0.0244, 0.02243, 0.02436, 0.0225, 0.02428, 0.02234, 0.02395, 0.02204, 0.02367, 0.02188, 0.00922],
        99.9: [ 2.1, 1.33717, 0.45253, 0.62473, 0.34438, 0.53689, 0.38276, 0.47881, 0.41504, 0.45709, 0.40318, 0.43915, 0.40399, 0.50801, 0.4181, 0.49432, 0.36656, 0.40223, 0.47122, 0.42522, 0.42015, 0.44103, 0.47387, 0.46495, 0.40014, 0.40319, 0.42792, 0.41605, 0.34219, 0.38321, 0.53986, 0.40391, 0.33984, 0.36552, 0.45424, 0.37131, 0.47991, 0.38843, 0.39676, 0.41815, 0.40629, 0.44786, 0.45162, 0.38871, 0.53148, 0.83622],
        0.01: [ 0.95433, 0.60182, 0.00662, 0.13649, 0.0067, 0.05307, 0.00682, 0.02316, 0.00689, 0.01486, 0.00686, 0.01171, 0.00684, 0.01017, 0.00683, 0.0093, 0.00683, 0.00878, 0.0069, 0.00851, 0.00705, 0.00843, 0.00727, 0.00843, 0.00739, 0.00833, 0.00733, 0.0081, 0.00713, 0.00785, 0.00702, 0.0077, 0.00691, 0.0076, 0.00688, 0.00754, 0.00689, 0.00753, 0.0069, 0.0075, 0.00686, 0.0074, 0.00676, 0.00731, 0.00671, 0.00089],
        99.99: [ 2.1, 1.33717, 0.45253, 0.62473, 0.34438, 0.53689, 0.38276, 0.47881, 0.41504, 0.45709, 0.40318, 0.43915, 0.40399, 0.50801, 0.4181, 0.49432, 0.36656, 0.40223, 0.47122, 0.42522, 0.42015, 0.44103, 0.47387, 0.46495, 0.40014, 0.40319, 0.42792, 0.41605, 0.34219, 0.38321, 0.53986, 0.40391, 0.33984, 0.36552, 0.45424, 0.37131, 0.47991, 0.38843, 0.39676, 0.41815, 0.40629, 0.44786, 0.45162, 0.38871, 0.53148, 0.83622],
        0.001: [ 0.88456, 0.55547, 0.002, 0.09079, 0.00202, 0.02308, 0.00207, 0.00758, 0.00208, 0.00473, 0.00209, 0.00371, 0.00208, 0.00322, 0.00207, 0.00292, 0.00207, 0.00277, 0.00211, 0.00269, 0.00215, 0.00266, 0.0022, 0.00266, 0.00225, 0.00263, 0.00222, 0.00255, 0.00216, 0.00247, 0.00214, 0.00243, 0.00209, 0.00239, 0.00208, 0.00237, 0.0021, 0.00237, 0.00209, 0.00237, 0.00209, 0.00233, 0.00206, 0.0023, 0.00203, 0.00011],
        99.999: [ 2.1, 1.33717, 0.45253, 0.62473, 0.34438, 0.53689, 0.38276, 0.47881, 0.41504, 0.45709, 0.40318, 0.43915, 0.40399, 0.50801, 0.4181, 0.49432, 0.36656, 0.40223, 0.47122, 0.42522, 0.42015, 0.44103, 0.47387, 0.46495, 0.40014, 0.40319, 0.42792, 0.41605, 0.34219, 0.38321, 0.53986, 0.40391, 0.33984, 0.36552, 0.45424, 0.37131, 0.47991, 0.38843, 0.39676, 0.41815, 0.40629, 0.44786, 0.45162, 0.38871, 0.53148, 0.83622],
        "std_dev": [ 0.09498, 0.06409, 0.03446, 0.05078, 0.03208, 0.0476, 0.03236, 0.04648, 0.0326, 0.04418, 0.03248, 0.04198, 0.03229, 0.04027, 0.03226, 0.0389, 0.03226, 0.03799, 0.03259, 0.03779, 0.03342, 0.03803, 0.03446, 0.03845, 0.03523, 0.03869, 0.03529, 0.03745, 0.03388, 0.03633, 0.03322, 0.03554, 0.03271, 0.03503, 0.03247, 0.0348, 0.03254, 0.03486, 0.03273, 0.03486, 0.03239, 0.03437, 0.032, 0.03405, 0.03196, 0.04437],
        "mean": [ 1.19887, 0.76515, 0.06175, 0.25943, 0.06054, 0.16131, 0.06135, 0.12206, 0.06192, 0.10192, 0.06162, 0.09023, 0.06134, 0.08301, 0.06123, 0.07824, 0.06127, 0.07516, 0.0618, 0.07378, 0.06329, 0.07367, 0.06539, 0.07403, 0.06663, 0.07378, 0.06635, 0.07162, 0.06419, 0.06941, 0.06295, 0.06803, 0.06207, 0.06712, 0.06172, 0.06663, 0.06169, 0.0666, 0.06199, 0.06646, 0.06146, 0.06556, 0.06065, 0.06488, 0.06037, 0.05854]
    },
    100: {
        5: [ 1.0351, 0.65587, 0.01423, 0.17838, 0.01452, 0.08649, 0.01457, 0.04882, 0.01487, 0.03302, 0.0148, 0.02606, 0.01475, 0.02248, 0.01474, 0.02042, 0.01471, 0.01914, 0.01472, 0.01832, 0.01485, 0.01789, 0.01514, 0.01777, 0.01559, 0.01781, 0.01591, 0.01772, 0.01591, 0.01741, 0.01562, 0.01691, 0.01524, 0.01653, 0.01501, 0.01627, 0.01489, 0.01612, 0.01482, 0.01601, 0.01482, 0.01598, 0.01487, 0.01597, 0.01478, 0.01578, 0.01458, 0.0156, 0.01442, 0.01554, 0.004],
        95: [ 1.3426, 0.86222, 0.11938, 0.33751, 0.11407, 0.23522, 0.11333, 0.19521, 0.11517, 0.17309, 0.11514, 0.15879, 0.11452, 0.14891, 0.11415, 0.14232, 0.11395, 0.13737, 0.11408, 0.13416, 0.11502, 0.13323, 0.11767, 0.13371, 0.12109, 0.13524, 0.12365, 0.13602, 0.12545, 0.13445, 0.1218, 0.12945, 0.11853, 0.12677, 0.11669, 0.12454, 0.11529, 0.12322, 0.11475, 0.1225, 0.11478, 0.1227, 0.11552, 0.12281, 0.11475, 0.12123, 0.11336, 0.11996, 0.11244, 0.11974, 0.1275],
        0.1: [ 1.0835, 0.68714, 0.02048, 0.19734, 0.02087, 0.10256, 0.02092, 0.06361, 0.02135, 0.0456, 0.02126, 0.03677, 0.02117, 0.03197, 0.02116, 0.02914, 0.02111, 0.02736, 0.02113, 0.02621, 0.02131, 0.02561, 0.02173, 0.02545, 0.02238, 0.02553, 0.02284, 0.02539, 0.02285, 0.02496, 0.02243, 0.02423, 0.02188, 0.02369, 0.02155, 0.02332, 0.02137, 0.0231, 0.02127, 0.02295, 0.02126, 0.02291, 0.02134, 0.02288, 0.02123, 0.02261, 0.02095, 0.02237, 0.02071, 0.02228, 0.0081],
        99.9: [ 2.1, 1.33712, 0.46871, 0.65234, 0.33846, 0.48984, 0.33331, 0.43732, 0.38868, 0.46428, 0.39102, 0.44554, 0.35724, 0.41225, 0.4327, 0.43611, 0.39247, 0.48531, 0.35597, 0.3849, 0.46406, 0.37749, 0.41995, 0.44824, 0.36314, 0.52068, 0.35245, 0.39, 0.39194, 0.44082, 0.39116, 0.37114, 0.37101, 0.51552, 0.44089, 0.3274, 0.33655, 0.43072, 0.41903, 0.33933, 0.47251, 0.36552, 0.37748, 0.468, 0.38346, 0.36143, 0.47198, 0.39357, 0.36987, 0.54482, 0.825],
        0.01: [ 0.9606, 0.60623, 0.00624, 0.14008, 0.00639, 0.0568, 0.00641, 0.02474, 0.00654, 0.01518, 0.00652, 0.0117, 0.00649, 0.01001, 0.00649, 0.00907, 0.00648, 0.00849, 0.00648, 0.00812, 0.00655, 0.00793, 0.00667, 0.00787, 0.00687, 0.00789, 0.007, 0.00785, 0.00701, 0.00771, 0.00688, 0.00749, 0.00671, 0.00731, 0.00661, 0.0072, 0.00656, 0.00713, 0.00652, 0.00709, 0.00653, 0.00707, 0.00654, 0.00707, 0.0065, 0.00699, 0.00641, 0.0069, 0.00633, 0.00688, 0.0008],
        99.99: [ 2.1, 1.33712, 0.46871, 0.65234, 0.33846, 0.48984, 0.33331, 0.43732, 0.38868, 0.46428, 0.39102, 0.44554, 0.35724, 0.41225, 0.4327, 0.43611, 0.39247, 0.48531, 0.35597, 0.3849, 0.46406, 0.37749, 0.41995, 0.44824, 0.36314, 0.52068, 0.35245, 0.39, 0.39194, 0.44082, 0.39116, 0.37114, 0.37101, 0.51552, 0.44089, 0.3274, 0.33655, 0.43072, 0.41903, 0.33933, 0.47251, 0.36552, 0.37748, 0.468, 0.38346, 0.36143, 0.47198, 0.39357, 0.36987, 0.54482, 0.825],
        0.001: [ 0.8937, 0.56169, 0.00188, 0.09537, 0.00193, 0.0264, 0.00193, 0.00821, 0.00197, 0.00484, 0.00198, 0.00371, 0.00196, 0.00316, 0.00196, 0.00286, 0.00196, 0.00269, 0.00195, 0.00256, 0.002, 0.0025, 0.00201, 0.00248, 0.00208, 0.0025, 0.00211, 0.00247, 0.00212, 0.00243, 0.00208, 0.00236, 0.00203, 0.00231, 0.00199, 0.00227, 0.00198, 0.00225, 0.00197, 0.00223, 0.002, 0.00223, 0.00197, 0.00223, 0.00196, 0.0022, 0.00193, 0.00218, 0.00191, 0.00217, 0.0001],
        99.999: [ 2.1, 1.33712, 0.46871, 0.65234, 0.33846, 0.48984, 0.33331, 0.43732, 0.38868, 0.46428, 0.39102, 0.44554, 0.35724, 0.41225, 0.4327, 0.43611, 0.39247, 0.48531, 0.35597, 0.3849, 0.46406, 0.37749, 0.41995, 0.44824, 0.36314, 0.52068, 0.35245, 0.39, 0.39194, 0.44082, 0.39116, 0.37114, 0.37101, 0.51552, 0.44089, 0.3274, 0.33655, 0.43072, 0.41903, 0.33933, 0.47251, 0.36552, 0.37748, 0.468, 0.38346, 0.36143, 0.47198, 0.39357, 0.36987, 0.54482, 0.825],
        "std_dev": [ 0.09205, 0.06196, 0.03278, 0.04853, 0.03072, 0.04517, 0.03045, 0.04428, 0.03091, 0.04243, 0.03092, 0.0404, 0.03075, 0.03863, 0.03065, 0.03735, 0.03058, 0.03628, 0.03063, 0.03558, 0.0309, 0.03547, 0.03162, 0.03566, 0.03251, 0.03615, 0.03319, 0.03642, 0.03378, 0.03607, 0.03273, 0.03464, 0.03192, 0.03398, 0.03136, 0.03334, 0.03094, 0.03298, 0.0308, 0.0328, 0.03088, 0.03287, 0.03102, 0.03295, 0.03081, 0.03249, 0.03045, 0.03215, 0.03023, 0.03214, 0.0394],
        "mean": [ 1.19846, 0.76473, 0.05842, 0.25876, 0.05778, 0.16021, 0.05766, 0.12046, 0.05872, 0.10009, 0.05858, 0.08812, 0.05831, 0.08054, 0.05817, 0.07563, 0.05803, 0.07221, 0.05813, 0.06994, 0.0586, 0.06897, 0.05986, 0.0689, 0.06166, 0.06944, 0.06292, 0.06946, 0.0634, 0.06848, 0.06188, 0.06622, 0.0603, 0.06481, 0.0594, 0.06372, 0.05876, 0.06309, 0.0585, 0.0627, 0.05849, 0.06271, 0.05878, 0.06269, 0.05843, 0.06193, 0.05769, 0.06125, 0.05717, 0.06109, 0.05165]
    }
}


def create_constellation(
        aa_vec: np.ndarray,
        window_size=WINDOW_SIZE,
        n_peaks=N_PEAKS,
        window=WINDOW_TYPE,
        overlap=OVERLAP
        ) -> ConstellationMap:
    """
    The function carries out a windowed fast fourier transformation,
    often called short time FFT (STFT), on the given vector and creates a list
    of frequencies, found by the STFT, and the index of the window they are found

    ...

    Parameter
    ---------
    aa_vec : np.ndarray
        An array of floats representing an amino acid sequence
    window_size : int
        size of a window, related to amino acids -> size 2 means two amino acids
    n_peaks : int, optional
        the number of frequency peaks that are selected from the STFT result,
        if 0 then all are selected (defaults to 0)
    window : str
        The window type for the STFT, read https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.get_window.html#scipy.signal.get_window
        for supported ones (defaults to "boxcar")
    overlap : int
        the overlap between two windows during doing the STFT (defaults to half the window size)

    Returns
    -------
     A ConstellationMap, a list of coordinates, meaning the pairs of
     window index and prominent frequency peak of the window
    """

    # adjust window size and overlap if invalid
    if overlap >= window_size:
        overlap = window_size - 1

    if len(aa_vec) < window_size:
        aa_vec = np.pad(aa_vec, (0, window_size - len(aa_vec)))
        assert len(aa_vec) == window_size

    # executing the STFT
    stft_result = signal.stft(
        aa_vec,
        nperseg=window_size,
        noverlap=overlap,
        window=window
    )

    return stft_to_constellation(*stft_result, n_peaks)


def stft_to_constellation(
        frequencies: np.ndarray,
        window_indexes: np.ndarray,
        stft: np.ndarray,
        n_peaks: int
        ) -> ConstellationMap:
    constellation_map: ConstellationMap = []

    # find and collect the most prominent frequencies from STFT per window
    for amplitudes in stft.T:

        # get rid of complex values to make them comparable
        spectrum: np.ndarray = abs(amplitudes)

        peaks: List[Tuple[int, int]] = find_peaks(spectrum, n_peaks)

        constellation_map.append(tuple((int(freq_idx), float(spectrum[freq_idx]), quantile) for freq_idx, quantile in peaks))

    return constellation_map


def find_peaks(spectrum: np.ndarray, n_peaks: int) -> List[Tuple[int, int]]:
    lower, upper = sorted((SIGNIFICANCE, 100 - SIGNIFICANCE))
    peaks = []

    # statistical normalization of spectrum
    def statistical_norm(amplitudes: np.ndarray):
        return (amplitudes - AMPL_QUANTILES[WINDOW_SIZE]["mean"]) / AMPL_QUANTILES[WINDOW_SIZE]["std_dev"]
    norm_spectrum = statistical_norm(spectrum)

    for quantile, (tail, op) in enumerate(((upper, gt), (lower, lt))):
        tail_idx = np.argwhere(op(spectrum, AMPL_QUANTILES[WINDOW_SIZE][tail])).flatten()

        if not len(tail_idx):
            continue

        norm_learned = statistical_norm(np.array(AMPL_QUANTILES[WINDOW_SIZE][tail]))
        quantile_deviations = abs(norm_spectrum - norm_learned)

        selected_deviations = quantile_deviations[tail_idx]
        peaks += [*zip(selected_deviations.tolist(), tail_idx.tolist(), [quantile] * len(selected_deviations))]

    peaks = sorted(peaks, reverse=True)
    if n_peaks:
        peaks = peaks[:n_peaks]

    return [(idx, q) for _, idx, q in peaks]
