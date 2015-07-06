from numpy import array, asarray, ceil
from models import LoadSampler

def case75(flex_level = 'MEDIUM'):

    ppc = {"version": '2'}

    ## system MVA base
    ppc["baseMVA"] = 1.0

    ## Bus data
    ppc["bus"] = array([
        [1000, 3, 0.0, 0.0, 0, 0, 1, 1, 0, 33, 1, 1.1, 0.9],
        [1100, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1101, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1102, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1103, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1104, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1105, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1106, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1107, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1108, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1109, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1110, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1111, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1112, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1113, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1114, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1115, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1116, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1117, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1118, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1119, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1120, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1121, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1122, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1123, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1124, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1125, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1126, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1127, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1128, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1129, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1130, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1131, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1132, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1133, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1134, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1135, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1136, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1137, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1138, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1139, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1140, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1141, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1142, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1143, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1144, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1145, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1146, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1147, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1148, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1149, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1150, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1151, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1152, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1153, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1154, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1155, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1156, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1157, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1158, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1159, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1160, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1161, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1162, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1163, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1164, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1165, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1166, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1167, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1168, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1169, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1170, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1171, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1172, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1173, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1174, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95],
        [1175, 1, 0.0, 0.0, 0, 0, 1, 1, 0, 11, 1, 1.05, 0.95]
    ])


    ## Gen data
    ppc["gen"] = array([
        [1000, 0.0, 0.0, 100.0, -100.0, 1.02, 100.0, 1, 100.0, -100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1101, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1103, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1104, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1106, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1107, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1109, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1110, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1111, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1113, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1114, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1115, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1117, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1120, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1122, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1123, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1124, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1125, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1126, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1128, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1131, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1133, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1134, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1135, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1136, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1137, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1139, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1142, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1144, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1146, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1147, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1148, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1149, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1150, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1151, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1153, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1154, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1156, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1157, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1158, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1160, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1161, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1163, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1164, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1165, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1167, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1168, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1169, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1170, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1171, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1172, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1173, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1174, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1175, -1.0, -0.2031, 0.0, 0.0, 0.0, 100.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1101.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1103.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.32, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1104.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1106.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.32, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1107.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1109.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.32, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1110.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1111.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1113.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1114.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.31, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1115.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.27, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1117.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.27, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1120.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.27, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1122.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1123.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1124.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1125.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1126.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.27, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1128.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.27, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1131.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.27, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1133.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1134.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1135.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1136.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1137.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.27, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1139.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.28, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1142.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.28, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1144.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.28, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1146.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.41, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1147.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.41, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1148.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.41, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1149.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.41, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1150.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.41, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1151.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1153.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1154.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1156.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1157.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1158.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1160.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1161.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1163.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1164.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1165.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1167.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.35, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1168.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1169.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.36, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1170.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1171.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.36, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1172.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1173.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.36, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1174.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1175.0, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1.0, 0.36, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1102, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1, 4.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1105, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1, 4.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1108, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1, 4.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1116, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1, 4.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1127, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1, 4.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1152, 1.0, 0.2031, 0.0, 0.0, 0.0, 100.0, 1, 4.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    ])

    ## Branch data
    ppc["branch"] = array([
        [1000.0, 1100.0, 0.0, 0.004667, 0.0, 25.0, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1100.0, 1101.0, 0.002245, 0.001164, 0.0, 6.82, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1101.0, 1102.0, 0.002245, 0.001164, 0.0, 6.82, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1102.0, 1103.0, 0.000688, 0.000187, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1100.0, 1104.0, 0.002245, 0.001164, 0.0, 6.82, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1104.0, 1105.0, 0.002245, 0.001164, 0.0, 6.82, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1105.0, 1106.0, 0.000688, 0.000187, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1100.0, 1107.0, 0.002245, 0.001164, 0.0, 6.82, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1107.0, 1108.0, 0.002245, 0.001164, 0.0, 6.82, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1108.0, 1109.0, 0.000688, 0.000187, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1100.0, 1110.0, 0.002931, 0.001519, 0.0, 6.82, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1110.0, 1111.0, 0.002931, 0.001519, 0.0, 6.82, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1111.0, 1112.0, 0.002931, 0.001519, 0.0, 6.82, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1111.0, 1113.0, 0.000731, 0.000199, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1112.0, 1114.0, 0.000731, 0.000199, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1100.0, 1115.0, 0.000821, 0.000632, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1115.0, 1116.0, 0.000821, 0.000632, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1116.0, 1117.0, 0.000821, 0.000632, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1117.0, 1118.0, 0.000821, 0.000632, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1118.0, 1119.0, 0.000821, 0.000632, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1119.0, 1120.0, 0.000821, 0.000632, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1120.0, 1121.0, 0.000821, 0.000632, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1116.0, 1122.0, 0.000597, 0.000161, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1118.0, 1123.0, 0.000597, 0.000161, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1119.0, 1124.0, 0.000597, 0.000161, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1121.0, 1125.0, 0.000597, 0.000161, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1100.0, 1126.0, 0.000821, 0.000632, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1126.0, 1127.0, 0.000821, 0.000632, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1127.0, 1128.0, 0.000821, 0.000632, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1128.0, 1129.0, 0.000821, 0.000632, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1129.0, 1130.0, 0.000821, 0.000632, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1130.0, 1131.0, 0.000821, 0.000632, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1131.0, 1132.0, 0.000821, 0.000632, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1127.0, 1133.0, 0.000597, 0.000161, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1129.0, 1134.0, 0.000597, 0.000161, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1130.0, 1135.0, 0.000597, 0.000161, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1132.0, 1136.0, 0.000597, 0.000161, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1100.0, 1137.0, 0.001011, 0.000777, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1137.0, 1138.0, 0.001011, 0.000777, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1138.0, 1139.0, 0.001011, 0.000777, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1139.0, 1140.0, 0.001011, 0.000777, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1140.0, 1141.0, 0.001011, 0.000777, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1141.0, 1142.0, 0.001011, 0.000777, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1142.0, 1143.0, 0.001011, 0.000777, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1143.0, 1144.0, 0.001011, 0.000777, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1144.0, 1145.0, 0.001011, 0.000777, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1138.0, 1146.0, 0.000629, 0.000171, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1140.0, 1147.0, 0.000629, 0.000171, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1141.0, 1148.0, 0.000629, 0.000171, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1143.0, 1149.0, 0.000629, 0.000171, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1145.0, 1150.0, 0.000629, 0.000171, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1100.0, 1151.0, 0.000733, 0.000564, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1151.0, 1152.0, 0.000733, 0.000564, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1152.0, 1153.0, 0.000733, 0.000564, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1153.0, 1154.0, 0.000733, 0.000564, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1154.0, 1155.0, 0.000733, 0.000564, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1155.0, 1156.0, 0.000733, 0.000564, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1156.0, 1157.0, 0.000733, 0.000564, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1157.0, 1158.0, 0.000733, 0.000564, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1158.0, 1159.0, 0.000733, 0.000564, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1159.0, 1160.0, 0.000733, 0.000564, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1160.0, 1161.0, 0.000733, 0.000564, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1161.0, 1162.0, 0.000733, 0.000564, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1162.0, 1163.0, 0.000733, 0.000564, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1163.0, 1164.0, 0.000733, 0.000564, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1164.0, 1165.0, 0.000733, 0.000564, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1165.0, 1166.0, 0.000733, 0.000564, 0.0, 8.86, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1152.0, 1167.0, 0.000803, 0.000219, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1154.0, 1168.0, 0.000803, 0.000219, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1155.0, 1169.0, 0.000803, 0.000219, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1157.0, 1170.0, 0.000803, 0.000219, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1159.0, 1171.0, 0.000803, 0.000219, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1161.0, 1172.0, 0.000803, 0.000219, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1162.0, 1173.0, 0.000803, 0.000219, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1164.0, 1174.0, 0.000803, 0.000219, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0],
        [1166.0, 1175.0, 0.000803, 0.000219, 0.0, 4.84, 0.0, 0.0, 0.0, 0.0, 1.0, -360.0, 360.0]
    ])

    ## List of all potential flexible services of loads.
    # Type, device, price, signal
    flex_conso = [
        [1, 1, 1.96, [0.0404, 0.0925, 0.1244, 0.1292, 0.1058, 0.0594, -0.0, -0.0594, -0.1058, -0.1292, -0.1244, -0.0925, -0.0404]],
        [1, 2, 2.54, [-0.0523, -0.0915, -0.1247, -0.1497, -0.1648, -0.1691, -0.1623, -0.1449, -0.1179, -0.0832, -0.043, 0.0, 0.043, 0.0832, 0.1179, 0.1449, 0.1623, 0.1691, 0.1648, 0.1497, 0.1247, 0.0915, 0.0523]],
        [1, 3, 1.96, [0.0404, 0.081, 0.1115, 0.1283, 0.1292, 0.1141, 0.0849, 0.0453, -0.0, -0.0453, -0.0849, -0.1141, -0.1292, -0.1283, -0.1115, -0.081, -0.0404]],
        [1, 4, 2.54, [-0.0523, -0.1304, -0.1679, -0.1531, -0.0907, 0.0, 0.0907, 0.1531, 0.1679, 0.1304, 0.0523]],
        [1, 5, 1.96, [0.0404, 0.0769, 0.1058, 0.1244, 0.1308, 0.1244, 0.1058, 0.0769, 0.0404, -0.0, -0.0404, -0.0769, -0.1058, -0.1244, -0.1308, -0.1244, -0.1058, -0.0769, -0.0404]],
        [1, 6, 2.54, [-0.0523, -0.102, -0.1405, -0.1636, -0.1687, -0.1554, -0.125, -0.081, -0.028, 0.028, 0.081, 0.125, 0.1554, 0.1687, 0.1636, 0.1405, 0.102, 0.0523]],
        [1, 7, 1.97, [0.0406, 0.0812, 0.1119, 0.1287, 0.1296, 0.1145, 0.0852, 0.0454, -0.0, -0.0454, -0.0852, -0.1145, -0.1296, -0.1287, -0.1119, -0.0812, -0.0406]],
        [1, 8, 1.97, [-0.0406, -0.0863, -0.1182, -0.1311, -0.1228, -0.0948, -0.0516, 0.0, 0.0516, 0.0948, 0.1228, 0.1311, 0.1182, 0.0863, 0.0406]],
        [1, 9, 0.5, [0.0103, 0.0235, 0.0316, 0.0328, 0.0269, 0.0151, -0.0, -0.0151, -0.0269, -0.0328, -0.0316, -0.0235, -0.0103]],
        [1, 10, 2.49, [-0.0513, -0.1174, -0.1579, -0.164, -0.1343, -0.0754, 0.0, 0.0754, 0.1343, 0.164, 0.1579, 0.1174, 0.0513]],
        [1, 11, 2.13, [0.0439, 0.1094, 0.1409, 0.1285, 0.0761, -0.0, -0.0761, -0.1285, -0.1409, -0.1094, -0.0439]],
        [1, 12, 2.13, [-0.0439, -0.1046, -0.1383, -0.1362, -0.099, -0.0361, 0.0361, 0.099, 0.1362, 0.1383, 0.1046, 0.0439]],
        [1, 13, 2.13, [0.0439, 0.1409, 0.0761, -0.0761, -0.1409, -0.0439]],
        [1, 14, 3.19, [-0.0658, -0.1283, -0.1767, -0.2057, -0.2122, -0.1955, -0.1573, -0.1018, -0.0352, 0.0352, 0.1018, 0.1573, 0.1955, 0.2122, 0.2057, 0.1767, 0.1283, 0.0658]],
        [1, 15, 3.19, [0.0658, 0.1568, 0.2073, 0.2042, 0.1483, 0.0541, -0.0541, -0.1483, -0.2042, -0.2073, -0.1568, -0.0658]],
        [1, 16, 3.20, [-0.0659, -0.1359, -0.1868, -0.2115, -0.2065, -0.1725, -0.1142, -0.0399, 0.0399, 0.1142, 0.1725, 0.2065, 0.2115, 0.1868, 0.1359, 0.0659]],
        [1, 17, 3.20, [0.0659, 0.1571, 0.2077, 0.2046, 0.1486, 0.0542, -0.0542, -0.1486, -0.2046, -0.2077, -0.1571, -0.0659]],
        [1, 18, 2.13, [-0.0439, -0.0934, -0.1279, -0.1419, -0.1329, -0.1026, -0.0558, 0.0, 0.0558, 0.1026, 0.1329, 0.1419, 0.1279, 0.0934, 0.0439]],
        [1, 19, 2.13, [0.0439, 0.0856, 0.1179, 0.1373, 0.1416, 0.1304, 0.1049, 0.0679, 0.0235, -0.0235, -0.0679, -0.1049, -0.1304, -0.1416, -0.1373, -0.1179, -0.0856, -0.0439]],
        [1, 20, 2.13, [-0.0439, -0.0905, -0.1244, -0.1409, -0.1375, -0.1149, -0.0761, -0.0266, 0.0266, 0.0761, 0.1149, 0.1375, 0.1409, 0.1244, 0.0905, 0.0439]],
        [1, 21, 3.19, [0.0658, 0.1449, 0.1971, 0.2126, 0.1884, 0.1292, 0.0459, -0.0459, -0.1292, -0.1884, -0.2126, -0.1971, -0.1449, -0.0658]],
        [1, 22, 3.19, [-0.0658, -0.1917, -0.1992, -0.0836, 0.0836, 0.1992, 0.1917, 0.0658]],
        [1, 23, 3.20, [0.0659, 0.1452, 0.1975, 0.213, 0.1888, 0.1294, 0.046, -0.046, -0.1294, -0.1888, -0.213, -0.1975, -0.1452, -0.0659]],
        [1, 24, 3.20, [-0.0659, -0.1359, -0.1868, -0.2115, -0.2065, -0.1725, -0.1142, -0.0399, 0.0399, 0.1142, 0.1725, 0.2065, 0.2115, 0.1868, 0.1359, 0.0659]],
        [1, 25, 2.17, [0.0448, 0.1067, 0.141, 0.1389, 0.1009, 0.0368, -0.0368, -0.1009, -0.1389, -0.141, -0.1067, -0.0448]],
        [1, 26, 2.18, [-0.0449, -0.1175, -0.1452, -0.1175, -0.0449, 0.0449, 0.1175, 0.1452, 0.1175, 0.0449]],
        [1, 27, 2.18, [0.0449, 0.0955, 0.1308, 0.1451, 0.1359, 0.1049, 0.0571, -0.0, -0.0571, -0.1049, -0.1359, -0.1451, -0.1308, -0.0955, -0.0449]],
        [1, 28, 2.18, [-0.0449, -0.0785, -0.107, -0.1284, -0.1414, -0.1451, -0.1393, -0.1243, -0.1012, -0.0714, -0.0369, 0.0, 0.0369, 0.0714, 0.1012, 0.1243, 0.1393, 0.1451, 0.1414, 0.1284, 0.107, 0.0785, 0.0449]],
        [1, 29, 3.25, [0.067, 0.1533, 0.2062, 0.2141, 0.1754, 0.0984, -0.0, -0.0984, -0.1754, -0.2141, -0.2062, -0.1533, -0.067]],
        [1, 30, 3.27, [-0.0674, -0.1158, -0.1573, -0.1893, -0.2099, -0.2179, -0.2128, -0.1949, -0.1653, -0.1257, -0.0786, -0.0267, 0.0267, 0.0786, 0.1257, 0.1653, 0.1949, 0.2128, 0.2179, 0.2099, 0.1893, 0.1573, 0.1158, 0.0674]],
        [1, 31, 3.27, [0.0674, 0.135, 0.1859, 0.2138, 0.2153, 0.1902, 0.1416, 0.0755, -0.0, -0.0755, -0.1416, -0.1902, -0.2153, -0.2138, -0.1859, -0.135, -0.0674]],
        [1, 32, 3.27, [-0.0674, -0.1485, -0.2019, -0.2177, -0.193, -0.1324, -0.047, 0.047, 0.1324, 0.193, 0.2177, 0.2019, 0.1485, 0.0674]],
        [1, 33, 3.27, [0.0674, 0.1485, 0.2019, 0.2177, 0.193, 0.1324, 0.047, -0.047, -0.1324, -0.193, -0.2177, -0.2019, -0.1485, -0.0674]],
        [1, 34, 1.71, [-0.0352, -0.0776, -0.1056, -0.1139, -0.1009, -0.0692, -0.0246, 0.0246, 0.0692, 0.1009, 0.1139, 0.1056, 0.0776, 0.0352]],
        [1, 35, 1.72, [0.0355, 0.0755, 0.1034, 0.1147, 0.1075, 0.083, 0.0451, -0.0, -0.0451, -0.083, -0.1075, -0.1147, -0.1034, -0.0755, -0.0355]],
        [1, 36, 1.72, [-0.0355, -0.0979, -0.1134, -0.0746, 0.0, 0.0746, 0.1134, 0.0979, 0.0355]],
        [1, 37, 1.72, [0.0355, 0.0846, 0.1118, 0.1101, 0.08, 0.0292, -0.0292, -0.08, -0.1101, -0.1118, -0.0846, -0.0355]],
        [1, 38, 1.72, [-0.0355, -0.0929, -0.1148, -0.0929, -0.0355, 0.0355, 0.0929, 0.1148, 0.0929, 0.0355]],
        [1, 39, 1.72, [0.0355, 0.0782, 0.1063, 0.1147, 0.1016, 0.0697, 0.0248, -0.0248, -0.0697, -0.1016, -0.1147, -0.1063, -0.0782, -0.0355]],
        [1, 40, 1.72, [-0.0355, -0.0755, -0.1034, -0.1147, -0.1075, -0.083, -0.0451, 0.0, 0.0451, 0.083, 0.1075, 0.1147, 0.1034, 0.0755, 0.0355]],
        [1, 41, 1.72, [0.0355, 0.0711, 0.0979, 0.1126, 0.1134, 0.1002, 0.0746, 0.0397, -0.0, -0.0397, -0.0746, -0.1002, -0.1134, -0.1126, -0.0979, -0.0711, -0.0355]],
        [1, 42, 1.72, [-0.0355, -0.0782, -0.1063, -0.1147, -0.1016, -0.0697, -0.0248, 0.0248, 0.0697, 0.1016, 0.1147, 0.1063, 0.0782, 0.0355]],
        [1, 43, 1.72, [0.0355, 0.0782, 0.1063, 0.1147, 0.1016, 0.0697, 0.0248, -0.0248, -0.0697, -0.1016, -0.1147, -0.1063, -0.0782, -0.0355]],
        [1, 44, 1.72, [-0.0355, -0.1092, -0.0929, 0.0, 0.0929, 0.1092, 0.0355]],
        [1, 45, 2.82, [0.0581, 0.108, 0.1484, 0.1757, 0.1876, 0.183, 0.1623, 0.1273, 0.0812, 0.0279, -0.0279, -0.0812, -0.1273, -0.1623, -0.183, -0.1876, -0.1757, -0.1484, -0.108, -0.0581]],
        [1, 46, 1.11, [-0.0229, -0.0704, -0.0599, 0.0, 0.0599, 0.0704, 0.0229]],
        [1, 47, 2.84, [0.0585, 0.1005, 0.1365, 0.1643, 0.1822, 0.1891, 0.1847, 0.1692, 0.1434, 0.1091, 0.0682, 0.0232, -0.0232, -0.0682, -0.1091, -0.1434, -0.1692, -0.1847, -0.1891, -0.1822, -0.1643, -0.1365, -0.1005, -0.0585]],
        [1, 48, 1.12, [-0.0231, -0.0551, -0.0729, -0.0718, -0.0521, -0.019, 0.019, 0.0521, 0.0718, 0.0729, 0.0551, 0.0231]],
        [1, 49, 2.839, [0.0585, 0.1531, 0.1892, 0.1531, 0.0585, -0.0585, -0.1531, -0.1892, -0.1531, -0.0585]],
        [1, 50, 1.12, [-0.0231, -0.0463, -0.0638, -0.0734, -0.0739, -0.0653, -0.0486, -0.0259, 0.0, 0.0259, 0.0486, 0.0653, 0.0739, 0.0734, 0.0638, 0.0463, 0.0231]],
        [1, 51, 2.839, [0.0585, 0.1087, 0.1493, 0.1768, 0.1888, 0.1842, 0.1633, 0.1281, 0.0817, 0.0281, -0.0281, -0.0817, -0.1281, -0.1633, -0.1842, -0.1888, -0.1768, -0.1493, -0.1087, -0.0585]],
        [1, 52, 1.12, [-0.0231, -0.044, -0.0605, -0.0711, -0.0748, -0.0711, -0.0605, -0.044, -0.0231, 0.0, 0.0231, 0.044, 0.0605, 0.0711, 0.0748, 0.0711, 0.0605, 0.044, 0.0231]],
        [1, 53, 2.839, [0.0585, 0.1171, 0.1613, 0.1856, 0.1869, 0.1651, 0.1229, 0.0655, -0.0, -0.0655, -0.1229, -0.1651, -0.1869, -0.1856, -0.1613, -0.1171, -0.0585]]
    ]

    # Sets of generators and loads from the set of all devices.
    gens = (ppc["gen"][:,1]>0.0).nonzero()[0]
    loads = (ppc["gen"][:,1]<0.0).nonzero()[0]

    # Determine subset of considered flexible services of loads given the flex lvl parameter.
    N_flex = len(flex_conso)
    flex_steps = {
        'NONE':[], # 0
        'LOW':range(0,N_flex,5), # 11
        'MEDIUM':range(0,N_flex,5)+range(1,N_flex,5), # 22
        'HIGH':range(0,N_flex,5)+range(1,N_flex,5)+range(2,N_flex,5), # 33
        'VERY HIGH': range(0,N_flex,5)+range(1,N_flex,5)+range(2,N_flex,5)+range(3,N_flex,5), # 43
        'FULL':range(53)
    }
    ppc["flex"] = []
    for k in flex_steps[flex_level]:
        ppc["flex"] += [flex_conso[k]]

    # Define the curtailment services for wind generators.
    curt_price = array([45.,  38.,  35.,  31.,  30.,  34.,  41.,  48.,  53.,  55.,  56.,
                        57.,  54.,  51.,  47.,  46.,  47.,  52.,  59.,  60.,  54.,  52.,
                        51.,  49.])
    for gen in gens[-6:]:
        ppc["flex"] += [[2, gen, lambda real_time, price_time: curt_price[ceil((price_time-1)/4)]]]

    # Define the stochastic models describing the consumption of loads and the production of generators.
    load_scales = array([0.4285,  0.4505,  0.3325,  0.3995,  0.4155,  0.4395,  0.3585,
                         0.405 ,  0.083 ,  0.443 ,  0.36  ,  0.3905,  0.4515,  0.575 ,
                         0.592 ,  0.5625,  0.4605,  0.3445,  0.3775,  0.357 ,  0.578 ,
                         0.404 ,  0.6275,  0.5865,  0.4015,  0.383 ,  0.4625,  0.395 ,
                         0.584 ,  0.4965,  0.6195,  0.5825,  0.5695,  0.279 ,  0.299 ,
                         0.3145,  0.302 ,  0.26  ,  0.493 ,  0.2825,  0.3235,  0.3045,
                         0.2835,  0.3005,  0.5495,  0.198 ,  0.5225,  0.194 ,  0.506 ,
                         0.17  ,  0.54  ,  0.244 ,  0.443 ])
    pv_scales = 2.e-6*array([159.,  206.,  159.,  206.,  159.,  206.,  160.,  160.,   40.,
                             202.,  173.,  173.,  173.,  259.,  259.,  260.,  260.,  173.,
                             173.,  173.,  259.,  259.,  260.,  260.,  176.,  177.,  177.,
                             177.,  264.,  266.,  266.,  266.,  266.,  139.,  140.,  140.,
                             140.,  140.,  140.,  140.,  140.,  140.,  140.,  140.,  229.,
                              90.,  231.,   91.,  231.,   91.,  231.,   91.,  231.])
    ppc["power"] = []
    for load, scale in zip(loads,load_scales):
        ppc["power"] += [[load, LoadSampler(scale)]] # Dev id, callable
    for gen, scale in zip(gens[:-6],pv_scales):
        ppc["power"] += [[gen, lambda ir, ws: scale*ir]] # Dev id, callable
    for gen in gens[-6:]:
        ppc["power"] += [[gen, lambda ir, ws: (9.0/6.0)*min(0.015*(ws-3.5)**3,3.0)*(ws>=3.5)*(ws<=25)]] # Dev id, callable

    return ppc
