from numpy import array, ceil
from models import LoadSampler

def case5(flex_level="MEDIUM"):
    case = {"version": "ANM"}
    case["baseMVA"] = 100.0
    # bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
    case["bus"] = array([
        [1, 3, 0, 0, 0, 0, 1, 1, 0, 13.8, 1, 1.1, 0.9],
        [2, 1, 0, 0, 0, 0, 1, 1, 0, 13.8, 1, 1.1, 0.9],
        [3, 1, 0, 0, 0, 0, 1, 1, 0, 13.8, 1, 1.1, 0.9],
        [4, 1, 0, 0, 0, 0, 1, 1, 0, 13.8, 1, 1.1, 0.9],
        [5, 1, 0, 0, 0, 0, 1, 1, 0, 13.8, 1, 1.1, 0.9]
    ])
    # bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin
    case["gen"] = array([
        [1,  0.,     0., 100, -100, 1.04, 100, 1, 100.0, -100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [3, -1., -0.267,   0,    0,  0.0, 100, 1, 100.0, -100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [4, -1., -0.254,   0,    0,  0.0, 100, 1, 100.0, -100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [5, -1., -0.233,   0,    0,  0.0, 100, 1, 100.0, -100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [2,  1.,    0.1,   5,   -5,  0.0, 100, 1, 20.0, 0.0, 7.5, 20.0, -5., 5., -2., 2., 0.0, 0.0, 0.0, 0.0, 0.0]
    ])
    # fbus tbus r x b rateA rateB rateC ratio angle status
    case["branch"] = array([
        [1, 2, 0.035, 0.02, 0., 14, 14, 14, 0, 0, 1, -360.0, 360.0],
        [2, 3, 0.5,    0.4, 0., 14, 14, 14, 0, 0, 1, -360.0, 360.0],
        [3, 4, 0.5,    0.4, 0.,  7,  7,  7, 0, 0, 1, -360.0, 360.0],
        [3, 5, 0.5,    0.4, 0.,  7,  7,  7, 0, 0, 1, -360.0, 360.0]
    ])

    # Sets of loads from the set of all devices.
    loads = (case["gen"][:,1]<0.0).nonzero()[0]

    # Scaling of load models
    load_scales = array([3.6, 3.9, 4.5])

    # Power data definition
    case["power"] = [[4, lambda ir, ws: (20./3.) *min(0.015*(ws-3.5)**3,3.0)*(ws>=3.5)*(ws<=25)]]
    for load, scale in zip(loads,load_scales):
        case["power"] += [[load, LoadSampler(scale)]]

    # Definition of flexibility services

    # Flexible generator
    curt_price = array([45.,  38.,  35.,  31.,  30.,  34.,  41.,  48.,  53.,  55.,  56.,
                        57.,  54.,  51.,  47.,  46.,  47.,  52.,  59.,  60.,  54.,  52.,
                        51.,  49.])
    case["flex"] = [[2, 4, lambda real_time, price_time: curt_price[ceil((price_time-1)/4)]]]
    if flex_level == "NONE":
        return case

    # Flexible loads
    case["flex"] += [[1, 1, 0.5, [0.15, 0.3, 0.15, 0., -0.15, -0.3, -0.15]]]
    if flex_level == "LOW":
        return case

    case["flex"] += [[1, 2, 0.3, [0.15, 0.3, 0.15, 0., -0.15, -0.3, -0.15]]]
    if flex_level == "MEDIUM":
        return case

    case["flex"] += [[1, 3, 0.2, [0.15, 0.3, 0.15, 0., -0.15, -0.3, -0.15]]]
    if flex_level == "HIGH":
        return case

    raise ValueError("Requested flexibility level (i.e. \"%s\") is not available." % flex_level)