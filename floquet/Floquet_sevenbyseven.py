import argparse

from utils import *
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as LA
import more_itertools
from scipy.optimize import minimize
import matplotlib.gridspec as gridspec

plt.rcParams["font.size"] = 18
plt.rcParams["figure.figsize"] = (12, 16)
#print(plt.rcParams)

def ansatz(qc, qnodes, vs, D1, D2):
    for i in range(D1):
        qc.rx(vs[6*i], qnodes[3])
        qc.rz(vs[6*i+1], qnodes[3])
        qc.rx(vs[6*i+2], qnodes[4])
        qc.rz(vs[6*i+3], qnodes[4])
        qc.rx(vs[6*i+4], qnodes[5])
        qc.rz(vs[6*i+5], qnodes[5])
    qc.cz(qnodes[3], qnodes[4])
    qc.cz(qnodes[4], qnodes[5])
    vs = vs[6*D1:]
    for j in range(D2):
        qc.rx(vs[12*j], qnodes[0])
        qc.rz(vs[12*j+1], qnodes[0])
        qc.rx(vs[12*j+2], qnodes[1])
        qc.rz(vs[12*j+3], qnodes[1])
        qc.rx(vs[12*j+4], qnodes[2])
        qc.rz(vs[12*j+5], qnodes[2])
        qc.rx(vs[12*j+6], qnodes[3])
        qc.rz(vs[12*j+7], qnodes[3])
        qc.rx(vs[12*j+8], qnodes[4])
        qc.rz(vs[12*j+9], qnodes[4])
        qc.rx(vs[12*j+10], qnodes[5])
        qc.rz(vs[12*j+11], qnodes[5])
    qc.cz(qnodes[0], qnodes[1])
    qc.cz(qnodes[1], qnodes[2])
    qc.cz(qnodes[2], qnodes[3])
    qc.cz(qnodes[3], qnodes[4])
    qc.cz(qnodes[4], qnodes[5])
    vs2 = vs[12*D2:]
    qc.rx(vs[0], qnodes[0])
    qc.rz(vs[1], qnodes[0])
    qc.rx(vs[2], qnodes[1])
    qc.rz(vs[3], qnodes[1])
    qc.rx(vs[4], qnodes[2])
    qc.rz(vs[5], qnodes[2])
    qc.rx(vs[6], qnodes[3])
    qc.rz(vs[7], qnodes[3])
    qc.rx(vs[8], qnodes[4])
    qc.rz(vs[9], qnodes[4])
    qc.rx(vs[10], qnodes[5])
    qc.rz(vs[11], qnodes[5])

def x(vs, initial_state, D1, D2):
    qnodes = QuantumRegister(6,'qc')
    cnodes = ClassicalRegister(6,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    ansatz(qc, qnodes, vs, D1, D2)
    qc.h(qnodes[5])
    qc.measure(qnodes, cnodes)
    return qc

def all_z(vs, initial_state, D1, D2):
    qnodes = QuantumRegister(6,'qc')
    cnodes = ClassicalRegister(6,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    ansatz(qc, qnodes, vs, D1, D2)
    qc.measure(qnodes, cnodes)
    return qc

def xx(vs, initial_state, D1, D2):
    qnodes = QuantumRegister(6,'qc')
    cnodes = ClassicalRegister(6,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    ansatz(qc, qnodes, vs, D1, D2)
    qc.h(qnodes[4])
    qc.h(qnodes[5])
    qc.measure(qnodes, cnodes)
    return qc

def yy(vs, initial_state, D1, D2):
    qnodes = QuantumRegister(6,'qc')
    cnodes = ClassicalRegister(6,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    ansatz(qc, qnodes, vs, D1, D2)
    qc.rx(np.pi/2, qnodes[4])
    qc.rx(np.pi/2, qnodes[5])
    qc.measure(qnodes, cnodes)
    return qc

def xxx(vs, initial_state, D1, D2):
    qnodes = QuantumRegister(6,'qc')
    cnodes = ClassicalRegister(6,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    ansatz(qc, qnodes, vs, D1, D2)
    qc.h(qnodes[3])
    qc.h(qnodes[4])
    qc.h(qnodes[5])
    qc.measure(qnodes, cnodes)
    return qc

def xyy(vs, initial_state, D1, D2):
    qnodes = QuantumRegister(6,'qc')
    cnodes = ClassicalRegister(6,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    ansatz(qc, qnodes, vs, D1, D2)
    qc.h(qnodes[3])
    qc.rx(np.pi/2, qnodes[4])
    qc.rx(np.pi/2, qnodes[5])
    qc.measure(qnodes, cnodes)
    return qc

def yyx(vs, initial_state, D1, D2):
    qnodes = QuantumRegister(6,'qc')
    cnodes = ClassicalRegister(6,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    ansatz(qc, qnodes, vs, D1, D2)
    qc.ry(np.pi/2, qnodes[3])
    qc.ry(np.pi/2, qnodes[4])
    qc.h(qnodes[5])
    qc.measure(qnodes, cnodes)
    return qc

def yxy(vs, initial_state, D1, D2):
    qnodes = QuantumRegister(6,'qc')
    cnodes = ClassicalRegister(6,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    ansatz(qc, qnodes, vs, D1, D2)
    qc.rx(np.pi/2, qnodes[3])
    qc.h(qnodes[4])
    qc.rx(np.pi/2, qnodes[5])
    qc.measure(qnodes, cnodes)
    return qc

def get_exp_x_000_001(res_x, shots):
    comb = more_itertools.powerset([0,1,2,3,4])
    exps = []
    for ind in comb:
        value = 0
        coef = 1/32
        for k,v in res_x.items():
            target_key = k[0]
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def get_exp_x_011_010(res_x, shots):
    comb = more_itertools.powerset([0,1,2,3,4])
    exps = []
    for ind in comb:
        value = 0
        if 4 in ind:
            coef = -1/32
        else:
            coef = 1/32
        for k,v in res_x.items():
            target_key = k[0]
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def get_exp_x_100_101(res_x, shots):
    comb = more_itertools.powerset([0,1,2,3,4])
    exps = []
    for ind in comb:
        value = 0
        if 3 in ind:
            coef = -1/32
        else:
            coef = 1/32
        for k,v in res_x.items():
            target_key = k[0]
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def get_exp_xx_001_010(res_xx, shots):
    comb = more_itertools.powerset([0,1,2,3])
    exps = []
    for ind in comb:
        value = 0
        coef = 1/32
        for k,v in res_xx.items():
            target_key = k[0]
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def get_exp_yy_001_010(res_yy, shots):
    comb = more_itertools.powerset([0,1,2,3])
    exps = []
    for ind in comb:
        value = 0
        coef = 1/32
        for k,v in res_yy.items():
            target_key = k[0]
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def get_exp_xx_101_110(res_xx, shots):
    comb = more_itertools.powerset([0,1,2,3])
    exps = []
    for ind in comb:
        value = 0
        if 3 in ind:
            coef = -1/32
        else:
            coef = 1/32
        for k,v in res_xx.items():
            target_key = k[0]
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def get_exp_yy_101_110(res_yy, shots):
    comb = more_itertools.powerset([0,1,2,3])
    exps = []
    for ind in comb:
        value = 0
        if 3 in ind:
            coef = -1/32
        else:
            coef = 1/32
        for k,v in res_yy.items():
            target_key = k[0]
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def get_exp_xxx_011_100(res_xxx, shots):
    comb = more_itertools.powerset([0,1,2])
    exps = []
    for ind in comb:
        value = 0
        coef = 1/32
        for k,v in res_xxx.items():
            target_key = k[0]
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def get_exp_yxy_011_100(res_yxy, shots):
    comb = more_itertools.powerset([0,1,2])
    exps = []
    for ind in comb:
        value = 0
        coef = 1/32
        for k,v in res_yxy.items():
            target_key = k[0]
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def get_exp_xyy_011_100(res_xyy, shots):
    comb = more_itertools.powerset([0,1,2])
    exps = []
    for ind in comb:
        value = 0
        coef = -1/32
        for k,v in res_xyy.items():
            target_key = k[0]
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def get_exp_yyx_011_100(res_yyx, shots):
    comb = more_itertools.powerset([0,1,2])
    exps = []
    for ind in comb:
        value = 0
        coef = 1/32
        for k,v in res_yyx.items():
            target_key = k[0]
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def get_exp_000(res, shots):
    comb = more_itertools.powerset([0,1,2,3,4,5])
    exps = []
    for ind in comb:
        value = 0
        coef = 1/64
        for k,v in res.items():
            target_key = ""
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def get_exp_001(res, shots):
    comb = more_itertools.powerset([0,1,2,3,4,5])
    exps = []
    for ind in comb:
        value = 0
        if 5 in ind:
            coef = -1/64
        else:
            coef = 1/64
        for k,v in res.items():
            target_key = ""
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def get_exp_010(res, shots):
    comb = more_itertools.powerset([0,1,2,3,4,5])
    exps = []
    for ind in comb:
        value = 0
        if 4 in ind:
            coef = -1/64
        else:
            coef = 1/64
        for k,v in res.items():
            target_key = ""
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def get_exp_011(res, shots):
    comb = more_itertools.powerset([0,1,2,3,4,5])
    exps = []
    for ind in comb:
        value = 0
        if 4 in ind and 5 in ind:
            coef = 1/64
        elif 4 in ind or 5 in ind:
            coef = -1/64
        else:
            coef = 1/64
        for k,v in res.items():
            target_key = ""
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def get_exp_100(res, shots):
    comb = more_itertools.powerset([0,1,2,3,4,5])
    exps = []
    for ind in comb:
        value = 0
        if 3 in ind:
            coef = -1/64
        else:
            coef = 1/64
        for k,v in res.items():
            target_key = ""
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def get_exp_101(res, shots):
    comb = more_itertools.powerset([0,1,2,3,4,5])
    exps = []
    for ind in comb:
        value = 0
        if 3 in ind and 5 in ind:
            coef = 1/64
        elif 3 in ind or 5 in ind:
            coef = -1/64
        else:
            coef = 1/64
        for k,v in res.items():
            target_key = ""
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def get_exp_110(res, shots):
    comb = more_itertools.powerset([0,1,2,3,4,5])
    exps = []
    for ind in comb:
        value = 0
        if 3 in ind and 4 in ind:
            coef = 1/64
        elif 3 in ind or 4 in ind:
            coef = -1/64
        else:
            coef = 1/64
        for k,v in res.items():
            target_key = ""
            for i in ind:
                target_key += k[5-i]
            if int(target_key.count('1'))%2 == 0:
                value += v
            else:
                value -= v
        exps.append(coef*(value/shots))
    return sum(exps)

def ssvqe_cost_function(thetas, t, width, shots, mu, E, Eu, Eg, dist, D1, D2):
    levels = []
    weights = [10, 9, 8, 7, 6, 5, 4, 3]
    for initial_state in more_itertools.powerset([3,4,5]):
        qc_all_z = all_z(thetas, initial_state, D1, D2)
        qc_x = x(thetas, initial_state, D1, D2)
        qc_xx = xx(thetas, initial_state, D1, D2)
        qc_yy = yy(thetas, initial_state, D1, D2)
        qc_xxx = xxx(thetas, initial_state, D1, D2)
        qc_yxy = yxy(thetas, initial_state, D1, D2)
        qc_xyy = xyy(thetas, initial_state, D1, D2)
        qc_yyx = yyx(thetas, initial_state, D1, D2)
        ret_z = ibmsim(qc_all_z, shots, 'qasm')
        ret_x = ibmsim(qc_x, shots, 'qasm')
        ret_xx = ibmsim(qc_xx, shots, 'qasm')
        ret_yy = ibmsim(qc_yy, shots, 'qasm')
        ret_xxx = ibmsim(qc_xxx, shots, 'qasm')
        ret_yxy = ibmsim(qc_yxy, shots, 'qasm')
        ret_xyy = ibmsim(qc_xyy, shots, 'qasm')
        ret_yyx = ibmsim(qc_yyx, shots, 'qasm')
        offs = mu*E*amplitude(t, width)*get_exp_x_000_001(ret_x, shots)+mu*E*amplitude(t, width)*get_exp_x_011_010(ret_x, shots)+mu*E*amplitude(t, width)*get_exp_x_100_101(ret_x, shots)+mu*E*amplitude(t, width)*get_exp_xx_001_010(ret_xx, shots)+mu*E*amplitude(t, width)*get_exp_xx_101_110(ret_xx, shots)+mu*E*amplitude(t, width)*get_exp_yy_010_011(ret_yy, shots)+mu*E*amplitude(t, width)*get_exp_yy_101_110(ret_yy, shots)+mu*E*amplitude(t, width)*get_exp_xx_101_110(ret_xx, shots)+mu*E*amplitude(t, width)*get_exp_xxx_011_100(ret_xxx, shots)+mu*E*amplitude(t, width)*get_exp_xyy_011_100(ret_xyy, shots)+mu*E*amplitude(t, width)*get_exp_yxy_011_100(ret_yxy, shots)+mu*E*amplitude(t, width)*get_exp_yyx_011_100(ret_yyx, shots)
        diags = (Eu-3)*get_exp_000(ret_z, shots)+(Eg-2)*get_exp_001(ret_z, shots)+(Eu-1)*get_exp_010(ret_z, shots)+Eg*get_exp_011(ret_z, shots)+(Eu+1)*get_exp_100(ret_z,shots)+(Eg+2)*get_exp_101(ret_z,shots)+(Eu+3)*get_exp_110(ret_z,shots)
        levels.append(offs+diags)
    return np.dot(levels, weights)

def ssvqe_optimized_levels(thetas, t, width, shots, mu, E, Eu, Eg, dist, D1, D2):
    levels = []
    for initial_state in more_itertools.powerset([3,4,5]):
        qc_all_z = all_z(thetas, initial_state, D1, D2)
        qc_x = x(thetas, initial_state, D1, D2)
        qc_xx = xx(thetas, initial_state, D1, D2)
        qc_yy = yy(thetas, initial_state, D1, D2)
        qc_xxx = xxx(thetas, initial_state, D1, D2)
        qc_yxy = yxy(thetas, initial_state, D1, D2)
        qc_xyy = xyy(thetas, initial_state, D1, D2)
        qc_yyx = yyx(thetas, initial_state, D1, D2)
        ret_z = ibmsim(qc_all_z, shots, 'qasm')
        ret_x = ibmsim(qc_x, shots, 'qasm')
        ret_xx = ibmsim(qc_xx, shots, 'qasm')
        ret_yy = ibmsim(qc_yy, shots, 'qasm')
        ret_xxx = ibmsim(qc_xxx, shots, 'qasm')
        ret_yxy = ibmsim(qc_yxy, shots, 'qasm')
        ret_xyy = ibmsim(qc_xyy, shots, 'qasm')
        ret_yyx = ibmsim(qc_yyx, shots, 'qasm')
        offs = mu*E*amplitude(t, width)*get_exp_x_000_001(ret_x, shots)+mu*E*amplitude(t, width)*get_exp_x_011_010(ret_x, shots)+mu*E*amplitude(t, width)*get_exp_x_100_101(ret_x, shots)+mu*E*amplitude(t, width)*get_exp_xx_001_010(ret_xx, shots)+mu*E*amplitude(t, width)*get_exp_xx_101_110(ret_xx, shots)+mu*E*amplitude(t, width)*get_exp_yy_010_011(ret_yy, shots)+mu*E*amplitude(t, width)*get_exp_yy_101_110(ret_yy, shots)+mu*E*amplitude(t, width)*get_exp_xx_101_110(ret_xx, shots)+mu*E*amplitude(t, width)*get_exp_xxx_011_100(ret_xxx, shots)+mu*E*amplitude(t, width)*get_exp_xyy_011_100(ret_xyy, shots)+mu*E*amplitude(t, width)*get_exp_yxy_011_100(ret_yxy, shots)+mu*E*amplitude(t, width)*get_exp_yyx_011_100(ret_yyx, shots)
        diags = (Eu-3)*get_exp_000(ret_z, shots)+(Eg-2)*get_exp_001(ret_z, shots)+(Eu-1)*get_exp_010(ret_z, shots)+Eg*get_exp_011(ret_z, shots)+(Eu+1)*get_exp_100(ret_z,shots)+(Eg+2)*get_exp_101(ret_z,shots)+(Eu+3)*get_exp_110(ret_z,shots)
        levels.append(offs+diags)
    return levels

def get_graph(res, exacts, gpath):
    fig = plt.figure(figsize=(12,16))
    gs = gridspec.GridSpec(2,2)

    axLU = plt.subplot(gs[0,0])
    axRU = plt.subplot(gs[0,1])
    axLB = plt.subplot(gs[1,0])
    axRB = plt.subplot(gs[1,1])
    axLU.set_xlim([0,6])
    axLU.set_ylim([-0.5,4])
    axLU.plot(distance, res[0,:,:].T, "k", label="ibm_Kawasaki")
    axLU.plot(distance, exacts[0,:,:].T, "k--", label="classical")
    #axLU.legend()
    axRU.set_xlim([0,6])
    axRU.set_ylim([-0.5,4])
    axRU.plot(distance, res[1,:,:].T, "k", label="ibm_Kawasaki")
    axRU.plot(distance, exacts[1,:,:].T, "k--", label="classical")
    #axRU.legend()
    axLB.set_xlim([0,6])
    axLB.set_ylim([-0.5,4])
    axLB.plot(distance, res[2,:,:].T, "k", label="ibm_Kawasaki")
    axLB.plot(distance, exacts[2,:,:].T, "k--", label="classical")
    #axLB.legend()
    axRB.set_xlim([0,6])
    axRB.set_ylim([-0.5,4])
    axRB.plot(distance, res[3,:,:].T, "k", label="ibm_Kawasaki")
    axRB.plot(distance, exacts[3,:,:].T, "k--", label="classical")
    #axRB.legend()

    plt.savefig(gpath)

def main(distance, ts, mu, E, width, D, alpha, r_0, t_0, n_state, _type, shots, gpath):
    res = np.zeros([len(ts), 8, len(distance)])
    exacts = np.zeros([len(ts), 8, len(distance)])
    init_theta_list = np.random.rand(D1*6+D2*12+12)*1e-1
    bounds = np.tile((-np.pi, np.pi), (D1*6+D2*12+12, 1))
    for i, t in enumerate(ts):
        for m, dist in enumerate(distance):
            Eg = morse_pot(dist, D, alpha, r_0)
            Eu = morse_pot_u(dist, D, alpha, r_0, t_0)
            if _type == "simulator":
                opt = minimize(ssvqe_cost_function, init_theta_list,
                        args=(t, shots, 10000, mu, E, Eu, Eg, dist, D1, D2),
                        method="L-BFGS-B",
                        #tol=1e-3,
                        bounds=bounds)
            #elif _type == "noise_simulator":
            #    opt = minimize(noisesim_cost_function, (0,0),
            #            args=(shots, t, mu, E, width, dist, D, alpha, r_0, t_0),
            #            method="SLSQP",
            #            #tol=1e-3,
            #            bounds=((-np.pi, np.pi), (-np.pi, np.pi)))
            #elif _type == "real":
            #    opt = minimize(real_cost_function, (0,0),
            #            args=(shots, t, mu, E, width, dist, D, alpha, r_0, t_0),
            #            method="SLSQP",
            #            #tol=1e-3,
            #            bounds=((-np.pi, np.pi), (-np.pi, np.pi)))
            else:
                print("please choose a suitable calculation type")
            res[i, :, m] = np.sort(np.array(ssvqe_optimized_levels(opt.x, t, width, shots, mu, E, Eu, Eg, dist, D1, D2)))
            fm = get_quasi_floquet_matrix(t, mu, E, morse_pot(dist, D, alpha, r_0), \
                                               morse_pot_u(dist, D, alpha, r_0, t_0), width, n_state)
            fm[7, -2:] = [0,0]
            exacts[i, :, m] = np.sort(LA.eig(fm)[0])
    np.save("res_floquet_{}".format(_type), res)
    get_graph(res, exacts, gpath)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get Floquet energies using Quantum Computer')
    parser.add_argument('type', type=str, help='an integer for the accumulator')
    args = parser.parse_args()

    mu = 1
    E = 0.8 #5.337e-3
    delta = 1
    omega = 1
    width = 100
    n_states = 3
    ts = [0, 150, 200, 300]
    distance = np.linspace(0.45, 6, 300)
    D = 2.79
    alpha = 0.72/0.7
    r_0 = 2*0.7
    t_0 = -0.15
    shots = 10000
    D1=1
    D2=1
    gpath = "floquet_result_{}.png".format(args.type)

    main(distance, ts, mu, E, width, D, alpha, r_0, t_0, n_states, args.type, shots, gpath)