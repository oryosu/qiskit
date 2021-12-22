from utils import *
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as LA
from scipy.optimize import minimize
import matplotlib.gridspec as gridspec

plt.rcParams["font.size"] = 18
plt.rcParams["figure.figsize"] = (12, 16)
#print(plt.rcParams)

def iix(vs, initial_state):
    qnodes = QuantumRegister(3,'qc')
    cnodes = ClassicalRegister(3,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    qc.cx(qnodes[0], qnodes[1])
    qc.ry(-vs[0], qnodes[0])
    qc.cx(qnodes[1], qnodes[0])
    qc.ry(vs[0], qnodes[0])
    qc.cx(qnodes[0],qnodes[1])
    qc.cx(qnodes[1], qnodes[2])
    qc.ry(-vs[1], qnodes[1])
    qc.cx(qnodes[2], qnodes[1])
    qc.ry(vs[1], qnodes[1])
    qc.cx(qnodes[1],qnodes[2])
    qc.h(qnodes[2])
    qc.measure(qnodes, cnodes)
    return qc

def zzx(vs, initial_state):
    qnodes = QuantumRegister(3,'qc')
    cnodes = ClassicalRegister(3,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    qc.cx(qnodes[0], qnodes[1])
    qc.ry(-vs[0], qnodes[0])
    qc.cx(qnodes[1], qnodes[0])
    qc.ry(vs[0], qnodes[0])
    qc.cx(qnodes[0],qnodes[1])
    qc.cx(qnodes[1], qnodes[2])
    qc.ry(-vs[1], qnodes[1])
    qc.cx(qnodes[2], qnodes[1])
    qc.ry(vs[1], qnodes[1])
    qc.cx(qnodes[1],qnodes[2])
    qc.h(qnodes[2])
    qc.measure(qnodes, cnodes)
    return qc

def zix(vs, initial_state):
    qnodes = QuantumRegister(3,'qc')
    cnodes = ClassicalRegister(3,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    qc.cx(qnodes[0], qnodes[1])
    qc.ry(-vs[0], qnodes[0])
    qc.cx(qnodes[1], qnodes[0])
    qc.ry(vs[0], qnodes[0])
    qc.cx(qnodes[0],qnodes[1])
    qc.cx(qnodes[1], qnodes[2])
    qc.ry(-vs[1], qnodes[1])
    qc.cx(qnodes[2], qnodes[1])
    qc.ry(vs[1], qnodes[1])
    qc.cx(qnodes[1],qnodes[2])
    qc.h(qnodes[2])
    qc.measure(qnodes, cnodes)
    return qc

def izx(vs, initial_state):
    qnodes = QuantumRegister(3,'qc')
    cnodes = ClassicalRegister(3,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    qc.cx(qnodes[0], qnodes[1])
    qc.ry(-vs[0], qnodes[0])
    qc.cx(qnodes[1], qnodes[0])
    qc.ry(vs[0], qnodes[0])
    qc.cx(qnodes[0],qnodes[1])
    qc.cx(qnodes[1], qnodes[2])
    qc.ry(-vs[1], qnodes[1])
    qc.cx(qnodes[2], qnodes[1])
    qc.ry(vs[1], qnodes[1])
    qc.cx(qnodes[1],qnodes[2])
    qc.h(qnodes[2])
    qc.measure(qnodes, cnodes)
    return qc

def ixx(vs, initial_state):
    qnodes = QuantumRegister(3,'qc')
    cnodes = ClassicalRegister(3,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    qc.cx(qnodes[0], qnodes[1])
    qc.ry(-vs[0], qnodes[0])
    qc.cx(qnodes[1], qnodes[0])
    qc.ry(vs[0], qnodes[0])
    qc.cx(qnodes[0],qnodes[1])
    qc.cx(qnodes[1], qnodes[2])
    qc.ry(-vs[1], qnodes[1])
    qc.cx(qnodes[2], qnodes[1])
    qc.ry(vs[1], qnodes[1])
    qc.cx(qnodes[1],qnodes[2])
    qc.h(qnodes[1])
    qc.h(qnodes[2])
    qc.measure(qnodes[1], cnodes[1])
    qc.measure(qnodes[2], cnodes[2])
    return qc

def iyy(vs, initial_state):
    qnodes = QuantumRegister(3,'qc')
    cnodes = ClassicalRegister(3,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    qc.cx(qnodes[0], qnodes[1])
    qc.ry(-vs[0], qnodes[0])
    qc.cx(qnodes[1], qnodes[0])
    qc.ry(vs[0], qnodes[0])
    qc.cx(qnodes[0],qnodes[1])
    qc.cx(qnodes[1], qnodes[2])
    qc.ry(-vs[1], qnodes[1])
    qc.cx(qnodes[2], qnodes[1])
    qc.ry(vs[1], qnodes[1])
    qc.cx(qnodes[1],qnodes[2])
    qc.sdg(qnodes[1])
    qc.sdg(qnodes[2])
    qc.measure(qnodes[1], cnodes[1])
    qc.measure(qnodes[2], cnodes[2])
    return qc

def xxx(vs, initial_state):
    qnodes = QuantumRegister(3,'qc')
    cnodes = ClassicalRegister(3,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    qc.cx(qnodes[0], qnodes[1])
    qc.ry(-vs[0], qnodes[0])
    qc.cx(qnodes[1], qnodes[0])
    qc.ry(vs[0], qnodes[0])
    qc.cx(qnodes[0],qnodes[1])
    qc.cx(qnodes[1], qnodes[2])
    qc.ry(-vs[1], qnodes[1])
    qc.cx(qnodes[2], qnodes[1])
    qc.ry(vs[1], qnodes[1])
    qc.cx(qnodes[1],qnodes[2])
    qc.h(qnodes[0])
    qc.h(qnodes[1])
    qc.h(qnodes[2])
    qc.measure(qnodes, cnodes)
    return qc

def yyx(vs, initial_state):
    qnodes = QuantumRegister(3,'qc')
    cnodes = ClassicalRegister(3,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    qc.cx(qnodes[0], qnodes[1])
    qc.ry(-vs[0], qnodes[0])
    qc.cx(qnodes[1], qnodes[0])
    qc.ry(vs[0], qnodes[0])
    qc.cx(qnodes[0],qnodes[1])
    qc.cx(qnodes[1], qnodes[2])
    qc.ry(-vs[1], qnodes[1])
    qc.cx(qnodes[2], qnodes[1])
    qc.ry(vs[1], qnodes[1])
    qc.cx(qnodes[1],qnodes[2])
    qc.sdg(qnodes[0])
    qc.sdg(qnodes[1])
    qc.h(qnodes[2])
    qc.measure(qnodes, cnodes)
    return qc

def yxy(vs, initial_state):
    qnodes = QuantumRegister(3,'qc')
    cnodes = ClassicalRegister(3,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    qc.cx(qnodes[0], qnodes[1])
    qc.ry(-vs[0], qnodes[0])
    qc.cx(qnodes[1], qnodes[0])
    qc.ry(vs[0], qnodes[0])
    qc.cx(qnodes[0],qnodes[1])
    qc.cx(qnodes[1], qnodes[2])
    qc.ry(-vs[1], qnodes[1])
    qc.cx(qnodes[2], qnodes[1])
    qc.ry(vs[1], qnodes[1])
    qc.cx(qnodes[1],qnodes[2])
    qc.sdg(qnodes[0])
    qc.h(qnodes[1])
    qc.sdg(qnodes[2])
    qc.measure(qnodes, cnodes)
    return qc

def xyy(vs, initial_state):
    qnodes = QuantumRegister(3,'qc')
    cnodes = ClassicalRegister(3,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    qc.cx(qnodes[0], qnodes[1])
    qc.ry(-vs[0], qnodes[0])
    qc.cx(qnodes[1], qnodes[0])
    qc.ry(vs[0], qnodes[0])
    qc.cx(qnodes[0],qnodes[1])
    qc.cx(qnodes[1], qnodes[2])
    qc.ry(-vs[1], qnodes[1])
    qc.cx(qnodes[2], qnodes[1])
    qc.ry(vs[1], qnodes[1])
    qc.cx(qnodes[1],qnodes[2])
    qc.h(qnodes[0])
    qc.sdg(qnodes[1])
    qc.sdg(qnodes[2])
    qc.measure(qnodes, cnodes)
    return qc

def zzz(vs, initial_state):
    qnodes = QuantumRegister(3,'qc')
    cnodes = ClassicalRegister(3,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if initial_state:
        for i in initial_state:
            qc.x(qnodes[i])
    qc.cx(qnodes[0], qnodes[1])
    qc.ry(-vs[0], qnodes[0])
    qc.cx(qnodes[1], qnodes[0])
    qc.ry(vs[0], qnodes[0])
    qc.cx(qnodes[0],qnodes[1])
    qc.cx(qnodes[1], qnodes[2])
    qc.ry(-vs[1], qnodes[1])
    qc.cx(qnodes[2], qnodes[1])
    qc.ry(vs[1], qnodes[1])
    qc.cx(qnodes[1],qnodes[2])
    qc.measure(qnodes, cnodes)
    return qc

def get_all_qcs(thetas):
    qcs = []
    qcs_z = []
    initial_states = [[],[2],[1],[1,2],[0],[0,2],[0,1]]
    for init in initial_states:
        qcs.append(iix(thetas, init))
        qcs.append(zzx(thetas, init))
        qcs.append(zix(thetas, init))
        qcs.append(izx(thetas, init))
        qcs.append(ixx(thetas, init))
        qcs.append(iyy(thetas, init))
        qcs.append(xxx(thetas, init))
        qcs.append(yyx(thetas, init))
        qcs.append(yxy(thetas, init))
        qcs.append(xyy(thetas, init))
        qcs_z.append(zzz(thetas, init))
    return qcs, qcs_z

def res_process_exp(res, shots):
    e = 0
    o = 0
    for k, v in res.items():
        if k.count("1")%2 == 0:
            e += v
        else:
            o += v
    return (e - o)/shots

def res_process_exp_z(res, shots):
    exps = [1]
    #each one qubit
    for i in range(3):
        e = 0
        o = 0
        for k, v in res.items():
            if int(k[i].count("1"))%2 == 0:
                e += v
            else:
                o += v
        exps.append((e-o)/shots)
    #two qubits
    for i in range(3):
        e = 0
        o = 0
        for k, v in res.items():
            if (int(k.count("1"))-int(k[i].count("1")))%2 == 0:
                e += v
            else:
                o += v
        exps.append((e-o)/shots)
    for k, v in res.items():
        if int(k.count("1"))%2 == 0:
            e += v
        else:
            o += v
    exps.append((e-o)/shots)
    return exps

def offdiagonal(res, t, mu, E, width):
    return mu*E*amplitude(t, width)*(0.75*res[0]-0.25*res[1]+0.25*res[2]+0.25*res[3]+0.5*res[4]+0.5*res[5]+0.25*res[6]+0.25*res[7]+0.25*res[8]-0.25*res[9])

def diagonal(res, Eg, Eu):
    ### [III, ZII, IZI, IIZ, iZZ, ZIZ, ZZI, ZZZ]
    coefs = [[1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8],
            [1/8, 1/8, 1/8, -1/8, -1/8, -1/8, 1/8, -1/8],
            [1/8, 1/8, -1/8, 1/8, -1/8, 1/8, -1/8, -1/8],
            [1/8, 1/8, -1/8, -1/8, 1/8, -1/8, -1/8, 1/8],
            [1/8, -1/8, 1/8, 1/8, 1/8, -1/8, -1/8, -1/8],
            [1/8, -1/8, 1/8, -1/8, -1/8, 1/8, -1/8, 1/8],
            [1/8, -1/8, -1/8, 1/8, -1/8, -1/8, 1/8, 1/8]]
    exps = [np.dot(c, res) for c in coefs]
    coefs_mat = [Eu-3, Eg-2, Eu-1, Eg, Eu+1, Eg+2, Eu+3]
    return np.dot(exps, coefs_mat)

def real_cost_function(thetas, shots, t, mu, E, width, dist, D, alpha, r_0, t_0, backend):
    qcs, qcs_z = get_all_qcs(thetas)
    Eg = morse_pot(dist, D, alpha, r_0)
    Eu = morse_pot_u(dist, D, alpha, r_0, t_0)
    off_exps = [res_process_exp(r, shots) for r in ibmsim_real(qcs, backend, shots)]
    dia_exps = [res_process_exp_z(r, shots) for r in ibmsim(qcs_z, backend, shots)]
    states_exps_ = []
    for i in range(len(dia_exps)):
        states_exps_.append(offdiagonal(off_exps[i:(i+1)*10], t, mu, E, width)+diagonal(dia_exps[i], Eg, Eu))
    ssvqe_coefs = [1, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02]
    return np.dot(states_exps_, ssvqe_coefs)

global states_exps
states_exps = []
def callback(thetas):
    qcs, qcs_z = get_all_qcs(thetas)
    Eg = morse_pot(dist,  2.79, 0.72/0.7, 2*0.7)
    Eu = morse_pot_u(dist,  2.79, 0.72/0.7, 2*0.7, -0.15)
    off_exps = [res_process_exp(r, 10000) for r in ibmsim(qcs, 10000)]
    dia_exps = [res_process_exp_z(r, 10000) for r in ibmsim(qcs_z, 10000)]
    for i in range(len(dia_exps)):
        states_exps.append(offdiagonal(off_exps[i:(i+1)*10], 0, mu, E, width)+diagonal(dia_exps[i], Eg, Eu))

def main(distance, ts, mu, E, width, shots, backend, gpath):
    res = np.zeros([len(ts), 7, len(distance)])
    exacts = np.zeros([len(ts), 7, len(distance)])
    for i, t in enumerate(ts):
        for m, dist in enumerate(distance):
        #dist = 1.2
            states_exps = []
            opt = minimize(real_cost_function, (0,0),
                    args=(shots, t, mu, E, width, dist, 2.79, 0.72/0.7, 2*0.7, -0.15, provider.get_backend(backend)),
                    method="SLSQP",
                    callback=callback,
                    #tol=1e-3,
                    bounds=((-np.pi, np.pi), (-np.pi, np.pi)))
            res[i, :, m] = np.sort(np.array(states_exps).reshape(-1, 7)[-1])
            exacts[i, :, m] = np.sort(LA.eig(get_quasi_floquet_matrix(t, mu, E, \
                                                    morse_pot(dist, 2.79, 0.72/0.7, 2*0.7), \
                                                    morse_pot_u(dist, 2.79, 0.72/0.7, 2*0.7, -0.15), width, 3))[0])
    get_graph(res, exacts, gpath)


def get_graph(res, exacts, gpath):
    fig = plt.figure(figsize=(12,16))
    gs = gridspec.GridSpec(2,2)

    axLU = plt.subplot(gs[0,0])
    axRU = plt.subplot(gs[0,1])
    axLB = plt.subplot(gs[1,0])
    axRB = plt.subplot(gs[1,1])
    axLU.set_xlim([0,6])
    axLU.set_ylim([-0.5,4])
    axLU.plot(distance, res[0,:,:], "k", label="ibm_Kawasaki")
    axLU.plot(distance, exacts[0,:,:], "k--", label="classical")
    axLU.legend()
    axRU.set_xlim([0,6])
    axRU.set_ylim([-0.5,4])
    axRU.plot(distance, res[1,:,:], "k", label="ibm_Kawasaki")
    axRU.plot(distance, exacts[1,:,:], "k--", label="classical")
    axRU.legend()
    axLB.set_xlim([0,6])
    axLB.set_ylim([-0.5,4])
    axLB.plot(distance, res[2,:,:], "k", label="ibm_Kawasaki")
    axLB.plot(distance, exacts[2,:,:], "k--", label="classical")
    axLB.legend()
    axRB.set_xlim([0,6])
    axRB.set_ylim([-0.5,4])
    axRB.plot(distance, res[3,:,:], "k", label="ibm_Kawasaki")
    axRB.plot(distance, exacts[3,:,:], "k--", label="classical")
    axRB.legend()

    plt.savefig(gpath)


if __name__ == "__main__":
    mu = 1
    E = 0.8 #5.337e-3
    delta = 1
    omega = 1
    width = 100
    n_states = 3
    ts = [0, 150, 200, 300]
    distance = np.linspace(0.45, 6, 70)
    shots = 10000
    backend = "ibm_kawasaki"
    gpath = "floquet_77_result.png"

    main(distance, ts, mu, E, width, shots, backend, gpath)