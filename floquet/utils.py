import numpy as np

import qiskit
from qiskit.providers.aer.noise import NoiseModel
from qiskit import IBMQ, QuantumCircuit, ClassicalRegister, QuantumRegister, execute, Aer

from qiskit.providers.ibmq import least_busy

IBMQ.load_account()
provider = IBMQ.get_provider(hub='ibm-q-utokyo', group='internal', project='yamanouchi-lab')
provider.backends()

def Hf_X(theta):
    qnodes = QuantumRegister(1,'qc')
    cnodes = ClassicalRegister(1,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    qc.ry(theta, qnodes[0])
    qc.x(qnodes[0])
    qc.measure(qnodes, cnodes)
    return qc

def Hf_X_ortho(theta):
    qnodes = QuantumRegister(1,'qc')
    cnodes = ClassicalRegister(1,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    qc.x(qnodes[0])
    qc.ry(theta, qnodes[0])
    qc.x(qnodes[0])
    qc.measure(qnodes, cnodes)
    return qc

def Hf_I(theta):
    qnodes = QuantumRegister(1,'qc')
    cnodes = ClassicalRegister(1,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    qc.ry(theta, qnodes[0])
    qc.measure(qnodes, cnodes)
    return qc

def Hf_I_ortho(theta):
    qnodes = QuantumRegister(1,'qc')
    cnodes = ClassicalRegister(1,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    qc.x(qnodes[0])
    qc.ry(theta, qnodes[0])
    qc.measure(qnodes, cnodes)
    return qc

def Hf_H(theta):
    qnodes = QuantumRegister(1,'qc')
    cnodes = ClassicalRegister(1,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    qc.ry(theta, qnodes[0])
    qc.h(qnodes[0])
    qc.measure(qnodes, cnodes)
    return qc

def Hf_H_ortho(theta):
    qnodes = QuantumRegister(1,'qc')
    cnodes = ClassicalRegister(1,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    qc.x(qnodes[0])
    qc.ry(theta, qnodes[0])
    qc.h(qnodes[0])
    qc.measure(qnodes, cnodes)
    return qc

def add_count(qc_result):
    for r in qc_result:
        if not "1" in r.keys():
            r["1"] = 0
        if not "0" in r.keys():
            r["0"] = 0
        r = dict(sorted(r.items()))
    return qc_result

def process_result(re, shots):
    vecs = []
    re = add_count(re)
    for r in re:
        vec = [np.sqrt(r['0']/shots), np.sqrt(r['1']/shots)]
        vecs.append(vec)
    return vecs

def get_eigenvalues(re, shots, coef, eigenstate):
    vecs = process_result(re, shots)
    eigenvalues = []
    for vec in vecs:
        eigenvalues.append(coef*np.dot(vec, eigenstate))
    return eigenvalues

def ibmsim(circ, shots, _type):
    if _type == "qasm":
        simulator = Aer.get_backend('qasm_simulator')
    elif _type == "statevector":
        simulator = Aer.get_backend('statevector_simulator')
    else:
        print("Please choose simulator type; qasm or statevector.")
    return execute(circ, simulator, shots=shots).result().get_counts()

def ibmKawasaki_sim(circ, simulator, noise_model, shots):
    return execute(circ, simulator, shots=shots, noise_model=noise_model).result().get_counts()

def ibmsim_real(circ, backend, shots):
    res = execute(circ, backend, shots=shots).result().get_counts()
    return res

def get_eigens(qc_results_HI, shots):
    eigenstates_list = []
    for n in range(int(len(qc_results_HI)/4)):
        re_H = qc_results_HI[4*n]
        re_H_ortho = qc_results_HI[4*n+1]
        re_I = qc_results_HI[4*n+2]
        re_I_ortho = qc_results_HI[4*n+3]
        eigenstates = np.zeros([2,2])
        e = 0
        for h, i in [(re_H, re_I), (re_H_ortho, re_I_ortho)]:
            h = dict(sorted(h.items()))
            i = dict(sorted(i.items()))
            if detect_sign(h, i, shots) > 0:
                eigenstates[e, :] = np.sqrt(np.array(list(i.values()))/shots)
            elif detect_sign(h, i, shots) < 0:
                eigenstates[e, :] = np.array([np.sqrt(int(list(i.values())[0])/shots), -np.sqrt(int(list(i.values())[1])/shots)])
            e += 1
        eigenstates_list.append(eigenstates)
    return eigenstates_list

def detect_sign(H_result, I_result, shots):
    H_result = add_count([H_result])[0]
    I_result = add_count([I_result])[0]
    return (H_result["0"]-0.5*I_result["0"]-0.5*I_result["1"])/shots

def pulse(t, omega, width):
    return np.exp(-((t-300)/width)**2)*np.exp(-1.j*omega*(t-300))

def amplitude(t, width):
    return np.exp(-((t-300)/width)**2)

def get_quasi_floquet_matrix(t, mu, E, Eg, Eu, width, n_states):
    fm = np.zeros([8, 8], dtype=np.complex)
    for i in range(8):
        if (i-n_states)%2 == 0:
            fm[i,i] = Eg + (i-n_states)
        else:
            fm[i,i] = Eu + (i-n_states)
    for j in range(8):
        fm[j,j-1] = mu*E*amplitude(t, width)
        fm[j-1,j] = mu*E*amplitude(t, width)
    fm[0,7] = 0
    fm[7,0] = 0
    return fm

def time_propagate_population(t, coef_ps, eigenvalues, init, omega, n_states):
    #coef_ps = eigens[1][:, 2*(n_states+1)]
    #eigenvalues = eigens[0]
    #vec = [coef_p_a/init, coef_q_a/init]
    #coef_p, coef_q = vec #/np.linalg.norm(vec)
    norm = np.linalg.norm(coef_ps, ord=2)
    pop_a = 0
    for i, (c, e) in enumerate(zip(coef_ps, eigenvalues)):
        pop_a += (c/norm)*np.exp(-1.j*(e+(i-n_states)*omega)*t)
    #pop_b = np.abs(b_norm*(coef_p_b*np.exp(-1.j*(eigenvalue_p+1)*t)+coef_q_b*np.exp(-1.j*(eigenvalue_q+1)*t)))**2
    return np.abs(pop_a)**2

def time_propagate_population_initial(t, eigens, init, omega):
    coef_ps= eigens[1][0, :]
    eigenvalues = eigens[0]
    #vec = [coef_p_a/init, coef_q_a/init]
    #coef_p, coef_q = vec #/np.linalg.norm(vec)
    #pop_a = 0
    #for i, (c, e) in enumerate(zip(coef_ps, eigenvalues)):
    pop_a  = coef_ps[0]*coef_ps[0]*np.exp(-1.j*eigenvalues[0]*t)+coef_ps[1]*coef_ps[1]*np.exp(-1.j*eigenvalues[1]*t)
    #pop_b = np.abs(b_norm*(coef_p_b*np.exp(-1.j*(eigenvalue_p+1)*t)+coef_q_b*np.exp(-1.j*(eigenvalue_q+1)*t)))**2
    return np.abs(pop_a)**2

def get_eigens_process(qc_params, modeling_params_list):
    """
    qc_params: [shots, qc_type, device=None]
    modeling_params: [mu, E, delta]
    """
    thetas = np.linspace(0, 2*np.pi, 100)
    minimized_thetas, eigenvalues_list = get_ssvqe_exp(thetas, *qc_params, modeling_params_list)
    #TODO Find the minimized theta by more sophiscated way
    #res_for_sign_detection = get_qc_results_for_eigenstates_with_sign(minimized_thetas, qc_params)
    #eigenstates_list = get_eigens(res_for_sign_detection, qc_params[0])
    #return [(eigenvalues, eigenstates) for eigenvalues, eigenstates in zip(eigenvalues_list, eigenstates_list)], minimized_thetas
    return eigenvalues_list

def get_qc_results_for_eigenstates_with_sign(minimized_thetas, qc_params):
    qcs = []
    for minimized_theta in minimized_thetas:
        qcs.extend([Hf_H(minimized_theta), Hf_H_ortho(minimized_theta), \
                    Hf_I(minimized_theta), Hf_I_ortho(minimized_theta)])

    shots, qc_type, device = qc_params
    if qc_type == "simulator":
        results = ibmsim(qcs, shots)
    elif qc_type == "noise_simulator":
        print("using the noise config of the Kawasaki machine")
        device = provider.get_backend("ibm_kawasaki")
        noise_model = NoiseModel.from_backend(device)
        simulator = Aer.get_backend('qasm_simulator')
        results = ibmKawasaki_sim(qcs, simulator, noise_model, shots)
    elif qc_type == "real":
        if device == "least_busy":
            device = least_busy(provider.backends(simulator=False, operational=True))
        elif device == "kawasaki":
            device = provider.get_backend("ibm_kawasaki")
        elif device == "jakarta":
            device = provider.get_backend("ibmq_jakarta")
        else:
            print("set 'device' param which will be used for calculation.")
        print("Calculation in {}...".format(device))
        if len(qcs) > 300:
            results = []
            for i in range(0,len(qcs),300):
                results.extend(ibmsim_real(qcs[i:i+300], device, shots))
        else:
            results = ibmsim_real(qcs, device, shots)
    else:
        print("set 'qc_type' param.")
    print("Second calculation on QC has been done.")
    #results = [results[n:n+4] for n in range(int(len(results)/4))]
    return results

def measure_x(theta, ortho=False):
    qnodes = QuantumRegister(1,'qc')
    cnodes = ClassicalRegister(1,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if ortho:
        qc.x(qnodes[0])
    qc.ry(theta, qnodes[0])
    qc.h(qnodes[0])
    qc.measure(qnodes, cnodes)
    return qc

def measure_z(theta, ortho=False):
    qnodes = QuantumRegister(1,'qc')
    cnodes = ClassicalRegister(1,'cr')
    qc = QuantumCircuit(qnodes, cnodes)
    if ortho:
        qc.x(qnodes[0])
    qc.ry(theta, qnodes[0])
    qc.measure(qnodes, cnodes)
    return qc

def get_exp(ret_x, ret_z, shots, mu, E, Eg, Eu, omega, width, n_states):
    ret_x = add_count([ret_x])[0]
    ret_z = add_count([ret_z])[0]
    return (Eg+Eu)/2+mu*E*(ret_x["0"]-ret_x["1"])/shots-((Eg-Eu)/2)*(ret_z["0"]-ret_z["1"])/shots
    
def get_ssvqe_exp(thetas, shots, qc_type, device, params_list):
    qcs = []
    for theta in thetas:
        qc_x = measure_x(theta)
        qc_z = measure_z(theta)
        qc_x_ortho = measure_x(theta, ortho=True)
        qc_z_ortho = measure_z(theta, ortho=True)
        qcs += [qc_x, qc_z, qc_x_ortho, qc_z_ortho]
    print("Circuits are collected.")

    if qc_type == "simulator":
        results = ibmsim(qcs, shots)
    elif qc_type == "noise_simulator":
        print("using the noise config of the Kawasaki machine")
        device = provider.get_backend("ibm_kawasaki")
        noise_model = NoiseModel.from_backend(device)
        simulator = Aer.get_backend('qasm_simulator')
        results = ibmKawasaki_sim(qcs, simulator, noise_model, shots)
    elif qc_type == "real":
        if device == "least_busy":
            device = least_busy(provider.backends(simulator=False, operational=True))
        elif device == "kawasaki":
            device = provider.get_backend("ibm_kawasaki")
        elif device == "jakarta":
            device = provider.get_backend("ibmq_jakarta")
        else:
            print("set 'device' param which will be used for calculation.")
        print("Calculation {} quantum circuits in {}...".format(len(qcs), device))
        if len(qcs) > 300:
            results = []
            for i in range(0,len(qcs),300):
                results.extend(ibmsim_real(qcs[i:i+300], device, shots))
        else:
            results = ibmsim_real(qcs, device, shots)
    else:
        print("set 'qc_type' param.")
    print("First calculation on QC has been done.")

    minimized_thetas = []
    eigenvalues_list = []
    exps = []
    for params in params_list:
        exps = [(get_exp(results[4*t], results[4*t+1], shots, *params), get_exp(results[4*t+2], results[4*t+3], shots, *params)) for t in range(len(thetas))]
        ssvqe_exps = [ex[0]+0.5*ex[1] for ex in exps]
        minimized_thetas.append(thetas[list(ssvqe_exps).index(min(ssvqe_exps))])
        eigenvalues_list.append(exps[list(ssvqe_exps).index(min(ssvqe_exps))])
    return minimized_thetas, eigenvalues_list

def morse_pot(dist, D, alpha, r_0):
    return D*(np.exp(-2*alpha*(dist-r_0))-2*np.exp(-alpha*(dist-r_0)))

def morse_pot_u(dist, D, alpha, r_0, t):
    return D*(np.exp(-2*alpha*(dist-r_0))-2*t*np.exp(-alpha*(dist-r_0)))