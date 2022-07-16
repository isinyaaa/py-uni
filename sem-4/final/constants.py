from enum import Enum

import numpy as np


# CONSTANTS
AVOGADRO_CONSTANT   = 6.02214086e23
SYSTEM_VOLUME       = 1e-15


class ReactionModel(Enum):
    Lotka_Volterra = 1
    Michaelis_Menten = 2
    auto_regulating_genetic_network = 3
    lac_operon = 4
    comparison = 5


# Michaelis-Menten
def MM_constants():
    X0 = np.array([301, 120, 0, 0]) # Initial number of molecules
    st_mat = np.array([
        [-1, -1,  1, 0], # Binding
        [ 1,  1, -1, 0], # Dissociation
        [ 0,  1, -1, 1]  # Conversion
    ])

    # Binding
    def b(x):
        return 1.66e-3 * x[0] * x[1]

    # Dissociation
    def d(x):
        return 1e-4 * x[2]

    # Conversion
    def c(x):
        return 0.1 * x[2]

    p_functions = [b, d, c]

    tf = 50
    tstep = 250

    # initial_state, stoichiometry_matrix, propensity_functions, final_time, time_step, system_volume, intervention
    return X0, st_mat, p_functions, tf, tstep, SYSTEM_VOLUME, None


# Lotka-Volterra
def LV_constants():
    X0 = np.array([301, 0])
    st_mat = np.array([
        [-2,  1], # Dimerisation
        [ 2, -1]  # Dissociation
    ])

    # Dimerisation
    def dim(x):
        return 1.66e-3 * x[0] * (x[0] - 1)/2

    # Dissociation
    def diss(x):
        return 0.2 * x[1]

    p_functions = [dim, diss]

    tf = 10
    tstep = 500

    # initial_state, stoichiometry_matrix, propensity_functions, final_time, time_step, system_volume, intervention
    return X0, st_mat, p_functions, tf, tstep, SYSTEM_VOLUME, None


# Prokariotic Gene auto-regulation
def PGAR_constants():
    X0 = [10, 0, 0, 0, 0] # Initial state
    st_mat = np.array([
        [-1,  1,  0,  0, -1], # Repression binding
        [ 1, -1,  0,  0,  1], # Reverse repression binding
        [ 0,  0,  1,  0,  0], # Transcription
        [ 0,  0,  0,  1,  0], # Translation
        [ 0,  0,  0, -2,  1], # Dimerisation
        [ 0,  0,  0,  2, -1], # Dissociation
        [ 0,  0, -1,  0,  0], # RNA degeneration
        [ 0,  0,  0, -1,  0]  # Protein degeneration
    ])

    # Repression binding:
    # Gene + P2 -> P2.Gene
    def rb(x):
        return 1 * x[0] * x[-1]

    # Reverse repression binding:
    # P2.Gene -> Gene + P2
    def rrb(x):
        return 10 * x[1]

    # Transcription:
    # Gene -> Gene + RNA
    def trsc(x):
        return 0.01 * x[0]

    # Translation:
    # RNA -> RNA + P
    def trans(x):
        return 10 * x[2]

    # Dimerisation:
    # P + P -> P2
    def dim(x):
        return 1 * 0.5 * x[-2] * (x[-2] - 1)

    # Dissociation:
    # P2 -> P + P
    def diss(x):
        return 1 * x[-1]

    # RNA degeneration:
    # RNA -> nothing
    def rnadeg(x):
        return 0.1 * x[2]

    # Protein degeneration:
    # P -> nothing
    def pdeg(x):
        return 0.01 * x[-2]

    p_functions = [rb, rrb, trsc, trans, dim, diss, rnadeg, pdeg]

    # initial_state, stoichiometry_matrix, propensity_functions, final_time, time_step, system_volume, intervention
    return X0, st_mat, p_functions, None, None, None, None


# lac operon
def LO_constants():
    X0 = np.array([1, 0, 50, 1, 100, 0, 0, 20, 0, 0, 0])
    st_mat = np.array([
        [0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # Inhibitor transcription
        [0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0], # Inhibitor translation
        [0,  0, -1,  0,  0,  0,  0, -1,  1,  0,  0], # Lactose inhibitor binding
        [0,  0,  1,  0,  0,  0,  0,  1, -1,  0,  0], # Lactose inhibitor dissociation
        [0,  0, -1, -1,  0,  0,  0,  0,  0,  1,  0], # Inhibitor binding
        [0,  0,  1,  1,  0,  0,  0,  0,  0, -1,  0], # Inhibitor dissociation
        [0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  1], # RNAp binding
        [0,  0,  0,  1,  1,  0,  0,  0,  0,  0, -1], # RNAp dissociation
        [0,  0,  0,  1,  1,  1,  0,  0,  0,  0, -1], # Transcription
        [0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0], # Translation
        [0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0], # Conversion
        [0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # Inhibitor RNA degradation
        [0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0], # Inhibitor degradation
        [0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0], # Lactose inhibitor degradation
        [0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0], # RNA degradation
        [0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0]  # Z degradation
    ])

    # Inhibitor transcription:
    # IDNA -> IDNA + IRNA
    def in_trsc(x):
        return 0.02 * x[0]

    # Inhibitor translation:
    # IRNA -> IRNA + I
    def in_trans(x):
        return 0.1 * x[1]

    # Lactose inhibitor binding:
    # I + Lactose -> ILactose
    def lac_in_bin(x):
        return 0.005 * x[2] * x[7]

    # Lactose inhibitor dissociation:
    # ILactose -> I + Lactose
    def lac_in_diss(x):
        return 0.1 * x[8]

    # Inhibitor binding:
    # I + Op -> IOp
    def in_bin(x):
        return 1 * x[2] * x[3]

    # Inhibitor dissociation:
    # IOp -> I + Op
    def in_diss(x):
        return 0.01 * x[9]

    # RNAp binding:
    # Op + RNAp -> RNApOp
    def rnap_bin(x):
        return 0.1 * x[3] * x[4]

    # RNAp dissociation:
    # RNApOp -> Op + RNAp
    def rnap_diss(x):
        return 0.01 * x[10]

    # Transcription:
    # RNApOp -> Op + RNAp + RNA
    def trans(x):
        return 0.03 * x[10]

    # Translation:
    # RNA -> RNA + Z
    def transl(x):
        return 0.1 * x[5]

    # Conversion:
    # Lactose + Z -> Z
    def conv(x):
        return 1e-5 * x[6] * x[7]

    # Inhibitor RNA degradation:
    # IRNA -> nothing
    def in_rna_deg(x):
        return 0.01 * x[1]

    # Inhibitor degradation:
    # I -> nothing
    def in_deg(x):
        return 0.002 * x[2]

    # Lactose inhibitor degradation:
    # ILactose -> Lactose
    def lac_in_deg(x):
        return 0.002 * x[8]

    # RNA degradation:
    # RNA -> nothing
    def rna_deg(x):
        return 0.01 * x[5]

    # Z degradation:
    # Z -> nothing
    def z_deg(x):
        return 0.001 * x[6]

    p_functions = [
        in_trsc, in_trans,   lac_in_bin, lac_in_diss,
        in_bin,  in_diss,    rnap_bin,   rnap_diss,
        trans,   transl,     conv,       in_rna_deg,
        in_deg,  lac_in_deg, rna_deg,    z_deg
    ]

    # Adds 10000 lactose molecules to the current state
    def intervention(t: float, x: list[int]) -> bool:
        if t >= 20_000:
            x[7] += 10_000
            return False
        return True

    tf = 50_000

    # initial_state, stoichiometry_matrix, propensity_functions, final_time, time_step, system_volume, intervention
    return X0, st_mat, p_functions, tf, None, None, intervention
