import matplotlib.pyplot as plt

from constants import AVOGADRO_CONSTANT, LV_constants, MM_constants, \
    PGAR_constants, LO_constants
from solve import stochastic_simulation_algorithm, euler_maruyama_algorithm

# Font configuration
font = { "family": "monospace", "size": 15 }
_font = { "family": "monospace", "size": 10 }
sub_font = { "family": "monospace", "size": 10, "weight": "bold" }
title_font = { "family" : "monospace", "size": 20, "weight": "bold" }

# Colors for CLI text
RUN = "\033[1;32m::\033[0m"
QST = "\033[1;35m=>\033[0m"
PLOT = "\033[1;34m->\033[0m"
INF = "\033[1;34m[#]\033[0m"
TITL = "\033[1;35m"
RES = "\033[0m"
ERR = "\033[1;31m"

def plot(ts: list[float], _xs: list[list[int]], _ys: list[list[float]], title:
        Optional[str], labels: Optional[list[str]]):
    """Visualization of the data:
    `ts`: time points
    `_xs`: record of the number of molecules of the system
    `_ys`: record of the concentration of each molecule
    `title`: title for the final plot
    `labels`: list of labels for each differen species of molecules
    """

    # Font config
    font = { "family": "monospace", "size": 15 }
    _font = { "family": "monospace", "size": 10 }
    title_font = { "family" : "monospace", "size": 20, "weight": "bold" }

    # Collection of colors to be used for each species of molecules
    c = [
        "tab:blue", "tab:orange", "tab:purple", "black",
        "magenta",  "darkorange", "tab:gray",   "forestgreen"
    ]

    xs = np.array(_xs) # Number of molecules
    kinds = len(xs[0]) # Number of different species

    # Concentration
    ys = np.array(_ys)
    if not _ys is None:
        # Only start a subplot if the concentration was specified
        plt.subplot(211)

    # Plot the evolution of each species
    for k in range(kinds):
       plt.plot(ts, xs[:, k], color=c[k])
    plt.xlabel("Time", font)
    plt.ylabel("Nr. molecules", font)
    plt.legend(labels, prop=_font)

    if _ys != None:
        plt.subplot(212)
        for k in range(kinds):
            plt.plot(ts, ys[:, k], color=c[k])
        plt.xlabel("Time", font)
        plt.ylabel("Concentration (M)", font)
        plt.legend(labels, prop=_font)
        plt.suptitle(title, fontproperties=title_font)
    else:
        plt.title(title, font)
    plt.show()


def michaelis_menten(X_1, X_2, X_3, X_4):
    print(f"\n{TITL}{' Michael-Menten ':-^79}{RES}\n")

    X0, st_mat, p_functions, vol = MM_constants()

    print(f"{RUN} Calculating SSA simulation and plotting...\n")
    ts, Xs, ys = stochastic_simulation_algorithm(X0, st_mat, p_functions, vol=vol, record=True)

    title = "Michael-Menten: SSA"
    ls = ["Substrate", "Enzyme", "Complex", "Product"]
    plot(ts, Xs, ys, title, ls)

    print(f"{RUN} Calculating Euler-Maruyama simulation and plotting...")
    ts, Xs, ys = euler_maruyama_algorithm(X0, st_mat, p_functions, record=True)

    title="Michael-Menten: Euler-Maruyama"
    plot(ts, Xs, ys, title, ls)


def dimerization_kinetics():
    print(f"\n{TITL}{' Dimerisation kinetics ':-^79}{RES}\n")
    ls = ["$P$", "$P_2$"]

    X0, st_mat, p_functions, tf, tstep, vol, _ = LV_constants()

    print(f"{RUN} Calculating and plotting the evolution over 10s (SSA)\n")
    ts, Xs, ys = stochastic_simulation_algorithm(X0, st_mat, p_functions, tf, vol, record=True)

    title = "Dimerisation kinetics: SSA"
    plot(ts, Xs, ys, title, ls)

    print(f"{RUN} Calculating and plotting the evolution"
          " over 10s (Euler-Maruyama)\n")
    ts, Xs, ys = euler_maruyama_algorithm(X0, st_mat, [dim, diss], vol, tf, tstep)

    title = "Dimerisation kinetics: Euler-Maruyama"
    plot(ts, Xs, ys, title, ls),


    print(f"{RUN} Running 20 simulation with SSA...")
    for _ in range(20):
        ts, Xs, _ = stochastic_simulation_algorithm(x0, r, [dim, diss], vol, tspan=10, _plot=False)
        plt.plot(ts, [X[0] for X in Xs], color="black")
    print(f"{PLOT} Plotting the overlay of 20 runs for the protein...\n")
    plt.title("Protein evolution: 20 runs", title_font)
    plt.ylim(0, x0[0])
    plt.xlabel("Time", font)
    plt.ylabel("Nr. molecules", font)
    plt.legend(["$P$"], prop=_font)
    plt.show()

    runs = 10000
    print(f"{RUN} Running {runs} stochastic simulations...")
    ps = []
    for _ in range(runs):
        _, Xs, _ = stochastic_simulation_algorithm(x0, r, [dim, diss], vol, tspan=10, _plot=False)
        ps.append(Xs[-1][0])

    print(f"{PLOT} Plotting the density histogram of the protein"
          " at time t = 10...\n")
    plt.hist(ps, bins=int(185/5), density=True,
             edgecolor="black", color="orange")
    plt.title("Density histogram of $P$: SSA", title_font)
    plt.xlabel("$P(10)$", font)
    plt.ylabel("Density", font)
    plt.show()

    runs = 1000
    print(f"{RUN} Running {runs} Euler-Maruyama simulations...")
    _L = 500
    _tspan = 10
    tau = _tspan / _L
    ts = np.array([n * tau for n in range(_L)])
    pem = np.zeros(_L)
    for _ in range(runs):
        _, xs, _ = euler_maruyama_algorithm(x0, r, [dim, diss], vol,
                          tspan=_tspan, L=_L, _plot=False)
        xs = np.array(xs)
        for t in range(_L):
            pem[t] += xs[t, 0]
    pem *= 1/runs # Take mean

    # Calculate standard deviation
    std_run = 3
    rec = np.zeros((std_run, _L))
    for run in range(std_run):
        _, xs, _ = euler_maruyama_algorithm(x0, r, [dim, diss], vol,
                          tspan=_tspan, L=_L, _plot=False)
        xs = np.array(xs)
        for t in range(_L):
            rec[run, t] = xs[t, 0]
    std = np.zeros(_L)
    for t in range(_L):
        std[t] = np.std(rec[:, t])

    print(f"{PLOT} Plotting sample mean for the protein evolution...")
    plt.plot(ts, pem, color="orange")
    plt.plot(ts, pem + std, color="magenta")
    plt.plot(ts, pem - std, color="green")
    plt.ylim(0, x0[0])
    plt.ylabel("Nr. molecules", font)
    plt.xlabel("Time", font)
    plt.legend(
        ["Mean of P",
         f"Mean of P + {std_run} samples std. deviation",
         f"Mean of P - {std_run} samples std. deviation"],
        prop=_font
    )
    plt.show()


def auto_regulatory_genetic_network():
    print(f"\n{TITL}{' Auto-regulatory network ':-^79}{RES}\n")
    ls = ["Gene", "$P_2 \\cdot$Gene", "RNA", "$P$", "$P_2$"]

    X0, st_mat, p_functions, *_ = PGAR_constants()

    print(f"{RUN} Calculating evolution of the network through"
          " a time span of 5000s (SSA)...")
    ts, Xs, _ = stochastic_simulation_algorithm(X0, st_mat, p_functions, 5000, record=True)
    plt.figure()

    print(f"{PLOT} Plotting evolution for time in [0, 5000]...")
    plt.subplot(321)
    plt.plot(ts, xs[:, 2], color="orange")
    plt.xlabel("Time", _font)
    plt.ylabel("RNA", _font)

    plt.subplot(323)
    plt.plot(ts, xs[:, 3], color="purple")
    plt.xlabel("Time", _font)
    plt.ylabel("$P$", _font)

    plt.subplot(325)
    plt.plot(ts, xs[:, 4], color="blue")
    plt.xlabel("Time", _font)
    plt.ylabel("$P_2$", _font)

    print(f"{PLOT} Plotting evolution for time in [0, 250]...\n")
    _ts = ts[ts <= 250]

    plt.subplot(322)
    plt.plot(_ts, xs[:len(_ts), 2], color="orange")
    plt.xlabel("Time", _font)
    plt.ylabel("RNA", _font)

    plt.subplot(324)
    plt.plot(_ts, xs[:len(_ts), 3], color="purple")
    plt.xlabel("Time", _font)
    plt.ylabel("$P$", _font)

    plt.subplot(326)
    plt.plot(_ts, xs[:len(_ts), 4], color="blue")
    plt.xlabel("Time", _font)
    plt.ylabel("$P_2$", _font)

    plt.suptitle("Auto-regulatory genetic network over 5000s (SSA)",
                 fontproperties=title_font)
    plt.show()


    print(f"{RUN} Calculating evolution of the network through"
          " a time span of 10s (SSA)...")
    ts, xs, _ = ssa(x0, r, a, tspan=10, _plot=False)
    xs = np.array(xs)
    print(f"{PLOT} Plotting evolution of P over 10s...\n")
    plt.plot(ts, xs[:, 3], color="purple")
    plt.title("Evolution of $P$ over 10s", title_font)
    plt.ylabel("Nr. molecules", font)
    plt.xlabel("Time", font)
    plt.show()

    def density_P10(runs: int, _c2: Optional[float], index: int):
        str_c2 = f"k2 = {_c2}"
        if _c2 == None:
            str_c2 = "k2 uniformly in [0.005, 0.03)"
        print(f"{RUN} Running {runs} simulations over 10s (SSA)"
              f" for {str_c2} ...")

        unif = False
        if _c2 == None:
            unif = True

        ps = []
        for _ in range(runs):
            if unif:
                _c2 = rand.uniform(0.005, 0.03)
            def _trsc(x):
                """Transcription with altered constant `c2`
                Gene -> Gene + RNA
                """
                return _c2 * x[0]

            _a = [rb, rrb, _trsc, trans, dim, diss, rnadeg, pdeg]

            _, xs, _ = ssa(x0, r, _a, tspan=10, labels=ls, _plot=False)
            ps.append(xs[-1][-2])

        print(f"{PLOT} Plotting density histogram of P over 10s...\n")
        plt.subplot(3, 1, index)
        plt.hist(ps, bins=int(185/5), density=True,
                edgecolor="black", color="orange")
        plt.xlabel("$P(10)$", _font)
        plt.ylabel("Density", _font)
        if not unif:
            plt.title(f"Density for $k_2 = {_c2}$", sub_font)
        else:
            plt.title(f"$k_2$ uniformly choosen in [0.005, 0.03)", sub_font)

    plt.subplots(constrained_layout=True)
    density_P10(1000, 0.01, 1)
    density_P10(1000, 0.02, 2)
    density_P10(1000, _c2=None, index=3)
    plt.suptitle("Density histogram of $P$ over 10s", fontproperties=title_font)
    plt.show()


def lac_operon():
    print(f"\n{TITL}{' Lac-operon model ':-^79}{RES}\n")

    ls = ["Gene", "RNA", "P", "P_2"]

    X0, st_mat, p_functions, tf, *_, intervention = PGAR_constants()

    print(f"{RUN} Running (SSA) simulation of the Lac-operon model over 50,000s...")
    ts, Xs, _ = stochastic_simulation_algorithm(X0, st_mat, p_functions, tf, record=True)

    print(f"{PLOT} Plotting results...")
    plt.subplot(311)
    plt.plot(ts, xs[:, 7], color="orange")
    plt.xlabel("Time", _font)
    plt.ylabel("Lactose", _font)

    plt.subplot(312)
    plt.plot(ts, xs[:, 5], color="purple")
    plt.xlabel("Time", _font)
    plt.ylabel("RNA", _font)

    plt.subplot(313)
    plt.plot(ts, xs[:, 6], color="tab:blue")
    plt.xlabel("Time", _font)
    plt.ylabel("Z", _font)

    plt.suptitle("Lac-Operon model for 50,000s (SSA)", fontproperties=title_font)
    plt.show()


def compare():
    """Comparison of the algorithms SSA and Euler-Maruyama"""
    print(f"\n{TITL}{' Comparison of SSA and Euler Maruyama ':-^79}{RES}\n")
    runs = 10_000

    print(
        f"{INF} In order to compare the algorithms SSA and Euler-Maruyama we'll use\n"
         "    the dimerisation kinetics model --- of which we'll analyse the evolution\n"
        f"    of the protein P throughout a time-span of 10s and {runs} independent\n"
         "    simulations.\n"
    )

    X0, st_mat, p_functions, tf, _, vol, _ = LV_constants()

    print(
        f"{RUN} Running {runs} simulations using both SSA and Euler-Maruyama,\n"
         "   this may take some time..."
    )

    print(
f"""{INF} It should be noted that the stochastic algorithm does not allow us to
    correctly calculate the point-wise mean:
    From its nature, each simulation gets a randomly generated list of
    time-points. In order to calculate the mean and standard deviation of each
    point through the simulations, we randomly choose a standard time-point
    array, by running a single simulation of SSA, that will be used as our basis
    for the record of the results throughout the simulations.
"""
    )

    print(
f"""{INF} In order to have a well behaved mean and standard deviation for SSA,
    we have to take care of two possibilities. In what follows, let `N` be the
    length of the standard time-point array and `M` be the length of the new
    time-point array:
    * If `N < M`, we choose to ignore the last `M - N` protein-states of the new
      array and only store the first `N` points.
    * If `M <= N`, we record the new protein-states as-is and choose to repeat,
      for the last `N - M` elements, the last obtained protein number. This is
      done in order to avoid bad behaviour of the mean and standard deviation.
"""
    )

    # Since the time points vary over the simulations using SSA, I'll determine
    # a standard array of time points `ts_ssa` as a support to calculate the
    # sample mean and standard deviation, this is chosen randomly by running a
    # single simulation of `ssa`
    ts_ssa, Xs0_ssa, _ = stochastic_simulation_algorithm(X0, st_mat, p_functions, tf, record=True)
    tstep = len(ts_ssa)

    # Record state of the protein for each simulation for SSA
    rec_ssa = np.zeros((runs, tstep))
    rec_ssa[0] = Xs0_ssa[:, 0]

    # `tau` is the fixed time-step used in `elmaru` and `ts_em` represents the
    # array of time-points that will be used by `elmaru` throughout each simulation
    tau = tf / tstep
    ts_em = np.array([n * tau for n in range(L)], dtype=float)

    # Record state of the protein for each simulation for Euler-Maruyama
    rec_em = np.zeros((runs, L))

    # In order to use the same simulation-loop as `ssa`, we run a single
    # simulation of `elmaru`
    _, xs0_em, _ = elmaru(x0, r, a, vol, tspan=tspan, L=L, _plot=False)
    xs0_em = np.array(xs0_em)
    rec_em[0] = xs0_em[:, 0]

    for run in range(1, runs):
        # `ssa` simulation: A certain care needs to go with the process of
        # recording the simulation, since the length of the time-points vary
        # through the simulations
        tnew, xs_ssa, _ = ssa(x0, r, a, tspan=tspan, _plot=False)
        xs_ssa, N = np.array(xs_ssa), len(tnew)
        M = min(N, L)
        rec_ssa[run, :M] = xs_ssa[:M, 0]
        if N < L:
            # Repeat the last value obtained in the simulation for the
            # protein. This is determined in order to bring stability to the
            # standard deviation at the end of the simulation, that is, near 10s
            last_value = xs_ssa[-1, 0]
            rec_ssa[run, N:] = np.repeat(last_value, L - N)

        # `elmaru` simulation
        _, xs_em, _ = elmaru(x0, r, a, vol, tspan=tspan, L=L, _plot=False)
        xs_em = np.array(xs_em)
        rec_em[run] = xs_em[:, 0]

    # Record of protein number at 10s
    ps_ssa = rec_ssa[:, -1]
    ps_em = rec_em[:, -1]

    # Point-wise mean
    mean_ssa = np.sum(rec_ssa, axis=0) / runs
    mean_em = np.sum(rec_em, axis=0) / runs

    # Point-wise standard deviation
    std_ssa = np.std(rec_ssa, axis=0)
    std_em = np.std(rec_em, axis=0)

    print(f"{PLOT} Plotting results...\n")
    print(
f"""{INF} The obtained results show that the stability of the Euler-Maruyama algorithm is
    far greater than that of SSA. We should not be deceived by the convergence
    of the SSA mean to the results of Euler-Maruyama near the 10s, since this
    region is where we got the least amount of protein-states throughout the
    simulation by the way we arranged the recording of the SSA simulations and
    therefore is not really reliable.
"""
    )

    plt.plot(ts_ssa, mean_ssa + 3 * std_ssa, color="tab:blue", linestyle=":")
    plt.plot(ts_ssa, mean_ssa, color="tab:blue")
    plt.plot(ts_ssa, mean_ssa - 3 * std_ssa, color="tab:blue", linestyle="-.")

    plt.plot(ts_em, mean_em + 3 * std_em, color="tab:purple", linestyle=":")
    plt.plot(ts_em, mean_em, color="tab:purple")
    plt.plot(ts_em, mean_em - 3 * std_em, color="tab:purple", linestyle="-.")

    plt.ylim(0, x0[0])
    plt.ylabel("Nr. molecules", font)
    plt.xlabel("Time", font)
    plt.legend([
        f"Mean of $P + 3 \\sigma$ (SSA)",
         "Mean of $P$     (SSA)",
        f"Mean of $P - 3 \\sigma$ (SSA)",
        f"Mean of $P + 3 \\sigma$ (Euler-Maruyama)",
         "Mean of $P$     (Euler-Maruyama)",
        f"Mean of $P - 3 \\sigma$ (Euler-Maruyama)"],
        prop=_font)
    plt.title("Comparison of SSA and Euler Maruyama", title_font)
    plt.show()

    print(f"{PLOT} Plotting density histogram for the protein at 10s...\n")

    print(
f"""{INF} Notice that the density histogram of SSA is way sparser when compared to that of
    Euler-Maruyama algorithm. This only reinforces that, if the user wants
    reliability and stability through their simulation, the Euler-Maruyama algorithm
    should be the prefered choice.
"""
    )

    plt.subplot(211)
    plt.hist(ps_ssa, bins=int(185/5), density=True,
             edgecolor="black", color="tab:blue")
    plt.ylabel("Density")
    plt.xlabel("$P(10)$")

    plt.subplot(212)
    plt.hist(ps_em, bins=int(185/5), density=True,
             edgecolor="black", color="tab:orange")
    plt.ylabel("Density")
    plt.xlabel("$P(10)$")

    plt.suptitle("Density histogram of $P(10)$: SSA and Euler-Maruyama",
                 fontproperties=title_font)
    plt.show()


def main(reaction_model: ReactionModel):
    match reaction_model:
        case ReactionModel.Lotka_Volterra:
            dimerisation_kinetics()
        case ReactionModel.Michaelis_Menten:
            michaelis_menten()
        case ReactionModel.auto_regulating_genetic_network:
            auto_regulating_genetic_network()
        case ReactionModel.lac_operon:
            lac_operon()
        case ReactionModel.comparison:
            compare()
        case _:
            raise ValueError('Please input one of the following:\n'
                            '1 - Dimerisation kinetics;\n'
                            '2 - Michael-Menten;\n'
                            '3 - Auto-regulating genetic network;\n'
                            '4 - Lac-operon model.\n'
                            '5 - Compare CME (SSA) and CLE (Euler-Maruyama).\n')


if __name__ == '__main__':
    from sys import argv
    if len(argv) < 2:
            raise ValueError('Please input one of the following:\n'
                            '1 - Dimerisation kinetics;\n'
                            '2 - Michael-Menten;\n'
                            '3 - Auto-regulating genetic network;\n'
                            '4 - Lac-operon model.\n'
                            '5 - Compare CME (SSA) and CLE (Euler-Maruyama).\n')

    main(ReactionModel(argv[1]))
