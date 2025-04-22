from backend import *


def gera_grafico(reator: ReatorBatelada, gera_dados: list[Substancia] | None = None):
    result = reator.simular()

    if gera_dados is not None:
        # Gera dict de pontos 'experimentais'
        dict_global = {}
        for subst in gera_dados:
            idx = reator.matriz_substancias.index(subst)

            pts = [0, 99, 199, 299, 399, 499, 599, 699, 799, 899, 999]
            pts_t = [float(result.t[i]) for i in pts]
            pts_y = [float(result.y[idx][i]) for i in pts]

            dict_global[subst] = dict(zip(pts_t, pts_y))
        
        print(dict_global)

    # Cria a figura e o eixo
    fig, ax = plt.subplots(figsize=(8, 6))

    # Para cada substância, plota a concentração versus o tempo
    for subst, conc in zip(reator.matriz_substancias, result.y):
        ax.plot(result.t, conc, label=subst.name)

    ax.set_xlabel("Tempo (s)")
    ax.set_ylabel("Concentração (mol/L)")
    ax.set_title("Evolução das concentrações no reator batelada")
    ax.legend()
    ax.grid(True)

    plt.show()


def ex_1_simulacao():
    pass


def ex_1_otimizacao():
    pass


def ex_2_simulacao():
    pass


def ex_2_otimizacao():
    pass


"""
4 reações diferentes
 1. A -> B + C; {"phi": 0.09, "Ep0": 1, "epsilon": 1, "Abs": 1, "l": 1, "l_reator": 1}
 2. A + B -> D; k= 0,34, ordem 1
 3. A + C -> E; k= 0,7, ordem 1
 """
def ex_3_simulacao():
    subst_a = Substancia("A", "A")
    subst_b = Substancia("B", "B")
    subst_c = Substancia("C", "C")
    subst_d = Substancia("D", "D")
    subst_e = Substancia("E", "E")

    r1 = Reacao("R1", 
                "Reação fotólise", 
                {subst_a: 1}, 
                {subst_b: 1, subst_c: 1},
                TaxaReacaoFotolise(subst_a, {"phi": 0.09, "Ep0": 1, "epsilon": 1, 
                                             "Abs": 1, "l": 1, "l_reator": 1}),
                True
    )

    r2 = Reacao(
        "R2",
        "Reação 2",
        {subst_a: 1, subst_b: 1},
        {subst_d: 1},
        TaxaReacaoPowerlaw({subst_a: 1, subst_b: 1}, 0.34),
        False
    )

    r3 = Reacao(
        "R3",
        "Reação 3",
        {subst_a: 1, subst_c: 1},
        {subst_e: 1},
        TaxaReacaoPowerlaw({subst_a: 1, subst_c: 1}, 0.7),
        False
    )

    reator = ReatorBatelada(
        {subst_a: 1},
        [r1, r2, r3],
        OpcoesDeSimulacao((0,300), t_eval=np.linspace(0, 20, 1000))
    )

    gera_grafico(reator, gera_dados=[subst_a, subst_e])


def ex_3_otimizacao():
    subst_a = Substancia("A", "A")
    subst_b = Substancia("B", "B")
    subst_c = Substancia("C", "C")
    subst_d = Substancia("D", "D")
    subst_e = Substancia("E", "E")

    r1 = Reacao("R1", 
                "Reação fotólise", 
                {subst_a: 1}, 
                {subst_b: 1, subst_c: 1},
                TaxaReacaoFotolise(subst_a, {"phi": 0.09, "Ep0": 1, "epsilon": 1, 
                                             "Abs": 1, "l": None, "l_reator": 1}),
                True
    )

    r2 = Reacao(
        "R2",
        "Reação 2",
        {subst_a: 1, subst_b: 1},
        {subst_d: 1},
        TaxaReacaoPowerlaw({subst_a: 1, subst_b: 1}, None),
        False
    )

    r3 = Reacao(
        "R3",
        "Reação 3",
        {subst_a: 1, subst_d: 1},
        {subst_e: 1},
        TaxaReacaoPowerlaw({subst_a: 1, subst_d: 1}, None),
        False
    )

    reator = ReatorBatelada(
        {subst_a: 1},
        [r1, r2, r3],
        OpcoesDeSimulacao((0,20), max_step=0.05)
    )

    pts_exp = {subst_e: {0.0: 0.0, 1.981981981981982: 0.01202595984628, 3.983983983983984: 0.04944958425057544, 5.985985985985986: 0.09104332089427424, 7.987987987987988: 0.12542680111801383, 9.98998998998999: 0.15085852570404154, 11.991991991991991: 0.16854592850623157, 13.993993993993994: 0.1805207385991573, 15.995995995995996: 0.1884811323167884, 17.997997997998: 0.1937587478081469, 20.0: 0.1972047565025434}}
    lst_metodos = [None, 'Nelder-Mead', 'Powell', 'CG', 'BFGS', 'Newton-CG', 'L-BFGS-B', 'TNC', 'COBYQA', 'SLSQP', 'trust-constr', 'dogleg', 'trust-ncg', 'trust-exact', 'trust-krylov']
    resultados = {}
    for opt in lst_metodos:
        try:
            print(f"Tentando {opt}")
            otimizador = Otimizador(
                reator,
                pts_exp,
                OpcoesDeOtimizacao(None, [r1, r2, r3], None, [0.4, 0.4, 0.4], 0.000000001)
            )

            result = otimizador.otimizar()

            resultados[opt] = otimizador.get_dados_para_grafico()[subst_a]

            print(f"Método: {opt}")
            print(result)
            print()
        except Exception as e:
            print(e)

    # Cria a figura e o eixo
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plota dados exp
    t_exp = list(pts_exp[subst_e].keys())
    y_exp = [pts_exp[subst_e][t] for t in t_exp]
    ax.plot(t_exp, y_exp, label="Dados Experimentais")

    # Plota evolução pra cada método
    for metodo, result in resultados.items():
        metodo = "Auto" if metodo is None else metodo
        #print(result)
        t = list(result.keys())
        y = [result[tempo] for tempo in t]
        ax.plot(t, y, label=metodo)

    ax.set_xlabel("Tempo (s)")
    ax.set_ylabel("Concentração (mol/L)")
    ax.set_title("Evolução das concentrações no reator batelada")
    ax.legend()
    ax.grid(True)

    plt.show()


if __name__ == "__main__":
    ex_3_otimizacao()