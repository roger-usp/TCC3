from dataclasses import dataclass
from abc import ABC, abstractmethod
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import copy as copy
import pandas as pd

import scipy.interpolate
import scipy.optimize


"""
2025-03-17
V3: essa versão está braba demais

To-do: integrar c front end

2025-04-05
Pré-Integração feita

V4: Output de dados de gráfico

V5: Integração feita a princípio
"""


class InputError(Exception):
    pass


@dataclass(frozen=True)
class Substancia:
    name: str
    id: str
    molar_mass: float | None = None  # g/mol


class TaxaReacao(ABC):
    @abstractmethod
    def calcula_taxa(self):
        pass


@dataclass
class TaxaReacaoPowerlaw(TaxaReacao):
    # Ordem da reação com relação a cada reagente
    # Exemplo: {A: 1, B: 2} significa que r = k * [A]^1 * [B]^2
    ordem: dict[Substancia: float]
    k: float | None = None  # Constante cinética k da reação na temperatura especificada

    def calcula_taxa(self, concentracoes: dict[Substancia: float]):  # Cálculo de r = k * Π (C_substancia^ordem_substancia)
        if self.k is None:
            raise ValueError("Constante de reação (k) não definida.")

        r = self.k
        for substancia, n in self.ordem.items():
            n = 0 if n is None else n
            r = r * concentracoes[substancia] ** n
        
        return r
    

@dataclass
class TaxaReacaoFotolise(TaxaReacao):
    reagente: Substancia
    dados_fotolise: dict[str: float | None]

    def calcula_taxa(self, concentracoes: dict[Substancia: float]) -> float:  
        parametros_fotolise = len(self.dados_fotolise.keys())
        parametros_especificados = [val is not None for val in self.dados_fotolise.values()].count(True)
        parametros_n_definidos = parametros_fotolise - parametros_especificados

        if parametros_n_definidos == 0:  # Calcula a taxa
            if not isinstance(self.reagente, Substancia):
                raise ValueError(f"O reagente da reação de fotólise '{self.name}' não foi identificado. Verifique se a reação possui apenas um reagente").__traceback__()
            
            phi, Ep0, epsilon, absorbancia, l, l_reator = (self.dados_fotolise[v] for v in ["phi", "Ep0", "epsilon", "Abs", "l", "l_reator"])
            conc = concentracoes[self.reagente]
            return phi * Ep0 * (epsilon * conc * l / absorbancia) * (1 - 10 ** (-absorbancia * l_reator / l))
        
        elif parametros_n_definidos > 0:
            raise ValueError(f"""Não foi possível calcular a taxa da reação '{self.name}'. Verifique se todos os parâmetros estão preenchidos: 
                             {phi=}, {Ep0=}, {epsilon=}, {absorbancia=}, {l=}, {l_reator=}, {conc=}""")


@dataclass
class Reacao:
    id: str
    name: str
    reagentes: dict[Substancia: float]  # Estequiometria dos reagentes
    produtos: dict[Substancia: float]  # Estequiometria dos produtos
    taxa: TaxaReacao | None = None
    is_fotolise: bool = False

    def calcula_taxa(self, concentracoes: dict[Substancia: float]) -> float:  # Chama o método dentro da taxa de reação para calcular r
        return self.taxa.calcula_taxa(concentracoes)





"""
Classes auxiliares
"""
class OpcoesDeSimulacao:  #TODO: Mudar o nome dps? achei meio genérico
    """Encapsula todas as opções que possam ser usadas em um simulador aqui"""
    def __init__(self, intervalo_t, **ivp_kwargs):
        self.intervalo_t = intervalo_t
        self.ivp_kwargs = ivp_kwargs


class OpcoesDeOtimizacao:
    def __init__(self, str_metodo_minimizacao: str | None, 
                 lista_reacoes_a_otimizar: list[Reacao], 
                 bounds: list[tuple] | None, 
                 chute_inicial: list[float], 
                 tolerancia: float):
        
        self.str_metodo_minimizacao = None if str_metodo_minimizacao == "Automático" else str_metodo_minimizacao
        self.lista_reacoes_a_otimizar = lista_reacoes_a_otimizar
        self.bounds = bounds
        self.chute_inicial = 0 if chute_inicial is None else chute_inicial
        self.tolerancia = 10**-5 if tolerancia is None else tolerancia




"""
Classes que representam reatores
"""
class ReatorBatelada:
    def __init__(self,
                dict_concentracoes_iniciais: dict[Substancia: float],  # mol/L
                reacoes: list[Reacao],
                opcoes_simulacao: OpcoesDeSimulacao):
        self.dict_concentracoes_iniciais = dict_concentracoes_iniciais
        self.reacoes = reacoes
        self.opcoes_simulacao = opcoes_simulacao

        # Adiciona os compostos não presentes nas concentrações iniciais
        # TODO: Não sei se essa verificação deveria ser feita no front end, ou no back end
        for reacao in self.reacoes:  # Itera sobre os componentes de cada reação
            reagentes = list(reacao.reagentes.keys())
            produtos = list(reacao.produtos.keys())

            for substancia in reagentes + produtos:
                if substancia not in self.dict_concentracoes_iniciais.keys():  # Caso não esteja presente, coloca concentração inicial como 0
                    self.dict_concentracoes_iniciais[substancia] = 0

        # Cria as matrizes de substância e concentração para cálculo das EDOs
        self.matriz_substancias = list(self.dict_concentracoes_iniciais.keys())
        self.concentracoes_iniciais = [self.dict_concentracoes_iniciais[i] for i in self.matriz_substancias]

    
    # Retorna uma matriz da derivada das concentrações
    # TODO: Colocar checks de erro para garantir que a matriz contem todas as concentrações
    def sistema_derivada_conc(self, t, y):  # A variável 't' não é usada aqui, mas o simulador precisa que ela seja declarada
        dydt_dict = {subst: 0 for subst in self.matriz_substancias}
        dict_conc = dict(zip(self.matriz_substancias, y))

        for reacao in self.reacoes:
            taxa_reacao = reacao.calcula_taxa(dict_conc)  # Calcula a taxa da reação

            # Subtrai o consumo de cada componente que é reagente
            for comp in reacao.reagentes.keys():
                dydt_dict[comp] -= taxa_reacao * reacao.reagentes[comp]
            
            # Adiciona a geração de cada componente que é produto
            for comp in reacao.produtos.keys():
                dydt_dict[comp] += taxa_reacao * reacao.produtos[comp]

        # Converte o dicionário de derivadas em uma lista
        dydt = [dydt_dict[i] for i in self.matriz_substancias]
        return dydt
    
    # Simula a concentração do reator ao longo do tempo
    def simular(self, **kwargs):
        combined_kwargs = {**self.opcoes_simulacao.ivp_kwargs, **kwargs}  # Combina os KWargs dando prioridade pros passados no método diretamente
        return scipy.integrate.solve_ivp(
                                        self.sistema_derivada_conc, 
                                        self.opcoes_simulacao.intervalo_t, 
                                        self.concentracoes_iniciais, 
                                        **combined_kwargs)


class Otimizador:
    def __init__(self, reator:ReatorBatelada, 
                 dados_exp:dict[Substancia: dict[float:float]], 
                 opcoes_otimizacao: OpcoesDeOtimizacao):
        """dados_exp: dict[tempo:concentracao]"""
        self.reator = reator
        self.reacao_a_otimizar = opcoes_otimizacao.lista_reacoes_a_otimizar
        self.dados_exp = dados_exp
        self.substancia_interesse = list(dados_exp.keys())
        self.opcoes_otimizacao = opcoes_otimizacao
        self.ultima_iteracao = None

    def simula_k_novo(self, k: list[float]):
        # Cria cópias do reator e da reação, pra não mudar os dados originais (é um pouco pedante isso, mas é melhor ter e n precisar do que precisar e n ter)
        reator = copy.deepcopy(self.reator)  # Possível fonte de erro
        reacao_a_otimizar = copy.deepcopy(self.reacao_a_otimizar)

        # Substitui as constantes k nas reações a otimizar
        for i in range(len(reacao_a_otimizar)):
            if reacao_a_otimizar[i].is_fotolise:  # Substitui constante na reação de Fotólise
                # Encontra qual fator é 'None' e substitui ele
                for key, value in reacao_a_otimizar[i].taxa.dados_fotolise.items():
                    if value is None:
                        reacao_a_otimizar[i].taxa.dados_fotolise[key] = k[i]
                        break

            else:  # Substitui constante na reação Powerlaw
                reacao_a_otimizar[i].taxa.k = k[i]
        
            # Coloca as reações modificadas no reator
            for j in range(len(reator.reacoes)):
                reacao = reator.reacoes[j]
                if reacao.name == reacao_a_otimizar[i].name:
                    reator.reacoes[j] = reacao_a_otimizar[i]

        # Simula o reator
        resultado = reator.simular()
        return resultado

    # Retorna o residual de uma simulação em relação aos dados experimentais
    def residual(self, resultado) -> float:
        # Pega os dados experimentais de interesse
        dados_experimentais = self.dados_exp
        dados_experimentais = {subst: dados if subst in self.substancia_interesse else None for subst, dados in dados_experimentais.items()}

        # Pega os resultados correspondentes aos dados experimentais de interesse
        tempos_resultado = resultado.t
        concentracoes_resultado = resultado.y

        matriz_substancias = self.reator.matriz_substancias

        # Cria um dicionário com as concentrações a partir da matriz de substancias do reator
        concentracoes_resultado = {matriz_substancias[i]: concentracoes_resultado[i] for i in range(len(concentracoes_resultado))}

        # Converte a lista de concentrações em um dicionário {tempo:concentração}
        for subst in concentracoes_resultado.keys():
            dict_tempo_conc = {tempos_resultado[i]: concentracoes_resultado[subst][i] for i in range(len(tempos_resultado))}
            concentracoes_resultado[subst] = dict_tempo_conc

        # Filtra o resultado para as substancias de interesse
        concentracoes_resultado = {subst: dados if subst in self.substancia_interesse else None for subst, dados in concentracoes_resultado.items()}


        # Calcula os residuais
        residual_quadrado = 0

        for substancia in self.substancia_interesse: 
            resultado_experimental = dados_experimentais[substancia]
            dados_modelo = concentracoes_resultado[substancia]
            x = list(dados_modelo.keys())
            y = [dados_modelo[i] for i in x]
            # Criando interpolação dos resultados e ajustando para os tempos do experimento
            interpolador = scipy.interpolate.interp1d(x=x, y=[y])
            resultado_modelo = interpolador(list(resultado_experimental.keys()))  # Resultado do modelo interpolado em t experimental

            # Calcula os residuais
            for i in range(len(resultado_experimental.keys())):
                diferenca = list(resultado_experimental.values())[i] - resultado_modelo[0][i]
                residual_quadrado += diferenca ** 2

        return residual_quadrado

    # Recebe o valor k da reação a ser otimizada, e retorna o erro em relação aos dados experimentais
    # Na prática, ele só junta as funções em um único lugar pra passar para o algoritmo de otimização
    def funcao_de_otimizacao(self, k) -> float:  
        resultado = self.simula_k_novo(k)

        self.ultima_iteracao = resultado  # Salva a ultima iteracao pra plotar o gráfico
        residual = self.residual(resultado)

        return residual
    
    def otimizar(self, **kwargs):
        x0 = self.opcoes_otimizacao.chute_inicial
        return scipy.optimize.minimize(
            fun=self.funcao_de_otimizacao, 
            x0=x0,
            method=self.opcoes_otimizacao.str_metodo_minimizacao,
            bounds=self.opcoes_otimizacao.bounds,
            tol=self.opcoes_otimizacao.tolerancia, 
            **kwargs) 
        
    def get_dados_para_grafico(self):
        matriz_substancias = self.reator.matriz_substancias
        dados_modelo = {}

        for subst in self.substancia_interesse:
            index_subst = matriz_substancias.index(subst)
            dados_modelo[subst] = dict(zip(self.ultima_iteracao.t, self.ultima_iteracao.y[index_subst]))

        return dados_modelo


"""
Integracao com o front end
"""
def cria_lista_substancia(df_substancias: pd.DataFrame) -> list[Substancia]:
    # Tratamento de dados faltantes
    df_substancias = df_substancias.replace('', None)
    df_substancias = df_substancias.dropna(how='all')

    # Verifica se tem substâncias no dataframe
    if df_substancias.empty:
        raise InputError("Nenhuma substância foi adicionada")
    
    # Verifica se tem linhas com id None ou vazio
    df_linhas_semi_preenchidas = df_substancias.loc[df_substancias["id"].isna()]

    if not df_linhas_semi_preenchidas.empty:  # Se tiver, sobe erro
        raise InputError(f"Todas as substâncias precisam ter um ID")

    # Verifica se algum ID está sendo usado por mais de uma substância
    df_linhas_repetidas = df_substancias[df_substancias.duplicated(subset=["id"], keep=False)]

    if not df_linhas_repetidas.empty:
        raise InputError(f"Existem {len(df_linhas_repetidas)} substâncias com o ID duplicado")

    list_substancias = []
    for idx, row in df_substancias.iterrows():
        subst = Substancia(row["name"], row["id"], row["molar_mass"])
        list_substancias.append(subst)

    return list_substancias


def cria_lista_reacoes(reactions_df: pd.DataFrame, list_substancias: list[Substancia]) -> list[Reacao]:
    """
    # Tratamento de dados faltantes
    reactions_df = reactions_df.replace('', None)
    reactions_df = reactions_df.dropna(how='all')
    
    # Verifica se tem reações no dataframe
    if reactions_df.empty:
        raise InputError("Nenhuma reação foi adicionada")
    
    # Verifica se tem linhas com id None ou vazio
    df_linhas_semi_preenchidas = reactions_df.loc[reactions_df["reaction_id"].isna()]

    if not df_linhas_semi_preenchidas.empty:  # Se tiver, sobe erro
        raise InputError(f"Todas as reações precisam ter um ID")

    # Verifica se algum ID está sendo usado por mais de uma substância
    df_linhas_repetidas = reactions_df[reactions_df.duplicated(subset=["reaction_id"], keep=False)]

    if not df_linhas_repetidas.empty:
        raise InputError(f"Existem {len(df_linhas_repetidas)} reações com o ID duplicado")
    
    # Verifica se pelo menos uma reação está ligada
    df_linhas_ativas = reactions_df.loc[reactions_df["selected"] is True]

    if df_linhas_ativas.empty:
        raise InputError(f"Não existe nenhuma reação ativa")
    """

    list_reacoes = []
    for idx, row in reactions_df.iterrows():
        if row["selected"] is False:  # Retira reações não selecionadas
            continue

        de_para_substancias = {subst.id: subst for subst in list_substancias}  # Dict pra encontrar o objeto substancia a partir do ID
        reagentes = {de_para_substancias[comp_id]: abs(coef) for comp_id, coef in row["stoic_coefs"].items() if coef < 0}
        produtos = {de_para_substancias[comp_id]: abs(coef) for comp_id, coef in row["stoic_coefs"].items() if coef > 0}

        # Declara a taxa de reação
        if row["kinectic_type"] == "Fotocatalítica":
            taxa_reacao = TaxaReacaoFotolise(
                reagente=list(reagentes.keys())[0],
                dados_fotolise = {value: row[value] for value in ["phi", "Ep0", "epsilon", "Abs", "l", "l_reator"]}
            )

        if row["kinectic_type"] == "Powerlaw":
            taxa_reacao = TaxaReacaoPowerlaw(
                ordem={de_para_substancias[comp_id]: coef for comp_id, coef in row["exponents"].items()},
                k = row["powerlaw_k"]
            )

        # Declara a reação
        reacao = Reacao(
            id = row["reaction_id"],
            name = row["reaction_str"],
            reagentes = reagentes,
            produtos = produtos,
            is_fotolise = (row["kinectic_type"] == "Fotocatalítica"),
            taxa = taxa_reacao
        )

        list_reacoes.append(reacao)

    return list_reacoes
            

def cria_opcoes_simulacao(str_metodo: str, passo_max: float | None, time_bound: list[float, float]):
    #print("time bound = ", time_bound)
    passo_max = np.inf if passo_max is None else passo_max
    
    return OpcoesDeSimulacao(time_bound, max_step=passo_max, method=str_metodo)


def cria_opcoes_otimizacao(str_metodo: str, df_reacoes: pd.DataFrame, tolerancia, lista_todas_reacoes):
    # Colunas: "reaction_str", "initial_guess", "bound_min", "bound_max"
    de_para_reacoes = {react.name: react for react in lista_todas_reacoes}  # Dict pra encontrar o objeto reacao a partir do ID
    
    lista_reacoes_otimizar = []
    lista_chutes_iniciais = []
    lista_bounds=[]

    for idx, row in df_reacoes.iterrows():
        reacao = de_para_reacoes[row["reaction_str"]]
        lista_reacoes_otimizar.append(reacao)

        k0 = row["initial_guess"]
        lista_chutes_iniciais.append(k0)

        bounds = (row["bound_min"], row["bound_max"])
        lista_bounds.append(bounds)

    return OpcoesDeOtimizacao(
        str_metodo_minimizacao=str_metodo,
        lista_reacoes_a_otimizar=lista_reacoes_otimizar,
        bounds=lista_bounds,
        chute_inicial=lista_chutes_iniciais,
        tolerancia=tolerancia
    )


def cria_otimizador(reator:ReatorBatelada, 
                df_dados_exp:pd.DataFrame, 
                opcoes_otimizacao: OpcoesDeOtimizacao, 
                molar_or_mass_base: str, 
                volume_reator: float, 
                list_substancias: list[Substancia]):
    #colunas: time, conc, comp_id

    # Tratamento de dados experimentais faltantes
    de_para_substancias = {subst.id: subst for subst in list_substancias}
    df_dados_exp = df_dados_exp.replace('', None)
    df_dados_exp = df_dados_exp.dropna(how='all')

    # Verifica se alguma das colunas está com um valor faltante, e se tiver sobe erro
    df_linhas_semi_preenchidas = df_dados_exp[df_dados_exp.isna().sum(axis=1).between(1, 2)]

    if not df_linhas_semi_preenchidas.empty:  # Se tiver, sobe erro
        raise ValueError(f"Dados experimentais incompletos detectados")


    # Se a base for mássica, converte pra molar
    if molar_or_mass_base == "Base mássica":
        for idx, row in df_dados_exp.iterrows():
            row["conc"] = row["conc"] / de_para_substancias[row["comp_id"]].molar_mass


    # Muda o dataframe para um dicionário de dicionários
    dict[Substancia: dict[float: float]]
    dados_exp = {}
    if list(df_dados_exp.comp_id.unique()) == [None]:
        raise ValueError("Os dados experimentais estão sem componente especificado")

    for comp_id in df_dados_exp.comp_id.unique():
        dados_exp_comp = df_dados_exp.loc[df_dados_exp["comp_id"] == comp_id]

        dados_exp_comp_dict = {}
        for idx, row in dados_exp_comp.iterrows():
            dados_exp_comp_dict[float(row["time"])] = float(row["conc"])
    
        dados_exp[de_para_substancias[comp_id]] = dados_exp_comp_dict
    
    return Otimizador(reator, dados_exp, opcoes_otimizacao)


def cria_reator_batelada(df_dados_experimentais: pd.DataFrame, molar_or_mass_base: str, volume_reator: float | None, list_reacoes: list[Reacao], 
                         list_substancias: list[Substancia], opcoes_simulacao: OpcoesDeSimulacao) -> ReatorBatelada:
    
    de_para_substancias = {subst.id: subst for subst in list_substancias}  # Dict pra encontrar o objeto substancia a partir do ID

    if molar_or_mass_base == "Base mássica":
        # TODO: Converter as concentrações de base mássica pra base molar
        for idx, row in df_dados_experimentais.iterrows():
            row["conc"] = row["conc"] / de_para_substancias[row["comp_id"]].molar_mass

    # Monta o dicionário de concentrações iniciais
    dict_conc_iniciais = {}

    for idx, row in df_dados_experimentais.iterrows():
        if row["time"] == 0:
            obj_substancia = de_para_substancias[row["comp_id"]]
            dict_conc_iniciais[obj_substancia] = row["conc"]

    # Gera o reator batelada
    reator = ReatorBatelada(
        dict_concentracoes_iniciais=dict_conc_iniciais,
        reacoes=list_reacoes,
        opcoes_simulacao=opcoes_simulacao
    )

    return reator
