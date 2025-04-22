import streamlit as st
from streamlit import session_state as sst

import numpy as np
import pandas as pd
import uuid

import json

from .backend import *

# Initialize
if "reactions_df" not in sst:
    sst.reactions_df = pd.DataFrame(columns=["selected", "reaction_str", "stoic_coefs", "exponents", "kinectic_type", "powerlaw_k", "phi", "Ep0", "epsilon", "Abs", "l", "l_reator", "reaction_id"])

if "reaction_id" not in sst:
    sst.reaction_id = None

reset_keys_list = ["create_reaction_df_key","powerlaw_k_key","phi_key","Ep0_key","epsilon_key","Abs_key","l_key","l_reator_key", "reactions_to_be_deleted_key"]
for reset_key in reset_keys_list:
    if reset_key not in sst:
        sst[reset_key] = str(uuid.uuid4())


def get_stoic_dict():
    # PS: só é possível colocar valores positivos de coeficiente, por isso tem a opção "produto/reagente"
    stoic_dict = dict()
    for idx, row in sst.create_reaction_df.iterrows():
        if row["comp_id"] not in stoic_dict.keys():
            if row["product_reagent"] == "Reagente" and row["coefficient"] is not None:
                stoic_dict[row["comp_id"]] = - row["coefficient"]
            elif row["product_reagent"] == "Produto" and row["coefficient"] is not None:
                stoic_dict[row["comp_id"]] = row["coefficient"]
    return stoic_dict


def get_exp_dict():
    exp_dict = dict()
    for idx, row in sst.create_reaction_df.iterrows():
        if row["comp_id"] not in exp_dict.keys():
            exp_dict[row["comp_id"]] = row["exponent"]
    return exp_dict

def get_reaction_str(stoic_dict):
    reagents = []
    products = []
    for comp_id, coef in stoic_dict.items():
        if coef < 0:
            reagents.append(f"{-coef} {comp_id}")
        elif coef > 0:
            products.append(f"{coef} {comp_id}")
    
    return "+".join(reagents) + " → " + "+".join(products)

# Callbacks
def add_reaction_clicked():
    # add to reaction_df
    new_data = {
        "stoic_coefs": sst.stoic_dict,
        "exponents": sst.exp_dict,
        "kinectic_type": sst.kinect_type,
        "selected":True,
        "reaction_str": get_reaction_str(sst.stoic_dict),
        "reaction_id": sst.reaction_id
    }
    
    if sst.kinect_type == "Powerlaw":
        new_data["powerlaw_k"] = sst.powerlaw_k
        for const in ["phi", "Ep0", "epsilon", "Abs", "l", "l_reator"]:
            new_data[const] = None
    
    elif sst.kinect_type == "Fotocatalítica":
        new_data["powerlaw_k"] = None
        for const in ["phi", "Ep0", "epsilon", "Abs", "l", "l_reator"]:
            new_data[const] = sst[const]
    
    sst.reactions_df = pd.concat([sst.reactions_df, pd.DataFrame([new_data])], ignore_index=True)


    # reset add_reaction data_editor
    sst.create_reaction_df_key = str(uuid.uuid4())
    
    # reset powerlaw
    sst.powerlaw_k_key =str(uuid.uuid4())

    # reset fotocatalitica
    sst.phi_key = str(uuid.uuid4())
    sst.Ep0_key = str(uuid.uuid4())
    sst.epsilon_key = str(uuid.uuid4())
    sst.Abs_key = str(uuid.uuid4())
    sst.l_key = str(uuid.uuid4())
    sst.l_reator_key = str(uuid.uuid4())

    sst.reaction_id = None


def delete_reaction_clicked():
    sst.reactions_df = sst.reactions_df.loc[sst.reactions_df.reaction_id != reaction_id_to_delete]
    sst.reactions_to_be_deleted_key = uuid.uuid4()



# Components
st.header("Componentes")

components_df = st.data_editor(
    pd.DataFrame(columns=["id", "name", "molar_mass"]),
    column_config={
        "id": st.column_config.TextColumn("ID"),
        "name": st.column_config.TextColumn("Nome"),
        "molar_mass": st.column_config.NumberColumn("Massa molar (g/mol)")
    },
    num_rows="dynamic",
    hide_index=True
)

# Reactions
st.header("Reações")
st.write("Monte a reação")

sst.reaction_id = st.text_input("ID Reação", value=sst.reaction_id)

# Isso aqui deixa uma brecha pro cara selecionar o mesmo componente 2 vezes
sst.create_reaction_df = st.data_editor(
    pd.DataFrame(columns=["comp_id", "product_reagent", "coefficient", "exponent"]),
    column_config={
        "comp_id": st.column_config.SelectboxColumn("ID componente", options=components_df["id"].unique(), required=True),
        "product_reagent": st.column_config.SelectboxColumn("Produto/Reagente", options=["Reagente","Produto"], required=True),
        "coefficient": st.column_config.NumberColumn("Coeficiente estequiométrico", min_value=0, required=True),
        "exponent": st.column_config.NumberColumn("Expoente") 
    },
    num_rows="dynamic",
    hide_index=True,
    key = sst.create_reaction_df_key
)

sst.stoic_dict = get_stoic_dict()
#print("stoic_dict", sst.stoic_dict)


sst.exp_dict = get_exp_dict()
sst.kinect_type = st.selectbox(
    "Selecione o tipo da cinética",
    options=["Powerlaw", "Fotocatalítica"],
    index=None,
    placeholder="Selecione",
    key="kinect_type_selectbox"
)

st.write("*Constante não preenchida indica que ela deverá ser encontrada com regressão dos dados experimentais")


# Reaction ID


# Reactions - Powerlaw
if sst.kinect_type == "Powerlaw":
    comp_exp_dict = dict() # {comp_id: exponent}

    for idx, row in sst.create_reaction_df.iterrows():
        if row["exponent"] is None:
            continue
        comp_exp_dict[row["comp_id"]] = row["exponent"]
    
    powerlaw_col1, powerlaw_col2 = st.columns(2)

    with powerlaw_col1:
        sst.powerlaw_k = st.number_input(
            "Constante k",
            value=None,
            placeholder=None,
            key=sst.powerlaw_k_key
        )
    
    with powerlaw_col2:
        # equation display
        equation_str = "k"
        for comp_id, exp in comp_exp_dict.items():
            if exp != 0:
                equation_str += "(C_{" + str(comp_id) + "})^{" + str(exp) + "}"
        
        if equation_str != "k":
            # Tem alguma equação pra mostrar
            st.latex(equation_str)




# Reactions - fotocatalítica
elif sst.kinect_type == "Fotocatalítica":
    st.write("Essa cinética permite apenas um reagente")
    comp_id = "C"
    latex_equation = fr"""
        \phi \cdot E_{{p,0}} \frac{{\varepsilon \cdot {comp_id} \cdot l}}{{Abs}} 
        \left(1 - 10^{{ \frac{{Abs}}{{l}} \cdot l_{{reator}} }}\right)
        """
    st.latex(latex_equation)

    fotocat_col1, fotocat_col2 = st.columns(2)
    with fotocat_col1:
        sst.phi = st.number_input("$\phi$",value=None,placeholder=None,key=sst.phi_key)
        sst.Ep0 = st.number_input("$E_{p,0}$",value=None,placeholder=None,key=sst.Ep0_key)
        sst.epsilon = st.number_input("$\\varepsilon$",value=None,placeholder=None,key=sst.epsilon_key)
    
    with fotocat_col2:
        sst.Abs = st.number_input("$Abs$",value=None,placeholder=None,key=sst.Abs_key)
        sst.l = st.number_input("$l$",value=None,placeholder=None,key=sst.l_key)
        sst.l_reator = st.number_input("$l_{reator}$",value=None,placeholder=None,key=sst.l_reator_key)

    


add_reaction = st.button("Adicionar Reação", on_click=add_reaction_clicked)


st.divider()

reactions_df_col, delete_reaction_col = st.columns(2)
with delete_reaction_col:
    reaction_id_to_delete = st.selectbox(
        "Selecione uma reação para deletar",
        options=sst.reactions_df.reaction_id,
        key=sst.reactions_to_be_deleted_key
    )
    delete_reaction_btn = st.button("Deletar reação selecionada", on_click=delete_reaction_clicked)
    
    if delete_reaction_btn:
        delete_reaction_clicked()



with reactions_df_col:
    st.dataframe(
        sst.reactions_df,
        column_order=["reaction_id", "reaction_str", "kinectic_type"],
        column_config={},
        hide_index=True
    )

sst.active_reactions = st.multiselect(
    "Reações ativas",
    options=sst.reactions_df.reaction_id,
    default=sst.reactions_df.reaction_id
)


#sst.reactions_df["selected"] = sst.reactions_df["selected"].apply(lambda x: x)


# Reactor
st.header("Reator")
reactor_col1, reactor_col2 = st.columns(2)
with reactor_col1:
    molar_or_mass_base = st.selectbox("Base", options=["Base molar", "Base mássica"])

with reactor_col2:
    reactor_volume = st.number_input("Volume do reator (L)")


st.write("É possível colar os dados experimentais direto do Excel, Google Sheets ou similar.")
st.write("O ID do componente deve ser o mesmo fornecido na primeira tabela")
st.write("A concentração do componente no tempo zero deve ser fornecida")

st.subheader("Dados experimentais")
sst.exp_data = st.data_editor(
    pd.DataFrame(columns=["time", "conc", "comp_id"]),
    hide_index=True,
    column_config={
        "time": "Tempo (min)",
        "conc": "Concentração (mmol/L)" if molar_or_mass_base == "Base molar" else "Concentração (mg/L)",
        "comp_id": "ID componente"
    },
    num_rows="dynamic"
)


sst.exp_data_float = sst.exp_data.copy()
for exp_data_col in ["time", "conc"]:
    # coloca o separador como ponto ao invés de vírgula
    sst.exp_data_float[exp_data_col] = pd.to_numeric(
        sst.exp_data_float[exp_data_col].apply(lambda s: str(s).replace(",", ".")),
        errors="coerce"
    )

sst.exp_data_float["comp_id"] = sst.exp_data["comp_id"]
#print("exp_data_float:\n",exp_data_float)



st.subheader("Parâmetros de simulação")
# st.write("Perguntar ao capitão em caso de dúvidas")
simul_ivp_method = st.selectbox("Método de resolução da EDO", ["RK45", "RK23", "DOP853", "Radau", "BDF", "LSODA"])
max_step = st.number_input("Passo máximo de integração", value=None, format="%.10g")


# tempos que o usuário passar [0,0] se nada for preenchido
exp_data_time_bound = [0, sst.exp_data_float.copy().time.fillna(0).astype(float).max()]


st.subheader("Parâmetros de otimização")
sst.initial_guess_bound_df = st.data_editor(
sst.reactions_df.copy().assign(initial_guess=None, bound_min=None, bound_max=None),
    column_order=["reaction_str", "initial_guess", "bound_min", "bound_max"],
    column_config={
        "reaction_str": st.column_config.Column("Reação",disabled=True),
        "initial_guess": st.column_config.NumberColumn("Estimativa inicial"),
        "bound_min": st.column_config.NumberColumn("Valor mínimo"),
        "bound_max": st.column_config.NumberColumn("Valor máximo")
    }
)

optimization_method = st.selectbox("Método de otimização", ["Automático", 'Nelder-Mead', 'Powell', 'CG', 'BFGS', 'Newton-CG', 'L-BFGS-B', 'TNC', 'COBYLA', 'COBYQA', 'SLSQP', 'trust-constr', 'dogleg', 'trust-ncg', 'trust-exact', 'trust-krylov'])

tolerance = st.number_input("Tolerância", value=None, format="%.10g")
max_iter = st.number_input("Número máximo de iterações", value=None, format="%.10g")





def calculate_btn_callback():
    try:
        list_substancias = cria_lista_substancia(components_df)
        #print(f"Lista de substâncias: {list_substancias}")
        list_reacoes = cria_lista_reacoes(sst.reactions_df, list_substancias)
        #print(f"Lista de reações: {list_reacoes}")

        #print("reactions_df\n", sst.reactions_df)
        opcoes_simulacao = cria_opcoes_simulacao(simul_ivp_method, max_step, exp_data_time_bound)
        #print(f"Opções de simulação: {opcoes_simulacao.__dict__}")

        opcoes_otimizacao = cria_opcoes_otimizacao(optimization_method, sst.initial_guess_bound_df, tolerance, list_reacoes)
        #print(f"Opções de otimização: {opcoes_otimizacao.__dict__}")
        
        reator = cria_reator_batelada(sst.exp_data_float, molar_or_mass_base, reactor_volume, list_reacoes, list_substancias, opcoes_simulacao)
        #print(f"Reator: {reator.__dict__}")
        
        sst.otimizador = cria_otimizador(reator, sst.exp_data, opcoes_otimizacao, molar_or_mass_base, reactor_volume, list_substancias)
        #print(f"Otimizador: {sst.otimizador.__dict__}")
        
        sst.optimization_result = sst.otimizador.otimizar()

        sst.raw_graph_data = sst.otimizador.get_dados_para_grafico() # {subst: {time:conc}}
        #print(sst.graph_data.keys())
    except InputError as e:
        st.warning(str(e))



calculate_btn = st.button("Calcular constantes", on_click=calculate_btn_callback)



if "optimization_result" in sst and "raw_graph_data" in sst:
    st.subheader("Resultados")

    # Resultados (k)
    st.write(f"Número de iterações: {sst.optimization_result['nit']}")
    st.write(f"Mean squared error (MSE): {sst.optimization_result['fun']:.2e}")
    st.dataframe(
        pd.DataFrame(
            {"k": [f"{result:.3E}" for result in sst.optimization_result["x"]],
             "Reação": [obj_reac.name for obj_reac in sst.otimizador.reacao_a_otimizar],
             "Determinante da Metriz Hessiana": [f"{result:.3E}" for result in sst.optimization_result["jac"]]
            }
        )
    )


    # Gráfico
    sst.optimized_subst_ids = [subst.id for subst in sst.raw_graph_data.keys()]
    sst.graph_data = {subst.id: sst.raw_graph_data[subst] for subst in sst.raw_graph_data.keys()}
    sst.graph_data = {subst_id: pd.DataFrame({"time": list(sst.graph_data[subst_id].keys()), "conc": list(sst.graph_data[subst_id].values())}) for subst_id in sst.graph_data.keys()}
    
    sst.graph_data = {subst_id: sst.graph_data[subst_id].sort_values(by="time", ascending=True) for subst_id in sst.graph_data.keys()}
    
    sst.graph_selected_substs = st.multiselect("Substâncias selecionadas", options=sst.optimized_subst_ids)

    sst.to_plot_dict = {subst_id: sst.graph_data[subst_id] for subst_id in sst.graph_selected_substs}
    
    if sst.to_plot_dict:
        sst.fig, sst.ax = plt.subplots()
        for subst_id, df in sst.to_plot_dict.items():
            sst.ax.plot("time", "conc", data=df, label=f"Curva {subst_id}")
            sst.ax.scatter(
                "time",
                "conc",
                data=sst.exp_data_float.loc[sst.exp_data_float["comp_id"] == subst_id],
                label=f"Dados exp. {subst_id}"
            )
        sst.ax.legend()
        sst.ax.set_xlabel("Tempo")
        sst.ax.set_ylabel("Concentração")

        st.pyplot(sst.fig)




