import streamlit as st
import pandas as pd
import numpy as np
from scipy.optimize import fsolve
import plotly.graph_objects as go
import plotly.express as px

def calc_srk(T,P):
    T = T[0]
    R = 8.31446261815324
    m = 0.48508 + 1.55171*omega - 0.15613*omega**2
    alpha = (1+m*(1-np.sqrt(T/Tc)))**2
    #a = 0.42748*R**2*Tc**2/(Pc*10**5)
    a = a_fac*R**2*Tc**2/(Pc*10**5)
    #b = 0.086649*R*Tc/(Pc*10**5)
    b = b_fac*R*Tc/(Pc*10**5)
    V_start = R*T/(P*10**5)*10**3
    # sol hat Einheite dm³/mol = m³/kmol
    sol = fsolve(lambda V: (R*T/(V*10**-3-b)-a*alpha/(V*10**-3*(V*10**-3+b))-P*10**5)**2, V_start)
    if EOS == "Ideales Gas":
        return 1/V_start
    if EOS == "SRK":
        return 1/sol[0]
    #return [1/sol[0], 1/V_start]

def prop_database():
    
    if component == "Wasser":
        Tc = 647.1 # K
        Pc = 220.64 # bar
        omega = 0.3443
    if component == "Wasserstoff":
        Tc = 33.145 # K
        Pc = 12.964 # bar
        omega = -0.219 
    if component == "Ethanol":
        Tc = 241.56+273.15 # K
        Pc = 62.68 # bar
        omega = 0.646
    if component == "CO2":
        Tc= 30.978 + 273.15 # K
        Pc = 73.773 # bar
        omega = 0.22394
    return Tc, Pc, omega

def density_plot():
    data = load_data(component)
    
    #data["model_results"], data["ideal_Gas"] = zip(*data.apply(lambda x:calc_srk(x.T+273.15,x.P),axis=1))
    data["model_results"]= data.apply(lambda x:calc_srk(x.T+273.15,x.P),axis=1)
    deviation = np.sqrt(sum((data["Density"]-data["model_results"])**2))

    fig = go.Figure()
    fig.add_scatter(
        x = data["T"],
        y = data["Density"],
        name = "Experimentelle Daten"
    )
    fig.add_scatter(
        x = data["T"],
        y = data["model_results"],
        name = "Modellierte Daten",
        mode = "lines",
        
    )
    fig.update_layout(        
        width=600,
        height=600,
        title="Dichte als Funktion der Temperatur",
        xaxis_title="Temperatur / °C",
        yaxis_title="Dichte / kmol/m³",
        legend_title="Legend",
        font=dict(            
            size=12,            
        )
    )
    return data, deviation, fig
def vle(T,x_1,P):
    x_2 = 1-x_1
    
    gamma1, gamma2, plv1, plv2, x_1, y_1 = gamma_calc(T,x_1)

    Pcalc = x_1*gamma1*plv1+x_2*gamma2*plv2
    return (Pcalc-P)**2

def gamma_calc(T,x_1):
    #1: water
    x_2 = 1-x_1

    #T = 20+273.15
    #P = 1

    alpha = alpha_fac
    a_11 = 0
    a_22 = 0
    b_11 = 0
    b_22 = 0
    a_12 = 3.4578
    a_21 = -0.8009
    b_12 = -586.0809
    b_21 = 246.18

    tau_11 = a_11+b_11/T
    tau_22 = a_22+b_22/T
    tau_12 = a_12+b_12/T
    tau_21 = a_21+b_21/T

    G_11 = np.exp(-alpha*tau_11)
    G_21 = np.exp(-alpha*tau_21)
    G_12 = np.exp(-alpha*tau_12)
    G_22 = np.exp(-alpha*tau_22)

    ln_gamma1=x_2**2*(tau_21*(G_21/(x_1+x_2*G_21))**2+(tau_12*G_12)/(x_2+x_1*G_12)**2)
    ln_gamma2=x_1**2*(tau_12*(G_12/(x_2+x_1*G_12))**2+(tau_21*G_21)/(x_1+x_2*G_21)**2)
    gamma1 = np.exp(ln_gamma1)
    gamma2 = np.exp(ln_gamma2)

    c_11 = 62.13607454
    c_12 = -7258.2
    c_15 = -7.3037
    c_16 = 4.1653E-06
    c_17 = 2
    c_21 = 61.79107454
    c_22 = -7122.3
    c_25 = -7.1424
    c_26 = 2.8853E-06
    c_27 = 2

    ln_plv1 = c_11+c_12/(T)+c_15*np.log(T)+c_16*(T)**c_17
    ln_plv2 = c_21+c_22/(T)+c_25*np.log(T)+c_26*(T)**c_27
    plv1 = np.exp(ln_plv1)
    plv2 = np.exp(ln_plv2)

    Pcalc = x_1*gamma1*plv1+x_2*gamma2*plv2
    y_1 = x_1*gamma1*plv1/Pcalc
    y_2 = x_2*gamma1*plv2/Pcalc
    return gamma1, gamma2, plv1, plv2, x_1, y_1
def vle_plot():
    color_list = px.colors.qualitative.Plotly

    P = 1
    T_start = 105 + 273.15
    x_val = np.linspace(0,1,250)
    vle_model_results = pd.DataFrame(x_val,columns=['x'])
    #vle_model_results['T'] = T_start
    vle_model_results['T']= vle_model_results.apply(lambda x:fsolve(vle,T_start,args=(x.x,P)),axis=1)
    vle_model_results['T'] = vle_model_results.apply(lambda x:x["T"].item(),axis=1)
    vle_model_results['gamma1'],vle_model_results['gamma2'],vle_model_results['plv1'],vle_model_results['plv2'],vle_model_results['x_1'],vle_model_results['y_1'] = gamma_calc(vle_model_results['T'].values,vle_model_results['x'].values)
    lit_data = load_vle_data()
    lit_data['T_calc']= lit_data.apply(lambda x:fsolve(vle,T_start,args=(x.x,P)),axis=1)
    lit_data['T_calc'] = lit_data.apply(lambda x:x["T_calc"].item(),axis=1)
    deviation = np.sqrt(sum((lit_data["T"]-lit_data["T_calc"])**2))
    fig = go.Figure()
    fig.add_scatter(
        x = 1 - vle_model_results['x_1'],
        y = vle_model_results['T'],
        mode = "lines",
        name = "Modellierte Siedelinie",
        #line=dict(color=color_list[1])
    )
    fig.add_scatter(
        x = 1 - vle_model_results['y_1'],
        y = vle_model_results['T'],
        name = "Modellierte Taulinie",
        mode = "lines",
        #line=dict(color=color_list[1])
    )
    fig.add_scatter(
        x = 1 - lit_data['x'],
        y = lit_data['T'],
        mode = "markers",
        name = "Experimentelle Daten Siedelinie",
        line=dict(color=color_list[1])
    )
    fig.add_scatter(
        x = 1 - lit_data['y'],
        y = lit_data['T'],
        name = "Experimentelle Daten Taulinie",
        mode = "markers",
        line=dict(color=color_list[1])
    )
    fig.update_layout(        
        width=600,
        height=600,
        title="VLE Ethanol - Wasser bei 1 bar",
        xaxis_title="Konzentration Ethanol",
        yaxis_title="Tempearatur / K",
        legend_title="Legend",
        font=dict(            
            size=12,            
        )
    )
    # for i, _x_val in enu.merate(x_val):
    #     T_solve = fsolve(vle,T_start,args=(_x_val,P))
    #     vle_model_results[i,'T'] = T_solve[0]
    return vle_model_results, fig, deviation
        
def reset_slider():
    st.session_state.a_slider = 0.42748
    st.session_state.b_slider = 0.086649

@st.cache_data
def load_data(component):
    if component == "Wasser":
        df = pd.read_excel("dens_data_h2o.xlsx")
    if component == "Wasserstoff":
        df = pd.read_excel("dens_data_h2.xlsx")
    if component == "CO2":
        df = pd.read_excel("dens_data_CO2.xlsx")
    return df
@st.cache_data
def load_vle_data():
    df = pd.read_excel("vle_data_h2o_etoh.xlsx")
    return df

st.set_page_config(layout="wide") 
st.title('Hands-On Übung zum Thema Stoffdatenmodelle')
a_fac = 0.4
b_fac = 0.4

with st.sidebar:
    
    st.header("Einstellungen")
    modus = st.selectbox("Berechnungsmodus",('Dichte eines Reinstoffes', 'VLE'))
    if modus == "Dichte eines Reinstoffes":
        component = st.selectbox("Komponente auswählen",('Wasser', 'Wasserstoff','CO2'), on_change=reset_slider)
        EOS = st.selectbox("EOS auswählen",('Ideales Gas', 'SRK'), on_change=reset_slider)
        if EOS == "SRK":
            a_fac = st.number_input("Parameter in Term für a",0.000,1.000,step=0.005, key="a_slider",format="%.4f")
            b_fac = st.number_input("Parameter in Term für b",-1.0,1.0,step=0.005, key="b_slider",format="%.4f")
        #T=st.slider("Temperatur / K",10.0,1000.0, 493.15)
        #P=st.slider("Druck / bar",1.0,100.0, 20.86)
        Tc, Pc, omega = prop_database()
    if modus == "VLE":
        alpha_fac = st.number_input("Parameter in Term für a",-1.000,1.000,value=0.2,step=0.05, key="alpha_slider",format="%.4f")
    #V = calc_srk()
if modus == "Dichte eines Reinstoffes":
    model_data, model_quality, fig = density_plot()
    st.write("Mittlere quadratische Abweichung zwischen Modellergebnissen und experimentellen Daten:",round(model_quality,10))
    st.plotly_chart(fig, theme="streamlit")
    st.dataframe(model_data)
if modus == "VLE":
    
    model_data, fig, model_quality =vle_plot()
    st.write("Mittlere quadratische Abweichung zwischen Modellergebnissen und experimentellen Daten:",round(model_quality,10))
    st.plotly_chart(fig, theme="streamlit")
    st.dataframe(model_data)
