#ESTA VERSIÓN CREA UNA INTERFAZ PARA SUBIR PDFS DE DATOS Y CALCULAR PROP. DE MIE

#Objetivos:
# 1. Lectura de archivo (habrá que fijar el formato)
# 2. Guardar en vectores y posible interpolación de valores
# 3. Cálculo de propiedades de Mie usando PyMieSim

import panel as pn
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
pn.extension('plotly')
pn.extension('tabulator')
import matplotlib.pyplot as plt
import io
from io import StringIO, BytesIO

import sys
from pathlib import Path

base_path = Path(__file__).resolve().parents[1]  # Esto apunta a PyMieSim_local
sys.path.insert(0, str(base_path))

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.units import nanometer, meter, degree, watt, AU, RIU, micrometer

file_input = pn.widgets.FileInput(accept='.csv', name='Select file')
tabla1 = pn.widgets.Tabulator(pd.DataFrame(), name='Uploaded values', height=400, show_index=False)
boton_mostrar = pn.widgets.Button(name='Show table', button_type='primary')
boton_mostrar2 = pn.widgets.Button(name='Calculate efficiency factors', button_type='primary')
boton_mostrar3 = pn.widgets.Button(name='Plot efficiency factors', button_type='primary')
grafica_todo = pn.pane.Plotly(width=600)


#--------------------------------------------------------------------------------------------------
# 1. Widget para subir el archivo
#--------------------------------------------------------------------------------------------------
def leer_archivo(event):
    wavelength = []
    real_n = []
    kappa = []
    if file_input.value is not None:
        df = pd.read_csv(io.BytesIO(file_input.value))
        df.columns = df.columns.str.strip()

        required_cols = ['Wavelength (micrometers)', 'n', 'k']
        if all(col in df.columns for col in required_cols):
            df = df[required_cols].copy()
            df.columns = ['Wavelength', 'n', 'kappa']
            wavelength = df['Wavelength'].to_numpy()
            real_n = df['n'].to_numpy()
            kappa = df['kappa'].to_numpy()
    
    return wavelength, real_n, kappa


file_input.param.watch(leer_archivo, 'value')


def datos_imprimir(event):
    wl, rn, kp = leer_archivo(event)
    df = pd.DataFrame({
        'Wavelength' : wl,
        'n' : rn,
        'k' : kp,
    })

    return df


#---------------------------------------------
# Siguientes funciones: n_complejo y PMS
#---------------------------------------------
def construir_n_complejo(n,k):
    n_complejo = []
    
    for indice_n, indice_k in zip(n, k):
        n_complejo.append(indice_n + indice_k * 1j)

    return n_complejo


def solo_qsca(event,n,k):#, diametro, pol, inten, AN):
    wavelength, cabra, oveja = leer_archivo(event)
    n_complejo = construir_n_complejo(n,k)
    qsca_values = []

    polarización = 0.
    intensidad = 1e-3
    apertura = 0.2
    diametro = 1000

    for longitud_de_onda, complex_n in zip(wavelength, n_complejo):
        source = Gaussian(
            wavelength=longitud_de_onda * micrometer,
            polarization=polarización * degree,
            optical_power=intensidad * watt,
            NA=apertura * AU
        )
        
        scatterer = Sphere(
            diameter=diametro * nanometer,
            property= complex_n * RIU,
            medium_property=1. * RIU,
            source=source
        )
        
        experiment = Setup(scatterer=scatterer, source=source)
        qsca = experiment.get('Qsca', as_numpy=True)
        qsca_values.append(qsca.tolist())
    return wavelength, qsca_values


def solo_qabs(event, n, k):#, diametro, pol, inten, AN):
    wavelength, cabra, oveja = leer_archivo(event)
    n_complejo = construir_n_complejo(n,k)
    qabs_values = []

    polarización = 0.
    intensidad = 1e-3
    apertura = 0.2
    diametro = 1000

    for longitud_de_onda, complex_n in zip(wavelength, n_complejo):
        source = Gaussian(
            wavelength=longitud_de_onda * micrometer,
            polarization=polarización * degree,
            optical_power=intensidad * watt,
            NA=apertura * AU
        )
        
        scatterer = Sphere(
            diameter=diametro * nanometer,
            property= complex_n * RIU,
            medium_property=1. * RIU,
            source=source
        )
        
        experiment = Setup(scatterer=scatterer, source=source)
        qabs = experiment.get('Qabs', as_numpy=True)
        qabs_values.append(qabs.tolist())
    return wavelength, qabs_values


def solo_qext(event, n, k):#, diametro, pol, inten, AN):
    wavelength, cabra, oveja = leer_archivo(event)
    n_complejo = construir_n_complejo(n,k)
    qext_values = []

    polarización = 0.
    intensidad = 1e-3
    apertura = 0.2
    diametro = 1000

    for longitud_de_onda, complex_n in zip(wavelength, n_complejo):
        source = Gaussian(
            wavelength=longitud_de_onda * micrometer,
            polarization=polarización * degree,
            optical_power=intensidad * watt,
            NA=apertura * AU
        )
        
        scatterer = Sphere(
            diameter=diametro * nanometer,
            property= complex_n * RIU,
            medium_property=1. * RIU,
            source=source
        )
        
        experiment = Setup(scatterer=scatterer, source=source)
        qext = experiment.get('Qext', as_numpy=True)
        qext_values.append(qext.tolist())
    return wavelength, qext_values


def datos_imprimir2(event, n, k):
    long, qsca = solo_qsca(event, n, k)
    cabra, qabs = solo_qabs(event, n, k)
    oveja, qext = solo_qext(event, n, k)

    df = pd.DataFrame({
        'Wavelength': long,
        'Qsca' : qsca,
        'Qabs' : qabs,
        'Qext' : qext,
    }
    ).reset_index(drop=True)

    return df

#-----------------------------------------------
# Representación de factores de eficiencia y .csv descargable
#-----------------------------------------------

def plotear_todo(event, n, k):
    long, qsca_values = solo_qsca(event, n, k)
    long, qabs_values = solo_qabs(event, n, k)
    long, qext_values = solo_qext(event, n, k)

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=long, y=qsca_values,
        mode='lines',
        name='Qsca'
    ))

    fig.add_trace(go.Scatter(
        x=long, y=qabs_values,
        mode='lines',
        name='Qabs'
    ))
    
    fig.add_trace(go.Scatter(
        x=long, y=qext_values,
        mode='lines',
        name='Qext'
    ))

    fig.update_layout(
        title='Efficiency factors',
        xaxis_title='λ (μm)',
        yaxis_title='',
        legend_title=''
    )

    return fig


def generar_csv2():
    try:
        wavelength, n, k = leer_archivo(event)  
        df = datos_imprimir2(event, n, k)     
    except Exception as e:
        df = pd.DataFrame({'Error': [str(e)]})

    sio = StringIO()
    df.to_csv(sio, index=False)
    sio.seek(0)

    bio = BytesIO()
    bio.write(sio.getvalue().encode('utf-8'))
    bio.seek(0)
    return bio


archivo_descargable = pn.widgets.FileDownload(callback=generar_csv2, filename='data2.csv', embed=False)

#--------------------------------------------------------
# Actualización de valores de tablas
#--------------------------------------------------------

def actualizar_tabla(event):
    try:
        df = datos_imprimir(event)
        tabla1.value = df
    except Exception as e:
        tabla1.value = pd.DataFrame({'Error': [str(e)]})

boton_mostrar.on_click(actualizar_tabla)


def actualizar_n_complejo(event):
    wavelength, n, k = leer_archivo(event)
    n_complejo = construir_n_complejo(n, k)
    global vector_n_complejo
    vector_n_complejo = n_complejo

file_input.param.watch(actualizar_n_complejo, 'value')


def actualizar_tabla2(event):
    try:
        wavelength, n, k = leer_archivo(event)  
        df = datos_imprimir2(event, n, k)
        tabla2.value = df     
    except Exception as e:
        tabla2.value = pd.DataFrame({'Error': [str(e)]})

boton_mostrar2.on_click(actualizar_tabla2)

tabla2 = pn.widgets.Tabulator(pd.DataFrame(), name='Efficiency factors', height=400, show_index=False)


def update_graficas_mie(event=None):
    try:
        wavelength, n, k = leer_archivo(event)  
        df = datos_imprimir2(event, n, k)     
        tabla2.value = df
    except Exception as e:
        tabla2.value = pd.DataFrame({'Error': [str(e)]})
    
    grafica_todo.object = plotear_todo(event, n, k)
    
boton_mostrar3.on_click(update_graficas_mie)

#---------------------------------------------
# GUI ANTOLAKETA
#---------------------------------------------
dashboard = pn.Column(
    pn.pane.Markdown("# GUI Variation"),
    pn.pane.Markdown('### Select the .csv file you would like to use to calculate the efficiency factors. The file must store 3 columns with the following names: Wavelength (micrometers), n, k'),
    pn.Row(file_input, boton_mostrar, boton_mostrar2, boton_mostrar3, archivo_descargable),
    pn.Row(tabla1, tabla2, grafica_todo),
)

dashboard.servable()
pn.serve(dashboard, port=8050)
