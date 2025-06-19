#ESTA VERSIÓN CREA LA INTERFAZ USANDO PANEL

#Objetivo: conectar la base de RefractiveIndex con la interfaz web

#Actualizaciones respecto versión anterior: pretende incluir diámetro, intensidad, polarización
#y AN como parámetros de entrada de las funciones solo_qsca y similar

import panel as pn
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
pn.extension('plotly')
pn.extension('tabulator')
import matplotlib.pyplot as plt
from io import StringIO, BytesIO

import sys
from pathlib import Path

base_path = Path(__file__).resolve().parents[1]  # Esto apunta a PyMieSim_local
sys.path.insert(0, str(base_path))

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.units import nanometer, meter, degree, watt, AU, RIU, micrometer

from refractiveindexmaster.refractiveindex2.refractivesqlite.dboperations import Database

db_operations = Database('C:/Users/usuario/Desktop/TFG/venv/PyMieSim_local/refractiveindexmaster/refractiveindex2/refractive.db')

# Funciones definidas en "db_operations":
#shelves = db_operations.get_all_shelves() 
#books = db_operations.get_books_in_shelf('shelf') 
#pages = db_operations.get_pages('shelf','book')
#todo = db_operations.search_pages()


#---------------------------------------------------------------------------------------------
# FUNCIONES OBTENCIÓN DE LOS DATOS DE LA BASE DE SQLITE (REFRACTIVE INDEX INFO)
#---------------------------------------------------------------------------------------------

def id_unico(tipo,material,author):
    nombre = str(tipo) + "\\" + str(material) + "\\" + str(author)
    list_of_pages = db_operations.search_pages(nombre)

    # Comprobar que hay resultados
    if not list_of_pages:
        return []

    # Suponiendo que cada fila es una lista o tupla, y que el primer elemento es el ID
    ids = [fila[0] for fila in list_of_pages]
    unico = ids[0]

    return unico


def info_completa_n(tipo, material, author):
    id = id_unico(tipo, material, author)
    info_n = db_operations.get_material_n_numpy(id)
    if info_n is None or len(info_n) == 0:
        return [], []
    long = [tupla[0] for tupla in info_n]
    index = [tupla[1] for tupla in info_n]

    return long, index


def info_completa_k(tipo, material, author):
    id = id_unico(tipo, material, author)
    info_k = db_operations.get_material_k_numpy(id)
    if info_k is None or len(info_k) == 0:
        return [], []
    long = [tupla[0] for tupla in info_k]
    kappa = [tupla[1] for tupla in info_k]

    return long, kappa



#---------------------------------------------------------------------------------------------
# OTRAS FUNCIONES - REPRESENTACIÓN DE DATOS DE N Y KAPPA
#---------------------------------------------------------------------------------------------

def datos_imprimir(tipo, material, author):
    long, index = info_completa_n(tipo, material, author)
    _, kappa = info_completa_k(tipo, material, author)
    df = pd.DataFrame({
        'wavelength' : long,
        'n' : index,
        'k' : kappa,
    }
    ).reset_index(drop=True)

    return df


def plotear_n_k(tipo,material,author):
    df = datos_imprimir(tipo,material,author)
    df.columns = df.columns.str.strip()

    long = df['wavelength']
    index = df['n']
    kappa = df['k']

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=long, y=index,
        mode='lines',
        name='Real refractive index (n)'
    ))

    fig.add_trace(go.Scatter(
        x=long, y=kappa,
        mode='lines',
        name='Extinction coefficient (k)'
    ))

    fig.update_layout(
        title='Optical properties',
        xaxis_title='λ (μm)',
        yaxis_title='',
        legend_title=''
    )

    return fig



#---------------------------------------------------------------------------------------------
# FUNCIONES PARA CÁLCULOS DE PROPIEDADES A PARTIR DE ÍNDICES DE MIE (PYMIESIM)
#---------------------------------------------------------------------------------------------

def construir_n_complejo(tipo, material, author):
    n_complejo = []
    df = datos_imprimir(tipo,material,author)
    df.columns = df.columns.str.strip()

    index = df['n']
    kappa = df['k']
    
    for indice_n, indice_k in zip(index, kappa):
        # Crear un número complejo y agregarlo a la lista
        n_complejo.append(indice_n + indice_k * 1j)

    return n_complejo


def solo_qsca(tipo, material, author, diametro, pol, inten, AN):
    df = datos_imprimir(tipo, material, author)
    long = df['wavelength'].tolist()
    #index = df['n'].tolist()
    n_complejo = construir_n_complejo(tipo, material, author)
    qsca_values = []

    # Parámetros de la fuente:
    if pol=='' or pol is None:
        polarización = [0.]
    else:
        polarización = pol
    if inten=='' or inten is None:
        intensidad = [1e-3]
    else:
        intensidad = inten
    if AN=='' or AN is None:
        apertura = [0.2]
    else:
        apertura = AN
    #Parámetros de la partícula:
    if diametro=='' or diametro is None:
        diametro = [1000]

    for longitud_de_onda, complex_n in zip(long, n_complejo):
        #print(f"Processing wavelength: {longitud_de_onda}, n_complejo: {complex_n}")  # Debugging line
        source = Gaussian(
            wavelength=longitud_de_onda * micrometer,
            polarization=float(polarización) * degree,
            optical_power=float(intensidad) * watt,
            NA=float(apertura) * AU
        )
        
        scatterer = Sphere(
            diameter=diametro * nanometer,
            property= complex_n * RIU,
            medium_property=1. * RIU,
            source=source
        )
        
        experiment = Setup(scatterer=scatterer, source=source)
        qsca = experiment.get('Qsca', as_numpy=True)
        #print(f"Qsca: {qsca}")  # Verifica si Qsca se está calculando correctamente
        qsca_values.append(qsca.tolist())
    return long, qsca_values


def solo_qabs(tipo, material, author, diametro, pol, inten, AN):
    df = datos_imprimir(tipo, material, author)
    long = df['wavelength'].tolist()
    #index = df['n'].tolist()
    n_complejo = construir_n_complejo(tipo, material, author)
    qabs_values = []

    # Parámetros de la fuente:
    if pol=='' or pol is None:
        polarización = [0.]
    else:
        polarización = pol
    if inten=='' or inten is None:
        intensidad = [1e-3]
    else:
        intensidad = inten
    if AN=='' or AN is None:
        apertura = [0.2]
    else:
        apertura = AN
    #Parámetros de la partícula:
    if diametro=='' or diametro is None:
        diametro = [1000]

    for longitud_de_onda, complex_n in zip(long, n_complejo):
        #print(f"Processing wavelength: {longitud_de_onda}, n_complejo: {complex_n}")  # Debugging line
        source = Gaussian(
            wavelength=longitud_de_onda * micrometer,
            polarization=float(polarización) * degree,
            optical_power=float(intensidad) * watt,
            NA=float(apertura) * AU
        )
        
        scatterer = Sphere(
            diameter=diametro * nanometer,
            property= complex_n * RIU,
            medium_property=1. * RIU,
            source=source
        )
        
        experiment = Setup(scatterer=scatterer, source=source)
        qabs = experiment.get('Qabs', as_numpy=True)
        #print(f"Qabs: {qabs}")  # Verifica si Qabs se está calculando correctamente
        qabs_values.append(qabs.tolist())
    return long, qabs_values


def solo_qext(tipo, material, author, diametro, pol, inten, AN):
    df = datos_imprimir(tipo, material, author)
    long = df['wavelength'].tolist()
    #index = df['n'].tolist()
    n_complejo = construir_n_complejo(tipo, material, author)
    qext_values = []

    # Parámetros de la fuente:
    if pol=='' or pol is None:
        polarización = [0.]
    else:
        polarización = pol
    if inten=='' or inten is None:
        intensidad = [1e-3]
    else:
        intensidad = inten
    if AN=='' or AN is None:
        apertura = [0.2]
    else:
        apertura = AN
    #Parámetros de la partícula:
    if diametro=='' or diametro is None:
        diametro = [1000]

    for longitud_de_onda, complex_n in zip(long, n_complejo):
        #print(f"Processing wavelength: {longitud_de_onda}, n_complejo: {complex_n}")  # Debugging line
        source = Gaussian(
            wavelength=longitud_de_onda * micrometer,
            polarization=float(polarización) * degree,
            optical_power=float(intensidad) * watt,
            NA=float(apertura) * AU
        )
        
        scatterer = Sphere(
            diameter=1000. * nanometer,
            property= complex_n * RIU,
            medium_property=1. * RIU,
            source=source
        )
        
        experiment = Setup(scatterer=scatterer, source=source)
        qext = experiment.get('Qext', as_numpy=True)
        #print(f"Qext: {qext}")  # Verifica si Qext se está calculando correctamente
        qext_values.append(qext.tolist())
    return long, qext_values



#---------------------------------------------------------------------------------------------
# FUNCIONES PARA REPRESENTACIÓN DE PROPIEDADES DERIVADAS DE ÍNDICES DE MIE
#---------------------------------------------------------------------------------------------

def plotear_todo(tipo,material,author, diametro, pol, inten, AN):
    long, qsca_values = solo_qsca(tipo, material, author, diametro, pol, inten, AN)
    long, qabs_values = solo_qabs(tipo, material, author, diametro, pol, inten, AN)
    long, qext_values = solo_qext(tipo, material, author, diametro, pol, inten, AN)

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



#---------------------------------------------------------------------------------------------
# FUNCIONES PARA CREAR TABLA DE VALORES E IMPRIMIR DATOS EN FORMATO .CSV
#---------------------------------------------------------------------------------------------

def datos_imprimir2(tipo,material,author, diametro, pol, inten, AN):
    long, qsca = solo_qsca(tipo,material,author, diametro, pol, inten, AN)
    cabra, qabs = solo_qabs(tipo,material,author, diametro, pol, inten, AN)
    oveja, qext = solo_qext(tipo,material,author, diametro, pol, inten, AN)

    df = pd.DataFrame({
        'Wavelength': long,
        'Qsca' : qsca,
        'Qabs' : qabs,
        'Qext' : qext,
    }
    ).reset_index(drop=True)

    return df



#---------------------------------------------------------------------------------------------
# WIDGETS 1 - SELECCIÓN DE UN TIPO, MATERIAL Y AUTOR CONCRETO
#---------------------------------------------------------------------------------------------

# SELECCIÓN DE LOS DATOS POR EL USUARIO + BOTÓN PARA SIGUIENTE FASE:
material1_dropdown = pn.widgets.Select(name='Material type', options=db_operations.get_all_shelves(), width=340)
material2_dropdown = pn.widgets.Select(name='Material', options=[], width=340)
material3_dropdown = pn.widgets.Select(name='Reference', options=[], width=340)
boton_mostrar = pn.widgets.Button(name='Show table', button_type='primary')


#---------------------------------------------------------------------------------------------
# FUNCIONES UPDATE - SELECCIÓN DE UN TIPO, MATERIAL Y AUTOR CONCRETO
#---------------------------------------------------------------------------------------------

def update_material_dropdown(event=None):
    lista_materiales = []
    tipo = material1_dropdown.value
    lista_materiales = db_operations.get_books_in_shelf(tipo)
    material2_dropdown.options = lista_materiales

material1_dropdown.param.watch(update_material_dropdown, 'value')
update_material_dropdown()


def update_author_dropdown(event=None):
    lista_autores = []
    tipo = material1_dropdown.value
    mat = material2_dropdown.value
    lista_autores = db_operations.get_author(tipo, mat)
    material3_dropdown.options = lista_autores

material1_dropdown.param.watch(update_author_dropdown, 'value')
material2_dropdown.param.watch(update_author_dropdown, 'value')
update_author_dropdown()



#---------------------------------------------------------------------------------------------
# FUNCIONES UPDATE - MOSTRAR DATAFRAME, GRÁFICA n-k Y ARCHIVO .CSV DESCARGABLE
#---------------------------------------------------------------------------------------------

def actualizar_tabla(event):
    tipo = material1_dropdown.value
    material = material2_dropdown.value
    autor = material3_dropdown.value

    try:
        df = datos_imprimir(tipo, material, autor)
        datos_material.value = df
    except Exception as e:
        datos_material.value = pd.DataFrame({'Error': [str(e)]})

boton_mostrar.on_click(actualizar_tabla)


def update_graficas(event=None):
    tipo = material1_dropdown.value
    material = material2_dropdown.value
    author = material3_dropdown.value
    if author == None:
        return
    grafica1.object = plotear_n_k(tipo,material,author)

boton_mostrar.on_click(update_graficas)


def generar_csv():
    tipo = material1_dropdown.value
    material = material2_dropdown.value
    author = material3_dropdown.value
    if author is None:
        return None

    df = datos_imprimir(tipo, material, author)
    sio = StringIO()
    df.to_csv(sio, index=False)
    sio.seek(0)

    # Convertimos el contenido a BytesIO para el widget:
    # El atributo callback de pn.widgets.FileDownload debe ser una función que devuelva un objeto 
    # de tipo BytesIO o similar, no el contenido directamente.
    bio = BytesIO()
    bio.write(sio.getvalue().encode('utf-8'))
    bio.seek(0)
    return bio



#---------------------------------------------------------------------------------------------
# WIDGETS 2 - MOSTRAR DATAFRAME, GRÁFICA n-k Y ARCHIVO .CSV DESCARGABLE
#---------------------------------------------------------------------------------------------

grafica1 = pn.pane.Plotly()
datos_material = pn.widgets.Tabulator(pd.DataFrame(), name='n and kappa', height=400, show_index=False)
archivo_descargable = pn.widgets.FileDownload(callback=generar_csv, filename='data.csv', embed=False)
boton_mostrar2 = pn.widgets.Button(name='Calculate efficiency factors', button_type='primary')
# OJO!!!!! callback espera una función que pueda llamar luego, no el resultado de esa función!
# NO PARÉNTESIS



#---------------------------------------------------------------------------------------------
# WIDGETS 3 - PARÁMETROS DE ENTRADA EXTRA
#---------------------------------------------------------------------------------------------

d_input = pn.widgets.FloatInput(name='Diameter in nm', value=1000, width=260)
pol_input = pn.widgets.FloatInput(name='Polarization in sexagesimal degrees', value=0, width=260)
int_input = pn.widgets.FloatInput(name='Incident intensity in watt', value=1e-3, width=260)
AN_input = pn.widgets.FloatInput(name='Numerical Aperture', value=0.2, width=260)

#---------------------------------------------------------------------------------------------
# REPRESENTAR QSCA, QABS, QEXT Y DF + .CSV DESCARGABLE
#---------------------------------------------------------------------------------------------

def update_graficas_mie(event=None):
    tipo = material1_dropdown.value
    material = material2_dropdown.value
    author = material3_dropdown.value
    diametro = d_input.value
    pol = pol_input.value
    inten = int_input.value
    AN = AN_input.value

    if author == None:
        return
    grafica_todo.object = plotear_todo(tipo,material,author, diametro, pol, inten, AN)
    
boton_mostrar2.on_click(update_graficas_mie)


def actualizar_tabla2(event):
    tipo = material1_dropdown.value
    material = material2_dropdown.value
    autor = material3_dropdown.value
    diametro = d_input.value
    pol = pol_input.value
    inten = int_input.value
    AN = AN_input.value

    try:
        df = datos_imprimir2(tipo, material, autor, diametro, pol, inten, AN)
        datos_material2.value = df
    except Exception as e:
        datos_material2.value = pd.DataFrame({'Error': [str(e)]})

boton_mostrar2.on_click(actualizar_tabla2)


def generar_csv2():
    tipo = material1_dropdown.value
    material = material2_dropdown.value
    author = material3_dropdown.value
    diametro = d_input.value
    pol = pol_input.value
    inten = int_input.value
    AN = AN_input.value

    if author is None:
        return None

    df = datos_imprimir2(tipo, material, author, diametro, pol, inten, AN)
    sio = StringIO()
    df.to_csv(sio, index=False)
    sio.seek(0)

    # Convertimos el contenido a BytesIO para el widget:
    # El atributo callback de pn.widgets.FileDownload debe ser una función que devuelva un objeto 
    # de tipo BytesIO o similar, no el contenido directamente.
    bio = BytesIO()
    bio.write(sio.getvalue().encode('utf-8'))
    bio.seek(0)
    return bio



#---------------------------------------------------------------------------------------------
# WIDGETS 3 - REPRESENTAR QSCA
#---------------------------------------------------------------------------------------------

grafica_todo = pn.pane.Plotly(width=600)
datos_material2 = pn.widgets.Tabulator(pd.DataFrame(), name='Efficiency factors', height=400, show_index=False)
archivo_descargable2 = pn.widgets.FileDownload(callback=generar_csv2, filename='data2.csv', embed=False)


#---------------------------------------------------------------------------------------------
# LAYOUT - ORGANIZACIÓN DE LOS WIDGETS EN LA GUI
#---------------------------------------------------------------------------------------------

dashboard = pn.Column(
    pn.pane.Markdown("# Graphical User Interface"),
    pn.pane.Markdown("## Refractive index data extraction and representation - From Refractive Index Info DataBase"),
    pn.Row(material1_dropdown, material2_dropdown, material3_dropdown),
    boton_mostrar,
    pn.Row(grafica1, datos_material, archivo_descargable),
    pn.pane.Markdown("## Efficiency factor calculation and representation - Using PyMieSim calculation software"),
    pn.Row(d_input, pol_input, int_input, AN_input),
    boton_mostrar2,
    pn.Row(grafica_todo, datos_material2, archivo_descargable2),
)


dashboard.servable()
pn.serve(dashboard, port=8050)



