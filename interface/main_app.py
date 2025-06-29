#ESTA VERSIÓN CREA LA INTERFAZ USANDO PANEL

#OBJETIVO: CORREGIR ERRORES PARA DATOS DE N^2; K=0
#Objetivo: corregir funcionalidad para glass, organic

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
aleaciones = ['Au-Ag', 'Cu-Zn', 'Nb-Sn', 'Ni-Fe', 'V-Ga']
aleaciones_semicond = ['AlAs-GaAs', 'AlSb-GaSb', 'AlN-Al2O3', 'GaAs-InAs', 'GaP-InP', 'GaAs-InAs-GaP-InP', 'Si-Ge']
cristales_mixtos = ['HfO2-Y2O3', 'In2O3-SnO2', 'TlBr-TlCl', 'TlBr-TlI', 'ZrO2-Y2O3']
cristales_dopados = ['Mg-LiTaO3', 'Nb-RbTiOPO4', 'Al-ZnO']
orga = ['(C8H8)n-(C3H3N)m', 'C6H4S4-C12H4N4', 'PTB7-PC71BM', 'P3HT-PC61BM']
orga_inorga = ['Cu-C12H4N4', 'Li-C12H4N4']
plasticos = ['CR-39', 'EVASKY_S87', 'EVASKY_S88', 'NAS-21', 'Optorez1330', 'ZeonexE48R']
cristales_liquidos = ['5CB', '5PCH', 'E7', 'E44', 'MLC-6241-000', 'MLC-6608', 'MLC-9200-000', 'MLC-9200-100', 'TL-216']
resis = ['pmma_resists', 'copolymer_resists', 'negative_tone_photoresists']
cuerpo = ['blood', 'liver', 'colon', 'DNA']
soot = ['acetylene_soot', 'propane_soot', 'diesel_soot']
exot = ['metamaterials']
clay = ['montmorillonite', 'kaolinite', 'illite']
seda = ['Bombyx_mori', 'Antheraea_mylitta', 'Samia_ricini', 'Antheraea_assamensis']

def id_unico(tipo,material,author):
    unico = ()
    if tipo=='main':
        if material=='BiB3O6':
            if author=='Umemura-Î±':
                nombre = str(tipo) + "\\" + str(material) + "\\Umemura-alpha"
            elif author=='Umemura-Î²':
                nombre = str(tipo) + "\\" + str(material) + "\\Umemura-beta"
            elif author=='Umemura-Î³':
                nombre = str(tipo) + "\\" + str(material) + "\\Umemura-gamma"
        elif material=='LiB3O5':
            if author=='Chen-Î±':
                nombre = str(tipo) + "\\" + str(material) + "\\Chen-alpha"
            elif author=='Chen-Î²':
                nombre = str(tipo) + "\\" + str(material) + "\\Chen-beta"
            elif author=='Chen-Î³':
                nombre = str(tipo) + "\\" + str(material) + "\\Chen-gamma"
        elif material=='AgBr':
            nombre = str(tipo) + "\\" + str(material) + "\\Schroter"
        elif material=='Co':
            if author=='Smith':
                nombre = str(tipo) + "\\" + str(material) + "1\\" + str(author)
            else:
                nombre = str(tipo) + "\\" + str(material) + "\\" + str(author)
        elif material=='KNbO3':
            if author=='Umemura-Î±':
                nombre = str(tipo) + "\\" + str(material) + "\\Umemura-alpha"
            elif author=='Umemura-Î²':
                nombre = str(tipo) + "\\" + str(material) + "\\Umemura-beta"
            elif author=='Umemura-Î³':
                nombre = str(tipo) + "\\" + str(material) + "\\Umemura-gamma"
        elif material=='H2O':
            if author=='Asfar-H2O':
                nombre = str(tipo) + "\\" + str(material) + "\\Asfar"
            elif author=='Kedenburg-D2O':
                nombre = str(tipo) + "\\D2O\\Kedenburg"
            elif author=='Asfar-D2O':
                nombre = str(tipo) + "\\D2O\\Asfar"
            else:
                nombre = str(tipo) + "\\" + str(material) + "\\" + str(author)
        #elif material=='MoO3':
        #    if author=='Lajaunie-Î±':
        #        nombre = str(tipo) + "\\" + str(material) + "\\Lajaunie-alpha"
        #    elif author=='Lajaunie-Î²':
        #        nombre = str(tipo) + "\\" + str(material) + "\\Lajaunie-beta"
        #    elif author=='Lajaunie-Î³':
        #        nombre = str(tipo) + "\\" + str(material) + "\\Lajaunie-gamma"
        elif material=='KTiOPO4':
            if author=='Kato-Î±':
                nombre = str(tipo) + "\\" + str(material) + "\\Kato-alpha"
            elif author=='Kato-Î²':
                nombre = str(tipo) + "\\" + str(material) + "\\Kato-beta"
            elif author=='Kato-Î³':
                nombre = str(tipo) + "\\" + str(material) + "\\Kato-gamma"
        elif material=='RbTiOPO4':
            if author=='Carvajal-Î±':
                nombre = str(tipo) + "\\" + str(material) + "\\Carvajal-alpha"
            elif author=='Carvajal-Î²':
                nombre = str(tipo) + "\\" + str(material) + "\\Carvajal-beta"
            elif author=='Carvajal-Î³':
                nombre = str(tipo) + "\\" + str(material) + "\\Carvajal-gamma"
        elif material=='CaSO4':
            if author=='Querry-Î±':
                nombre = str(tipo) + "\\" + str(material) + "\\Querry-alpha"
            elif author=='Querry-Î²':
                nombre = str(tipo) + "\\" + str(material) + "\\Querry-beta"
            elif author=='Querry-Î³':
                nombre = str(tipo) + "\\" + str(material) + "\\Querry-gamma"
        else:
            nombre = str(tipo) + "\\" + str(material) + "\\" + str(author)
    elif tipo=='other':
        autor = author.lower()
        if material=='air':
            nombre = "mixed gases\\" + str(material) + "\\" + str(author)
        elif material in aleaciones:
            nombre = "alloys\\" + str(material) + "\\" + str(author)
        elif material in aleaciones_semicond:
            nombre = "semiconductor alloys\\" + str(material) + "\\" + str(author)
        elif material=="AuAl2":
            nombre = "intermetallics\\" + str(material) + "\\" + str(author)
        elif material in cristales_mixtos:
            nombre = "mixed crystals\\" + str(material) + "\\" + str(author)
        elif material in cristales_dopados:
            if author=="Carvajal-Î±":
                nombre = "doped crystals\\" + str(material) + "\\Carvajal-alpha"
            elif author=="Carvajal-Î²":
                nombre = "doped crystals\\" + str(material) + "\\Carvajal-beta"
            elif author=='Carvajal-Î³':
                nombre = "doped crystals\\" + str(material) + "\\Carvajal-gamma"
            else:
                nombre = "doped crystals\\" + str(material) + "\\" + str(author)
        elif material in orga:
            nombre = "mixed organic\\" + str(material) + "\\" + str(author)
        elif material in orga_inorga:
            nombre = "mixed organic-inorganic\\" + str(material) + "\\" + str(author)
        elif material in plasticos:
            nombre = "commercial plastics\\" + str(material) + "\\" + str(author)
        elif material in cristales_liquidos:
            nombre = "liquid crystals\\" + str(material) + "\\" + str(author)
        elif material=='SiN':
            nombre = "anti-reflective coatings\\" + str(material) + "\\" + str(author)
        elif material in cuerpo:
            nombre = "human body\\" + str(material) + "\\" + str(author)
        elif material in exot:
            nombre = "exotic\\" + str(material) + "\\" + str(author)
        elif material=='BK7_matching_liquid':
            nombre = "index-matching liquids\\" + str(autor) + "\\" + str(material)
        elif material=='acrylic_matching_liquid':
            nombre = "index-matching liquids\\" + str(autor) + "\\" + str(material)
        elif material in resis:
            alpha = "other\\resists\\"
            if material=='pmma_resists':
                if author=='Microchem495':
                    nombre = alpha + "Microchem 495"
                elif author=='Microchem950':
                    nombre = alpha + "Microchem 950"
            elif material=='copolymer_resists':
                nombre = alpha + "Microchem 8.5 mEL"
            elif material=='negative_tone_photoresists':
                if author=='MicroresistEpoCore':
                    nombre = alpha + "Micro resist EpoCore"
                elif author=='MicroresistEpoClad':
                    nombre = alpha + "Micro resist EpoClad"
                elif author=='Microchem_SU8_2000':
                    nombre = alpha + "Microchem SU-8 2000"
                elif author=='Microchem_SU8_3000':
                    nombre = alpha + "Microchem SU-8 3000"
        elif material=='CH3NH3PbI3':
            nombre = "perovskite\\" + str(material) + "\\" + str(author)
        elif material in soot:
            alpha = "soots\\"
            if material=='acetylene_soot':
                nombre = alpha + "acetylene soot\\Dalzell"
            elif material=='propane_soot':
                nombre = alpha + "propane soot\\Dalzell"
            elif material=='diesel_soot':
                nombre = alpha + "diesel soot\\" + str(author)
        elif material in clay:
            alpha = "clays\\"
            nombre = alpha + str(material) + "\\" + str(author)
        elif material in seda:
            alpha = "silk\\"
            if material=='Bombyx_mori':
                nombre = alpha + "Bombyx mori\\" + str(author)
            elif material=='Antheraea_mylitta':
                nombre = alpha + "Antheraea mylitta\\" + str(author)
            elif material=='Samia_ricini':
                nombre = alpha + "Samia ricini\\" + str(author)
            elif material=='Antheraea_assamensis':
                nombre = alpha + "Antheraea assamensis\\" + str(author)
    elif tipo=='glass':
        autor = author.lower()
        if author=='SCHOTT':
            if material=='FK51':
                nombre = str(tipo) + "\\" + str(autor) + "\\N-" + str(material) + "A"
            else:
                nombre = str(tipo) + "\\" + str(autor) + "\\N-" + str(material) 
        elif author=='OHARA':
            if material=='BK7':
                nombre = str(tipo) + "\\" + str(autor) + "\\S-BSL7"
            elif material=='BAF10':
                nombre = str(tipo) + "\\" + str(autor) + "\\S-BAH10"
            elif material=='BAK1':
                nombre = str(tipo) + "\\" + str(autor) + "\\S-BAL11"
        elif author=='HIKARI':
            if material=='BK7':
                nombre = str(tipo) + "\\" + str(autor) + "\\J-BK7A"
            else:
                nombre = str(tipo) + "\\" + str(autor) + "\\J-" + str(material)
        elif author=='CDGM':
            if material=='BK7':
                nombre = str(tipo) + "\\" + str(autor) + "\\H-K9L"
            elif material=='BAF10':
                nombre = str(tipo) + "\\" + str(autor) + "\\H-ZBAF52"
            elif material=='BAK1':
                nombre = str(tipo) + "\\" + str(autor) + "\\H-BAK8"
            elif material=='SF5':
                nombre = str(tipo) + "\\" + str(autor) + "\\ZF2"
            elif material=='SF10':
                nombre = str(tipo) + "\\" + str(autor) + "\\ZF4"
            elif material=='SF11':
                nombre = str(tipo) + "\\" + str(autor) + "\\ZF13"
        elif author=="HOYA":
            if material=='BK7':
                nombre = str(tipo) + "\\" + str(autor) + "\\BSC7"
            elif material=='BAF10':
                nombre = str(tipo) + "\\" + str(autor) + "\\BAF10"
            elif material=='SF5':
                nombre = str(tipo) + "\\" + str(autor) + "\\E-FD5"
            elif material=='SF10':
                nombre = str(tipo) + "\\" + str(autor) + "\\E-FD10"
        elif author=="SUMITA":
            if material=='BK7':
                nombre = str(tipo) + "\\" + str(autor) + "\\K-BK7"
            elif material=='SF11':
                nombre = str(tipo) + "\\" + str(autor) + "\\K-SFLD11"
        elif author=="LZOS":
            nombre = str(tipo) + "\\" + str(autor) + "\\K8"
    elif tipo=='organic':
        if material=='methane':
            nombre = str(tipo) + "\\CH4 - " + str(material) + "\\" + str(author)
        elif material=='ethane':
            nombre = str(tipo) + "\\C2H6 - " + str(material) + "\\" + str(author)
        elif material=='pentane':
            nombre = str(tipo) + "\\C5H12 - " + str(material) + "\\" + str(author)
        elif material=='hexane':
            nombre = str(tipo) + "\\C6H14 - " + str(material) + "\\" + str(author)
        elif material=='heptane':
            nombre = str(tipo) + "\\C7H16 - " + str(material) + "\\" + str(author)
        elif material=='octane':
            nombre = str(tipo) + "\\C8H18 - " + str(material) + "\\" + str(author)
        elif material=='acetylene':
            nombre = str(tipo) + "\\C2H2 - " + str(material) + "\\" + str(author)
        elif material=='ethylene':
            nombre = str(tipo) + "\\C2H4 - " + str(material) + "\\" + str(author)
        elif material=='methanol':
            nombre = str(tipo) + "\\CH4O - " + str(material) + "\\" + str(author)
        elif material=='ethanol':
            nombre = str(tipo) + "\\C2H6O - " + str(material) + "\\" + str(author)
        elif material=='propanol':
            nombre = str(tipo) + "\\C3H8O - " + str(material) + "\\" + str(author)
        elif material=='butanol':
            nombre = str(tipo) + "\\C4H10O - " + str(material) + "\\" + str(author)
        elif material=='pentanol':
            nombre = str(tipo) + "\\C5H12O - " + str(material) + "\\" + str(author)
        elif material=='octanol':
            nombre = str(tipo) + "\\C8H18O - " + str(material) + "\\" + str(author)
        elif material=='glycerol':
            nombre = str(tipo) + "\\C3H8O3 - " + str(material) + "\\" + str(author)
        elif material=='ethylene_glycol':
            nombre = str(tipo) + "\\C2H6O2 - ethylene glycol\\" + str(author)
        elif material=='propylene_glycol':
            nombre = str(tipo) + "\\C3H8O2 - propylene glycol\\" + str(author)
        elif material=='ethyl_acetate':
            nombre = str(tipo) + "\\C4H8O2 - ethyl acetate\\" + str(author)
        elif material=='methyl_salicylate':
            nombre = str(tipo) + "\\C8H8O3 - methyl salicylate\\" + str(author)
        elif material=='ethyl_salicylate':
            nombre = str(tipo) + "\\C9H10O3 - ethyl salicylate\\" + str(author)
        elif material=='ethyl_cinnamate':
            nombre = str(tipo) + "\\C11H12O2 - ethyl cinnamate\\" + str(author)
        elif material=='diethyl_phthalate':
            nombre = str(tipo) + "\\C12H14O4 - diethyl phthalate\\" + str(author)
        elif material=='cyclohexane':
            nombre = str(tipo) + "\\C6H12 - cyclohexane\\" + str(author)
        elif material=='benzene':
            nombre = str(tipo) + "\\C6H6 - benzene\\" + str(author)
        elif material=='styrene':
            nombre = str(tipo) + "\\C8H8 - styrene\\" + str(author)
        elif material=='toluene':
            nombre = str(tipo) + "\\C7H8 - toluene\\" + str(author)
        elif material=='nitrobenzene':
            nombre = str(tipo) + "\\C6H5NO2 - nitrobenzene\\" + str(author)
        elif material=='dioxane':
            nombre = str(tipo) + "\\C4H8O2 - dioxane\\" + str(author)
        elif material=='acetic_acid':
            nombre = str(tipo) + "\\C2H4O2 - acetic acid\\" + str(author)
        elif material=='pentanediol':
            nombre = str(tipo) + "\\C5H12O2 - pentanediol\\" + str(author)
        elif material=='acetone':
            nombre = str(tipo) + "\\C3H6O - acetone\\" + str(author)
        elif material=='bromoform':
            nombre = str(tipo) + "\\CHBr3 - bromoform\\" + str(author)
        elif material=='dichloromethane':
            nombre = str(tipo) + "\\CH2Cl2 - dichloromethane\\" + str(author)
        elif material=='chloroform':
            nombre = str(tipo) + "\\CHCl3 - chloroform\\" + str(author)
        elif material=='carbon_tetrachloride':
            nombre = str(tipo) + "\\CCl4 - carbon tetrachloride\\" + str(author)
        elif material=='acetonitrile':
            nombre = str(tipo) + "\\C2H3N - acetonitrile\\" + str(author)
        elif material=='urea':
            nombre = str(tipo) + "\\CH4N2O - urea\\" + str(author)
        elif material=='dimethyl_sulfoxide':
            nombre = str(tipo) + "\\C2H6OS - dimethyl sulfoxide\\" + str(author)
        elif material=='dimethyl_methylphosphonate':
            nombre = str(tipo) + "\\C3H9O3P - dimethyl methylphosphonate\\" + str(author)
        elif material=='diethyl_sulfite':
            nombre = str(tipo) + "\\C4H10O3S - diethyl sulfite\\" + str(author)
        elif material=='diisopropyl_methylphosphonate':
            nombre = str(tipo) + "\\C7H17O3P -  diisopropyl methylphosphonate\\" + str(author)
        elif material=='potassium_hydrogen_phthalate':
            nombre = str(tipo) + "\\C8H5KO4 - potassium hydrogen phthalate\\" + str(author)
        elif material=='cinnamaldehyde':
            nombre = str(tipo) + "\\C9H8O - cinnamaldehyde\\" + str(author)
        elif material=='polyvinyl_alcohol':
            nombre = str(tipo) + "\\(C2H4O)n - polyvinyl alcohol\\" + str(author)
        elif material=='polydimethylsiloxane':
            nombre = str(tipo) + "\\(C2H6OSi)n - polydimethylsiloxane\\" + str(author)
        elif material=='polylactic_acid':
            nombre = str(tipo) + "\\(C3H4O2)n - polylactic acid\\" + str(author)
        elif material=='poly(methyl_methacrylate)':
            nombre = str(tipo) + "\\(C5H8O2)n - poly(methyl methacrylate)\\" + str(author)
        elif material=='polyvinylpyrrolidone':
            nombre = str(tipo) + "\\(C6H9NO)n - polyvinylpyrrolidone\\" + str(author)
        elif material=='cellulose':
            nombre = str(tipo) + "\\(C6H10O5)n - cellulose\\" + str(author)
        elif material=='poly(N-isopropylacrylamide)':
            nombre = str(tipo) + "\\(C6H11NO)n - poly(N-isopropylacrylamide)\\" + str(author)
        elif material=='polystyren':
            nombre = str(tipo) + "\\(C8H8)n - polystyren\\" + str(author)
        elif material=='polycarbonate':
            nombre = str(tipo) + "\\(C16H14O3)n - polycarbonate\\" + str(author)
    list_of_pages = db_operations.search_pages(nombre)

    # Comprobar que hay resultados
    if not list_of_pages:
        return None

    # Suponiendo que cada fila es una lista o tupla, y que el primer elemento es el ID
    ids = [fila[0] for fila in list_of_pages]
    unico = ids[0]
    print(f"[DEBUG] Nombre buscado: {nombre}")
    print(f"[DEBUG] Resultados: {list_of_pages}")
    return unico


def info_completa_n(tipo, material, author):
    id = id_unico(tipo, material, author)
    print(f"[DEBUG] pageid: {id}, tipo: {type(id)}")
    info_n = db_operations.get_material_n_numpy(id)
    if info_n is None or len(info_n) == 0:
        return [], []
    long = [tupla[0] for tupla in info_n]
    index = [tupla[1] for tupla in info_n]

    return long, index


def info_completa_k(tipo, material, author):
    id = id_unico(tipo, material, author)
    print(f"[DEBUG] pageid: {id}, tipo: {type(id)}")
    info_k = db_operations.get_material_k_numpy(id)
    if info_k is None or len(info_k) == 0:
        info_k_var = db_operations.get_material_n_numpy(id)
        long = [tupla[0] for tupla in info_k_var]
        kappa = [0] * len(long)
    else:
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
        mode='lines+markers',
        name='Qsca',
        marker=dict(symbol='circle', size=8)
    ))

    fig.add_trace(go.Scatter(
        x=long, y=qabs_values,
        mode='lines+markers',
        name='Qabs',
        marker=dict(symbol='triangle-up', size=8)
    ))

    fig.add_trace(go.Scatter(
        x=long, y=qext_values,
        mode='lines+markers',
        name='Qext',
        marker=dict(symbol='square', size=8)
    ))

    fig.update_layout(
        title='Efficiency factors',
        xaxis_title='λ (μm)',
        yaxis_title='',
        legend_title='',
        template='plotly_white'  # Fondo blanco para mejor contraste
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
shelves = db_operations.get_all_shelves()
shelves_filtered = [shelf for shelf in shelves if shelf != '3d']

# SELECCIÓN DE LOS DATOS POR EL USUARIO + BOTÓN PARA SIGUIENTE FASE:
material1_dropdown = pn.widgets.Select(name='Material type', options=shelves_filtered, width=340)
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

    if tipo=='organic':
        valores_a_descartar = ['ethylene_glycol', 'propylene_glycol', 'diisopropyl_methylphosphonate', 'potassium_hydrogen_phthalate'] 
        filtered_list = [material for material in lista_materiales if material not in valores_a_descartar]
    elif tipo=='glass':
        valores_a_incluir = ['BK7', 'BAF10', 'BAK1', 'FK51', 'LASF9', 'SF5', 'SF10', 'SF11']
        filtered_list = [material for material in lista_materiales if material in valores_a_incluir]
    elif tipo=='other':
        valores_a_descartar = ['fused_silica_matching_liquid', 'TherminolVP-1', 'Optical adhesives', 'Leica_Type_F', 'Olympus_IMMOIL-F30CC', 'Sigma_Aldrich_M5904']
        filtered_list = [material for material in lista_materiales if material not in valores_a_descartar]
    elif tipo=='main':
        valores_a_descartar = ['B', 'B4C', 'Ce', 'Er', 'Eu', 'GaSe', 'Ho', 'Lu', 'Mg', 'SiO', 'Pr', 'Sr', 'Tm', 'Yb', 'HfO2', 'ZrO2', 'GaSb']
        filtered_list = [material for material in lista_materiales if material not in valores_a_descartar]
    else:
        filtered_list = lista_materiales

    material2_dropdown.options = filtered_list

material1_dropdown.param.watch(update_material_dropdown, 'value')
update_material_dropdown()


def update_author_dropdown(event=None):
    lista_autores = []
    tipo = material1_dropdown.value
    mat = material2_dropdown.value
    lista_autores = db_operations.get_author(tipo, mat)

    if tipo=='main':
        if mat=='H2O':
            valores_a_descartar = ['Wang', 'Wang-D2O', 'Segelstein']
            filtered_list = [autor for autor in lista_autores if autor not in valores_a_descartar]
        elif mat=='Ag':
            valores_a_descartar = ['Hagemann']
            filtered_list = [autor for autor in lista_autores if autor not in valores_a_descartar]
        elif mat=='Al':
            valores_a_descartar = ['Hagemann']
            filtered_list = [autor for autor in lista_autores if autor not in valores_a_descartar]
        elif mat=='Au':
            valores_a_descartar = ['Hagemann', 'Hagemann']
            filtered_list = [autor for autor in lista_autores if autor not in valores_a_descartar]
        elif mat=='Bi':
            valores_a_descartar = ['Hagemann']
            filtered_list = [autor for autor in lista_autores if autor not in valores_a_descartar]
        elif mat=='C':
            valores_a_descartar = ['Hagemann']
            filtered_list = [autor for autor in lista_autores if autor not in valores_a_descartar]
        elif mat=='Ca':
            valores_a_descartar = ['Rodriguez-de_Marcos']
            filtered_list = [autor for autor in lista_autores if autor not in valores_a_descartar]
        elif mat=='Cu':
            valores_a_descartar = ['Hagemann']
            filtered_list = [autor for autor in lista_autores if autor not in valores_a_descartar]
        elif mat=='Al2O3':
            valores_a_descartar = ['Hagemann']
            filtered_list = [autor for autor in lista_autores if autor not in valores_a_descartar]
        elif mat=='MoO3':
            valores_a_descartar = ['Lajaunie-Î±', 'Lajaunie-Î²', 'Lajaunie-Î³']
            filtered_list = [autor for autor in lista_autores if autor not in valores_a_descartar]
        elif mat=='SiO':
            valores_a_descartar = ['Hass']
            filtered_list = [autor for autor in lista_autores if autor not in valores_a_descartar]
        elif mat=='MoS2':
            valores_a_descartar = ['Yim-2nm', 'Yim-3nm', 'Yim-20nm']
            filtered_list = [autor for autor in lista_autores if autor not in valores_a_descartar]
        elif mat=='Si':
            valores_a_descartar = ['Green-1995', 'Daub']
            filtered_list = [autor for autor in lista_autores if autor not in valores_a_descartar]
        else:
            filtered_list = lista_autores

    material3_dropdown.options = filtered_list

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



