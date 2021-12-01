
# Needed packages

from PyBioMed import Pyprotein
import pandas as pd
from pandas import DataFrame
from Bio import SeqIO
import numpy as np

import streamlit as st
from pathlib import Path
import base64
from PIL import Image
import io

# Streamlit config

#%%

#---------------------------------#
# Page layout
## Page expands to full width
st.set_page_config(page_title='LIDEB Tools - Druggability',
    layout='wide')

######
# Function to put a picture as header   
def img_to_bytes(img_path):
    img_bytes = Path(img_path).read_bytes()
    encoded = base64.b64encode(img_bytes).decode()
    return encoded


image = Image.open('cropped-header.png')
st.image(image)

st.markdown("![Twitter Follow](https://img.shields.io/twitter/follow/LIDeB_UNLP?style=social)")
st.subheader(":pushpin:" "About Us")
st.markdown("We are a drug discovery team with an interest in the development of publicly available open-source customizable cheminformatics tools to be used in computer-assisted drug discovery. We belong to the Laboratory of Bioactive Research and Development (LIDeB) of the National University of La Plata (UNLP), Argentina. Our research group is focused on computer-guided drug repurposing and rational discovery of new drug candidates to treat epilepsy and neglected tropical diseases.")
st.markdown(":computer:""**Web Site** " "<https://lideb.biol.unlp.edu.ar>")


# Introduction
#---------------------------------#

st.write("""
# LIDeB Tools

**It is a free web-application for Druggability Prediction**

The tool uses the following packages [RDKIT](https://www.rdkit.org/docs/index.html), [PyBioMed](https://github.com/gadsbyfly/PyBioMed/), [Biopython](https://biopython.org/), [Plotly](https://plotly.com/python/)

The next workflow summarizes the steps performed by this method:
    
    
""")


# image = Image.open('clustering_metodo_manu.png')
# st.image(image, caption='Clustering Workflow')


########### OPTIONS #######
# SIDEBAR

# Loading file
st.sidebar.header('Upload your fasta')

archivo = st.sidebar.file_uploader("Upload a fasta file with sequences", type=["fasta"])


##########################################DE ACÁ PARA ABAJO NO HAY QUE TOCAR NADA################################################

def druggability_app(archivo):
    seq=[]
    proteinas = []
    for record in SeqIO.parse(archivo,'fasta'):
        seq.append(record.id)
        sequence = str(record.seq)
        sequence = sequence.strip()
        proteinas.append(sequence)
    
    df=pd.DataFrame()
    df2=pd.DataFrame()
    valores=pd.DataFrame()
    a=Pyprotein.PyProtein(str(proteinas[0]))
    a=a.GetCTD()
    b=a.keys()

    for i in b:
    	claves = pd.DataFrame(columns=[i])
    	df = pd.concat([df,claves],sort=False)
    
    for i in proteinas:
    	valor = pd.DataFrame(columns=[i])
    	df2 = pd.concat([df2,valor],sort=False)
    
    df2=df2.T
    
    for i in range(len(proteinas)):
    	protein_class = Pyprotein.PyProtein(str(proteinas[i]))
    	CTD = protein_class.GetCTD()
    	valor = pd.DataFrame(CTD.values())
    	valor = valor.T
    	valores=pd.concat([valores,valor])
    	valores=valores.rename(index={0:str(seq[i])})
    
    b=list(b)
    for i in range(len(b)):
    	valores=valores.rename(columns={i:str(b[i])})
        
    valores=DataFrame(valores[['_SolventAccessibilityD1100','_SecondaryStrD1025','_ChargeD1100','_SolventAccessibilityD1050','_HydrophobicityD3025',
    			  '_SecondaryStrD3001','_ChargeD2075','_PolarityD1001','_PolarizabilityD2050','_ChargeD3100',
                              '_SolventAccessibilityD1100','_ChargeD2075','_NormalizedVDWVD3025','_NormalizedVDWVD3075','_PolarizabilityD3100',
                              '_SolventAccessibilityD1100','_ChargeD2050','_NormalizedVDWVD3075','_SolventAccessibilityD1075','_PolarizabilityD3025']])
    coeficientes=[0.160,0.038,-0.027,-0.037,0.025,
                  -0.109,-0.096,-0.089,0.028,-0.019,
                  0.174,-0.091,0.026,-0.034,0.044,
                  0.160,-0.054,-0.027,-0.040,0.016]
    descriptores=[]
    cuentas=[]
    
    for j in range(len(valores.columns)):
    	v=valores.iloc[:,j].tolist()
    	descriptores.append(v)
    	v=[numero*coeficientes[j] for numero in v]
    	cuentas.append(v)
    
    M1=[]
    M2=[]
    M3=[]
    M4=[]
    
    z=0
    while z<len(valores.index):
    	m1=cuentas[0][z]+cuentas[1][z]+cuentas[2][z]+cuentas[3][z]+cuentas[4][z]-12.232
    	m2=cuentas[5][z]+cuentas[6][z]+cuentas[7][z]+cuentas[8][z]+cuentas[9][z]+8.400
    	m3=cuentas[10][z]+cuentas[11][z]+cuentas[12][z]+cuentas[13][z]+cuentas[14][z]-12.245
    	m4=cuentas[15][z]+cuentas[16][z]+cuentas[17][z]+cuentas[18][z]+cuentas[19][z]-8.019
    	
    	M1.append(m1)
    	M2.append(m2)
    	M3.append(m3)
    	M4.append(m4)
    	z=z+1
    
    corte=0.5
    scores=(M1,M2,M3,M4)
    		
    predicciones=[] #Lista de listas, cada elemento de la lista es uno de los modelos. Cada uno de los modelos tiene las train predicciones binarizadas.
    for models in scores:
      for j in models:#Transformo a valor de clase a cada una de las predicciones de los modelos.
          float(j)
          if j<corte:
            j='Non-druggable'
            predicciones.append(j)
          else: 
            j='Druggable'
            predicciones.append(j)	
    
    
    predicciones=[predicciones[j:j+len(valores.index)] for j in range (0,len(predicciones),len(valores.index))]#Separo la lista original en (n°modelos) listas
    
    resultados=DataFrame(valores.index,columns=['Protein_ID'])
    l=1
    for x in predicciones:
    	g=DataFrame(x,columns=['PREDICT_' + str(l)])
    	resultados=pd.concat([resultados,g],axis=1)
    	l=l+1
    
    train=pd.read_csv('train_todo.csv',sep=',')
    
    train=DataFrame(train[['_SolventAccessibilityD1100','_SecondaryStrD1025','_ChargeD1100','_SolventAccessibilityD1050','_HydrophobicityD3025',
    			  '_SecondaryStrD3001','_ChargeD2075','_PolarityD1001','_PolarizabilityD2050','_ChargeD3100',
                              '_SolventAccessibilityD1100','_ChargeD2075','_NormalizedVDWVD3025','_NormalizedVDWVD3075','_PolarizabilityD3100',
                              '_SolventAccessibilityD1100','_ChargeD2050','_NormalizedVDWVD3075','_SolventAccessibilityD1075','_PolarizabilityD3025']])
    
    h=[]
    
    for i in range(4):
    	m=train.iloc[:,i*5:i*5+5]
    	x=m.to_numpy()
    	xt=np.transpose(x)
    	y=np.matmul(xt,x)
    	z=np.linalg.inv(y)
    	for f in range(len(valores.index)):
    		m_=valores.iloc[f,i*5:i*5+5]
    		x_=m_.to_numpy()
    		xt_=np.transpose(x_)
    		y_=np.matmul(x_,z)
    		z_=np.matmul(y_,xt_)
    		hi=z_
    		h.append(hi)
    h2=10/156 # 2k/n
    h3=15/156 # 3k/n
    
    dominio=[]
    
    for g in h:
    	if g<h2:
    		w='YES**'
    		dominio.append(w)
    	else:
    		if g<h3:
    			w='YES*'
    			dominio.append(w)
    		else:
    			w='NO'
    			dominio.append(w)
    
    dominio=[dominio[j:j+len(valores.index)] for j in range (0,len(predicciones*len(valores.index)),len(valores.index))]#Separo la lista original en (n°modelos) listas
    v=0
    dominio=DataFrame(dominio).T
    for y in predicciones:
    	dominio=dominio.rename(columns={v:'DA_'+str(v+1)})
    	v=v+1
    resultados_final = pd.concat([resultados,dominio],axis=1)
    resultados_final['name'] = seq
    columnas=["Protein_ID","PREDICT_1","DA_1","PREDICT_2","DA_2","PREDICT_3","DA_3","PREDICT_4","DA_4"]
    resultados_final = resultados_final[columnas]
    st.write(resultados_final)
    
    return resultados_final

def druggability_heatmap(resultados_final):
    resultados_final['-'] = 5
    resultados_final['--'] = 5
    resultados_final['---'] = 5
    resultados_final['----'] = 5

    cols = ["Protein_ID","-","PREDICT_1","DA_1","--","PREDICT_2","DA_2","---","PREDICT_3","DA_3","----","PREDICT_4","DA_4"]
    resultados_final = resultados_final[cols]
    resultados_final.replace("Non-druggable","0",inplace=True)
    resultados_final1 = resultados_final.replace("Druggable","10")
    resultados_final1 = resultados_final1.replace("YES*","7")
    resultados_final1 = resultados_final1.replace("YES**","10")
    resultados_final1 = resultados_final1.replace("NO","0")
    resultados_final1.set_index('Protein_ID',inplace=True)
    resultados_final = resultados_final1.apply(pd.to_numeric, errors = 'coerce')
    import seaborn as sns
    import matplotlib.pylab as plt
    
    # st.write(resultados_final)
    cmap = sns.diverging_palette(20,130,center='light',as_cmap=True,s=100,l=50,sep=100)
    # cmap = sns.diverging_palette(10,230,sep=1,center='light',n=4)
    
    # plt.figure(figsize=(45,4))
    graph = sns.heatmap(resultados_final[0:], annot=False, square=True, linewidths=1.0, cmap=cmap, fmt=".1f",robust=True,
                        annot_kws={'fontsize':5}, cbar_kws={"orientation": "vertical"}, cbar=False, vmin=0, vmax=10)
    # bottom, top = graph.get_ylim()
    # graph.set_ylim(bottom, top)
    
    # here set the labelsize by 20
    
    plt.yticks(fontsize=5, rotation=0);plt.xticks(fontsize=5, rotation=90)
    plt.ylabel('Targets',fontsize=5);plt.xlabel("",fontsize=5)
    
    # cbar = graph.collections[0].colorbar
    
    # cbar.ax.tick_params(labelsize=10)
    
    plt.tight_layout()
    st.pyplot(plt)


if archivo is not None:
    run = st.sidebar.button("PREDICT")
    if run == True:
        file_name_st = archivo.name
        archivo_ok = io.TextIOWrapper(archivo)
        resultados_final = druggability_app(archivo_ok)
        druggability_heatmap(resultados_final)

st.markdown("""
         **To cite the application, please reference XXXXXXXXX**
         """)
