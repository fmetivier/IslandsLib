#
#
# Extraction de la Désirade
#
#


#
#librairies
#
import sys 
sys.path.append("..")

import cartopy.io.shapereader as shapereader
import numpy as np
import matplotlib.pyplot as plt

import IslandsLib as il



#
# 1 extraction des contours dans un dictionnaire en choisissant le niveau des plus hautes eaux
#
fname="/home/metivier/Nextcloud/Recherche/GroundWater/islands/data/Guadeloupe/LimiteTerreMer_GLP.shp"
reader = shapereader.Reader(fname)

features = reader.records()
lg = {}
i = 0
for feature in features:
    # print(feature.geometry)
    print(feature.attributes)
    print("===================== next test is specific and shoud be changed")
    if "hautes" in feature.attributes["NiveauLTM"]:
        lg[feature.attributes["CdOH"]] = il.linestring_to_list( str(feature.geometry) )


#
# Pour la suite je procède par essais et erreurs
#

#on dirait que c'est beau...
fig, ax = plt.subplots(1)
for key, val in lg.items():
    ax.plot(val[0], val[1], color='k')
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.axis("equal")
plt.savefig('../docs/source/figures/fig1.svg', bbox_inches='tight')

#Mais en fait...
fig, ax = plt.subplots(1)
for key, val in lg.items():
    ax.plot(val[0], val[1])
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.axis("equal")
plt.savefig('../docs/source/figures/fig2.svg', bbox_inches='tight')


#on ne garde que La Désirade (facile)
fig, ax = plt.subplots(1)
for key, val in lg.items():
    if min(val[0]) > 7e5 and min(val[1])>1.8e6:
        ax.plot(val[0], val[1])
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.axis("equal")
plt.savefig('../docs/source/figures/fig3.svg', bbox_inches='tight')

# on enlève tous les petites iles fermées 
# Attention cela enlève aussi les lacs et lagunes fermées à vérifier ultérieurement si besoin
#
fig, ax = plt.subplots(1)
i = 0
c = ['C0','C1']
for key, val in lg.items():
    if min(val[0]) > 7e5 and min(val[1])>1.8e6:
        if val[0][0] != val[0][-1] and val[1][0] != val[1][-1]:
            ax.plot(val[0], val[1], color = c[i%2])
            ax.text(val[0][0], val[1][0], key, color = c[i % 2])
            i += 1
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.axis("equal")
plt.savefig('../docs/source/figures/fig4.svg', bbox_inches='tight')

DesiradeList = ['07L0000001169872212','07L0000001169872209']
ex={}
for key, val in lg.items():     
    if key in DesiradeList:
           ex[key] = [val[0][0], val[1][0], val[0][-1], val[1][-1]]

x,y = il.IGN_to_segment_v2(lg, DesiradeList, ex, plot_res=False)
fig, ax = plt.subplots(1)
ax.plot(x,y,'-', color='C1')
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.axis("equal")
plt.savefig('../docs/source/figures/fig5.svg', bbox_inches='tight')

fname = "../data/Contours/Guadeloupe/Desirade.txt"
entete = "#Guadeloupe\n#Désirade\n#RGAF09 / UTM zone 20N\n#Source IGN\n#anticlockwise"
il.output_to_file(fname,entete,x,y)


plt.show()
