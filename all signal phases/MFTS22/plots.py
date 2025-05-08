import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img


font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22,}

plt.rc('font', **font)
plt.rcParams['font.size'] = '22'
plt.rcParams['figure.figsize'] = [12, 7]
# Rendering math equations using TeX
# plt.rcParams['text.usetex'] = True
plt.rcParams.update({
  "text.usetex": True,
  "font.family": "Times"
})

data = []
with open('output_mfts.txt','r') as f:    
    lines = f.readlines()
    for line in lines:
        elements = line.strip().split("\t")
        data.append(list(map(float, elements)))   
f.close()
data_idm = []
with open('output_idm_mfts.txt','r') as f:    
    lines = f.readlines()
    for line in lines:
        elements = line.strip().split("\t")
        data_idm.append(list(map(float, elements)))   
f.close()
data_deter = []
with open('output_deter_mfts.txt','r') as f:    
    lines = f.readlines()
    for line in lines:
        elements = line.strip().split("\t")
        data_deter.append(list(map(float, elements)))   
f.close()


pos = [row[1] for row in data]
spd = [row[2] for row in data]
acc = [row[3] for row in data]
time = [i for i in range(len(pos))]

pos_idm = [row[0] for row in data_idm]
spd_idm = [row[1] for row in data_idm]
acc_idm = [row[2] for row in data_idm]
timeIDM = [i for i in range(len(pos_idm))]

posDeter = [row[1] for row in data_deter]
spdDeter = [row[2] for row in data_deter]
accDeter = [row[3] for row in data_deter]
timeDeter = [row[0] for row in data_deter]

# continuation of acc and accDet in plots
acc[len(time)-1] = accDeter[0]

real_cost = [row[10] for row in data]
real_escape_cost = [row[11] for row in data]

kg_min = data[0][12]
kg_max = data[0][13]
kr_min = data[0][14]
kr_max = data[0][15]
tr_switch = data[0][16]

arrb = [row[17] for row in data]
escape_arrb = [row[18] for row in data]

total_arrb = []
for i in range(len(arrb)):
    total_arrb.append(arrb[i] + escape_arrb[i])

# PLOTS
plt.plot(time, pos, linewidth = 3, color = '#1f77b4')
plt.plot(timeIDM, pos_idm, linewidth = 3, color = 'magenta')
plt.plot([0, kg_min], [300, 300], 'y', linewidth = 3, color = 'green')
plt.plot([kg_min, kg_max], [300, 300], 'y', linewidth = 3, color = 'green', linestyle='dashed')
plt.plot([kg_max, kr_min], [300, 300], 'y', linewidth = 3, color = 'red')
plt.plot([kr_min, kr_max], [300, 300], 'y', linewidth = 3, color = 'red', linestyle='dashed')
plt.plot([0, kr_max], [370, 370], 'y', linewidth = 3, color = 'black')
plt.plot(timeDeter, posDeter, linewidth = 3, color = '#1f77b4')
plt.plot([tr_switch, tr_switch], [0, 370], 'y', linewidth = 1, color = 'black', linestyle='dashed')
plt.plot([60, 60], [0, 370], 'y', linewidth = 1, color = 'black', linestyle='dashed')
plt.legend([r'\bf{SDP}', r'\bf{IDM}', r'\bf{Geen Phase}', r'\bf{Uncertain Geen Phase}', r'\bf{Red Phase}', r'\bf{Uncertain Red Phase}'],  loc = 4)
plt.xlabel(r"\bf{Time Steps}", fontweight = 'bold')
plt.ylabel(r"\bf{Position (in m)}", fontweight = 'bold')
plt.savefig(str(pos[0])+'_'+str(spd[0])+'_'+str(tr_switch)+'_'+'position.png', dpi = 300)
plt.show()

plt.plot(time, spd, linewidth = 3, color = '#1f77b4')
plt.plot(timeIDM, spd_idm, linewidth = 3, color = 'magenta')
plt.plot(timeDeter, spdDeter, linewidth = 3, color = '#1f77b4')
plt.legend([r'\bf{SDP}', r'\bf{IDM}'])
plt.xlabel(r"\bf{Time Steps}", fontweight = 'bold')
plt.ylabel(r"\bf{Speed (in m/s)}", fontweight = 'bold')
plt.savefig(str(pos[0])+'_'+str(spd[0])+'_'+str(tr_switch)+'_'+'speed.png', dpi = 300)
plt.show()

plt.plot(time, acc, linewidth = 3, color = '#1f77b4')
plt.plot(timeIDM, acc_idm, linewidth = 3, color = 'magenta')
plt.plot(timeDeter, accDeter, linewidth = 3, color = '#1f77b4')
plt.legend([r'\bf{SDP}', r'\bf{IDM}'])
plt.xlabel(r"\bf{Time Steps}", fontweight = 'bold')
plt.ylabel(r"\bf{Acceleration (in m/s$^2$)}", fontweight = 'bold')
plt.savefig(str(pos[0])+'_'+str(spd[0])+'_'+str(tr_switch)+'_'+'acceleration.png', dpi = 300)
plt.show()


# Traffic signal phases example
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
from scipy import interpolate, ndimage
# read image
img_file = "test.png"
logo = img.imread(img_file)
# Convert image into numpy array
image_arr = np.array(logo)
# Rotate image
angle = 45
image_arr = ndimage.rotate(image_arr, angle, reshape=True)

fig, ax = plt.subplots(figsize = (12, 8))
ax.plot([0, kg_min], [300, 300], 'y', linewidth = 3, color = 'green')
ax.plot([kg_min, kg_max], [300, 300], 'y', linewidth = 3, color = 'green', linestyle='dashed', dashes = (2,3))
ax.plot([kg_max, kr_min], [300, 300], 'y', linewidth = 3, color = 'red')
ax.plot([kr_min, kr_max], [300, 300], 'y', linewidth = 3, color = 'red', linestyle='dashed', dashes = (2,3))
ax.plot([0, kr_max], [370, 370], 'y', linewidth = 3, color = 'black')


plt.ylim([0, 400])
ax.legend([r'\bf{Certain Geen Phase}', r'\bf{Uncertain Geen Extension Phase}', r'\bf{Certain Red Phase}', r'\bf{Uncertain Red Extension Phase}'], loc = 4)
plt.xlabel(r"\bf{Time}", fontweight = 'bold')
plt.ylabel(r"\bf{Position}", fontweight = 'bold')

ax.set_xticks([0, kr_min, kr_max, kg_min, kg_max])
ax.set_yticks([300, 370])

labelsx = [item.get_text() for item in ax.get_xticklabels()]
labelsy = [item.get_text() for item in ax.get_yticklabels()]

labelsy = [r'$x_1$', r'$x_e$']
labelsx = [0, r'$k^R_{min}$', r'$k^R_{max}$', r'$k^G_{min}$', r'$k^G_{max}$']

ax.set_xticklabels(labelsx)
ax.set_yticklabels(labelsy)

ax.xaxis.grid(True)

imagebox = OffsetImage(image_arr, zoom = 0.03)
ab = AnnotationBbox(imagebox, (0, 25), frameon = False)
ax.add_artist(ab)

# ax.annotate("", xy = (6.0, 75), xytext = (3, 50), arrowprops = dict(arrowstyle="fancy"))
# ax.annotate("", xy = (8.0, 72.5), xytext = (5, 50))

# plt.savefig('signal_example.png', dpi = 300)
# plt.show()