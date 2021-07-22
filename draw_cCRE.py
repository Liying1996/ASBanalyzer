import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse

# plt.switch_backend('agg')

parser = argparse.ArgumentParser()
parser.add_argument("-file", help="please give the file path", required=True)
parser.add_argument("-name", help="the title of the figure", required=True)
args = parser.parse_args()

path='/'.join(args.file.split('/')[:-1])
mpl.rcParams['font.size'] = 7.0 
 

data = pd.read_csv(args.file, header=None, sep='\t')
data.replace(',CTCF-bound', '',regex=True, inplace=True)

as_data = data[data.iloc[:,13] == "significant"]
nonas_data = data[data.iloc[:,13] == "not significant"]

x = data.iloc[:,-1].unique().tolist()
x1 = [0 for i in range(len(x))]
for i in as_data.iloc[:,-1]:
    for j in range(len(x)):
        if i == x[j]:
            x1[j] += 1
            continue
x2 =  [0 for i in range(len(x))]
for i in nonas_data.iloc[:,-1]:
    for j in range(len(x)):
        if i == x[j]:
            x2[j] += 1
            continue

x_0 = [1,0,0,0] 

labels = x
color_dict = {'dELS':'#FFCD00', 'pELS':'#FFA700', 'PLS':'#FF0000', 'DNase-H3K4me3':'#ffaaaa', 'CTCF-only':'#00B0F0', 'Unclassified':'#8c8c8c'}
colors = [color_dict[i] for i in labels]

font1 = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 15,
        }
      
plt.figure(figsize=(20,20))
fig, ax = plt.subplots()

# background
rect = fig.patch
rect.set_facecolor('white')

# labels
def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        if val/total < 0.01:
            return ''
        return '{p:.2f}% ({v:d})'.format(p=pct,v=val)
    return my_autopct

# pie
pie_1 = ax.pie(x1,startangle = 90,radius=1.8,pctdistance = 1.0,colors=colors,autopct=make_autopct(x1),textprops = {'fontsize':13, 'color':'k'})
pie_2 = ax.pie(x2,startangle = 90,radius=1.3,pctdistance = 0.9,colors=colors,autopct=make_autopct(x2),textprops = {'fontsize':13, 'color':'k'})
pie_0 = ax.pie(x_0, radius=0.8,colors = 'w')

# title
ax.text(0.1, 2.3, args.name, fontsize=24, style='oblique', ha='center',va='top',wrap=True)

# pie color
for pie_wedge in pie_0[0][:1]:
    pie_wedge.set_edgecolor('gray')
for pie_wedge in pie_1[0]:
    pie_wedge.set_edgecolor('gray')
for pie_wedge in pie_2[0]:
    pie_wedge.set_edgecolor('gray')

ax.legend(labels, bbox_to_anchor=(1.3,1.0), loc='center left', prop=font1)
ax.set(aspect="equal")

fig.savefig(path+'/'+args.name+'_cCREs.png',dpi=200,bbox_inches='tight',facecolor=fig.get_facecolor(), transparent=True)
plt.show()


