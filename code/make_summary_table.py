"""Write ../results/clark_evans_summary.csv  (Table 1 of the paper)."""
import os, csv, numpy as np
import stats_common as sc
HERE=os.path.dirname(os.path.abspath(__file__)); RES=os.path.join(HERE,'..','results')
lon=[];lat=[]
for row in csv.DictReader(open(os.path.join(HERE,'..','data','gulf_salt_diapirs.csv'))):
    lon.append(float(row['lon']));lat.append(float(row['lat']))
lon=np.array(lon);lat=np.array(lat)
def gsub(a,b,c,d):
    m=(lon>=a)&(lon<=b)&(lat>=c)&(lat<=d);return sc.project_lonlat_km(lon[m],lat[m])
X=[];Y=[];IT=[]
for row in csv.DictReader(open(os.path.join(HERE,'..','data','ngb_salt_structures.csv'))):
    X.append(float(row['easting_km']));Y.append(float(row['northing_km']));IT.append(row['internal_type'])
X=np.array(X);Y=np.array(Y);IT=np.array(IT);Pall=np.c_[X,Y]
COMP='Salzdiapir kompressiv ueberpraegt';HAL=np.isin(IT,['Salzdiapir','Salzkissen','Doppelsalinar',COMP])
def box(P,a,b,c,d):
    m=(P[:,0]>=a)&(P[:,0]<=b)&(P[:,1]>=c)&(P[:,1]<=d);return P[m]
rng=np.random.default_rng(101)
rows=[('Gulf of Mexico','Whole basin (all pooled)',gsub(-100,-84,15.5,33)),
 ('Gulf of Mexico','N. Gulf shelf/slope + coast',gsub(-97,-88,26.5,30)),
 ('Gulf of Mexico','Texas-Louisiana shelf core',gsub(-96,-90,27,29)),
 ('Gulf of Mexico','Bay of Campeche',gsub(-96,-90,16,22)),
 ('Gulf of Mexico','Onshore interior basins',gsub(-96,-88,30,33)),
 ('North German Basin','All structures',Pall),
 ('North German Basin','All halokinetic',Pall[HAL]),
 ('North German Basin','Central basin halokinetic',box(Pall[HAL],350,650,5850,6090)),
 ('North German Basin','Glueckstadt Graben halokinetic',box(Pall[HAL],460,600,5900,6070)),
 ('North German Basin','Salt diapirs only (mature)',Pall[IT=="Salzdiapir"]),
 ('North German Basin','Salt pillows only',Pall[IT=="Salzkissen"])]
out=[['basin','subset','N','median_NND_km','clark_evans_R','p_dispersion','verdict']]
for basin,name,P in rows:
    R,p,_,med,_=sc.clark_evans(P,S=499,rng=rng)
    v='clustered' if R<0.95 else('regular' if R>1.05 else 'random')
    out.append([basin,name,len(P),f'{med:.1f}',f'{R:.3f}',f'{p:.3f}',v])
with open(os.path.join(RES,'clark_evans_summary.csv'),'w',newline='') as f:
    csv.writer(f).writerows(out)
print('wrote ../results/clark_evans_summary.csv')
for r in out: print(r)
