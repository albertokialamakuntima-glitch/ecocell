import json,re,glob,numpy as np,pandas as pd,shapely,pycountry
import pycountry_convert as pc
from shapely.geometry import shape,mapping,MultiPolygon,box
from shapely.ops import transform
from shapely.strtree import STRtree
PT={'Algeria': 'Argélia', 'Benin': 'Benim', 'Botswana': 'Botsuana', 'Cameroon': 'Camarões', 'Central African Rep.': 'Rep. Centro-Africana', 'Chad': 'Chade', 'Comoros': 'Comores', 'Dem. Rep. Congo': 'RD Congo', 'Djibouti': 'Jibuti', 'Egypt': 'Egipto', 'Eq. Guinea': 'Guiné Equatorial', 'Eritrea': 'Eritreia', 'eSwatini': 'Essuatíni', 'Ethiopia': 'Etiópia', 'Gabon': 'Gabão', 'Gambia': 'Gâmbia', 'Ghana': 'Gana', 'Guinea': 'Guiné', 'Guinea-Bissau': 'Guiné-Bissau', "Côte d'Ivoire": 'Costa do Marfim', 'Kenya': 'Quénia', 'Lesotho': 'Lesoto', 'Liberia': 'Libéria', 'Libya': 'Líbia', 'Madagascar': 'Madagáscar', 'Malawi': 'Maláui', 'Mauritania': 'Mauritânia', 'Mauritius': 'Maurícia', 'Morocco': 'Marrocos', 'Mozambique': 'Moçambique', 'Namibia': 'Namíbia', 'Niger': 'Níger', 'Nigeria': 'Nigéria', 'Rwanda': 'Ruanda', 'São Tomé and Principe': 'São Tomé e Príncipe', 'Seychelles': 'Seicheles', 'Sierra Leone': 'Serra Leoa', 'Somalia': 'Somália', 'Somaliland': 'Somalilândia', 'South Africa': 'África do Sul', 'S. Sudan': 'Sudão do Sul', 'Sudan': 'Sudão', 'Tanzania': 'Tanzânia', 'Tunisia': 'Tunísia', 'W. Sahara': 'Saara Ocidental', 'Zambia': 'Zâmbia', 'Zimbabwe': 'Zimbábue'}
d=pd.read_csv('cells.csv').drop_duplicates(['lat','lon']).reset_index(drop=True)
M=json.load(open('meta.json'));from datetime import date;DAYS=(date.fromisoformat(M['end'])-date.fromisoformat(M['start'])).days
if 'recent' not in d: d['recent']=np.nan
print('células fundidas',len(d))
def fix(s):
    parts=list(s.geoms) if s.geom_type=='MultiPolygon' else [s]; out=[]
    for p in parts:
        b=p.bounds
        if b[2]-b[0]<=180: out.append(p); continue
        q=transform(lambda x,y,z=None:(np.where(np.asarray(x)<0,np.asarray(x)+360,x),y),p).buffer(0)
        for bx,sh in ((box(0,-90,180,90),0),(box(180,-90,360,90),-360)):
            r=q.intersection(bx)
            if r.is_empty: continue
            r=transform(lambda x,y,z=None:(np.asarray(x)+sh,y),r)
            out+=[x for x in (r.geoms if hasattr(r,'geoms') else [r]) if x.geom_type=='Polygon' and x.area>0]
    return MultiPolygon(out)
g=json.load(open('world.json')); feats=[f for f in g['features'] if f['properties']['name']!='Antarctica']
shp=[fix(shape(f['geometry'])) for f in feats]; names=[f['properties']['name'] for f in feats]
# continentes
CONT=[('África',[-18,-35,52,38],8),('Europa',[-25,34,45,72],3),('Ásia',[25,-11,180,78],30),('América do Norte',[-170,7,-50,83],30),('América do Sul',[-92,-56,-34,13],10),('Oceânia',[110,-48,180,-1],20)]
M={'AF':0,'EU':1,'AS':2,'NA':3,'SA':4,'OC':5}; ASIA={'RU','KZ','GE','AM','AZ','CY','TR'}
def cont_of(f,s):
    try:
        a2=pycountry.countries.get(numeric=str(f.get('id')).zfill(3)).alpha_2
        return 2 if a2 in ASIA else M[pc.country_alpha2_to_continent_code(a2)]
    except Exception:
        rp=s.representative_point()
        for k in [0,1,2,3,4,5]:
            b=CONT[k][1]
            if b[0]<=rp.x<=b[2] and b[1]<=rp.y<=b[3]: return k
        return 2
ARID=set("Algeria,Libya,Egypt,Niger,Mali,Chad,Mauritania,W. Sahara,Sudan,Saudi Arabia,Iraq,Iran,Afghanistan,Pakistan,Oman,Yemen,United Arab Emirates,Kuwait,Qatar,Jordan,Syria,Israel,Palestine,Turkmenistan,Uzbekistan,Mongolia,Kazakhstan,Australia,Namibia,Botswana,Somalia,Somaliland,Djibouti,Eritrea,Morocco,Tunisia,Bahrain,Lebanon".split(','))
SEMI=set("China,India,South Africa,Kenya,Ethiopia,Argentina,Turkey,Spain,Burkina Faso,S. Sudan,Zimbabwe,Lesotho,Tajikistan,Kyrgyzstan,Azerbaijan,Armenia,Greece,Portugal,Chile,Paraguay,Uganda,Tanzania,Zambia,Angola,Mozambique,Madagascar,Senegal,Gambia,Guinea-Bissau,Benin,Togo,Nigeria".split(','))
FOREST=set("Brazil,Colombia,Venezuela,Peru,Ecuador,Bolivia,Suriname,Guyana,Dem. Rep. Congo,Congo,Gabon,Cameroon,Central African Rep.,Eq. Guinea,Indonesia,Malaysia,Papua New Guinea,Myanmar,Laos,Cambodia,Vietnam,Thailand,Russia,Sweden,Finland,Norway,Canada,Liberia,Sierra Leone,Guinea,Côte d'Ivoire,Ghana,Japan,Philippines,Panama,Costa Rica,Nicaragua,Honduras,Guatemala,Belize".split(','))
def vclass(n): return 4 if n=='Greenland' else 0 if n in ARID else 1 if n in SEMI else 3 if n in FOREST else 2
K=[cont_of(f,s) for f,s in zip(feats,shp)]
# atribuir células
pts=shapely.points(d.lon.values,d.lat.values); tree=STRtree(shp)
ci=np.full(len(d),-1); sea=np.zeros(len(d),int)
a,b=tree.query(pts,predicate='within'); ci[a[::-1]]=b[::-1]
rest=np.where(ci<0)[0]
a,b=tree.query_nearest(pts[rest],max_distance=2,all_matches=False); ci[rest[a]]=b; sea[rest[a]]=1
d['ci']=ci; d['sea']=sea
d['an']=d['mean']-d.groupby('lat')['mean'].transform('median')
d['chg']=d.slope*DAYS/d['mean']*100
d=d[d.ci>=0].reset_index(drop=True)
d['ch']=d.chg.where(d.chg.abs()<=3)
z=lambda s:(s-s.mean())/s.std()
d['sc']=z(d.an)+z(d.ch).fillna(0); d['rk']=d.sc.rank(ascending=False).astype(int)
q1,q2=d.an.quantile([1/3,2/3]); d['lv']=np.where(d.an<=q1,0,np.where(d.an<=q2,1,2))
def mainb(s):
    ps=list(s.geoms) if s.geom_type=='MultiPolygon' else [s]; wa=lambda p:p.area*np.cos(np.radians(p.centroid.y)); m=max(wa(p) for p in ps)
    bs=np.array([p.bounds for p in ps if wa(p)>=.25*m])
    return [max(-179.9,round(float(bs[:,0].min())-.5,1)),max(-84,round(float(bs[:,1].min())-.5,1)),min(179.9,round(float(bs[:,2].max())+.5,1)),min(84,round(float(bs[:,3].max())+.5,1))]
def expected(s):
    x0,y0,x1,y1=s.bounds; X,Y=np.meshgrid(np.arange(np.floor(x0)+.5,x1,1.0),np.arange(np.ceil(y0),y1,1.0))
    return int(shapely.contains_xy(s,X.ravel(),Y.ravel()).sum())
cs=[]
for i,n in enumerate(names):
    s=d[d.ci==i]; l=s[s.sea==0]
    cs.append({'n':PT.get(n,n),'k':K[i],'v':vclass(n),'b':mainb(shp[i]),'ne':expected(shp[i]) if len(d[d.ci==i]) else 0,'a':round(shp[i].area,1),'nc':len(s),'nl':len(l),
      'ma':None if l.empty else round(float(l.an.mean()),1),'hi':int(round(100*(s.lv==2).mean())) if len(s) else 0})
top=[];cc={}
for ix in d.sort_values('rk').index:
    c=int(d.ci[ix])
    if pd.isna(d.ch[ix]) or abs(d.lat[ix])>60: continue
    if cc.get(c,0)<2: top.append(int(ix)); cc[c]=cc.get(c,0)+1
    if len(top)==30: break
cells=[[r.lat,r.lon,round(r.mean,1),None if pd.isna(r.ch) else round(r.ch,2),round(r.an,1),round(r.sc,2),int(r.lv),int(r.sea),int(r.ci),int(r.rk),None if pd.isna(r.recent) else round(r.recent-r.mean,1)] for r in d.itertuples()]
def rd(x): return [rd(i) for i in x] if isinstance(x,(list,tuple)) else round(x,2)
gf=[]
for i,(s,n) in enumerate(zip(shp,names)):
    m=mapping(s.simplify(0.03 if cs[i]['nc'] and K[i]==0 else 0.1)); big=max(s.geoms,key=lambda p:p.area); rp=big.representative_point()
    gf.append({'type':'Feature','properties':{'i':i,'lon':round(rp.x,2),'lat':round(rp.y,2)},'geometry':{'type':m['type'],'coordinates':rd(m['coordinates'])}})
data={'cells':cells,'countries':cs,'geo':{'type':'FeatureCollection','features':gf},'top':top,'thr':[round(q1,1),round(q2,1)],
      'conts':[{'n':n,'b':b,'lt':lt} for n,b,lt in CONT]}
h=open('template.html').read().replace('__DATA__',json.dumps(data,ensure_ascii=False,separators=(',',':')))
open('index.html','w').write(h)
print(len(cells),'células',len(h)//1024,'KB','países com dados',sum(c['nc']>0 for c in cs),'mar',int(d.sea.sum()),'fiável',round(d.ch.notna().mean(),2),q1,q2)
print(d.sort_values('rk').head(10).assign(c=lambda x:x.ci.map(lambda i:names[i]))[['c','lat','lon','sea','an','ch','sc']].round(2).to_string())
print('sem dados (grandes):',[(names[i],round(c['a'])) for i,c in enumerate(cs) if c['nc']==0 and c['a']>40])
