"""Re-render every original comparison-report figure on exact square canvases.

Uses saved 4000 W results only; no CFD/tracer calculations are rerun.
"""
from pathlib import Path
import json,csv,hashlib,zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PIL import Image
from visualize_compare import Field,boxes,vertices,baffle,FACES

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'deliverables/figures_1to1'
OUT.mkdir(parents=True,exist_ok=True)
CASES=['both_open','both_closed','left_open','right_open']
LAYOUTS=['corner','center']
LABEL={'both_open':'兩門全開','both_closed':'兩門全關','left_open':'左開右關','right_open':'右開左關','corner':'靠左方案','center':'置中方案'}
DOOR_COLORS=['#308bb1','#bf633b','#62946b','#815bad']
LAYOUT_COLORS=['#2589a1','#df9450']
plt.rcParams.update({'font.family':'Noto Sans CJK JP','font.size':13,
                    'axes.unicode_minus':False,'pdf.fonttype':42,
                    'svg.fonttype':'path','figure.facecolor':'white'})
CONFIG={k:json.loads((ROOT/k/'config.json').read_text()) for k in LAYOUTS}
INDEX=[]
COMBINED=OUT/'all_figures_1to1.pdf'

def save(fig,name,title,sources,pdf):
    fig.set_size_inches(12,12)
    for suffix in ['png','pdf','svg']:
        fig.savefig(OUT/(name+'.'+suffix),dpi=300,bbox_inches=None,facecolor='white')
    pdf.savefig(fig,bbox_inches=None)
    INDEX.append(dict(name=name,title=title,sources=sources,canvas_pixels=[3600,3600]))
    plt.close(fig)

def geometry(pdf):
    fig=plt.figure(figsize=(12,12))
    for i,layout in enumerate(LAYOUTS):
        cfg=CONFIG[layout];ax=fig.add_subplot(2,1,i+1,projection='3d')
        for label,b,color in boxes(cfg):
            p=vertices(b)
            ax.add_collection3d(Poly3DCollection(p[FACES],facecolor=color,
                         edgecolor='#607d8b',linewidth=.65,alpha=.34))
        bp=baffle(cfg)
        fs=np.array([[0,1,2],[0,2,3],[4,6,5],[4,7,6],[0,4,5],[0,5,1],
                     [1,5,6],[1,6,2],[2,6,7],[2,7,3],[3,7,4],[3,4,0]])
        ax.add_collection3d(Poly3DCollection(bp[fs],facecolor='#e5a037',
                            edgecolor='#a66a16',alpha=.95,linewidth=.7))
        burner=np.array(cfg['burners'])
        ax.scatter(burner[:,0],burner[:,1],burner[:,2]+.012,s=45,
                   color='#d35532',depthshade=False)
        ax.set(xlim=(0,3.3),ylim=(0,1.65),zlim=(0,2.2),
               xlabel='x（m）',ylabel='y（m）',zlabel='z（m）')
        ax.set_xticks([0,.7,1.5,2.5,3.3]);ax.set_yticks([.6,1.05,1.65]);ax.set_zticks([0,.9,1.6,2.2])
        ax.tick_params(labelsize=10,pad=1)
        ax.set_box_aspect((3.3,1.65,2.2));ax.view_init(elev=20,azim=-62)
        ax.set_title(f"{LABEL[layout]}｜煙機與雙爐整體中心 x＝{cfg['layout_center_x_m']:.2f} m",fontsize=15,pad=4)
    fig.subplots_adjust(left=.03,right=.95,bottom=.05,top=.89,hspace=.08)
    fig.suptitle('廚房幾何比較｜靠左與置中\n頂吸煙機90×30 cm；前擋板10 cm、向下45°',fontsize=19,y=.97)
    fig.text(.5,.022,'橙色擋板只在煙機90 cm寬度內，隨煙機與爐頭一起移動；各爐對空氣放熱2000 W。',ha='center',fontsize=11)
    save(fig,'01_geometry_comparison_1to1','幾何模型：靠左／置中',
         [f'{k}/config.json' for k in LAYOUTS],pdf)

def spatial(pdf):
    u=np.linspace(0,3.3,65);v=np.linspace(0,1.65,39);xx,yy=np.meshgrid(u,v)
    q=np.c_[xx.ravel(),yy.ravel(),np.full(xx.size,1.5)]
    for n,case in enumerate(CASES,2):
        fig,axes=plt.subplots(2,1,figsize=(12,12))
        sources=[]
        for ax,layout in zip(axes,LAYOUTS):
            path=ROOT/layout/'results'/(case+'_refined_tracer.npz')
            with np.load(path) as d:field=Field(dict(d))
            ci=int(np.flatnonzero(np.isclose(field.times,120))[0])
            values,_=field.at(q);c=values[:,6+ci].reshape(xx.shape)
            image=ax.pcolormesh(u,v,np.maximum(c,1e-7),cmap='magma',
                               norm=LogNorm(1e-7,10),shading='auto',rasterized=True)
            ax.add_patch(Rectangle((0,0),.7,.6,facecolor='#becbd2',edgecolor='#69818d'))
            hood=CONFIG[layout]['hood']
            ax.add_patch(Rectangle((hood['x0'],0),.9,.3,fill=False,
                                  edgecolor='#54cad1',linestyle='--',linewidth=1.8))
            for bx,by,bz in CONFIG[layout]['burners']:
                ax.plot(bx,by,'o',color='#5bd5d1',markersize=5)
            metric=json.loads((ROOT/layout/'results'/(case+'_refined_tracer_metrics.json')).read_text())
            ax.set_title(f"{LABEL[layout]}｜{LABEL[case]}｜120 s滯留 {100*metric['retained_fraction_at_120s']:.2f}%",fontsize=15,pad=12)
            ax.set(xlim=(0,3.3),ylim=(0,1.65),xlabel='x（m）',ylabel='y（m）')
            ax.set_aspect('equal');sources.append(str(path.relative_to(ROOT)))
        fig.subplots_adjust(left=.09,right=.83,bottom=.145,top=.85,hspace=.35)
        cb=fig.colorbar(image,cax=fig.add_axes([.88,.15,.025,.65]))
        cb.set_label('歸一化示蹤物濃度 c/R（s/m³）',labelpad=13)
        fig.suptitle(f'油煙分布比較｜{LABEL[case]}\n釋放120 s；水平剖面 z＝1.50 m',fontsize=19,y=.965)
        fig.text(.46,.027,'兩位置使用相同色階（對數）；雙爐合計4000 W。\n數值按總排放率歸一化，並非實測油煙濃度。',ha='center',fontsize=12)
        save(fig,f'{n:02d}_smoke_{case}_1to1',f'油煙空間分布：{LABEL[case]}',sources,pdf)

def retention(pdf):
    fig,axes=plt.subplots(2,1,figsize=(12,12));x=np.arange(4);sources=[]
    for ax,(key,title,factor) in zip(axes,[
        ('steady_equivalent_residence_time_s','全室平均滯留時間 τ（s）',1),
        ('occupied_zone_normalized_concentration_at_120s','走道濃度指數（c/R ×1000）',1000)]):
        for i,layout in enumerate(LAYOUTS):
            vals=[]
            for case in CASES:
                path=ROOT/layout/'results'/(case+'_refined_tracer_metrics.json')
                vals.append(json.loads(path.read_text())[key]*factor)
                sources.append(str(path.relative_to(ROOT)))
            ax.bar(x+(i-.5)*.36,vals,.36,label=LABEL[layout],color=LAYOUT_COLORS[i])
        ax.set(xticks=x,xticklabels=[LABEL[c] for c in CASES],ylabel=title)
        ax.set_title(title,fontsize=16,pad=12)
        ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True);ax.legend(loc='upper left',fontsize=12)
        ax.margins(y=.17)
    fig.subplots_adjust(left=.13,right=.95,bottom=.10,top=.86,hspace=.38)
    fig.suptitle('兩種位置的油煙滯留指標比較\n雙爐合計4000 W；相同油煙源強',fontsize=19,y=.965)
    fig.text(.5,.035,'全室滯留與走道濃度是不同指標；走道濃度指數不是實際mg/m³。',ha='center',fontsize=12)
    save(fig,'06_retention_exposure_comparison_1to1','全室滯留／走道濃度指標',sorted(set(sources)),pdf)

def decay(pdf):
    fig,axes=plt.subplots(2,2,figsize=(12,12));sources=[]
    for col,layout in enumerate(LAYOUTS):
        for row,(index,end,ylim,ylabel) in enumerate([
            (1,180,(1e-5,1),'全室存量／停止排放時存量'),
            (7,420,(1e-4,1),'走道平均濃度／停止排放時濃度')]):
            ax=axes[row,col]
            for case,color in zip(CASES,DOOR_COLORS):
                path=ROOT/layout/'results'/(case+'_refined_tracer_history.csv')
                a=np.loadtxt(path,delimiter=',',skiprows=1)
                at120=a[np.argmin(np.abs(a[:,0]-120)),index]
                assert at120>0
                sel=(a[:,0]>=120)&(a[:,0]<=end)
                ax.semilogy(a[sel,0]-120,np.maximum(a[sel,index]/at120,1e-8),
                            color=color,label=LABEL[case],linewidth=1.8)
                sources.append(str(path.relative_to(ROOT)))
            ax.axhline(.05,color='#777777',linestyle='--',linewidth=1.2)
            ax.set(xlim=(0,end-120),ylim=ylim,xlabel='停止示蹤物排放後時間（s）',ylabel=ylabel)
            kind='全室存量衰減' if row==0 else '走道尾端衰減'
            ax.set_title(f'{LABEL[layout]}｜{kind}',fontsize=14,pad=12)
            ax.grid(alpha=.2,which='major');ax.tick_params(labelsize=11)
    handles,labels=axes[0,0].get_legend_handles_labels()
    from matplotlib.lines import Line2D
    handles.append(Line2D([0],[0],color='#777777',linestyle='--'));labels.append('剩餘5%')
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.895),ncol=5,fontsize=11,frameon=False)
    fig.subplots_adjust(left=.115,right=.96,bottom=.115,top=.825,hspace=.42,wspace=.33)
    fig.suptitle('全室存量與走道尾端衰減｜靠左 vs 置中\n四幅合併；各位置均比較四種門狀態',fontsize=19,y=.965)
    fig.text(.5,.032,'120 s後只停止示蹤物排放，風機與4000 W熱負荷保持不變。\n走道尾端是走道平均濃度的長尾衰減；並非走廊末端單點，亦非實際關火清除秒數。',ha='center',fontsize=11)
    save(fig,'07_clearance_four_panels_1to1','全室存量＋走道尾端：四幅合一',sorted(set(sources)),pdf)

def main():
    with PdfPages(COMBINED) as pdf:
        geometry(pdf);spatial(pdf);retention(pdf);decay(pdf)
    for item in INDEX:
        with Image.open(OUT/(item['name']+'.png')) as im:
            assert im.size==(3600,3600)
    from pypdf import PdfReader
    pages=PdfReader(COMBINED).pages
    assert len(pages)==7
    assert all(abs(float(p.mediabox.width)-float(p.mediabox.height))<1e-8 for p in pages)
    summary=dict(original_report=str(ROOT/'deliverables/kitchen_smoke_comparison_report.pdf'),
                 count=7,canvas_ratio='1:1',png_pixels=[3600,3600],dpi=300,
                 pdf_page_points=[864,864],source_heat_total_W=4000,
                 new_simulations_run=False,figures=INDEX,
                 decay_layout=dict(columns=['靠左','置中'],rows=['全室存量','走道平均濃度長尾']),
                 spatial_layout='each door state: corner above, center below',
                 validation=dict(all_png_exact_square=True,all_pdf_pages_square=True,pdf_pages=7))
    (OUT/'figure_index.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
    (OUT/'README_zh.md').write_text('''# 中文優劣分析PDF：全部圖表正方形重輸出

共7張，畫布比例1:1，PNG為3600×3600像素、300 dpi；各圖另有PDF及SVG。all_figures_1to1.pdf含全部7張正方形圖頁。

01幾何比較；02–05為四種門狀態的油煙分布，上靠左、下置中；06全室平均滯留與走道濃度；07全室存量與走道尾端四幅合一，上排全室、下排走道，左靠左、右置中。

只重新绘製已存4000 W的加密網格結果，沒有重新計算流場或示蹤物。保持原色階與座標範圍，房間剖面長寬比不拉伸。走道尾端指原報告走道平均濃度的長尾衰減，並非走廊末端單點濃度。

全室曲線仍顯示停止排放後0–60 s；走道曲線0–300 s，四種門狀態均保留。只停止示蹤物，風機與4000 W熱負荷不變，不代表熄火後實際秒數。

圖表來源映射與輸出規格見figure_index.json。源程式export_square_figures.py保留於工作專案，依賴原計算資料與visualize_compare.py。
''')
    (OUT/'figure_gallery.html').write_text('<!doctype html><html lang="zh-Hant"><meta charset="utf-8"><title>正方形比較圖</title><style>body{font-family:system-ui;max-width:1050px;margin:35px auto;padding:0 20px;color:#183747}img{width:100%}section{margin:35px 0}a{color:#16788d}</style><h1>中文優劣分析：全部7張1:1圖</h1><p>3600×3600 px・300 dpi・沿用原4000 W計算。</p>'+''.join(f'<section><h2>{x["title"]}</h2><p><a href="{x["name"]}.png">PNG</a> · <a href="{x["name"]}.pdf">PDF</a> · <a href="{x["name"]}.svg">SVG</a></p><img src="{x["name"]}.png"></section>' for x in INDEX)+'</html>')
    manifest={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(OUT.iterdir()) if p.is_file() and p.suffix!='.zip' and p.name!='manifest.json'}
    (OUT/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
    archive=OUT.parent/'kitchen_comparison_figures_1to1.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(OUT.iterdir()):
            if p.is_file():z.write(p,p.name)
        z.write(Path(__file__),'export_square_figures.py')
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for name,item in manifest.items():
            assert hashlib.sha256(z.read(name)).hexdigest()==item['sha256']
    print(json.dumps(dict(count=7,archive=str(archive),archive_MiB=archive.stat().st_size/1024**2,validation=summary['validation']),ensure_ascii=False,indent=2))

if __name__=='__main__':main()
