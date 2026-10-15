from pathlib import Path
import json, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve, average_precision_score, auc

def root_from_file(file): return Path(file).resolve().parents[1]
def load_config(root): return json.loads((Path(root)/'config/analysis_config.json').read_text())
def pair_key_vec(a,b):
    a=pd.Series(a,dtype='string').astype(str); b=pd.Series(b,dtype='string').astype(str)
    return pd.Series(np.where(a.to_numpy()<=b.to_numpy(), a.to_numpy()+'||'+b.to_numpy(), b.to_numpy()+'||'+a.to_numpy()))
def score_for_auc(x,direction):
    x=np.asarray(x,float); return x if direction=='higher' else -x
def threshold_pred(x,t,rule):
    x=np.asarray(x,float)
    if rule=='>=': return x>=t
    if rule=='<=': return x<=t
    raise ValueError(rule)
def safe_div(a,b): return float(a/b) if b else np.nan
def confusion_metrics(y,p):
    y=np.asarray(y,int); p=np.asarray(p,bool)
    tp=int(((y==1)&p).sum()); fn=int(((y==1)&~p).sum()); tn=int(((y==0)&~p).sum()); fp=int(((y==0)&p).sum())
    sens=safe_div(tp,tp+fn); spec=safe_div(tn,tn+fp); ppv=safe_div(tp,tp+fp); npv=safe_div(tn,tn+fn)
    bal=np.nanmean([sens,spec]); f1=safe_div(2*ppv*sens,ppv+sens) if np.isfinite(ppv) and np.isfinite(sens) and ppv+sens else np.nan
    return dict(tp=tp,fn=fn,tn=tn,fp=fp,sensitivity=sens,specificity=spec,precision=ppv,npv=npv,balanced_accuracy=bal,f1=f1)
def youden_threshold(y,raw,direction,sample_weight=None):
    s=score_for_auc(raw,direction); fpr,tpr,thr=roc_curve(y,s,sample_weight=sample_weight); j=tpr-fpr
    finite=np.isfinite(thr); idx=np.where(finite)[0][np.argmax(j[finite])] if finite.any() else int(np.argmax(j))
    t=thr[idx]; return float(t if direction=='higher' else -t), float(j[idx])
def roc_pr(y,raw,direction):
    s=score_for_auc(raw,direction); roc=float(roc_auc_score(y,s)); ap=float(average_precision_score(y,s)); p,r,_=precision_recall_curve(y,s); pr=float(auc(r,p))
    return roc,ap,pr
def bootstrap_confusion_by_group(genus_table,B=5000,seed=20260914):
    # genus_table: one row/genus with tp/fn/tn/fp
    if len(genus_table)<2: return {}
    a=genus_table[['tp','fn','tn','fp']].to_numpy(float); rng=np.random.default_rng(seed); n=len(a)
    vals={k:[] for k in ['sensitivity','specificity','precision','npv','balanced_accuracy','f1']}
    for _ in range(B):
        idx=rng.integers(0,n,n); tp,fn,tn,fp=a[idx].sum(axis=0)
        sens=safe_div(tp,tp+fn); spec=safe_div(tn,tn+fp); ppv=safe_div(tp,tp+fp); npv=safe_div(tn,tn+fn); bal=np.nanmean([sens,spec]); f1=safe_div(2*ppv*sens,ppv+sens) if np.isfinite(ppv) and np.isfinite(sens) and ppv+sens else np.nan
        for k,v in [('sensitivity',sens),('specificity',spec),('precision',ppv),('npv',npv),('balanced_accuracy',bal),('f1',f1)]:
            if np.isfinite(v): vals[k].append(v)
    out={}
    for k,v in vals.items():
        if v:
            v=np.asarray(v); out[k+'_ci_low']=float(np.quantile(v,.025)); out[k+'_ci_high']=float(np.quantile(v,.975)); out[k+'_bootstrap_n']=len(v)
    return out
def genus_confusion_table(df,pred_col='pred',group_col='genus_eval',y_col='y'):
    rows=[]
    for g,z in df.groupby(group_col,sort=False):
        y=z[y_col].to_numpy(int); p=z[pred_col].to_numpy(bool)
        m=confusion_metrics(y,p); rows.append({'genus':g,'n':len(z),**{k:m[k] for k in ['tp','fn','tn','fp']}})
    return pd.DataFrame(rows)
def bootstrap_macro_mean(values,B=5000,seed=20260914):
    v=np.asarray(pd.Series(values).dropna(),float)
    if len(v)<2: return (np.nan,np.nan,np.nan,0)
    rng=np.random.default_rng(seed); means=np.empty(B,float)
    for i in range(B): means[i]=rng.choice(v,size=len(v),replace=True).mean()
    return float(v.mean()),float(np.quantile(means,.025)),float(np.quantile(means,.975)),len(v)
