#!/usr/bin/env python3
"""Symbolize samples.txt / maps.txt written by sigprof_sampler.rs.

Both files are written to the directory the sampler was run in; run this
script from that directory.

    python3 symbolize_samples.py /path/to/profiled/binary > perf.txt
    python3 perf_to_svg.py out.svg "title" < perf.txt

Frame addresses inside the binary are resolved with addr2line (inline frames
expanded); addresses in shared libraries (libm, libc) are labelled with the
library only, since their internal functions have no exported symbol. Output
is `perf script -F ip,sym` format, leaf frame first. Standard library only.
"""
import subprocess,re,bisect,sys
maps=[]
for l in open('maps.txt'):
    p=l.split(None,5)
    if len(p)<6 or 'x' not in p[1]: continue
    a,b=[int(x,16) for x in p[0].split('-')]
    maps.append((a,b,int(p[2],16),p[5].strip()))
import os
exe=os.path.realpath(sys.argv[1])
exebase=min(int(l.split('-')[0],16) for l in open('maps.txt') if l.strip().endswith(exe))
L=open('samples.txt').read().split('\n')
stacks=[[int(x,16) for x in l.split()] for l in L[1:] if l.strip()]
stacks=[[s[0]]+[a-1 for a in s[1:]] for s in stacks]
uniq=sorted({a for s in stacks for a in s})
def find(a):
    for m in maps:
        if m[0]<=a<m[1]: return m
syms={}
def dsyms(path):
    if path not in syms:
        o=subprocess.run(['nm','-D','--defined-only','-C',path],capture_output=True,text=True).stdout
        v=sorted((int(l.split()[0],16),l.split(None,2)[2]) for l in o.split('\n') if len(l.split(None,2))==3)
        syms[path]=([x[0] for x in v],[x[1] for x in v])
    return syms[path]
exe_addrs=[];res={}
for a in uniq:
    m=find(a)
    if m is None: res[a]=['[unknown]']
    elif m[3]==exe: exe_addrs.append(a)
    else:
        res[a]=['['+m[3].split('/')[-1]+']']
inp="\n".join(f"{a-exebase:x}" for a in exe_addrs)
out=subprocess.run(['addr2line','-i','-f','-C','-a','-e',exe],input=inp,capture_output=True,text=True).stdout.split('\n')
cur=None;i=0
while i<len(out):
    l=out[i]
    if l.startswith('0x'): cur=int(l,16)+exebase;res[cur]=[];i+=1
    elif l=='': i+=1
    else: res[cur].append(l);i+=2
def simp(n):
    # drop generic args ::<...> and truncate
    prev=None
    while prev!=n:
        prev=n; n=re.sub(r'::<[^<>]*>','',n)
    n=re.sub(r'<(\w[\w:]*)[^<>]*>(?=::)',r'\1',n) if len(n)>110 else n
    return n if len(n)<=110 else n[:107]+'...'
f=sys.stdout
if True:
    for s in stacks:
        for a in s:
            for fn in res.get(a,['[unknown]']):
                f.write(f"\t{a:x} {simp(fn)}\n")
        f.write("\n")
