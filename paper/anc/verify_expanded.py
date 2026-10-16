#!/usr/bin/env python3
"""Exact checks for the central extension and printed data; no quotient searches.
Requires only sibling verify.py and Python's standard library.
"""
from verify import red, inv, mul, conj, comm, subst, determinant, matmul
from fractions import Fraction as F
import json

R = {'R':'xyXYaxyXYA','T':'aaxAAX','X':'AYaayAxaYAAyay',
     'Y':'AYaayxyaYAAyaXY','S':'AYaayAxyaaYA'}
rows = [('X',1,''),('Y',1,'Y'),('X',-1,'xAYaayAXaYAAya'),
('Y',-1,'xyAYaayAYXYAAya'),('T',-1,'xyXYaxAAX'),('R',-1,'xyXYaxAAXYa'),
('R',1,'xyXYaxAAYX'),('T',1,'xyXYaxAAYX'),('T',1,'xyXYaxAAYaayA'),
('S',1,'xyXYaxA'),('T',-1,'xyXYaxyAAYaayXYAAyaaxAAYX'),
('R',-1,'xyXYaxyAAYaayXYAAyaaxAAYX'),('R',1,'xyXYaxyAAYaayXYAAyaaxAAXYa'),
('T',1,'xyXYaxyAAYaayXYAAyaaxAAX'),('S',-1,'xyXYaxyXAAYaaXAYAAya'),
('T',-1,'xyXYaxyXAAYaaXAA'),('R',1,''),('T',1,'AYaayA'),('R',-1,'AYaayAxa')]
assert mul(*(conj(v,R[r] if e==1 else inv(R[r])) for r,e,v in rows)) == ''
assert [sum(e for r,e,v in rows if r==k) for k in 'RTXYS']==[0,1,0,0,0]
C=[[0,1,1],[1,0,0],[1,1,0]]
L=[[0,-1,0],[0,1,-1],[-1,-1,1]]
assert matmul(C,L)==[[-1,0,0],[0,-1,0],[0,0,-1]]
J=comm('x','y')
assert comm(mul(J,'a'),R['R'])==conj(R['R'],comm('aa',J))
assert mul(J,'axm',J,inv('axm'))==mul(R['R'],conj('a',comm(inv(J),'xm')))
F0=mul(J,comm('r','s')); k=conj('R','m')
assert mul(comm('x','My'),comm('r',mul(inv(k),'s')))==mul(
    conj('xMX',F0),conj('xM',inv(comm('m',inv('rx')))))
assert mul(comm('x','my'),comm('r',mul(k,'s')))==mul(
    conj('xmX',F0),conj('x',comm('m',inv('rx'))))
psi={'x':'Y','y':'yx','r':'r','s':'s'}
for c in 'xy':
    w=c
    for _ in range(6): w=subst(w,psi)
    assert w==conj(inv(J),c)

# The automorphism and exact arithmetic in F(a,b,c,d) semidirect <u>.
phi={'a':'b','b':'c','c':'d','d':'CAcc'}
iphi={'a':'bbDB','b':'a','c':'b','d':'c'}
for c in 'abcd':
    assert subst(phi[c],iphi)==c and subst(iphi[c],phi)==c

def ph(w,n):
    for _ in range(abs(n)): w=subst(w,phi if n>0 else iphi)
    return w

def prod(*pairs):
    w,n='',0
    for v,j in pairs: w,n=mul(w,ph(v,n)),n+j
    return w,n

def inverse(p): return ph(inv(p[0]),-p[1]),-p[1]
def cj(a,b): return prod(a,b,inverse(a))
def cm(a,b): return prod(a,b,inverse(a),inverse(b))
def evaluate(word,images):
    return prod(*(images[c] if c.islower() else inverse(images[c.lower()]) for c in word))
u=('',1); ss=('c',1); yy=cj(inverse(u),ss); zz=prod(u,u)
NN=cm(ss,inverse(zz)); BB=cj(inverse(ss),NN); kk=cm(ss,inverse(BB))
rr=prod(inverse(u),inverse(kk)); xx=cj(inverse(u),rr); MM=cj(rr,kk)
TT=cm(zz,xx); XX=prod(cj(BB,xx),yy)
YY=prod(cj(BB,yy),inverse(prod(cj(BB,MM),yy,xx)))
assert prod(cj(BB,TT),XX)==('',0)
f='AcaCA'
sig=mul(ph(f,1),'B',inv(f),'b')
assert cj(('BAc',1),TT)==(sig,0)
kf='aCbdB'; eta='AcaCbdBCacbDBcAAcDB'
assert kk==(kf,0) and YY==(inv(eta),1)
x0=prod((inv(kf),0),inverse(u))
assert x0==prod(TT,xx)
j0=mul(inv(kf),'a',ph(kf,1),'B')
assert cm(x0,yy)==(j0,0)
assert mul(f,inv(subst(conj('a',f),{'a':'b','c':'d'})))==conj('b',inv(sig))

# Printed labelled action of C2*C3: alpha=p, beta=q; free labels e1=e,e2=f,e3=g.
p=[1,0,17,4,3,8,7,6,5,10,9,27,13,12,28,16,15,2,24,20,19,29,23,22,18,26,25,11,14,21]
q=[25,2,3,1,5,6,4,0,9,16,12,10,11,15,13,14,8,19,17,18,22,20,21,24,26,7,23,28,29,27]
lp={12:'E',13:'e',22:'f',23:'F'}; lq={4:'g',6:'G'}
def action(word,state=0):
    label=''
    for c in word:
        perm,labels=(p,lp) if c.lower()=='p' else (q,lq)
        if c.islower(): label,state=mul(label,labels.get(state,'')),perm[state]
        else:
            state=perm.index(state); label=mul(label,inv(labels.get(state,'')))
    return label,state
for i in range(30):
    assert action('pp',i)==('',i) and action('qqq',i)==('',i)
h='Qp'; v='Pqq'; aa=mul(h,inv(v)); cc=mul(v,h,inv(v),inv(v)); z=mul(v,v,inv(h))
ff=subst(f,{'a':aa,'c':cc}); gg=conj(aa,ff)
assert action(z)==('g',0) and action(ff)==('gE',0) and action(gg)==('fg',0)
syllables=['AcaC','acA','Cac','aCA','cAAc','CA']
assert [action(subst(w,{'a':aa,'c':cc}))[1] for w in syllables]==[3,16,18,21,15,2]
assert action(aa)[1]==3
assert determinant([[-1,1,0,0],[0,0,-1,1],[0,-1,1,0],[-1,0,1,-1]]) in (1,-1)

# Whitehead graph and all 96 cut differences, not a search for quotients.
words=['xyAXYaaxyAXY','aaxAAX','YaayAxaaYAAyyA','YaayxyaYAAyXY','YaayAxyaaYA']
images={'a':'a','x':'ax','y':'yA'}
def cyc(w):
    w=red(w)
    while w and w[0]==inv(w[-1]): w=w[1:-1]
    return w
for old,new in zip(R.values(),words):
    w=cyc(subst(old,images)); assert len(w)==len(new) and new in w+w
for w in words:
    assert all(w != w[:d]*(len(w)//d) for d in range(1,len(w)) if len(w)%d==0)
assert len(set(map(len,words)))==5
edges=[(w[i],inv(w[(i+1)%len(w)])) for w in words for i in range(len(w))]
expected=[
[0,8,15,17,15,19,10,8,2,10,9,11,15,19,2,0],
[0,17,21,26,23,28,24,17,11,24,24,25,18,19,11,0],
[0,7,4,9,2,1,6,3,11,14,7,8,11,6,7,0],
[0,7,6,11,8,7,14,11,3,6,1,2,9,4,7,0],
[0,11,19,18,25,24,24,11,17,24,28,23,26,21,17,0],
[0,2,19,15,11,9,10,2,8,10,19,15,17,15,8,0]]
for v,expect in zip('YXAaxy',expected):
    others=[c for c in 'YXAaxy' if c not in [v,inv(v)]]
    differences=[]
    for mask in range(16):
        cut={v}|{c for j,c in enumerate(others) if mask&(1<<j)}
        differences.append(sum((a in cut)!=(b in cut) for a,b in edges)-sum(a==v or b==v for a,b in edges))
    assert differences==expect
for a in 'YXx':
    for b in 'Aay': assert (a,b) in edges or (b,a) in edges

# Exterior algebra: indices Q1,Q2,Y1,Y2; verify the correction for symbolic f.
def wedge(a,b):
    out={}
    for I,c in a.items():
        for J,d in b.items():
            if set(I)&set(J): continue
            sign=(-1)**sum(i>j for i in I for j in J)
            key=tuple(sorted(I+J)); out[key]=out.get(key,0)+c*d*sign
    return {k:v for k,v in out.items() if v}
AF={(0,2):F(1),(1,2):F(-1,2)}; AB={(1,3):F(1),(1,2):F(1,2)}; normal={(2,3):F(1)}
assert wedge(AF,AF)==wedge(AB,AB)=={}
assert wedge(AF,normal)==wedge(AB,normal)==wedge(normal,normal)=={}
assert wedge(AF,AB)
# Fold explicit subgroup generators in C2*C3, completing only finite-factor
# polygons. This verifies a fixed intersection certificate, not group quotients.
def modular(word):
    stack=[]
    for c in word:
        factor=c.lower(); exponent=1 if c.islower() else -1
        if stack and stack[-1][0]==factor: exponent+=stack.pop()[1]
        exponent%=2 if factor=='p' else 3
        if exponent: stack.append((factor,exponent))
    return ''.join(c if n==1 else 'Q' for c,n in stack)

def fold(words):
    parent=[0]; orbits={'p':[], 'q':[]}
    def new(): parent.append(len(parent)); return len(parent)-1
    def root(i):
        while parent[i]!=i: i=parent[i]
        return i
    for word in words:
        word=modular(word); state=0
        for k,c in enumerate(word):
            end=0 if k==len(word)-1 else new()
            if c=='p': orbits['p'].append((state,end))
            else:
                third=new(); orbits['q'].append((state,end,third) if c=='q' else (end,state,third))
            state=end
    while True:
        conflict=None
        for key,cycles in orbits.items():
            seen={}
            for cycle in cycles:
                cycle=[root(i) for i in cycle]
                for j,a in enumerate(cycle):
                    b=cycle[(j+1)%len(cycle)]
                    if a in seen and seen[a]!=b: conflict=(seen[a],b); break
                    seen[a]=b
                if conflict: break
            if conflict: break
        if not conflict: break
        a,b=sorted(map(root,conflict)); parent[b]=a
    reps=sorted({root(i) for i in range(len(parent))}); indices={v:i for i,v in enumerate(reps)}
    out={key:{} for key in orbits}
    for key,cycles in orbits.items():
        for cycle in cycles:
            cycle=[indices[root(i)] for i in cycle]
            for j,a in enumerate(cycle): out[key][a]=cycle[(j+1)%len(cycle)]
    out['Q']={v:k for k,v in out['q'].items()}
    return out
left=fold([z,ff,gg]); v0=subst('CaaC',{'a':aa,'c':cc})
right=fold([conj(v0,w) for w in [z,ff,gg]])
expected={'p':{},'q':{},'Q':{}}
for a,b in [(0,1),(2,4),(5,7)]: expected['p'][a]=b; expected['p'][b]=a
for cyc3 in [(1,3,2),(4,6,5),(7,9,8)]:
    for i,a in enumerate(cyc3): expected['q'][a]=cyc3[(i+1)%3]
expected['Q']={v:k for k,v in expected['q'].items()}
match={(0,0):0}; queue=[(0,0)]
for a,b in queue:
    for key in ['p','q','Q']:
        exists=a in left[key] and b in right[key]
        assert exists==(match[a,b] in expected[key])
        if exists:
            nxt=left[key][a],right[key][b]; target=expected[key][match[a,b]]
            if nxt in match: assert match[nxt]==target
            else: assert target not in match.values(); match[nxt]=target; queue.append(nxt)
assert len(match)==10

print(json.dumps(dict(verification='PASS',central_identity_factors=19,
  central_signed_counts=[0,1,0,0,0],central_lift_matrix=True,
  crossing_and_boundary_identities=True,trefoil_dictionary=True,
  corrected_mixed_sign_conjugator='b^-1 a^-1 c u',modular_action_states=30,
  whitehead_differences=96,pair_exterior_algebra=True,
  trefoil_intersection_product_vertices=10,
  geometric_arguments_verified_by_program=False,full_group_decided=False),indent=2))
