# скольжение + косинусовое (пример)
from collections import Counter
import math

# скольжение и оптимальный выбор
def sliding_details(target, candidate):
    n,m = len(target), len(candidate)
    rows=[]
    for offset in range(-m+1, n):
        i0 = max(0, offset); i1 = min(n, offset+m)
        overlap = i1-i0
        if overlap<=0:
            rows.append((offset,"","",overlap,0,0.0)); continue
        j0 = i0-offset
        tseg = target[i0:i1]; sseg = candidate[j0:j0+overlap]
        matches = sum(1 for a,b in zip(tseg,sseg) if a==b)
        rows.append((offset,tseg,sseg,overlap,matches,matches/overlap))
    return rows

def best_sliding(target, candidate, min_overlap=1):
    rows = sliding_details(target, candidate)
    valid = [r for r in rows if r[3]>=min_overlap]
    if not valid: return None
    # выбор про %, затем по совпадениям, затем по перекрытию
    valid.sort(key=lambda x: (x[5], x[4], x[3]), reverse=True)
    return valid[0]

def cosine(a,b):
    v1 = Counter(a); v2 = Counter(b)
    keys = sorted(set(v1)|set(v2))
    vec1 = [v1[k] for k in keys]; vec2 = [v2[k] for k in keys]
    dot = sum(x*y for x,y in zip(vec1,vec2))
    na = math.sqrt(sum(x*x for x in vec1)); nb = math.sqrt(sum(y*y for y in vec2))
    return dot/(na*nb) if na>0 and nb>0 else 0.0

# Пример
t="ABCDEF"; s="BCDA"
print("Sliding table:")
for r in sliding_details(t,s):
    print(r)
best=best_sliding(t,s,min_overlap=1)
print("Best:", best)

A='ababcbbc'; B='bc'
print("Cosine:", cosine(A,B))

# комбинированный пример
w1,w2 = 0.6,0.4
best_percent = best[5] if best else 0.0
combined = w1*best_percent + w2*cosine(A,B)
print("Combined:", combined)
