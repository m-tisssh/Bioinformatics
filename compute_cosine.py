import argparse, csv
from collections import Counter

# чтение предобработанного TSV 
def read_minimal(path):
    rows = []
    with open(path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for r in reader:
            rows.append(r)
    return rows

# косинусное сходство между строками 
def cosine(a,b):
    va = Counter(a); vb = Counter(b)
    keys = set(va) | set(vb)
    num = sum(va[k]*vb[k] for k in keys)
    norm_a = sum(v*v for v in va.values()) ** 0.5
    norm_b = sum(v*v for v in vb.values()) ** 0.5
    if norm_a==0 or norm_b==0:
        return 0.0
    return num / (norm_a * norm_b)

# приведение % идентичности
def to_float_percent(pstr):
    if pstr is None or pstr == '':
        return 0.0
    try:
        v = float(pstr)
    except:
        return 0.0
    if v > 1.0 + 1e-9:
        v = v / 100.0
    return v

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--preprocessed', required=True)
    p.add_argument('--top', type=int, default=10)
    p.add_argument('--w-match', type=float, default=0.6)
    p.add_argument('--w-cos', type=float, default=0.4)
    p.add_argument('--out', default=None, help='optional CSV output')
    args = p.parse_args()

    # нормализация весов (сумма = 1)
    s = args.w_match + args.w_cos
    if s <= 0:
        raise SystemExit("non-positive weight sum")
    args.w_match /= s; args.w_cos /= s

    rows = read_minimal(args.preprocessed)
    if not rows:
        raise SystemExit("no rows read")

    # находим таргет (отмечен offset == TARGET)
    target = None
    for r in rows:
        if r.get('offset','') == 'TARGET':
            target = r; break
    if target is None:
        target = rows[0]
    target_seq = (target.get('seq_clean') or '').upper()

    results = []
    for r in rows:
        h = r.get('header','')
        if h == target.get('header'): continue
        seq = (r.get('seq_clean') or '').upper()
        perc = to_float_percent(r.get('percent_identity',''))
        cos = cosine(target_seq, seq)
        combined = args.w_match * perc + args.w_cos * cos
        results.append({
            'header': h,
            'len': len(seq),
            'percent_identity': round(perc,6),
            'cosine': round(cos,6),
            'combined': round(combined,6)
        })

    results.sort(key=lambda x: x['combined'], reverse=True)

    # печать топ-N
    print(f"Target: {target.get('header','<unknown>')} (len={len(target_seq)})")
    print(f"{'rank':>4} {'combined':>8} {'%id':>8} {'cosine':>8}  header")
    for i, r in enumerate(results[:args.top], 1):
        print(f"{i:4d} {r['combined']:8.6f} {r['percent_identity']:8.6f} {r['cosine']:8.6f}  {r['header']}")

    # опционально сохраняем в CSV
    if args.out:
        with open(args.out, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=['rank','header','combined','percent_identity','cosine','len'])
            writer.writeheader()
            for i, r in enumerate(results, 1):
                out = r.copy(); out['rank']=i
                writer.writerow(out)
        print("wrote:", args.out)

if __name__ == '__main__':
    main()

# шаг 1: предобработка
# python preprocess.py --fasta uniprot_sprot.fasta --target-header "sp|Q6GZX4|001R_FRG3G" --out test_4.preproc.tsv
    
# шаг 2: расчёт метрик
# python compute_cosine.py --preprocessed test_4.preproc.tsv --top 10 --w-match 0.6 --w-cos 0.4 --out results.csv
