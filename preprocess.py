import argparse, csv, re
from pathlib import Path
from collections import Counter
import numpy as np

# чтение .fastа  
def parse_fasta(path):
    recs = []
    h = None
    buf = []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line: continue
            if line.startswith('>'): # новая запись
                if h is not None:
                    recs.append((h, ''.join(buf)))
                h = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip()) # добавляем последовательность
        if h is not None:
            recs.append((h, ''.join(buf)))
    return recs

# очистка последовательности (только буквы) 
def clean(s):
    return re.sub(r'[^A-Za-z]', '', s).upper()

# выбор длины 
def next_len(n):
    try:
        return np.fft.next_fast_len(n)
    except AttributeError:
        p = 1
        while p < n:
            p <<= 1
        return p

# построение бинарных массивов для каждой буквы 
def build_target_arrays(target_seq, alphabet):
    return {ch: np.fromiter((1.0 if c==ch else 0.0 for c in target_seq), dtype=np.float32) for ch in alphabet}

# Метод скольжения 
def best_by_fft(target_arrays, target_seq, cand_seq, min_overlap=1, fft_cache=None):
    n = len(target_seq); m = len(cand_seq)
    if n==0 or m==0:
        return {"matches":0,"overlap":0,"percent":0.0,"offset":None}

    conv_len = n + m - 1
    nfft = next_len(conv_len)

    # для ускорения
    if fft_cache is None:
        fft_cache = {}
    if nfft not in fft_cache:
        fft_cache[nfft] = {ch: np.fft.rfft(np.pad(arr, (0, nfft - len(arr)))) for ch, arr in target_arrays.items()}
    target_ffts = fft_cache[nfft]

    # суммируем свёртки для совпадающих букв
    sum_freq = None
    cand_counts = Counter(cand_seq)
    for ch in cand_counts:
        if ch not in target_ffts:
            continue
        b = np.zeros(nfft, dtype=np.float32)
        for i, c in enumerate(cand_seq):
            if c == ch:
                b[m-1-i] = 1.0
        B = np.fft.rfft(b)
        prod = target_ffts[ch] * B
        sum_freq = prod if sum_freq is None else sum_freq + prod

    if sum_freq is None:
        return {"matches":0,"overlap":0,"percent":0.0,"offset":None}

    # восстанавливаем результаты 
    conv = np.fft.irfft(sum_freq, nfft)[:conv_len]
    offsets = np.arange(-m+1, n)
    idx = offsets + (m - 1)
    matches_vec = np.rint(conv[idx]).astype(int)
    overlap_vec = np.minimum(n, offsets + m) - np.maximum(0, offsets)

    valid = overlap_vec >= min_overlap
    if not valid.any():
        return {"matches":0,"overlap":0,"percent":0.0,"offset":None}
    percent = np.full_like(matches_vec, -1.0, dtype=float)
    percent[valid] = matches_vec[valid] / overlap_vec[valid]

    # выбираем лучшее совпадение
    best_idx = int(np.argmax(percent))
    cand_idx = np.where(percent == percent[best_idx])[0]
    if cand_idx.size > 1:
        mm = matches_vec[cand_idx]
        j = int(cand_idx[np.argmax(mm)])
        best_idx = j

    return {
        "matches": int(matches_vec[best_idx]),
        "overlap": int(overlap_vec[best_idx]),
        "percent": float(matches_vec[best_idx] / overlap_vec[best_idx]) if overlap_vec[best_idx]>0 else 0.0,
        "offset": int(offsets[best_idx])
    }

def main():
    # аргументы запуска
    p = argparse.ArgumentParser()
    p.add_argument('--fasta', required=True)
    p.add_argument('--out', default='preprocessed_minimal.tsv')
    p.add_argument('--target-header', required=True, help="Header of the sequence to use as target (exact or substring match)")
    p.add_argument('--min-overlap', type=int, default=3)
    args = p.parse_args()

    recs = parse_fasta(args.fasta)
    if not recs:
        raise SystemExit("no records")

    # выбираем таргет только по заголовку
    found = [r for r in recs if r[0] == args.target_header or args.target_header in r[0]]
    if not found:
        raise SystemExit("target-header not found")
    target_header, target_raw = found[0]

    target_seq = clean(target_raw)
    if not target_seq:
        raise SystemExit("empty target after cleaning")

    # алфавит = все символы
    observed = set(target_seq)
    for _, raw in recs:
        observed.update(clean(raw))
    alphabet = sorted(observed)

    target_arrays = build_target_arrays(target_seq, alphabet)
    fft_cache = {}

    # сохраняем результат
    out_path = Path(args.out)
    with out_path.open('w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['header','seq_clean','percent_identity','matches','overlap','offset'])
        w.writerow([target_header, target_seq, '', '', '', 'TARGET'])
        for h, raw in recs:
            if h == target_header: continue
            s = clean(raw)
            if not s:
                w.writerow([h, '', 0.0, 0, 0, ''])
                continue
            best = best_by_fft(target_arrays, target_seq, s, min_overlap=args.min_overlap, fft_cache=fft_cache)
            w.writerow([h, s, f"{best['percent']:.6f}", best['matches'], best['overlap'], best['offset'] if best['offset'] is not None else ''])
    print("done. wrote:", out_path)

if __name__ == '__main__':
    main()

# запуск:
# python preprocess.py --fasta uniprot_sprot.fasta --target-header "sp|Q6GZX4|001R_FRG3G" --out test.preproc.tsv
