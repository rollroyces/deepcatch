#!/usr/bin/env python3
"""Full pipeline comparison: Classical MDS vs Neural Model"""
import torch, numpy as np, json, time
from model import *
from model.tokenizer import DNATokenizer, DNA_BASES, SPECIAL_TOKENS
from collections import Counter

torch.manual_seed(42); np.random.seed(42)

# ── STEP 1: Train tokenizer with max merge length constraint ──
corpus_seqs = [
    'ACGT'*400 + 'CGCG'*200 + 'ATAT'*100 + 'GCGC'*300 + 'CCGG'*100,
    'TTTT'*300 + 'AAAA'*300 + 'CGCG'*150 + 'ACGT'*50 + 'GGCC'*200,
    'GCGC'*450 + 'CCGG'*200 + 'ATAT'*100 + 'NNNN'*50,
    'ACGT'*600 + 'GGCC'*150 + 'TTTT'*100 + 'AAAA'*50 + 'CG'*300,
    'CCGG'*300 + 'CGCG'*300 + 'ATAT'*200 + 'ACGT'*100 + 'G'*500,
]

VOCAB_SIZE = 128
MAX_MERGE_LEN = 8

vocab = {}
for i, tok in enumerate(SPECIAL_TOKENS + DNA_BASES):
    vocab[tok] = i

token_lists = [[c for c in s.upper() if c in set(DNA_BASES)] for s in corpus_seqs if s]
merges = []
for step in range(VOCAB_SIZE - len(vocab)):
    pair_counts = Counter()
    for tokens in token_lists:
        for i in range(len(tokens) - 1):
            pair_counts[(tokens[i], tokens[i + 1])] += 1
    if not pair_counts:
        break
    best_pair, best_count = pair_counts.most_common(1)[0]
    if best_count < 2:
        break
    merged_token = best_pair[0] + best_pair[1]
    if len(merged_token) > MAX_MERGE_LEN:
        break
    merges.append(best_pair)
    vocab[merged_token] = len(vocab)
    for idx, tokens in enumerate(token_lists):
        merged = []
        i = 0
        while i < len(tokens):
            if (i+1 < len(tokens) and tokens[i] == best_pair[0]
                    and tokens[i+1] == best_pair[1]):
                merged.append(merged_token)
                i += 2
            else:
                merged.append(tokens[i])
                i += 1
        token_lists[idx] = merged

tok = DNATokenizer(vocab=vocab, merges=merges)
tok.save('model/dna_bpe_tokenizer.json')
print(f'Tokenizer: vocab_size={tok.vocab_size}, merges={len(merges)}')
top_tokens = sorted(vocab.items(), key=lambda x: x[1])[-10:]
for t, i in top_tokens:
    display = t if len(t) <= 20 else t[:17] + "..."
    print(f'  token[{i:3d}] = "{display}" (len={len(t)})')

# ── STEP 2: Generate training data ──
cancer_signatures = ['CG', 'GC', 'CGC', 'GCG', 'CCG', 'CGG', 'GGC', 'GCC']
healthy_signatures = ['ACGT', 'AC', 'GT', 'TG', 'CA', 'AG', 'CT', 'GA']

def gen_seq(signatures, noise=200):
    parts = [s * np.random.randint(15, 60) for s in signatures]
    parts.append(''.join(np.random.choice(['A','C','G','T'], noise)))
    np.random.shuffle(parts)
    return ''.join(parts)

def gen_cancer():
    return gen_seq(cancer_signatures)

def gen_healthy():
    return gen_seq(healthy_signatures)

n_train, n_test = 500, 200
X_tr, y_tr = [], []
for _ in range(n_train):
    is_cancer = np.random.random() < 0.5
    seq = gen_cancer() if is_cancer else gen_healthy()
    X_tr.append(tok.count_frequencies(seq))
    y_tr.append(1.0 if is_cancer else 0.0)

X_te, y_te = [], []
for _ in range(n_test):
    is_cancer = np.random.random() < 0.5
    seq = gen_cancer() if is_cancer else gen_healthy()
    X_te.append(tok.count_frequencies(seq))
    y_te.append(1.0 if is_cancer else 0.0)

X_tr = np.array(X_tr); y_tr = np.array(y_tr)
X_te = np.array(X_te); y_te = np.array(y_te)
print(f'\nTrain: {len(X_tr)} | Test: {len(X_te)}')
print(f'Train balance: cancer={y_tr.sum():.0f} healthy={(1-y_tr).sum():.0f}')

# ── STEP 3: Train ──
model = MotifDiversityModel(vocab_size=tok.vocab_size, d_model=64,
                            n_heads=4, n_layers=2)
opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-5)
crit = torch.nn.BCELoss()

X_t = torch.tensor(X_tr, dtype=torch.float32)
y_t = torch.tensor(y_tr, dtype=torch.float32).unsqueeze(1)
X_v = torch.tensor(X_te, dtype=torch.float32)
y_v = torch.tensor(y_te, dtype=torch.float32).unsqueeze(1)

best_acc = 0
for epoch in range(60):
    model.train()
    perm = torch.randperm(len(X_t))
    for i in range(0, len(X_t), 32):
        idx = perm[i:i+32]
        opt.zero_grad()
        loss = crit(model(X_t[idx]), y_t[idx])
        loss.backward()
        opt.step()
    model.eval()
    with torch.no_grad():
        vp = model(X_v)
        vacc = ((vp > 0.5).float() == y_v).float().mean().item()
    if vacc > best_acc:
        best_acc = vacc
        torch.save({
            'model_state_dict': model.state_dict(),
            'config': {'vocab_size': tok.vocab_size, 'd_model': 64,
                       'n_heads': 4, 'n_layers': 2},
        }, 'model/motif_model_checkpoint.pt')
    if epoch % 15 == 0:
        print(f'  epoch {epoch:3d}: val_acc={vacc:.4f}')
print(f'Best val acc: {best_acc:.4f}')

# ── STEP 4: COMPARISON TABLE ──
print()
print('='*78)
print('FINAL COMPARISON: Classical MDS vs Neural Model')
print('='*78)
header = f'{"Sample":<20s} {"MDS(Simp)":>10s} {"MDS(Shan)":>10s} {"Neural":>8s} {"Pred":>8s} {"Truth":>8s}'
print(header)
print('-'*78)

predictor = MotifPredictor(model, device='cpu')

# Collect predictions
all_classical_preds = []
all_neural_preds = []

for i in range(min(20, len(X_te))):
    freqs = X_te[i]
    est_tokens = int(np.sum(freqs > 0) * 100)
    # Reconstruct approx counts for MDS
    total_tok = max(1, int(np.sum(freqs * 10000)))
    counts = (freqs * total_tok).astype(np.int64)
    mds_s = compute_classical_mds(counts, 'simpson')
    mds_h = compute_classical_mds(counts, 'shannon')
    r = predictor.predict(freqs)
    
    all_classical_preds.append(1.0 if mds_s > 0.7 else 0.0)
    all_neural_preds.append(1.0 if r['probability'] >= 0.5 else 0.0)
    
    label = f'Sample {i}'
    cancer_str = "cancer" if r['probability'] >= 0.5 else "healthy"
    truth_str = "cancer" if y_te[i] else "healthy"
    print(f'{label:<20s} {mds_s:10.4f} {mds_h:10.4f} {r["probability"]:8.4f} {cancer_str:>8s} {truth_str:>8s}')

print('-'*78)
cp_arr = np.array(all_classical_preds)
np_arr = np.array(all_neural_preds)
yt_arr = y_te[:20]
print(f'{"Classical MDS accuracy":<20s} {(cp_arr == yt_arr).mean():>10.4f}')
print(f'{"Neural model accuracy":<20s} {(np_arr == yt_arr).mean():>10.4f}')

# Full test set
full_classical = []
for i in range(len(X_te)):
    total_tok = max(1, int(np.sum(X_te[i] * 10000)))
    counts = (X_te[i] * total_tok).astype(np.int64)
    mds_s = compute_classical_mds(counts, 'simpson')
    full_classical.append(1.0 if mds_s > 0.7 else 0.0)
full_classical = np.array(full_classical)

model.eval()
with torch.no_grad():
    full_neural = (model(X_v).numpy().flatten() > 0.5).astype(float)

print(f'{"Full Classical MDS":<20s} {(full_classical == y_te).mean():>10.4f}')
print(f'{"Full Neural Model":<20s} {(full_neural == y_te).mean():>10.4f}')

# ── STEP 5: Explain one prediction ──
print()
print('='*78)
print('INTEGRATED GRADIENTS — What drives a cancer prediction?')
print('='*78)
c_seq = gen_cancer()
c_freqs = tok.count_frequencies(c_seq)
expl = explain_prediction(c_freqs, model, tokenizer=tok, top_k=8)
r = predictor.predict(c_freqs)
print(f'Prediction: prob={r["probability"]:.4f} -> {r["prediction"]}')
print(f'Classical MDS (baseline): {compute_classical_mds((c_freqs*len(tok.encode(c_seq))).astype(np.int64), "simpson"):.4f}')
print(f'Convergence delta: {expl.get("convergence_delta", "N/A")}')
print(f'Top motifs (increase cancer probability):')
for m in expl['top_motifs'][:5]:
    print(f'  "{m["token"]}"  ->  {m["attribution"]:+.6f}')
print(f'Bottom motifs (decrease cancer probability):')
for m in expl['bottom_motifs'][:5]:
    print(f'  "{m["token"]}"  ->  {m["attribution"]:+.6f}')

print()
print(f'Model saved: model/motif_model_checkpoint.pt')
print(f'Tokenizer saved: model/dna_bpe_tokenizer.json')
