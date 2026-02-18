import random
from pathlib import Path

random.seed(42)
base_dir = Path('/Users/aotoole/repositories/raccoon/examples/large_demo')
inputs_dir = base_dir / 'seq_qc_inputs'
inputs_dir.mkdir(parents=True, exist_ok=True)

n_seqs = 30
length = 1000
bases = ['A', 'C', 'G', 'T']

ref = ''.join(random.choice(bases) for _ in range(length))

seqs = {}
issues = {}
for i in range(1, n_seqs + 1):
    sid = f"KBX{i:03d}"
    seq = list(ref)
    issues[sid] = []

    if i % 7 == 0:
        start = random.randint(200, 240)
        for j in range(start, start + 8):
            seq[j] = random.choice([b for b in bases if b != seq[j]])
        issues[sid].append('clustered_snps')

    if i % 6 == 0:
        for j in range(400, 500):
            seq[j] = 'N'
        issues[sid].append('high_N')

    if i % 10 == 0:
        seq = seq[:700]
        issues[sid].append('short')

    if i % 5 == 0:
        for j in random.sample(range(50, 950), 5):
            if j < len(seq):
                seq[j] = random.choice([b for b in bases if b != seq[j]])
        issues[sid].append('scattered_snps')

    if i % 4 == 0:
        for j in range(800, 820):
            if j < len(seq):
                seq[j] = 'N'
        issues[sid].append('N_block')

    seqs[sid] = ''.join(seq)

alignment_path = base_dir / 'alignment.fasta'
with alignment_path.open('w') as handle:
    for sid, seq in seqs.items():
        if len(seq) < length:
            seq = seq + ('N' * (length - len(seq)))
        handle.write(f">{sid}\n{seq}\n")

set_a = list(seqs.items())[:10]
set_b = list(seqs.items())[10:20]
set_c = list(seqs.items())[20:]
for name, subset in [('set_a.fasta', set_a), ('set_b.fasta', set_b), ('set_c.fasta', set_c)]:
    with (inputs_dir / name).open('w') as handle:
        for sid, seq in subset:
            handle.write(f">{sid}\n{seq}\n")

meta_path = base_dir / 'metadata.csv'
locations = ['SiteA', 'SiteB', 'SiteC', 'SiteD', 'SiteE']
with meta_path.open('w') as handle:
    handle.write('id,location,date,host\n')
    for i, sid in enumerate(seqs.keys(), start=1):
        loc = locations[i % len(locations)]
        date = f"2024-{(i % 12) + 1:02d}-{(i % 27) + 1:02d}"
        host = 'human' if i % 3 else 'animal'
        handle.write(f"{sid},{loc},{date},{host}\n")

issues_path = base_dir / 'issues_summary.csv'
with issues_path.open('w') as handle:
    handle.write('id,issues\n')
    for sid, iss in issues.items():
        handle.write(f"{sid},{';'.join(iss)}\n")

print(f"Wrote large demo to {base_dir}")
