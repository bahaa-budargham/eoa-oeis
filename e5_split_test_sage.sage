# EOA English E5 13-13 split null + small-core null (SageMath version)
# Replaces e5_split_test.py / verify_coreclaim.py numpy approximations with
# exact Sage factorization over QQ. Author: Bahaa BuDargham (adapted by EOA).
# Usage:  sage e5_split_test_sage.sage  [n_alphabets]
import random, sys
from collections import Counter
set_random_seed(42)
R.<x> = PolynomialRing(QQ)

N = int(sys.argv[1]) if len(sys.argv) > 1 else 250
equal_split_count = 0
core_size_counts = Counter()
lin_hist = Counter()

for i in range(N):
    if (i + 1) % 50 == 0:
        print("  ... alphabet %d / %d" % (i + 1, N)); sys.stdout.flush()
    M = matrix(ZZ, 26, 26, 0)
    for row in range(26):
        for _ in range(random.randint(2, 5)):
            M[row, random.randint(0, 25)] += 1
    sccs = DiGraph(M).strongly_connected_components()
    core_size_counts[max(len(c) for c in sccs)] += 1
    degs = [f.degree() for f, m in M.charpoly().factor() for _ in range(m)]
    lin_hist[degs.count(1)] += 1
    if degs.count(1) == 13 and degs.count(13) == 1:
        equal_split_count += 1

print("\nDominant core sizes (26x26, name lengths 2-5):")
for sz, c in sorted(core_size_counts.items()):
    print("  size %2d: %3d (%.1f%%)" % (sz, c, 100.0 * c / N))
print("\n13-13 equal splits: %d of %d" % (equal_split_count, N))
print("linear-factor histogram:", dict(sorted(lin_hist.items())))
