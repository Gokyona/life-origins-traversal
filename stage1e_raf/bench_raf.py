"""Stage-1e benchmark: the stage-1b harness (../stage1b/bench1b.py) run on the RAF model.

Reuses bench1b unchanged: its methods, map-maker and scoring look the model up as
the module-level name `model`, so swapping in raf.py is enough. All bench1b
options work, e.g.

  python bench_raf.py --device xpu --surrogate es --fast --budget 64000 --reps 5 \
      --n-test 4000 --gt-seeds 64 --methods random_1seed,sobol_1seed,engine3,engine4_soft --out results_raf.json

One RAF instance costs ~1 ms on the Arc B70 (vs ~1 us for model10), so ground
truth and budgets are smaller than in stage 1b.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
STAGE1B = os.environ.get("LOT_STAGE1B", os.path.join(HERE, "..", "stage1b"))
sys.path.insert(0, HERE)
sys.path.insert(1, STAGE1B)

import bench1b  # noqa: E402
import raf  # noqa: E402

bench1b.model = raf

if __name__ == "__main__":
    bench1b.main()
