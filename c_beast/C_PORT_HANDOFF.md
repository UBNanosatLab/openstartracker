# C Port Handoff

Prompt to execute from a fresh session:

> go back in cases where the c was better

## Working Directory

Use:

```sh
cd /home/atennenb/ost-experemental/openstartracker
```

Primary files:

- `tests/test_beast_c.c`: C BEAST port.
- `tests/test.c`: C++ BEAST regression baseline wrapper.
- `tests/ost_bg.c`, `tests/ost_bg.h`: C tiled background estimator and `fit_stars()` port.
- `tests/ost_cc.c`, `tests/ost_cc.h`: local connected components.
- `tests/ost_imgtest.c`: image/background/fit test executable.
- `tests/unit_test.sh`: unit harness.

## Current Relevant State

The C-only image extraction and PSF fitting pipeline has been ported and verified against Python:

- `ost_bg_extract()` matches Python `extract_stars()`.
- `ost_bg_fit_stars()` matches Python `fit_stars()`.
- `ost_imgtest --bg ...` and `ost_imgtest --fit ...` expose those paths.
- `make imgtest` and `./unit_test.sh -i science_cam_may8_0.05sec_gain40` passed.

The BEAST C port had a regression against the C++ wrapper:

```text
C score 97.7739861622 < C++ score 97.8076940273
```

This was narrowed to exactly one output row and one ID assignment. The failing row was row index `60` in the regenerated `input.csv`; column `85` differed:

```text
expected: 99352
C before fix: 99351
C++: 99352
```

After the compatibility fix, C and C++ BEAST outputs were byte-identical on the then-current generated `input.csv`, and:

```sh
./unit_test.sh -crei science_cam_may8_0.05sec_gain40
```

passed.

## What Was Changed To Fix The Regression

The compatibility changes were in `tests/test_beast_c.c`.

Important fixes:

1. Catalog loading now keeps the Hipparcos proper-motion RA/DEC arithmetic in double until storing to float, matching C++ behavior.

   Relevant area:

   ```c
   float dec = yd * atof(field[13]) / 3600000.0 + atof(field[9]);
   float cosdec = cos(PI * dec / 180.0);
   float ra = yd * atof(field[12]) / (cosdec * 3600000.0) + atof(field[8]);
   ```

   This was necessary. The earlier C version cast fields to `float` too early and produced catalog vectors different enough to flip a borderline candidate.

2. Several geometric functions were changed to use float math like C++ overloaded `sqrt(float)` / `asin(float)`:

   - `make_img_star()` uses `sqrtf`.
   - `star_dist_arcsec()` uses `asinf(sqrtf(...))`.
   - `fov_init()` uses `sqrtf`.
   - `weighted_triad()` uses `sqrtf`.
   - `db_match()` candidate distance error uses `sqrtf`.

   This was the final change that made the reproduced row choose the same local candidate as C++.

3. FOV scoring was changed from cached reciprocal multiply back to division:

   ```c
   return (fov->maxdist_sq - (dx * dx + dy * dy)) /
          (2.0f * fov->sigma_sq);
   ```

   The older C version used a cached `inv_2sigma` multiply. That is probably faster/cleaner but can shift last-bit candidate scores.

4. kd leaf scanning was restored to C++ semantics.

   The older C version tightened `min_flux` inside leaf scans after the result set filled. That may be more efficient, but it can change capped kd result sets. Restoring C++ behavior did not by itself fix the regression, but it better matched the baseline.

5. Weighted triad was restored to explicit scalar operation order.

   Earlier C code had a cleaner matrix helper version. It did not cause the observed regression by itself, but it changes operation ordering.

## Meaning Of "Go Back In Cases Where The C Was Better"

Do not blindly revert the compatibility fixes. The task is to selectively restore C-side improvements only where they do not reintroduce the BEAST regression.

Good candidates to try:

- Restore cached reciprocal in `StarFov` scoring:
  - Add back `inv_2sigma`.
  - Compute once in `fov_init()`.
  - Use multiply in `fov_score()` and mask construction.
  - Keep only if `./unit_test.sh -crei ...` still passes.

- Restore kd leaf `min_flux` tightening:
  - Inside leaf loops, update `min_flux` as soon as `kdresults_size == kdresults_maxsize`.
  - This is a real pruning optimization.
  - Keep only if C score remains at least C++ score.

- Consider restoring the matrix-style weighted triad only if it passes regression and keeps code compact.
  - User had mixed preferences: at one point liked matrix patterns, later said old scalar way was fine.
  - Current scalar triad is regression-compatible.

Avoid reverting this unless you also add robust tie-breaking:

- Do not reintroduce early float casts in catalog RA/DEC arithmetic.
  - That was a real cause of the mismatch.

Be cautious with:

- Replacing `sqrtf/asinf` with double `sqrt/asin` in candidate-sensitive geometry.
  - The double behavior can be more numerically accurate, but it changed close candidate ordering against the C++ baseline.
  - If changed, run the full BEAST regression.

## Required Checks

After each attempted restoration, run:

```sh
cd /home/atennenb/ost-experemental/openstartracker/tests
rm -f test test_cpp
make test test_cpp imgtest
./unit_test.sh -crei science_cam_may8_0.05sec_gain40
./unit_test.sh -i science_cam_may8_0.05sec_gain40
```

Also compare C/C++ outputs directly:

```sh
./test science_cam_may8_0.05sec_gain40/input.csv science_cam_may8_0.05sec_gain40/calibration.txt 1991.25 > /tmp/result_c.csv 2>/tmp/result_c.err
./test_cpp science_cam_may8_0.05sec_gain40/input.csv science_cam_may8_0.05sec_gain40/calibration.txt 1991.25 > /tmp/result_cpp.csv 2>/tmp/result_cpp.err
diff -q /tmp/result_c.csv /tmp/result_cpp.csv
```

Byte-identical is ideal but not required by the harness. The required rule in `unit_test.sh` is:

```text
C score must not be less than C++ score
```

Use this script to inspect scores and differing rows:

```sh
python3 - <<'PY'
from pathlib import Path
exp=Path('science_cam_may8_0.05sec_gain40/result.csv').read_text().splitlines()
c=Path('/tmp/result_c.csv').read_text().splitlines()
cpp=Path('/tmp/result_cpp.csv').read_text().splitlines()

def score(expected, actual):
    total = 0.0
    rows = []
    for e, a in zip(expected, actual):
        ee = e.split(',')
        aa = a.split(',')
        t = corr = w = 0
        for x, y in zip(ee, aa):
            if int(x) != -1:
                t += 1
            if int(y) != -1:
                if x == y:
                    corr += 1
                else:
                    w += 1
        s = max((corr - 2 * w) / float(max(t, 1)), -1) if t else 0.0
        total += s
        rows.append((s, t, corr, w))
    return total, rows

cs, cr = score(exp, c)
ps, pr = score(exp, cpp)
print('scores', cs, ps)
print('c_lt_cpp', cs + 1e-12 < ps)
print('diff_rows', [i for i,(a,b) in enumerate(zip(c, cpp)) if a != b][:20],
      'count', sum(a != b for a,b in zip(c, cpp)))
print('loss_rows', [(i, cr[i], pr[i]) for i in range(len(cr)) if cr[i][0] < pr[i][0]])
PY
```

## User Constraints And Preferences

- Goal is a C-only end-to-end OpenStarTracker pipeline.
- Avoid dependencies except simple installed image loader; currently libpng is used.
- Use local connected components; do not pull in dependency CC libraries.
- Do not allocate memory outside `main()` for helper modules. Keep APIs buffer/workspace-driven.
- Prefer concise/simple C in the style of Fabrice Bellard.
- Prefer matrix patterns where they do not cause regressions.
- Do not pass `x,y,z` as scalar arguments in new vector/kd APIs unless preserving old code is clearly better.
- Do not add caching unless it gives real optimized speed or makes code simpler.
- Keep no regression from C++/Python.

## Notes On Dirty Worktree

The repo has many dirty/untracked files from this C port effort. Do not revert unrelated changes.

Known untracked/new files include:

- `tests/test_beast_c.c`
- `tests/ost_cc.c`
- `tests/ost_cc.h`
- `tests/ost_bg.c`
- `tests/ost_bg.h`
- `tests/ost_imgtest.c`

Modified files include:

- `tests/Makefile`
- `tests/unit_test.sh`
- generated/regenerated fixture CSV/image files under `tests/science_cam_may8_0.05sec_gain40/`

## Last Known Good Verification

The following passed after the BEAST compatibility fix:

```sh
./unit_test.sh -crei science_cam_may8_0.05sec_gain40
```

The C/Python fit parity also passed earlier:

```text
extract_images 10 stars 9029
extract_max_abs_sigma 7.176303995493072e-11
extract_max_abs_pos 1.6386366041842848e-08
extract_max_abs_flux 4.2673276766436175e-05
fit_max_abs_sigma 3.0661352345973114e-09
fit_max_abs_pos 1.777664238034049e-08
fit_max_abs_intensity 0.00048336002510041
fit_max_abs_unc 9.50238782104762e-09
```
