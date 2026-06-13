# c_beast

`c_beast.h` is a single-header C implementation of the current OpenStarTracker C path:

```text
PNG loader supplied by caller
RGBA -> gray16 -> tiled histogram background -> direct thresholded CC runs
    -> weighted connected components -> PSF fit -> BEAST star ID
```

Use it by including the header normally for declarations:

```c
#include "c_beast/c_beast.h"
```

In exactly one C file, define `C_BEAST_IMPLEMENTATION` before including it:

```c
#define C_BEAST_IMPLEMENTATION
#include "c_beast/c_beast.h"
```

The library does not allocate inside the image-processing or tracker helpers.
Callers allocate the tracker, connected-components workspace, background maps,
fit workspace, and result buffers, then pass them in.

The current test consumers are:

- `tests/test_beast_c.c`: CSV unit-test BEAST runner.
- `tests/ost_imgtest.c`: image extraction and fit test runner.
- `tests/ost_pipeline.c`: full C image-to-IDs pipeline.

Run the regression suite from `tests/`:

```sh
make all
./unit_test.sh -ei science_cam_may8_0.05sec_gain40
```
