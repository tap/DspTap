# DspTap

Shared DSP primitives for the **Tap** family of audio libraries. Header-only,
plain portable C++ (C++20, standard library only), no Max/Min or framework
dependency — consumed as a git submodule by the individual libraries.

Today it holds one primitive:

## `dsptap::real_fft` — real FFT with a fixed numeric contract

`include/dsptap/fft.h` wraps the vendored [Ooura split-radix real
FFT](third_party/ooura/) behind a small, well-specified interface, with
**optional, mutually-exclusive float32 backends** that re-present the *exact*
same numeric contract for speed on specific hardware:

| Backend | Build option | Target | Notes |
|---------|-------------|--------|-------|
| Ooura (default) | — | everywhere | the golden model; double is **always** Ooura |
| CMSIS-DSP Helium | `DSPTAP_FFT_CMSIS` | bare-metal Cortex-M55 (MVE) | ~3× fewer instructions/transform |
| Apple vDSP | `DSPTAP_FFT_ACCELERATE` | macOS / Apple Silicon | ~3× faster/transform |

The two float32 backends conjugate imaginary bins and rescale so every
intermediate spectrum matches the Ooura build to single-precision rounding —
so the whole double-precision test battery stays a valid oracle for the
accelerated float paths. `tests/test_fft_backend.cpp` pins each backend to
Ooura's `rdft_f` bin-for-bin at the certified geometries (N = 512, 2048).

```cpp
#include "dsptap/fft.h"

dsptap::real_fft   fft(1024);   // double, the desktop/golden profile
dsptap::real_fft32 fft32(1024); // float,  the embedded / accelerated profile

std::vector<double> x(1024, 0.0);
fft.forward_inplace(x.data());        // packed spectrum, W = exp(+2πi/N)
fft.inverse(x.data(), x.data());      // out-of-place inverse, normalized (2/N)
```

Key contract points (full detail in the header docstring):

- **Packing** (N/2 + 1 bins): `data[0]` = DC real, `data[1]` = Nyquist real,
  `data[2k]`/`data[2k+1]` = bin *k* real/imag for `1 ≤ k < N/2`.
- **Sign convention** `W = exp(+2πi/N)` — imaginary parts are *conjugated*
  relative to the engineering-convention DFT. Consistent across every operand,
  so spectral products (fast convolution, adaptive-filter regressors) are
  unaffected; conjugate only when importing spectra computed elsewhere.
- **Normalization**: `*_inplace` inverse is unnormalized (multiply by `2/N`);
  the out-of-place `inverse()` applies the `2/N` for you.
- Transforms are `noexcept` and allocation-free after construction — real-time
  safe. Size must be a power of two, `≥ 4`, fixed at construction.

## Build

Standalone (builds the Ooura path, plus vDSP on macOS, and runs the tests):

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

The CMSIS-DSP Helium backend is compile-verified under the bare-metal M55
toolchain (`cmake/arm-cortex-m55-mps3.cmake`); running its parity suite under
QEMU is done in the consuming library's embedded harness.

### As a submodule

```cmake
add_subdirectory(submodules/dsptap)   # or however it is pinned
target_link_libraries(my_dsp PRIVATE DspTap::DspTap)
```

`DspTap::DspTap` is an INTERFACE target (the headers + the compiled Ooura
static lib `DspTap::fft`); it does not build the tests when added as a
subdirectory (`DSPTAP_BUILD_TESTS` defaults OFF unless top-level). The
per-platform float32 backend defaults follow the target: vDSP on Apple, CMSIS
on the bare-metal M55 profile, Ooura elsewhere — override with
`-DDSPTAP_FFT_ACCELERATE=OFF` etc.

## Provenance

This code was carried, byte-for-byte in its vendored Ooura sources, inside both
**MuTap** (adaptive filtering) and **AmbiTap** (ambisonics / binaural
convolution). The C++ wrappers had begun to diverge — MuTap grew the templated
`basic_real_fft<Sample>` and the CMSIS/vDSP backends; AmbiTap kept an older
double-engine wrapper with no backends — so a bug fix or a new backend in one
would silently miss the other. DspTap is the consolidation: one wrapper, one
contract, one home for the next backend. The unified wrapper is MuTap's
backend-capable `basic_real_fft`, generalized to the `dsptap` namespace.

See [`third_party/ooura/readme.txt`](third_party/ooura/readme.txt) and
[`third_party/cmsis-dsp/VENDOR.md`](third_party/cmsis-dsp/VENDOR.md) for the
vendored-code provenance and licenses.

## License

The DspTap wrapper is MIT (`LICENSE`). Vendored third-party code keeps its own
license: Ooura FFT (permissive, see its readme), CMSIS-DSP / CMSIS-Core
(Apache-2.0, SPDX headers retained in every file). See [`NOTICE.md`](NOTICE.md).
