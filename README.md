# DspTap

Shared DSP primitives for the **Tap** family of audio libraries. Header-only,
plain portable C++ (C++20, standard library only), no Max/Min or framework
dependency — consumed as a git submodule by the individual libraries.

Today it holds four primitives:

## `tap::dsp::real_fft` — real FFT with a fixed numeric contract

`include/tap/dsp/fft.h` wraps the vendored [Ooura split-radix real
FFT](third_party/ooura/) behind a small, well-specified interface, with
**optional, mutually-exclusive float32 backends** that re-present the *exact*
same numeric contract for speed on specific hardware:

| Backend | Build option | Target | Notes |
|---------|-------------|--------|-------|
| Ooura (default) | — | everywhere | the golden model; double is **always** Ooura |
| CMSIS-DSP Helium | `TAP_DSP_FFT_CMSIS` | bare-metal Cortex-M55 (MVE) | ~3× fewer instructions/transform |
| Apple vDSP | `TAP_DSP_FFT_ACCELERATE` | macOS / Apple Silicon | ~3× faster/transform |

The two float32 backends conjugate imaginary bins and rescale so every
intermediate spectrum matches the Ooura build to single-precision rounding —
so the whole double-precision test battery stays a valid oracle for the
accelerated float paths. `tests/test_fft_backend.cpp` pins each backend to
Ooura's `rdft_f` bin-for-bin at the certified geometries (N = 512, 2048).

```cpp
#include "tap/dsp/fft.h"

tap::dsp::real_fft   fft(1024);   // double, the desktop/golden profile
tap::dsp::real_fft32 fft32(1024); // float,  the embedded / accelerated profile

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

## `tap::dsp::yin` — YIN pitch detector

`include/tap/dsp/yin.h` implements the time-domain YIN estimator (de Cheveigné
& Kawahara 2002, steps 1–5): squared-difference function, cumulative-mean
normalization, absolute threshold with local-minimum descent, and parabolic
interpolation for sub-sample period precision. Header-only, allocation-free
after construction, `noexcept` on the analysis path.

```cpp
#include "tap/dsp/yin.h"

tap::dsp::yin   det(800, 20, 800);   // window, tau_min, tau_max — double golden model
tap::dsp::yin32 det32(800, 20, 800); // float, the embedded profile

const auto r = det.analyze(frame);   // frame_size() == window + tau_max samples
if (r.voiced()) {
    const double freq = sample_rate / r.period; // fractional-sample period
}
```

Key contract points (full detail in the header docstring):

- **Geometry** fixed at construction: integration `window`, searched lag range
  `[tau_min, tau_max]` (bound from your frequency range as `τ = sr / f`), with
  `window ≥ tau_max`; `analyze()` reads `window + tau_max` samples, oldest first.
- **Result**: fractional period in samples (0 = unvoiced) plus the normalized
  aperiodicity at the chosen lag (global minimum when unvoiced).
- **Threshold**: the paper's absolute threshold on the normalized difference,
  default 0.1, settable at runtime.
- Both precisions run the identical algorithm; the test battery pins sine/
  sawtooth accuracy (sub-cent in double), octave robustness, unvoiced
  rejection, and float/double agreement. The difference-function inner loop is
  the designated Helium-MVE / HVX backend candidate behind this same contract,
  mirroring the FFT's golden-model-plus-backends pattern.

## `tap::dsp::psola` — pitch-synchronous overlap-add shifter

`include/tap/dsp/psola.h` is the real-time TD-PSOLA resynthesis stage:
Hann grains two source periods long, extracted at period-spaced analysis marks
and overlap-added at period/ratio-spaced synthesis marks with sub-sample
(Hermite) placement. Detection-agnostic — the caller supplies the period (from
`tap::dsp::yin` or any tracker); fixed latency of `2 * max_period + 2` samples.

```cpp
tap::dsp::psola shifter(900);            // deepest period it will be given
double y = shifter.process(x, period, ratio); // per sample; ratio 2 = octave up
```

Know what PSOLA is: it resamples the source's **spectral envelope** at the new
harmonic spacing — which is why it preserves formants on voice, and why a pure
tone shifted far from any new harmonic thins toward silence. Feed it
harmonic-rich material; both behaviors are pinned by the tests.

## `tap::dsp::pvoc` — phase-vocoder pitch shifter

`include/tap/dsp/pvoc.h` is an STFT pitch shifter (Hann, 4× overlap, built on
`tap::dsp::real_fft` so the float profile rides the vDSP/CMSIS backends) using
Laroche–Dolson-style **peak-region shifting**: each spectral peak's region is
translated rigidly by an integer bin offset and rotated by one accumulated
residual phase, so phase relationships across the peak stay intact. At
ratio 1 the output reconstructs the input's waveform delayed by exactly one
FFT frame (pinned by the tests). Latency = the FFT size (1024 default).

```cpp
tap::dsp::pvoc shifter(1024);
shifter.set_formant(true);               // optional LPC formant preservation
double y = shifter.process(x, ratio);    // per sample; ratio sampled per hop
```

Optional **formant preservation** (`set_formant`) uses the classic source-filter
method: an LPC spectral envelope (autocorrelation + Levinson–Durbin, order 48)
per analysis frame, with every relocated bin rescaled by
`envelope(target)/envelope(source)` — the excitation moves, the envelope stays.
At ratio 1 the correction is exactly unity, so the identity contract holds
either way.

## Notebooks

`notebooks/pitchshift.ipynb` measures the three pitch primitives — driving the
**actual shipping C++** through the C ABI in `tools/capi/` (ctypes bridge:
`notebooks/dsptap_py.py`, which builds `build_capi/` on first import). It
documents the two findings from the primitives' development: PSOLA's
envelope-resampling nature (why it preserves formants *and* why a pure tone
shifted an octave thins out), and the measured level collapse of naive
phase-vocoder bin remapping vs the shipping peak-locked design — plus the LPC
formant-preservation demo.

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
target_link_libraries(my_dsp PRIVATE tap::dsp)
```

`tap::dsp` is an INTERFACE target (the headers + the compiled Ooura
static lib `tap::dsp_fft`); it does not build the tests when added as a
subdirectory (`TAP_DSP_BUILD_TESTS` defaults OFF unless top-level). The
per-platform float32 backend defaults follow the target: vDSP on Apple, CMSIS
on the bare-metal M55 profile, Ooura elsewhere — override with
`-DTAP_DSP_FFT_ACCELERATE=OFF` etc.

## Provenance

This code was carried, byte-for-byte in its vendored Ooura sources, inside both
**MuTap** (adaptive filtering) and **AmbiTap** (ambisonics / binaural
convolution). The C++ wrappers had begun to diverge — MuTap grew the templated
`basic_real_fft<Sample>` and the CMSIS/vDSP backends; AmbiTap kept an older
double-engine wrapper with no backends — so a bug fix or a new backend in one
would silently miss the other. DspTap is the consolidation: one wrapper, one
contract, one home for the next backend. The unified wrapper is MuTap's
backend-capable `basic_real_fft`, generalized to the `tap::dsp` namespace.

See [`third_party/ooura/readme.txt`](third_party/ooura/readme.txt) and
[`third_party/cmsis-dsp/VENDOR.md`](third_party/cmsis-dsp/VENDOR.md) for the
vendored-code provenance and licenses.

## License

The DspTap wrapper is MIT (`LICENSE`). Vendored third-party code keeps its own
license: Ooura FFT (permissive, see its readme), CMSIS-DSP / CMSIS-Core
(Apache-2.0, SPDX headers retained in every file). See [`NOTICE.md`](NOTICE.md).
