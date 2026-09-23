/// @file dsptap_capi.h
/// @brief Minimal C ABI over the DspTap primitives (fft, yin, psola, pvoc, log_mel, decimate), for language
///        bindings and the verification notebooks (notebooks/ drive it via ctypes).
///
///        Conventions: plain C types only; the caller owns all arrays and sizes them. Handle-based
///        functions return 0 on success and -1 on any error (bad argument, bad handle); a value
///        a function has to hand back (an exponent, a count) goes through an out-parameter or is
///        documented as a non-negative count where that predates this note (dsptap_yin_track,
///        dsptap_log_mel_process, dsptap_decimator_process). No global state. Everything runs the
///        double-precision golden profile — the notebooks verify the same code the consuming
///        libraries compile — except where a function documents a profile selector (the FFT), so
///        the notebooks can measure the embedded profiles as well.
///
///        Handles are typed opaque pointers, one struct type per primitive, so a C or C++
///        compiler rejects a yin handle passed to a pvoc function. No exception crosses this
///        boundary, from any entry point: every function is noexcept to a C++ caller; the create
///        functions (and the log_mel setters, which rebuild) catch allocation failure and return
///        NULL (-1 for a setter, which then leaves the handle as it was); every other function
///        performs no allocation, so there is nothing to throw. Geometry — every buffer a process
///        call touches — is fixed at create.
///
///        Changelog (contract-visible): 2026-09-23, Stage 3c — DSPTAP_FFT_PROFILE_Q15 / _Q31
///        un-reserved, _Q15_BFP / _Q31_BFP appended, dsptap_fft_fixed_scaling_exponent and the
///        four *_exp transforms added. dsptap_fft_forward_inplace_raw / _inverse_inplace_raw keep
///        their 0 / -1 contract and return -1 for a fixed-point handle, whose exponent only the
///        *_exp variants carry.
///        Stage 6 (docs/audit-fft-and-code-smells.md, Part 2 "capi"): the six handle typedefs are
///        pointers to distinct incomplete structs instead of void* (binary-compatible: every
///        handle is still one pointer, and the ctypes bridge's c_void_p is unchanged;
///        source-compatible for C; a C++ caller that stored a handle as void* now needs the
///        typedef); every entry point is DSPTAP_NOEXCEPT; dsptap_pvoc_create returns NULL above
///        fft_size 2^28 and dsptap_psola_create at max_period 2^26 or more (the headers'
///        preconditions, previously unchecked here); dsptap_decimator_process no longer
///        allocates; DSPTAP_API is dllexport only while building the library (DSPTAP_BUILDING)
///        and dllimport for a consumer on Windows.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/// Symbol visibility. On Windows the library's own build defines DSPTAP_BUILDING (tools/capi/
/// CMakeLists.txt, PRIVATE) and exports; every consumer that includes this header imports. On
/// other platforms the library builds with hidden visibility and these are the exported set.
#if defined(_WIN32)
#if defined(DSPTAP_BUILDING)
#define DSPTAP_API __declspec(dllexport)
#else
#define DSPTAP_API __declspec(dllimport)
#endif
#else
#define DSPTAP_API __attribute__((visibility("default")))
#endif

/// Every entry point is noexcept to a C++ caller (a C translation unit sees nothing). The create
/// functions and the log_mel setters catch and report failure (NULL / -1); the rest allocate
/// nothing and call only noexcept header functions and plain loops, so they cannot throw (if that
/// ever changed they would std::terminate rather than unwind into a C frame).
#ifdef __cplusplus
#define DSPTAP_NOEXCEPT noexcept
#else
#define DSPTAP_NOEXCEPT
#endif

/// Opaque handles: a pointer to a distinct incomplete type per primitive, defined only inside the
/// library. NULL is the invalid handle every function checks for.
typedef struct dsptap_fft_s*       dsptap_fft;
typedef struct dsptap_yin_s*       dsptap_yin;
typedef struct dsptap_psola_s*     dsptap_psola;
typedef struct dsptap_pvoc_s*      dsptap_pvoc;
typedef struct dsptap_log_mel_s*   dsptap_log_mel;
typedef struct dsptap_decimator_s* dsptap_decimator;

/// -- fft ------------------------------------------------------------------------------------

/// Numeric profile of a real FFT handle: which tap::dsp::basic_real_fft<Sample, Scaling>
/// instantiation it runs. DOUBLE is the golden model; FLOAT is the embedded floating profile
/// (the split-radix engine, or the platform backend the build selected — see
/// dsptap_fft_backend). Q15 and Q31
/// are the fixed-point profiles (std::int16_t in Q0.15, std::int32_t in Q0.31) under
/// tap::dsp::scaling::fixed; Q15_BFP and Q31_BFP the same two under
/// tap::dsp::scaling::block_floating (Stage 3b/3c of docs/audit-fft-and-code-smells.md). The
/// values are part of the ABI: existing ones never move, new ones are appended.
///
/// THE EXPONENT (fft.h, the fixed-point contract, restated here because every fixed-point entry
/// point below is defined in its terms). A fixed-point transform leaves its buffer scaled by
/// 2^-e and returns e: read the buffer as fractions of full scale (x / 2^15 for Q15, x / 2^31
/// for Q31) and let G be the DOUBLE profile on the same input read the same way — after a
/// forward, G's forward result == buffer * 2^e; after an inverse, G's UNNORMALIZED inverse
/// result == buffer * 2^e (up to the kernel's rounding noise; same packing, same exp(+i) sign).
/// Under fixed scaling e is the constant dsptap_fft_fixed_scaling_exponent(): log2 N for Q15
/// (forward output exactly X / N) and log2 N + 1 for Q31 (X / 2N); under block floating point
/// 0 <= e <= that constant, chosen per block from the data's headroom. The fixed-point inverse
/// applies NO 2/N (the exponent carries the scale): a raw round trip reconstructs
/// x == out * 2^(e_fwd + e_inv + 1 - log2 N) under both policies. The floating profiles have
/// no exponent and report e == 0 wherever one is returned. The fixed-point kernel is the
/// portable int32 kernel on every host (fft/fixed_point.h): its output for a given input is
/// the same bit pattern everywhere, whatever dsptap_fft_backend() says about float32.
///
/// Size: the floating profiles take any power of two >= 4; the four fixed-point profiles take
/// a power of two in [4, 65536] (fft/fixed_point.h, k_min_size / k_max_size, asserted by the
/// header's constructor), and dsptap_fft_create returns NULL above that bound rather than
/// hand out a handle the header does not support.
#define DSPTAP_FFT_PROFILE_DOUBLE 0
#define DSPTAP_FFT_PROFILE_FLOAT 1
#define DSPTAP_FFT_PROFILE_Q15 2
#define DSPTAP_FFT_PROFILE_Q31 3
#define DSPTAP_FFT_PROFILE_Q15_BFP 4
#define DSPTAP_FFT_PROFILE_Q31_BFP 5

/// Create a real FFT of `size` points in the given profile. All workspace is allocated here; the
/// transforms below are allocation-free.
/// @param size     a power of two >= 4; for the fixed-point profiles at most 65536 (see the
///                 profile codes)
/// @param profile  one of the DSPTAP_FFT_PROFILE_* codes
/// @return the handle, or NULL on a bad size, an unknown profile, or allocation failure
DSPTAP_API dsptap_fft dsptap_fft_create(int size, int profile) DSPTAP_NOEXCEPT;
DSPTAP_API void       dsptap_fft_destroy(dsptap_fft h) DSPTAP_NOEXCEPT;
DSPTAP_API int        dsptap_fft_size(dsptap_fft h) DSPTAP_NOEXCEPT;
/// size/2 + 1 (DC .. Nyquist).
DSPTAP_API int dsptap_fft_num_bins(dsptap_fft h) DSPTAP_NOEXCEPT;
/// The DSPTAP_FFT_PROFILE_* the handle was created with.
DSPTAP_API int dsptap_fft_profile(dsptap_fft h) DSPTAP_NOEXCEPT;
/// sizeof the profile's native sample (8 for DOUBLE, 4 for FLOAT, 2 for Q15 / Q15_BFP, 4 for
/// Q31 / Q31_BFP) — the element width of the buffers the *_inplace_raw transforms take.
DSPTAP_API int dsptap_fft_sample_bytes(dsptap_fft h) DSPTAP_NOEXCEPT;
/// The float32 engine this build compiled: "split_radix" (the C++20 engine of
/// fft/split_radix.h, bit-identical to the Ooura C it replaced; the string was "ooura" until
/// Stage 2c), "accelerate" (Apple vDSP) or "cmsis" (CMSIS-DSP Helium). The double profile is
/// always the split-radix engine; the four fixed-point profiles are always the portable int32
/// kernel (fft/fixed_point.h), on every host.
DSPTAP_API const char* dsptap_fft_backend(void) DSPTAP_NOEXCEPT;
/// The exponent contract's constant for this handle's size: basic_real_fft<Sample,
/// Scaling>::fixed_scaling_exponent(size) — the exponent every transform of a Q15 / Q31 handle
/// returns (log2 size, log2 size + 1) and the upper bound of what a Q15_BFP / Q31_BFP handle
/// returns.
/// @param h  the handle
/// @return the constant; 0 for DOUBLE and FLOAT (no exponent); -1 on a NULL handle
DSPTAP_API int dsptap_fft_fixed_scaling_exponent(dsptap_fft h) DSPTAP_NOEXCEPT;

/// Forward transform of size() real samples into the packed spectrum, in double regardless of
/// profile. `out` may alias `in`. Packing and sign convention are the header's contract
/// (fft.h): out[0] = DC, out[1] = Nyquist (both real), out[2k] + i*out[2k+1] = bin k for
/// 1 <= k < size/2, with W = exp(+2*pi*i/size) — CONJUGATE to the engineering convention
/// numpy.fft.rfft uses.
///
/// Per profile, at this double boundary:
///   - DOUBLE: the golden model, no conversion.
///   - FLOAT:  `in` is rounded to float32 first and the result widened back.
///   - Q15 / Q31 / Q15_BFP / Q31_BFP: `in` is read as fractions of full scale and QUANTIZED to
///     the profile's Q format by the substrate's round_sat (sample_traits.h: round half away
///     from zero, saturating at +-full scale, so |in| <= 1 - LSB is the saturation-free input
///     range); the transform runs on that native buffer and returns e; `out` is the native
///     result read as fractions and multiplied by 2^e. So `out` is directly comparable to the
///     DOUBLE profile's output on the same quantized input: out == G.forward(x_q) up to the
///     kernel's noise, in the DOUBLE profile's units, whatever the policy and whatever e was.
///     dsptap_fft_forward_exp is the same call with e written out.
DSPTAP_API int dsptap_fft_forward(dsptap_fft h, const double* in, double* out) DSPTAP_NOEXCEPT;
/// Inverse at the double boundary. `out` may alias `in`.
///   - DOUBLE / FLOAT: scaled by 2/size like basic_real_fft::inverse(), so forward -> inverse
///     reproduces the input.
///   - Q15 / Q31 / Q15_BFP / Q31_BFP: `in` (a packed spectrum in fractions of full scale) is
///     quantized to the Q format as for the forward, the fixed-point inverse runs, and `out` is
///     the native result read as fractions times 2^e. That is the UNNORMALIZED inverse: NO 2/N
///     is applied (fft.h: the fixed-point inverse() carries no 2/N; the exponent is the scale),
///     so `out` == G.inverse_inplace(a_q), the DOUBLE profile's raw inverse of the same quantized
///     spectrum. Note that a forward result scaled by 2^e can exceed full scale, so the
///     forward -> inverse round trip in double is not this function's job: scale the spectrum
///     into [-1, 1) yourself (e.g. by 2^-e_fwd) or use the raw path, whose round trip identity
///     is stated above. dsptap_fft_inverse_exp is the same call with e written out.
DSPTAP_API int dsptap_fft_inverse(dsptap_fft h, const double* in, double* out) DSPTAP_NOEXCEPT;
/// dsptap_fft_forward / dsptap_fft_inverse with the transform's exponent written out.
/// @param h         the handle
/// @param in        size() doubles (samples for the forward, a packed spectrum for the inverse)
/// @param out       size() doubles; may alias `in`
/// @param exponent  receives e: 0 for DOUBLE and FLOAT, the fixed-point e otherwise
/// @return 0, or -1 on a NULL handle, array or exponent pointer
DSPTAP_API int dsptap_fft_forward_exp(dsptap_fft h, const double* in, double* out, int* exponent) DSPTAP_NOEXCEPT;
DSPTAP_API int dsptap_fft_inverse_exp(dsptap_fft h, const double* in, double* out, int* exponent) DSPTAP_NOEXCEPT;

/// In-place transforms on the profile's NATIVE sample type: `data` points at size() samples of
/// dsptap_fft_sample_bytes() each (double for DOUBLE, float for FLOAT, int16_t for Q15 /
/// Q15_BFP, int32_t for Q31 / Q31_BFP), in the same packing, and must be aligned for that type
/// (a pointer the caller obtained as such; not an offset into a byte buffer). These are the
/// header's forward_inplace()/inverse_inplace() with no conversion at all, so the notebooks
/// measure the embedded profile's own arithmetic rather than a double round trip. The inverse
/// is UNSCALED, exactly as in fft.h: multiply by 2/size for a floating round trip; for a
/// fixed-point round trip apply the exponent identity stated at the profile codes.
///
/// These two are the ABI's original raw transforms and keep the file-level convention: for a
/// DOUBLE or FLOAT handle they transform `data` and return 0, exactly as before the fixed-point
/// profiles existed. For a fixed-point handle they do nothing and return -1: the transform's
/// exponent is part of its result and these signatures have nowhere to put it (fft.h: under
/// block floating point a discarded exponent is a silent scale error), so a fixed-point caller
/// must use the *_exp variants below. A caller that follows the 0 / -1 convention therefore can
/// never receive a 2^-e-scaled buffer and a success status.
/// @param h     the handle
/// @param data  size() native samples, aligned for the type; transformed in place
/// @return 0 on success; -1 on a NULL handle or buffer, and for any fixed-point handle
DSPTAP_API int dsptap_fft_forward_inplace_raw(dsptap_fft h, void* data) DSPTAP_NOEXCEPT;
DSPTAP_API int dsptap_fft_inverse_inplace_raw(dsptap_fft h, void* data) DSPTAP_NOEXCEPT;
/// The raw transforms for every profile, with the exponent written out: the path a binding uses.
/// @param h         the handle
/// @param data      size() native samples, aligned for the type; transformed in place
/// @param exponent  receives e: 0 for DOUBLE and FLOAT, the fixed-point e otherwise
/// @return 0, or -1 on a NULL handle, buffer or exponent pointer
DSPTAP_API int dsptap_fft_forward_inplace_raw_exp(dsptap_fft h, void* data, int* exponent) DSPTAP_NOEXCEPT;
DSPTAP_API int dsptap_fft_inverse_inplace_raw_exp(dsptap_fft h, void* data, int* exponent) DSPTAP_NOEXCEPT;

/// -- yin ------------------------------------------------------------------------------------

/// Create a detector (window, tau_min, tau_max as in tap::dsp::basic_yin), or NULL on bad geometry
/// (the header's 2 <= tau_min < tau_max <= window) or allocation failure.
DSPTAP_API dsptap_yin dsptap_yin_create(int window, int tau_min, int tau_max) DSPTAP_NOEXCEPT;
DSPTAP_API void       dsptap_yin_destroy(dsptap_yin h) DSPTAP_NOEXCEPT;
DSPTAP_API int        dsptap_yin_set_threshold(dsptap_yin h, double threshold) DSPTAP_NOEXCEPT;
DSPTAP_API int        dsptap_yin_frame_size(dsptap_yin h) DSPTAP_NOEXCEPT;

/// Analyze one frame of frame_size() samples (oldest first). Writes the fractional period in
/// samples (0 = unvoiced) and the normalized aperiodicity. Returns 0, or -1 on a bad handle.
DSPTAP_API int dsptap_yin_analyze(dsptap_yin h, const double* frame, double* period,
                                  double* aperiodicity) DSPTAP_NOEXCEPT;

/// Track a whole signal: analyze every `hop` samples and write up to max_out periods (0 where
/// unvoiced). Returns the number of analyses written, or -1 on error.
DSPTAP_API int dsptap_yin_track(dsptap_yin h, const double* x, int n, int hop, double* periods,
                                int max_out) DSPTAP_NOEXCEPT;

/// -- psola ----------------------------------------------------------------------------------

/// Create a TD-PSOLA shifter (tap::dsp::basic_psola<double>), or NULL on a max_period outside the
/// header's [16, 2^26) or on allocation failure.
DSPTAP_API dsptap_psola dsptap_psola_create(int max_period) DSPTAP_NOEXCEPT;
DSPTAP_API void         dsptap_psola_destroy(dsptap_psola h) DSPTAP_NOEXCEPT;
DSPTAP_API int          dsptap_psola_latency(dsptap_psola h) DSPTAP_NOEXCEPT;
DSPTAP_API int          dsptap_psola_clear(dsptap_psola h) DSPTAP_NOEXCEPT;

/// Shift n samples at a fixed source period and ratio (state persists across calls).
DSPTAP_API int dsptap_psola_process(dsptap_psola h, const double* in, double* out, int n, double period,
                                    double ratio) DSPTAP_NOEXCEPT;

/// -- pvoc -----------------------------------------------------------------------------------

/// Create a phase-vocoder shifter (tap::dsp::basic_pvoc<double>), or NULL unless fft_size is a
/// power of two in the header's [64, 2^28], or on allocation failure.
DSPTAP_API dsptap_pvoc dsptap_pvoc_create(int fft_size) DSPTAP_NOEXCEPT;
DSPTAP_API void        dsptap_pvoc_destroy(dsptap_pvoc h) DSPTAP_NOEXCEPT;
DSPTAP_API int         dsptap_pvoc_latency(dsptap_pvoc h) DSPTAP_NOEXCEPT;
DSPTAP_API int         dsptap_pvoc_set_formant(dsptap_pvoc h, int on) DSPTAP_NOEXCEPT;
DSPTAP_API int         dsptap_pvoc_clear(dsptap_pvoc h) DSPTAP_NOEXCEPT;

/// Shift n samples at a fixed ratio (state persists across calls).
DSPTAP_API int dsptap_pvoc_process(dsptap_pvoc h, const double* in, double* out, int n, double ratio) DSPTAP_NOEXCEPT;

/// -- log_mel --------------------------------------------------------------------------------

/// Create a log-mel front end at the given geometry (tap::dsp::log_mel_geometry; sqrt_window
/// selects the sqrt-Hann window, preemphasis 0 = off), or NULL on invalid geometry. The log
/// constants and PCEN default per the header; the setters below rebuild the object (geometry is
/// fixed at construction), which also resets it; a setter that fails (invalid values, or
/// allocation failure while rebuilding) returns -1 and leaves the handle exactly as it was.
DSPTAP_API dsptap_log_mel dsptap_log_mel_create(double sample_rate, int frame, int hop, int fft_size, int bands,
                                                double fmin_hz, double fmax_hz, int sqrt_window,
                                                double preemphasis) DSPTAP_NOEXCEPT;
DSPTAP_API void           dsptap_log_mel_destroy(dsptap_log_mel h) DSPTAP_NOEXCEPT;
DSPTAP_API int dsptap_log_mel_set_log(dsptap_log_mel h, double floor, double shift, double scale) DSPTAP_NOEXCEPT;
DSPTAP_API int dsptap_log_mel_set_pcen(dsptap_log_mel h, int enabled, double smoother, double alpha, double delta,
                                       double power, double epsilon) DSPTAP_NOEXCEPT;
DSPTAP_API int dsptap_log_mel_reset(dsptap_log_mel h) DSPTAP_NOEXCEPT;
DSPTAP_API int dsptap_log_mel_bands(dsptap_log_mel h) DSPTAP_NOEXCEPT;
DSPTAP_API int dsptap_log_mel_latency(dsptap_log_mel h) DSPTAP_NOEXCEPT;
/// The formula-level contract version of log_mel.h (trained models record it).
DSPTAP_API int dsptap_log_mel_contract_version(void) DSPTAP_NOEXCEPT;

/// Stream n samples; writes up to max_frames frames of bands() features (row-major) and returns the
/// number written, or -1 on error. Partial hops carry over to the next call.
DSPTAP_API int dsptap_log_mel_process(dsptap_log_mel h, const double* x, int n, double* features,
                                      int max_frames) DSPTAP_NOEXCEPT;

/// -- decimate -------------------------------------------------------------------------------

/// Create a decimator by ratio 2, 3 or 6 (tap::dsp::basic_decimator<float, M>; transparent selects
/// that profile over economy), or NULL on a bad ratio or allocation failure. The float profile is
/// the deployed profile on the FIR substrate (decimate.h's double golden model is not exposed
/// here); the double arrays are converted at the boundary through float staging buffers
/// allocated here, so dsptap_decimator_process allocates nothing for any n.
DSPTAP_API dsptap_decimator dsptap_decimator_create(int ratio, int transparent) DSPTAP_NOEXCEPT;
DSPTAP_API void             dsptap_decimator_destroy(dsptap_decimator h) DSPTAP_NOEXCEPT;
DSPTAP_API int              dsptap_decimator_taps(dsptap_decimator h) DSPTAP_NOEXCEPT;
DSPTAP_API int              dsptap_decimator_latency(dsptap_decimator h) DSPTAP_NOEXCEPT;
DSPTAP_API int              dsptap_decimator_reset(dsptap_decimator h) DSPTAP_NOEXCEPT;
/// Outputs the next call with n inputs will produce.
DSPTAP_API int dsptap_decimator_outputs_for(dsptap_decimator h, int n) DSPTAP_NOEXCEPT;

/// Decimate n samples; writes up to max_out outputs and returns the number written (all of
/// outputs_for(n) when max_out allows; excess outputs are dropped, and all n inputs are consumed
/// either way), or -1 on error. Any n: the input is staged in fixed-size blocks, and the result
/// is bit-identical to one call over the whole stream (decimate.h: chunking-invariant).
DSPTAP_API int dsptap_decimator_process(dsptap_decimator h, const double* in, int n, double* out,
                                        int max_out) DSPTAP_NOEXCEPT;

#ifdef __cplusplus
}
#endif
