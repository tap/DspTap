/* Single-precision instantiation of Ooura's fftsg.c.
 *
 * fftsg.c is written for double throughout. Embedded targets in MuTap's
 * real-time profile (e.g. Cortex-M55, whose FPU is single-precision only)
 * need a float FFT: this file compiles the same, unmodified fftsg.c source
 * a second time with `double` mapped to `float` and every file-scope symbol
 * suffixed `_f` so both precisions can link into one binary.
 *
 * <math.h> must be included BEFORE the `double` remap so the libm
 * declarations keep their real prototypes; fftsg.c's own include then
 * becomes a no-op via the header guard.
 *
 * The underlying algorithm and code remain Takuya Ooura's; see readme.txt
 * for the license/attribution terms.
 */

#include <math.h>

#define double float

/* Public API */
#define cdft cdft_f
#define rdft rdft_f
#define ddct ddct_f
#define ddst ddst_f
#define dfct dfct_f
#define dfst dfst_f

/* File-scope helpers (fftsg.c defines these as externally visible) */
#define makewt makewt_f
#define makeipt makeipt_f
#define makect makect_f
#define cftfsub cftfsub_f
#define cftbsub cftbsub_f
#define bitrv2 bitrv2_f
#define bitrv2conj bitrv2conj_f
#define bitrv216 bitrv216_f
#define bitrv216neg bitrv216neg_f
#define bitrv208 bitrv208_f
#define bitrv208neg bitrv208neg_f
#define cftf1st cftf1st_f
#define cftb1st cftb1st_f
#define cftrec4_th cftrec4_th_f
#define cftrec4 cftrec4_f
#define cfttree cfttree_f
#define cftleaf cftleaf_f
#define cftmdl1 cftmdl1_f
#define cftmdl2 cftmdl2_f
#define cftfx41 cftfx41_f
#define cftf161 cftf161_f
#define cftf162 cftf162_f
#define cftf081 cftf081_f
#define cftf082 cftf082_f
#define cftf040 cftf040_f
#define cftb040 cftb040_f
#define cftx020 cftx020_f
#define rftfsub rftfsub_f
#define rftbsub rftbsub_f
#define dctsub dctsub_f
#define dstsub dstsub_f

#include "fftsg.c"
