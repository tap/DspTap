/* Minimal QEMU TCG plugin: counts executed guest instructions and prints
 * one line at exit:
 *
 *   TAP_DSP_INSN_COUNT <n>
 *
 * Used by scripts/icount.py for the deterministic performance ratchet
 * (bench/README.md). Counting uses the inline-add fast path; the single
 * counter is exact for our single-vCPU deterministic workloads.
 *
 * Build (qemu-plugin.h fetched for the matching QEMU 8.2.x, whose header
 * defines QEMU_PLUGIN_VERSION 1 — the plugin API version is the header's, and
 * the pinned v8.2.2 header says 1, not 2 as MuTap's copy of this comment had it):
 *   gcc -shared -fPIC $(pkg-config --cflags glib-2.0) \
 *       -I<dir with qemu-plugin.h> insn_count.c -o libinsncount.so
 *
 * Licensing: qemu-plugin.h is QEMU's, SPDX GPL-2.0-or-later. It is fetched at
 * CI time (digest-verified), never vendored into this repo, and this file is
 * compiled against it only to build a test tool that runs in CI and is not
 * shipped; nothing in the DspTap tree or in what consumers link is GPL.
 * This file itself is MIT.
 *
 * Provenance: copied from MuTap's tools/qemu_insn_plugin/insn_count.c (MIT,
 * MuTap contributors) with only the output marker renamed; the same file
 * pattern lives in SampleRateTap and RatioTap. A taphouse-style
 * consolidation is the eventual home (docs/audit-fft-and-code-smells.md,
 * Part 11, "Sharing").
 */
// SPDX-License-Identifier: MIT
// Copyright 2026 MuTap contributors
// Copyright 2026 Timothy Place and the DspTap contributors.
#include <glib.h>
#include <inttypes.h>
#include <qemu-plugin.h>
#include <stdint.h>

QEMU_PLUGIN_EXPORT int qemu_plugin_version = QEMU_PLUGIN_VERSION;

static uint64_t insn_count;

static void tb_trans(qemu_plugin_id_t id, struct qemu_plugin_tb* tb) {
    (void)id;
    size_t n = qemu_plugin_tb_n_insns(tb);
    for (size_t i = 0; i < n; i++) {
        struct qemu_plugin_insn* insn = qemu_plugin_tb_get_insn(tb, i);
        qemu_plugin_register_vcpu_insn_exec_inline(insn, QEMU_PLUGIN_INLINE_ADD_U64, &insn_count, 1);
    }
}

static void at_exit(qemu_plugin_id_t id, void* userdata) {
    (void)id;
    (void)userdata;
    g_autofree gchar* msg = g_strdup_printf("TAP_DSP_INSN_COUNT %" PRIu64 "\n", insn_count);
    qemu_plugin_outs(msg);
}

QEMU_PLUGIN_EXPORT int qemu_plugin_install(qemu_plugin_id_t id, const qemu_info_t* info, int argc, char** argv) {
    (void)info;
    (void)argc;
    (void)argv;
    qemu_plugin_register_vcpu_tb_trans_cb(id, tb_trans);
    qemu_plugin_register_atexit_cb(id, at_exit, NULL);
    return 0;
}
