/* Minimal QEMU TCG plugin: counts executed guest instructions and prints
 * one line at exit:
 *
 *   TAP_DSP_INSN_COUNT <n>
 *
 * Used by scripts/icount.py for the deterministic performance ratchet
 * (bench/README.md). Counting uses the inline-add fast path; the single
 * counter is exact for our single-vCPU deterministic workloads.
 *
 * A workload that brackets calls with bench/call_marker.h's
 * tap_dsp_call_begin() / tap_dsp_call_end() (the layout workloads,
 * bench/layout/) gets a second line:
 *
 *   TAP_DSP_CALL_MAX <max> call=<k> calls=<n> total=<sum>
 *
 * the most instructions executed between one begin and the next end, which
 * of the n bracketed calls it was (1-based), and their sum. The markers are
 * found by symbol (qemu_plugin_insn_symbol, the guest image's symbol
 * table), at their entry: the address of the first marker instruction in
 * the first translation block that enters it, latched, so only that
 * address gets the callback in every later translation. Without the latch
 * a marker straddling a page boundary (two translation blocks) fired twice
 * per call (SampleRateTap's plugin, fixed the same way; reproduced on
 * mps2-an505). Only the two entry instructions get a callback, so a
 * workload without the markers (every ratchet scenario) counts exactly as
 * before. The mechanism is
 * SampleRateTap's (its tools/qemu_insn_plugin, the roadmap's A6), with the
 * marker names and the output line renamed.
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
#include <string.h>

QEMU_PLUGIN_EXPORT int qemu_plugin_version = QEMU_PLUGIN_VERSION;

static uint64_t insn_count;

/* The bracketed calls (the file comment). */
static uint64_t call_begin_at;
static uint64_t call_max;
static uint64_t call_max_at;
static uint64_t call_n;
static uint64_t call_total;

static void on_call_begin(unsigned int vcpu, void* userdata) {
    (void)vcpu;
    (void)userdata;
    call_begin_at = insn_count;
}

static void on_call_end(unsigned int vcpu, void* userdata) {
    (void)vcpu;
    (void)userdata;
    const uint64_t d = insn_count - call_begin_at;
    ++call_n;
    call_total += d;
    if (d > call_max) {
        call_max    = d;
        call_max_at = call_n;
    }
}

/* Each marker's entry address, latched at its first translation (the file
 * comment); 0 until then. Translation is serialized; the workloads are
 * single-vCPU. */
static uint64_t begin_entry;
static uint64_t end_entry;

static int is_entry(uint64_t* entry, uint64_t vaddr) {
    if (*entry == 0) {
        *entry = vaddr;
    }
    return *entry == vaddr;
}

static void tb_trans(qemu_plugin_id_t id, struct qemu_plugin_tb* tb) {
    (void)id;
    size_t      n    = qemu_plugin_tb_n_insns(tb);
    const char* prev = NULL;
    for (size_t i = 0; i < n; i++) {
        struct qemu_plugin_insn* insn = qemu_plugin_tb_get_insn(tb, i);
        qemu_plugin_register_vcpu_insn_exec_inline(insn, QEMU_PLUGIN_INLINE_ADD_U64, &insn_count, 1);
        const char* sym = qemu_plugin_insn_symbol(insn);
        if (sym != NULL && (prev == NULL || strcmp(sym, prev) != 0)) {
            const uint64_t vaddr = qemu_plugin_insn_vaddr(insn);
            if (strcmp(sym, "tap_dsp_call_begin") == 0 && is_entry(&begin_entry, vaddr)) {
                qemu_plugin_register_vcpu_insn_exec_cb(insn, on_call_begin, QEMU_PLUGIN_CB_NO_REGS, NULL);
            }
            else if (strcmp(sym, "tap_dsp_call_end") == 0 && is_entry(&end_entry, vaddr)) {
                qemu_plugin_register_vcpu_insn_exec_cb(insn, on_call_end, QEMU_PLUGIN_CB_NO_REGS, NULL);
            }
        }
        prev = sym;
    }
}

static void at_exit(qemu_plugin_id_t id, void* userdata) {
    (void)id;
    (void)userdata;
    g_autofree gchar* msg = g_strdup_printf("TAP_DSP_INSN_COUNT %" PRIu64 "\n", insn_count);
    qemu_plugin_outs(msg);
    if (call_n != 0) {
        g_autofree gchar* calls =
            g_strdup_printf("TAP_DSP_CALL_MAX %" PRIu64 " call=%" PRIu64 " calls=%" PRIu64 " total=%" PRIu64 "\n",
                            call_max, call_max_at, call_n, call_total);
        qemu_plugin_outs(calls);
    }
}

QEMU_PLUGIN_EXPORT int qemu_plugin_install(qemu_plugin_id_t id, const qemu_info_t* info, int argc, char** argv) {
    (void)info;
    (void)argc;
    (void)argv;
    qemu_plugin_register_vcpu_tb_trans_cb(id, tb_trans);
    qemu_plugin_register_atexit_cb(id, at_exit, NULL);
    return 0;
}
