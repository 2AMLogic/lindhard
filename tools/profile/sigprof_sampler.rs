// SIGPROF sampler used for docs/img/flamegraph-b-si.svg (see
// docs/benchmarks.md, "Profile"). It is not part of the workspace: build it
// as its own crate with `lindhard` (path), `libc` and `rayon` as
// dependencies, `RUSTFLAGS="-C force-frame-pointers=yes"`, release profile
// with `debug = "line-tables-only"` and `strip = "none"`. Written for this
// repository (MIT); needed because `perf_event_paranoid = 4` blocked `perf`.
// x86_64 Linux only: it walks the rbp chain from the signal context.
#![allow(dead_code)]
#[path = "<repo>/lindhard/benches/common/mod.rs"] // set <repo> to the checkout
mod common;
use common::*;
use lindhard::ion::bca::{Bca, SummaryTally};
use lindhard::ion::stopping::lindhard_scharff::LindhardScharff;
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};

const MAXD: usize = 48;
const MAXS: usize = 200_000;
static mut BUF: [[u64; MAXD]; MAXS] = [[0; MAXD]; MAXS];
static NS: AtomicUsize = AtomicUsize::new(0);

extern "C" fn handler(_: i32, _: *mut libc::siginfo_t, ctx: *mut libc::c_void) {
    unsafe {
        let uc = &*(ctx as *const libc::ucontext_t);
        let rip = uc.uc_mcontext.gregs[libc::REG_RIP as usize] as u64;
        let mut rbp = uc.uc_mcontext.gregs[libc::REG_RBP as usize] as u64;
        let i = NS.fetch_add(1, Relaxed);
        if i >= MAXS { return; }
        let b = &mut (*std::ptr::addr_of_mut!(BUF))[i];
        b[0] = rip;
        let mut d = 1;
        while d < MAXD && rbp != 0 && rbp % 8 == 0 {
            let p = rbp as *const u64;
            let ret = *p.add(1);
            let next = *p;
            if ret == 0 { break; }
            b[d] = ret; d += 1;
            if next <= rbp || next - rbp > (1 << 24) { break; }
            rbp = next;
        }
    }
}

fn main() {
    let reps: usize = std::env::args().nth(1).map(|s| s.parse().unwrap()).unwrap_or(5);
    let ls = LindhardScharff::new();
    let p = &PROBLEMS[0];
    let st = stack(p);
    let bca = Bca::new(beam(p, 10_000), &st, config(1), &ls, table()).unwrap();
    let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
    // table built; start sampling
    unsafe {
        let mut sa: libc::sigaction = std::mem::zeroed();
        sa.sa_sigaction = handler as usize;
        sa.sa_flags = libc::SA_SIGINFO | libc::SA_RESTART;
        libc::sigaction(libc::SIGPROF, &sa, std::ptr::null_mut());
        let it = libc::itimerval {
            it_interval: libc::timeval { tv_sec: 0, tv_usec: 1000 },
            it_value: libc::timeval { tv_sec: 0, tv_usec: 1000 },
        };
        libc::setitimer(libc::ITIMER_PROF, &it, std::ptr::null_mut());
    }
    let t = std::time::Instant::now();
    let mut h = 0;
    for _ in 0..reps {
        h += pool.install(|| bca.run(|| SummaryTally::new(1e-9, 100)).unwrap().histories);
    }
    let el = t.elapsed().as_secs_f64();
    unsafe { libc::setitimer(libc::ITIMER_PROF, &std::mem::zeroed(), std::ptr::null_mut()); }
    eprintln!("{h} ions in {el:.3} s = {:.0} ions/s, {} samples", h as f64 / el, NS.load(Relaxed));
    // base of the main binary mapping
    let maps = std::fs::read_to_string("/proc/self/maps").unwrap();
    let exe = std::env::current_exe().unwrap();
    let base = maps.lines().filter(|l| l.ends_with(exe.to_str().unwrap())).map(|l| u64::from_str_radix(l.split('-').next().unwrap(), 16).unwrap()).min().unwrap();
    std::fs::write("maps.txt", &maps).unwrap();
    let mut out = std::io::BufWriter::new(std::fs::File::create("samples.txt").unwrap());
    writeln!(out, "base {base:x}").unwrap();
    let n = NS.load(Relaxed).min(MAXS);
    for i in 0..n {
        let b = unsafe { &(*std::ptr::addr_of!(BUF))[i] };
        let s: Vec<String> = b.iter().take_while(|&&a| a != 0).map(|a| format!("{:x}", a)).collect();
        writeln!(out, "{}", s.join(" ")).unwrap();
    }
}
