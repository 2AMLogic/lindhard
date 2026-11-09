//! GPU feasibility prototype (#425): the hot kernel of the BCA engine, a
//! scattering-table lookup plus a Lindhard-Scharff stopping evaluation, on
//! the CPU (the `lindhard` library, f64, one thread) and as a wgpu compute
//! shader (f32). Prints a table of throughput and the f32-vs-f64 deviation.
//! Method and results: docs/gpu-feasibility.md.

use lindhard::ion::potential::{Potential, Screening};
use lindhard::ion::scattering::{ScatteringTable, TableSpec};
use std::time::Instant;
use wgpu::util::DeviceExt;

const K_L: f64 = 0.1; // order-of-magnitude reduced LS coefficient; the value is irrelevant to throughput

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
    n_eps: u32,
    n_beta: u32,
    n: u32,
    pad: u32,
    ln_eps_min: f32,
    ln_beta_min: f32,
    inv_d: f32,
    k_l: f32,
}

fn sample_points(n: usize) -> Vec<(f64, f64)> {
    // Same generator and ranges as lindhard/benches/scattering.rs.
    let mut st = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = || {
        st = st
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (st >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..n)
        .map(|_| {
            (
                1e-3 * 10f64.powf(5.0 * next()),
                1e-3 * 10f64.powf(4.5 * next()),
            )
        })
        .collect()
}

fn cpu(table: &ScatteringTable, pts: &[(f64, f64)]) -> Vec<(f64, f64)> {
    pts.iter()
        .map(|&(e, b)| (table.half_angle_tan(e, b).unwrap(), K_L * e.sqrt()))
        .collect()
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn main() {
    let pot = Potential::new(Screening::ZblUniversal, 14.0, 14.0);
    let spec = TableSpec {
        eps_min: 1e-3,
        eps_max: 1e2,
        beta_min: 1e-3,
        beta_max: 1e2,
        per_decade: 32,
    };
    let table = ScatteringTable::build(&pot, &spec);
    let (n_eps, n_beta) = table.shape();
    let d = std::f64::consts::LN_10 / spec.per_decade as f64;
    // The stored quantity is ln tan(theta/2) at the grid nodes; interpolation
    // at a node returns it exactly, so read it back through the public API.
    let mut y = Vec::with_capacity(n_eps * n_beta);
    for i in 0..n_eps {
        let eps = (spec.eps_min * (i as f64 * d).exp()).min(spec.eps_max);
        for j in 0..n_beta {
            let beta = spec.beta_min * (j as f64 * d).exp();
            let t = table.half_angle_tan(eps, beta.min(spec.beta_max * 0.999_999));
            y.push(t.unwrap().ln() as f32);
        }
    }
    println!(
        "cpu: {} threads available; table {n_eps}x{n_beta}",
        std::thread::available_parallelism().map_or(1, |n| n.get())
    );

    let sizes = [1usize << 14, 1 << 18, 1 << 21];
    let inst = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
    let adapter = pollster::block_on(inst.request_adapter(&wgpu::RequestAdapterOptions {
        power_preference: wgpu::PowerPreference::HighPerformance,
        ..Default::default()
    }));
    let gpu = adapter.ok().map(|a| {
        let info = a.get_info();
        if info.device_type != wgpu::DeviceType::DiscreteGpu
            && info.device_type != wgpu::DeviceType::IntegratedGpu
        {
            println!("WARNING: adapter is not a GPU; the wgpu columns are NOT GPU numbers");
        }
        println!(
            "adapter: {} ({:?}, {:?})",
            info.name, info.device_type, info.backend
        );
        let (dev, q) = pollster::block_on(a.request_device(&wgpu::DeviceDescriptor::default()))
            .expect("device");
        (info, dev, q)
    });
    if gpu.is_none() {
        println!("adapter: none (no GPU or software rasteriser visible)");
    }

    println!(
        "{:>10} {:>14} {:>14} {:>14} {:>12}",
        "points", "cpu 1T Mpt/s", "gpu kern Mpt/s", "gpu e2e Mpt/s", "max rel dev"
    );
    for &n in &sizes {
        let pts = sample_points(n);
        let reps = (1usize << 24) / n + 1;
        let ts: Vec<f64> = (0..reps.max(3))
            .map(|_| {
                let t = Instant::now();
                std::hint::black_box(cpu(&table, std::hint::black_box(&pts)));
                t.elapsed().as_secs_f64()
            })
            .collect();
        let cpu_rate = n as f64 / median(ts) / 1e6;
        let reference = cpu(&table, &pts);

        let (mut kern, mut e2e, mut dev_max) = (f64::NAN, f64::NAN, f64::NAN);
        if let Some((_, dev, queue)) = &gpu {
            let r = run_gpu(dev, queue, &spec, n_eps, n_beta, d, &y, &pts);
            kern = n as f64 / r.0 / 1e6;
            e2e = n as f64 / r.1 / 1e6;
            dev_max =
                r.2.iter()
                    .zip(&reference)
                    .map(|(g, c)| {
                        ((f64::from(g.0) - c.0).abs() / c.0).max((f64::from(g.1) - c.1).abs() / c.1)
                    })
                    .fold(0.0, f64::max);
        }
        println!("{n:>10} {cpu_rate:>14.1} {kern:>14.1} {e2e:>14.1} {dev_max:>12.2e}");
    }
}

/// Returns (median kernel-only seconds, median upload+kernel+readback seconds, outputs).
#[allow(clippy::too_many_arguments)]
fn run_gpu(
    dev: &wgpu::Device,
    queue: &wgpu::Queue,
    spec: &TableSpec,
    n_eps: usize,
    n_beta: usize,
    d: f64,
    y: &[f32],
    pts: &[(f64, f64)],
) -> (f64, f64, Vec<(f32, f32)>) {
    let n = pts.len();
    let params = Params {
        n_eps: n_eps as u32,
        n_beta: n_beta as u32,
        n: n as u32,
        pad: 0,
        ln_eps_min: spec.eps_min.ln() as f32,
        ln_beta_min: spec.beta_min.ln() as f32,
        inv_d: (1.0 / d) as f32,
        k_l: K_L as f32,
    };
    let inp: Vec<[f32; 2]> = pts.iter().map(|&(e, b)| [e as f32, b as f32]).collect();
    let buf = |label, contents: &[u8], usage| {
        dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some(label),
            contents,
            usage,
        })
    };
    let pbuf = buf(
        "p",
        bytemuck::bytes_of(&params),
        wgpu::BufferUsages::UNIFORM,
    );
    let tbuf = buf("t", bytemuck::cast_slice(y), wgpu::BufferUsages::STORAGE);
    let ibuf = buf(
        "i",
        bytemuck::cast_slice(&inp),
        wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
    );
    let bytes = (n * 8) as u64;
    let obuf = dev.create_buffer(&wgpu::BufferDescriptor {
        label: Some("o"),
        size: bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });
    let rbuf = dev.create_buffer(&wgpu::BufferDescriptor {
        label: Some("r"),
        size: bytes,
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let module = dev.create_shader_module(wgpu::include_wgsl!("shader.wgsl"));
    let pipe = dev.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: None,
        layout: None,
        module: &module,
        entry_point: Some("main"),
        compilation_options: Default::default(),
        cache: None,
    });
    let bg = dev.create_bind_group(&wgpu::BindGroupDescriptor {
        label: None,
        layout: &pipe.get_bind_group_layout(0),
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: pbuf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: tbuf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: ibuf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 3,
                resource: obuf.as_entire_binding(),
            },
        ],
    });
    let dispatch = |readback: bool| {
        let mut enc = dev.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            pass.set_pipeline(&pipe);
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups((n as u32).div_ceil(64), 1, 1);
        }
        if readback {
            enc.copy_buffer_to_buffer(&obuf, 0, &rbuf, 0, bytes);
        }
        queue.submit([enc.finish()]);
        let _ = dev.poll(wgpu::PollType::wait_indefinitely());
    };
    // NB: n <= 2^21 keeps the dispatch within the 65535 workgroup limit in x, so
    // the shader indexes on gid.x only.
    assert!(n.div_ceil(64) <= 65535);
    dispatch(false); // warm-up / pipeline compile
    let mut kt = Vec::new();
    let mut et = Vec::new();
    for _ in 0..7 {
        let t = Instant::now();
        dispatch(false);
        kt.push(t.elapsed().as_secs_f64());
        let t = Instant::now();
        queue.write_buffer(&ibuf, 0, bytemuck::cast_slice(&inp));
        dispatch(true);
        let s = rbuf.slice(..);
        s.map_async(wgpu::MapMode::Read, |r| r.unwrap());
        let _ = dev.poll(wgpu::PollType::wait_indefinitely());
        et.push(t.elapsed().as_secs_f64());
        rbuf.unmap();
    }
    dispatch(true);
    let s = rbuf.slice(..);
    s.map_async(wgpu::MapMode::Read, |r| r.unwrap());
    let _ = dev.poll(wgpu::PollType::wait_indefinitely());
    let out: Vec<(f32, f32)> = bytemuck::cast_slice::<u8, [f32; 2]>(&s.get_mapped_range().unwrap())
        .iter()
        .map(|v| (v[0], v[1]))
        .collect();
    rbuf.unmap();
    (median(kt), median(et), out)
}
