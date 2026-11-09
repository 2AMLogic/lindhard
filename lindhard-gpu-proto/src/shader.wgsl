// Prototype kernel: bilinear lookup in the log-tan(theta/2) scattering table
// (same scheme as lindhard::ion::scattering::ScatteringTable::half_angle_tan)
// plus Lindhard-Scharff reduced stopping k_L * sqrt(eps). f32 throughout:
// WGSL has no portable f64.
struct Params {
    n_eps: u32,
    n_beta: u32,
    n: u32,
    _pad: u32,
    ln_eps_min: f32,
    ln_beta_min: f32,
    inv_d: f32,
    k_l: f32,
};
@group(0) @binding(0) var<uniform> p: Params;
@group(0) @binding(1) var<storage, read> table: array<f32>;
@group(0) @binding(2) var<storage, read> inp: array<vec2<f32>>; // (eps, beta)
@group(0) @binding(3) var<storage, read_write> out: array<vec2<f32>>; // (tan(theta/2), S_e)

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    if (k >= p.n) { return; }
    let q = inp[k];
    let fi = (log(q.x) - p.ln_eps_min) * p.inv_d;
    let fj = (log(q.y) - p.ln_beta_min) * p.inv_d;
    let i = min(u32(floor(fi)), p.n_eps - 2u);
    let j = min(u32(floor(fj)), p.n_beta - 2u);
    let ti = fi - f32(i);
    let tj = fj - f32(j);
    let a = table[i * p.n_beta + j];
    let b = table[i * p.n_beta + j + 1u];
    let c = table[(i + 1u) * p.n_beta + j];
    let d = table[(i + 1u) * p.n_beta + j + 1u];
    let y = (1.0 - ti) * ((1.0 - tj) * a + tj * b) + ti * ((1.0 - tj) * c + tj * d);
    out[k] = vec2<f32>(exp(y), p.k_l * sqrt(q.x));
}
