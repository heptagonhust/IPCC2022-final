use rayon::slice::ParallelSliceMut;
use ::safer_ffi::prelude::*;
#[derive_ReprC]
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Edge {
  a: i32,
  b: i32,
  weight: f64, 
  origin_weight: f64,
  lca: i32,
}

#[ffi_export]
pub fn parallel_sort_edges(edges: c_slice::Mut<Edge>) {
  edges.as_slice().par_sort_by(|a, b| {a.weight.partial_cmp(&b.weight).unwrap().reverse()});
}

#[ffi_export]
pub fn sort_edges(edges: c_slice::Mut<Edge>) {
  edges.as_slice().sort_by(|a, b| {a.weight.partial_cmp(&b.weight).unwrap().reverse()});
}

// The following function is only necessary for the header generation.
#[cfg(feature = "headers")] // c.f. the `Cargo.toml` section
pub fn generate_headers() -> ::std::io::Result<()> {
    ::safer_ffi::headers::builder()
        .to_file("rust_parasort.h")?
        .generate()
}
