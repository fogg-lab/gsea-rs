use wasm_bindgen::prelude::*;
use std::collections::HashMap;
use serde_wasm_bindgen;

mod algorithm;
mod stats;
mod utils;
use stats::GSEAResult;

#[wasm_bindgen]
pub fn prerank_rs(
    genes: Vec<String>,
    metric: Vec<f64>,
    gene_sets: JsValue,
    weight: f64,
    min_size: usize,
    max_size: usize,
    nperm: usize,
    seed: u64,
) -> Result<JsValue, JsValue> {
    let gene_sets: HashMap<String, Vec<String>> = serde_wasm_bindgen::from_value(gene_sets)?;
    let mut gmt = HashMap::<&str, &[String]>::new();
    for (k, v) in gene_sets.iter() {
        gmt.insert(k.as_str(), v.as_slice());
    }
    let mut gsea = GSEAResult::new(weight, max_size, min_size, nperm, seed);
    gsea.prerank(&genes, &metric, &gmt);
    serde_wasm_bindgen::to_value(&gsea).map_err(|e| JsValue::from_str(&e.to_string()))
}
