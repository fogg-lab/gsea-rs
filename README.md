# gsea-rs

WebAssembly-ready Rust implementation of preranked Gene Set Enrichment Analysis (GSEA) distilled from [GSEApy](https://github.com/zqfang/GSEApy).

## Setup

### Prerequisites
- Rust
- wasm-pack (`cargo install wasm-pack`)

### Installation

1. Clone repo:
```
git clone https://github.com/wigginno/gsea-rs.git
cd gsea-rs
```

2. Build the WebAssembly module:
```
wasm-pack build --target web
```

3. Serve the demo `index.html` file e.g. by running `npx http-server` in the project directory then open `http://localhost:8080` in your browser.

## Demo usage

Options:
- Genes (comma-separated)
- Metric (comma-separated)
- Gene Sets (in JSON format). This could easily be adapted to load .json gene sets from MSigDB.
- Weight
- Min Size
- Max Size
- Number of Permutations
- Random Seed

After filling in the form, click the "Run GSEA" button to perform the analysis. The results will be displayed in a table below the form.

## API

The main function exposed by the WebAssembly module is `prerank_rs`:

```rust
pub fn prerank_rs(
    genes: Vec<String>,
    metric: Vec<f64>,
    gene_sets: Vec<Vec<String>>,
    weight: f64,
    min_size: u32,
    max_size: u32,
    num_permutations: u32,
    random_seed: u32,
) -> Vec<GseaResult>
```

This function takes the following inputs:

- `genes`: A vector of gene names
- `metric`: A vector of metric values corresponding to the genes
- `gene_sets`: A vector of vectors, where each inner vector represents a gene set
- `weight`: The weight parameter for the GSEA algorithm
- `min_size`: The minimum size of the gene sets to consider
- `max_size`: The maximum size of the gene sets to consider
- `num_permutations`: The number of permutations to perform for the GSEA algorithm
- `random_seed`: A random seed for the permutation generation

The function returns a vector of `GseaResult` structs, which contain the following fields:

- `gene_set`: The name of the gene set
- `es`: The enrichment score for the gene set
- `nes`: The normalized enrichment score for the gene set
- `p_val`: The p-value for the gene set
- `fdr`: The false discovery rate for the gene set
