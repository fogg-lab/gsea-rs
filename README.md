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
    gene_sets: JsValue,
    weight: f64,
    min_size: usize,
    max_size: usize,
    nperm: usize,
    seed: u64,
) -> Result<JsValue, JsValue>
```

This function takes the following inputs:

- `genes`: A vector of gene names (as strings)
- `metric`: A vector of metric values corresponding to the genes (as floating-point numbers)
- `gene_sets`: A JavaScript object representing gene sets, where keys are set names and values are arrays of gene names
- `weight`: The weight parameter for the GSEA algorithm (as a floating-point number)
- `min_size`: The minimum size of the gene sets to consider (as an integer)
- `max_size`: The maximum size of the gene sets to consider (as an integer)
- `nperm`: The number of permutations to perform for the GSEA algorithm (as an integer)
- `seed`: A random seed for the permutation generation (as an unsigned 64-bit integer)

The function returns a `Result` containing a `JsValue` on success, which represents an array of summary objects with the following structure:

```javascript
[
  {
    term: string,
    es: number,
    nes: number,
    pval: number,
    fwerp: number,
    fdr: number,
    run_es: number[],
    hits: number[],
    esnull: number[],
    index: number | null
  },
  // ... more summary objects
]
```

Each summary object in the array contains:

- `term`: The name of the gene set
- `es`: The enrichment score for the gene set
- `nes`: The normalized enrichment score for the gene set
- `pval`: The p-value for the gene set
- `fwerp`: The family-wise error rate p-value for the gene set
- `fdr`: The false discovery rate for the gene set
- `run_es`: An array of running enrichment scores
- `hits`: An array of indices representing the positions of hits in the ranked list
- `esnull`: An array of null enrichment scores
- `index`: An optional index value (may be null)

To use this function in JavaScript, you'll need to use a Web Worker to run the WebAssembly module. See `gsea-worker.js` and `gsea-app.js` for implementation details.
