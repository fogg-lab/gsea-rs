async function init_wasm_in_worker() {
  try {
    const wasm = await import('./pkg/gsea_rs.js');
    await wasm.default();

    self.onmessage = async function (event) {
      const { genes, metric, geneSets, weight, minSize, maxSize, nperm, seed } = event.data;

      try {
        const result = wasm.prerank_rs(
          genes,
          metric,
          geneSets,
          weight,
          minSize,
          maxSize,
          nperm,
          BigInt(seed)
        );
        self.postMessage({ status: 'success', result });
      } catch (error) {
        self.postMessage({ status: 'error', error: error.toString() });
      }
    };

    self.postMessage({ status: 'ready' });
  } catch (error) {
    self.postMessage({ status: 'error', error: `Failed to initialize WebAssembly: ${error.toString()}` });
  }
}

init_wasm_in_worker().catch(err => {
  self.postMessage({ status: 'error', error: `Failed to initialize WebAssembly: ${err.toString()}` });
});
