let worker;

function initWorker() {
    worker = new Worker('gsea-worker.js');
    worker.onmessage = handleWorkerMessage;
    worker.onerror = handleWorkerError;
}

function handleWorkerMessage(event) {
    console.log('Worker message received:', event.data);
    const { status, result, error } = event.data;

    if (status === 'ready') {
        console.log('Worker is ready');
        document.getElementById('run-button').disabled = false;
    } else if (status === 'success') {
        displayResults(result);
    } else if (status === 'error') {
        displayError(error);
    }
}

function handleWorkerError(error) {
    console.error('Worker error:', error);
    displayError('An error occurred in the worker: ' + error.message);
}

function runGSEA(event) {
    event.preventDefault();

    const genes = document.getElementById('genes').value.split(',').map(g => g.trim());
    const metric = document.getElementById('metric').value.split(',').map(m => parseFloat(m.trim()));
    const geneSets = JSON.parse(document.getElementById('gene-sets').value);
    const weight = parseFloat(document.getElementById('weight').value);
    const minSize = parseInt(document.getElementById('min-size').value);
    const maxSize = parseInt(document.getElementById('max-size').value);
    const nperm = parseInt(document.getElementById('nperm').value);
    const seed = parseInt(document.getElementById('seed').value);

    worker.postMessage({
        genes,
        metric,
        geneSets,
        weight,
        minSize,
        maxSize,
        nperm,
        seed
    });

    document.getElementById('run-button').disabled = true;
    document.getElementById('results').innerHTML = 'Calculating...';
}

function displayResults(result) {
    console.log('Result received:', result);
    const resultsDiv = document.getElementById('results');
    resultsDiv.innerHTML = '<h2>Results:</h2>';

    if (!result || typeof result !== 'object') {
        displayError('Invalid result structure');
        return;
    }

    const summaries = result.summaries;

    if (!Array.isArray(summaries)) {
        displayError('Invalid summaries structure');
        return;
    }

    const table = document.createElement('table');
    table.innerHTML = `
        <tr>
            <th>Term</th>
            <th>ES</th>
            <th>NES</th>
            <th>P-value</th>
            <th>FDR</th>
        </tr>
    `;

    summaries.forEach(summary => {
        const row = document.createElement('tr');
        row.innerHTML = `
            <td>${summary.term}</td>
            <td>${summary.es.toFixed(4)}</td>
            <td>${summary.nes.toFixed(4)}</td>
            <td>${summary.pval.toExponential(2)}</td>
            <td>${summary.fdr.toExponential(2)}</td>
        `;
        table.appendChild(row);
    });

    resultsDiv.appendChild(table);
    document.getElementById('run-button').disabled = false;
}

function displayError(error) {
    const resultsDiv = document.getElementById('results');
    resultsDiv.innerHTML = `<h2>Error:</h2><p>${error}</p>`;
    document.getElementById('run-button').disabled = false;
}

document.addEventListener('DOMContentLoaded', () => {
    initWorker();
    document.getElementById('gsea-form').addEventListener('submit', runGSEA);
});
