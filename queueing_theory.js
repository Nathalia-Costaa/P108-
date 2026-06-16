// queueing_theory.js
// Utility functions and core calculations for queueing theory models

function getNumber(id, defaultValue) {
    let val = parseFloat(document.getElementById(id).value);
    return isNaN(val) ? defaultValue : val;
}

function factorial(n) {
    if (n < 0) return 1;
    if (n === 0 || n === 1) return 1;
    let res = 1;
    for (let i = 2; i <= n; i++) res *= i;
    return res;
}

function powSafe(base, exp) {
    return Math.pow(base, exp);
}

function formatValue(x) {
    return (typeof x === 'number' && isFinite(x)) ? x.toFixed(6) : '—';
}

function updateDynamicFields() {
    const dynamicDiv = document.getElementById('dynamicParams');
    const model = document.getElementById('modelType').value;
    dynamicDiv.innerHTML = '';
    
    if (model === 'mms') {
        dynamicDiv.innerHTML = `<label>🧑‍🤝‍🧑 Número de servidores (s)</label>
                                <input type="number" id="servers_s" step="1" value="2" placeholder="ex: 2">`;
    } 
    else if (model === 'mm1k') {
        dynamicDiv.innerHTML = `<label>📦 Capacidade máxima do sistema (K)</label>
                                <input type="number" id="capacityK" step="1" value="10" placeholder="máx clientes (inclui serviço)">`;
    }
    else if (model === 'mmsk') {
        dynamicDiv.innerHTML = `<div class="param-row">
                                    <div><label>🧑‍🤝‍🧑 Servidores (s)</label><input type="number" id="servers_s_msk" step="1" value="2"></div>
                                    <div><label>📦 Capacidade K (total)</label><input type="number" id="capacityK_msk" step="1" value="8"></div>
                                </div>`;
    }
    else if (model === 'mg1') {
        dynamicDiv.innerHTML = `<label>📉 Variância do tempo de serviço (σ²)</label>
                                <input type="number" id="variance_mg1" step="any" value="0.5" placeholder="ex: 0.25 ou 0.5">`;
    }
    else if (model === 'finitePop') {
        dynamicDiv.innerHTML = `<label>👥 População total finita (N)</label>
                                <input type="number" id="popN" step="1" value="15" placeholder="N clientes no universo">`;
    }
}

function computeMM1(lambda, mu) {
    const rho = lambda / mu;
    if (rho >= 1) return { stable: false };
    
    const P0 = 1 - rho;
    const L = rho / (1 - rho);
    const Lq = (rho * rho) / (1 - rho);
    const W = 1 / (mu - lambda);
    const Wq = lambda / (mu * (mu - lambda));
    
    return { stable: true, rho, P0, L, Lq, W, Wq };
}

function computeMMS(lambda, mu, s) {
    const rho = lambda / (s * mu);
    if (rho >= 1) return { stable: false, rho };
    
    let sumP0 = 0;
    for (let n = 0; n < s; n++) {
        sumP0 += powSafe(lambda/mu, n) / factorial(n);
    }
    const termoS = powSafe(lambda/mu, s) / factorial(s) * (1 / (1 - rho));
    const P0 = 1 / (sumP0 + termoS);
    const Lq = (P0 * powSafe(lambda/mu, s) * rho) / (factorial(s) * (1 - rho) * (1 - rho));
    const L = Lq + (lambda / mu);
    const W = L / lambda;
    const Wq = Lq / lambda;
    
    return { stable: true, rho, P0, L, Lq, W, Wq };
}

function computeMM1K(lambda, mu, K) {
    const rhoVal = lambda / mu;
    let P0, L, Lq, W, Wq, lambda_bar;
    
    if (Math.abs(rhoVal - 1) < 1e-9) {
        P0 = 1 / (K + 1);
        L = K / 2;
        lambda_bar = lambda * (1 - P0);
        Lq = L - (1 - P0);
        W = L / lambda_bar;
        Wq = Lq / lambda_bar;
    } else {
        const denom = (1 - Math.pow(rhoVal, K + 1));
        P0 = (1 - rhoVal) / denom;
        L = rhoVal / (1 - rhoVal) - ((K + 1) * Math.pow(rhoVal, K + 1)) / denom;
        lambda_bar = lambda * (1 - (P0 * Math.pow(rhoVal, K)));
        Lq = L - (lambda_bar / mu);
        W = L / lambda_bar;
        Wq = Lq / lambda_bar;
    }
    
    return { stable: true, rho: rhoVal, P0, L, Lq, W, Wq, lambda_bar, K };
}

function computeMMSK(lambda, mu, s, K) {
    let sumP0 = 0;
    for (let n = 0; n <= s; n++) {
        sumP0 += powSafe(lambda/mu, n) / factorial(n);
    }
    
    let secondSum = 0;
    for (let n = s + 1; n <= K; n++) {
        secondSum += powSafe(lambda/mu, n) / (factorial(s) * Math.pow(s, n - s));
    }
    
    const P0 = 1 / (sumP0 + secondSum);
    
    let Lq = 0;
    for (let n = s + 1; n <= K; n++) {
        const pn = (powSafe(lambda/mu, n) / (factorial(s) * Math.pow(s, n - s))) * P0;
        Lq += (n - s) * pn;
    }
    
    const lastTerm = (powSafe(lambda/mu, K) / (factorial(s) * Math.pow(s, K - s))) * P0;
    const lambda_bar = lambda * (1 - lastTerm);
    
    let Lsum = 0;
    for (let n = 0; n <= K; n++) {
        let pn;
        if (n <= s) {
            pn = (powSafe(lambda/mu, n) / factorial(n)) * P0;
        } else {
            pn = (powSafe(lambda/mu, n) / (factorial(s) * Math.pow(s, n - s))) * P0;
        }
        Lsum += n * pn;
    }
    
    const L = Lsum;
    const W = L / lambda_bar;
    const Wq = Lq / lambda_bar;
    const rho = lambda / (s * mu);
    
    return { stable: true, rho, P0, L, Lq, W, Wq, lambda_bar, K, s };
}

function computeMG1(lambda, mu, variance) {
    const rho = lambda / mu;
    if (rho >= 1) return { stable: false, rho };
    
    const P0 = 1 - rho;
    const Lq = ((lambda * lambda * variance) + (rho * rho)) / (2 * (1 - rho));
    const L = rho + Lq;
    const W = L / lambda;
    const Wq = Lq / lambda;
    
    return { stable: true, rho, P0, L, Lq, W, Wq };
}

function computeFinitePopulation(lambda, mu, N) {
    let sumInv = 0;
    for (let i = 0; i <= N; i++) {
        sumInv += factorial(N) / factorial(N - i) * Math.pow(lambda/mu, i);
    }
    
    const P0 = 1 / sumInv;
    const L = N - (mu/lambda) * (1 - P0);
    const lambda_bar = lambda * (N - L);
    const Pc = P0 * (factorial(N) / factorial(N - N)) * Math.pow(lambda/mu, N);
    let Lq = N - ((lambda + mu)/lambda) * (1 - Pc);
    if (isNaN(Lq) || Lq < 0) Lq = Math.max(0, N - (lambda+mu)/lambda * (1 - Pc));
    
    const W = L / lambda_bar;
    const Wq = Lq / lambda_bar;
    const rho = lambda / mu;
    
    return { stable: true, rho, P0, L, Lq, W, Wq, lambda_bar, N, Pc };
}

function displayResults(results, model, extraParams = {}) {
    document.getElementById('stabilityBadge').innerHTML = results.stable ? '✅ Sistema estável' : '⚠️ Instável (ρ ≥ 1)';
    
    if (!results.stable) {
        document.getElementById('rhoVal').innerHTML = results.rho.toFixed(4);
        document.getElementById('Lval').innerHTML = '∞ (instável)';
        document.getElementById('Lqval').innerHTML = '∞';
        document.getElementById('Wval').innerHTML = '∞';
        document.getElementById('Wqval').innerHTML = '∞';
        document.getElementById('P0val').innerHTML = '0 (saturado)';
        document.getElementById('extraMetrics').innerHTML = '<div class="note">⚠️ Aumente μ ou reduza λ para estabilidade.</div>';
        return;
    }
    
    document.getElementById('rhoVal').innerHTML = formatValue(results.rho);
    document.getElementById('Lval').innerHTML = formatValue(results.L);
    document.getElementById('Lqval').innerHTML = formatValue(results.Lq);
    document.getElementById('Wval').innerHTML = formatValue(results.W);
    document.getElementById('Wqval').innerHTML = formatValue(results.Wq);
    document.getElementById('P0val').innerHTML = formatValue(results.P0);
    
    let extraHtml = '';
    
    if (model === 'mmsk') {
        extraHtml = `<div class="stat-card"><div class="stat-title">Capacidade K / s</div><div class="stat-value">K=${extraParams.K}, s=${extraParams.s}</div></div>
                     <div class="stat-card"><div class="stat-title">Taxa efetiva λ̅</div><div class="stat-value">${formatValue(results.lambda_bar)}</div></div>`;
    } else if (model === 'mm1k') {
        extraHtml = `<div class="stat-card"><div class="stat-title">Capacidade máxima K</div><div class="stat-value">${extraParams.K}</div></div>
                     <div class="stat-card"><div class="stat-title">λ efetiva</div><div class="stat-value">${formatValue(results.lambda_bar)}</div></div>`;
    } else if (model === 'mg1') {
        extraHtml = `<div class="stat-card"><div class="stat-title">Variância σ²</div><div class="stat-value">${extraParams.variance.toFixed(4)}</div></div>
                     <div class="stat-card"><div class="stat-title">Fórmula Pollaczek-Khintchine</div><div class="stat-value">✓ M/G/1</div></div>`;
    } else if (model === 'mms') {
        extraHtml = `<div class="stat-card"><div class="stat-title">Servidores s</div><div class="stat-value">${extraParams.s}</div></div>`;
    } else if (model === 'finitePop') {
        extraHtml = `<div class="stat-card"><div class="stat-title">População finita N</div><div class="stat-value">${extraParams.N}</div></div>
                     <div class="stat-card"><div class="stat-title">λ̅ (taxa efetiva)</div><div class="stat-value">${formatValue(results.lambda_bar)}</div></div>`;
    } else {
        extraHtml = `<div class="stat-card"><div class="stat-title">Regime M/M/1</div><div class="stat-value">λ < μ</div></div>`;
    }
    
    document.getElementById('extraMetrics').innerHTML = extraHtml;
}

function computeAll() {
    const model = document.getElementById('modelType').value;
    let lambda = getNumber('lambda', 1.0);
    let mu = getNumber('mu', 1.0);
    
    if (lambda <= 0) lambda = 0.001;
    if (mu <= 0) mu = 0.001;
    
    let results;
    
    switch(model) {
        case 'mm1':
            results = computeMM1(lambda, mu);
            displayResults(results, model);
            break;
            
        case 'mms':
            const s_mms = Math.max(1, Math.floor(getNumber('servers_s', 2)));
            results = computeMMS(lambda, mu, s_mms);
            displayResults(results, model, { s: s_mms });
            break;
            
        case 'mm1k':
            const K_mm1k = Math.max(1, Math.floor(getNumber('capacityK', 10)));
            results = computeMM1K(lambda, mu, K_mm1k);
            displayResults(results, model, { K: K_mm1k });
            break;
            
        case 'mmsk':
            const s_msk = Math.max(1, Math.floor(getNumber('servers_s_msk', 2)));
            const K_msk = Math.max(s_msk, Math.floor(getNumber('capacityK_msk', 8)));
            results = computeMMSK(lambda, mu, s_msk, K_msk);
            displayResults(results, model, { s: s_msk, K: K_msk });
            break;
            
        case 'mg1':
            const variance = Math.max(0, getNumber('variance_mg1', 0.5));
            results = computeMG1(lambda, mu, variance);
            displayResults(results, model, { variance: variance });
            break;
            
        case 'finitePop':
            const N_pop = Math.max(1, Math.floor(getNumber('popN', 15)));
            results = computeFinitePopulation(lambda, mu, N_pop);
            displayResults(results, model, { N: N_pop });
            break;
            
        default:
            results = computeMM1(lambda, mu);
            displayResults(results, model);
    }
}

// Event Listeners and Initialization
document.getElementById('modelType').addEventListener('change', () => {
    updateDynamicFields();
    computeAll();
});

document.getElementById('calcBtn').addEventListener('click', computeAll);

// Add listeners to all inputs for real-time calculation
const allInputs = document.querySelectorAll('input, select');
allInputs.forEach(el => {
    el.addEventListener('change', computeAll);
    el.addEventListener('input', computeAll);
});

// Initialize
updateDynamicFields();
computeAll();