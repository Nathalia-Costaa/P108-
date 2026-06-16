// queueing_theory.js
// Sistema completo de filas com cálculos de P(n), P(n>r), P(W>t), P(Wq>t)

// ==================== ESTADO GLOBAL ====================
let currentState = {
    model: 'mm1',
    lambda: 3.0,
    mu: 4.0,
    s: 1,
    K: Infinity,
    Npop: Infinity,
    variance: 0.5,
    rho: 0.75,
    P0: 0.25,
    lambda_efetiva: 3.0,
    isStable: true
};

// ==================== UTILITÁRIOS ====================
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

function formatValue(x) {
    return (typeof x === 'number' && isFinite(x)) ? x.toFixed(6) : '—';
}

// ==================== FUNÇÕES P(n) POR MODELO ====================

function getPn(n) {
    const model = currentState.model;
    const lambda = currentState.lambda;
    const mu = currentState.mu;
    const rho = currentState.rho;
    const P0 = currentState.P0;
    const s = currentState.s;
    const K = currentState.K;
    
    if (n < 0) return 0;
    
    if (model === 'mm1') {
        return P0 * Math.pow(rho, n);
    }
    else if (model === 'mms') {
        const a = lambda / mu;
        if (n <= s) {
            return (Math.pow(a, n) / factorial(n)) * P0;
        } else {
            return (Math.pow(a, n) / (factorial(s) * Math.pow(s, n - s))) * P0;
        }
    }
    else if (model === 'mm1k') {
        if (n > K) return 0;
        if (Math.abs(rho - 1) < 1e-9) return P0;
        return P0 * Math.pow(rho, n);
    }
    else if (model === 'mmsk') {
        if (n > K) return 0;
        const a = lambda / mu;
        if (n <= s) {
            return (Math.pow(a, n) / factorial(n)) * P0;
        } else {
            return (Math.pow(a, n) / (factorial(s) * Math.pow(s, n - s))) * P0;
        }
    }
    else if (model === 'mg1') {
        return P0 * Math.pow(rho, n);
    }
    else if (model === 'finitePop') {
        const N = K;
        if (n > N) return 0;
        return P0 * (factorial(N) / factorial(N - n)) * Math.pow(lambda / mu, n);
    }
    return 0;
}

// ==================== FUNÇÕES P(n > r) ====================

function getPnGr(r) {
    if (r < 0) return 1;
    
    if (currentState.model === 'mm1') {
        return Math.pow(currentState.rho, r + 1);
    }
    
    let sum = 0;
    let maxN = 1000;
    
    for (let n = r + 1; n <= maxN; n++) {
        let pn = getPn(n);
        if (pn < 1e-15 && n > 100) break;
        sum += pn;
        if (sum > 0.999999) break;
    }
    return Math.min(sum, 1);
}

// ==================== FUNÇÕES P(W > t) ====================

function getPWt(t) {
    if (t < 0) return 1;
    const model = currentState.model;
    const mu = currentState.mu;
    const rho = currentState.rho;
    const lambda = currentState.lambda;
    const s = currentState.s;
    const P0 = currentState.P0;
    
    if (model === 'mm1') {
        return Math.exp(-mu * (1 - rho) * t);
    }
    else if (model === 'mms') {
        const a = lambda / mu;
        const term1 = (P0 * Math.pow(a, s)) / (factorial(s) * (1 - rho));
        const denom = (s - 1 - a);
        let secondPart;
        if (Math.abs(denom) < 1e-9) {
            secondPart = term1 * mu * t;
        } else {
            secondPart = term1 * (1 - Math.exp(-mu * t * denom)) / denom;
        }
        return Math.exp(-mu * t) * (1 + secondPart);
    }
    else {
        return Math.exp(-mu * (1 - rho) * t);
    }
}

// ==================== FUNÇÕES P(Wq > t) ====================

function getPWqt(t) {
    if (t < 0) return 1;
    const model = currentState.model;
    const mu = currentState.mu;
    const rho = currentState.rho;
    const lambda = currentState.lambda;
    const s = currentState.s;
    const P0 = currentState.P0;
    const K = currentState.K;
    
    if (model === 'mm1') {
        return rho * Math.exp(-mu * (1 - rho) * t);
    }
    else if (model === 'mms') {
        let sumPn0toSminus1 = 0;
        for (let n = 0; n < s; n++) {
            sumPn0toSminus1 += getPn(n);
        }
        const pWqGtZero = 1 - Math.min(sumPn0toSminus1, 1);
        return pWqGtZero * Math.exp(-s * mu * (1 - rho) * t);
    }
    else if (model === 'mm1k') {
        const PK = getPn(K);
        return PK * Math.exp(-mu * (1 - rho) * t);
    }
    else if (model === 'mmsk') {
        let sumPn0toSminus1 = 0;
        for (let n = 0; n < s; n++) {
            sumPn0toSminus1 += getPn(n);
        }
        const pWqGtZero = 1 - Math.min(sumPn0toSminus1, 1);
        return pWqGtZero * Math.exp(-s * mu * (1 - rho) * t);
    }
    else if (model === 'mg1') {
        const variance = currentState.variance;
        const c_squared = variance * mu * mu;
        const theta = (2 * mu * (1 - rho)) / (1 + c_squared);
        return rho * Math.exp(-theta * t);
    }
    else if (model === 'mmsFinitePop') {
        const N = K; // K será usado como população total
        const s = currentState.s;
        const a = lambda / mu;
        if (n > N) return 0;
        
        if (n <= s) {
            return (factorial(N) / factorial(N - n) / factorial(n)) * Math.pow(a, n) * P0;
        } else {
            return (factorial(N) / factorial(N - n) / factorial(s) / Math.pow(s, n - s)) * Math.pow(a, n) * P0;
        }
    }
    else {
        return rho * Math.exp(-mu * (1 - rho) * t);
    }
}

// ==================== CÁLCULO DAS MÉTRICAS ====================

function updateDynamicFields() {
    const model = document.getElementById('modelType').value;
    const dynamicDiv = document.getElementById('dynamicParams');
    dynamicDiv.innerHTML = '';
    
    if (model === 'mms') {
        dynamicDiv.innerHTML = `<label>🧑‍🤝‍🧑 Número de servidores (s)</label>
                                <input type="number" id="servers_s" step="1" value="2">`;
    } 
    else if (model === 'mm1k') {
        dynamicDiv.innerHTML = `<label>📦 Capacidade máxima K</label>
                                <input type="number" id="capacityK" step="1" value="10">`;
    }
    else if (model === 'mmsk') {
        dynamicDiv.innerHTML = `<div class="param-row">
                                    <label>🧑‍🤝‍🧑 Servidores (s)</label>
                                    <input type="number" id="servers_s_msk" step="1" value="2">
                                    <label>📦 Capacidade K</label>
                                    <input type="number" id="capacityK_msk" step="1" value="8">
                                </div>`;
    }
    else if (model === 'mg1') {
        dynamicDiv.innerHTML = `<label>📉 Variância do tempo de serviço (σ²)</label>
                                <input type="number" id="variance_mg1" step="any" value="0.5">`;
    }
    else if (model === 'finitePop') {
        dynamicDiv.innerHTML = `<label>👥 População total N</label>
                                <input type="number" id="popN" step="1" value="15">`;
    }
    else if (model === 'mmsFinitePop') {
        dynamicDiv.innerHTML = `
            <label>🧑‍🤝 Número de servidores (s)</label>
            <input type="number" id="servers_s_finite" step="1" value="2">
            <label> População total N</label>
            <input type="number" id="popN_finite" step="1" value="15">
        `;
    }
    else if (model === 'mmsFinitePop') {
        extra = `<div class="stat-card"><div class="stat-title">Servidores s / População N</div>
                <div class="stat-value">${s} / ${K}</div></div>
                <div class="stat-card"><div class="stat-title">λ efetiva</div>
                <div class="stat-value">${formatValue(lambda_efetiva)}</div></div>`;
    }
}

function computeAll() {
    const model = document.getElementById('modelType').value;
    let lambda = getNumber('lambda', 3.0);
    let mu = getNumber('mu', 4.0);
    
    if (lambda <= 0) lambda = 0.001;
    if (mu <= 0) mu = 0.001;
    
    currentState.model = model;
    currentState.lambda = lambda;
    currentState.mu = mu;
    
    let rho = lambda / mu;
    let P0 = 1 - rho;
    let L = rho / (1 - rho);
    let Lq = (rho * rho) / (1 - rho);
    let W = 1 / (mu - lambda);
    let Wq = lambda / (mu * (mu - lambda));
    let lambda_efetiva = lambda;
    let s = 1;
    let K = Infinity;
    
    // Calcular conforme modelo
    if (model === 'mm1') {
        if (lambda >= mu) currentState.isStable = false;
        else currentState.isStable = true;
    }
    else if (model === 'mms') {
        s = Math.max(1, Math.floor(getNumber('servers_s', 2)));
        rho = lambda / (s * mu);
        if (rho >= 1) {
            currentState.isStable = false;
        } else {
            currentState.isStable = true;
            let sumP0 = 0;
            for (let n = 0; n < s; n++) {
                sumP0 += Math.pow(lambda/mu, n) / factorial(n);
            }
            const termoS = Math.pow(lambda/mu, s) / factorial(s) * (1 / (1 - rho));
            P0 = 1 / (sumP0 + termoS);
            Lq = (P0 * Math.pow(lambda/mu, s) * rho) / (factorial(s) * (1 - rho) * (1 - rho));
            L = Lq + (lambda / mu);
            W = L / lambda;
            Wq = Lq / lambda;
        }
    }
    else if (model === 'mm1k') {
        K = Math.max(1, Math.floor(getNumber('capacityK', 10)));
        rho = lambda / mu;
        currentState.isStable = true;
        if (Math.abs(rho - 1) < 1e-9) {
            P0 = 1 / (K + 1);
            L = K / 2;
            lambda_efetiva = lambda * (1 - P0);
            Lq = L - (1 - P0);
            W = L / lambda_efetiva;
            Wq = Lq / lambda_efetiva;
        } else {
            const denom = (1 - Math.pow(rho, K + 1));
            P0 = (1 - rho) / denom;
            L = rho / (1 - rho) - ((K + 1) * Math.pow(rho, K + 1)) / denom;
            lambda_efetiva = lambda * (1 - (P0 * Math.pow(rho, K)));
            Lq = L - (lambda_efetiva / mu);
            W = L / lambda_efetiva;
            Wq = Lq / lambda_efetiva;
        }
    }
    else if (model === 'mmsk') {
        s = Math.max(1, Math.floor(getNumber('servers_s_msk', 2)));
        K = Math.max(s, Math.floor(getNumber('capacityK_msk', 8)));
        rho = lambda / (s * mu);
        currentState.isStable = true;
        
        let sumP0 = 0;
        for (let n = 0; n <= s; n++) {
            sumP0 += Math.pow(lambda/mu, n) / factorial(n);
        }
        let secondSum = 0;
        for (let n = s + 1; n <= K; n++) {
            secondSum += Math.pow(lambda/mu, n) / (factorial(s) * Math.pow(s, n - s));
        }
        P0 = 1 / (sumP0 + secondSum);
        
        Lq = 0;
        for (let n = s + 1; n <= K; n++) {
            const pn = (Math.pow(lambda/mu, n) / (factorial(s) * Math.pow(s, n - s))) * P0;
            Lq += (n - s) * pn;
        }
        
        const lastTerm = (Math.pow(lambda/mu, K) / (factorial(s) * Math.pow(s, K - s))) * P0;
        lambda_efetiva = lambda * (1 - lastTerm);
        
        let Lsum = 0;
        for (let n = 0; n <= K; n++) {
            let pn;
            if (n <= s) {
                pn = (Math.pow(lambda/mu, n) / factorial(n)) * P0;
            } else {
                pn = (Math.pow(lambda/mu, n) / (factorial(s) * Math.pow(s, n - s))) * P0;
            }
            Lsum += n * pn;
        }
        L = Lsum;
        W = L / lambda_efetiva;
        Wq = Lq / lambda_efetiva;
    }
    else if (model === 'mg1') {
        const variance = Math.max(0, getNumber('variance_mg1', 0.5));
        currentState.variance = variance;
        rho = lambda / mu;
        if (rho >= 1) {
            currentState.isStable = false;
        } else {
            currentState.isStable = true;
            P0 = 1 - rho;
            Lq = ((lambda * lambda * variance) + (rho * rho)) / (2 * (1 - rho));
            L = rho + Lq;
            W = L / lambda;
            Wq = Lq / lambda;
        }
    }
    else if (model === 'finitePop') {
        const N = Math.max(1, Math.floor(getNumber('popN', 15)));
        K = N;
        rho = lambda / mu;
        currentState.isStable = true;
        
        let sumInv = 0;
        for (let i = 0; i <= N; i++) {
            sumInv += factorial(N) / factorial(N - i) * Math.pow(lambda / mu, i);
        }
        P0 = 1 / sumInv;
        L = N - (mu / lambda) * (1 - P0);
        lambda_efetiva = lambda * (N - L);
        const Pc = P0 * (factorial(N) / factorial(N - N)) * Math.pow(lambda / mu, N);
        Lq = N - ((lambda + mu) / lambda) * (1 - Pc);
        if (isNaN(Lq) || Lq < 0) Lq = Math.max(0, N - (lambda + mu) / lambda * (1 - Pc));
        W = L / lambda_efetiva;
        Wq = Lq / lambda_efetiva;
    }
    else if (model === 'mmsFinitePop') {
        const N = Math.max(1, Math.floor(getNumber('popN_finite', 15)));
        const s = Math.max(1, Math.floor(getNumber('servers_s_finite', 2)));
        const a = lambda / mu;
        
        K = N;
        currentState.s = s;
        currentState.isStable = true; // População finita sempre estável
        
        // --- Cálculo de P₀ ---
        let sumP0 = 0;
        
        // Primeiro somatório: n = 0 até s-1
        for (let n = 0; n < s; n++) {
            sumP0 += (factorial(N) / factorial(N - n) / factorial(n)) * Math.pow(a, n);
        }
        
        // Segundo somatório: n = s até N
        for (let n = s; n <= N; n++) {
            sumP0 += (factorial(N) / factorial(N - n) / factorial(s) / Math.pow(s, n - s)) * Math.pow(a, n);
        }
        
        P0 = 1 / sumP0;
        
        // --- Cálculo de L (número médio no sistema) ---
        let L = 0;
        for (let n = 1; n <= N; n++) {
            let pn;
            if (n <= s) {
                pn = (factorial(N) / factorial(N - n) / factorial(n)) * Math.pow(a, n) * P0;
            } else {
                pn = (factorial(N) / factorial(N - n) / factorial(s) / Math.pow(s, n - s)) * Math.pow(a, n) * P0;
            }
            L += n * pn;
        }
        
        // --- Taxa efetiva de chegada ---
        lambda_efetiva = lambda * (N - L);
        
        // --- Lq (fila) ---
        Lq = L - (a) * (N - L);
        if (Lq < 0) Lq = 0;
        
        // --- W e Wq ---
        W = L / lambda_efetiva;
        Wq = Lq / lambda_efetiva;
}
    
    // Atualizar estado
    currentState.rho = rho;
    currentState.P0 = P0;
    currentState.lambda_efetiva = lambda_efetiva;
    currentState.s = s;
    currentState.K = K;
    
    // Mostrar resultados
    const stable = currentState.isStable && (rho < 1);
    document.getElementById('stabilityBadge').innerHTML = stable ? '✅ Sistema estável' : '⚠️ Instável (ρ ≥ 1)';
    
    if (!stable) {
        document.getElementById('rhoVal').innerHTML = rho.toFixed(4);
        document.getElementById('P0val').innerHTML = '0';
        document.getElementById('Lval').innerHTML = '∞';
        document.getElementById('Lqval').innerHTML = '∞';
        document.getElementById('Wval').innerHTML = '∞';
        document.getElementById('Wqval').innerHTML = '∞';
        document.getElementById('extraMetrics').innerHTML = '<div class="note">⚠️ Sistema instável! Aumente μ ou reduza λ.</div>';
        return;
    }
    
    document.getElementById('rhoVal').innerHTML = formatValue(rho);
    document.getElementById('P0val').innerHTML = formatValue(P0);
    document.getElementById('Lval').innerHTML = formatValue(L);
    document.getElementById('Lqval').innerHTML = formatValue(Lq);
    document.getElementById('Wval').innerHTML = formatValue(W);
    document.getElementById('Wqval').innerHTML = formatValue(Wq);
    
    let extra = '';
    if (model === 'mms') {
        extra = `<div class="stat-card"><div class="stat-title">Servidores s</div><div class="stat-value">${s}</div></div>`;
    } else if (model === 'mm1k') {
        extra = `<div class="stat-card"><div class="stat-title">Capacidade K</div><div class="stat-value">${K}</div></div>
                 <div class="stat-card"><div class="stat-title">λ efetiva</div><div class="stat-value">${formatValue(lambda_efetiva)}</div></div>`;
    } else if (model === 'mmsk') {
        extra = `<div class="stat-card"><div class="stat-title">Servidores s / Capacidade K</div><div class="stat-value">${s} / ${K}</div></div>
                 <div class="stat-card"><div class="stat-title">λ efetiva</div><div class="stat-value">${formatValue(lambda_efetiva)}</div></div>`;
    } else if (model === 'mg1') {
        extra = `<div class="stat-card"><div class="stat-title">Variância σ²</div><div class="stat-value">${currentState.variance.toFixed(4)}</div></div>`;
    } else if (model === 'finitePop') {
        extra = `<div class="stat-card"><div class="stat-title">População N</div><div class="stat-value">${K}</div></div>
                 <div class="stat-card"><div class="stat-title">λ efetiva</div><div class="stat-value">${formatValue(lambda_efetiva)}</div></div>`;
    }
    
    document.getElementById('extraMetrics').innerHTML = extra;
}

// ==================== ATUALIZAR PROBABILIDADES ====================

function updatePn() {
    const n = getNumber('nValue', 2);
    const result = getPn(n);
    document.getElementById('PnResult').innerHTML = formatValue(result);
}

function updatePnGr() {
    const r = getNumber('rValue', 3);
    const result = getPnGr(r);
    document.getElementById('PnGrResult').innerHTML = formatValue(result);
}

function updatePWt() {
    const t = getNumber('tValueW', 0.5);
    const result = getPWt(t);
    document.getElementById('PWtResult').innerHTML = formatValue(result);
}

function updatePWqt() {
    const t = getNumber('tValueWq', 0.5);
    const result = getPWqt(t);
    document.getElementById('PWqtResult').innerHTML = formatValue(result);
}

// ==================== EVENTOS ====================

document.getElementById('modelType').addEventListener('change', () => {
    updateDynamicFields();
    computeAll();
    updatePn();
    updatePnGr();
    updatePWt();
    updatePWqt();
});

document.getElementById('calcBtn').addEventListener('click', () => {
    computeAll();
    updatePn();
    updatePnGr();
    updatePWt();
    updatePWqt();
});

document.getElementById('calcPnBtn').addEventListener('click', updatePn);
document.getElementById('calcPnGrBtn').addEventListener('click', updatePnGr);
document.getElementById('calcPWtBtn').addEventListener('click', updatePWt);
document.getElementById('CalcPWqtBtn').addEventListener('click', updatePWqt);

// Listeners para inputs para atualizar automaticamente
document.getElementById('lambda').addEventListener('input', () => {
    computeAll();
    updatePn();
    updatePnGr();
    updatePWt();
    updatePWqt();
});

document.getElementById('mu').addEventListener('input', () => {
    computeAll();
    updatePn();
    updatePnGr();
    updatePWt();
    updatePWqt();
});

// ==================== INICIALIZAÇÃO ====================
updateDynamicFields();
computeAll();
updatePn();
updatePnGr();
updatePWt();
updatePWqt();