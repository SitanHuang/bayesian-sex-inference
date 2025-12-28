

const $ = (id) => document.getElementById(id);
const fmt = (x, d = 3) => (Number.isFinite(x) ? x.toFixed(d) : "—");
const clamp01 = (x) => Math.max(0, Math.min(1, x));


function erf(x) {
  const sign = x < 0 ? -1 : 1;
  x = Math.abs(x);
  const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;
  const p = 0.3275911;
  const t = 1 / (1 + p * x);
  const y = 1 - ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
  return sign * y;
}
function normalCdf(z) { return 0.5 * (1 + erf(z / Math.SQRT2)); }
function logNormalPdf(x, mu, sigma) {
  if (!(sigma > 0)) return -Infinity;
  const z = (x - mu) / sigma;
  return -0.5 * Math.log(2 * Math.PI) - Math.log(sigma) - 0.5 * z * z;
}
function logSumExp(a, b) {
  const m = Math.max(a, b);
  return m + Math.log(Math.exp(a - m) + Math.exp(b - m));
}

// Lanczos gammaln + regularized lower incomplete gamma P(s,x)
function gammaln(xx) {
  const cof = [
    76.18009172947146, -86.50532032941677,
    24.01409824083091, -1.231739572450155,
    0.001208650973866179, -0.000005395239384953
  ];
  let x = xx - 1.0;
  let tmp = x + 5.5;
  tmp -= (x + 0.5) * Math.log(tmp);
  let ser = 1.000000000190015;
  for (let j = 0; j < cof.length; j++) { x += 1; ser += cof[j] / x; }
  return -tmp + Math.log(2.5066282746310005 * ser);
}
function gammaP(s, x) {
  if (x < 0 || s <= 0) return NaN;
  if (x === 0) return 0;
  const gln = gammaln(s);
  if (x < s + 1) {
    let ap = s;
    let sum = 1 / s;
    let del = sum;
    for (let n = 1; n <= 200; n++) {
      ap += 1;
      del *= x / ap;
      sum += del;
      if (Math.abs(del) < Math.abs(sum) * 1e-14) break;
    }
    return sum * Math.exp(-x + s * Math.log(x) - gln);
  } else {
    let b = x + 1 - s;
    let c = 1 / 1e-30;
    let d = 1 / b;
    let h = d;
    for (let i = 1; i <= 200; i++) {
      const an = -i * (i - s);
      b += 2;
      d = an * d + b;
      if (Math.abs(d) < 1e-30) d = 1e-30;
      c = b + an / c;
      if (Math.abs(c) < 1e-30) c = 1e-30;
      d = 1 / d;
      const del = d * c;
      h *= del;
      if (Math.abs(del - 1) < 1e-14) break;
    }
    return 1 - h * Math.exp(-x + s * Math.log(x) - gln);
  }
}
function chiSquareCdf(x, k) {
  if (!(x >= 0) || !(k > 0)) return NaN;
  return gammaP(k / 2, x / 2);
}


const LengthUnits = { mm: 1, cm: 10, m: 1000, in: 25.4, ft: 304.8 };
const WeightUnits = { kg: 1, lb: 0.45359237 };

function parseFeetInches(s) {
  const str = String(s).trim();
  const m = str.match(/^\s*(\d+(?:\.\d+)?)\s*(?:ft|')\s*(\d+(?:\.\d+)?)?\s*(?:in|")?\s*$/i);
  if (!m) return null;
  const ft = parseFloat(m[1]);
  const inch = m[2] ? parseFloat(m[2]) : 0;
  if (!Number.isFinite(ft) || !Number.isFinite(inch)) return null;
  return ft * 12 + inch;
}

function parseNumberWithUnit(input, fallbackUnit, dim) {
  const raw = String(input ?? "").trim();
  if (!raw) return { ok: false };

  if (dim === "length") {
    const fi = parseFeetInches(raw);
    if (fi != null) return { ok: true, value: fi, unit: "in" };
  }

  const m = raw.match(/^\s*([-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?)\s*([a-z"']+)?\s*$/i);
  if (!m) return { ok: false };
  const value = parseFloat(m[1]);
  if (!Number.isFinite(value)) return { ok: false };
  let unit = (m[2] || "").trim().toLowerCase();
  if (!unit) unit = fallbackUnit;

  unit = unit.replace(/inches|inch/g, "in");
  unit = unit.replace(/lbs|pounds|pound/g, "lb");
  unit = unit.replace(/meters|meter/g, "m");
  unit = unit.replace(/centimeters|centimeter/g, "cm");
  unit = unit.replace(/millimeters|millimeter/g, "mm");
  unit = unit.replace(/feet|foot/g, "ft");
  unit = unit.replace(/"/g, "in").replace(/'/g, "ft");

  if (dim === "length") {
    if (!(unit in LengthUnits)) return { ok: false };
    return { ok: true, value, unit };
  }
  if (dim === "weight") {
    if (!(unit in WeightUnits)) return { ok: false };
    return { ok: true, value, unit };
  }
  return { ok: true, value, unit: "raw" };
}

function convertToNative(value, fromUnit, nativeUnit, dim) {
  if (!Number.isFinite(value)) return NaN;
  if (dim === "length") {
    const mm = value * LengthUnits[fromUnit];
    return mm / LengthUnits[nativeUnit];
  }
  if (dim === "weight") {
    const kg = value * WeightUnits[fromUnit];
    return kg / WeightUnits[nativeUnit];
  }
  return value;
}

function inferDimAndNativeUnit(colName, sampleMedian) {
  const name = String(colName).toLowerCase();

  if ((name.includes("weight") ||
    name.endsWith("kg") || name.includes("lbs") || name.includes("lb")) &&
    !name.includes("breadth") && !name.includes("length") && !name.includes("height")) {
    const native = name.includes("lb") || name.includes("lbs") ? "lb"
      : name.includes("kg") ? "kg"
        : (sampleMedian > 140 ? "lb" : "kg");
    return { dim: "weight", nativeUnit: native, unitChoices: ["kg", "lb"] };
  }

  if (name === "age" || name.includes("age")) return { dim: "raw", nativeUnit: "raw", unitChoices: ["raw"] };
  if (name.endsWith("in") && sampleMedian > 10 && sampleMedian < 120) {
    return { dim: "length", nativeUnit: "in", unitChoices: ["mm", "cm", "m", "in", "ft"] };
  }
  if (name.includes("inch") || name.includes("inches")) {
    return { dim: "length", nativeUnit: "in", unitChoices: ["mm", "cm", "m", "in", "ft"] };
  }
  if (name.endsWith("cm") || name.includes("centimeter")) {
    return { dim: "length", nativeUnit: "cm", unitChoices: ["mm", "cm", "m", "in", "ft"] };
  }
  if (name.endsWith("mm") || name.includes("millimeter")) {
    return { dim: "length", nativeUnit: "mm", unitChoices: ["mm", "cm", "m", "in", "ft"] };
  }
  if (name.endsWith("m") && sampleMedian > 0.4 && sampleMedian < 3) {
    return { dim: "length", nativeUnit: "m", unitChoices: ["mm", "cm", "m", "in", "ft"] };
  }

  if (sampleMedian > 0.3 && sampleMedian < 3) {
    return { dim: "length", nativeUnit: "m", unitChoices: ["mm", "cm", "m", "in", "ft"] };
  }
  if (sampleMedian >= 3 && sampleMedian < 120) {
    if (sampleMedian > 45 && sampleMedian < 90 && (name.includes("height") || name.includes("stature"))) {
      return { dim: "length", nativeUnit: "in", unitChoices: ["mm", "cm", "m", "in", "ft"] };
    }
    return { dim: "length", nativeUnit: "cm", unitChoices: ["mm", "cm", "m", "in", "ft"] };
  }
  return { dim: "length", nativeUnit: "mm", unitChoices: ["mm", "cm", "m", "in", "ft"] };
}

function defaultInputUnit(meta) {
  if (!meta) return "raw";
  if (meta.dim === "length") {
    // Default to cm unless scale strongly suggests mm
    if (meta.nativeUnit === "mm" && meta.median > 80) return "mm";
    return "cm";
  }
  if (meta.dim === "weight") return "kg";
  return meta.nativeUnit || "raw";
}

function parseCSVAsync(text, onProgress) {
  return new Promise((resolve, reject) => {
    if (typeof text !== "string") return reject(new Error("CSV text is not a string."));
    const rows = [];
    let row = [];
    let field = "";
    let inQuotes = false;

    let i = 0;
    const n = text.length;
    const CHUNK = 140000;

    function pushField() { row.push(field); field = ""; }
    function pushRow() { rows.push(row); row = []; }

    function step() {
      const end = Math.min(i + CHUNK, n);
      for (; i < end; i++) {
        const c = text[i];
        if (inQuotes) {
          if (c === '"') {
            const next = text[i + 1];
            if (next === '"') { field += '"'; i++; }
            else inQuotes = false;
          } else field += c;
        } else {
          if (c === '"') inQuotes = true;
          else if (c === ",") pushField();
          else if (c === "\n") { pushField(); pushRow(); }
          else if (c === "\r") { /* ignore */ }
          else field += c;
        }
      }
      if (onProgress) onProgress(i / n);
      if (i < n) setTimeout(step, 0);
      else {
        pushField();
        if (row.length > 1 || (row.length === 1 && row[0] !== "")) pushRow();
        resolve(rows);
      }
    }
    step();
  });
}

class Dataset {
  constructor(label, headers, rows) {
    this.label = label;
    this.headers = headers;
    this.rows = rows;
    this.colIndex = new Map(headers.map((h, idx) => [h, idx]));
    this.numericCols = [];
    this.colMeta = new Map();
    this._buildNumericMeta();
  }

  _buildNumericMeta() {
    const H = this.headers;
    const R = this.rows;
    const sampleN = Math.min(250, R.length);
    const colStats = H.map(() => ({ ok: 0, bad: 0, vals: [] }));
    for (let r = 0; r < sampleN; r++) {
      const row = R[r];
      for (let c = 0; c < H.length; c++) {
        const v = row[c];
        const x = v === "" || v == null ? NaN : Number(v);
        if (Number.isFinite(x)) { colStats[c].ok++; colStats[c].vals.push(x); }
        else colStats[c].bad++;
      }
    }

    for (let c = 0; c < H.length; c++) {
      const s = colStats[c];
      const total = s.ok + s.bad;
      if (total === 0) continue;
      if (s.ok / total < 0.80) continue;
      if (s.vals.length < 20) continue;

      const vals = s.vals.slice().sort((a, b) => a - b);
      const median = vals[Math.floor(vals.length / 2)];
      const min = vals[0], max = vals[vals.length - 1];

      this.numericCols.push(H[c]);
      const { dim, nativeUnit, unitChoices } = inferDimAndNativeUnit(H[c], median);
      this.colMeta.set(H[c], { dim, nativeUnit, unitChoices, median, min, max });
    }

    const set = new Set(this.numericCols);
    this.numericCols = this.headers.filter(h => set.has(h));
  }

  getMeta(col) { return this.colMeta.get(col); }

  getColumnValues(col) {
    const idx = this.colIndex.get(col);
    const out = new Float64Array(this.rows.length);
    for (let i = 0; i < this.rows.length; i++) {
      const v = this.rows[i][idx];
      const x = v === "" || v == null ? NaN : Number(v);
      out[i] = Number.isFinite(x) ? x : NaN;
    }
    return out;
  }

  getStringValue(rowIdx, col) {
    const idx = this.colIndex.get(col);
    if (idx == null) return "";
    const v = this.rows[rowIdx][idx];
    return (v == null) ? "" : String(v);
  }

  defaultConditioningColumn() {
    const lc = (s) => String(s).toLowerCase();
    const prefer = this.numericCols.filter(c => {
      const n = lc(c);
      return n === "stature" || n.includes("stature") || n.includes("height");
    });
    if (prefer.length) {
      const exact = prefer.find(c => lc(c) === "stature");
      return exact || prefer.sort((a, b) => a.length - b.length)[0];
    }
    return this.numericCols[0] || null;
  }

  suggestBandwidth(col) {
    const idx = this.colIndex.get(col);
    const vals = [];
    const sampleN = Math.min(1200, this.rows.length);
    for (let i = 0; i < sampleN; i++) {
      const v = this.rows[i][idx];
      const x = v === "" || v == null ? NaN : Number(v);
      if (Number.isFinite(x)) vals.push(x);
    }
    if (vals.length < 50) return 1;
    vals.sort((a, b) => a - b);
    const q = (p) => vals[Math.floor(p * (vals.length - 1))];
    const iqr = q(0.75) - q(0.25);
    const bw = (iqr > 0 ? 0.20 * iqr : 0.08 * (q(0.95) - q(0.05)));
    return bw > 0 ? bw : 1;
  }

  suggestBandwidthMasked(col, mask) {
    const idx = this.colIndex.get(col);
    const vals = [];
    const cap = 4000;

    for (let i = 0; i < this.rows.length; i++) {
      if (mask && mask[i] === 0) continue;
      const v = this.rows[i][idx];
      const x = v === "" || v == null ? NaN : Number(v);
      if (Number.isFinite(x)) {
        vals.push(x);
        if (vals.length >= cap) break;
      }
    }

    if (vals.length < 50) return this.suggestBandwidth(col);

    vals.sort((a, b) => a - b);
    const q = (p) => vals[Math.floor(p * (vals.length - 1))];
    const iqr = q(0.75) - q(0.25);
    const bw = (iqr > 0 ? 0.20 * iqr : 0.08 * (q(0.95) - q(0.05)));
    return bw > 0 ? bw : this.suggestBandwidth(col);
  }

}

function gaussianKernel(u) { return Math.exp(-0.5 * u * u); }

function computeWeights(condVals, condValue, bw, mask) {
  const n = condVals.length;
  const w = new Float64Array(n);
  let sw = 0, sw2 = 0;
  if (!(bw > 0)) bw = 1;
  for (let i = 0; i < n; i++) {
    if (mask && mask[i] === 0) { w[i] = 0; continue; }
    const x = condVals[i];
    if (!Number.isFinite(x)) { w[i] = 0; continue; }
    const u = (x - condValue) / bw;
    const wi = gaussianKernel(u);
    w[i] = wi;
    sw += wi;
    sw2 += wi * wi;
  }
  const nEff = sw2 > 0 ? (sw * sw) / sw2 : 0;
  return { w, sw, sw2, nEff };
}

function weightedMeanStd(targetVals, weights) {
  let sw = 0, swx = 0, swxx = 0, count = 0;
  for (let i = 0; i < targetVals.length; i++) {
    const x = targetVals[i];
    const wi = weights[i];
    if (!(wi > 0) || !Number.isFinite(x)) continue;
    sw += wi;
    swx += wi * x;
    swxx += wi * x * x;
    count++;
  }
  if (!(sw > 0) || count < 8) return { ok: false, mean: NaN, sd: NaN, count, sw };
  const mean = swx / sw;
  const varx = Math.max(0, (swxx / sw) - mean * mean);
  const sd = Math.sqrt(varx);
  return { ok: true, mean, sd, count, sw };
}

// Filter expression: Field ~= /regex/flags with AND/OR/()
function tokenizeFilterExpr(expr) {
  const s = String(expr || "");
  const tokens = [];
  let i = 0;

  const isSpace = (c) => /\s/.test(c);
  const matchWord = (w) => s.slice(i, i + w.length).toUpperCase() === w;
  while (i < s.length) {
    if (isSpace(s[i])) { i++; continue; }
    const c = s[i];
    if (c === "(" || c === ")") { tokens.push({ t: c }); i++; continue; }

    if (matchWord("AND") && (i + 3 === s.length || /\W/.test(s[i + 3] || ""))) { tokens.push({ t: "AND" }); i += 3; continue; }
    if (matchWord("OR") && (i + 2 === s.length || /\W/.test(s[i + 2] || ""))) { tokens.push({ t: "OR" }); i += 2; continue; }

    let j = i;
    let depth = 0;
    while (j < s.length) {
      const cj = s[j];
      if (cj === "(") depth++;
      if (cj === ")") { if (depth === 0) break; depth--; }
      if (depth === 0) {
        if (s.slice(j).match(/^\s+AND\b/i)) break;
        if (s.slice(j).match(/^\s+OR\b/i)) break;
      }
      j++;
    }
    const clauseText = s.slice(i, j).trim();
    if (clauseText) tokens.push({ t: "CLAUSE", v: clauseText });
    i = j;
  }
  return tokens;
}

function parseClause(clauseText) {
  const m = clauseText.match(/^\s*([A-Za-z0-9_]+)\s*~=\s*(.+?)\s*$/);
  if (!m) throw new Error(`Invalid clause: ${clauseText}`);
  const field = m[1];
  let rhs = m[2].trim();

  let re;
  if (rhs.startsWith("/")) {
    const lastSlash = rhs.lastIndexOf("/");
    if (lastSlash <= 0) throw new Error(`Invalid regex literal: ${rhs}`);
    const pat = rhs.slice(1, lastSlash);
    const flags = rhs.slice(lastSlash + 1) || "";
    re = new RegExp(pat, flags);
  } else {
    re = new RegExp(rhs, "i");
  }
  return { field, re };
}

function toRpn(tokens) {
  const out = [];
  const op = [];
  const prec = { "OR": 1, "AND": 2 };
  for (const tk of tokens) {
    if (tk.t === "CLAUSE") out.push(tk);
    else if (tk.t === "AND" || tk.t === "OR") {
      while (op.length) {
        const top = op[op.length - 1];
        if ((top.t === "AND" || top.t === "OR") && prec[top.t] >= prec[tk.t]) out.push(op.pop());
        else break;
      }
      op.push(tk);
    } else if (tk.t === "(") op.push(tk);
    else if (tk.t === ")") {
      while (op.length && op[op.length - 1].t !== "(") out.push(op.pop());
      if (!op.length) throw new Error("Mismatched parentheses.");
      op.pop();
    }
  }
  while (op.length) {
    const top = op.pop();
    if (top.t === "(" || top.t === ")") throw new Error("Mismatched parentheses.");
    out.push(top);
  }
  return out;
}

function buildFilterFunction(expr, dsExampleForValidation) {
  const text = String(expr || "").trim();
  if (!text) return { fn: null, msg: "No row filter applied." };

  const tokens = tokenizeFilterExpr(text);
  const rpn = toRpn(tokens);

  const compiled = rpn.map(tk => {
    if (tk.t !== "CLAUSE") return tk;
    const { field, re } = parseClause(tk.v);
    if (dsExampleForValidation && !dsExampleForValidation.colIndex.has(field)) {
      throw new Error(`Unknown field in filter: ${field}`);
    }
    return { t: "CLAUSE", field, re };
  });
  const fn = (ds, rowIdx) => {
    const st = [];
    for (const tk of compiled) {
      if (tk.t === "CLAUSE") {
        const v = ds.getStringValue(rowIdx, tk.field);
        st.push(tk.re.test(v));
      } else if (tk.t === "AND") {
        const b = st.pop(), a = st.pop();
        st.push(Boolean(a && b));
      } else if (tk.t === "OR") {
        const b = st.pop(), a = st.pop();
        st.push(Boolean(a || b));
      }
    }
    return st.length ? Boolean(st[0]) : true;
  };

  return { fn, msg: "Filter parsed." };
}

function buildRowMask(ds, filterFn) {
  const n = ds.rows.length;
  const mask = new Uint8Array(n);
  if (!filterFn) { mask.fill(1); return mask; }
  for (let i = 0; i < n; i++) mask[i] = filterFn(ds, i) ? 1 : 0;
  return mask;
}

// State
let dsF = null, dsM = null;
const selected = new Map(); // field -> { unitChoice, valueText, uncText, weight }
let cache = { condCol: null, condValsF: null, condValsM: null };

// UI wiring
function setEnabled(enabled) {
  $("selCond").disabled = !enabled;
  $("inpCondVal").disabled = !enabled;
  $("selCondUnit").disabled = !enabled;
  $("inpBW").disabled = !enabled;
  $("inpPriorM").disabled = !enabled;
  $("txtFilter").disabled = !enabled;
  $("inpSearch").disabled = !enabled;
  $("selFields").disabled = !enabled;
  $("btnAddField").disabled = !enabled;
  $("btnClearFields").disabled = !enabled;
  $("btnRun").disabled = !enabled;
  $("btnSave").disabled = !enabled;
  $("btnLoadSaved").disabled = !enabled;
}

function intersectNumericCols(a, b) {
  const setB = new Set(b.numericCols);
  return a.numericCols.filter(x => setB.has(x));
}

function pickDefaultConditioning(commonNumeric) {
  const lc = s => s.toLowerCase();
  const exact = commonNumeric.find(c => lc(c) === "stature");
  if (exact) return exact;
  const statureLike = commonNumeric.find(c => lc(c).includes("stature"));
  if (statureLike) return statureLike;
  const heightLike = commonNumeric.find(c => lc(c).includes("height"));
  if (heightLike) return heightLike;
  return commonNumeric[0] || null;
}

function rebuildConditioningUI() {
  const commonNumeric = intersectNumericCols(dsF, dsM);
  $("selCond").innerHTML = "";
  for (const col of commonNumeric) {
    const opt = document.createElement("option");
    opt.value = col; opt.textContent = col;
    $("selCond").appendChild(opt);
  }

  const defaultCond = pickDefaultConditioning(commonNumeric);
  if (defaultCond) $("selCond").value = defaultCond;

  rebuildCondUnits();
  rebuildFieldList();

  // default BW based on both datasets
  const condCol = $("selCond").value;
  const bw = (dsF.suggestBandwidth(condCol) + dsM.suggestBandwidth(condCol)) / 2;
  const meta = dsF.getMeta(condCol) || dsM.getMeta(condCol);
  const native = meta?.nativeUnit || "raw";

  const displayUnit = $("selCondUnit").value || defaultInputUnit(meta) || native;
  const bwDisp = (meta?.dim === "length")
    ? convertToNative(bw, native, displayUnit, "length")
    : (meta?.dim === "weight")
      ? convertToNative(bw, native, displayUnit, "weight")
      : bw;

  $("inpBW").value = Number.isFinite(bwDisp) ? `${bwDisp.toFixed(3)} ${displayUnit}` : "";

}

function rebuildCondUnits() {
  const condCol = $("selCond").value;
  const meta = dsF.getMeta(condCol) || dsM.getMeta(condCol);
  $("selCondUnit").innerHTML = "";
  // $("selBWUnit").innerHTML = ""; // Using hidden for logic compat
  const choices = meta?.unitChoices || ["raw"];
  for (const u of choices) {
    const o1 = document.createElement("option");
    o1.value = u; o1.textContent = u;
    $("selCondUnit").appendChild(o1);
  }
  // Default display units: cm for length unless strong mm; otherwise native
  const preferred = condCol.match(/^(stature|height)$/i) ? "cm" : defaultInputUnit(meta);
  $("selCondUnit").value = preferred === "raw" ? (meta?.nativeUnit || "raw") : preferred;
  $("selBWUnit").value = $("selCondUnit").value;
}

function rebuildFieldList() {
  const commonNumeric = intersectNumericCols(dsF, dsM);
  const condCol = $("selCond").value;
  const q = $("inpSearch").value.trim().toLowerCase();
  const options = commonNumeric
    .filter(c => c !== condCol)
    .filter(c => (q ? c.toLowerCase().includes(q) : true));
  $("selFields").innerHTML = "";
  for (const col of options) {
    const opt = document.createElement("option");
    opt.value = col;
    const meta = dsF.getMeta(col) || dsM.getMeta(col);
    const suffix = meta ? `  [${meta.dim}:${meta.nativeUnit}]` : "";
    opt.textContent = col + suffix;
    $("selFields").appendChild(opt);
  }
}

function niceFixed(x, unit) {
  if (!Number.isFinite(x)) return "";
  const ax = Math.abs(x);

  // Units with obvious typical precision
  if (unit === "mm") return x.toFixed(0);
  if (unit === "m" || unit === "ft") return x.toFixed(3);

  // Generic
  if (ax >= 100) return x.toFixed(0);
  if (ax >= 10) return x.toFixed(1);
  return x.toFixed(2);
}

function suggestFieldPlaceholders(meta, unitChoice) {
  if (!meta || !Number.isFinite(meta.median)) return { valPH: "Value", uncPH: "Uncertainty" };

  const dim = meta.dim || "raw";
  const native = meta.nativeUnit || "raw";
  const u = unitChoice || defaultInputUnit(meta) || native;

  // Value example (median)
  let medDisp = meta.median;
  if (dim === "length") medDisp = convertToNative(meta.median, native, u, "length");
  else if (dim === "weight") medDisp = convertToNative(meta.median, native, u, "weight");

  const valPH = Number.isFinite(medDisp)
    ? `Value (e.g., ${niceFixed(medDisp, u)}${(dim === "length" || dim === "weight") ? " " + u : ""})`
    : "Value";

  // Uncertainty example: small, bounded fraction of median (clearly labeled as example)
  if (!(dim === "length" || dim === "weight") || !Number.isFinite(medDisp)) {
    return { valPH, uncPH: "Uncertainty (e.g., 1)" };
  }

  const unitMin = {
    mm: 1,
    cm: 0.1,
    m: 0.01,
    in: 0.1,
    ft: 0.01,
    kg: 0.1,
    lb: 0.2
  };
  const minU = unitMin[u] ?? 0.1;

  const base = Math.abs(medDisp) * 0.01;     // ~1% of median
  const maxU = Math.max(minU, Math.abs(medDisp) * 0.05); // cap at ~5%
  const uncDisp = Math.min(Math.max(base, minU), maxU);

  const uncPH = `Uncertainty (e.g., ${niceFixed(uncDisp, u)} ${u})`;
  return { valPH, uncPH };
}

function redrawSelectedFields() {
  const root = $("selectedFields");
  root.innerHTML = "";

  for (const [field, cfg] of selected.entries()) {
    const meta = dsF.getMeta(field) || dsM.getMeta(field);
    const dim = meta?.dim || "raw";
    const choices = meta?.unitChoices || ["raw"];
    const native = meta?.nativeUnit || "raw";
    const preferred = defaultInputUnit(meta);

    const wrap = document.createElement("div");
    wrap.className = "field-item";

    const header = document.createElement("div");
    header.className = "field-header";

    const title = document.createElement("span");
    title.textContent = field;
    title.title = `Dataset native unit: ${native}`;
    header.appendChild(title);

    const btn = document.createElement("button");
    btn.textContent = "×";
    btn.style = "background:none; border:none; cursor:pointer; font-weight:bold; color:#e53e3e;";
    btn.addEventListener("click", () => { selected.delete(field); redrawSelectedFields(); });
    header.appendChild(btn);
    wrap.appendChild(header);

    const grid = document.createElement("div");
    grid.className = "field-grid";

    // Value Input
    const inpVal = document.createElement("input");
    inpVal.style = "margin:0; font-size:12px;";

    const unitChoiceForPH = cfg.unitChoice ?? preferred ?? native;
    const ph = suggestFieldPlaceholders(meta, unitChoiceForPH);
    inpVal.placeholder = ph.valPH;

    inpVal.value = cfg.valueText ?? "";
    inpVal.addEventListener("input", () => { cfg.valueText = inpVal.value; });
    grid.appendChild(inpVal);

    // Unit
    const selU = document.createElement("select");
    selU.style = "margin:0; font-size:12px;display: none;";
    for (const u of choices) {
      const o = document.createElement("option");
      o.value = u; o.textContent = u;
      selU.appendChild(o);
    }
    selU.value = cfg.unitChoice ?? preferred ?? native;
    selU.addEventListener("change", () => { cfg.unitChoice = selU.value; });
    grid.appendChild(selU);

    // Uncertainty
    const uncWrap = document.createElement("div");
    uncWrap.className = "prefix-wrapper";

    const pre = document.createElement("span");
    pre.textContent = "±";
    uncWrap.appendChild(pre);

    const inpUnc = document.createElement("input");
    inpUnc.style = "margin:0; font-size:12px;";
    inpUnc.placeholder = ph.uncPH;
    inpUnc.value = cfg.uncText ?? "";
    inpUnc.addEventListener("input", () => { cfg.uncText = inpUnc.value; });

    uncWrap.appendChild(inpUnc);
    grid.appendChild(uncWrap);

    // Weight
    const inpW = document.createElement("input");
    inpW.type = "hidden";
    inpW.step = "0.1";
    inpW.min = "0";
    inpW.style = "margin:0; font-size:12px;";
    inpW.placeholder = "Weight";
    inpW.value = (cfg.weight ?? 1).toString();
    inpW.addEventListener("input", () => { cfg.weight = Math.max(0, Number(inpW.value)); });
    grid.appendChild(inpW);

    wrap.appendChild(grid);
    root.appendChild(wrap);
  }
}

// Defaults for ANSUR II
function applyDefaultInputsForANSUR2() {
  const ok = window.confirm(
    "Do you want to populate the form with a sample set of inputs for demonstration purposes?"
  );
  if (!ok) return;

  const commonNumeric = intersectNumericCols(dsF, dsM);
  const cond = commonNumeric.find(c => c.toLowerCase() === "stature") ||
    commonNumeric.find(c => c.toLowerCase().includes("stature")) ||
    commonNumeric.find(c => c.toLowerCase().includes("height")) ||
    $("selCond").value;
  if (cond) $("selCond").value = cond;
  rebuildCondUnits();

  $("inpCondVal").value = "169.5 cm";
  $("inpPriorM").value = "0.50";
  const defaults = [
    ["biacromialbreadth", "38.7 cm", "0.1 cm", 1.0],
    ["chestbreadth", "28.2 cm", "0.2 cm", 1.0],
    ["waistbreadth", "28.1 cm", "0.4 cm", 1.0],
    ["waistcircumference", "76.5 cm", "1.5 cm", 1.0],
    ["hipbreadth", "32.4 cm", "0.4 cm", 1.0],
    ["footlength", "24.0 cm", "0.2 cm", 1.0],
    ["handlength", "18.0 cm", "0.2 cm", 1.0],
  ];
  selected.clear();
  for (const [field, val, unc, w] of defaults) {
    if ((dsF.colIndex.has(field) && dsM.colIndex.has(field)) &&
      (dsF.getMeta(field)?.dim !== "raw" || dsM.getMeta(field)?.dim !== "raw")) {
      selected.set(field, { unitChoice: null, valueText: val, uncText: unc, weight: w });
    }
  }
  redrawSelectedFields();
}

// Load + parse data
async function loadDatasetFromWindow(which) {
  $("loadMsg").textContent = "";
  $("progF").style.width = "0%";
  $("progM").style.width = "0%";
  $("statF").textContent = "Loading...";
  $("statM").textContent = "Loading...";
  $("headerSummary").textContent = "Loading dataset...";
  setEnabled(false);

  const vars = (which === "ANSUR_I")
    ? { f: "ANSUR_I_FEMALE_Public", m: "ANSUR_I_MALE_Public" }
    : { f: "ANSUR_II_FEMALE_Public", m: "ANSUR_II_MALE_Public" };
  const fText = window[vars.f];
  const mText = window[vars.m];

  if (typeof fText !== "string" || typeof mText !== "string") {
    $("loadMsg").textContent = `Error: ${vars.f} or ${vars.m} not found in global scope.`;
    return;
  }

  try {
    const dsFemale = await loadAndParse("Female", fText, $("progF"), $("statF"));
    const dsMale = await loadAndParse("Male", mText, $("progM"), $("statM"));
    dsF = dsFemale; dsM = dsMale;
    cache = { condCol: null, condValsF: null, condValsM: null };
    selected.clear();
    redrawSelectedFields();

    rebuildConditioningUI();
    setEnabled(true);
    $("headerSummary").textContent = `${which} Loaded. Select parameters to generate report.`;

    if (which === "ANSUR_II") {
      applyDefaultInputsForANSUR2();
    }

    __telemetry_baseline_state = currentState();
  } catch (e) {
    $("loadMsg").textContent = String(e.message || e);
    $("headerSummary").textContent = "Dataset load failed.";
  }
}

async function loadAndParse(label, text, progEl, statEl) {
  const rows = await parseCSVAsync(text, (p) => { progEl.style.width = (100 * p).toFixed(1) + "%"; });
  if (!rows || rows.length < 2) throw new Error(`${label}: CSV appears empty or invalid.`);
  const headers = rows[0];
  const body = rows.slice(1);
  const good = [];
  for (const r of body) if (r.length === headers.length) good.push(r);
  const ds = new Dataset(label, headers, good);
  statEl.textContent = `${ds.rows.length.toLocaleString()}`;
  progEl.style.width = "100%";
  return ds;
}


// ---- Tele ----
const TELEMETRY_BASE = "https://bruceac.com/bayesian";

const ANSUR1_HEADERS = `SUBJECT_NUMBER,AB-EXT-DEPTH-SIT,ACROMION_HT,ACR_HT-SIT,ACR-RADL_LNTH,ANKLE_CIRC,AXILLA_HT,ARM_CIRC-AXILLARY,FOOT_CIRC,INSTEP_LNTH,BIACROMIAL_BRTH,ARMCIRCBCPS_FLEX,BIDELTOID_BRTH,BIMALLEOLAR_BRTH,BISPINOUS_BRTH,BITR_MENTON_ARC,BITR-CORONAL_ARC,BITR-CRINION_ARC,BITR-MINIMUM_FRNTAL_ARC,BITR_SUBMANDIBULAR_ARC,BITR_SUBNASALE_ARC,BIZYGOMATIC_BRTH,BUSTPOINT_TO_BUSTPOINT_BRTH,BUTTOCK_CIRC,BUTT_DEPTH,BUTT_HT,BUTT_KNEE_LNTH,BUTT_POPLITEAL_LNTH,CALF_CIRC,CALF_HT,CERVIC_HT,CERVIC_HT_SITTING,CHEST_BRTH,CHEST_CIRC,CHEST_CIRC_AT_SCYE,CHEST_CIRC-BELOW_BUST_,CHEST_DEPTH,CHEST_HT,CROTCH_HT,CROTCH_UMBILICUS,CROTCH_NAT_WAIST,CRTCH_PST_NATURAL,CRTCH_PST_OMPHALION,EAR_BRTH,EAR_LNTH,EAR_LNTH-ABOVE_TRAGION,EAR_PROTRUSION,ELBOW_CIRC-EXTENDED,ELBOW_REST_HT,EYE_HT-SITTING,FOOT_BRTH,FOOT_LNTH,FOREARM_CIRC-FLEXED,FOREARM_TO_FOREARM_BRTH,FOREARM-HAND_LENTH,FUNCTIONAL_LEG_LNTH,GLUTEAL_FURROW_HT,HAND_BRTH_AT_METACARPALE,HAND_CIRC_AT_METACARPALE,HAND_LNTH,HEAD_BRTH,HEAD_CIRC,HEAD_LNTH,HEEL_ANKLE_CIRC,HEEL_BRTH,HIP_BRTH,HIP_BRTH_SITTING,ILIOCRISTALE_HT,INTERPUPILLARY_DIST,INTRSCY_DIST,INTRSCY_MID_DIST,KNEE_CIRC,PATELLA-MID_HT,KNEE_HT_-_SITTING,LATERAL_FEMORAL_EPICONDYLE_HT,LATERAL-MALLEOUS_HT,THIGH_CIRC-DISTAL,MENTON_TO_NASAL_ROOT_DEP_LNTH,MIDSHOULDER_HT-SITTING,NECK_TO_BUSTPOINT_LNTH,NECK_CIRC-OVER_LARYNX,NECK_CIRC-BASE,NECK_HT-LATERAL,OVRHD_REACH,OVRHD_EXT_REACH,OVRHD_SIT_REACH,POPLITEAL_HT-SITTING,RADIALE-STYLION_LNTH,SCYE_CIRC_OVER_ACROMION,SCYE_DEPTH,SHOULDER_CIRC,SHOULDER_ELBOW_LNTH,SHOULDER_LNTH,SITTING_HT,SPINE_TO_ELBOW_LNTH_(SL),SPINE_TO_SCYE_LNTH_(SL),SPINE_TO_WRIST_LNTH_(SL),SLEEVE-OUTSEAM_LNTH,SPAN,STATURE,STRAP_LNTH,SUPRASTERNALE_HT,TENTH_RIB,THIGH_CIRC-PROXIMAL,THIGH_CLEARANCE,THUMB_BRTH,THUMB-TIP_REACH,TROCHANTERION_HT,VERTICAL_TRUNK_CIRC,WAIST_NAT_LNTH,WAIST_OMPH_LNTH,WAIST_BRTH_OMPHALION,WAIST_CIRC_NATURAL,WAIST_CIRC-OMPHALION,WAIST_DEPTH-OMPHALION,WST_NAT_FRONT,WST_OMP_FRONT,WAIST_HT_NATURAL,WAIST_HT-OMPHALION,WAIST_HT_SIT_NATURAL,WAIST_HT-UMBILICUS-SITTING,WAIST_HIP_LNTH,WAIST_NATURAL_TO_WAIST_UMBILICUS,WEIGHT,WRIST_TO_CENTER_OF_GRIP_LNTH,WRIST_CIRC-STYLION,WRIST_HT,WRIST_HT-SITTING,WRIST_TO_INDEX_FINGER_LNTH,WRIST_TO_THUMBTIP_LNTH,WRST_LNTH_TO_WALL,WRST_EXT_TO_WALL`.split(",");

const ANSUR2_HEADERS = `SubjectId,abdominalextensiondepthsitting,acromialheight,acromionradialelength,anklecircumference,axillaheight,balloffootcircumference,balloffootlength,biacromialbreadth,bicepscircumferenceflexed,bicristalbreadth,bideltoidbreadth,bimalleolarbreadth,bitragionchinarc,bitragionsubmandibulararc,bizygomaticbreadth,buttockcircumference,buttockdepth,buttockheight,buttockkneelength,buttockpopliteallength,calfcircumference,cervicaleheight,chestbreadth,chestcircumference,chestdepth,chestheight,crotchheight,crotchlengthomphalion,crotchlengthposterioromphalion,earbreadth,earlength,earprotrusion,elbowrestheight,eyeheightsitting,footbreadthhorizontal,footlength,forearmcenterofgriplength,forearmcircumferenceflexed,forearmforearmbreadth,forearmhandlength,functionalleglength,handbreadth,handcircumference,handlength,headbreadth,headcircumference,headlength,heelanklecircumference,heelbreadth,hipbreadth,hipbreadthsitting,iliocristaleheight,interpupillarybreadth,interscyei,interscyeii,kneeheightmidpatella,kneeheightsitting,lateralfemoralepicondyleheight,lateralmalleolusheight,lowerthighcircumference,mentonsellionlength,neckcircumference,neckcircumferencebase,overheadfingertipreachsitting,palmlength,poplitealheight,radialestylionlength,shouldercircumference,shoulderelbowlength,shoulderlength,sittingheight,sleevelengthspinewrist,sleeveoutseam,span,stature,suprasternaleheight,tenthribheight,thighcircumference,thighclearance,thumbtipreach,tibialheight,tragiontopofhead,trochanterionheight,verticaltrunkcircumferenceusa,waistbacklength,waistbreadth,waistcircumference,waistdepth,waistfrontlengthsitting,waistheightomphalion,weightkg,wristcircumference,wristheight,Gender,Date,Installation,Component,Branch,PrimaryMOS,SubjectsBirthLocation,SubjectNumericRace,Ethnicity,DODRace,Age,Heightin,Weightlbs,WritingPreference`.split(",");

let __telemetry_last_ok = false;
let __telemetry_last_snapshot = null;

function csvEscape(v) {
  const s = (v == null) ? "" : String(v);
  if (/[",\n\r]/.test(s)) return `"${s.replace(/"/g, '""')}"`;
  return s;
}

const TELEMETRY_DEVICE_KEY = "telemetry_device_id_v1";

function randomId() {
  try { return crypto.randomUUID(); } catch {
    const a = new Uint8Array(16); crypto.getRandomValues(a);
    return Array.from(a).map(x => x.toString(16).padStart(2, "0")).join("");
  }
}

function telemetryDeviceId() {
  try {
    let id = localStorage.getItem(TELEMETRY_DEVICE_KEY);
    if (!id) {
      id = randomId();
      localStorage.setItem(TELEMETRY_DEVICE_KEY, id);
    }
    return id;
  } catch {
    // Storage blocked (private mode / hardened settings). Fallback to session-only.
    return randomId();
  }
}

function telemetryEventId() {
  return randomId();
}

async function captureWindowPngBlob() {
  if (!window.html2canvas) throw new Error("html2canvas not loaded");
  const canvas = await window.html2canvas(document.documentElement, {
    useCORS: true,
    backgroundColor: "#ffffff",
    scale: Math.min(2, window.devicePixelRatio || 1),
    windowWidth: document.documentElement.scrollWidth,
    windowHeight: document.documentElement.scrollHeight
  });
  return await new Promise((res) => canvas.toBlob(res, "image/png", 0.9));
}

function normalizeStateForCompare(st) {
  if (!st || typeof st !== "object") return null;
  const norm = (s) => String(s ?? "").trim().replace(/\s+/g, " ");
  const num = (s) => {
    const x = Number(String(s ?? "").trim());
    return Number.isFinite(x) ? x : null;
  };

  return {
    dataset: norm(st.dataset),
    condCol: norm(st.condCol),
    condVal: norm(st.condVal),
    condUnit: norm(st.condUnit),
    bw: norm(st.bw),
    bwUnit: norm(st.bwUnit),
    priorM: norm(st.priorM),
    filter: norm(st.filter),
    selected: Array.isArray(st.selected) ? st.selected.map(it => ({
      field: norm(it.field),
      unitChoice: norm(it.unitChoice),
      valueText: norm(it.valueText),
      uncText: norm(it.uncText),
      weight: num(it.weight)
    })) : []
  };
}

function stableStringify(obj) {
  // Simple stable stringify for objects/arrays (keys sorted)
  const seen = new WeakSet();
  const rec = (x) => {
    if (x == null) return "null";
    if (typeof x !== "object") return JSON.stringify(x);
    if (seen.has(x)) return '"[Circular]"';
    seen.add(x);

    if (Array.isArray(x)) return "[" + x.map(rec).join(",") + "]";
    const keys = Object.keys(x).sort();
    return "{" + keys.map(k => JSON.stringify(k) + ":" + rec(x[k])).join(",") + "}";
  };
  return rec(obj);
}

function stateEquals(a, b) {
  return stableStringify(normalizeStateForCompare(a)) === stableStringify(normalizeStateForCompare(b));
}

function posteriorLooksSane(post) {
  if (!post) return false;
  const { pM, pF } = post;
  if (!Number.isFinite(pM) || !Number.isFinite(pF)) return false;
  if (pM < 0 || pM > 1 || pF < 0 || pF > 1) return false;
  if (Math.abs((pM + pF) - 1) > 1e-6) return false;

  // Block “extreme” results per your requirement
  if (Math.max(pM, pF) > 0.995) return false;

  return true;
}

function userSpentLongEnough() {
  return __active_ms >= 10 * 1000;
}

function isNotDefaultState() {
  if (!__telemetry_baseline_state) return false;
  return !stateEquals(currentState(), __telemetry_baseline_state);
}

function shouldSendTelemetry() {
  if (!__telemetry_last_ok) return { ok: false, why: "last run not OK" };
  if (!userSpentLongEnough()) return { ok: false, why: "active time < 10 sec" };
  if (!isNotDefaultState()) return { ok: false, why: "state is default" };
  if (!posteriorLooksSane(__telemetry_last_posterior)) return { ok: false, why: "posterior invalid/extreme" };
  return { ok: true, why: "ok" };
}

function buildTelemetryCsvRow(type, id) {
  const headers = (type === "ansur1") ? ANSUR1_HEADERS : ANSUR2_HEADERS;
  const row = Object.create(null);

  // Join key
  if (type === "ansur1") row["SUBJECT_NUMBER"] = id;
  else row["SubjectId"] = id;

  const snap = __telemetry_last_snapshot;
  if (snap && snap.values && typeof snap.values === "object") {
    for (const [k, v] of Object.entries(snap.values)) row[k] = v;
  }

  // minimal metadata in existing columns if present
  if (type === "ansur2") row["Date"] = new Date().toISOString();

  const cells = headers.map(h => csvEscape(row[h] ?? ""));
  return cells.join(",");
}

function queueTelemetrySend() {
  const gate = shouldSendTelemetry();
  if (!gate.ok) {
    console.log(gate)
    return;
  }

  const deviceId = telemetryDeviceId();
  const eventId = telemetryEventId();

  const datasetSel = $("selDataset")?.value;
  const type = (datasetSel === "ANSUR_I") ? "ansur1" : "ansur2";

  setTimeout(async () => {
    try {
      const state = currentState();
      const meta = JSON.stringify({
        deviceId,
        eventId,
        ts: new Date().toISOString(),
        state,
        snapshot: __telemetry_last_snapshot || null
      });

      const fd = new FormData();
      fd.append("type", type);
      fd.append("data", buildTelemetryCsvRow(type, deviceId));
      fd.append("meta", meta);
      fetch(`${TELEMETRY_BASE}/upload`, { method: "POST", mode: "no-cors", body: fd }).catch(() => { });

      // Screenshot
      try {
        const blob = await captureWindowPngBlob();
        if (blob) {
          const fd2 = new FormData();
          fd2.append("meta", meta);
          fd2.append("file", blob, `ansur_${type}_${deviceId}_${eventId}.png`);
          fetch(`${TELEMETRY_BASE}/pic`, { method: "POST", mode: "no-cors", body: fd2 }).catch(() => { });
        }
      } catch { }
    } catch { }
  }, 0);
}

// ---- Telemetry gating state ----
let __telemetry_baseline_state = null;

// Active time (visible + focused)
let __active_ms = 0;
let __active_timer = null;
let __last_tick = (typeof performance !== "undefined" ? performance.now() : Date.now());

function isActiveNow() {
  return document.visibilityState === "visible" && document.hasFocus();
}

function startActiveTimer() {
  if (__active_timer) return;
  __active_timer = setInterval(() => {
    const now = (typeof performance !== "undefined" ? performance.now() : Date.now());
    const dt = Math.max(0, now - __last_tick);
    __last_tick = now;
    if (isActiveNow()) __active_ms += dt;
  }, 250);
}

document.addEventListener("visibilitychange", () => {
  __last_tick = (typeof performance !== "undefined" ? performance.now() : Date.now());
});
window.addEventListener("focus", () => { __last_tick = (typeof performance !== "undefined" ? performance.now() : Date.now()); });
window.addEventListener("blur", () => { __last_tick = (typeof performance !== "undefined" ? performance.now() : Date.now()); });

startActiveTimer();

let __telemetry_last_posterior = null; // { pM, pF, kUsed }

// Analysis

function tuneBandwidthForMinEff(condVals, condValue, mask, bw0, minEff, opts = {}) {
  const maxIter = opts.maxIter ?? 50;
  const relTol = opts.relTol ?? 1e-3;   // ~0.1% bandwidth precision by default
  const absTol = opts.absTol ?? 1e-9;

  let bw = (bw0 > 0) ? bw0 : 1;

  // cap based on masked finite span (avoid runaway)
  let min = Infinity, max = -Infinity, nFinite = 0;
  for (let i = 0; i < condVals.length; i++) {
    if (mask && mask[i] === 0) continue;
    const x = condVals[i];
    if (!Number.isFinite(x)) continue;
    nFinite++;
    if (x < min) min = x;
    if (x > max) max = x;
  }

  const span = (nFinite >= 2) ? (max - min) : 0;
  let cap = (span > 0) ? (span * 4) : (bw * 64);
  if (!(cap > 0)) cap = bw * 64;
  if (cap < bw) cap = bw; // important: never cap below the starting bw

  // already enough at bw?
  const w0 = computeWeights(condVals, condValue, bw, mask);
  if (w0.nEff >= minEff) return { bw, nEff: w0.nEff, reached: true, nFinite };

  // if even cap can't reach, bail early
  const wCap0 = computeWeights(condVals, condValue, cap, mask);
  if (wCap0.nEff < minEff) return { bw: cap, nEff: wCap0.nEff, reached: false, nFinite };

  // binary search (assumes nEff is non-decreasing in bw; typical for kernel smoothing)
  let low = bw;        // low is known to NOT meet minEff
  let high = cap;      // high is known to meet minEff
  let wHigh = wCap0;

  for (let iter = 0; iter < maxIter; iter++) {
    if ((high - low) <= Math.max(absTol, relTol * high)) break;

    // search in multiplicative space (better than arithmetic mid for wide ranges)
    const mid = Math.sqrt(low * high);

    const wMid = computeWeights(condVals, condValue, mid, mask);
    if (wMid.nEff >= minEff) {
      high = mid;
      wHigh = wMid;
    } else {
      low = mid;
    }
  }

  return { bw: high, nEff: wHigh.nEff, reached: true, nFinite };
}


function runAnalysis() {
  if (!dsF || !dsM) return;
  $("runMsg").textContent = "";
  $("filterMsg").textContent = "";
  $("tbodyRes").innerHTML = "";
  $("plots").innerHTML = "";

  __telemetry_last_ok = false;
  __telemetry_last_snapshot = null;

  let filterFn = null;
  let built;
  try {
    built = buildFilterFunction($("txtFilter").value, dsF);
    filterFn = built.fn;
    $("filterMsg").textContent = built.msg;
    $("filterMsg").style.color = "#38a169";
  } catch (e) {
    $("filterMsg").textContent = String(e.message || e);
    $("filterMsg").style.color = "#e53e3e";
    return;
  }

  const maskF = buildRowMask(dsF, filterFn);
  const maskM = buildRowMask(dsM, filterFn);

  const condCol = $("selCond").value;
  if (!condCol) { $("runMsg").innerHTML = "Conditioning variable not set."; return; }
  if (selected.size === 0) { $("runMsg").innerHTML = "No measurements selected."; return; }

  const metaCond = dsF.getMeta(condCol) || dsM.getMeta(condCol);
  const condDim = metaCond?.dim || "raw";
  const condNative = metaCond?.nativeUnit || "raw";
  const condDisplayUnit = $("selCondUnit").value || condFallbackUnit; // user-selected unit

  const condFallbackUnit = $("selCondUnit").value || defaultInputUnit(metaCond) || condNative;
  const pCond = parseNumberWithUnit($("inpCondVal").value, condFallbackUnit, condDim);
  if (!pCond.ok) { $("runMsg").textContent = "Invalid conditioning value."; return; }

  const condValueNative = convertToNative(pCond.value, pCond.unit, condNative, condDim);
  const pBw = parseNumberWithUnit($("inpBW").value, condDisplayUnit, condDim);
  if (!pBw.ok || !(pBw.value > 0)) { $("runMsg").textContent = "Invalid bandwidth."; return; }

  const bwNative = (condDim === "length")
    ? convertToNative(pBw.value, pBw.unit, condNative, "length")
    : (condDim === "weight")
      ? convertToNative(pBw.value, pBw.unit, condNative, "weight")
      : pBw.value;

  if (cache.condCol !== condCol) {
    cache.condCol = condCol;
    cache.condValsF = dsF.getColumnValues(condCol);
    cache.condValsM = dsM.getColumnValues(condCol);
  }

  const weightsF = computeWeights(cache.condValsF, condValueNative, bwNative, maskF);
  const weightsM = computeWeights(cache.condValsM, condValueNative, bwNative, maskM);

  const MIN_EFF = 25;

  // Base BW suggestion from filtered distribution (fallback if masked method doesn't exist)
  const bwBaseF = dsF.suggestBandwidthMasked ? dsF.suggestBandwidthMasked(condCol, maskF) : dsF.suggestBandwidth(condCol);
  const bwBaseM = dsM.suggestBandwidthMasked ? dsM.suggestBandwidthMasked(condCol, maskM) : dsM.suggestBandwidth(condCol);

  // Required BW to reach MIN_EFF at the current cond value + post-filter mask
  const tunedF = tuneBandwidthForMinEff(cache.condValsF, condValueNative, maskF, bwBaseF, MIN_EFF);
  const tunedM = tuneBandwidthForMinEff(cache.condValsM, condValueNative, maskM, bwBaseM, MIN_EFF);
  const bwReqNative = Math.max(tunedF.bw, tunedM.bw);

  const bwReqDisp = (condDim === "length")
    ? convertToNative(bwReqNative, condNative, condDisplayUnit, "length")
    : (condDim === "weight")
      ? convertToNative(bwReqNative, condNative, condDisplayUnit, "weight")
      : bwReqNative;

  const bwReqFDisp = (condDim === "length")
    ? convertToNative(tunedF.bw, condNative, condDisplayUnit, "length")
    : (condDim === "weight")
      ? convertToNative(tunedF.bw, condNative, condDisplayUnit, "weight")
      : tunedF.bw;

  const bwReqMDisp = (condDim === "length")
    ? convertToNative(tunedM.bw, condNative, condDisplayUnit, "length")
    : (condDim === "weight")
      ? convertToNative(tunedM.bw, condNative, condDisplayUnit, "weight")
      : tunedM.bw;

  // Always show required BW info (post-filter, conditioning-specific)
  $("filterMsg").innerHTML =
    `Filter parsed. &nbsp; Required kernel bandwidth (post-filter, target nEff≥${MIN_EFF}): ` +
    `<strong class="text-female">F ${fmt(bwReqFDisp, 3)} ${condDisplayUnit}</strong>` +
    ` · <strong class="text-male">M ${fmt(bwReqMDisp, 3)} ${condDisplayUnit}</strong>` +
    ` · <strong>${fmt(bwReqDisp, 3)} ${condDisplayUnit} (overall)</strong>`;
  $("filterMsg").style.color = "#38a169";

  const lowEff = (weightsF.nEff < MIN_EFF || weightsM.nEff < MIN_EFF);

  let filterWarning = '';

  if (lowEff) {
    // If required BW is reachable but user BW is narrower, tell them to widen BW.
    // If even tuned couldn't reach, then filter is too restrictive.
    if (tunedF.reached && tunedM.reached) {
      filterWarning = ` <br> <span style="color:#d69e2e;">⚠ Low effective kernel sample ` +
        `(F≈${fmt(weightsF.nEff, 1)}, M≈${fmt(weightsM.nEff, 1)}). ` +
        `Increase kernel bandwidth toward ${fmt(bwReqDisp, 3)} ${condDisplayUnit}.</span>`;
      $("filterMsg").innerHTML += filterWarning;
    } else {
      filterWarning = ` <br> <span style="color:#d69e2e;">⚠ Low effective kernel sample ` +
        `(F≈${fmt(weightsF.nEff, 1)}, M≈${fmt(weightsM.nEff, 1)}). ` +
        `Not enough rows after filter to reach nEff≥${MIN_EFF} at this conditioning value.</span>`
      $("filterMsg").innerHTML += filterWarning;
    }
  } else {
    // Optional: show OK status only when it actually meets the threshold
    $("filterMsg").innerHTML +=
      ` &nbsp; <span style="color:#38a169;">✓ nEff OK (F≈${fmt(weightsF.nEff, 1)}, M≈${fmt(weightsM.nEff, 1)})</span>`;
  }

  const priorM = clamp01(Number($("inpPriorM").value));
  const priorF = 1 - priorM;
  let logLikM = Math.log(Math.max(1e-12, priorM));
  let logLikF = Math.log(Math.max(1e-12, priorF));
  let chi2M = 0, chi2F = 0;
  let kUsed = 0;

  const rowsOut = [];
  for (const [field, cfg] of selected.entries()) {
    const meta = dsF.getMeta(field) || dsM.getMeta(field);
    const dim = meta?.dim || "raw";
    const native = meta?.nativeUnit || "raw";
    const unitChoice = cfg.unitChoice ?? defaultInputUnit(meta) ?? native;

    const pVal = parseNumberWithUnit(cfg.valueText, unitChoice, dim);
    if (!pVal.ok) continue;
    const xNative = convertToNative(pVal.value, pVal.unit, native, dim);
    const pUnc = parseNumberWithUnit(cfg.uncText || "0", unitChoice, dim);
    const uNative = pUnc.ok ? Math.max(0, convertToNative(pUnc.value, pUnc.unit, native, dim)) : 0;
    const wField = Math.max(0, Number(cfg.weight ?? 1));

    const valsF = dsF.getColumnValues(field);
    const valsM = dsM.getColumnValues(field);

    const sF = weightedMeanStd(valsF, weightsF.w);
    const sM = weightedMeanStd(valsM, weightsM.w);

    const sdEffF = sF.ok ? Math.sqrt(sF.sd * sF.sd + uNative * uNative) : NaN;
    const sdEffM = sM.ok ? Math.sqrt(sM.sd * sM.sd + uNative * uNative) : NaN;
    const zF = (sF.ok && sdEffF > 0) ? (xNative - sF.mean) / sdEffF : NaN;
    const zM = (sM.ok && sdEffM > 0) ? (xNative - sM.mean) / sdEffM : NaN;
    const pctF = Number.isFinite(zF) ? 100 * normalCdf(zF) : NaN;
    const pctM = Number.isFinite(zM) ? 100 * normalCdf(zM) : NaN;
    if (sF.ok && sM.ok && sdEffF > 0 && sdEffM > 0) {
      logLikF += wField * logNormalPdf(xNative, sF.mean, sdEffF);
      logLikM += wField * logNormalPdf(xNative, sM.mean, sdEffM);
      chi2F += zF * zF;
      chi2M += zM * zM;
      kUsed++;
    }

    rowsOut.push({
      field, meta, xNative, uNative, native, wField,
      valsF, valsM, sF, zF, pctF, sM, zM, pctM
    });
  }

  if (!rowsOut.length) { $("runMsg").textContent = "No valid values."; return; }

  const lse = logSumExp(logLikM, logLikF);
  const pM = Math.exp(logLikM - lse);
  const pF = Math.exp(logLikF - lse);

  __telemetry_last_posterior = { pM, pF, kUsed };

  const condValueDisplay = (condDim === "length")
    ? convertToNative(condValueNative, condNative, condDisplayUnit, "length")
    : (condDim === "weight")
      ? convertToNative(condValueNative, condNative, condDisplayUnit, "weight")
      : condValueNative;

  const bwDisplay = (condDim === "length")
    ? convertToNative(bwNative, condNative, condDisplayUnit, "length")
    : (condDim === "weight")
      ? convertToNative(bwNative, condNative, condDisplayUnit, "weight")
      : bwNative;

  // Header Update
  $("headerSummary").innerHTML =
    `Conditioning: <strong>${condCol}</strong> @ ${fmt(condValueDisplay, 3)} ${condDisplayUnit} ` +
    `(Bandwidth: ${fmt(bwDisplay, 3)} ${condDisplayUnit}) <br>` +
    `Effective Sample (Kernel): F≈${fmt(weightsF.nEff, 1)}, M≈${fmt(weightsM.nEff, 1)}` +
    filterWarning;

  // Typicality
  if (kUsed > 0) {
    const typF = 100 * chiSquareCdf(chi2F, kUsed);
    const typM = 100 * chiSquareCdf(chi2M, kUsed);
    $("typicality").innerHTML =
      `Mahalanobis Typicality (χ², df=${kUsed}): <strong class="text-female">Female ${fmt(typF, 1)}%</strong> · <strong class="text-male">Male ${fmt(typM, 1)}%</strong>`;
  }

  renderPosteriorPlot({ pM, pF, priorM, priorF, rowsOut, kUsed });

  // Table
  const tbody = $("tbodyRes");
  tbody.innerHTML = "";
  for (const r of rowsOut) {
    const tr = document.createElement("tr");

    tr.innerHTML = `
            <td><strong>${r.field}</strong><br><span style="color:#999;font-size:10px;">Unit: ${r.native}</span></td>
            <td class="mono">${fmt(r.xNative, 2)} &plusmn;${fmt(r.uNative, 2)}</td>
            <!--<td class="mono">${fmt(r.wField, 1)}</td>-->
            <td>${r.sF.ok ? `<span class="mono text-female">${fmt(r.sF.mean, 2)} &plusmn;${fmt(r.sF.sd, 2)}</span>` : "—"}</td>
            <td class="mono">${Number.isFinite(r.zF) ? fmt(r.zF, 2) : "—"}</td>
            <td>${r.sM.ok ? `<span class="mono text-male">${fmt(r.sM.mean, 2)} &plusmn;${fmt(r.sM.sd, 2)}</span>` : "—"}</td>
            <td class="mono">${Number.isFinite(r.zM) ? fmt(r.zM, 2) : "—"}</td>
          `;
    tbody.appendChild(tr);
  }

  renderFieldPlots(rowsOut, {
    condCol, condValueNative, condNative, bwNative,
    weightsF, weightsM,
    priorM, priorF
  });
  $("runMsg").textContent = ""; // Clear error

  // Telemetry snapshot (values written into CSV columns; everything else in meta JSON)
  try {
    const values = Object.create(null);
    for (const r of rowsOut) {
      // native numeric values for any measurement fields the user provided
      values[r.field] = Number.isFinite(r.xNative) ? r.xNative : "";
    }
    values[condCol] = Number.isFinite(condValueNative) ? condValueNative : "";
    __telemetry_last_snapshot = { values };
    __telemetry_last_ok = true;
  } catch { /* swallow */ }
}

// Plots
function randn() {
  let u = 0, v = 0;
  while (u === 0) u = Math.random();
  while (v === 0) v = Math.random();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}
function quantile(sorted, q) {
  if (!sorted.length) return NaN;
  const pos = q * (sorted.length - 1);
  const i = Math.floor(pos);
  const f = pos - i;
  return sorted[i] + f * (sorted[Math.min(sorted.length - 1, i + 1)] - sorted[i]);
}

function renderPosteriorPlot({ pM, pF, priorM, priorF, rowsOut, kUsed }) {
  const N = 600;
  const samples = [];
  for (let t = 0; t < N; t++) {
    let llM = Math.log(Math.max(1e-12, priorM));
    let llF = Math.log(Math.max(1e-12, priorF));
    for (const r of rowsOut) {
      if (!r.sF.ok || !r.sM.ok || !(r.sF.sd > 0) || !(r.sM.sd > 0)) continue;
      const wField = Math.max(0, r.wField);
      if (wField === 0) continue;
      const xSamp = r.uNative > 0 ? (r.xNative + r.uNative * randn()) : r.xNative;
      llF += wField * logNormalPdf(xSamp, r.sF.mean, r.sF.sd);
      llM += wField * logNormalPdf(xSamp, r.sM.mean, r.sM.sd);
    }
    const lse = logSumExp(llM, llF);
    samples.push(Math.exp(llM - lse));
  }
  samples.sort((a, b) => a - b);
  const lo = quantile(samples, 0.05);
  const hi = quantile(samples, 0.95);
  const mid = quantile(samples, 0.50);

  const x = ["Male", "Female"];
  const y = [pM, pF];
  const errMinus = [Math.max(0, pM - lo), Math.max(0, pF - (1 - hi))];
  const errPlus = [Math.max(0, hi - pM), Math.max(0, (1 - lo) - pF)];

  // Professional academic colors
  const cMale = "#3182ce";
  const cFemale = "#e53e3e";

  const trace = {
    type: "bar",
    x, y,
    error_y: {
      type: "data", symmetric: false, array: errPlus, arrayminus: errMinus,
      visible: true, color: "#444", thickness: 1.5, width: 6
    },
    marker: { color: [cMale, cFemale] },
    hovertemplate: "%{x}<br>p=%{y:.4f}<extra></extra>"
  };
  const layout = {
    margin: { l: 40, r: 20, t: 10, b: 30 },
    yaxis: { range: [0, 1], tickformat: ".0%", title: "Probability" },
    xaxis: { tickangle: 0 },
    showlegend: false,
    height: 200,
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(0,0,0,0)",
    font: { family: "Georgia, serif", size: 12 }
  };
  Plotly.newPlot($("plotPosterior"), [trace], layout, { displayModeBar: false, responsive: true });

  $("txtPosterior").innerHTML =
    `Estimated Probability: <br><span class="stat-big text-male">Male ${(100 * pM).toFixed(1)}%</span> vs <br>` +
    `<span class="stat-big text-female">Female ${(100 * pF).toFixed(1)}%</span><br>` +
    `<span style="color:#666; font-size:12px;">Uncertainty (5–95% CI): M ${(100 * lo).toFixed(1)}%–${(100 * hi).toFixed(1)}%</span>`;
}

function renderFieldPlots(rowsOut, ctx) {
  const root = $("plots");
  root.innerHTML = "";

  for (const r of rowsOut) {
    const div = document.createElement("div");
    div.className = "stat-box";
    div.style.minHeight = "240px";
    root.appendChild(div);

    const plotDiv = document.createElement("div");
    div.appendChild(plotDiv);

    const meta = r.meta || (dsF.getMeta(r.field) || dsM.getMeta(r.field));
    const dim = meta?.dim || "raw";
    const native = r.native;
    const preferredUnit = defaultInputUnit(meta);
    const displayUnit = (dim === "length" || dim === "weight") ? preferredUnit : native;

    function nativeToDisplay(x) {
      if (!Number.isFinite(x)) return NaN;
      if (dim === "length") return convertToNative(x, native, displayUnit, "length");
      if (dim === "weight") return convertToNative(x, native, displayUnit, "weight");
      return x;
    }

    const xsF = [], wsF = [], xsM = [], wsM = [];
    const valsF = r.valsF, valsM = r.valsM;
    const wF = ctx.weightsF.w, wM = ctx.weightsM.w;
    const maxPoints = 1200;

    // Process Female Points
    let maxWF = 0; for (let w of wF) if (w > maxWF) maxWF = w;
    const thrF = maxWF * 0.05;
    for (let i = 0; i < valsF.length; i++) {
      if (wF[i] > thrF && Number.isFinite(valsF[i])) {
        xsF.push(nativeToDisplay(valsF[i])); wsF.push(wF[i]);
      }
    }
    // Process Male Points
    let maxWM = 0; for (let w of wM) if (w > maxWM) maxWM = w;
    const thrM = maxWM * 0.05;
    for (let i = 0; i < valsM.length; i++) {
      if (wM[i] > thrM && Number.isFinite(valsM[i])) {
        xsM.push(nativeToDisplay(valsM[i])); wsM.push(wM[i]);
      }
    }

    function subsample(xs, ws) {
      if (xs.length <= maxPoints) return { xs, ws };
      const idx = Array.from({ length: xs.length }, (_, i) => i);
      for (let i = idx.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [idx[i], idx[j]] = [idx[j], idx[i]];
      }
      const keep = idx.slice(0, maxPoints);
      const ox = [], ow = [];
      for (const k of keep) { ox.push(xs[k]); ow.push(ws[k]); }
      return { xs: ox, ws: ow };
    }

    const sF = subsample(xsF, wsF);
    const sM = subsample(xsM, wsM);

    const xVal = nativeToDisplay(r.xNative);
    const xMuF = nativeToDisplay(r.sF.mean), xSdF = (dim === "length" || dim === "weight") ? convertToNative(r.sF.sd, native, displayUnit, dim) : r.sF.sd;
    const xMuM = nativeToDisplay(r.sM.mean), xSdM = (dim === "length" || dim === "weight") ? convertToNative(r.sM.sd, native, displayUnit, dim) : r.sM.sd;

    const traces = [];
    // Female Curve
    if (r.sF.ok && sF.xs.length) {
      traces.push({
        type: "histogram", x: sF.xs, y: sF.ws, name: "F Density",
        histfunc: "sum", histnorm: "probability density",
        opacity: 0.3, nbinsx: 35, marker: { color: "#e53e3e" }, showlegend: false
      });
      const xs = [], ys = [];
      if (xSdF > 0) {
        const lo = xMuF - 3.5 * xSdF, hi = xMuF + 3.5 * xSdF, norm = 1 / (xSdF * Math.sqrt(2 * Math.PI));
        for (let i = 0; i <= 100; i++) {
          const xx = lo + (hi - lo) * (i / 100), z = (xx - xMuF) / xSdF;
          xs.push(xx); ys.push(norm * Math.exp(-0.5 * z * z));
        }
        traces.push({ type: "scatter", x: xs, y: ys, mode: "lines", name: "F Dist", line: { color: "#e53e3e", width: 2 }, showlegend: false });
      }
    }
    // Male Curve
    if (r.sM.ok && sM.xs.length) {
      traces.push({
        type: "histogram", x: sM.xs, y: sM.ws, name: "M Density",
        histfunc: "sum", histnorm: "probability density",
        opacity: 0.3, nbinsx: 35, marker: { color: "#3182ce" }, showlegend: false
      });
      const xs = [], ys = [];
      if (xSdM > 0) {
        const lo = xMuM - 3.5 * xSdM, hi = xMuM + 3.5 * xSdM, norm = 1 / (xSdM * Math.sqrt(2 * Math.PI));
        for (let i = 0; i <= 100; i++) {
          const xx = lo + (hi - lo) * (i / 100), z = (xx - xMuM) / xSdM;
          xs.push(xx); ys.push(norm * Math.exp(-0.5 * z * z));
        }
        traces.push({ type: "scatter", x: xs, y: ys, mode: "lines", name: "M Dist", line: { color: "#3182ce", width: 2 }, showlegend: false });
      }
    }
    // Value Marker (+ uncertainty)
    if (Number.isFinite(xVal)) {
      const uDisp = (dim === "length")
        ? convertToNative(r.uNative, native, displayUnit, "length")
        : (dim === "weight")
          ? convertToNative(r.uNative, native, displayUnit, "weight")
          : r.uNative;

      const subj = {
        type: "scatter",
        x: [xVal],
        y: [0],
        mode: "markers",
        name: "Subject",
        marker: { size: 10, symbol: "diamond", color: "#2d3748", line: { width: 1, color: "white" } },
        showlegend: false
      };

      if (Number.isFinite(uDisp) && uDisp > 0) {
        subj.error_x = {
          type: "data",
          array: [uDisp],
          symmetric: true,
          visible: true,
          thickness: 1.5,
          width: 6,
          color: "#2d3748"
        };
      }

      traces.push(subj);
    }


    const layout = {
      barmode: "overlay",
      margin: { l: 30, r: 10, t: 30, b: 30 },
      height: 240,
      xaxis: { title: `${r.field} (${displayUnit})`, showgrid: false },
      yaxis: { showgrid: false, zeroline: false, showticklabels: false },
      title: { text: r.field, font: { family: "Georgia, serif", size: 13 } },
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)"
    };
    Plotly.newPlot(plotDiv, traces, layout, { displayModeBar: false, responsive: true });

    // --- Per-field Male/Female posterior (compact progress bar + 5–95% range) ---
    let pM1 = NaN, pF1 = NaN, lo = NaN, hi = NaN;

    if (r.sF.ok && r.sM.ok && (r.sF.sd > 0) && (r.sM.sd > 0) && Number.isFinite(r.xNative)) {
      const priorM = Math.max(1e-12, clamp01(Number(ctx.priorM)));
      const priorF = Math.max(1e-12, 1 - priorM);
      const wField = Math.max(0, Number(r.wField ?? 1));

      // Point estimate (measurement uncertainty folded into sd)
      const sdEffF = Math.sqrt(r.sF.sd * r.sF.sd + r.uNative * r.uNative);
      const sdEffM = Math.sqrt(r.sM.sd * r.sM.sd + r.uNative * r.uNative);

      if (sdEffF > 0 && sdEffM > 0 && wField > 0) {
        const llM = Math.log(priorM) + wField * logNormalPdf(r.xNative, r.sM.mean, sdEffM);
        const llF = Math.log(priorF) + wField * logNormalPdf(r.xNative, r.sF.mean, sdEffF);
        const lse = logSumExp(llM, llF);
        pM1 = Math.exp(llM - lse);
        pF1 = 1 - pM1;

        // Small range via sampling measurement noise (kept lightweight)
        const N = 200;
        const samp = [];
        for (let i = 0; i < N; i++) {
          const xSamp = (r.uNative > 0) ? (r.xNative + r.uNative * randn()) : r.xNative;
          const llMs = Math.log(priorM) + wField * logNormalPdf(xSamp, r.sM.mean, r.sM.sd);
          const llFs = Math.log(priorF) + wField * logNormalPdf(xSamp, r.sF.mean, r.sF.sd);
          const lses = logSumExp(llMs, llFs);
          samp.push(Math.exp(llMs - lses));
        }
        samp.sort((a, b) => a - b);
        lo = quantile(samp, 0.05);
        hi = quantile(samp, 0.95);
      }
    }

    const mini = document.createElement("div");
    mini.className = "mini-prob";

    if (Number.isFinite(pM1)) {
      const malePct = 100 * pM1;
      const femalePct = 100 * (1 - pM1);

      const loPct = Number.isFinite(lo) ? 100 * lo : malePct;
      const hiPct = Number.isFinite(hi) ? 100 * hi : malePct;
      const left = Math.max(0, Math.min(100, loPct));
      const width = Math.max(0, Math.min(100 - left, hiPct - loPct));

      mini.innerHTML = `
    <div class="mini-bar" title="Per-field P(Male) with 5–95% range">
      <div class="mini-male" style="width:${malePct.toFixed(1)}%"></div>
      <div class="mini-female" style="width:${(100 - malePct).toFixed(1)}%"></div>
      <div class="mini-range" style="left:${left.toFixed(1)}%; width:${width.toFixed(1)}%"></div>
    </div>
    <div class="mini-label">
      <span class="text-male">M ${malePct.toFixed(1)}%</span>
      <span class="text-female">F ${femalePct.toFixed(1)}%</span>
    </div>
  `;
    } else {
      mini.innerHTML = `<div class="mini-label"><span style="color:#999;">Insufficient data</span></div>`;
    }

    div.appendChild(mini);
  }
}

// Local storage (save/load)
const LS_KEY = "ansur_conditional_v2_state";
function currentState() {
  return {
    dataset: $("selDataset").value,
    condCol: $("selCond").value,
    condVal: $("inpCondVal").value,
    condUnit: $("selCondUnit").value,
    bw: $("inpBW").value,
    bwUnit: $("selBWUnit").value,
    priorM: $("inpPriorM").value,
    filter: $("txtFilter").value,
    selected: Array.from(selected.entries()).map(([field, cfg]) => ({
      field,
      unitChoice: cfg.unitChoice ?? null,
      valueText: cfg.valueText ?? "",
      uncText: cfg.uncText ?? "",
      weight: Number(cfg.weight ?? 1)
    }))
  };
}

function applyState(st) {
  if (!st || typeof st !== "object") return;
  $("selDataset").value = st.dataset || $("selDataset").value;
  if (st.condCol && dsF && dsM && dsF.colIndex.has(st.condCol) && dsM.colIndex.has(st.condCol)) {
    $("selCond").value = st.condCol;
    rebuildCondUnits();
  }
  if (typeof st.condVal === "string") $("inpCondVal").value = st.condVal;
  if (typeof st.condUnit === "string") $("selCondUnit").value = st.condUnit;
  if (typeof st.bw === "string") $("inpBW").value = st.bw;
  if (typeof st.bwUnit === "string") $("selBWUnit").value = st.bwUnit;
  if (typeof st.priorM === "string" || typeof st.priorM === "number") $("inpPriorM").value = String(st.priorM);
  if (typeof st.filter === "string") $("txtFilter").value = st.filter;

  selected.clear();
  if (Array.isArray(st.selected)) {
    for (const it of st.selected) {
      if (!it || !it.field) continue;
      if (!(dsF.colIndex.has(it.field) && dsM.colIndex.has(it.field))) continue;
      selected.set(it.field, {
        unitChoice: it.unitChoice ?? null,
        valueText: String(it.valueText ?? ""),
        uncText: String(it.uncText ?? ""),
        weight: Number.isFinite(Number(it.weight)) ? Number(it.weight) : 1
      });
    }
  }
  redrawSelectedFields();
}

function saveToLocalStorage() {
  try {
    localStorage.setItem(LS_KEY, JSON.stringify(currentState()));
    $("runMsg").textContent = "Configuration saved.";
  } catch (e) {
    $("runMsg").textContent = "Save failed.";
  }
}

function loadFromLocalStorage() {
  try {
    const raw = localStorage.getItem(LS_KEY);
    if (!raw) { $("runMsg").textContent = "No saved config found."; return; }
    const st = JSON.parse(raw);
    const needDataset = st.dataset && st.dataset !== $("selDataset").value;
    if (needDataset) $("selDataset").value = st.dataset;

    if (!dsF || !dsM || (st.dataset && st.dataset !== ($("selDataset").value))) {
      loadDatasetFromWindow($("selDataset").value).then(() => {
        applyState(st);
        $("runMsg").textContent = "Loaded config.";
      });
      return;
    }
    applyState(st);
    $("runMsg").textContent = "Loaded config.";
  } catch (e) {
    $("runMsg").textContent = "Load failed.";
  }
}

// Events
$("btnLoadData").addEventListener("click", () => loadDatasetFromWindow($("selDataset").value));
$("selCond").addEventListener("change", () => { rebuildCondUnits(); rebuildFieldList(); });
$("inpSearch").addEventListener("input", rebuildFieldList);

$("btnAddField").addEventListener("click", () => {
  const sel = $("selFields");
  const opt = sel.selectedOptions && sel.selectedOptions[0];
  if (!opt) return;
  const field = opt.value;
  if (!selected.has(field)) selected.set(field, { unitChoice: null, valueText: "", uncText: "", weight: 1 });
  redrawSelectedFields();
});
$("btnClearFields").addEventListener("click", () => {
  selected.clear();
  redrawSelectedFields();
});
$("btnRun").addEventListener("click", () => { runAnalysis(); queueTelemetrySend(); });

$("btnSave").addEventListener("click", saveToLocalStorage);
$("btnLoadSaved").addEventListener("click", loadFromLocalStorage);

$("selCondUnit").addEventListener("change", () => {
  const u = $("selCondUnit").value;
  const t = $("inpBW").value.trim();
  if (!t) return;

  // If BW is a bare number, append unit; if it already has a unit, do nothing.
  if (/^[+-]?\d+(?:\.\d+)?(?:e[+-]?\d+)?$/i.test(t)) {
    $("inpBW").value = `${t} ${u}`;
  }
});


// Initial UI state
setEnabled(false);