"use strict";

/**
 * Parametros do NACA 4412 (serie de 4 digitos):
 * m = 0.04 (camber maximo), p = 0.40 (posicao do camber maximo), t = 0.12 (espessura relativa).
 */
const NACA_4412 = Object.freeze({
  m: 0.04,
  p: 0.40,
  t: 0.12,
});

let airfoilChart = null;

/**
 * Gera pontos x no intervalo [0, c] com espacamento em cosseno.
 * O clustering proximo ao bordo de ataque melhora a precisao geometrica.
 */
function cosineSpacing(pointCount, chord = 1) {
  if (!Number.isInteger(pointCount) || pointCount < 2) {
    throw new Error("pointCount deve ser inteiro >= 2.");
  }
  if (!(chord > 0)) {
    throw new Error("chord deve ser positivo.");
  }

  const betaStep = Math.PI / (pointCount - 1);
  const xCoordinates = [];

  for (let i = 0; i < pointCount; i += 1) {
    const beta = i * betaStep;
    const x = 0.5 * chord * (1 - Math.cos(beta));
    xCoordinates.push(x);
  }

  return xCoordinates;
}

/**
 * Distribuicao de espessura (y_t) para perfis NACA 4-digitos.
 * Foi usado -0.1036 no termo de x^4 para fechar o bordo de fuga.
 */
function thicknessDistribution(x, chord, t) {
  const xOverC = x / chord;
  const yT =
    5 *
    t *
    chord *
    (0.2969 * Math.sqrt(xOverC) -
      0.1260 * xOverC -
      0.3516 * xOverC ** 2 +
      0.2843 * xOverC ** 3 -
      0.1036 * xOverC ** 4);

  return yT;
}

/**
 * Calcula linha media (y_c) e inclinacao (dy_c/dx) no ponto x.
 */
function camberLineAndSlope(x, chord, m, p) {
  const xOverC = x / chord;

  if (xOverC < p) {
    const yC = (m / (p ** 2)) * (2 * p * xOverC - xOverC ** 2) * chord;
    const dyCDx = (2 * m / (p ** 2)) * (p - xOverC);
    return { yC, dyCDx };
  }

  const yC = (m / ((1 - p) ** 2)) * ((1 - 2 * p) + 2 * p * xOverC - xOverC ** 2) * chord;
  const dyCDx = (2 * m / ((1 - p) ** 2)) * (p - xOverC);
  return { yC, dyCDx };
}

/**
 * Gera coordenadas discretizadas do NACA 4412.
 * pointCount representa o numero de pontos por superficie (superior e inferior).
 */
function generateNACA4412Coordinates(pointCount = 100, chord = 1) {
  const xList = cosineSpacing(pointCount, chord);
  const upper = [];
  const lower = [];

  for (const x of xList) {
    const { yC, dyCDx } = camberLineAndSlope(x, chord, NACA_4412.m, NACA_4412.p);
    const yT = thicknessDistribution(x, chord, NACA_4412.t);
    const theta = Math.atan(dyCDx);

    // Rotacao local da espessura ao redor da linha media.
    const xUpper = x - yT * Math.sin(theta);
    const yUpper = yC + yT * Math.cos(theta);
    const xLower = x + yT * Math.sin(theta);
    const yLower = yC - yT * Math.cos(theta);

    upper.push({ x: xUpper, y: yUpper });
    lower.push({ x: xLower, y: yLower });
  }

  /**
   * Construcao de poligono fechado:
   * 1) superficie superior de TE para LE
   * 2) superficie inferior de LE para TE (evitando repetir o ponto LE)
   */
  const polygon = [...upper.slice().reverse(), ...lower.slice(1)];
  return { upper, lower, polygon };
}

/**
 * Soma do termo cruzado (2 * area assinada).
 */
function signedAreaTwice(vertices) {
  let area2 = 0;

  for (let i = 0; i < vertices.length; i += 1) {
    const a = vertices[i];
    const b = vertices[(i + 1) % vertices.length];
    area2 += a.x * b.y - b.x * a.y;
  }

  return area2;
}

/**
 * Garante orientacao anti-horaria para manter area positiva nas formulas.
 */
function ensureCounterClockwise(vertices) {
  const area2 = signedAreaTwice(vertices);
  return area2 < 0 ? vertices.slice().reverse() : vertices.slice();
}

/**
 * Aplica o Teorema de Green em um poligono simples e fechado.
 * Retorna area, centroide e inercia em relacao ao centroide.
 */
function polygonSectionProperties(vertices) {
  if (!Array.isArray(vertices) || vertices.length < 3) {
    throw new Error("Sao necessarios pelo menos 3 vertices para formar um poligono.");
  }

  const v = ensureCounterClockwise(vertices);

  let area2 = 0;
  let cxAcc = 0;
  let cyAcc = 0;
  let ixxAcc = 0;
  let iyyAcc = 0;
  let ixyAcc = 0;

  for (let i = 0; i < v.length; i += 1) {
    const p0 = v[i];
    const p1 = v[(i + 1) % v.length];

    const cross = p0.x * p1.y - p1.x * p0.y;
    area2 += cross;

    // Termos do centroide via formula de poligono.
    cxAcc += (p0.x + p1.x) * cross;
    cyAcc += (p0.y + p1.y) * cross;

    // Segundos momentos em relacao a origem global.
    ixxAcc += (p0.y ** 2 + p0.y * p1.y + p1.y ** 2) * cross;
    iyyAcc += (p0.x ** 2 + p0.x * p1.x + p1.x ** 2) * cross;

    // Produto de inercia em relacao a origem global.
    ixyAcc += (p0.x * p1.y + 2 * p0.x * p0.y + 2 * p1.x * p1.y + p1.x * p0.y) * cross;
  }

  const area = area2 / 2;
  if (Math.abs(area) < Number.EPSILON) {
    throw new Error("Area numerica nula. Verifique os vertices do poligono.");
  }

  const cx = cxAcc / (3 * area2);
  const cy = cyAcc / (3 * area2);

  const ixxOrigin = ixxAcc / 12;
  const iyyOrigin = iyyAcc / 12;
  const ixyOrigin = ixyAcc / 24;

  // Teorema dos eixos paralelos para mover da origem para o centroide.
  const ixxCentroid = ixxOrigin - area * cy ** 2;
  const iyyCentroid = iyyOrigin - area * cx ** 2;
  const ixyCentroid = ixyOrigin - area * cx * cy;

  return {
    area,
    centroid: { x: cx, y: cy },
    inertiaAtOrigin: { ixx: ixxOrigin, iyy: iyyOrigin, ixy: ixyOrigin },
    inertiaAtCentroid: { ixx: ixxCentroid, iyy: iyyCentroid, ixy: ixyCentroid },
  };
}

/**
 * Pipeline completo para o NACA 4412.
 */
function analyzeNACA4412(pointCount = 100, chord = 1) {
  const geometry = generateNACA4412Coordinates(pointCount, chord);
  const properties = polygonSectionProperties(geometry.polygon);

  return {
    input: { pointCount, chord },
    geometry,
    properties,
  };
}

function formatNumber(value) {
  return Number(value).toExponential(8);
}

/**
 * Atualiza metricas tanto por id fixo quanto por data-metric para permitir
 * reuso do mesmo valor em cards e tabelas sem duplicar logica.
 */
function setMetricValue(metricName, formattedValue) {
  const metricIdMap = {
    area: "areaValue",
    cx: "cxValue",
    cy: "cyValue",
    ixx: "ixxValue",
    iyy: "iyyValue",
    ixy: "ixyValue",
  };

  const metricId = metricIdMap[metricName];
  if (metricId) {
    const byId = document.getElementById(metricId);
    if (byId) {
      byId.textContent = formattedValue;
    }
  }

  const byAttribute = document.querySelectorAll(`[data-metric="${metricName}"]`);
  byAttribute.forEach((node) => {
    node.textContent = formattedValue;
  });
}

function getCoordinatePreview(points, limit = 6) {
  return points.slice(0, limit).map((pt, index) => {
    return `${index + 1}. x=${pt.x.toFixed(6)}, y=${pt.y.toFixed(6)}`;
  });
}

function renderResults(result) {
  const { area, centroid, inertiaAtCentroid } = result.properties;

  setMetricValue("area", formatNumber(area));
  setMetricValue("cx", formatNumber(centroid.x));
  setMetricValue("cy", formatNumber(centroid.y));
  setMetricValue("ixx", formatNumber(inertiaAtCentroid.ixx));
  setMetricValue("iyy", formatNumber(inertiaAtCentroid.iyy));
  setMetricValue("ixy", formatNumber(inertiaAtCentroid.ixy));

  const upperPreview = getCoordinatePreview(result.geometry.upper).join("\n");
  const lowerPreview = getCoordinatePreview(result.geometry.lower).join("\n");
  const sampleText = `SUPERFICIE SUPERIOR (primeiros pontos)\n${upperPreview}\n\nSUPERFICIE INFERIOR (primeiros pontos)\n${lowerPreview}`;
  document.getElementById("sampleCoordinates").textContent = sampleText;
}

/**
 * Renderiza o perfil aerodinamico no canvas usando Chart.js.
 */
function renderAirfoilChart(geometry) {
  const canvas = document.getElementById("airfoilChart");
  if (!canvas || typeof Chart === "undefined") {
    return;
  }

  const gridStep = 0.1;
  const tickFormatter = (value) => {
    const safeValue = Math.abs(value) < 1e-12 ? 0 : value;
    return safeValue.toFixed(1);
  };

  const upperData = geometry.upper.map((point) => ({ x: point.x, y: point.y }));
  const lowerData = geometry.lower.map((point) => ({ x: point.x, y: point.y }));

  if (airfoilChart) {
    airfoilChart.destroy();
  }

  airfoilChart = new Chart(canvas, {
    type: "scatter",
    data: {
      datasets: [
        {
          label: "Extradorso",
          data: upperData,
          showLine: true,
          pointRadius: 0,
          borderWidth: 2,
          borderColor: "#38bdf8",
        },
        {
          label: "Intradorso",
          data: lowerData,
          showLine: true,
          pointRadius: 0,
          borderWidth: 2,
          borderColor: "#f59e0b",
        },
      ],
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      animation: {
        duration: 650,
      },
      plugins: {
        legend: {
          labels: {
            color: "#e2e8f0",
          },
        },
      },
      scales: {
        x: {
          type: "linear",
          min: 0,
          max: 1.0,
          title: {
            display: true,
            text: "X",
            color: "#e2e8f0",
          },
          grid: {
            color: "rgba(148, 163, 184, 0.18)",
            borderDash: [2, 3],
          },
          ticks: {
            color: "#e2e8f0",
            stepSize: gridStep,
            callback: tickFormatter,
          },
        },
        y: {
          min: -0.4,
          max: 0.4,
          title: {
            display: true,
            text: "Y",
            color: "#e2e8f0",
          },
          grid: {
            color: (context) => {
              return context.tick.value === 0 ? "rgba(226, 232, 240, 0.85)" : "rgba(148, 163, 184, 0.18)";
            },
            borderDash: (context) => {
              return context.tick.value === 0 ? [] : [2, 3];
            },
            lineWidth: (context) => {
              return context.tick.value === 0 ? 1.6 : 1;
            },
          },
          ticks: {
            color: "#e2e8f0",
            stepSize: gridStep,
            callback: tickFormatter,
          },
        },
      },
    },
  });
}

function runAnalysis() {
  const result = analyzeNACA4412(100, 1);
  renderResults(result);
  renderAirfoilChart(result.geometry);

  // Facilita comparacao com CAD: dados ficam acessiveis no console.
  window.naca4412Result = result;
  return result;
}

document.addEventListener("DOMContentLoaded", () => {
  runAnalysis();

  const button = document.getElementById("recalculateBtn");
  if (button) {
    button.addEventListener("click", () => {
      runAnalysis();
    });
  }
});
