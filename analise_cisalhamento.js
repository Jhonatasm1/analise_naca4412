"use strict";

/**
 * Analise do Centro de Cisalhamento (x_SC, y_SC) para o NACA 4412
 * tratado como casca de parede fina fechada.
 *
 * Convencoes:
 * - Origem em (0,0): bordo de ataque.
 * - Eixo x: ao longo da corda (m).
 * - Eixo y: vertical (m).
 * - Fluxo de cisalhamento q (N/m), para cargas Vx e Vy (N).
 */

const INPUT = Object.freeze({
  chord: 1.0,
  shellThickness: 0.002, // 2 mm
  pointCount: 100,
});

const NACA_4412 = Object.freeze({
  m: 0.04,
  p: 0.40,
  t: 0.12,
});

/**
 * Propriedades da Parte 1 (usadas diretamente no modelo de cisalhamento).
 */
const SECTION_PROPERTIES = Object.freeze({
  centroid: Object.freeze({
    x: 4.17406377e-1,
    y: 3.12234339e-2,
  }),
  inertia: Object.freeze({
    Ixx: 7.50460784e-5,
    Iyy: 4.47113694e-3,
    Ixy: 1.85235787e-5,
  }),
});

function cosineSpacing(pointCount, chord = 1) {
  if (!Number.isInteger(pointCount) || pointCount < 2) {
    throw new Error("pointCount deve ser inteiro >= 2.");
  }

  const xs = [];
  const step = Math.PI / (pointCount - 1);
  for (let i = 0; i < pointCount; i += 1) {
    const beta = i * step;
    xs.push(0.5 * chord * (1 - Math.cos(beta)));
  }

  return xs;
}

function thicknessDistribution(x, chord, t) {
  const xr = x / chord;
  return (
    5 *
    t *
    chord *
    (0.2969 * Math.sqrt(xr) -
      0.1260 * xr -
      0.3516 * xr ** 2 +
      0.2843 * xr ** 3 -
      0.1036 * xr ** 4)
  );
}

function camberLineAndSlope(x, chord, m, p) {
  const xr = x / chord;

  if (xr < p) {
    return {
      yC: (m / (p ** 2)) * (2 * p * xr - xr ** 2) * chord,
      dyCDx: (2 * m / (p ** 2)) * (p - xr),
    };
  }

  return {
    yC: (m / ((1 - p) ** 2)) * ((1 - 2 * p) + 2 * p * xr - xr ** 2) * chord,
    dyCDx: (2 * m / ((1 - p) ** 2)) * (p - xr),
  };
}

function signedAreaTwice(vertices) {
  let a2 = 0;
  for (let i = 0; i < vertices.length; i += 1) {
    const p0 = vertices[i];
    const p1 = vertices[(i + 1) % vertices.length];
    a2 += p0.x * p1.y - p1.x * p0.y;
  }
  return a2;
}

function ensureCounterClockwise(vertices) {
  return signedAreaTwice(vertices) < 0 ? vertices.slice().reverse() : vertices.slice();
}

/**
 * Gera o contorno fechado da secao com a mesma logica da Parte 1:
 * 100 pontos por superficie + fechamento via ultimo segmento.
 */
function generateNACA4412Contour(pointCount = 100, chord = 1) {
  const xs = cosineSpacing(pointCount, chord);
  const upper = [];
  const lower = [];

  for (const x of xs) {
    const { yC, dyCDx } = camberLineAndSlope(x, chord, NACA_4412.m, NACA_4412.p);
    const yT = thicknessDistribution(x, chord, NACA_4412.t);
    const theta = Math.atan(dyCDx);

    upper.push({
      x: x - yT * Math.sin(theta),
      y: yC + yT * Math.cos(theta),
    });

    lower.push({
      x: x + yT * Math.sin(theta),
      y: yC - yT * Math.cos(theta),
    });
  }

  const contour = [...upper.slice().reverse(), ...lower.slice(1)];
  return ensureCounterClockwise(contour);
}

function buildSegments(vertices) {
  const segments = [];

  for (let i = 0; i < vertices.length; i += 1) {
    const p0 = vertices[i];
    const p1 = vertices[(i + 1) % vertices.length];
    const dx = p1.x - p0.x;
    const dy = p1.y - p0.y;
    const ds = Math.hypot(dx, dy);

    if (ds <= 1e-14) {
      continue;
    }

    segments.push({
      p0,
      p1,
      dx,
      dy,
      ds,
      xMid: 0.5 * (p0.x + p1.x),
      yMid: 0.5 * (p0.y + p1.y),
    });
  }

  return segments;
}

function computePerimeter(segments) {
  return segments.reduce((acc, seg) => acc + seg.ds, 0);
}

/**
 * Calcula q_b, q_s0 e o torque interno para um par (Vx, Vy).
 *
 * q_b = -A*integral(y~ t_s ds) - B*integral(x~ t_s ds)
 * A = (Vy*Iyy - Vx*Ixy)/D
 * B = (Vx*Ixx - Vy*Ixy)/D
 * D = Ixx*Iyy - Ixy^2
 */
function solveShearFlowForLoad(segments, properties, thickness, Vx, Vy) {
  const { Ixx, Iyy, Ixy } = properties.inertia;
  const { x: cx, y: cy } = properties.centroid;

  const D = Ixx * Iyy - Ixy ** 2;
  if (Math.abs(D) < Number.EPSILON) {
    throw new Error("Denominador D ~ 0. Verifique Ixx, Iyy e Ixy.");
  }

  const A = (Vy * Iyy - Vx * Ixy) / D;
  const B = (Vx * Ixx - Vy * Ixy) / D;

  let intY = 0;
  let intX = 0;
  const qb = [];

  // Integracao incremental ao longo da linha media da casca.
  for (const seg of segments) {
    const xTildeMid = seg.xMid - cx;
    const yTildeMid = seg.yMid - cy;

    const dIntY = yTildeMid * thickness * seg.ds;
    const dIntX = xTildeMid * thickness * seg.ds;

    const intYMid = intY + 0.5 * dIntY;
    const intXMid = intX + 0.5 * dIntX;

    const qBasic = -A * intYMid - B * intXMid;
    qb.push(qBasic);

    intY += dIntY;
    intX += dIntX;
  }

  const perimeter = computePerimeter(segments);

  // Fechamento da celula por taxa de torcao nula.
  let numerator = 0;
  let denominator = 0;
  for (let i = 0; i < segments.length; i += 1) {
    numerator += (qb[i] / thickness) * segments[i].ds;
    denominator += (1 / thickness) * segments[i].ds;
  }

  const qS0 = -numerator / denominator;

  let torqueZ = 0;
  for (let i = 0; i < segments.length; i += 1) {
    const seg = segments[i];
    const qTotal = qb[i] + qS0;

    // Braço vetorial: r x ds => (x*dy - y*dx)
    const crossFactor = seg.xMid * seg.dy - seg.yMid * seg.dx;
    torqueZ += qTotal * crossFactor;
  }

  return {
    D,
    A,
    B,
    qb,
    qS0,
    torqueZ,
    perimeter,
  };
}

function computeShearCenter(segments, properties, thickness) {
  const unitVy = solveShearFlowForLoad(segments, properties, thickness, 0, 1);
  const unitVx = solveShearFlowForLoad(segments, properties, thickness, 1, 0);

  // T = x_SC*Vy - y_SC*Vx (momento externo equivalente sobre a origem).
  const xSC = unitVy.torqueZ;
  const ySC = -unitVx.torqueZ;

  return {
    xSC,
    ySC,
    unitVy,
    unitVx,
  };
}

function fmt(value, digits = 10) {
  return Number(value).toExponential(digits);
}

function printResults(summary, properties) {
  const { Ixx, Iyy, Ixy } = properties.inertia;
  const D = summary.unitVy.D;

  // Constantes de Vy e Vx nos dois termos de q_b.
  const AcoefVy = Iyy / D;
  const AcoefVx = -Ixy / D;
  const BcoefVx = Ixx / D;
  const BcoefVy = -Ixy / D;

  console.log("\n=== ANALISE DO CENTRO DE CISALHAMENTO | NACA 4412 ===");
  console.log("Parametros:");
  console.log(`  c = ${INPUT.chord}`);
  console.log(`  t_s = ${INPUT.shellThickness} m`);
  console.log(`  pontos por superficie = ${INPUT.pointCount}`);

  console.log("\n[1] Constantes multiplicadoras de V_y e V_x na equacao de q_b:");
  console.log(`  D = Ixx*Iyy - Ixy^2 = ${fmt(D)}`);
  console.log("  q_b = -A*integral(y~ t_s ds) - B*integral(x~ t_s ds)");
  console.log(`  A = (${fmt(AcoefVy)})*V_y + (${fmt(AcoefVx)})*V_x`);
  console.log(`  B = (${fmt(BcoefVx)})*V_x + (${fmt(BcoefVy)})*V_y`);

  console.log("\n[2] Perimetro total da celula:");
  console.log(`  L = ${fmt(summary.unitVy.perimeter)} m`);

  console.log("\n[3] Fluxo de fechamento q_s0:");
  console.log(`  q_s0(V_y) = (${fmt(summary.unitVy.qS0)}) * V_y  [com V_x = 0]`);
  console.log(
    `  referencia completa: q_s0 = (${fmt(summary.unitVy.qS0)})*V_y + (${fmt(summary.unitVx.qS0)})*V_x`
  );

  console.log("\n[4] Centro de Cisalhamento (em relacao ao bordo de ataque):");
  console.log(`  x_SC = ${fmt(summary.xSC)} m`);
  console.log(`  y_SC = ${fmt(summary.ySC)} m`);

  console.log("\n======================================================\n");
}

function runShearCenterAnalysis() {
  const contour = generateNACA4412Contour(INPUT.pointCount, INPUT.chord);
  const segments = buildSegments(contour);

  const summary = computeShearCenter(segments, SECTION_PROPERTIES, INPUT.shellThickness);
  printResults(summary, SECTION_PROPERTIES);

  return {
    input: INPUT,
    properties: SECTION_PROPERTIES,
    result: summary,
  };
}

// Execucao direta em Node ou no browser.
if (typeof window !== "undefined") {
  window.resultadoCisalhamentoNACA4412 = runShearCenterAnalysis();
} else {
  runShearCenterAnalysis();
}
