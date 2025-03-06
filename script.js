// Helper functions for vector operations
const vectorSubtract = (v1, v2) => [v1[0] - v2[0], v1[1] - v2[1]];
const vectorDot = (v1, v2) => v1[0] * v2[0] + v1[1] * v2[1];
const vectorCross = (v1, v2) => v1[0] * v2[1] - v1[1] * v2[0];
const vectorNormal = (v) => [-v[1], v[0]]; // 90-degree rotation (clockwise)
const vectorNormalize = (v) => {
  const magnitude = Math.sqrt(v[0] * v[0] + v[1] * v[1]);
  return [v[0] / magnitude, v[1] / magnitude];
};

/**
 * Determines which edges could be top facets of a polygon
 * @param {number[][]} polygon - Array of vertices in clockwise order
 * @returns {number[][][]} - Array of edges that could be top facets
 */
function getTopFacets(polygon) {
  const n = polygon.length;
  const topFacets = [];

  for (let i = 0; i < n; i++) {
    const edge = [polygon[i], polygon[(i + 1) % n]];
    const edgeVector = vectorSubtract(edge[1], edge[0]);
    const normal = vectorNormal(edgeVector);

    // Check if all other vertices lie on the same side of this edge
    let allSameSide = true;
    let sideSign = null; // Will be set on first check

    for (let j = 0; j < n; j++) {
      if (j !== i && j !== (i + 1) % n) {
        const vertexToCheck = polygon[j];
        const checkVector = vectorSubtract(vertexToCheck, edge[0]);
        const crossProduct = vectorCross(edgeVector, checkVector);

        // Determine which side the vertex is on
        if (sideSign === null) {
          sideSign = Math.sign(crossProduct);
        } else if (Math.sign(crossProduct) !== sideSign && crossProduct !== 0) {
          allSameSide = false;
          break;
        }
      }
    }

    if (allSameSide) {
      topFacets.push(edge);
    }
  }

  return topFacets;
}

/**
 * Rotates a polygon so the top facet is horizontal and all other edges are below it
 * @param {number[][]} polygon - Array of vertices in clockwise order
 * @param {number[][]} topEdge - The edge to be positioned as the top facet
 * @returns {number[][]} - The rotated polygon
 */
function rotateToEdge(polygon, topEdge) {
  // Check if edge is in polygon
  let edgeFound = false;
  for (let i = 0; i < polygon.length; i++) {
    const edge = [polygon[i], polygon[(i + 1) % polygon.length]];
    if (
      (edge[0][0] === topEdge[0][0] &&
        edge[0][1] === topEdge[0][1] &&
        edge[1][0] === topEdge[1][0] &&
        edge[1][1] === topEdge[1][1]) ||
      (edge[0][0] === topEdge[1][0] &&
        edge[0][1] === topEdge[1][1] &&
        edge[1][0] === topEdge[0][0] &&
        edge[1][1] === topEdge[0][1])
    ) {
      edgeFound = true;
      break;
    }
  }

  if (!edgeFound) {
    throw new Error("Edge is not in polygon");
  }

  // Calculate the angle to rotate
  const edgeVector = vectorSubtract(topEdge[1], topEdge[0]);
  let angle = Math.atan2(edgeVector[1], edgeVector[0]);

  // We want the edge to be horizontal, so rotate by -angle
  angle = -angle;

  // Rotate all vertices
  const rotatedPolygon = polygon.map((vertex) => {
    const x = vertex[0] * Math.cos(angle) - vertex[1] * Math.sin(angle);
    const y = vertex[0] * Math.sin(angle) + vertex[1] * Math.cos(angle);
    return [x, y];
  });

  // Translate so that the leftmost vertex of the edge is at the origin
  const leftVertex = edgeVector[0] <= 0 ? topEdge[0] : topEdge[1];
  const rotatedLeftVertex = rotatedPolygon.find(
    (v) =>
      Math.abs(
        v[0] -
          (leftVertex[0] * Math.cos(angle) - leftVertex[1] * Math.sin(angle))
      ) < 1e-10 &&
      Math.abs(
        v[1] -
          (leftVertex[0] * Math.sin(angle) + leftVertex[1] * Math.cos(angle))
      ) < 1e-10
  );

  const tx = -rotatedLeftVertex[0];
  const ty = -rotatedLeftVertex[1];

  return rotatedPolygon.map((v) => [v[0] + tx, v[1] + ty]);
}

/**
 * Gets half-plane constraints for the rotation center
 * @param {number[][]} polygon - Array of vertices (after rotation)
 * @param {number[][]} topFacet - The top facet (horizontal edge)
 * @returns {Object[]} - Array of half-plane constraints
 */
function getRotationConstraints(polygon, topFacet) {
  const constraints = [];
  const n = polygon.length;

  // Find the indices of the top facet vertices
  let topIndex1 = -1,
    topIndex2 = -1;
  for (let i = 0; i < n; i++) {
    if (
      Math.abs(polygon[i][0] - topFacet[0][0]) < 1e-10 &&
      Math.abs(polygon[i][1] - topFacet[0][1]) < 1e-10
    ) {
      topIndex1 = i;
    }
    if (
      Math.abs(polygon[i][0] - topFacet[1][0]) < 1e-10 &&
      Math.abs(polygon[i][1] - topFacet[1][1]) < 1e-10
    ) {
      topIndex2 = i;
    }
  }

  if (topIndex1 === -1 || topIndex2 === -1) {
    throw new Error("Top facet vertices not found in polygon");
  }

  // For each edge except the top facet
  for (let i = 0; i < n; i++) {
    if (
      (i !== topIndex1 || (i + 1) % n !== topIndex2) &&
      (i !== topIndex2 || (i + 1) % n !== topIndex1)
    ) {
      const edge = [polygon[i], polygon[(i + 1) % n]];
      const edgeVector = vectorSubtract(edge[1], edge[0]);
      let normal = vectorNormal(edgeVector);

      // We want the normal to point away from the mold (towards polygon center)
      // Since the polygon is in clockwise order, we invert the normal
      normal = [-normal[0], -normal[1]];

      // Create constraint based on normal direction
      if (normal[0] > 0) {
        // Region below half plane passing through bottom vertex
        constraints.push({
          type: "below",
          normal: normal,
          point: edge[0][1] < edge[1][1] ? edge[0] : edge[1], // Choose the bottom vertex
        });
      } else if (normal[0] < 0) {
        // Region above half plane passing through top vertex
        constraints.push({
          type: "above",
          normal: normal,
          point: edge[0][1] > edge[1][1] ? edge[0] : edge[1], // Choose the top vertex
        });
      }
      // If normal[0] === 0, the edge is vertical and doesn't constrain the rotation center
    }
  }

  return constraints;
}

/**
 * Solves the linear programming problem to find the rotation center
 * @param {Object[]} constraints - Array of half-plane constraints
 * @returns {number[]|null} - Coordinates of rotation center or null if none exists
 */
function solveLP(constraints) {
  // Define a bounding box for our search
  // In practice, we'd need a more sophisticated approach
  const BOUND = 1000;
  const candidates = [];

  // Check if a point satisfies all constraints
  const satisfiesConstraints = (point) => {
    for (const constraint of constraints) {
      const normal = constraint.normal;
      const constraintPoint = constraint.point;
      const vectorToPoint = vectorSubtract(point, constraintPoint);
      const dotProduct = vectorDot(normal, vectorToPoint);

      if (
        (constraint.type === "below" && dotProduct < 0) ||
        (constraint.type === "above" && dotProduct > 0)
      ) {
        return false;
      }
    }
    return true;
  };

  // Find intersection points of constraint lines
  for (let i = 0; i < constraints.length; i++) {
    for (let j = i + 1; j < constraints.length; j++) {
      const c1 = constraints[i];
      const c2 = constraints[j];

      // Calculate line equations: a1x + b1y = c1, a2x + b2y = c2
      const a1 = c1.normal[0];
      const b1 = c1.normal[1];
      const c1val = a1 * c1.point[0] + b1 * c1.point[1];

      const a2 = c2.normal[0];
      const b2 = c2.normal[1];
      const c2val = a2 * c2.point[0] + b2 * c2.point[1];

      // Check if lines are parallel
      const det = a1 * b2 - a2 * b1;
      if (Math.abs(det) < 1e-10) continue;

      // Find intersection
      const x = (c1val * b2 - c2val * b1) / det;
      const y = (a1 * c2val - a2 * c1val) / det;

      // Check if the intersection point satisfies all constraints
      if (satisfiesConstraints([x, y])) {
        candidates.push([x, y]);
      }
    }
  }

  // If no intersections found, check the extreme points in each direction
  // (This is a simplification; a real implementation would be more thorough)
  const extremePoints = [
    [-BOUND, -BOUND],
    [BOUND, -BOUND],
    [BOUND, BOUND],
    [-BOUND, BOUND],
  ];

  for (const point of extremePoints) {
    if (satisfiesConstraints(point)) {
      candidates.push(point);
    }
  }

  if (candidates.length === 0) {
    return null; // No valid rotation center
  }

  // Find the lowest point among valid candidates
  candidates.sort((a, b) => a[1] - b[1]);
  return candidates[0];
}

/**
 * Determines if a polygon is castable by clockwise rotation
 * @param {number[][]} polygon - Array of vertices in clockwise order
 * @param {number[][]} topFacet - The top facet
 * @returns {number[]|null} - Center of rotation or null if not castable
 */
function castabilityByRotation(polygon, topFacet) {
  // First, rotate the polygon so the top facet is horizontal
  const rotatedPolygon = rotateToEdge(polygon, topFacet);

  // Get the rotated top facet
  let rotatedTopFacet = null;
  for (let i = 0; i < rotatedPolygon.length; i++) {
    const edge = [
      rotatedPolygon[i],
      rotatedPolygon[(i + 1) % rotatedPolygon.length],
    ];
    if (Math.abs(edge[0][1] - edge[1][1]) < 1e-10 && edge[0][1] >= 0) {
      rotatedTopFacet = edge;
      break;
    }
  }

  if (!rotatedTopFacet) {
    throw new Error("Rotated top facet not found");
  }

  // Get constraints for rotation center
  const constraints = getRotationConstraints(rotatedPolygon, rotatedTopFacet);

  // Solve the linear programming problem
  const center = solveLP(constraints);

  return center;
}

// Main function to check if a polygon is castable
function isPolygonCastable(polygon) {
  const topFacets = getTopFacets(polygon);

  if (topFacets.length === 0) {
    return { castable: false, message: "No valid top facets found" };
  }

  // Try each potential top facet
  for (const topFacet of topFacets) {
    const center = castabilityByRotation(polygon, topFacet);
    if (center) {
      return {
        castable: true,
        center: center,
        topFacet: topFacet,
      };
    }
  }

  return { castable: false, message: "Polygon is not castable by rotation" };
}

// Visualization functions
function createCanvas(containerId, width, height) {
  const canvas = document.createElement("canvas");
  canvas.width = width;
  canvas.height = height;
  document.getElementById(containerId).appendChild(canvas);
  return canvas;
}

function scalePolygon(polygon, scale, offsetX, offsetY) {
  return polygon.map((p) => [p[0] * scale + offsetX, p[1] * scale + offsetY]);
}

function drawPolygon(
  ctx,
  polygon,
  color = "#2196F3",
  fillColor = "rgba(33, 150, 243, 0.2)"
) {
  ctx.beginPath();
  ctx.moveTo(polygon[0][0], polygon[0][1]);

  for (let i = 1; i < polygon.length; i++) {
    ctx.lineTo(polygon[i][0], polygon[i][1]);
  }

  ctx.closePath();
  ctx.fillStyle = fillColor;
  ctx.fill();
  ctx.strokeStyle = color;
  ctx.lineWidth = 2;
  ctx.stroke();
}

function drawEdge(ctx, edge, color = "#4CAF50", lineWidth = 3) {
  ctx.beginPath();
  ctx.moveTo(edge[0][0], edge[0][1]);
  ctx.lineTo(edge[1][0], edge[1][1]);
  ctx.strokeStyle = color;
  ctx.lineWidth = lineWidth;
  ctx.stroke();
}

function drawRotationCenter(ctx, center, color = "#F44336") {
  ctx.beginPath();
  ctx.arc(center[0], center[1], 5, 0, Math.PI * 2);
  ctx.fillStyle = color;
  ctx.fill();
  ctx.strokeStyle = "#000";
  ctx.lineWidth = 1;
  ctx.stroke();
}

function drawRotationPath(
  ctx,
  polygon,
  center,
  color = "rgba(255, 152, 0, 0.3)"
) {
  for (const vertex of polygon) {
    // Calculate radius (distance from center to vertex)
    const dx = vertex[0] - center[0];
    const dy = vertex[1] - center[1];
    const radius = Math.sqrt(dx * dx + dy * dy);

    // Calculate start angle
    const startAngle = Math.atan2(dy, dx);

    // Draw arc (clockwise rotation)
    ctx.beginPath();
    ctx.arc(
      center[0],
      center[1],
      radius,
      startAngle,
      startAngle - Math.PI / 2,
      true
    );
    ctx.strokeStyle = color;
    ctx.lineWidth = 2;
    ctx.stroke();
  }
}

function visualizePolygon(testCase) {
  const containerDiv = document.createElement("div");
  containerDiv.className = "test-container";
  containerDiv.innerHTML = `<h2>${testCase.name}</h2>`;

  document.getElementById("testContainer").appendChild(containerDiv);

  // Create a canvas for this test case
  const canvasId = `canvas-${testCase.name.replace(/\s+/g, "-").toLowerCase()}`;
  containerDiv.innerHTML += `<canvas id="${canvasId}" width="600" height="400"></canvas>`;

  const canvas = document.getElementById(canvasId);
  const ctx = canvas.getContext("2d");

  // Clear canvas
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  // Calculate scale and offset to center the polygon
  let minX = Infinity,
    minY = Infinity,
    maxX = -Infinity,
    maxY = -Infinity;
  for (const point of testCase.polygon) {
    minX = Math.min(minX, point[0]);
    minY = Math.min(minY, point[1]);
    maxX = Math.max(maxX, point[0]);
    maxY = Math.max(maxY, point[1]);
  }

  const polygonWidth = maxX - minX;
  const polygonHeight = maxY - minY;
  const scale = Math.min(
    (canvas.width - 100) / polygonWidth,
    (canvas.height - 100) / polygonHeight
  );
  const offsetX = (canvas.width - polygonWidth * scale) / 2 - minX * scale;
  const offsetY = (canvas.height - polygonHeight * scale) / 2 - minY * scale;

  // Scale the polygon to fit the canvas
  const scaledPolygon = scalePolygon(testCase.polygon, scale, offsetX, offsetY);

  // Get castability result
  let result;
  try {
    result = isPolygonCastable(testCase.polygon);

    // Draw the polygon
    drawPolygon(ctx, scaledPolygon);

    // Display result info
    const resultInfo = document.createElement("div");
    resultInfo.className = "result-info";

    if (result.castable) {
      resultInfo.innerHTML = `<span class="castable">✓ Castable</span> - Rotation center: (${result.center[0].toFixed(
        2
      )}, ${result.center[1].toFixed(2)})`;

      // Scale the top facet
      const scaledTopFacet = [
        [
          result.topFacet[0][0] * scale + offsetX,
          result.topFacet[0][1] * scale + offsetY,
        ],
        [
          result.topFacet[1][0] * scale + offsetX,
          result.topFacet[1][1] * scale + offsetY,
        ],
      ];

      // Draw the top facet
      drawEdge(ctx, scaledTopFacet, "#4CAF50");

      // Scale the rotation center
      const scaledCenter = [
        result.center[0] * scale + offsetX,
        result.center[1] * scale + offsetY,
      ];

      // Draw the rotation center
      drawRotationCenter(ctx, scaledCenter);

      // Draw rotation paths
      drawRotationPath(ctx, scaledPolygon, scaledCenter);
    } else {
      resultInfo.innerHTML = `<span class="not-castable">✗ Not Castable</span> - ${result.message}`;
    }

    containerDiv.appendChild(resultInfo);

    // Add facets info
    const topFacets = getTopFacets(testCase.polygon);
    const facetsInfo = document.createElement("div");
    facetsInfo.className = "result-info";
    facetsInfo.textContent = `Found ${topFacets.length} potential top facet(s)`;
    containerDiv.appendChild(facetsInfo);
  } catch (error) {
    console.error(`Error visualizing ${testCase.name}:`, error.message);
    const errorInfo = document.createElement("div");
    errorInfo.className = "result-info not-castable";
    errorInfo.textContent = `Error: ${error.message}`;
    containerDiv.appendChild(errorInfo);
  }
}

// Test cases
const testCases = [
  {
    name: "Convex Pentagon",
    polygon: [
      [0, 0],
      [5, 0],
      [5, 3],
      [3, 5],
      [0, 3],
    ],
  },
  {
    name: "Square",
    polygon: [
      [0, 0],
      [5, 0],
      [5, 5],
      [0, 5],
    ],
  },
  {
    name: "Triangle",
    polygon: [
      [0, 0],
      [5, 0],
      [2, 4],
    ],
  },
  {
    name: "L-Shape",
    polygon: [
      [0, 0],
      [3, 0],
      [3, 1],
      [1, 1],
      [1, 3],
      [0, 3],
    ],
  },
  {
    name: "Star Shape",
    polygon: [
      [3, 0],
      [4, 2],
      [6, 2],
      [5, 3],
      [6, 5],
      [3, 4],
      [0, 5],
      [1, 3],
      [0, 2],
      [2, 2],
    ],
  },
  {
    name: "U-Shape",
    polygon: [
      [0, 0],
      [5, 0],
      [5, 5],
      [4, 5],
      [4, 1],
      [1, 1],
      [1, 5],
      [0, 5],
    ],
  },
  {
    name: "Arrow Shape",
    polygon: [
      [0, 0],
      [2, 0],
      [3, 1],
      [4, 0],
      [6, 0],
      [3, 3],
    ],
  },
];

// Visualize each test case
testCases.forEach(visualizePolygon);
