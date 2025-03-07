// ==========================================
// RANDOM POLYGON GENERATION
// ==========================================

/**
 * Checks if a line segment intersects with another line segment
 * @param {number[]} p1 First point of first segment
 * @param {number[]} p2 Second point of first segment
 * @param {number[]} p3 First point of second segment
 * @param {number[]} p4 Second point of second segment
 * @returns {boolean} True if the segments intersect
 */
function doSegmentsIntersect(p1, p2, p3, p4) {
  // Calculate the direction vectors
  const d1 = [p2[0] - p1[0], p2[1] - p1[1]];
  const d2 = [p4[0] - p3[0], p4[1] - p3[1]];
  const d3 = [p3[0] - p1[0], p3[1] - p1[1]];

  // Calculate the cross products
  const cross1 = d1[0] * d3[1] - d1[1] * d3[0]; // Cross product of d1 and d3
  const cross2 = d1[0] * d2[1] - d1[1] * d2[0]; // Cross product of d1 and d2

  if (Math.abs(cross2) < 1e-10) {
    // Lines are parallel or collinear
    return false;
  }

  // Calculate the parameter t for the intersection point on the second line segment
  const t1 = cross1 / cross2;

  if (t1 < 0 || t1 > 1) {
    // Intersection point is not on the second segment
    return false;
  }

  // Calculate the parameter s for the intersection point on the first line segment
  const cross3 = d3[0] * d2[1] - d3[1] * d2[0]; // Cross product of d3 and d2
  const s = cross3 / cross2;

  if (s < 0 || s > 1) {
    // Intersection point is not on the first segment
    return false;
  }

  return true;
}

/**
 * Checks if a polygon is simple (no self-intersections)
 * @param {number[][]} polygon Array of vertices
 * @returns {boolean} True if the polygon is simple
 */
function isSimplePolygon(polygon) {
  const n = polygon.length;

  // Check for intersections between non-adjacent edges
  for (let i = 0; i < n; i++) {
    const p1 = polygon[i];
    const p2 = polygon[(i + 1) % n];

    for (let j = i + 2; j < n; j++) {
      // Skip adjacent edges (they share a vertex)
      if (i === 0 && j === n - 1) continue;

      const p3 = polygon[j];
      const p4 = polygon[(j + 1) % n];

      if (doSegmentsIntersect(p1, p2, p3, p4)) {
        return false;
      }
    }
  }

  return true;
}

/**
 * Generate a random polygon
 * @param {number} numVertices Number of vertices
 * @param {number} minSize Minimum size
 * @param {number} maxSize Maximum size
 * @returns {number[][]} Array of vertices
 */
function generateRandomPolygon(numVertices, minSize = 50, maxSize = 150) {
  // Generate random convex polygon
  if (Math.random() < 0.5) {
    return generateRandomConvexPolygon(numVertices, minSize, maxSize);
  } else {
    return generateRandomNonConvexPolygon(numVertices, minSize, maxSize);
  }
}

/**
 * Generate a random convex polygon
 * @param {number} numVertices Number of vertices
 * @param {number} minSize Minimum size
 * @param {number} maxSize Maximum size
 * @returns {number[][]} Array of vertices
 */
function generateRandomConvexPolygon(numVertices, minSize = 50, maxSize = 150) {
  // Generate random points
  const angles = [];
  for (let i = 0; i < numVertices; i++) {
    angles.push(Math.random() * 2 * Math.PI);
  }

  // Sort angles to ensure convexity
  angles.sort();

  const centerX = 0;
  const centerY = 0;
  const radius = minSize + Math.random() * (maxSize - minSize);

  // Create vertices using polar coordinates
  const vertices = [];
  for (let i = 0; i < numVertices; i++) {
    const angle = angles[i];
    // Add some randomness to the radius for more varied shapes
    const r = radius * (0.5 + 0.5 * Math.random());
    const x = centerX + r * Math.cos(angle);
    const y = centerY + r * Math.sin(angle);
    vertices.push([x, y]);
  }

  // Convex polygons created this way are guaranteed to be simple
  return vertices;
}

/**
 * Generate a random non-convex polygon
 * @param {number} numVertices Number of vertices
 * @param {number} minSize Minimum size
 * @param {number} maxSize Maximum size
 * @returns {number[][]} Array of vertices
 */
function generateRandomNonConvexPolygon(
  numVertices,
  minSize = 50,
  maxSize = 150
) {
  // First generate a convex polygon
  const convexPolygon = generateRandomConvexPolygon(
    numVertices,
    minSize,
    maxSize
  );

  // Calculate centroid
  let centroidX = 0,
    centroidY = 0;
  for (const vertex of convexPolygon) {
    centroidX += vertex[0];
    centroidY += vertex[1];
  }
  centroidX /= numVertices;
  centroidY /= numVertices;

  // Make multiple attempts to create a simple non-convex polygon
  for (let attempt = 0; attempt < 10; attempt++) {
    // Make a copy of the convex polygon
    const nonConvexPolygon = [...convexPolygon];

    // Number of vertices to push inward (at least 1, at most half)
    // Use fewer vertices for more complex polygons to reduce chance of self-intersection
    const maxVerticesToModify = Math.max(1, Math.floor(numVertices / 3));
    const numToModify = 1 + Math.floor(Math.random() * maxVerticesToModify);

    // Randomly select vertices to modify, but avoid selecting adjacent vertices
    // to reduce the chance of self-intersections
    const indicesToModify = new Set();
    const avoidIndices = new Set();

    while (
      indicesToModify.size < numToModify &&
      indicesToModify.size + avoidIndices.size < numVertices
    ) {
      const index = Math.floor(Math.random() * numVertices);

      // Skip if this index or adjacent indices should be avoided
      if (avoidIndices.has(index)) continue;

      indicesToModify.add(index);

      // Avoid modifying adjacent vertices
      avoidIndices.add(index);
      avoidIndices.add((index + 1) % numVertices);
      avoidIndices.add((index - 1 + numVertices) % numVertices);
    }

    // Start with smaller move factors to reduce chance of self-intersection
    const maxMoveFactor = 0.4;

    // Push selected vertices inward with controlled movement
    for (const index of indicesToModify) {
      const vertex = nonConvexPolygon[index];
      const vectorToCentroid = [centroidX - vertex[0], centroidY - vertex[1]];
      const magnitude = Math.sqrt(
        vectorToCentroid[0] ** 2 + vectorToCentroid[1] ** 2
      );
      const normalizedVector = [
        vectorToCentroid[0] / magnitude,
        vectorToCentroid[1] / magnitude,
      ];

      // Use a smaller movement factor to reduce self-intersections
      // Start with small movements
      const moveFactor = 0.1 + Math.random() * maxMoveFactor;

      // Calculate the new position
      const newX = vertex[0] + normalizedVector[0] * magnitude * moveFactor;
      const newY = vertex[1] + normalizedVector[1] * magnitude * moveFactor;

      // Update the vertex
      nonConvexPolygon[index] = [newX, newY];
    }

    // Check if the resulting polygon is simple (no self-intersections)
    if (isSimplePolygon(nonConvexPolygon)) {
      return nonConvexPolygon;
    }
  }

  // If all attempts fail, fall back to the original convex polygon
  console.log(
    "Failed to generate simple non-convex polygon, falling back to convex"
  );
  return convexPolygon;
}

/**
 * Function to generate a set of random test cases
 * @param {number} count Number of test cases to generate
 * @param {number} minVertices Minimum number of vertices
 * @param {number} maxVertices Maximum number of vertices
 * @returns {Array} Array of test cases
 */
function generateRandomTestCases(count = 5, minVertices = 5, maxVertices = 10) {
  const randomTestCases = [];

  for (let i = 0; i < count; i++) {
    const numVertices =
      minVertices + Math.floor(Math.random() * (maxVertices - minVertices + 1));
    const polygon = generateRandomPolygon(numVertices);
    randomTestCases.push({
      name: `Random Polygon ${i + 1}`,
      polygon: polygon,
    });
  }

  return randomTestCases;
} // ==========================================
// VECTOR OPERATIONS AND HELPER FUNCTIONS
// ==========================================
const vectorSubtract = (v1, v2) => [v1[0] - v2[0], v1[1] - v2[1]];
const vectorDot = (v1, v2) => v1[0] * v2[0] + v1[1] * v2[1];
const vectorCross = (v1, v2) => v1[0] * v2[1] - v1[1] * v2[0];
const vectorNormal = (v) => [-v[1], v[0]]; // 90-degree rotation (clockwise)
const vectorNormalize = (v) => {
  const magnitude = Math.sqrt(v[0] * v[0] + v[1] * v[1]);
  return [v[0] / magnitude, v[1] / magnitude];
};

// Helper function to check if a polygon is in clockwise order
function isClockwise(polygon) {
  let sum = 0;
  for (let i = 0; i < polygon.length; i++) {
    const current = polygon[i];
    const next = polygon[(i + 1) % polygon.length];
    sum += (next[0] - current[0]) * (next[1] + current[1]);
  }
  return sum > 0;
}

// Ensure polygon is in clockwise order
function ensureClockwise(polygon) {
  if (!isClockwise(polygon)) {
    return polygon.slice().reverse();
  }
  return polygon;
}

// ==========================================
// CORE POLYGON CASTABILITY ALGORITHM
// ==========================================

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

  // Find the rightmost vertex in the polygon
  let rightmostX = -Infinity;
  let rightmostVertex = null;

  for (const vertex of polygon) {
    if (vertex[0] > rightmostX) {
      rightmostX = vertex[0];
      rightmostVertex = vertex;
    }
  }

  // Add initial constraint: rotation center must be to the right of the rightmost vertex
  constraints.push({
    type: "right",
    normal: [1, 0], // Points to the right
    point: rightmostVertex, // The rightmost vertex
  });

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
        (constraint.type === "above" && dotProduct > 0) ||
        (constraint.type === "right" && dotProduct <= 0)
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

// ==========================================
// RANDOM POLYGON GENERATION
// ==========================================

/**
 * Checks if a line segment intersects with another line segment
 * @param {number[]} p1 First point of first segment
 * @param {number[]} p2 Second point of first segment
 * @param {number[]} p3 First point of second segment
 * @param {number[]} p4 Second point of second segment
 * @returns {boolean} True if the segments intersect
 */
function doSegmentsIntersect(p1, p2, p3, p4) {
  // Calculate the direction vectors
  const d1 = [p2[0] - p1[0], p2[1] - p1[1]];
  const d2 = [p4[0] - p3[0], p4[1] - p3[1]];
  const d3 = [p3[0] - p1[0], p3[1] - p1[1]];

  // Calculate the cross products
  const cross1 = d1[0] * d3[1] - d1[1] * d3[0]; // Cross product of d1 and d3
  const cross2 = d1[0] * d2[1] - d1[1] * d2[0]; // Cross product of d1 and d2

  if (Math.abs(cross2) < 1e-10) {
    // Lines are parallel or collinear
    return false;
  }

  // Calculate the parameter t for the intersection point on the second line segment
  const t1 = cross1 / cross2;

  if (t1 < 0 || t1 > 1) {
    // Intersection point is not on the second segment
    return false;
  }

  // Calculate the parameter s for the intersection point on the first line segment
  const cross3 = d3[0] * d2[1] - d3[1] * d2[0]; // Cross product of d3 and d2
  const s = cross3 / cross2;

  if (s < 0 || s > 1) {
    // Intersection point is not on the first segment
    return false;
  }

  return true;
}

/**
 * Checks if a polygon is simple (no self-intersections)
 * @param {number[][]} polygon Array of vertices
 * @returns {boolean} True if the polygon is simple
 */
function isSimplePolygon(polygon) {
  const n = polygon.length;

  // Check for intersections between non-adjacent edges
  for (let i = 0; i < n; i++) {
    const p1 = polygon[i];
    const p2 = polygon[(i + 1) % n];

    for (let j = i + 2; j < n; j++) {
      // Skip adjacent edges (they share a vertex)
      if (i === 0 && j === n - 1) continue;

      const p3 = polygon[j];
      const p4 = polygon[(j + 1) % n];

      if (doSegmentsIntersect(p1, p2, p3, p4)) {
        return false;
      }
    }
  }

  return true;
}

// Function to generate a set of random test cases
function generateRandomTestCases(count = 5, minVertices = 5, maxVertices = 10) {
  const randomTestCases = [];

  for (let i = 0; i < count; i++) {
    const numVertices =
      minVertices + Math.floor(Math.random() * (maxVertices - minVertices + 1));
    const polygon = generateRandomPolygon(numVertices);
    randomTestCases.push({
      name: `Random Polygon ${i + 1}`,
      polygon: polygon,
    });
  }

  return randomTestCases;
}

// ==========================================
// VISUALIZATION FUNCTIONS
// ==========================================

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

    // Draw complete circular arc (full 360 degrees)
    ctx.beginPath();
    ctx.arc(center[0], center[1], radius, 0, Math.PI * 2);
    ctx.strokeStyle = color;
    ctx.lineWidth = 1;
    ctx.stroke();

    // Draw the initial vertex position
    ctx.beginPath();
    ctx.arc(vertex[0], vertex[1], 3, 0, Math.PI * 2);
    ctx.fillStyle = color;
    ctx.fill();

    // Draw multiple positions along the rotation path to illustrate movement
    const numPoints = 12; // Number of points to display along the path
    for (let i = 1; i <= numPoints; i++) {
      const angle = startAngle - (Math.PI * 2 * i) / numPoints;
      const x = center[0] + radius * Math.cos(angle);
      const y = center[1] + radius * Math.sin(angle);

      // Draw a small dot at each position
      ctx.beginPath();
      ctx.arc(x, y, 1.5, 0, Math.PI * 2);
      ctx.fillStyle =
        i === numPoints / 4
          ? "rgba(255, 100, 0, 0.7)"
          : "rgba(255, 152, 0, 0.5)";
      ctx.fill();
    }
  }
}

function visualizePolygon(testCase, containerId = "testContainer") {
  const containerDiv = document.createElement("div");
  containerDiv.className = "test-container";
  containerDiv.innerHTML = `<h2>${testCase.name}</h2>`;

  document.getElementById(containerId).appendChild(containerDiv);

  // Create a canvas for this test case
  const canvasId = `canvas-${testCase.name.replace(/\s+/g, "-").toLowerCase()}`;
  containerDiv.innerHTML += `<canvas id="${canvasId}" width="600" height="400"></canvas>`;

  const canvas = document.getElementById(canvasId);
  const ctx = canvas.getContext("2d");

  // Clear canvas
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  // Get castability result first to find the top facet
  let result;
  try {
    result = isPolygonCastable(testCase.polygon);

    let scaledPolygon;
    let scaledCenter = null;
    let scaledTopFacet = null;
    let rotatedPolygon = null; // Declare rotatedPolygon here

    // Calculate scale and offset to center the polygon
    let polygonToDisplay = [...testCase.polygon]; // Clone the polygon

    // If castable, properly orient the polygon so the top facet is at the top
    if (result.castable) {
      // Step 1: Rotate the polygon so that top facet is horizontal
      rotatedPolygon = rotateToEdge(testCase.polygon, result.topFacet);

      // Step 2: Find the y-coordinate of the top facet
      let topFacetY = null;
      for (let i = 0; i < rotatedPolygon.length; i++) {
        const v1 = rotatedPolygon[i];
        const v2 = rotatedPolygon[(i + 1) % rotatedPolygon.length];

        // Check if this edge corresponds to the rotated top facet
        // (it should be a horizontal line)
        if (Math.abs(v1[1] - v2[1]) < 1e-10) {
          topFacetY = v1[1];
          break;
        }
      }

      if (topFacetY === null) {
        console.error(
          "Could not find the top facet y-coordinate in rotated polygon"
        );
        topFacetY = 0;
      }

      // Step 3: Find if any vertex is above the top facet
      let anyVertexAbove = false;
      let minY = topFacetY;

      for (const vertex of rotatedPolygon) {
        if (vertex[1] < topFacetY - 1e-10) {
          // Vertex is above the top facet
          anyVertexAbove = true;
        }
        minY = Math.min(minY, vertex[1]);
      }

      // Step 4: Ensure top facet is at the top, and no other vertices are above it
      const offsetY = topFacetY - minY; // This moves the top facet to the highest point

      // Adjust the polygon so top facet is at the top
      polygonToDisplay = rotatedPolygon.map((vertex) => [
        vertex[0],
        // If any parts are above the top facet, we need to push the top facet higher
        // by offsetting all vertices
        anyVertexAbove ? vertex[1] - minY : vertex[1],
      ]);

      // Step 5: Invert the y-coordinates so top facet is at the top of the screen
      // (Canvas y-coord increases downward)
      polygonToDisplay = polygonToDisplay.map((vertex) => [
        vertex[0],
        -vertex[1],
      ]);
    }
    // Else - If not castable, we leave the polygon as is (no rotation needed)

    // Calculate bounds for the properly oriented polygon
    let minX = Infinity,
      minY = Infinity,
      maxX = -Infinity,
      maxY = -Infinity;
    for (const point of polygonToDisplay) {
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
    scaledPolygon = scalePolygon(polygonToDisplay, scale, offsetX, offsetY);

    // Draw the polygon
    drawPolygon(ctx, scaledPolygon);

    // Display result info
    const resultInfo = document.createElement("div");
    resultInfo.className = "result-info";

    if (result.castable && rotatedPolygon) {
      // For the rotation center, we need to apply the same transformations
      let center = result.center;

      // Step 1: Apply rotation using the same method as for the polygon
      const rotationMatrix = getRotationMatrix(result.topFacet);
      let rotatedCenter = applyTransformation(center, rotationMatrix);

      // Step 2: Apply the same y offset and inversion as the polygon
      // Find the y-coordinate of top facet after all transformations
      let topFacetY = null;
      for (let i = 0; i < rotatedPolygon.length; i++) {
        const v1 = rotatedPolygon[i];
        const v2 = rotatedPolygon[(i + 1) % rotatedPolygon.length];

        if (Math.abs(v1[1] - v2[1]) < 1e-10) {
          topFacetY = v1[1];
          break;
        }
      }

      if (topFacetY !== null) {
        // Find if any vertex is above the top facet
        let anyVertexAbove = false;
        let minY = topFacetY;

        for (const vertex of rotatedPolygon) {
          if (vertex[1] < topFacetY - 1e-10) {
            anyVertexAbove = true;
          }
          minY = Math.min(minY, vertex[1]);
        }

        // Apply the same transformations to the center as to the polygon
        if (anyVertexAbove) {
          rotatedCenter = [rotatedCenter[0], rotatedCenter[1] - minY];
        }
      }

      // Invert y-coordinate
      rotatedCenter = [rotatedCenter[0], -rotatedCenter[1]];

      // Scale and translate
      scaledCenter = [
        rotatedCenter[0] * scale + offsetX,
        rotatedCenter[1] * scale + offsetY,
      ];

      // Find the top-most horizontal edge after all transformations
      let topY = Infinity;
      for (let i = 0; i < scaledPolygon.length; i++) {
        const v1 = scaledPolygon[i];
        const v2 = scaledPolygon[(i + 1) % scaledPolygon.length];

        // If this is a horizontal edge (y coordinates are the same)
        if (Math.abs(v1[1] - v2[1]) < 1e-10) {
          if (v1[1] < topY) {
            topY = v1[1];
            scaledTopFacet = [v1, v2];
          }
        }
      }

      resultInfo.innerHTML = `<span class="castable">✓ Castable</span> - Rotation center: (${result.center[0].toFixed(
        2
      )}, ${result.center[1].toFixed(2)})`;

      // Draw the top facet
      if (scaledTopFacet) {
        drawEdge(ctx, scaledTopFacet, "#4CAF50");
      }

      // Draw the rotation center
      drawRotationCenter(ctx, scaledCenter);

      // Draw rotation paths
      drawRotationPath(ctx, scaledPolygon, scaledCenter);
    } else if (result.castable) {
      // Fallback if rotatedPolygon is not available but the shape is castable
      resultInfo.innerHTML = `<span class="castable">✓ Castable</span> - Rotation center: (${result.center[0].toFixed(
        2
      )}, ${result.center[1].toFixed(2)})`;

      // Just scale the original center and top facet directly
      scaledCenter = [
        result.center[0] * scale + offsetX,
        result.center[1] * scale + offsetY,
      ];

      scaledTopFacet = [
        [
          result.topFacet[0][0] * scale + offsetX,
          result.topFacet[0][1] * scale + offsetY,
        ],
        [
          result.topFacet[1][0] * scale + offsetX,
          result.topFacet[1][1] * scale + offsetY,
        ],
      ];

      if (scaledTopFacet) {
        drawEdge(ctx, scaledTopFacet, "#4CAF50");
      }

      drawRotationCenter(ctx, scaledCenter);
      drawRotationPath(ctx, scaledPolygon, scaledCenter);
    } else {
      resultInfo.innerHTML = `<span class="not-castable">✗ Not Castable</span> - ${result.message}`;
      // For non-castable shapes, we don't draw the top facet or rotation center
    }

    containerDiv.appendChild(resultInfo);

    // Add facets info
    const topFacets = getTopFacets(testCase.polygon);
    const facetsInfo = document.createElement("div");
    facetsInfo.className = "result-info";
    facetsInfo.textContent = `Found ${topFacets.length} potential top facet(s)`;
    containerDiv.appendChild(facetsInfo);

    // Add vertex count info
    const vertexInfo = document.createElement("div");
    vertexInfo.className = "result-info";
    vertexInfo.textContent = `Polygon has ${testCase.polygon.length} vertices`;
    containerDiv.appendChild(vertexInfo);
  } catch (error) {
    console.error(`Error visualizing ${testCase.name}:`, error.message);
    const errorInfo = document.createElement("div");
    errorInfo.className = "result-info not-castable";
    errorInfo.textContent = `Error: ${error.message}`;
    containerDiv.appendChild(errorInfo);
  }
}

// Helper function to get rotation matrix from top facet
function getRotationMatrix(topEdge) {
  const edgeVector = vectorSubtract(topEdge[1], topEdge[0]);
  const angle = -Math.atan2(edgeVector[1], edgeVector[0]);

  return {
    cos: Math.cos(angle),
    sin: Math.sin(angle),
    angle: angle,
  };
}

// Apply rotation matrix to a point
function applyTransformation(point, matrix) {
  return [
    point[0] * matrix.cos - point[1] * matrix.sin,
    point[0] * matrix.sin + point[1] * matrix.cos,
  ];
}

// ==========================================
// TEST CASES & DATA
// ==========================================

// Test cases - predefined shapes
const predefinedTestCases = [
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

// ==========================================
// INITIALIZATION AND EVENT HANDLING
// ==========================================

// Initialize the application when the DOM is loaded
document.addEventListener("DOMContentLoaded", function () {
  console.log("DOM loaded, initializing application...");

  // Set up event handler for the Generate button
  const generateButton = document.getElementById("generateButton");
  if (generateButton) {
    generateButton.addEventListener("click", handleGenerateButtonClick);
    console.log("Generate button handler attached");
  } else {
    console.error("Generate button not found in the DOM");
  }

  // Visualize predefined test cases
  predefinedTestCases.forEach((test) => visualizePolygon(test));
});

/**
 * Handle generate button click
 */
function handleGenerateButtonClick() {
  console.log("Generate button clicked!");

  // Add a loading message immediately for visual feedback
  const randomTestContainer = document.getElementById("randomTestContainer");
  randomTestContainer.innerHTML =
    '<div class="test-container"><h3>Generating random polygons...</h3></div>';

  // Use setTimeout to ensure the UI updates before we start the heavy computation
  setTimeout(function () {
    try {
      // Clear the loading message
      randomTestContainer.innerHTML = "";

      // Generate and visualize random test cases
      const randomTests = generateRandomTestCases(5);
      console.log("Generated random polygons:", randomTests);

      // Make sure random polygons are in clockwise order
      for (const test of randomTests) {
        test.polygon = ensureClockwise(test.polygon);
      }

      // Display random test case stats
      const statsDiv = document.createElement("div");
      statsDiv.className = "test-container";
      statsDiv.innerHTML = "<h3>Random Test Statistics</h3>";
      randomTestContainer.appendChild(statsDiv);

      let castableCount = 0;
      randomTests.forEach((test) => {
        try {
          const result = isPolygonCastable(test.polygon);
          if (result.castable) castableCount++;

          // Visualize each random polygon
          visualizePolygon(test, "randomTestContainer");
        } catch (error) {
          console.error(`Error testing ${test.name}:`, error);

          // Still display the polygon even if there was an error analyzing it
          const errorDiv = document.createElement("div");
          errorDiv.className = "test-container";
          errorDiv.innerHTML = `<h3>${test.name}</h3><div class="result-info not-castable">Error analyzing polygon: ${error.message}</div>`;
          randomTestContainer.appendChild(errorDiv);
        }
      });

      // Add statistics
      const resultInfo = document.createElement("div");
      resultInfo.className = "result-info";
      resultInfo.innerHTML = `<strong>Results:</strong> ${castableCount} out of ${
        randomTests.length
      } random polygons are castable (${(
        (castableCount / randomTests.length) *
        100
      ).toFixed(0)}%).`;
      statsDiv.appendChild(resultInfo);
    } catch (error) {
      console.error("Error in random polygon generation:", error);
      randomTestContainer.innerHTML = `<div class="test-container not-castable">Error generating random polygons: ${error.message}</div>`;
    }
  }, 10); // Short timeout to allow UI to update
}

// ==========================================
// EXPORT FUNCTIONS TO GLOBAL SCOPE
// ==========================================

// Make necessary functions available globally for use in inline scripts
window.isClockwise = isClockwise;
window.ensureClockwise = ensureClockwise;
window.generateRandomTestCases = generateRandomTestCases;
window.generateRandomPolygon = generateRandomPolygon;
window.generateRandomConvexPolygon = generateRandomConvexPolygon;
window.generateRandomNonConvexPolygon = generateRandomNonConvexPolygon;
window.isSimplePolygon = isSimplePolygon;
window.doSegmentsIntersect = doSegmentsIntersect;
window.visualizePolygon = visualizePolygon;
window.isPolygonCastable = isPolygonCastable;
window.getTopFacets = getTopFacets;
window.getRotationMatrix = getRotationMatrix;
window.applyTransformation = applyTransformation;
