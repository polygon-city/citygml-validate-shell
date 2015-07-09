// Note: QIE paper states that only LOD1–3 geometries are to be tested,
// therefore any cavities (interior shells) in solids can be ignored. This
// implies that a solid has exactly 1 shell representing its exterior surfaces.
//
// See: https://github.com/tudelft3d/CityGML-QIE-3Dvalidation/blob/master/errors/errors.md#shell

// TODO: Test against a CityGML dataset containing valid and invalid geometry
// TODO: Support inner polygons

var _ = require("lodash");
var async = require("async");
var sylvester = require("sylvester");
var polygonjs = require("polygon");

var citygmlBoundaries = require("citygml-boundaries");
var citygmlPoints = require("citygml-points");
var points3dto2d = require("points-3d-to-2d");

// Input: [polygon, polygon, ...]
var citygmlValidateShell = function(shell, callback) {
  // Validate shell
  // Validation errors are stored within results array
  async.series([
    GE_S_TOO_FEW_POLYGONS(shell),
    GE_S_NOT_CLOSED(shell),
    GE_S_NON_MANIFOLD_VERTEX(shell),
    GE_S_NON_MANIFOLD_EDGE(shell),
    GE_S_MULTIPLE_CONNECTED_COMPONENTS(shell),
    GE_S_SELF_INTERSECTION(shell),
    GE_S_POLYGON_WRONG_ORIENTATION(shell)
    // GE_S_ALL_POLYGONS_WRONG_ORIENTATION(shell)
  ], function(err, results) {
    // Remove passes (undefined)
    var failures = _.filter(results, function(result) {
      return (result && result[0]);
    });

    callback(err, failures);
  });
};

var GE_S_TOO_FEW_POLYGONS = function(shell) {
  // A shell should have at least 4 polygons - the simplest volumetric shape in
  // 3D is a tetrahedron.

  // Async pattern as per Pro JavaScript Frameworks book
  // Missing process.nextTick trick inside the function
  return function(callback) {
    if (shell.length > 3) {
      callback(null);
    } else {
      callback(null, [new Error("GE_S_TOO_FEW_POLYGONS: A shell should have at least 4 polygons"), shell.length]);
    }
  };
};

// IMPL: Look for edges that are only referenced by 1 polygon
// TODO: Store reference to edges
// TODO: Count number of times each edge is referenced
// TODO: Fail if edge count is 1
var GE_S_NOT_CLOSED = function(shell) {
  // The shell must not have 'holes', ie it must be 'watertight'. This refers
  // only to the topology of the shell, not to its geometry (see
  // GE_S_SELF_INTERSECTION).
  //
  // Example: https://github.com/tudelft3d/CityGML-QIE-3Dvalidation/raw/master/errors/figs/302.png

  return function(callback) {
    var polygonPoints = [];

    _.each(shell, function(polygonXML) {
      var boundaries = citygmlBoundaries(polygonXML);

      _.each(boundaries.exterior.concat(boundaries.interior), function(boundary) {
        var bPoints = citygmlPoints(boundary);
        polygonPoints.push(bPoints);
      });
    });

    var edges = findEdges(polygonPoints);
    var edgeCounts = countMatchingEdges(edges);

    var holes = [];

    _.each(edgeCounts, function(count, edgeId) {
      if (count < 2) {
        holes.push(edges[edgeId]);
      }
    });

    if (holes.length === 0) {
      callback(null);
    } else {
      callback(null, [new Error("GE_S_NOT_CLOSED: The shell must be watertight (ie. no holes)"), holes]);
    }
  };
};

// NOTE: Most non-manifold errors will fail GE_S_NOT_CLOSED before this test
// IMPL: Look for vertices that are referenced by 2 or more faces but don’t share an edge
// IMPL: Look for polygons that share vertices (all of them, otherwise fail GE_S_MULTIPLE_CONNECTED_COMPONENTS?)
// IMPL: Look for polygons that share edges
// IMPL: Fail if polygons that share a vertex don't also share an edge
var GE_S_NON_MANIFOLD_VERTEX = function(shell) {
  // Polygons shared by a vertex must also share an edge
  //
  // Example: https://github.com/tudelft3d/CityGML-QIE-3Dvalidation/raw/master/errors/figs/303.png

  return function(callback) {
    var polygonPoints = [];

    _.each(shell, function(polygonXML) {
      var boundaries = citygmlBoundaries(polygonXML);

      _.each(boundaries.exterior.concat(boundaries.interior), function(boundary) {
        var bPoints = citygmlPoints(boundary);
        polygonPoints.push(bPoints);
      });
    });

    var edges = findEdges(polygonPoints);

    var edgesByPolygon = findEdgesByPolygon(edges);
    var polygonsByVertex = {};

      _.each(polygonPoints, function(checkRingPoints, polygonIndex) {
        _.each(checkRingPoints, function(point) {
          var pointId = point.toString();

          if (!polygonsByVertex[pointId]) {
            polygonsByVertex[pointId] = [];
          }

          // Skip polygon if added already
          if (_.indexOf(polygonsByVertex[pointId], polygonIndex) < 0) {
            polygonsByVertex[pointId].push(polygonIndex);
          }
        });
      });
    // });

    var badVertices = [];

    _.each(polygonsByVertex, function(polygonIndexes, vertexId) {
      // Skip if only one polygon is using this vertex
      // This shouldn't fail here but will fail GE_S_NOT_CLOSED
      if (polygonIndexes.length < 2) {
        return;
      }

      var indexesCopy = _.clone(polygonIndexes);

      var polygonIndex;
      var polygonEdges;
      var checkEdges;
      var edgeMatch;

      var polygonEdgeIdSplit;
      var checkEdgeIdSplit;


      while (indexesCopy.length > 0) {
        polygonIndex = indexesCopy.shift();
        polygonEdges = edgesByPolygon[polygonIndex];

        // Find an edge shared by the current polygon and one of the others
        // sharing this vertex
        _.each(indexesCopy, function(checkIndex) {
          checkEdges = edgesByPolygon[checkIndex];

          edgeMatch = _.find(polygonEdges, function(polygonEdge) {
            polygonEdgeIdSplit = polygonEdge.split(":");

            return _.find(checkEdges, function(checkEdge) {
              checkEdgeIdSplit = checkEdge.split(":");

              return (polygonEdgeIdSplit[0] === checkEdgeIdSplit[1] && polygonEdgeIdSplit[1] === checkEdgeIdSplit[0]);
            });
          });

          // Exit the each loop early
          if (edgeMatch) {
            return false;
          }
        });

        // Exit the while loop early
        if (edgeMatch) {
          break;
        }
      }

      // Check for lack of required edges
      if (!edgeMatch) {
        badVertices.push(vertexId);
      }
    });

    if (badVertices.length === 0) {
      callback(null)
    } else {
      callback(null, [new Error("GE_S_NON_MANIFOLD_VERTEX: Polygons shared by a vertex must also share an edge"), badVertices]);
    }
  };
};

// IMPL: Look for edges that are referenced more than twice
var GE_S_NON_MANIFOLD_EDGE = function(shell) {
  // Each edge of a shell should have exactly 2 incident polygons.
  //
  // Example: https://github.com/tudelft3d/CityGML-QIE-3Dvalidation/raw/master/errors/figs/304.png

  return function(callback) {
    var polygonPoints = [];

    _.each(shell, function(polygonXML) {
      var boundaries = citygmlBoundaries(polygonXML);

      _.each(boundaries.exterior.concat(boundaries.interior), function(boundary) {
        var bPoints = citygmlPoints(boundary);
        polygonPoints.push(bPoints);
      });
    });

    var edges = findEdges(polygonPoints);
    var edgeCounts = countMatchingEdges(edges);

    var nonManifolds = [];

    _.each(edgeCounts, function(count, edgeId) {
      if (count !== 2) {
        nonManifolds.push(edges[edgeId]);
      }
    });

    if (nonManifolds.length === 0) {
      callback(null);
    } else {
      callback(null, [new Error("GE_S_NON_MANIFOLD_EDGE: Each edge of a shell should have exactly 2 incident polygons"), nonManifolds]);
    }
  };
};

// IMPL: http://stackoverflow.com/a/21901612/997339
var GE_S_MULTIPLE_CONNECTED_COMPONENTS = function(shell) {
  // Polygons that are not connected to the shell should be reported as an
  // error.
  //
  // Example: https://github.com/tudelft3d/CityGML-QIE-3Dvalidation/raw/master/errors/figs/305.png

  return function(callback) {
    var polygonPoints = [];

    _.each(shell, function(polygonXML) {
      var boundaries = citygmlBoundaries(polygonXML);

      _.each(boundaries.exterior.concat(boundaries.interior), function(boundary) {
        var bPoints = citygmlPoints(boundary);
        polygonPoints.push(bPoints);
      });
    });

    var edges = findEdges(polygonPoints);

    // List of edges as strings: [["0,0,0", "0,1,0"], [...]]
    var edgeList = _.map(edges, function(edge) {
      return [edge[0].toString(), edge[1].toString()];
    });

    var edgeListAdj = convert_edgelist_to_adjlist(edgeList);

    var groups = [];
    var visited = {};
    var v;

    for (v in edgeListAdj) {
      if (edgeListAdj.hasOwnProperty(v) && !visited[v]) {
        groups.push(bfs(v, edgeListAdj, visited));
      }
    }

    // Fail if groups is more than 1
    if (groups.length === 1) {
      callback(null);
    } else {
      callback(null, [new Error("GE_S_MULTIPLE_CONNECTED_COMPONENTS: All polygons must be connected to the shell"), groups]);
    }
  };
};

// Based on: https://github.com/mrdoob/three.js/blob/master/src/math/Plane.js
// TODO: Support polygon holes
// TODO: Support i_GE_S_SELF_INTERSECTION_3.gml failure
// TODO: Fix NaN intersection results in i_GE_S_SELF_INTERSECTION_4.gml failure
var GE_S_SELF_INTERSECTION = function(shell) {
  // If topology of the shell is correct and the shell is closed (thus no error
  // GE_S_TOO_FEW_POLYGONS / GE_S_NOT_CLOSED / GE_S_NON_MANIFOLD_VERTEX /
  // GE_S_NON_MANIFOLD_EDGE / GE_S_MULTIPLE_CONNECTED_COMPONENTS), it is
  // possible that the geometry introduces errors, eg intersections. For
  // instance, the topology of both these shells is identical, but the geometry
  // differs. The left shell is valid while the right one is invalid.
  //
  // Example: https://github.com/tudelft3d/CityGML-QIE-3Dvalidation/raw/master/errors/figs/306.png

  return function(callback) {
    var polygonPoints = [];

    _.each(shell, function(polygonXML) {
      var boundaries = citygmlBoundaries(polygonXML);

      _.each(boundaries.exterior.concat(boundaries.interior), function(boundary) {
        var bPoints = citygmlPoints(boundary);
        polygonPoints.push(bPoints);
      });
    });

    // var checkPolygons = _.clone(shell);
    // var edges = findEdges(shell);
    var edges = findEdges(polygonPoints);

    var intersections = [];

    _.each(polygonPoints, function(checkRingPoints, polygonIndex) {
      // var checkBoundaries = citygmlBoundaries(polygonXML);
      // var checkRingPoints = citygmlPoints(checkBoundaries.exterior[0]);

      // TODO: Should be abstracted into a Plane module as it's also used in
      // GE_P_NON_PLANAR_POLYGON_DISTANCE_PLANE
      var plane = {
        normal: undefined,
        constant: undefined
      };

      // Describe plane
      plane.normal = $V(normalUnit(checkRingPoints[0], checkRingPoints[1], checkRingPoints[2]));
      plane.constant = - $V(checkRingPoints[0]).dot(plane.normal);

      _.each(edges, function(edge) {
        // Ignore edges of the same polygon
        if (_.contains(edge[2], polygonIndex)) {
          return;
        }

        // Check each edge for intersection against the plane
        var intersection = lineIntersectsPlane(plane, edge);

        if (intersection) {
          intersections.push(intersection);
        }
      });
    });

    if (intersections.length === 0) {
      callback(null);
    } else {
      callback(null, [new Error("GE_S_SELF_INTERSECTION: Shell must not intersect itself"), intersections]);
    }
  };
};

// Note: This also tests for GE_S_ALL_POLYGONS_WRONG_ORIENTATION
//
// Note: See note at top about cavities being ignored in QIE testing. In short,
// you can assume that the given shell is outside and so all exterior boundaries
// should be anti-clockwise.
//
// Note: Wising depends entirely on which side you want to be viewing the
// polygon from.
//
// Note: This overlooks issues with collinear points (which should fail before
// anyway)
//
// IMPL: Find a polygon we know is facing up/down based on horizontal plane (eg.
// top-most must have normal angled above horizontal, bottom-most below). Seems
// impossible to have a valid shell with a "top" face that points at or below
// horizontal.
// IMPL: Look at all neighbouring polygons and check that orientation is opposite
var GE_S_POLYGON_WRONG_ORIENTATION = function(shell) {
  // If a polygon is used to construct a shell, its exterior ring must be
  // oriented in such as way that when viewed from outside the shell the points
  // are ordered counterclockwise.
  //
  // And GE_S_ALL_POLYGONS_WRONG_ORIENTATION, if all the polygons have the wrong
  // orientation (as defined in GE_S_POLYGON_WRONG_ORIENTATION), ie they all
  // point inwards.

  // TODO: This approach is way too complex – find a simpler way
  // TODO: Run more checks to ensure this approach is valid in all circumstances
  // - Overhanging top face /|
  // - Pyramid /\
  // - Wedge on top \ /
  return function(callback) {
    var polygonIndexes = [];
    var polygonPoints = [];
    var polygonWindings = [];

    _.each(shell, function(polygonXML, polygonIndex) {
      var boundaries = citygmlBoundaries(polygonXML);
      var exteriorPoints = citygmlPoints(boundaries.exterior[0]);

      polygonIndexes.push(polygonIndex);
      polygonPoints.push(exteriorPoints);
    });

    if (polygonPoints.length === 0) {
      callback(new Error("GE_S_POLYGON_WRONG_ORIENTATION: No polygons found"));
      return;
    }

    // Polygons that share the top-most Z position
    var topPolygons = {};
    var maxZ;

    _.each(polygonPoints, function(pPoints) {
      var points2d = points3dto2d(pPoints);
      var poly = polygonjs(points2d.points);

      // Find winding
      polygonWindings.push(poly.winding());

      // Find top-most point
      _.each(pPoints, function(point) {
        if (maxZ === undefined || point[2] > maxZ) {
          maxZ = point[2];
        }
      });
    });

    _.each(polygonPoints, function(pPoints, polygonIndex) {
      var top;

      var totalZ = 0;
      var avgZ;

      // Find polygons that use top-most point
      _.each(pPoints, function(point) {
        if (point[2] === maxZ) {
          top = true;
        }

        totalZ += point[2];
      });

      if (top) {
        avgZ = totalZ / pPoints.length;
        topPolygons[polygonIndex] = [pPoints, avgZ];
      }
    });

    var maxNormalZ;
    var outerPolygon;
    var outerPolygonIndex;
    var outerPolygonNormal;

    var flipped = [];
    var outerFlipped = false;

    // Pick outer polygon with the largest absolute Z normal as this will always     // be on top after the previous filters in a valid solid (has the shallowest
    // angle) and will filter out overhanging polygons
    _.each(topPolygons, function(oPolygon, oPolygonIndex) {
      var p0, p1, p2;
      var normal;
      var collinearThreshold = 0.01;

      // Find first sequence of points that aren't collinear
      _.each(oPolygon[0], function(point, pIndex) {
        // Exit if no more points are available
        if (pIndex === oPolygon[0].length - 2) {
          return false;
        }

        p0 = $V(point);
        p1 = $V(oPolygon[0][pIndex+1]);
        p2 = $V(oPolygon[0][pIndex+2]);

        // Colinear or near-colinear?
        var cross = p0.subtract(p1).cross(p0.subtract(p2));

        // Exit if non-collinear points are found
        if (Math.abs(cross.e(1)) > collinearThreshold || Math.abs(cross.e(2)) > collinearThreshold || Math.abs(cross.e(3)) > collinearThreshold) {
          normal = normalUnit(p0, p1, p2);
          return false;
        }
      });

      // Skip first point to avoid normal inconsistancies with polygons that
      // have identical start and end points.
      // var normal = normalUnit(oPolygon[0][0], oPolygon[0][1], oPolygon[0][2]);
      var absNormalZ = Math.abs(normal.e(3));

      if (maxNormalZ === undefined || absNormalZ >= maxNormalZ) {
        // It's possible for multiple polygons to share both the same top vertex
        // and have the same absolute Z normal. This is likely because one of
        // the polygons is invalid (eg. not planar) and so the polygon with the
        // most vertices shared with the top vertex is chosen.
        if (absNormalZ >= maxNormalZ) {
          // Compare average Z values to pick top-most polygon
          if (topPolygons[oPolygonIndex][1] < topPolygons[outerPolygonIndex]) {
            return;
          }
        }

        maxNormalZ = absNormalZ;
        outerPolygon = oPolygon;
        outerPolygonIndex = Number(oPolygonIndex);
        outerPolygonNormal = normal;
      }
    });

    if (!outerPolygon) {
      callback(new Error("GE_S_POLYGON_WRONG_ORIENTATION: Outer polygon not found"));
      return;
    }

    var checkPoints2d = points3dto2d(outerPolygon[0], false);
    var checkPolygon = polygonjs(checkPoints2d.points);
    var checkWinding = checkPolygon.winding();

    // If normal has a negative Z then fail as shallowest (top-most) polygon
    // will always have a positive Z normal
    if (outerPolygonNormal.e(3) < 0) {
      var insideOut = true;

      var indexes = [];

      // Fail GE_S_ALL_POLYGONS_WRONG_ORIENTATION if all polygons have the same,
      // incorrect winding
      _.each(polygonPoints, function(ePoints, eIndex) {
        var ePoints2d = points3dto2d(ePoints, false);
        var ePolygon = polygonjs(ePoints2d.points);
        var eWinding = ePolygon.winding();

        indexes.push(eIndex);

        if (eWinding != checkWinding) {
          insideOut = false;
          return false;
        }
      });

      if (insideOut) {
        callback(null, [new Error("GE_S_ALL_POLYGONS_WRONG_ORIENTATION: All the polygons have the wrong orientation"), indexes]);
        return;
      } else {
        // This will cause only flipped polygons to be valid, as they match up
        // with the flipped outer polygon.
        outerFlipped = true;
      }
    }

    // This polygon is pointing outward and can be used to walk the rest of
    // the solid to check and fix normals
    //
    // All other polygons should have the same winding as this one, assuming
    // the solid is manifold

    // Find edges
    var edges = findEdges(polygonPoints);
    var edgesByPolygon = findEdgesByPolygon(edges);
    var edgeCounts = countMatchingEdges(edges);

    // If an edge count is anything other than 2 then there has to be a polygon
    // that has opposite winding to those adjacent to it (ie. flipped normal)
    //
    // TODO: This isn't robust enough – it fails when polygons aren't properly
    // connected, which aren't valid either but aren't part of this test

    _.each(edgeCounts, function(count, edgeId) {
      if (count !== 2) {
        // At this point we know edges[edgeId][2] contains the indexes for all
        // the polygons which have incorrectly-matched edges
        _.each(edges[edgeId][2], function(pIndex) {
          var pEdges = edgesByPolygon[pIndex];

          var pFlipped = true;

          // If all edges have a count of 1 then this polygon is flipped,
          // otherwise it's valid as it shares edges with other polygons
          _.each(pEdges, function(pEdge) {
            if (edgeCounts[pEdge] === 2) {
              pFlipped = false;
            }
          });

          if (pFlipped) {
            if (!_.contains(flipped, pIndex)) {
              flipped.push(pIndex);
            }
          }

          // REMOVED: Winding checks can't be relied on – winding can be either
          // the same or opposite depending on how 2D projection happens
          //
          // var polyWinding = polygonWindings[pIndex];
          //
          // // Compare windings (different === failure)
          // if (polyWinding !== checkWinding) {
          //   if (!_.contains(flipped, pIndex)) {
          //     // Add polygon index to failures
          //     flipped.push(pIndex);
          //   }
          // }
        });
      }
    });

    // All edges match up, now check windings match
    // It's possible that a polygon can have matched edges yet still be pointing
    // the wrong way. Eg. A flipped polygon would behave this way.
    // We can't just check winding *or* edges as it's also possible for all
    // polygons to be wound correctly yet not share edges properly.
    if (flipped.length === 0) {
      _.each(polygonWindings, function(winding, index) {
        if (winding !== checkWinding) {
          if (!_.contains(flipped, index)) {
            flipped.push(index);
          }
        }
      });
    }

    // If outer polygon is flipped then everything in the flipped array is
    // actually valid (as it was compared against a flipped polygon). We can
    // diff the flipped array with the original indexes to get the actually
    // flipped polygons.
    if (outerFlipped) {
      flipped = _.difference(polygonIndexes, flipped);
    }

    callback(null, [new Error("GE_S_POLYGON_WRONG_ORIENTATION: When an exterior polygon is viewed from outside the shell the points must be ordered counterclockwise"), flipped]);
  };
};

var findEdges = function(polygons) {
  var edges = {};

  // REMOVED: CityGML-specific processing should be done externally
  // _.each(shell, function(polygonXML, polygonIndex) {
  //   var boundaries = citygmlBoundaries(polygonXML);
  //
  //   var checkRings = boundaries.exterior.concat(boundaries.interior);
  //   var checkRingXML;
  //   var checkRingPoints;
  //   var prevPoint;
  //
  //   _.each(checkRings, function(checkRingXML) {
  //     checkRingPoints = citygmlPoints(checkRingXML);

  _.each(polygons, function(polygon, polygonIndex) {
    var prevPoint;
    _.each(polygon, function(point) {
      if (!prevPoint) {
        prevPoint = point;
        return;
      }

      // Serialise edge
      var edgeId = prevPoint.toString() + ":" + point.toString();

      // First time create new edge definition, second time only add polygon
      // index to edge
      if (!edges[edgeId]) {
        edges[edgeId] = [prevPoint, point, [polygonIndex]];
      } else {
        edges[edgeId][2].push(polygonIndex);
      }

      prevPoint = point;
    });
  });

  //   });
  // });

  return edges;
};

var countMatchingEdges = function(edges) {
  var edgeCounts = {};

  _.each(edges, function(edge, edgeId) {
    var idParts = edgeId.split(":");

    if (edgeCounts[edgeId]) {
      edgeCounts[edgeId] += 1;
    } else {
      edgeCounts[edgeId] = 1;
    }

    // Look for reversed edge (ie. a shared one)
    var sharedEdgeId = idParts[1] + ":" + idParts[0];
    var sharedEdge = edges[sharedEdgeId];

    if (sharedEdge) {
      if (edgeCounts[sharedEdgeId]) {
        edgeCounts[sharedEdgeId] += 1;
      } else {
        edgeCounts[sharedEdgeId] = 1;
      }
    }
  });

  return edgeCounts;
};

var findEdgesByPolygon = function(edges) {
  var edgesByPolygon = {};

  _.each(edges, function(edge, edgeId) {
    _.each(edge[2], function(pIndex) {
      if (!edgesByPolygon[pIndex]) {
        edgesByPolygon[pIndex] = [];
      }

      edgesByPolygon[pIndex].push(edgeId);
    });
  });

  return edgesByPolygon;
};

// From: http://stackoverflow.com/a/21901612/997339
// Converts an edgelist to an adjacency list representation
// In this program, we use a dictionary as an adjacency list,
// where each key is a vertex, and each value is a list of all
// vertices adjacent to that vertex
var convert_edgelist_to_adjlist = function(edgelist) {
  var adjlist = {};
  var i, len, pair, u, v;
  for (i = 0, len = edgelist.length; i < len; i += 1) {
    pair = edgelist[i];
    u = pair[0];
    v = pair[1];
    if (adjlist[u]) {
      // append vertex v to edgelist of vertex u
      adjlist[u].push(v);
    } else {
      // vertex u is not in adjlist, create new adjacency list for it
      adjlist[u] = [v];
    }
    if (adjlist[v]) {
      adjlist[v].push(u);
    } else {
      adjlist[v] = [u];
    }
  }
  return adjlist;
};

// From: http://stackoverflow.com/a/21901612/997339
// Breadth First Search using adjacency list
var bfs = function(v, adjlist, visited) {
  var q = [];
  var current_group = [];
  var i, len, adjV, nextVertex;
  q.push(v);
  visited[v] = true;
  while (q.length > 0) {
    v = q.shift();
    current_group.push(v);
    // Go through adjacency list of vertex v, and push any unvisited
    // vertex onto the queue.
    // This is more efficient than our earlier approach of going
    // through an edge list.
    adjV = adjlist[v];
    for (i = 0, len = adjV.length; i < len; i += 1) {
      nextVertex = adjV[i];
      if (!visited[nextVertex]) {
        q.push(nextVertex);
        visited[nextVertex] = true;
      }
    }
  }
  return current_group;
};

var normalUnit = function(p1, p2, p3) {
  var v1 = $V(p1);
  var v2 = $V(p2);
  var v3 = $V(p3);

  var a = v2.subtract(v1);
  var b = v3.subtract(v1);

  return a.cross(b).toUnitVector();
};

// From: https://github.com/mrdoob/three.js/blob/master/src/math/Plane.js
// TODO: Fix NaN intersection results (eg. in i_GE_S_SELF_INTERSECTION_4.gml)
var lineIntersectsPlane = function(plane, line) {
  var startVector = $V(line[0]);
  var endVector = $V(line[1]);

  // End vector minus start vector
  var direction = endVector.subtract(startVector);

  var denominator = plane.normal.dot(direction);

  if (denominator === 0) {
    // Note: Removed as unsure we need this value returned
    // // Line is coplanar, return origin
    // if (this.distanceToPoint(line.start) === 0) {
    //   return result.copy(line.start);
    // }

    // Unsure if this is the correct method to handle this case.
    return undefined;
  }

  var t = - (startVector.dot(plane.normal) + plane.constant) / denominator;

  // TODO: Just less / more than, or equal to as well?
  // Leaving as just less / more than seems to fail even for valid shells (eg. a
  // simple box retuns t values of 0 and 1)
  if (t <= 0 || t >= 1) {
    return undefined;
  }

  var result = direction.dup().multiply(t).add(startVector);

  return result;
};

module.exports = citygmlValidateShell;
