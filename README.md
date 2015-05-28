# CityGML Validate Shell

Validate a CityGML shell against the [QIE suite](https://github.com/tudelft3d/CityGML-QIE-3Dvalidation/blob/master/errors/errors.md#shell)

## Usage

```javascript
var citygmlPolygons = require("citygml-polygons");
var citygmlValidateShell = require("citygml-validate-shell");

var xml = "..."; // Some CityGML
var polygons = citygmlPolygons(xml);

// Validate shell as a whole
citygmlValidateShell(polygons, function(err, results) {
  if (err) {
    console.log("Shell not valid:", err, results);
  } else {
    console.log("Shell valid");
  }
});
```
