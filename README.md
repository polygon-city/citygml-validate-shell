# CityGML Validate Shell

Validate a CityGML shell against the [QIE suite](https://github.com/tudelft3d/CityGML-QIE-3Dvalidation/blob/master/errors/errors.md#shell)

## Usage

```javascript
var _ = require("lodash");
var citygmlPolygons = require("citygml-polygons");
var citygmlValidateShell = require("citygml-validate-shell");

var xml = "..."; // Some CityGML
var polygons = citygmlPolygons(xml);

// Validate shell as a whole
citygmlValidateShell(polygons, function(err, results) {
  _.each(results, function(vError) {
    // Should always be an error, but check anyway
    if (!vError || !vError[0]) {
      return;
    }

    // Output validation error name
    console.log(vError[0].message);
  });
});
```
