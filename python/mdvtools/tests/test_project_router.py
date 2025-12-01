"""
Tests for ProjectBlueprint routing behavior.

These tests verify that the custom routing implementation behaves consistently
with Flask's standard routing, particularly around:
- Proper escaping of special regex characters (dots, etc.)
- Flask/Werkzeug converter types (string, path, int, float)
- Edge cases that caused bugs in production
"""

import pytest
from mdvtools.project_router import ProjectBlueprint


class TestRoutePatternCompilation:
    """Test that Flask-style routes compile to correct regex patterns."""
    
    def test_dot_escaping_in_routes(self):
        """
        REGRESSION TEST: Dots in routes should be literal, not regex wildcards.
        
        Bug: Files named '.geojson' were matching '/<file>.json' route as '.ge' + '.json'
        because the dot wasn't being escaped in the regex pattern.
        
        Expected: Only files actually ending in '.json' should match.
        """
        bp = ProjectBlueprint("test", "", "/test/")
        
        @bp.route("/<file>.json")
        def json_route(file):
            return f"json: {file}"
        
        # Get the compiled pattern
        patterns = list(bp.routes.keys())
        assert len(patterns) == 1
        pattern = patterns[0]
        
        # Should match files ending in .json
        assert pattern.match("/state.json"), "Should match 'state.json'"
        assert pattern.match("/config.json"), "Should match 'config.json'"
        
        # Should NOT match .geojson (the literal dot must be escaped)
        match = pattern.match("/data.geojson")
        assert match is None, "Should NOT match 'data.geojson' (dot should be literal, not wildcard)"
        
        # Verify the dot is actually escaped in the pattern
        assert r'\.' in pattern.pattern, f"Pattern should contain escaped dot: {pattern.pattern}"
        print(f"✓ Pattern correctly escapes dot: {pattern.pattern}")
    
    def test_multiple_dots_in_routes(self):
        """Routes with multiple dots should escape all of them."""
        bp = ProjectBlueprint("test", "", "/test/")
        
        @bp.route("/<file>.tar.gz")
        def archive_route(file):
            return f"archive: {file}"
        
        pattern = list(bp.routes.keys())[0]
        
        # Should match exact pattern
        assert pattern.match("/backup.tar.gz")
        
        # Should NOT match variations
        assert not pattern.match("/backup_tar.gz")
        assert not pattern.match("/backupXtarYgz")
        
        # Both dots should be escaped
        assert pattern.pattern.count(r'\.') == 2


class TestConverterTypes:
    """Test different Flask/Werkzeug converter types."""
    
    def test_string_converter_excludes_slashes(self):
        """String converter (default) should NOT match paths with slashes."""
        bp = ProjectBlueprint("test", "", "/test/")
        
        @bp.route("/file/<name>")
        def string_route(name):
            return f"file: {name}"
        
        pattern = list(bp.routes.keys())[0]
        
        # Should match simple names
        assert pattern.match("/file/document")
        assert pattern.match("/file/image.png")
        
        # Should NOT match paths with slashes
        assert not pattern.match("/file/nested/path")
        assert not pattern.match("/file/dir/file.txt")
    
    def test_path_converter_includes_slashes(self):
        """Path converter should match everything including slashes."""
        bp = ProjectBlueprint("test", "", "/test/")
        
        @bp.route("/images/<path:filepath>")
        def path_route(filepath):
            return f"path: {filepath}"
        
        pattern = list(bp.routes.keys())[0]
        
        # Should match simple names
        assert pattern.match("/images/photo.jpg")
        
        # Should match nested paths
        assert pattern.match("/images/nested/dir/photo.jpg")
        assert pattern.match("/images/a/b/c/d/file.png")
    
    def test_int_converter_only_matches_integers(self):
        """Int converter should only match digit sequences."""
        bp = ProjectBlueprint("test", "", "/test/")
        
        @bp.route("/user/<int:id>")
        def int_route(id):
            return f"user: {id}"
        
        pattern = list(bp.routes.keys())[0]
        
        # Should match integers
        assert pattern.match("/user/123")
        assert pattern.match("/user/999999")
        
        # Should NOT match non-integers
        assert not pattern.match("/user/abc")
        assert not pattern.match("/user/12.5")
        assert not pattern.match("/user/12a")
    
    def test_float_converter_matches_numbers(self):
        """Float converter should match integers and decimals."""
        bp = ProjectBlueprint("test", "", "/test/")
        
        @bp.route("/price/<float:amount>")
        def float_route(amount):
            return f"price: {amount}"
        
        pattern = list(bp.routes.keys())[0]
        
        # Should match integers
        assert pattern.match("/price/100")
        
        # Should match floats
        assert pattern.match("/price/99.99")
        assert pattern.match("/price/0.5")
        
        # Should NOT match non-numeric
        assert not pattern.match("/price/abc")


class TestRegressionBugs:
    """Tests for specific bugs that were found in production."""
    
    def test_geojson_bug(self):
        """
        CRITICAL REGRESSION TEST
        
        Bug description:
        - Files named 'data.geojson' in images folder were incorrectly matching
          the '/<file>.json' route
        - The dot in '.json' was being treated as regex wildcard
        - Pattern matched 'data.ge' as the filename, then '.json' as any char + 'json'
        - Result: 404 because 'data.ge.json' doesn't exist
        
        This happened ONLY with ProjectBlueprint, not with Flask's native routing,
        making it harder to catch in single-project mode.
        """
        bp = ProjectBlueprint("test", "", "/test/")
        
        @bp.route("/<file>.json")
        def json_route(file):
            return f"json: {file}"
        
        @bp.route("/images/<path:path>")
        def images_route(path):
            return f"image: {path}"
        
        # The .geojson file should ONLY match the images route
        json_pattern = [p for p in bp.routes.keys() if 'images' not in p.pattern][0]
        images_pattern = [p for p in bp.routes.keys() if 'images' in p.pattern][0]
        
        # .geojson should NOT match the .json route
        assert not json_pattern.match("/data.geojson"), \
            "BUG: .geojson incorrectly matched .json route (dot not escaped)"
        
        # .geojson SHOULD match the images route
        assert images_pattern.match("/images/data.geojson"), \
            ".geojson files should match images/<path> route"
        
        # .json files should still work
        assert json_pattern.match("/state.json"), \
            "Regular .json files should still match"


class TestFlaskParity:
    """
    Tests to ensure ProjectBlueprint behaves like Flask's Blueprint.
    
    Note: These are behavioral tests. Full parity isn't possible due to the
    custom implementation, but core routing behavior should match.
    """
    
    def test_route_ordering(self):
        """Routes should be tried in registration order."""
        bp = ProjectBlueprint("test", "", "/test/")
        
        results = []
        
        @bp.route("/<path:anything>")
        def catch_all(anything):
            results.append("catch_all")
            return "catch all"
        
        @bp.route("/specific")
        def specific():
            results.append("specific")
            return "specific"
        
        # Both patterns would match "/specific", but the first registered wins
        # (This is different from Werkzeug which uses specificity rules)
        for pattern, (handler, _) in bp.routes.items():
            if pattern.match("/specific"):
                handler("specific")
                break
        
        # Due to dict ordering (Python 3.7+), catch_all was registered first
        assert results[0] == "catch_all"
    
    def test_empty_route(self):
        """Empty route should match exactly '/'."""
        bp = ProjectBlueprint("test", "", "/test/")
        
        @bp.route("/")
        def index():
            return "index"
        
        pattern = list(bp.routes.keys())[0]
        
        assert pattern.match("/")
        assert not pattern.match("/other")


class TestEdgeCases:
    """Test edge cases and unusual patterns."""
    
    def test_consecutive_converters(self):
        """Multiple converters in sequence should work."""
        bp = ProjectBlueprint("test", "", "/test/")
        
        @bp.route("/<category>/<item>")
        def two_parts(category, item):
            return f"{category}/{item}"
        
        pattern = list(bp.routes.keys())[0]
        
        assert pattern.match("/books/fiction")
        # Neither should match slashes (default string converter)
        assert not pattern.match("/books/fiction/novel")
    
    def test_converter_with_dots(self):
        """Converter followed by literal dot should work."""
        bp = ProjectBlueprint("test", "", "/test/")
        
        @bp.route("/<name>.txt")
        def text_file(name):
            return f"text: {name}"
        
        pattern = list(bp.routes.keys())[0]
        
        assert pattern.match("/readme.txt")
        assert not pattern.match("/readme_txt")  # Dot must be literal
    
    def test_special_chars_in_static_parts(self):
        """Special regex chars in static route parts should be escaped."""
        bp = ProjectBlueprint("test", "", "/test/")
        
        @bp.route("/api/v1.0/users")
        def api_v1():
            return "api"
        
        pattern = list(bp.routes.keys())[0]
        
        # Should match exactly
        assert pattern.match("/api/v1.0/users")
        
        # Should NOT match with dot as wildcard
        assert not pattern.match("/api/v1X0/users")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

