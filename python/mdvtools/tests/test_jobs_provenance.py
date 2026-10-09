from mdvtools.jobs.provenance import content_hash


def test_content_hash_is_stable_and_param_sensitive():
    a = content_hash("concat_columns", {"x": 1, "y": 2}, None)
    assert a == content_hash(
        "concat_columns", {"y": 2, "x": 1}, None
    )  # key order independent
    assert a != content_hash(
        "concat_columns", {"x": 1, "y": 3}, None
    )  # params change → new hash
    assert a != content_hash(
        "concat_columns", {"x": 1, "y": 2}, "sub"
    )  # subset change → new hash
