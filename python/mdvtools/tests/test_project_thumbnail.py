import json
import os

from mdvtools import mdvproject as mp


def test_pick_thumbnail_follows_all_views_order_not_views_json_order():
    # views.json key order differs from the order the project lists its views
    views = {"second": {"viewImage": "img2"}, "first": {"viewImage": "img1"}}
    state = {"all_views": ["first", "second"]}
    assert mp.pick_thumbnail(state, views) == "img1"


def test_pick_thumbnail_skips_views_without_an_image():
    views = {"first": {"viewImage": ""}, "second": {}, "third": {"viewImage": "img3"}}
    state = {"all_views": ["first", "second", "third"]}
    assert mp.pick_thumbnail(state, views) == "img3"


def test_pick_thumbnail_falls_back_to_unlisted_views():
    views = {"listed": {}, "unlisted": {"viewImage": "img"}}
    assert mp.pick_thumbnail({"all_views": ["listed"]}, views) == "img"
    assert mp.pick_thumbnail({}, {}) is None


def _write(path, data, mtime):
    with open(path, "w") as f:
        json.dump(data, f)
    os.utime(path, (mtime, mtime))


def test_get_project_thumbnail_rereads_after_views_change(tmp_path):
    project_dir = str(tmp_path)
    statefile = os.path.join(project_dir, "state.json")
    viewsfile = os.path.join(project_dir, "views.json")
    _write(statefile, {"all_views": ["a", "b"]}, 1000)
    _write(viewsfile, {"b": {"viewImage": "img-b"}}, 1000)
    assert mp.get_project_thumbnail(project_dir) == "img-b"

    # the first view gets an image later
    _write(viewsfile, {"b": {"viewImage": "img-b"}, "a": {"viewImage": "img-a"}}, 2000)
    assert mp.get_project_thumbnail(project_dir) == "img-a"


def test_get_project_thumbnail_returns_none_for_missing_project(tmp_path):
    assert mp.get_project_thumbnail(str(tmp_path / "missing")) is None
