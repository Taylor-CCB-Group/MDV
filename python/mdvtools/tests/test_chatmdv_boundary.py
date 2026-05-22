from mdvtools.llm import check_chatmdv_boundary as boundary


def test_is_allowed_llm_prefix():
    assert boundary.is_allowed("python/mdvtools/llm/templates.py")


def test_is_allowed_explicit_test_files():
    for path in boundary.ALLOWED_FILES:
        assert boundary.is_allowed(path)


def test_is_allowed_test_chat_glob():
    assert boundary.is_allowed("python/mdvtools/tests/test_chat_cli_module.py")
    assert boundary.is_allowed("python/mdvtools/tests/test_chatmdv_roles_injection.py")
    assert boundary.is_allowed("python/mdvtools/tests/test_chatlog_handler.py")


def test_is_allowed_policy_and_preflight_globs():
    assert boundary.is_allowed(
        "python/mdvtools/tests/test_chat_first_text_table_policy.py"
    )
    assert boundary.is_allowed("python/mdvtools/tests/test_code_preflight.py")


def test_is_allowed_rejects_unrelated_tests():
    assert not boundary.is_allowed("python/mdvtools/tests/test_langchain.py")
    assert not boundary.is_allowed("python/mdvtools/mdvproject.py")
