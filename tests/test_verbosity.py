import chem.verbosity as v


def test_is_quiet_unset(monkeypatch):
    monkeypatch.delenv("CHEM_QUIETNESS", raising=False)
    assert v.is_quiet() is False


def test_is_quiet_falsy_values(monkeypatch):
    for val in ["0", "N", "n", "FALSE", "false", "False"]:
        monkeypatch.setenv("CHEM_QUIETNESS", val)
        assert v.is_quiet() is False, val


def test_is_quiet_truthy_values(monkeypatch):
    for val in ["1", "Y", "y", "TRUE", "true", "yes", "anything"]:
        monkeypatch.setenv("CHEM_QUIETNESS", val)
        assert v.is_quiet() is True, val


def test_logged_prints_call_by_default(monkeypatch, capsys):
    monkeypatch.delenv("CHEM_QUIETNESS", raising=False)

    @v.logged
    def add(a, b=1):
        return a + b

    assert add(2, b=3) == 5
    err = capsys.readouterr().err
    assert "add(a=2, b=3)" in err


def test_logged_silent_when_quiet(monkeypatch, capsys):
    monkeypatch.setenv("CHEM_QUIETNESS", "1")

    @v.logged
    def add(a, b=1):
        return a + b

    assert add(2, b=3) == 5
    err = capsys.readouterr().err
    assert err == ""
