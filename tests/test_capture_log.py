import logging
import warnings

import pytest

import swmmbuilder
from swmmbuilder import DataWarning, OptionWarning, capture_log

LOGGER = logging.getLogger("swmmbuilder")


def test_captures_swmmbuilder_records_at_level():
    with capture_log(level="INFO") as captured:
        logging.getLogger("swmmbuilder.gaps").info("Filled %d inverts", 3)
        logging.getLogger("swmmbuilder.network").debug("not captured")
        logging.getLogger("swmmbuilder.network").warning("Dropped 2 pipes")
    assert captured.text == (
        "INFO:swmmbuilder.gaps:Filled 3 inverts\nWARNING:swmmbuilder.network:Dropped 2 pipes\n"
    )


def test_level_as_number_and_unknown_name():
    with capture_log(level=logging.DEBUG) as captured:
        logging.getLogger("swmmbuilder.x").debug("detail")
    assert captured.text == "DEBUG:swmmbuilder.x:detail\n"
    with pytest.raises(ValueError, match="unknown level"):
        with capture_log(level="LOUD"):
            pass


def test_ignores_other_loggers():
    with capture_log() as captured:
        logging.getLogger("other").warning("elsewhere")
        logging.getLogger("swmmbuilderish").warning("elsewhere too")
        logging.getLogger().warning("root")
    assert captured.text == ""


def test_sequential_captures_are_independent():
    texts = []
    for i in range(2):
        with capture_log() as captured:
            logging.getLogger("swmmbuilder.swmm").info("build %d", i)
        texts.append(captured.text)
    assert texts == ["INFO:swmmbuilder.swmm:build 0\n", "INFO:swmmbuilder.swmm:build 1\n"]


def test_nested_captures_each_get_their_text():
    with capture_log() as outer:
        logging.getLogger("swmmbuilder.a").info("one")
        with capture_log(level="WARNING") as inner:
            logging.getLogger("swmmbuilder.a").info("two")
            logging.getLogger("swmmbuilder.a").warning("three")
    assert inner.text == "WARNING:swmmbuilder.a:three\n"
    assert outer.text == (
        "INFO:swmmbuilder.a:one\nINFO:swmmbuilder.a:two\nWARNING:swmmbuilder.a:three\n"
    )


def test_handler_removed_and_level_restored_after_exception():
    handlers = list(LOGGER.handlers)
    level = LOGGER.level
    showwarning = warnings.showwarning
    with pytest.raises(RuntimeError):
        with capture_log(level="DEBUG") as captured:
            assert LOGGER.getEffectiveLevel() == logging.DEBUG
            raise RuntimeError("build failed")
    assert LOGGER.handlers == handlers
    assert LOGGER.level == level
    assert warnings.showwarning is showwarning
    logging.getLogger("swmmbuilder.a").warning("after")
    assert captured.text == ""


def test_warnings_captured_and_still_emitted():
    with pytest.warns(UserWarning, match="not ours"):
        with pytest.warns(DataWarning, match="adverse slope"):
            with pytest.warns(OptionWarning):
                with capture_log() as captured:
                    warnings.warn("P7: adverse slope", DataWarning, stacklevel=1)
                    warnings.warn("keys name no node", OptionWarning, stacklevel=1)
                    warnings.warn("not ours", UserWarning, stacklevel=1)
    assert captured.text == (
        "WARNING:DataWarning:P7: adverse slope\nWARNING:OptionWarning:keys name no node\n"
    )


def test_warnings_below_level_not_captured():
    with pytest.warns(DataWarning):
        with capture_log(level="ERROR") as captured:
            warnings.warn("suspect", DataWarning, stacklevel=1)
    assert captured.text == ""


def test_path_written(tmp_path):
    path = tmp_path / "build_log.txt"
    with capture_log(path=path) as captured:
        logging.getLogger("swmmbuilder.a").info("café mètre")
    assert path.read_text(encoding="utf-8") == captured.text == "INFO:swmmbuilder.a:café mètre\n"


def test_exported():
    assert swmmbuilder.capture_log is capture_log
    assert "capture_log" in swmmbuilder.__all__
