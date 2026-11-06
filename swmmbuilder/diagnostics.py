"""How swmmbuilder reports problems.

Each problem goes to one of four channels, chosen by what the user can do about it:

``raise`` (``ValueError``)
    The output would be wrong, or SWMM would refuse it. Every such problem found in one check
    is listed in one error (see :class:`Problems`).

``warnings.warn`` (a :class:`SwmmBuilderWarning`)
    The result the user asked for is written, but suspect data or an option the user should act
    on affects it. Filterable by category, and testable with ``pytest.warns``:

    - :class:`DataWarning`: the data are suspect (an adverse slope, a pipe below its node, a rim
      not above the invert, an unknown shape code);
    - :class:`DesignWarning`: a design check failed (e.g. ``force_main_checks``) or a value
      breaks a design convention (a pattern whose mean is not 1);
    - :class:`OptionWarning`: an option was partly ignored or needs confirming (force mains not
      pumped, pump-station nodes not in ``pump_stations``, ``height_correction`` keys that name
      no node, plot values whose IDs are not in the network).

``logger.warning``
    Input that was dropped, repaired or ignored where the caller cannot act item by item, or
    which does not change the hydraulics: clean-up of whole-city data on loading, nodes left out
    of the map coordinates, cosmetic plot issues.

``logger.info``
    An audit trail: values filled or derived by rules, pipes removed by the filters, what was
    read and written.

:func:`capture_log` collects the log and the warnings of one build as text.
"""

import difflib
import io
import logging
import warnings
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TextIO


class SwmmBuilderWarning(UserWarning):
    """Base class of the warnings of swmmbuilder."""


class DataWarning(SwmmBuilderWarning):
    """The data are suspect; the model is written as the data say."""


class DesignWarning(SwmmBuilderWarning):
    """A design check failed, or a value breaks a design convention."""


class OptionWarning(SwmmBuilderWarning):
    """An option was partly ignored, or needs confirming."""


class Problems:
    """Problems collected by one check, raised together as one ``ValueError``::

    problems = Problems("Invalid pumps")
    problems.add("P1: curve must be ...")
    problems.raise_if_any()  # ValueError: Invalid pumps:\\n  P1: curve must be ...
    """

    def __init__(self, title, items=()):
        self.title = title
        self.items = list(items)

    def add(self, problem):
        self.items.append(problem)

    def extend(self, problems):
        self.items.extend(problems)

    def __bool__(self):
        return bool(self.items)

    def __len__(self):
        return len(self.items)

    def message(self):
        """The title and one indented line per problem."""
        return f"{self.title}:\n  " + "\n  ".join(self.items)

    def raise_if_any(self):
        """Raise ValueError with :meth:`message` if any problem was added."""
        if self.items:
            raise ValueError(self.message())


def check_options(given, known, where):
    """Raise TypeError if a keyword of ``given`` is not in ``known``, suggesting the closest
    known names (e.g. ``size_overide`` -> ``size_overrides``), instead of the TypeError of a
    function the keywords are passed on to.

    Parameters
    ----------
    given : iterable of str
        The keyword names passed.
    known : iterable of str
        The keyword names accepted.
    where : str
        What the keywords were passed to, for the message, e.g. ``'build_model'``.
    """
    known = sorted(set(known))
    unknown = sorted(set(given) - set(known))
    if not unknown:
        return
    hints = []
    for name in unknown:
        close = difflib.get_close_matches(name, known, n=2)
        hints.append(f"{name} (did you mean {' or '.join(close)}?)" if close else name)
    raise TypeError(
        f"{where}: unknown option{'s' if len(unknown) > 1 else ''} {', '.join(hints)}. "
        f"Options: {', '.join(known)}."
    )


class CapturedLog:
    """The text collected by :func:`capture_log`."""

    def __init__(self) -> None:
        self._buffer = io.StringIO()

    @property
    def text(self) -> str:
        """The messages captured so far, one per line (final once the ``with`` block exits)."""
        return self._buffer.getvalue()


def _level_number(level: int | str) -> int:
    if isinstance(level, int):
        return level
    number = logging.getLevelName(level.upper())
    if not isinstance(number, int):
        raise ValueError(f"capture_log: unknown level {level!r}")
    return number


@contextmanager
def capture_log(level: int | str = "INFO", path: str | Path | None = None) -> Iterator[CapturedLog]:
    """Collect swmmbuilder's log and warnings while the ``with`` block runs::

        with swmmbuilder.capture_log(level="INFO", path=out / "build_log.txt") as captured:
            model = sub.to_inp("model.inp", ...)
        print(captured.text)

    A handler is attached to the ``swmmbuilder`` logger (not the root logger) for the block and
    removed on exit, even if the block raises. If the logger's level would drop records at
    ``level``, it is lowered for the block and restored on exit. Unlike
    ``logging.basicConfig``, which does nothing once the root logger has a handler, each use gets
    its own text, so repeated builds in one process (a calibration loop, say) each get their log.
    Other loggers are not captured. Records still propagate to the root logger as usual, so a
    root handler set up by the user sees the records at ``level`` during the block too.

    Warnings of :class:`SwmmBuilderWarning` and its subclasses shown during the block are added
    to the text as well, and are still shown (or recorded by ``pytest.warns``) as usual. Python
    shows a given warning once per location by default, so a warning repeated at the same place
    appears once; a warning the user's filters ignore or turn into errors does not appear.
    Other warning categories are not captured.

    Each line has the format of ``logging.basicConfig``, ``LEVEL:logger:message``, e.g.
    ``INFO:swmmbuilder.gaps:Filled 3 inverts``. A warning has its category in place of the
    logger, e.g. ``WARNING:DataWarning:P7: adverse slope``.

    Parameters
    ----------
    level : int or str, default "INFO"
        The lowest level captured, as a ``logging`` level or its name (``"DEBUG"``, ``"INFO"``,
        ``"WARNING"``, ...). Warnings count as ``WARNING``.
    path : str or Path, optional
        A file the text is written to (UTF-8) on exit, including after an exception.

    Yields
    ------
    CapturedLog
        Its ``text`` holds the messages captured so far.

    Raises
    ------
    ValueError
        If ``level`` is not a known level name.
    """
    number = _level_number(level)
    captured = CapturedLog()
    handler = logging.StreamHandler(captured._buffer)
    handler.setLevel(number)
    handler.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
    logger = logging.getLogger("swmmbuilder")
    old_level = logger.level
    old_showwarning = warnings.showwarning

    def showwarning(
        message: Warning | str,
        category: type[Warning],
        filename: str,
        lineno: int,
        file: TextIO | None = None,
        line: str | None = None,
    ) -> None:
        if issubclass(category, SwmmBuilderWarning) and logging.WARNING >= number:
            record = logging.makeLogRecord(
                {
                    "name": category.__name__,
                    "levelno": logging.WARNING,
                    "levelname": "WARNING",
                    "msg": str(message),
                }
            )
            handler.handle(record)
        old_showwarning(message, category, filename, lineno, file, line)

    logger.addHandler(handler)
    if logger.getEffectiveLevel() > number:
        logger.setLevel(number)
    warnings.showwarning = showwarning
    try:
        yield captured
    finally:
        # Restore only our own hook, in case the block replaced it in turn.
        if warnings.showwarning is showwarning:
            warnings.showwarning = old_showwarning
        logger.setLevel(old_level)
        logger.removeHandler(handler)
        handler.close()
        if path is not None:
            Path(path).write_text(captured.text, encoding="utf-8")
