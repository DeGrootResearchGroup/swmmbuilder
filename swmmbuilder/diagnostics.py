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
"""

import difflib


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
